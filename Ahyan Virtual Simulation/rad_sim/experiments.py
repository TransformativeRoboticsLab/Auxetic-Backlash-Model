from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np

from .cell_geometry import RADHardwareProfile
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, LoadCase, SimulationResult
from .spring_hinge import solve_spring_hinge_3d


@dataclass(frozen=True)
class SourceCommand:
    cell: tuple[int, int]
    alpha: float = 0.0
    z: float = 0.0


@dataclass(frozen=True)
class ResponseCharacterization:
    commands: tuple[SourceCommand, ...]
    alpha_delta: np.ndarray
    height_delta: np.ndarray
    actuator_influence: np.ndarray
    z_residual: np.ndarray
    alpha_die_off: np.ndarray
    z_die_off: np.ndarray
    alpha_reach: int
    z_reach: int
    effective_alpha_die_off: int
    effective_z_die_off: int
    max_abs_alpha_delta: float
    max_abs_height_delta: float


@dataclass(frozen=True)
class ResponseDecayProfile:
    model: str
    alpha_ratio: float
    z_ratio: float
    alpha_length: float
    z_length: float
    alpha_shells: int
    z_shells: int
    alpha_reach: int
    z_reach: int
    alpha_first: float
    z_first: float
    alpha_last: float
    z_last: float


@dataclass(frozen=True)
class PairCharacterization:
    combined: ResponseCharacterization
    first: ResponseCharacterization
    second: ResponseCharacterization
    alpha_superposition_error: float
    height_superposition_error: float


@dataclass(frozen=True)
class OperatorPairInteraction:
    first_index: int
    second_index: int
    first: SourceCommand
    second: SourceCommand
    manhattan_distance: int
    alpha_superposition_error: float
    height_superposition_error: float
    tolerance: float

    @property
    def max_error(self) -> float:
        return max(self.alpha_superposition_error, self.height_superposition_error)

    @property
    def nonadditive(self) -> bool:
        return self.max_error > self.tolerance


@dataclass(frozen=True)
class OperatorInteractionGraph:
    commands: tuple[SourceCommand, ...]
    interactions: tuple[OperatorPairInteraction, ...]
    alpha_error_matrix: np.ndarray
    height_error_matrix: np.ndarray
    interaction_hotspot_map: np.ndarray
    interaction_degree_map: np.ndarray
    total_pair_count: int
    truncated: bool
    tolerance: float

    @property
    def evaluated_pair_count(self) -> int:
        return len(self.interactions)

    @property
    def nonadditive_pair_count(self) -> int:
        return sum(1 for interaction in self.interactions if interaction.nonadditive)

    @property
    def max_alpha_error(self) -> float:
        if self.alpha_error_matrix.size == 0:
            return 0.0
        return float(np.max(self.alpha_error_matrix))

    @property
    def max_height_error(self) -> float:
        if self.height_error_matrix.size == 0:
            return 0.0
        return float(np.max(self.height_error_matrix))

    @property
    def max_interaction_error(self) -> float:
        return max(self.max_alpha_error, self.max_height_error)

    @property
    def max_hotspot_error(self) -> float:
        if self.interaction_hotspot_map.size == 0:
            return 0.0
        return float(np.max(self.interaction_hotspot_map))

    @property
    def max_interaction_degree(self) -> int:
        if self.interaction_degree_map.size == 0:
            return 0
        return int(np.max(self.interaction_degree_map))

    @property
    def interaction_density(self) -> float:
        if self.evaluated_pair_count == 0:
            return 0.0
        return self.nonadditive_pair_count / self.evaluated_pair_count


@dataclass(frozen=True)
class PhysicalResponseComparison:
    commands: tuple[SourceCommand, ...]
    locked_cells: tuple[tuple[int, int], ...]
    kinematic: ResponseCharacterization
    physical_baseline: SimulationResult
    physical_result: SimulationResult
    physical_alpha_delta: np.ndarray
    physical_height_delta: np.ndarray
    physical_center_delta: np.ndarray
    alpha_rms_error: float
    height_rms_error: float
    center_rms_error: float
    max_abs_height_error: float
    max_abs_center_error: float
    physical_energy: float
    physical_success: bool


@dataclass(frozen=True)
class PhysicalPairComparison:
    combined: PhysicalResponseComparison
    first: PhysicalResponseComparison
    second: PhysicalResponseComparison
    physical_height_superposition_error: float
    physical_center_superposition_error: float


@dataclass(frozen=True)
class ResponseMatrix:
    commands: tuple[SourceCommand, ...]
    alpha: np.ndarray
    height: np.ndarray
    cell_shape: tuple[int, int]

    @property
    def alpha_rank(self) -> int:
        return int(np.linalg.matrix_rank(self.alpha))

    @property
    def height_rank(self) -> int:
        return int(np.linalg.matrix_rank(self.height))

    def reachable_alpha_cells(self, tolerance: float = 1e-9) -> int:
        return int(np.count_nonzero(np.any(np.abs(self.alpha) > tolerance, axis=1)))

    def reachable_height_cells(self, tolerance: float = 1e-9) -> int:
        return int(np.count_nonzero(np.any(np.abs(self.height) > tolerance, axis=1)))


PROTOCOL_MEASUREMENT_FIELDS: tuple[str, ...] = (
    "alpha_delta_grid",
    "height_delta_grid",
    "center_displacement_grid",
    "actuator_command",
    "lock_state",
    "pin_hole_slip_mm",
    "actuator_force_n",
)


@dataclass(frozen=True)
class CalibrationExperimentStep:
    id: str
    scope: Literal["single", "pair", "cluster", "lock"]
    commands: tuple[SourceCommand, ...]
    observation_cells: tuple[tuple[int, int], ...]
    locked_cells: tuple[tuple[int, int], ...]
    measurement_fields: tuple[str, ...]
    purpose: str
    expected_response: str
    repeat_count: int = 3

    def to_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "scope": self.scope,
            "commands": [
                {"row": row, "col": col, "alpha": command.alpha, "z": command.z}
                for command in self.commands
                for row, col in (command.cell,)
            ],
            "observationCells": [
                {"row": row, "col": col} for row, col in self.observation_cells
            ],
            "lockedCells": [
                {"row": row, "col": col} for row, col in self.locked_cells
            ],
            "measurementFields": list(self.measurement_fields),
            "purpose": self.purpose,
            "expectedResponse": self.expected_response,
            "repeatCount": self.repeat_count,
        }


@dataclass(frozen=True)
class CalibrationExperimentProtocol:
    config_shape: tuple[int, int]
    center_cell: tuple[int, int]
    steps: tuple[CalibrationExperimentStep, ...]
    hardware_profile_name: str = "paper-reference"
    schema: str = "rad-sim.calibration-experiment-protocol.v1"
    notes: str = (
        "Protocol defines repeatable simulator/bench measurements; it does not "
        "claim the current spring-hinge solver is calibrated."
    )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "hardwareProfile": self.hardware_profile_name,
            "grid": {"rows": self.config_shape[0], "cols": self.config_shape[1]},
            "centerCell": {"row": self.center_cell[0], "col": self.center_cell[1]},
            "measurementFields": list(PROTOCOL_MEASUREMENT_FIELDS),
            "notes": self.notes,
            "steps": [step.to_dict() for step in self.steps],
        }


@dataclass(frozen=True)
class CalibrationExperimentSimulation:
    step_id: str
    scope: str
    command_count: int
    alpha_reach: int
    z_reach: int
    max_abs_alpha_delta: float
    max_abs_height_delta: float
    physical_height_rms_error: float | None = None
    physical_center_rms_error: float | None = None
    physical_success: bool | None = None


@dataclass(frozen=True)
class CalibrationCellMeasurement:
    cell: tuple[int, int]
    alpha_delta: float | None = None
    height_delta: float | None = None
    center_delta: tuple[float, float, float] | None = None
    pin_hole_slip_mm: float | None = None
    actuator_force_n: float | None = None

    def to_dict(self) -> dict[str, object]:
        row, col = self.cell
        return {
            "row": row,
            "col": col,
            "alphaDelta": self.alpha_delta,
            "heightDelta": self.height_delta,
            "centerDelta": list(self.center_delta) if self.center_delta is not None else None,
            "pinHoleSlipMm": self.pin_hole_slip_mm,
            "actuatorForceN": self.actuator_force_n,
        }


@dataclass(frozen=True)
class CalibrationStepMeasurement:
    step_id: str
    repeat_index: int
    cells: tuple[CalibrationCellMeasurement, ...]
    notes: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "stepId": self.step_id,
            "repeatIndex": self.repeat_index,
            "cells": [cell.to_dict() for cell in self.cells],
            "notes": self.notes,
        }


@dataclass(frozen=True)
class CalibrationExperimentMeasurements:
    protocol_schema: str
    hardware_profile_name: str
    steps: tuple[CalibrationStepMeasurement, ...]
    schema: str = "rad-sim.calibration-experiment-results.v1"
    notes: str = "Fill optional measured fields with real bench measurements."

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "protocolSchema": self.protocol_schema,
            "hardwareProfile": self.hardware_profile_name,
            "notes": self.notes,
            "steps": [step.to_dict() for step in self.steps],
        }


@dataclass(frozen=True)
class CalibrationExperimentComparison:
    step_id: str
    repeat_index: int
    measured_cell_count: int
    missing_observation_count: int
    alpha_rmse: float | None
    height_rmse: float | None
    center_rmse: float | None
    max_abs_height_error: float | None
    mean_actuator_force_n: float | None
    mean_pin_hole_slip_mm: float | None
    mean_signed_alpha_error: float | None = None
    mean_signed_height_error: float | None = None
    mean_abs_alpha_error: float | None = None
    mean_abs_height_error: float | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "stepId": self.step_id,
            "repeatIndex": self.repeat_index,
            "measuredCellCount": self.measured_cell_count,
            "missingObservationCount": self.missing_observation_count,
            "alphaRmse": self.alpha_rmse,
            "heightRmse": self.height_rmse,
            "centerRmse": self.center_rmse,
            "maxAbsHeightError": self.max_abs_height_error,
            "meanSignedAlphaError": self.mean_signed_alpha_error,
            "meanSignedHeightError": self.mean_signed_height_error,
            "meanAbsAlphaError": self.mean_abs_alpha_error,
            "meanAbsHeightError": self.mean_abs_height_error,
            "meanActuatorForceN": self.mean_actuator_force_n,
            "meanPinHoleSlipMm": self.mean_pin_hole_slip_mm,
        }


def _state_with_commands(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    locked_cells: Iterable[tuple[int, int]] = (),
) -> LatticeState:
    state = LatticeState.uniform(config)
    for r, c in locked_cells:
        state.locked_mask[r, c] = True
    for command in commands:
        r, c = command.cell
        state.actuator_grid[r, c] += command.alpha
        state.z_actuator_grid[r, c] += command.z
    return state


def _effective_die_off(delta: np.ndarray, die_off: np.ndarray, tolerance: float) -> int:
    mask = np.abs(delta) > tolerance
    if not np.any(mask):
        return 0
    finite = die_off[mask & np.isfinite(die_off)]
    if finite.size == 0:
        return 0
    return int(np.max(finite))


def _reach(delta: np.ndarray, tolerance: float) -> int:
    return int(np.count_nonzero(np.abs(delta) > tolerance))


def _source_cells_from_commands(
    commands: Iterable[SourceCommand],
) -> tuple[tuple[int, int], ...]:
    cells: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for command in commands:
        cell = (int(command.cell[0]), int(command.cell[1]))
        if cell in seen:
            continue
        seen.add(cell)
        cells.append(cell)
    return tuple(cells)


def _nearest_source_distance(
    row: int,
    col: int,
    source_cells: tuple[tuple[int, int], ...],
) -> int:
    return min(abs(row - r) + abs(col - c) for r, c in source_cells)


def _clamp_cell(config: LatticeConfig, cell: tuple[int, int]) -> tuple[int, int]:
    return (
        min(max(int(cell[0]), 0), config.rows - 1),
        min(max(int(cell[1]), 0), config.cols - 1),
    )


def _unique_cells(cells: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    unique: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for cell in cells:
        normalized = (int(cell[0]), int(cell[1]))
        if normalized in seen:
            continue
        seen.add(normalized)
        unique.append(normalized)
    return tuple(unique)


def _neighbor_cells(config: LatticeConfig, center: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    row, col = center
    return _unique_cells(
        _clamp_cell(config, cell)
        for cell in (
            (row, col + 1),
            (row + 1, col),
            (row, col - 1),
            (row - 1, col),
        )
    )


def _fit_shell_decay(
    shells: dict[int, dict[str, float]],
    key: str,
    tolerance: float,
) -> tuple[float, float, int, int, float, float]:
    values = [
        (distance, metrics[key])
        for distance, metrics in sorted(shells.items())
        if metrics[key] > tolerance
    ]
    if not values:
        return 0.0, 0.0, 0, 0, 0.0, 0.0

    shell_count = len(values)
    reach = int(values[-1][0])
    first = float(values[0][1])
    last = float(values[-1][1])
    if shell_count < 2:
        return 0.0, 0.0, shell_count, reach, first, last

    distances = np.array([distance for distance, _ in values], dtype=float)
    magnitudes = np.array([value for _, value in values], dtype=float)
    slope = float(np.polyfit(distances, np.log(magnitudes), 1)[0])
    ratio = float(np.exp(slope))
    length = float(-1.0 / slope) if slope < -tolerance else 0.0
    return ratio, length, shell_count, reach, first, last


def response_decay_profile(
    response: ResponseCharacterization,
    source_cells: Iterable[tuple[int, int]] | None = None,
    tolerance: float = 1e-9,
) -> ResponseDecayProfile:
    """Fit shellwise response attenuation from active source cells.

    This mirrors the browser locality diagnostic: cells are grouped by Manhattan
    distance from the nearest source, each shell uses its maximum absolute
    response, and a log-linear fit estimates the per-shell decay ratio.
    """

    sources = (
        _source_cells_from_commands(response.commands)
        if source_cells is None
        else tuple((int(r), int(c)) for r, c in source_cells)
    )
    if not sources:
        return ResponseDecayProfile(
            model="log-linear shell max",
            alpha_ratio=0.0,
            z_ratio=0.0,
            alpha_length=0.0,
            z_length=0.0,
            alpha_shells=0,
            z_shells=0,
            alpha_reach=0,
            z_reach=0,
            alpha_first=0.0,
            z_first=0.0,
            alpha_last=0.0,
            z_last=0.0,
        )

    shells: dict[int, dict[str, float]] = {}
    rows, cols = response.alpha_delta.shape
    for row in range(rows):
        for col in range(cols):
            distance = _nearest_source_distance(row, col, sources)
            shell = shells.setdefault(distance, {"alpha": 0.0, "z": 0.0})
            shell["alpha"] = max(shell["alpha"], abs(float(response.alpha_delta[row, col])))
            shell["z"] = max(shell["z"], abs(float(response.height_delta[row, col])))

    (
        alpha_ratio,
        alpha_length,
        alpha_shells,
        alpha_reach,
        alpha_first,
        alpha_last,
    ) = _fit_shell_decay(shells, "alpha", tolerance)
    (
        z_ratio,
        z_length,
        z_shells,
        z_reach,
        z_first,
        z_last,
    ) = _fit_shell_decay(shells, "z", tolerance)
    return ResponseDecayProfile(
        model="log-linear shell max",
        alpha_ratio=alpha_ratio,
        z_ratio=z_ratio,
        alpha_length=alpha_length,
        z_length=z_length,
        alpha_shells=alpha_shells,
        z_shells=z_shells,
        alpha_reach=alpha_reach,
        z_reach=z_reach,
        alpha_first=alpha_first,
        z_first=z_first,
        alpha_last=alpha_last,
        z_last=z_last,
    )


def characterize_response(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    locked_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
) -> ResponseCharacterization:
    commands = tuple(commands)
    locked_cells = tuple(locked_cells)
    baseline = simulate_kinematic(config, _state_with_commands(config, (), locked_cells))
    result = simulate_kinematic(config, _state_with_commands(config, commands, locked_cells))
    alpha_delta = result.alpha - baseline.alpha
    height_delta = result.metadata["height"] - baseline.metadata["height"]
    alpha_die_off = result.metadata["die_off"]
    z_die_off = result.metadata["z_die_off"]
    return ResponseCharacterization(
        commands=commands,
        alpha_delta=alpha_delta,
        height_delta=height_delta,
        actuator_influence=result.metadata["actuator_influence"],
        z_residual=result.metadata["z_residual"],
        alpha_die_off=alpha_die_off,
        z_die_off=z_die_off,
        alpha_reach=_reach(alpha_delta, tolerance),
        z_reach=_reach(height_delta, tolerance),
        effective_alpha_die_off=_effective_die_off(alpha_delta, alpha_die_off, tolerance),
        effective_z_die_off=_effective_die_off(height_delta, z_die_off, tolerance),
        max_abs_alpha_delta=float(np.max(np.abs(alpha_delta))),
        max_abs_height_delta=float(np.max(np.abs(height_delta))),
    )


def characterize_single_cell(
    config: LatticeConfig,
    cell: tuple[int, int],
    alpha: float = 0.0,
    z: float = 0.0,
    locked_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
) -> ResponseCharacterization:
    return characterize_response(
        config,
        (SourceCommand(cell=cell, alpha=alpha, z=z),),
        locked_cells=locked_cells,
        tolerance=tolerance,
    )


def characterize_pair(
    config: LatticeConfig,
    first: SourceCommand,
    second: SourceCommand,
    locked_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
) -> PairCharacterization:
    first_response = characterize_response(config, (first,), locked_cells, tolerance)
    second_response = characterize_response(config, (second,), locked_cells, tolerance)
    combined = characterize_response(config, (first, second), locked_cells, tolerance)
    alpha_expected = first_response.alpha_delta + second_response.alpha_delta
    height_expected = first_response.height_delta + second_response.height_delta
    return PairCharacterization(
        combined=combined,
        first=first_response,
        second=second_response,
        alpha_superposition_error=float(np.max(np.abs(combined.alpha_delta - alpha_expected))),
        height_superposition_error=float(np.max(np.abs(combined.height_delta - height_expected))),
    )


def characterize_pairwise_interactions(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    locked_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
    max_pairs: int | None = None,
) -> OperatorInteractionGraph:
    command_tuple = tuple(commands)
    command_count = len(command_tuple)
    total_pair_count = command_count * (command_count - 1) // 2
    alpha_error_matrix = np.zeros((command_count, command_count), dtype=float)
    height_error_matrix = np.zeros((command_count, command_count), dtype=float)
    interaction_hotspot_map = np.zeros((config.rows, config.cols), dtype=float)
    interaction_degree_map = np.zeros((config.rows, config.cols), dtype=int)
    if command_count < 2:
        return OperatorInteractionGraph(
            commands=command_tuple,
            interactions=(),
            alpha_error_matrix=alpha_error_matrix,
            height_error_matrix=height_error_matrix,
            interaction_hotspot_map=interaction_hotspot_map,
            interaction_degree_map=interaction_degree_map,
            total_pair_count=total_pair_count,
            truncated=False,
            tolerance=tolerance,
        )

    max_pairs = total_pair_count if max_pairs is None else max(0, int(max_pairs))
    singles = [
        characterize_response(config, (command,), locked_cells, tolerance)
        for command in command_tuple
    ]
    interactions: list[OperatorPairInteraction] = []
    evaluated = 0
    for first_index in range(command_count):
        for second_index in range(first_index + 1, command_count):
            if evaluated >= max_pairs:
                return OperatorInteractionGraph(
                    commands=command_tuple,
                    interactions=tuple(interactions),
                    alpha_error_matrix=alpha_error_matrix,
                    height_error_matrix=height_error_matrix,
                    interaction_hotspot_map=interaction_hotspot_map,
                    interaction_degree_map=interaction_degree_map,
                    total_pair_count=total_pair_count,
                    truncated=True,
                    tolerance=tolerance,
                )
            first = command_tuple[first_index]
            second = command_tuple[second_index]
            combined = characterize_response(
                config,
                (first, second),
                locked_cells,
                tolerance,
            )
            expected_alpha = (
                singles[first_index].alpha_delta + singles[second_index].alpha_delta
            )
            expected_height = (
                singles[first_index].height_delta + singles[second_index].height_delta
            )
            alpha_error = float(np.max(np.abs(combined.alpha_delta - expected_alpha)))
            height_error = float(np.max(np.abs(combined.height_delta - expected_height)))
            alpha_error_matrix[first_index, second_index] = alpha_error
            alpha_error_matrix[second_index, first_index] = alpha_error
            height_error_matrix[first_index, second_index] = height_error
            height_error_matrix[second_index, first_index] = height_error
            max_error = max(alpha_error, height_error)
            for command in (first, second):
                row, col = command.cell
                if 0 <= row < config.rows and 0 <= col < config.cols:
                    interaction_hotspot_map[row, col] = max(
                        interaction_hotspot_map[row, col],
                        max_error,
                    )
            if max_error > tolerance:
                for row, col in {first.cell, second.cell}:
                    if 0 <= row < config.rows and 0 <= col < config.cols:
                        interaction_degree_map[row, col] += 1
            distance = abs(first.cell[0] - second.cell[0]) + abs(first.cell[1] - second.cell[1])
            interactions.append(
                OperatorPairInteraction(
                    first_index=first_index,
                    second_index=second_index,
                    first=first,
                    second=second,
                    manhattan_distance=int(distance),
                    alpha_superposition_error=alpha_error,
                    height_superposition_error=height_error,
                    tolerance=tolerance,
                )
            )
            evaluated += 1

    return OperatorInteractionGraph(
        commands=command_tuple,
        interactions=tuple(interactions),
        alpha_error_matrix=alpha_error_matrix,
        height_error_matrix=height_error_matrix,
        interaction_hotspot_map=interaction_hotspot_map,
        interaction_degree_map=interaction_degree_map,
        total_pair_count=total_pair_count,
        truncated=False,
        tolerance=tolerance,
    )


def characterize_cluster(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    locked_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
) -> ResponseCharacterization:
    return characterize_response(config, tuple(commands), locked_cells, tolerance)


def build_calibration_experiment_protocol(
    config: LatticeConfig,
    center_cell: tuple[int, int] | None = None,
    *,
    alpha_step: float = -0.25,
    z_step: float = 0.30,
    hardware_profile: RADHardwareProfile | None = None,
    repeat_count: int = 3,
    include_lock_control: bool = True,
) -> CalibrationExperimentProtocol:
    """Build a repeatable single/pair/cluster protocol for physical calibration."""

    center = _clamp_cell(
        config,
        center_cell if center_cell is not None else (config.rows // 2, config.cols // 2),
    )
    neighbors = tuple(cell for cell in _neighbor_cells(config, center) if cell != center)
    primary_neighbor = neighbors[0] if neighbors else center
    secondary_neighbor = neighbors[1] if len(neighbors) > 1 else primary_neighbor
    observation_pair = _unique_cells((center, primary_neighbor))
    observation_cluster = _unique_cells((center, primary_neighbor, secondary_neighbor, *neighbors))
    measurement_fields = PROTOCOL_MEASUREMENT_FIELDS
    alpha_expand = abs(float(alpha_step)) * 0.75
    z_step = float(z_step)
    alpha_step = float(alpha_step)

    steps: list[CalibrationExperimentStep] = [
        CalibrationExperimentStep(
            id="single_alpha_contract",
            scope="single",
            commands=(SourceCommand(center, alpha=alpha_step),),
            observation_cells=(center,),
            locked_cells=(),
            measurement_fields=measurement_fields,
            purpose="Measure the local rotating-square dilation response to contraction.",
            expected_response="Primary alpha change at the commanded cell with backlash-gated neighbor influence.",
            repeat_count=repeat_count,
        ),
        CalibrationExperimentStep(
            id="single_alpha_expand",
            scope="single",
            commands=(SourceCommand(center, alpha=alpha_expand),),
            observation_cells=(center,),
            locked_cells=(),
            measurement_fields=measurement_fields,
            purpose="Measure expansion-side travel and check for asymmetric backlash.",
            expected_response="Positive alpha response at the commanded cell with smaller expansion command.",
            repeat_count=repeat_count,
        ),
        CalibrationExperimentStep(
            id="single_z_lift",
            scope="single",
            commands=(SourceCommand(center, z=z_step),),
            observation_cells=observation_pair,
            locked_cells=(),
            measurement_fields=measurement_fields,
            purpose="Measure direct vertical actuation and residual neighbor lift.",
            expected_response="Commanded cell moves vertically; adjacent observation cell captures pin-hole residual coupling.",
            repeat_count=repeat_count,
        ),
        CalibrationExperimentStep(
            id="pair_z_residual",
            scope="pair",
            commands=(SourceCommand(center, z=z_step),),
            observation_cells=observation_pair,
            locked_cells=(),
            measurement_fields=measurement_fields,
            purpose="Quantify vertical die-off from one actuated cell into a neighboring cell.",
            expected_response="Neighbor height response should decay with clearance, backlash, and graph distance.",
            repeat_count=repeat_count,
        ),
        CalibrationExperimentStep(
            id="pair_superposition",
            scope="pair",
            commands=(
                SourceCommand(center, alpha=alpha_step, z=0.5 * z_step),
                SourceCommand(primary_neighbor, alpha=alpha_step, z=0.5 * z_step),
            ),
            observation_cells=observation_pair,
            locked_cells=(),
            measurement_fields=measurement_fields,
            purpose="Measure whether adjacent cell commands add linearly or interact through backlash.",
            expected_response="Any deviation from summed single-cell responses identifies a programmable-discontinuity interaction.",
            repeat_count=repeat_count,
        ),
        CalibrationExperimentStep(
            id="cluster_mixed_actuation",
            scope="cluster",
            commands=(
                SourceCommand(center, z=z_step),
                SourceCommand(primary_neighbor, alpha=alpha_step),
                SourceCommand(secondary_neighbor, alpha=0.5 * alpha_expand, z=-0.5 * z_step),
            ),
            observation_cells=observation_cluster,
            locked_cells=(),
            measurement_fields=measurement_fields,
            purpose="Measure collective response of mixed horizontal and vertical actuation.",
            expected_response="Cluster field should reveal multi-operator coupling, residual height spread, and reachable directions.",
            repeat_count=repeat_count,
        ),
    ]
    if include_lock_control:
        steps.append(
            CalibrationExperimentStep(
                id="locked_cell_control",
                scope="lock",
                commands=(SourceCommand(center, alpha=alpha_step, z=z_step),),
                observation_cells=observation_pair,
                locked_cells=(center,),
                measurement_fields=measurement_fields,
                purpose="Verify lock enforcement against commanded alpha and vertical motion.",
                expected_response="Locked cell should remain fixed while any neighbor residual exposes compliance leakage.",
                repeat_count=repeat_count,
            )
        )

    return CalibrationExperimentProtocol(
        config_shape=(config.rows, config.cols),
        center_cell=center,
        hardware_profile_name=(
            hardware_profile.name if hardware_profile is not None else "paper-reference"
        ),
        steps=tuple(steps),
    )


def export_calibration_experiment_protocol_json(
    protocol: CalibrationExperimentProtocol,
) -> str:
    return json.dumps(protocol.to_dict(), indent=2)


def calibration_experiment_results_template(
    protocol: CalibrationExperimentProtocol,
) -> CalibrationExperimentMeasurements:
    steps = tuple(
        CalibrationStepMeasurement(
            step_id=step.id,
            repeat_index=repeat_index,
            cells=tuple(
                CalibrationCellMeasurement(cell=cell)
                for cell in step.observation_cells
            ),
            notes="Replace null fields with measured bench data.",
        )
        for step in protocol.steps
        for repeat_index in range(1, step.repeat_count + 1)
    )
    return CalibrationExperimentMeasurements(
        protocol_schema=protocol.schema,
        hardware_profile_name=protocol.hardware_profile_name,
        steps=steps,
    )


def export_calibration_experiment_results_template_json(
    protocol: CalibrationExperimentProtocol,
) -> str:
    return json.dumps(calibration_experiment_results_template(protocol).to_dict(), indent=2)


def _optional_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _measurement_cell_from_dict(raw: dict[str, object]) -> CalibrationCellMeasurement:
    center_delta_raw = raw.get("centerDelta")
    center_delta: tuple[float, float, float] | None = None
    if isinstance(center_delta_raw, (list, tuple)) and len(center_delta_raw) == 3:
        center_delta = tuple(float(value) for value in center_delta_raw)  # type: ignore[assignment]
    return CalibrationCellMeasurement(
        cell=(int(raw["row"]), int(raw["col"])),
        alpha_delta=_optional_float(raw.get("alphaDelta")),
        height_delta=_optional_float(raw.get("heightDelta")),
        center_delta=center_delta,
        pin_hole_slip_mm=_optional_float(raw.get("pinHoleSlipMm")),
        actuator_force_n=_optional_float(raw.get("actuatorForceN")),
    )


def calibration_experiment_measurements_from_dict(
    raw: dict[str, object],
) -> CalibrationExperimentMeasurements:
    schema = str(raw.get("schema", ""))
    if schema != "rad-sim.calibration-experiment-results.v1":
        raise ValueError("unsupported calibration experiment results schema")
    steps_raw = raw.get("steps")
    if not isinstance(steps_raw, list):
        raise ValueError("calibration experiment results require a steps list")
    steps: list[CalibrationStepMeasurement] = []
    for step_raw in steps_raw:
        if not isinstance(step_raw, dict):
            raise ValueError("each calibration result step must be an object")
        cells_raw = step_raw.get("cells", [])
        if not isinstance(cells_raw, list):
            raise ValueError("each calibration result step requires a cells list")
        steps.append(
            CalibrationStepMeasurement(
                step_id=str(step_raw["stepId"]),
                repeat_index=int(step_raw.get("repeatIndex", 1)),
                cells=tuple(
                    _measurement_cell_from_dict(cell)
                    for cell in cells_raw
                    if isinstance(cell, dict)
                ),
                notes=str(step_raw.get("notes", "")),
            )
        )
    return CalibrationExperimentMeasurements(
        protocol_schema=str(raw.get("protocolSchema", "")),
        hardware_profile_name=str(raw.get("hardwareProfile", "")),
        steps=tuple(steps),
        notes=str(raw.get("notes", "")),
    )


def calibration_experiment_measurements_from_json(
    text: str,
) -> CalibrationExperimentMeasurements:
    return calibration_experiment_measurements_from_dict(json.loads(text))


def _rmse(values: list[float]) -> float | None:
    if not values:
        return None
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array**2)))


def _mean_optional(values: list[float]) -> float | None:
    if not values:
        return None
    return float(np.mean(np.asarray(values, dtype=float)))


CalibrationPair = tuple[int, int, float, float]


def _calibration_linear_fit(pairs: list[CalibrationPair]) -> dict[str, object]:
    sample_count = len(pairs)
    if not pairs:
        return {
            "sampleCount": 0,
            "gain": None,
            "bias": None,
            "suggestedGain": None,
            "suggestedBias": None,
            "gainIdentifiable": False,
            "rmsRawError": None,
            "rmsResidual": None,
            "meanPredicted": None,
            "meanMeasured": None,
            "predictedRange": None,
            "measuredRange": None,
        }
    predicted = np.asarray([pair[2] for pair in pairs], dtype=float)
    measured = np.asarray([pair[3] for pair in pairs], dtype=float)
    mean_predicted = float(np.mean(predicted))
    mean_measured = float(np.mean(measured))
    variance = float(np.sum((predicted - mean_predicted) ** 2))
    covariance = float(np.sum((predicted - mean_predicted) * (measured - mean_measured)))
    gain_identifiable = variance > 1e-12
    gain = covariance / variance if gain_identifiable else None
    suggested_gain = gain if gain_identifiable else 1.0
    bias = (
        mean_measured - suggested_gain * mean_predicted
        if gain_identifiable
        else float(np.mean(measured - predicted))
    )
    residuals = measured - (suggested_gain * predicted + bias)
    raw_errors = measured - predicted
    return {
        "sampleCount": sample_count,
        "gain": gain,
        "bias": bias,
        "suggestedGain": suggested_gain,
        "suggestedBias": bias,
        "gainIdentifiable": gain_identifiable,
        "rmsRawError": float(np.sqrt(np.mean(raw_errors**2))),
        "rmsResidual": float(np.sqrt(np.mean(residuals**2))),
        "meanPredicted": mean_predicted,
        "meanMeasured": mean_measured,
        "predictedRange": [float(np.min(predicted)), float(np.max(predicted))],
        "measuredRange": [float(np.min(measured)), float(np.max(measured))],
    }


def _calibration_error_field(
    config: LatticeConfig,
    alpha_pairs: list[CalibrationPair],
    height_pairs: list[CalibrationPair],
    *,
    fit: dict[str, object] | None = None,
) -> dict[str, object]:
    alpha_error = np.zeros((config.rows, config.cols), dtype=float)
    height_error = np.zeros((config.rows, config.cols), dtype=float)
    sample_count = np.zeros((config.rows, config.cols), dtype=int)
    alpha_sample_count = np.zeros((config.rows, config.cols), dtype=int)
    height_sample_count = np.zeros((config.rows, config.cols), dtype=int)
    touched: set[tuple[int, int]] = set()
    alpha_fit = fit.get("alpha", {}) if isinstance(fit, dict) else {}
    height_fit = fit.get("height", {}) if isinstance(fit, dict) else {}
    alpha_gain = float(alpha_fit.get("suggestedGain") or 1.0)
    alpha_bias = float(alpha_fit.get("suggestedBias") or 0.0)
    height_gain = float(height_fit.get("suggestedGain") or 1.0)
    height_bias = float(height_fit.get("suggestedBias") or 0.0)
    for row, col, predicted, measured in alpha_pairs:
        error = measured - (alpha_gain * predicted + alpha_bias)
        alpha_error[row, col] += error
        alpha_sample_count[row, col] += 1
        touched.add((row, col))
    for row, col, predicted, measured in height_pairs:
        error = measured - (height_gain * predicted + height_bias)
        height_error[row, col] += error
        height_sample_count[row, col] += 1
        touched.add((row, col))
    for row, col in touched:
        sample_count[row, col] += 1
    nonzero_alpha = alpha_sample_count > 0
    nonzero_height = height_sample_count > 0
    alpha_error[nonzero_alpha] = alpha_error[nonzero_alpha] / alpha_sample_count[nonzero_alpha]
    height_error[nonzero_height] = height_error[nonzero_height] / height_sample_count[nonzero_height]
    combined = np.hypot(alpha_error, height_error)
    max_combined = float(np.max(combined)) if combined.size else 0.0
    worst_cell = None
    top_cells: list[dict[str, object]] = []
    if np.any(sample_count > 0):
        measured_indices = np.argwhere(sample_count > 0)
        measured_errors = np.asarray(
            [combined[row, col] for row, col in measured_indices],
            dtype=float,
        )
        row, col = measured_indices[int(np.argmax(measured_errors))]
        worst_cell = {
            "row": int(row),
            "col": int(col),
            "alphaError": float(alpha_error[row, col]),
            "heightError": float(height_error[row, col]),
            "combinedError": float(combined[row, col]),
            "sampleCount": int(sample_count[row, col]),
            "alphaSampleCount": int(alpha_sample_count[row, col]),
            "heightSampleCount": int(height_sample_count[row, col]),
        }
        for row, col in np.argwhere(sample_count > 0):
            top_cells.append(
                {
                    "row": int(row),
                    "col": int(col),
                    "alphaError": float(alpha_error[row, col]),
                    "heightError": float(height_error[row, col]),
                    "combinedError": float(combined[row, col]),
                    "sampleCount": int(sample_count[row, col]),
                    "alphaSampleCount": int(alpha_sample_count[row, col]),
                    "heightSampleCount": int(height_sample_count[row, col]),
                }
            )
        top_cells.sort(
            key=lambda cell: (
                -float(cell["combinedError"]),
                int(cell["row"]),
                int(cell["col"]),
            )
        )
    return {
        "alphaError": alpha_error.tolist(),
        "heightError": height_error.tolist(),
        "combinedError": combined.tolist(),
        "sampleCount": sample_count.tolist(),
        "alphaSampleCount": alpha_sample_count.tolist(),
        "heightSampleCount": height_sample_count.tolist(),
        "maxAbsAlphaError": float(np.max(np.abs(alpha_error))) if alpha_error.size else 0.0,
        "maxAbsHeightError": float(np.max(np.abs(height_error))) if height_error.size else 0.0,
        "maxCombinedError": max_combined,
        "worstCell": worst_cell,
        "topCells": top_cells[:12],
    }


def _summarize_calibration_comparison(
    comparisons: tuple[CalibrationExperimentComparison, ...],
    field: dict[str, object],
    fit: dict[str, object],
    fit_residual_field: dict[str, object],
) -> dict[str, object]:
    height_rmse = [item.height_rmse for item in comparisons if item.height_rmse is not None]
    alpha_rmse = [item.alpha_rmse for item in comparisons if item.alpha_rmse is not None]
    center_rmse = [item.center_rmse for item in comparisons if item.center_rmse is not None]
    signed_alpha = [
        item.mean_signed_alpha_error
        for item in comparisons
        if item.mean_signed_alpha_error is not None
    ]
    signed_height = [
        item.mean_signed_height_error
        for item in comparisons
        if item.mean_signed_height_error is not None
    ]
    abs_alpha = [item.mean_abs_alpha_error for item in comparisons if item.mean_abs_alpha_error is not None]
    abs_height = [item.mean_abs_height_error for item in comparisons if item.mean_abs_height_error is not None]
    force = [item.mean_actuator_force_n for item in comparisons if item.mean_actuator_force_n is not None]
    slip = [item.mean_pin_hole_slip_mm for item in comparisons if item.mean_pin_hole_slip_mm is not None]
    worst = max(
        (item for item in comparisons if item.max_abs_height_error is not None),
        key=lambda item: abs(float(item.max_abs_height_error)),
        default=None,
    )
    return {
        "schema": "rad-sim.calibration-experiment-comparison-summary.v1",
        "stepCount": len(comparisons),
        "measuredCellCount": sum(item.measured_cell_count for item in comparisons),
        "missingObservationCount": sum(item.missing_observation_count for item in comparisons),
        "alphaRmseMean": _mean_optional(alpha_rmse),
        "heightRmseMean": _mean_optional(height_rmse),
        "centerRmseMean": _mean_optional(center_rmse),
        "meanSignedAlphaError": _mean_optional(signed_alpha),
        "meanSignedHeightError": _mean_optional(signed_height),
        "meanAbsAlphaError": _mean_optional(abs_alpha),
        "meanAbsHeightError": _mean_optional(abs_height),
        "fit": fit,
        "fitResidualMaxCombinedError": fit_residual_field.get("maxCombinedError"),
        "fitResidualWorstCell": fit_residual_field.get("worstCell"),
        "topCells": field.get("topCells", []),
        "fitResidualTopCells": fit_residual_field.get("topCells", []),
        "maxAbsHeightError": worst.max_abs_height_error if worst else None,
        "maxAbsAlphaError": field.get("maxAbsAlphaError"),
        "maxCombinedError": field.get("maxCombinedError"),
        "worstCell": field.get("worstCell"),
        "worstStepId": worst.step_id if worst else None,
        "meanActuatorForceN": _mean_optional(force),
        "meanPinHoleSlipMm": _mean_optional(slip),
    }


def compare_calibration_experiment_measurements(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    measurements: CalibrationExperimentMeasurements,
    *,
    tolerance: float = 1e-9,
) -> tuple[CalibrationExperimentComparison, ...]:
    """Compare imported bench measurements against current simulated protocol response."""

    protocol_steps = {step.id: step for step in protocol.steps}
    comparisons: list[CalibrationExperimentComparison] = []
    for measured_step in measurements.steps:
        step = protocol_steps.get(measured_step.step_id)
        if step is None:
            continue
        response = characterize_response(
            config,
            step.commands,
            locked_cells=step.locked_cells,
            tolerance=tolerance,
        )
        baseline = simulate_kinematic(
            config,
            _state_with_commands(config, (), step.locked_cells),
        )
        result = simulate_kinematic(
            config,
            _state_with_commands(config, step.commands, step.locked_cells),
        )
        center_delta = result.deformed_centers_3d - baseline.deformed_centers_3d
        observed = {measurement.cell: measurement for measurement in measured_step.cells}
        alpha_errors: list[float] = []
        height_errors: list[float] = []
        center_errors: list[float] = []
        slip_values: list[float] = []
        force_values: list[float] = []
        for measurement in measured_step.cells:
            row, col = measurement.cell
            if not (0 <= row < config.rows and 0 <= col < config.cols):
                continue
            if measurement.alpha_delta is not None:
                alpha_errors.append(
                    measurement.alpha_delta - float(response.alpha_delta[row, col])
                )
            if measurement.height_delta is not None:
                height_errors.append(
                    measurement.height_delta - float(response.height_delta[row, col])
                )
            if measurement.center_delta is not None:
                simulated_center = center_delta[row, col, :]
                center_errors.append(
                    float(
                        np.linalg.norm(
                            np.asarray(measurement.center_delta, dtype=float)
                            - simulated_center
                        )
                    )
                )
            if measurement.pin_hole_slip_mm is not None:
                slip_values.append(measurement.pin_hole_slip_mm)
            if measurement.actuator_force_n is not None:
                force_values.append(measurement.actuator_force_n)
        missing = sum(1 for cell in step.observation_cells if cell not in observed)
        comparisons.append(
            CalibrationExperimentComparison(
                step_id=measured_step.step_id,
                repeat_index=measured_step.repeat_index,
                measured_cell_count=len(observed),
                missing_observation_count=missing,
                alpha_rmse=_rmse(alpha_errors),
                height_rmse=_rmse(height_errors),
                center_rmse=_rmse(center_errors),
                max_abs_height_error=(
                    max(abs(value) for value in height_errors) if height_errors else None
                ),
                mean_actuator_force_n=_mean_optional(force_values),
                mean_pin_hole_slip_mm=_mean_optional(slip_values),
                mean_signed_alpha_error=_mean_optional(alpha_errors),
                mean_signed_height_error=_mean_optional(height_errors),
                mean_abs_alpha_error=_mean_optional([abs(value) for value in alpha_errors]),
                mean_abs_height_error=_mean_optional([abs(value) for value in height_errors]),
            )
        )
    return tuple(comparisons)


def calibration_experiment_comparison_report(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    measurements: CalibrationExperimentMeasurements,
    *,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Build a JSON-ready simulator-vs-bench calibration comparison report."""

    protocol_steps = {step.id: step for step in protocol.steps}
    comparisons: list[CalibrationExperimentComparison] = []
    alpha_pairs: list[CalibrationPair] = []
    height_pairs: list[CalibrationPair] = []
    for measured_step in measurements.steps:
        step = protocol_steps.get(measured_step.step_id)
        if step is None:
            continue
        response = characterize_response(
            config,
            step.commands,
            locked_cells=step.locked_cells,
            tolerance=tolerance,
        )
        baseline = simulate_kinematic(
            config,
            _state_with_commands(config, (), step.locked_cells),
        )
        result = simulate_kinematic(
            config,
            _state_with_commands(config, step.commands, step.locked_cells),
        )
        center_delta = result.deformed_centers_3d - baseline.deformed_centers_3d
        observed = {measurement.cell: measurement for measurement in measured_step.cells}
        alpha_errors: list[float] = []
        height_errors: list[float] = []
        center_errors: list[float] = []
        slip_values: list[float] = []
        force_values: list[float] = []
        for measurement in measured_step.cells:
            row, col = measurement.cell
            if not (0 <= row < config.rows and 0 <= col < config.cols):
                continue
            if measurement.alpha_delta is not None:
                predicted = float(response.alpha_delta[row, col])
                alpha_errors.append(measurement.alpha_delta - predicted)
                alpha_pairs.append((row, col, predicted, measurement.alpha_delta))
            if measurement.height_delta is not None:
                predicted = float(response.height_delta[row, col])
                height_errors.append(measurement.height_delta - predicted)
                height_pairs.append((row, col, predicted, measurement.height_delta))
            if measurement.center_delta is not None:
                simulated_center = center_delta[row, col, :]
                center_errors.append(
                    float(
                        np.linalg.norm(
                            np.asarray(measurement.center_delta, dtype=float)
                            - simulated_center
                        )
                    )
                )
            if measurement.pin_hole_slip_mm is not None:
                slip_values.append(measurement.pin_hole_slip_mm)
            if measurement.actuator_force_n is not None:
                force_values.append(measurement.actuator_force_n)
        missing = sum(1 for cell in step.observation_cells if cell not in observed)
        comparisons.append(
            CalibrationExperimentComparison(
                step_id=measured_step.step_id,
                repeat_index=measured_step.repeat_index,
                measured_cell_count=len(observed),
                missing_observation_count=missing,
                alpha_rmse=_rmse(alpha_errors),
                height_rmse=_rmse(height_errors),
                center_rmse=_rmse(center_errors),
                max_abs_height_error=(
                    max(abs(value) for value in height_errors) if height_errors else None
                ),
                mean_actuator_force_n=_mean_optional(force_values),
                mean_pin_hole_slip_mm=_mean_optional(slip_values),
                mean_signed_alpha_error=_mean_optional(alpha_errors),
                mean_signed_height_error=_mean_optional(height_errors),
                mean_abs_alpha_error=_mean_optional([abs(value) for value in alpha_errors]),
                mean_abs_height_error=_mean_optional([abs(value) for value in height_errors]),
            )
        )
    comparison_tuple = tuple(comparisons)
    fit = {
        "alpha": _calibration_linear_fit(alpha_pairs),
        "height": _calibration_linear_fit(height_pairs),
    }
    field = _calibration_error_field(config, alpha_pairs, height_pairs)
    fit_residual_field = _calibration_error_field(
        config,
        alpha_pairs,
        height_pairs,
        fit=fit,
    )
    summary = _summarize_calibration_comparison(
        comparison_tuple,
        field,
        fit,
        fit_residual_field,
    )
    hardware_profile = getattr(config, "hardware_profile", None)
    profile_name = hardware_profile.name if hardware_profile else None
    return {
        "schema": "rad-sim.calibration-comparison-report.v1",
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "cellSize": config.cell_size,
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
        },
        "hardwareProfile": profile_name,
        "sourceResultsSchema": measurements.schema,
        "summary": summary,
        "comparison": {
            "schema": "rad-sim.calibration-experiment-comparison.v1",
            "comparisons": [comparison.to_dict() for comparison in comparison_tuple],
            "field": field,
            "fit": fit,
            "fitResidualField": fit_residual_field,
        },
    }


def export_calibration_experiment_comparison_json(
    comparisons: Iterable[CalibrationExperimentComparison],
) -> str:
    return json.dumps(
        {
            "schema": "rad-sim.calibration-experiment-comparison.v1",
            "comparisons": [comparison.to_dict() for comparison in comparisons],
        },
        indent=2,
    )


def export_calibration_experiment_comparison_report_json(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    measurements: CalibrationExperimentMeasurements,
    *,
    tolerance: float = 1e-9,
) -> str:
    return json.dumps(
        calibration_experiment_comparison_report(
            config,
            protocol,
            measurements,
            tolerance=tolerance,
        ),
        indent=2,
    )


def run_calibration_experiment_protocol(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    *,
    physical: bool = False,
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
) -> tuple[CalibrationExperimentSimulation, ...]:
    """Run protocol steps through the current simulator for baseline expectations."""

    simulations: list[CalibrationExperimentSimulation] = []
    for step in protocol.steps:
        response = characterize_response(
            config,
            step.commands,
            locked_cells=step.locked_cells,
            tolerance=tolerance,
        )
        physical_height_rms_error: float | None = None
        physical_center_rms_error: float | None = None
        physical_success: bool | None = None
        if physical:
            comparison = compare_physical_response(
                config,
                step.commands,
                locked_cells=step.locked_cells,
                load_case=load_case,
                tolerance=tolerance,
            )
            physical_height_rms_error = comparison.height_rms_error
            physical_center_rms_error = comparison.center_rms_error
            physical_success = comparison.physical_success
        simulations.append(
            CalibrationExperimentSimulation(
                step_id=step.id,
                scope=step.scope,
                command_count=len(step.commands),
                alpha_reach=response.alpha_reach,
                z_reach=response.z_reach,
                max_abs_alpha_delta=response.max_abs_alpha_delta,
                max_abs_height_delta=response.max_abs_height_delta,
                physical_height_rms_error=physical_height_rms_error,
                physical_center_rms_error=physical_center_rms_error,
                physical_success=physical_success,
            )
        )
    return tuple(simulations)


def compare_physical_response(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    locked_cells: Iterable[tuple[int, int]] = (),
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
) -> PhysicalResponseComparison:
    """Compare kinematic backlash propagation against 3D spring-hinge relaxation."""

    commands = tuple(commands)
    locked_cells = tuple(locked_cells)
    kinematic = characterize_response(config, commands, locked_cells, tolerance)
    baseline_state = _state_with_commands(config, (), locked_cells)
    command_state = _state_with_commands(config, commands, locked_cells)
    physical_baseline = solve_spring_hinge_3d(config, baseline_state, load_case)
    physical_result = solve_spring_hinge_3d(config, command_state, load_case)
    physical_alpha_delta = physical_result.alpha - physical_baseline.alpha
    physical_center_delta = (
        physical_result.deformed_centers_3d - physical_baseline.deformed_centers_3d
    )
    physical_height_delta = physical_center_delta[..., 2]
    alpha_error = physical_alpha_delta - kinematic.alpha_delta
    height_error = physical_height_delta - kinematic.height_delta
    kinematic_center_delta = (
        simulate_kinematic(config, command_state).deformed_centers_3d
        - simulate_kinematic(config, baseline_state).deformed_centers_3d
    )
    center_error = physical_center_delta - kinematic_center_delta
    return PhysicalResponseComparison(
        commands=commands,
        locked_cells=locked_cells,
        kinematic=kinematic,
        physical_baseline=physical_baseline,
        physical_result=physical_result,
        physical_alpha_delta=physical_alpha_delta,
        physical_height_delta=physical_height_delta,
        physical_center_delta=physical_center_delta,
        alpha_rms_error=float(np.sqrt(np.mean(alpha_error**2))),
        height_rms_error=float(np.sqrt(np.mean(height_error**2))),
        center_rms_error=float(np.sqrt(np.mean(center_error**2))),
        max_abs_height_error=float(np.max(np.abs(height_error))),
        max_abs_center_error=float(np.max(np.linalg.norm(center_error, axis=2))),
        physical_energy=float(physical_result.metadata.get("energy", np.nan)),
        physical_success=bool(
            physical_baseline.metadata.get("success", False)
            and physical_result.metadata.get("success", False)
        ),
    )


def compare_physical_pair(
    config: LatticeConfig,
    first: SourceCommand,
    second: SourceCommand,
    locked_cells: Iterable[tuple[int, int]] = (),
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
) -> PhysicalPairComparison:
    first_response = compare_physical_response(
        config, (first,), locked_cells, load_case, tolerance
    )
    second_response = compare_physical_response(
        config, (second,), locked_cells, load_case, tolerance
    )
    combined = compare_physical_response(
        config, (first, second), locked_cells, load_case, tolerance
    )
    height_expected = first_response.physical_height_delta + second_response.physical_height_delta
    center_expected = first_response.physical_center_delta + second_response.physical_center_delta
    return PhysicalPairComparison(
        combined=combined,
        first=first_response,
        second=second_response,
        physical_height_superposition_error=float(
            np.max(np.abs(combined.physical_height_delta - height_expected))
        ),
        physical_center_superposition_error=float(
            np.max(np.linalg.norm(combined.physical_center_delta - center_expected, axis=2))
        ),
    )


def compare_physical_cluster(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    locked_cells: Iterable[tuple[int, int]] = (),
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
) -> PhysicalResponseComparison:
    return compare_physical_response(
        config,
        tuple(commands),
        locked_cells=locked_cells,
        load_case=load_case,
        tolerance=tolerance,
    )


def build_response_matrix(
    config: LatticeConfig,
    actuator_cells: Iterable[tuple[int, int]] | None = None,
    *,
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    include_alpha: bool = True,
    include_z: bool = True,
    locked_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
) -> ResponseMatrix:
    if not include_alpha and not include_z:
        raise ValueError("at least one command family must be included")
    cells = tuple(
        actuator_cells
        if actuator_cells is not None
        else ((r, c) for r in range(config.rows) for c in range(config.cols))
    )
    commands: list[SourceCommand] = []
    for cell in cells:
        if include_alpha:
            commands.append(SourceCommand(cell=cell, alpha=alpha_step))
        if include_z:
            commands.append(SourceCommand(cell=cell, z=z_step))

    alpha_columns = []
    height_columns = []
    for command in commands:
        response = characterize_response(config, (command,), locked_cells, tolerance)
        alpha_columns.append(response.alpha_delta.reshape(-1))
        height_columns.append(response.height_delta.reshape(-1))

    if commands:
        alpha_matrix = np.column_stack(alpha_columns)
        height_matrix = np.column_stack(height_columns)
    else:
        alpha_matrix = np.zeros((config.rows * config.cols, 0), dtype=float)
        height_matrix = np.zeros((config.rows * config.cols, 0), dtype=float)
    return ResponseMatrix(
        commands=tuple(commands),
        alpha=alpha_matrix,
        height=height_matrix,
        cell_shape=(config.rows, config.cols),
    )
