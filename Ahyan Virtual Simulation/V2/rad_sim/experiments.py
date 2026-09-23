from __future__ import annotations

import importlib
import importlib.util
import json
from copy import deepcopy
from dataclasses import dataclass, replace
from typing import Iterable, Literal

import numpy as np

from .cell_geometry import RADHardwareProfile
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, LoadCase, SimulationResult
from .operators import (
    ProgrammableDiscontinuityEvent,
    apply_programmable_event,
    lattice_topology_diagnostic,
)
try:
    from .spring_hinge import solve_spring_hinge_3d
except ModuleNotFoundError as exc:  # pragma: no cover - exercised only in minimal runtimes.
    if exc.name not in {"scipy", "scipy.optimize"}:
        raise
    _SPRING_HINGE_IMPORT_ERROR = exc

    def solve_spring_hinge_3d(*args, **kwargs):
        raise ImportError(
            "SciPy is required for spring-hinge mechanics. Install project dependencies "
            "before using solve_spring_hinge_3d."
        ) from _SPRING_HINGE_IMPORT_ERROR


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
    positive_z_reach: int
    negative_z_reach: int
    max_positive_height_delta: float
    max_negative_height_delta: float


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
        if self.alpha.size == 0 or self.alpha.shape[1] == 0:
            return 0
        return int(np.linalg.matrix_rank(self.alpha))

    @property
    def height_rank(self) -> int:
        if self.height.size == 0 or self.height.shape[1] == 0:
            return 0
        return int(np.linalg.matrix_rank(self.height))

    def reachable_alpha_cells(self, tolerance: float = 1e-9) -> int:
        return int(np.count_nonzero(np.any(np.abs(self.alpha) > tolerance, axis=1)))

    def reachable_height_cells(self, tolerance: float = 1e-9) -> int:
        return int(np.count_nonzero(np.any(np.abs(self.height) > tolerance, axis=1)))

    def alpha_underactuated_cells(self, tolerance: float = 1e-9) -> int:
        return max(0, self.cell_shape[0] * self.cell_shape[1] - self.reachable_alpha_cells(tolerance))

    def height_underactuated_cells(self, tolerance: float = 1e-9) -> int:
        return max(0, self.cell_shape[0] * self.cell_shape[1] - self.reachable_height_cells(tolerance))

    def to_dict(self, tolerance: float = 1e-9) -> dict[str, object]:
        rows, cols = self.cell_shape

        def command_family(command: SourceCommand) -> str:
            has_alpha = abs(command.alpha) > tolerance
            has_z = abs(command.z) > tolerance
            if has_alpha and has_z:
                return "combined"
            if has_alpha:
                return "alpha"
            if has_z:
                return "z"
            return "zero"

        return {
            "schema": "rad-sim.response-matrix.v1",
            "grid": {"rows": rows, "cols": cols},
            "cellOrder": [
                {"index": row * cols + col, "row": row, "col": col}
                for row in range(rows)
                for col in range(cols)
            ],
            "commands": [
                {
                    "index": index,
                    "row": command.cell[0],
                    "col": command.cell[1],
                    "family": command_family(command),
                    "alpha": command.alpha,
                    "z": command.z,
                }
                for index, command in enumerate(self.commands)
            ],
            "alpha": self.alpha.tolist(),
            "height": self.height.tolist(),
            "diagnostics": {
                "columnCount": len(self.commands),
                "alphaRank": self.alpha_rank,
                "heightRank": self.height_rank,
                "reachableAlphaCells": self.reachable_alpha_cells(tolerance),
                "reachableHeightCells": self.reachable_height_cells(tolerance),
                "alphaUnderactuatedCells": self.alpha_underactuated_cells(tolerance),
                "heightUnderactuatedCells": self.height_underactuated_cells(tolerance),
                "tolerance": tolerance,
            },
        }


@dataclass(frozen=True)
class RemovedTopologyReachabilityComparison:
    removed_cells: tuple[tuple[int, int], ...]
    actuator_cells: tuple[tuple[int, int], ...]
    locked_cells: tuple[tuple[int, int], ...]
    intact_matrix: ResponseMatrix
    removed_matrix: ResponseMatrix
    intact_topology: dict[str, object]
    removed_topology: dict[str, object]
    tolerance: float

    @property
    def alpha_reachable_cell_loss(self) -> int:
        return max(
            0,
            self.intact_matrix.reachable_alpha_cells(self.tolerance)
            - self.removed_matrix.reachable_alpha_cells(self.tolerance),
        )

    @property
    def height_reachable_cell_loss(self) -> int:
        return max(
            0,
            self.intact_matrix.reachable_height_cells(self.tolerance)
            - self.removed_matrix.reachable_height_cells(self.tolerance),
        )

    @property
    def alpha_rank_loss(self) -> int:
        return max(0, self.intact_matrix.alpha_rank - self.removed_matrix.alpha_rank)

    @property
    def height_rank_loss(self) -> int:
        return max(0, self.intact_matrix.height_rank - self.removed_matrix.height_rank)

    @property
    def component_count_delta(self) -> int:
        return int(self.removed_topology["component_count"]) - int(
            self.intact_topology["component_count"]
        )

    @property
    def deleted_edge_delta(self) -> int:
        return int(self.removed_topology["deleted_edge_count"]) - int(
            self.intact_topology["deleted_edge_count"]
        )

    def to_dict(self, tolerance: float | None = None) -> dict[str, object]:
        tol = self.tolerance if tolerance is None else tolerance
        return {
            "schema": "rad-sim.removed-topology-reachability-comparison.v1",
            "removedCells": [_cell_dict(cell) for cell in self.removed_cells],
            "actuatorCells": [_cell_dict(cell) for cell in self.actuator_cells],
            "lockedCells": [_cell_dict(cell) for cell in self.locked_cells],
            "summary": {
                "alphaReachableCellLoss": self.alpha_reachable_cell_loss,
                "heightReachableCellLoss": self.height_reachable_cell_loss,
                "alphaRankLoss": self.alpha_rank_loss,
                "heightRankLoss": self.height_rank_loss,
                "componentCountDelta": self.component_count_delta,
                "deletedEdgeDelta": self.deleted_edge_delta,
                "tolerance": tol,
            },
            "intact": {
                "responseMatrix": self.intact_matrix.to_dict(tolerance=tol)["diagnostics"],
                "topology": _json_ready_topology(self.intact_topology),
            },
            "removed": {
                "responseMatrix": self.removed_matrix.to_dict(tolerance=tol)["diagnostics"],
                "topology": _json_ready_topology(self.removed_topology),
            },
        }


@dataclass(frozen=True)
class TopologyExperimentScenario:
    name: str
    commands: tuple[SourceCommand, ...] = ()
    locked_cells: tuple[tuple[int, int], ...] = ()
    removed_cells: tuple[tuple[int, int], ...] = ()
    actuator_cells: tuple[tuple[int, int], ...] | None = None

    def normalized(self) -> "TopologyExperimentScenario":
        return TopologyExperimentScenario(
            name=self.name,
            commands=tuple(
                SourceCommand(
                    cell=(int(command.cell[0]), int(command.cell[1])),
                    alpha=float(command.alpha),
                    z=float(command.z),
                )
                for command in self.commands
            ),
            locked_cells=tuple((int(r), int(c)) for r, c in self.locked_cells),
            removed_cells=tuple((int(r), int(c)) for r, c in self.removed_cells),
            actuator_cells=(
                None
                if self.actuator_cells is None
                else tuple((int(r), int(c)) for r, c in self.actuator_cells)
            ),
        )


@dataclass(frozen=True)
class TopologyComponentResponseSummary:
    label: int
    cell_count: int
    actuator_cell_count: int
    command_column_count: int
    alpha_rank: int
    height_rank: int
    reachable_alpha_cells: int
    reachable_height_cells: int

    @property
    def blocked_alpha_cells(self) -> int:
        return max(0, self.cell_count - self.reachable_alpha_cells)

    @property
    def blocked_height_cells(self) -> int:
        return max(0, self.cell_count - self.reachable_height_cells)

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "cellCount": self.cell_count,
            "actuatorCellCount": self.actuator_cell_count,
            "commandColumnCount": self.command_column_count,
            "alphaRank": self.alpha_rank,
            "heightRank": self.height_rank,
            "reachableAlphaCells": self.reachable_alpha_cells,
            "reachableHeightCells": self.reachable_height_cells,
            "blockedAlphaCells": self.blocked_alpha_cells,
            "blockedHeightCells": self.blocked_height_cells,
        }


@dataclass(frozen=True)
class TopologyExperimentScenarioResult:
    scenario: TopologyExperimentScenario
    topology: dict[str, object]
    response_matrix: ResponseMatrix
    response: ResponseCharacterization | None
    component_reachable_cells: int
    component_blocked_cells: int
    positive_height_reachable_cells: int
    negative_height_reachable_cells: int
    component_summaries: tuple[TopologyComponentResponseSummary, ...]
    target_alpha_rms: float | None
    target_height_rms: float | None
    max_abs_target_alpha_residual: float | None
    max_abs_target_height_residual: float | None
    tolerance: float

    def to_dict(self, tolerance: float | None = None) -> dict[str, object]:
        tol = self.tolerance if tolerance is None else tolerance
        response = self.response
        return {
            "name": self.scenario.name,
            "commands": [_source_command_dict(command) for command in self.scenario.commands],
            "lockedCells": [_cell_dict(cell) for cell in self.scenario.locked_cells],
            "removedCells": [_cell_dict(cell) for cell in self.scenario.removed_cells],
            "actuatorCells": (
                None
                if self.scenario.actuator_cells is None
                else [_cell_dict(cell) for cell in self.scenario.actuator_cells]
            ),
            "topology": _json_ready_topology(self.topology),
            "responseMatrix": self.response_matrix.to_dict(tolerance=tol)["diagnostics"],
            "metrics": {
                "componentReachableCells": self.component_reachable_cells,
                "componentBlockedCells": self.component_blocked_cells,
                "componentResponseRank": [
                    summary.to_dict() for summary in self.component_summaries
                ],
                "positiveHeightReachableCells": self.positive_height_reachable_cells,
                "negativeHeightReachableCells": self.negative_height_reachable_cells,
                "alphaReach": 0 if response is None else response.alpha_reach,
                "zReach": 0 if response is None else response.z_reach,
                "effectiveAlphaDieOff": 0 if response is None else response.effective_alpha_die_off,
                "effectiveZDieOff": 0 if response is None else response.effective_z_die_off,
                "targetAlphaRms": self.target_alpha_rms,
                "targetHeightRms": self.target_height_rms,
                "maxAbsTargetAlphaResidual": self.max_abs_target_alpha_residual,
                "maxAbsTargetHeightResidual": self.max_abs_target_height_residual,
                "tolerance": tol,
            },
        }


@dataclass(frozen=True)
class TopologyExperimentReport:
    config_shape: tuple[int, int]
    scenarios: tuple[TopologyExperimentScenarioResult, ...]
    target_alpha: np.ndarray | None = None
    target_height: np.ndarray | None = None
    notes: str = (
        "Finite-response topology experiment for comparing programmable "
        "discontinuity states. This is a simulator diagnostic, not a calibrated "
        "physical validation."
    )

    def to_dict(self, tolerance: float | None = None) -> dict[str, object]:
        return {
            "schema": "rad-sim.topology-experiment-report.v1",
            "grid": {"rows": self.config_shape[0], "cols": self.config_shape[1]},
            "notes": self.notes,
            "target": {
                "alpha": None if self.target_alpha is None else self.target_alpha.tolist(),
                "height": None if self.target_height is None else self.target_height.tolist(),
            },
            "summary": {
                "scenarioCount": len(self.scenarios),
                "maxComponentCount": max(
                    (int(scenario.topology["component_count"]) for scenario in self.scenarios),
                    default=0,
                ),
                "maxDeletedEdges": max(
                    (int(scenario.topology["deleted_edge_count"]) for scenario in self.scenarios),
                    default=0,
                ),
                "maxComponentBlockedCells": max(
                    (scenario.component_blocked_cells for scenario in self.scenarios),
                    default=0,
                ),
                "minHeightReachableCells": min(
                    (
                        scenario.response_matrix.reachable_height_cells(
                            scenario.tolerance
                        )
                        for scenario in self.scenarios
                    ),
                    default=0,
                ),
                "maxComponentHeightRank": max(
                    (
                        component.height_rank
                        for scenario in self.scenarios
                        for component in scenario.component_summaries
                    ),
                    default=0,
                ),
                "maxComponentBlockedHeightCells": max(
                    (
                        component.blocked_height_cells
                        for scenario in self.scenarios
                        for component in scenario.component_summaries
                    ),
                    default=0,
                ),
            },
            "scenarioComparisons": _scenario_comparisons(self.scenarios),
            "scenarios": [
                scenario.to_dict(tolerance=tolerance) for scenario in self.scenarios
            ],
        }


PROTOCOL_MEASUREMENT_FIELDS: tuple[str, ...] = (
    "alpha_delta_grid",
    "height_delta_grid",
    "center_displacement_grid",
    "actuator_command",
    "lock_state",
    "pin_hole_slip_mm",
    "actuator_force_n",
)

REACHABLE_EQUILIBRIUM_BENCH_FIELDS: tuple[str, ...] = (
    *PROTOCOL_MEASUREMENT_FIELDS,
    "baseline_equilibrium_residual",
    "response_column_id",
    "target_reachability_flag",
    "topology_component_label",
    "topology_blocked_flag",
)

SWEEP_SENSITIVITY_METRICS: tuple[str, ...] = (
    "meanAlphaReach",
    "meanZReach",
    "maxObservedNeighborZResidual",
    "maxSuperpositionError",
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
class ResponseAtlasCell:
    cell: tuple[int, int]
    alpha_delta: float
    height_delta: float
    z_residual: float
    alpha_die_off: float | None
    z_die_off: float | None

    def to_dict(self) -> dict[str, object]:
        row, col = self.cell
        return {
            "row": row,
            "col": col,
            "alphaDelta": self.alpha_delta,
            "heightDelta": self.height_delta,
            "zResidual": self.z_residual,
            "alphaDieOff": self.alpha_die_off,
            "zDieOff": self.z_die_off,
        }


@dataclass(frozen=True)
class ResponseAtlasEntry:
    step_id: str
    scope: str
    purpose: str
    expected_response: str
    commands: tuple[SourceCommand, ...]
    locked_cells: tuple[tuple[int, int], ...]
    observation_cells: tuple[ResponseAtlasCell, ...]
    alpha_reach: int
    z_reach: int
    effective_alpha_die_off: int
    effective_z_die_off: int
    max_abs_alpha_delta: float
    max_abs_height_delta: float
    mean_abs_alpha_delta: float
    mean_abs_height_delta: float
    alpha_superposition_error: float
    height_superposition_error: float
    physical_height_rms_error: float | None = None
    physical_center_rms_error: float | None = None
    physical_success: bool | None = None
    physical_energy: float | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "stepId": self.step_id,
            "scope": self.scope,
            "purpose": self.purpose,
            "expectedResponse": self.expected_response,
            "commands": [_source_command_dict(command) for command in self.commands],
            "lockedCells": [_cell_dict(cell) for cell in self.locked_cells],
            "observationCells": [cell.to_dict() for cell in self.observation_cells],
            "alphaReach": self.alpha_reach,
            "zReach": self.z_reach,
            "effectiveAlphaDieOff": self.effective_alpha_die_off,
            "effectiveZDieOff": self.effective_z_die_off,
            "maxAbsAlphaDelta": self.max_abs_alpha_delta,
            "maxAbsHeightDelta": self.max_abs_height_delta,
            "meanAbsAlphaDelta": self.mean_abs_alpha_delta,
            "meanAbsHeightDelta": self.mean_abs_height_delta,
            "alphaSuperpositionError": self.alpha_superposition_error,
            "heightSuperpositionError": self.height_superposition_error,
            "physicalHeightRmsError": self.physical_height_rms_error,
            "physicalCenterRmsError": self.physical_center_rms_error,
            "physicalSuccess": self.physical_success,
            "physicalEnergy": self.physical_energy,
        }


@dataclass(frozen=True)
class ResponseAtlas:
    config_shape: tuple[int, int]
    center_cell: tuple[int, int]
    entries: tuple[ResponseAtlasEntry, ...]
    physical: bool = False
    schema: str = "rad-sim.response-atlas.v1"
    notes: str = (
        "Simulator-generated single, pair, cluster, and lock response atlas for "
        "comparing backlash, clearance, and physical-preview settings."
    )

    @property
    def scope_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for entry in self.entries:
            counts[entry.scope] = counts.get(entry.scope, 0) + 1
        return counts

    def to_dict(self) -> dict[str, object]:
        max_alpha_reach = max((entry.alpha_reach for entry in self.entries), default=0)
        max_z_reach = max((entry.z_reach for entry in self.entries), default=0)
        max_superposition_error = max(
            (
                max(entry.alpha_superposition_error, entry.height_superposition_error)
                for entry in self.entries
            ),
            default=0.0,
        )
        physical_entries = [
            entry for entry in self.entries if entry.physical_success is not None
        ]
        return {
            "schema": self.schema,
            "grid": {"rows": self.config_shape[0], "cols": self.config_shape[1]},
            "centerCell": _cell_dict(self.center_cell),
            "physical": self.physical,
            "notes": self.notes,
            "summary": {
                "entryCount": len(self.entries),
                "scopeCounts": self.scope_counts,
                "maxAlphaReach": max_alpha_reach,
                "maxZReach": max_z_reach,
                "maxSuperpositionError": max_superposition_error,
                "physicalEntryCount": len(physical_entries),
                "physicalSuccessCount": sum(
                    1 for entry in physical_entries if entry.physical_success
                ),
            },
            "assumptions": {
                "source": "Calibration protocol commands run through the current simulator.",
                "physical": (
                    "Optional spring-hinge pass is a model-disagreement diagnostic, not calibrated hardware truth."
                ),
            },
            "entries": [entry.to_dict() for entry in self.entries],
        }


@dataclass(frozen=True)
class ResponseAtlasSweepSample:
    backlash: float
    pin_hole_clearance: float
    pin_radius: float
    hole_radius: float
    atlas: ResponseAtlas

    @property
    def mean_alpha_reach(self) -> float:
        return _mean(entry.alpha_reach for entry in self.atlas.entries)

    @property
    def mean_z_reach(self) -> float:
        return _mean(entry.z_reach for entry in self.atlas.entries)

    @property
    def max_alpha_reach(self) -> int:
        return max((entry.alpha_reach for entry in self.atlas.entries), default=0)

    @property
    def max_z_reach(self) -> int:
        return max((entry.z_reach for entry in self.atlas.entries), default=0)

    @property
    def max_superposition_error(self) -> float:
        return max(
            (
                max(entry.alpha_superposition_error, entry.height_superposition_error)
                for entry in self.atlas.entries
            ),
            default=0.0,
        )

    @property
    def max_observed_z_residual(self) -> float:
        return max(
            (
                abs(cell.z_residual)
                for entry in self.atlas.entries
                for cell in entry.observation_cells
            ),
            default=0.0,
        )

    @property
    def max_observed_neighbor_z_residual(self) -> float:
        residuals: list[float] = []
        for entry in self.atlas.entries:
            command_cells = {command.cell for command in entry.commands}
            residuals.extend(
                abs(cell.z_residual)
                for cell in entry.observation_cells
                if cell.cell not in command_cells
            )
        return max(residuals, default=0.0)

    @property
    def mean_abs_alpha_delta(self) -> float:
        return _mean(entry.mean_abs_alpha_delta for entry in self.atlas.entries)

    @property
    def mean_abs_height_delta(self) -> float:
        return _mean(entry.mean_abs_height_delta for entry in self.atlas.entries)

    @property
    def physical_success_rate(self) -> float | None:
        physical_entries = [
            entry for entry in self.atlas.entries if entry.physical_success is not None
        ]
        if not physical_entries:
            return None
        return sum(1 for entry in physical_entries if entry.physical_success) / len(physical_entries)

    @property
    def mean_physical_height_rms_error(self) -> float | None:
        return _optional_mean(
            entry.physical_height_rms_error for entry in self.atlas.entries
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "settings": {
                "backlash": self.backlash,
                "pinHoleClearance": self.pin_hole_clearance,
                "pinRadius": self.pin_radius,
                "holeRadius": self.hole_radius,
            },
            "summary": {
                "meanAlphaReach": self.mean_alpha_reach,
                "meanZReach": self.mean_z_reach,
                "maxAlphaReach": self.max_alpha_reach,
                "maxZReach": self.max_z_reach,
                "maxSuperpositionError": self.max_superposition_error,
                "maxObservedZResidual": self.max_observed_z_residual,
                "maxObservedNeighborZResidual": self.max_observed_neighbor_z_residual,
                "meanAbsAlphaDelta": self.mean_abs_alpha_delta,
                "meanAbsHeightDelta": self.mean_abs_height_delta,
                "physicalSuccessRate": self.physical_success_rate,
                "meanPhysicalHeightRmsError": self.mean_physical_height_rms_error,
            },
            "atlas": self.atlas.to_dict(),
        }


@dataclass(frozen=True)
class ResponseAtlasSweep:
    config_shape: tuple[int, int]
    center_cell: tuple[int, int]
    backlash_values: tuple[float, ...]
    clearance_values: tuple[float, ...]
    samples: tuple[ResponseAtlasSweepSample, ...]
    physical: bool = False
    schema: str = "rad-sim.response-atlas-sweep.v1"
    notes: str = (
        "Simulator-generated sweep for comparing how backlash and pin-hole "
        "clearance change locality, residual vertical motion, and operator "
        "interaction before calibrated bench data exists."
    )

    def to_dict(self) -> dict[str, object]:
        by_backlash = _sweep_trend(self.samples, "backlash")
        by_clearance = _sweep_trend(self.samples, "pin_hole_clearance")
        trends = {
            "byBacklash": by_backlash,
            "byPinHoleClearance": by_clearance,
        }
        sensitivity = _sweep_sensitivity_from_trends(trends)
        return {
            "schema": self.schema,
            "grid": {"rows": self.config_shape[0], "cols": self.config_shape[1]},
            "centerCell": _cell_dict(self.center_cell),
            "physical": self.physical,
            "parameters": {
                "backlashValues": list(self.backlash_values),
                "pinHoleClearanceValues": list(self.clearance_values),
            },
            "notes": self.notes,
            "summary": {
                "sampleCount": len(self.samples),
                "maxAlphaReach": max(
                    (sample.max_alpha_reach for sample in self.samples),
                    default=0,
                ),
                "maxZReach": max((sample.max_z_reach for sample in self.samples), default=0),
                "maxObservedZResidual": max(
                    (sample.max_observed_z_residual for sample in self.samples),
                    default=0.0,
                ),
                "maxObservedNeighborZResidual": max(
                    (sample.max_observed_neighbor_z_residual for sample in self.samples),
                    default=0.0,
                ),
                "maxSuperpositionError": max(
                    (sample.max_superposition_error for sample in self.samples),
                    default=0.0,
                ),
            },
            "trends": trends,
            "sensitivity": sensitivity,
            "operatorLawCandidates": _sweep_operator_law_candidates(trends, sensitivity),
            "assumptions": {
                "source": "Each sample rebuilds the response atlas from the same protocol commands.",
                "interpretation": (
                    "Backlash and clearance are treated as programmable-discontinuity "
                    "dead-zone parameters; trend summaries are simulator diagnostics."
                ),
                "physical": (
                    "Optional spring-hinge metrics are model-disagreement diagnostics, "
                    "not calibrated hardware validation."
                ),
            },
            "samples": [sample.to_dict() for sample in self.samples],
        }


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
    dataset_id: str | None = None
    dataset_role: str = "unassigned"
    source_file_id: str | None = None
    collected_at: str | None = None
    operator: str | None = None
    profile_id: str | None = None
    profile_frozen_at: str | None = None
    provenance_notes: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "protocolSchema": self.protocol_schema,
            "hardwareProfile": self.hardware_profile_name,
            "provenance": {
                "schema": "rad-sim.calibration-dataset-provenance.v1",
                "datasetId": self.dataset_id,
                "datasetRole": self.dataset_role,
                "sourceFileId": self.source_file_id,
                "collectedAt": self.collected_at,
                "operator": self.operator,
                "profileId": self.profile_id,
                "profileFrozenAt": self.profile_frozen_at,
                "notes": self.provenance_notes,
            },
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
    removed_cells: Iterable[tuple[int, int]] = (),
) -> LatticeState:
    state = LatticeState.uniform(config)
    removed = {(int(r), int(c)) for r, c in removed_cells}
    for r, c in removed:
        state.removed_mask[r, c] = True
    for r, c in locked_cells:
        if (int(r), int(c)) in removed:
            continue
        state.locked_mask[r, c] = True
    for command in commands:
        r, c = command.cell
        if (int(r), int(c)) in removed:
            continue
        state.actuator_grid[r, c] += command.alpha
        state.z_actuator_grid[r, c] += command.z
    return state


def _cell_dict(cell: tuple[int, int]) -> dict[str, int]:
    return {"row": int(cell[0]), "col": int(cell[1])}


def _json_ready_topology(topology: dict[str, object]) -> dict[str, object]:
    ready: dict[str, object] = {}
    for key, value in topology.items():
        if isinstance(value, np.ndarray):
            ready[key] = value.tolist()
        elif isinstance(value, np.integer):
            ready[key] = int(value)
        elif isinstance(value, np.floating):
            ready[key] = float(value)
        else:
            ready[key] = value
    return ready


def _array_rms(value: np.ndarray | None) -> float | None:
    if value is None:
        return None
    if value.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(value**2)))


def _array_max_abs(value: np.ndarray | None) -> float | None:
    if value is None:
        return None
    if value.size == 0:
        return 0.0
    return float(np.max(np.abs(value)))


def _target_matrix(
    value: np.ndarray | Iterable[Iterable[float]] | None,
    shape: tuple[int, int],
    name: str,
) -> np.ndarray | None:
    if value is None:
        return None
    matrix = np.asarray(value, dtype=float)
    if matrix.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    return matrix


def _component_reach_counts(
    topology: dict[str, object],
    actuator_cells: Iterable[tuple[int, int]],
) -> tuple[int, int]:
    labels = np.asarray(topology["component_labels"], dtype=int)
    active_components = {
        int(labels[r, c])
        for r, c in actuator_cells
        if 0 <= r < labels.shape[0] and 0 <= c < labels.shape[1] and int(labels[r, c]) >= 0
    }
    present_mask = labels >= 0
    reachable = present_mask & np.isin(labels, list(active_components))
    return int(np.count_nonzero(reachable)), int(
        np.count_nonzero(present_mask & ~reachable)
    )


def _matrix_rank_or_zero(matrix: np.ndarray) -> int:
    if matrix.size == 0 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        return 0
    return int(np.linalg.matrix_rank(matrix))


def _reachable_rows(matrix: np.ndarray, tolerance: float) -> int:
    if matrix.size == 0 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        return 0
    return int(np.count_nonzero(np.any(np.abs(matrix) > tolerance, axis=1)))


def _component_response_summaries(
    topology: dict[str, object],
    matrix: ResponseMatrix,
    tolerance: float,
) -> tuple[TopologyComponentResponseSummary, ...]:
    labels = np.asarray(topology["component_labels"], dtype=int)
    summaries: list[TopologyComponentResponseSummary] = []
    component_count = int(topology["component_count"])
    command_labels = []
    for command in matrix.commands:
        r, c = command.cell
        command_labels.append(int(labels[r, c]) if 0 <= r < labels.shape[0] and 0 <= c < labels.shape[1] else -1)
    command_labels_array = np.asarray(command_labels, dtype=int)
    for label in range(component_count):
        row_indices = np.flatnonzero(labels.reshape(-1) == label)
        column_indices = np.flatnonzero(command_labels_array == label)
        alpha_component = matrix.alpha[np.ix_(row_indices, column_indices)]
        height_component = matrix.height[np.ix_(row_indices, column_indices)]
        actuator_cells = {
            matrix.commands[index].cell
            for index in column_indices
        }
        summaries.append(
            TopologyComponentResponseSummary(
                label=label,
                cell_count=int(row_indices.size),
                actuator_cell_count=len(actuator_cells),
                command_column_count=int(column_indices.size),
                alpha_rank=_matrix_rank_or_zero(alpha_component),
                height_rank=_matrix_rank_or_zero(height_component),
                reachable_alpha_cells=_reachable_rows(alpha_component, tolerance),
                reachable_height_cells=_reachable_rows(height_component, tolerance),
            )
        )
    return tuple(summaries)


def _scenario_comparisons(
    scenarios: tuple[TopologyExperimentScenarioResult, ...]
) -> list[dict[str, object]]:
    if not scenarios:
        return []
    baseline = scenarios[0]
    comparisons: list[dict[str, object]] = []
    for scenario in scenarios[1:]:
        comparisons.append(
            {
                "baseline": baseline.scenario.name,
                "scenario": scenario.scenario.name,
                "alphaRankDelta": scenario.response_matrix.alpha_rank - baseline.response_matrix.alpha_rank,
                "heightRankDelta": scenario.response_matrix.height_rank - baseline.response_matrix.height_rank,
                "reachableAlphaCellDelta": scenario.response_matrix.reachable_alpha_cells(scenario.tolerance)
                - baseline.response_matrix.reachable_alpha_cells(baseline.tolerance),
                "reachableHeightCellDelta": scenario.response_matrix.reachable_height_cells(scenario.tolerance)
                - baseline.response_matrix.reachable_height_cells(baseline.tolerance),
                "componentCountDelta": int(scenario.topology["component_count"])
                - int(baseline.topology["component_count"]),
                "componentBlockedCellDelta": scenario.component_blocked_cells
                - baseline.component_blocked_cells,
                "deletedEdgeDelta": int(scenario.topology["deleted_edge_count"])
                - int(baseline.topology["deleted_edge_count"]),
            }
        )
    return comparisons


def _signed_height_reach_counts(
    config: LatticeConfig,
    actuator_cells: Iterable[tuple[int, int]],
    locked_cells: Iterable[tuple[int, int]],
    removed_cells: Iterable[tuple[int, int]],
    z_step: float,
    tolerance: float,
) -> tuple[int, int]:
    positive = np.zeros((config.rows, config.cols), dtype=bool)
    negative = np.zeros((config.rows, config.cols), dtype=bool)
    step = abs(float(z_step))
    if step <= tolerance:
        return 0, 0
    locked = tuple(locked_cells)
    removed = tuple(removed_cells)
    for cell in actuator_cells:
        for signed_step in (step, -step):
            response = characterize_response(
                config,
                (SourceCommand(cell=cell, z=signed_step),),
                locked,
                tolerance,
                removed_cells=removed,
            )
            positive |= response.height_delta > tolerance
            negative |= response.height_delta < -tolerance
    return int(np.count_nonzero(positive)), int(np.count_nonzero(negative))


def _source_command_dict(command: SourceCommand) -> dict[str, object]:
    return {
        "row": int(command.cell[0]),
        "col": int(command.cell[1]),
        "alpha": float(command.alpha),
        "z": float(command.z),
    }


def _finite_or_none(value: object) -> float | None:
    try:
        numeric = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return numeric if np.isfinite(numeric) else None


def _mean(values: Iterable[float | int]) -> float:
    numeric = [float(value) for value in values]
    return float(np.mean(numeric)) if numeric else 0.0


def _optional_mean(values: Iterable[float | None]) -> float | None:
    numeric = [float(value) for value in values if value is not None]
    return float(np.mean(numeric)) if numeric else None


def _sweep_trend(
    samples: Iterable[ResponseAtlasSweepSample],
    parameter: Literal["backlash", "pin_hole_clearance"],
) -> list[dict[str, object]]:
    buckets: dict[float, list[ResponseAtlasSweepSample]] = {}
    for sample in samples:
        value = float(getattr(sample, parameter))
        buckets.setdefault(value, []).append(sample)

    trend: list[dict[str, object]] = []
    for value, bucket in sorted(buckets.items()):
        trend.append(
            {
                "value": value,
                "sampleCount": len(bucket),
                "meanAlphaReach": _mean(sample.mean_alpha_reach for sample in bucket),
                "meanZReach": _mean(sample.mean_z_reach for sample in bucket),
                "maxObservedZResidual": max(
                    (sample.max_observed_z_residual for sample in bucket),
                    default=0.0,
                ),
                "maxObservedNeighborZResidual": max(
                    (sample.max_observed_neighbor_z_residual for sample in bucket),
                    default=0.0,
                ),
                "maxSuperpositionError": max(
                    (sample.max_superposition_error for sample in bucket),
                    default=0.0,
                ),
            }
        )
    return trend


def _trend_endpoint_slope(trend: list[dict[str, object]], metric: str) -> float | None:
    if len(trend) < 2:
        return None
    first = trend[0]
    last = trend[-1]
    dx = float(last["value"]) - float(first["value"])
    if abs(dx) <= 1e-12:
        return None
    return (float(last[metric]) - float(first[metric])) / dx


def _sweep_sensitivity_from_trends(
    trends: dict[str, list[dict[str, object]]],
) -> dict[str, object]:
    metrics: dict[str, dict[str, float | None]] = {}
    dominant: dict[str, object] | None = None
    parameter_labels = {
        "byBacklash": "backlash",
        "byPinHoleClearance": "pinHoleClearance",
    }
    for trend_key, parameter in parameter_labels.items():
        trend = trends.get(trend_key, [])
        parameter_metrics: dict[str, float | None] = {}
        for metric in SWEEP_SENSITIVITY_METRICS:
            slope = _trend_endpoint_slope(trend, metric)
            parameter_metrics[metric] = slope
            if slope is None:
                continue
            candidate = {
                "parameter": parameter,
                "metric": metric,
                "slope": slope,
                "absSlope": abs(slope),
            }
            if dominant is None or float(candidate["absSlope"]) > float(dominant["absSlope"]):
                dominant = candidate
        metrics[parameter] = parameter_metrics
    return {
        "method": "endpoint finite difference over each parameter trend",
        "metrics": metrics,
        "dominant": dominant,
    }


def _trend_monotonicity(
    trend: list[dict[str, object]],
    metric: str,
    tolerance: float = 1e-12,
) -> str:
    if len(trend) < 2:
        return "insufficient"
    increases = 0
    decreases = 0
    flats = 0
    for previous, current in zip(trend, trend[1:]):
        delta = float(current[metric]) - float(previous[metric])
        if abs(delta) <= tolerance:
            flats += 1
        elif delta > 0:
            increases += 1
        else:
            decreases += 1
    if increases and not decreases:
        return "increasing"
    if decreases and not increases:
        return "decreasing"
    if flats and not increases and not decreases:
        return "flat"
    return "mixed"


def _operator_law_statement(parameter: str, metric: str, monotonicity: str) -> str:
    parameter_labels = {
        "backlash": "backlash dead-zone width",
        "pinHoleClearance": "pin-hole clearance",
    }
    metric_labels = {
        "meanAlphaReach": "mean dilation reach",
        "meanZReach": "mean vertical reach",
        "maxObservedNeighborZResidual": "neighbor vertical residual motion",
        "maxSuperpositionError": "non-additive superposition residual",
    }
    parameter_label = parameter_labels.get(parameter, parameter)
    metric_label = metric_labels.get(metric, metric)
    if monotonicity == "increasing":
        return f"Increasing {parameter_label} increases {metric_label} over the sampled simulator sweep."
    if monotonicity == "decreasing":
        return f"Increasing {parameter_label} decreases {metric_label} over the sampled simulator sweep."
    if monotonicity == "flat":
        return f"Changing {parameter_label} leaves {metric_label} approximately flat over the sampled simulator sweep."
    if monotonicity == "mixed":
        return f"{parameter_label} has a mixed sampled relationship with {metric_label}; no monotone candidate is supported."
    return f"{parameter_label} has insufficient sampled data to propose a {metric_label} law candidate."


def _sweep_operator_law_candidates(
    trends: dict[str, list[dict[str, object]]],
    sensitivity: dict[str, object],
) -> dict[str, object]:
    trend_map = {
        "backlash": trends.get("byBacklash", []),
        "pinHoleClearance": trends.get("byPinHoleClearance", []),
    }
    sensitivity_metrics = sensitivity.get("metrics", {})
    laws: list[dict[str, object]] = []
    for parameter, trend in trend_map.items():
        parameter_sensitivity = {}
        if isinstance(sensitivity_metrics, dict):
            parameter_sensitivity = sensitivity_metrics.get(parameter, {}) or {}
        for metric in SWEEP_SENSITIVITY_METRICS:
            monotonicity = _trend_monotonicity(trend, metric)
            slope = (
                parameter_sensitivity.get(metric)
                if isinstance(parameter_sensitivity, dict)
                else None
            )
            laws.append(
                {
                    "parameter": parameter,
                    "metric": metric,
                    "monotonicity": monotonicity,
                    "slope": slope,
                    "supportedBySweep": monotonicity in {"increasing", "decreasing", "flat"},
                    "status": "simulator-diagnostic",
                    "statement": _operator_law_statement(parameter, metric, monotonicity),
                }
            )
    return {
        "schema": "rad-sim.operator-law-candidates.v1",
        "method": "adjacent monotonicity over sampled parameter trend plus endpoint sensitivity",
        "laws": laws,
    }


def _atlas_observation_cells(
    response: ResponseCharacterization,
    cells: Iterable[tuple[int, int]],
) -> tuple[ResponseAtlasCell, ...]:
    observed: list[ResponseAtlasCell] = []
    for row, col in cells:
        observed.append(
            ResponseAtlasCell(
                cell=(int(row), int(col)),
                alpha_delta=float(response.alpha_delta[row, col]),
                height_delta=float(response.height_delta[row, col]),
                z_residual=float(response.z_residual[row, col]),
                alpha_die_off=_finite_or_none(float(response.alpha_die_off[row, col])),
                z_die_off=_finite_or_none(float(response.z_die_off[row, col])),
            )
        )
    return tuple(observed)


def _atlas_superposition_errors(
    config: LatticeConfig,
    commands: tuple[SourceCommand, ...],
    combined: ResponseCharacterization,
    locked_cells: tuple[tuple[int, int], ...],
    tolerance: float,
) -> tuple[float, float]:
    if len(commands) <= 1:
        return 0.0, 0.0
    alpha_expected = np.zeros_like(combined.alpha_delta)
    height_expected = np.zeros_like(combined.height_delta)
    for command in commands:
        single = characterize_response(
            config,
            (command,),
            locked_cells=locked_cells,
            tolerance=tolerance,
        )
        alpha_expected += single.alpha_delta
        height_expected += single.height_delta
    return (
        float(np.max(np.abs(combined.alpha_delta - alpha_expected))),
        float(np.max(np.abs(combined.height_delta - height_expected))),
    )


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
    removed_cells: Iterable[tuple[int, int]] = (),
) -> ResponseCharacterization:
    commands = tuple(commands)
    locked_cells = tuple(locked_cells)
    removed_cells = tuple(removed_cells)
    baseline = simulate_kinematic(
        config, _state_with_commands(config, (), locked_cells, removed_cells)
    )
    result = simulate_kinematic(
        config, _state_with_commands(config, commands, locked_cells, removed_cells)
    )
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
        positive_z_reach=int(np.count_nonzero(height_delta > tolerance)),
        negative_z_reach=int(np.count_nonzero(height_delta < -tolerance)),
        max_positive_height_delta=float(max(0.0, np.max(height_delta))),
        max_negative_height_delta=float(min(0.0, np.min(height_delta))),
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


def calibration_experiment_protocol_from_dict(
    raw: dict[str, object],
) -> CalibrationExperimentProtocol:
    schema = str(raw.get("schema", ""))
    if schema != "rad-sim.calibration-experiment-protocol.v1":
        raise ValueError("unsupported calibration experiment protocol schema")
    grid = raw.get("grid", {})
    if not isinstance(grid, dict):
        raise ValueError("calibration experiment protocol requires a grid object")
    center_raw = raw.get("centerCell", {})
    if not isinstance(center_raw, dict):
        raise ValueError("calibration experiment protocol requires a centerCell object")
    steps_raw = raw.get("steps", [])
    if not isinstance(steps_raw, list):
        raise ValueError("calibration experiment protocol requires a steps list")
    steps: list[CalibrationExperimentStep] = []
    for step_raw in steps_raw:
        if not isinstance(step_raw, dict):
            raise ValueError("each calibration experiment step must be an object")
        commands_raw = step_raw.get("commands", [])
        observation_raw = step_raw.get("observationCells", [])
        locked_raw = step_raw.get("lockedCells", [])
        if not isinstance(commands_raw, list):
            raise ValueError("each calibration experiment step requires commands")
        if not isinstance(observation_raw, list):
            raise ValueError("each calibration experiment step requires observationCells")
        if not isinstance(locked_raw, list):
            raise ValueError("each calibration experiment step requires lockedCells")

        def _cells(values: list[object]) -> tuple[tuple[int, int], ...]:
            cells: list[tuple[int, int]] = []
            for cell_raw in values:
                if not isinstance(cell_raw, dict):
                    continue
                cells.append((int(cell_raw["row"]), int(cell_raw["col"])))
            return tuple(cells)

        commands: list[SourceCommand] = []
        for command_raw in commands_raw:
            if not isinstance(command_raw, dict):
                continue
            commands.append(
                SourceCommand(
                    cell=(int(command_raw["row"]), int(command_raw["col"])),
                    alpha=float(command_raw.get("alpha", 0.0)),
                    z=float(command_raw.get("z", 0.0)),
                )
            )
        fields_raw = step_raw.get("measurementFields", PROTOCOL_MEASUREMENT_FIELDS)
        fields = (
            tuple(str(field) for field in fields_raw)
            if isinstance(fields_raw, list)
            else PROTOCOL_MEASUREMENT_FIELDS
        )
        steps.append(
            CalibrationExperimentStep(
                id=str(step_raw["id"]),
                scope=str(step_raw.get("scope", "single")),  # type: ignore[arg-type]
                commands=tuple(commands),
                observation_cells=_cells(observation_raw),
                locked_cells=_cells(locked_raw),
                measurement_fields=fields,
                purpose=str(step_raw.get("purpose", "")),
                expected_response=str(step_raw.get("expectedResponse", "")),
                repeat_count=int(step_raw.get("repeatCount", 1)),
            )
        )
    return CalibrationExperimentProtocol(
        config_shape=(int(grid["rows"]), int(grid["cols"])),
        center_cell=(int(center_raw["row"]), int(center_raw["col"])),
        steps=tuple(steps),
        hardware_profile_name=str(raw.get("hardwareProfile", "paper-reference")),
        notes=str(raw.get("notes", "")),
    )


def calibration_experiment_protocol_from_json(
    payload: str,
) -> CalibrationExperimentProtocol:
    raw = json.loads(payload)
    if not isinstance(raw, dict):
        raise ValueError("calibration experiment protocol JSON must be an object")
    return calibration_experiment_protocol_from_dict(raw)


def calibration_experiment_results_template(
    protocol: CalibrationExperimentProtocol,
    *,
    dataset_id: str | None = None,
    dataset_role: str = "unassigned",
    source_file_id: str | None = None,
    collected_at: str | None = None,
    operator: str | None = None,
    profile_id: str | None = None,
    profile_frozen_at: str | None = None,
    provenance_notes: str = "",
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
        dataset_id=dataset_id,
        dataset_role=dataset_role,
        source_file_id=source_file_id,
        collected_at=collected_at,
        operator=operator,
        profile_id=profile_id,
        profile_frozen_at=profile_frozen_at,
        provenance_notes=provenance_notes,
    )


def export_calibration_experiment_results_template_json(
    protocol: CalibrationExperimentProtocol,
    *,
    dataset_id: str | None = None,
    dataset_role: str = "unassigned",
    source_file_id: str | None = None,
    collected_at: str | None = None,
    operator: str | None = None,
    profile_id: str | None = None,
    profile_frozen_at: str | None = None,
    provenance_notes: str = "",
) -> str:
    return json.dumps(
        calibration_experiment_results_template(
            protocol,
            dataset_id=dataset_id,
            dataset_role=dataset_role,
            source_file_id=source_file_id,
            collected_at=collected_at,
            operator=operator,
            profile_id=profile_id,
            profile_frozen_at=profile_frozen_at,
            provenance_notes=provenance_notes,
        ).to_dict(),
        indent=2,
    )


def calibration_bench_notebook(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    *,
    fit_dataset_id: str = "fit-run-001",
    holdout_dataset_id: str = "holdout-run-001",
    profile_id: str = "calibration-profile-v1",
    profile_frozen_at: str | None = None,
) -> dict[str, object]:
    """Build a bench-facing protocol notebook for calibration data collection."""

    frozen_at = profile_frozen_at or "set-after-fit-before-holdout"
    fit_filename = f"{fit_dataset_id}.json"
    holdout_filename = f"{holdout_dataset_id}.json"
    datasets = {
        "fit": {
            "datasetId": fit_dataset_id,
            "role": "fit",
            "suggestedFilename": fit_filename,
            "sourceFileId": fit_filename,
            "profileId": None,
            "profileFrozenAt": None,
            "provenanceSchema": "rad-sim.calibration-dataset-provenance.v1",
            "purpose": "Estimate diagnostic response parameters and freeze a bounded model profile.",
        },
        "holdout": {
            "datasetId": holdout_dataset_id,
            "role": "holdout",
            "suggestedFilename": holdout_filename,
            "sourceFileId": holdout_filename,
            "profileId": profile_id,
            "profileFrozenAt": frozen_at,
            "provenanceSchema": "rad-sim.calibration-dataset-provenance.v1",
            "purpose": "Replay the frozen profile against data collected after the fit file is closed.",
        },
    }
    scenarios: list[dict[str, object]] = []
    for step in protocol.steps:
        command_cells = [_cell_dict(command.cell) for command in step.commands]
        scenarios.append(
            {
                "stepId": step.id,
                "scope": step.scope,
                "repeatCount": step.repeat_count,
                "commandCells": command_cells,
                "commands": [_source_command_dict(command) for command in step.commands],
                "lockedCells": [_cell_dict(cell) for cell in step.locked_cells],
                "observationCells": [_cell_dict(cell) for cell in step.observation_cells],
                "measurementFields": list(step.measurement_fields),
                "datasetRoles": ["fit", "holdout"],
                "purpose": step.purpose,
                "expectedResponse": step.expected_response,
                "claimLabel": "bench protocol step; no physical law claimed until measured",
            }
        )
    return {
        "schema": "rad-sim.calibration-bench-notebook.v1",
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "cellSize": config.cell_size,
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "pinHoleClearance": config.pin_hole_clearance,
        },
        "protocol": {
            "schema": protocol.schema,
            "hardwareProfile": protocol.hardware_profile_name,
            "centerCell": _cell_dict(protocol.center_cell),
            "stepCount": len(protocol.steps),
            "repeatCountTotal": sum(step.repeat_count for step in protocol.steps),
            "notes": protocol.notes,
        },
        "datasetPlan": datasets,
        "instruments": [
            {
                "id": "calipers-or-micrometer",
                "measures": ["pin radius", "hole radius", "plate thickness", "joint stack height"],
                "claimLabel": "required geometry measurement",
            },
            {
                "id": "motion-tracking-camera",
                "measures": ["cell center displacement", "height delta", "neighbor residual z"],
                "claimLabel": "required response measurement",
            },
            {
                "id": "force-gauge-or-load-cell",
                "measures": ["actuator force", "contact/load proxy"],
                "claimLabel": "required force measurement for contact calibration",
            },
            {
                "id": "fixture-and-clamps",
                "measures": ["locked cell enforcement", "fixed boundary repeatability"],
                "claimLabel": "required boundary-condition control",
            },
            {
                "id": "actuator-controller",
                "measures": ["commanded alpha", "commanded z", "command repeat timing"],
                "claimLabel": "required command provenance",
            },
            {
                "id": "scale-reference-marker",
                "measures": ["pixel-to-mm scale", "coordinate registration"],
                "claimLabel": "required unit-scale evidence",
            },
        ],
        "measurementColumns": [
            "datasetId",
            "datasetRole",
            "sourceFileId",
            "collectedAt",
            "operator",
            "profileId",
            "profileFrozenAt",
            *PROTOCOL_MEASUREMENT_FIELDS,
        ],
        "phases": [
            {
                "id": "setup",
                "description": "Measure hardware profile fields and register the cell coordinate frame.",
                "requiredEvidence": ["hardware profile", "scale marker", "fixture notes"],
            },
            {
                "id": "fit-collection",
                "description": "Collect the fit dataset without using any holdout measurements.",
                "requiredEvidence": [fit_dataset_id, "datasetRole=fit"],
            },
            {
                "id": "profile-freeze",
                "description": "Generate and freeze a model profile before holdout collection starts.",
                "requiredEvidence": [profile_id, frozen_at],
            },
            {
                "id": "holdout-collection",
                "description": "Collect a separate holdout file after the model profile is frozen.",
                "requiredEvidence": [holdout_dataset_id, "datasetRole=holdout"],
            },
            {
                "id": "validation",
                "description": "Run residual validation and independent split checks before making physical claims.",
                "requiredEvidence": [
                    "rad-sim.calibration-model-profile-holdout-validation.v1",
                    "rad-sim.calibration-train-holdout-split.v1",
                ],
            },
        ],
        "scenarios": scenarios,
        "outputs": [
            {
                "id": "protocol-json",
                "schema": protocol.schema,
                "function": "export_calibration_experiment_protocol_json",
                "browserAction": "Save Protocol",
            },
            {
                "id": "fit-results-template-json",
                "schema": "rad-sim.calibration-experiment-results.v1",
                "function": "export_calibration_experiment_results_template_json",
                "datasetRole": "fit",
                "browserAction": "Save Results Template",
            },
            {
                "id": "holdout-results-template-json",
                "schema": "rad-sim.calibration-experiment-results.v1",
                "function": "export_calibration_experiment_results_template_json",
                "datasetRole": "holdout",
                "browserAction": "Save Results Template",
            },
            {
                "id": "model-profile-json",
                "schema": "rad-sim.calibration-model-profile.v1",
                "function": "export_calibration_model_profile_json",
                "browserAction": "Save Model Profile",
            },
            {
                "id": "holdout-validation-json",
                "schema": "rad-sim.calibration-model-profile-holdout-validation.v1",
                "function": "export_calibration_model_profile_holdout_validation_json",
                "browserAction": "Save Holdout Check",
            },
            {
                "id": "holdout-validation-csv",
                "schema": "rad-sim.calibration-model-profile-holdout-validation.v1",
                "function": "export_calibration_model_profile_holdout_validation_csv",
                "browserAction": "Save Holdout CSV",
            },
            {
                "id": "bench-notebook-csv",
                "schema": "rad-sim.calibration-bench-notebook.v1",
                "function": "export_calibration_bench_notebook_csv",
                "browserAction": "Save Bench CSV",
            },
        ],
        "passFailCriteria": {
            "residualValidationPass": "fit and holdout residual scores do not increase after an applied profile update",
            "independentValidationPass": "residual pass plus complete train/holdout provenance metadata",
            "missingEvidence": "must be empty before claiming independent hardware validation",
            "physicalClaimLimit": "passing this protocol supports calibration bookkeeping only; contact, friction, stiffness, and material laws still require model-specific validation",
        },
        "formalization": {
            "targetId": "calibration_bench_protocol_coverage",
            "leanStructure": "Mechanics.CalibrationBenchProtocolCoverageNat",
            "leanPredicate": "calibrationBenchProtocolCoverageReadyNat",
            "schema": "rad-sim.calibration-bench-notebook.v1",
        },
        "claimLabels": {
            "notebook": "bench protocol artifact",
            "scenarioRows": "experimentally unvalidated physical procedure until performed",
            "formalCoverage": "Lean-proven finite coverage predicate, not physical accuracy",
        },
        "limitations": [
            "The notebook enforces artifact coverage and provenance fields but cannot prove the lab actually collected independent files.",
            "The protocol is normalized and must be mapped to measured hardware units before physical claims.",
            "Rigid-body contact, friction, stiffness, gravity sag, and actuator-force laws remain uncalibrated.",
        ],
    }


def export_calibration_bench_notebook_json(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    **kwargs: object,
) -> str:
    return json.dumps(
        calibration_bench_notebook(config, protocol, **kwargs),
        indent=2,
    )


def export_calibration_bench_notebook_csv(notebook: dict[str, object]) -> str:
    """Export one row per calibration protocol step for bench execution."""

    header = [
        "schema",
        "step_id",
        "scope",
        "repeat_count",
        "command_cells",
        "locked_cells",
        "observation_cells",
        "measurement_columns",
        "fit_dataset_id",
        "holdout_dataset_id",
        "purpose",
        "expected_response",
        "claim_label",
    ]
    dataset_plan = notebook.get("datasetPlan", {})
    fit = dataset_plan.get("fit", {}) if isinstance(dataset_plan, dict) else {}
    holdout = dataset_plan.get("holdout", {}) if isinstance(dataset_plan, dict) else {}
    rows: list[list[object]] = [header]
    for scenario in notebook.get("scenarios", []):
        if not isinstance(scenario, dict):
            continue
        rows.append(
            [
                notebook.get("schema", ""),
                scenario.get("stepId", ""),
                scenario.get("scope", ""),
                scenario.get("repeatCount", ""),
                json.dumps(scenario.get("commandCells", []), separators=(",", ":")),
                json.dumps(scenario.get("lockedCells", []), separators=(",", ":")),
                json.dumps(scenario.get("observationCells", []), separators=(",", ":")),
                ";".join(str(field) for field in scenario.get("measurementFields", []))
                if isinstance(scenario.get("measurementFields"), list)
                else "",
                fit.get("datasetId") if isinstance(fit, dict) else "",
                holdout.get("datasetId") if isinstance(holdout, dict) else "",
                scenario.get("purpose", ""),
                scenario.get("expectedResponse", ""),
                scenario.get("claimLabel", ""),
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row) for row in rows
    )


def calibration_bench_packet(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    *,
    fit_dataset_id: str = "fit-run-001",
    holdout_dataset_id: str = "holdout-run-001",
    profile_id: str = "calibration-profile-v1",
    profile_frozen_at: str | None = None,
) -> dict[str, object]:
    """Bundle the calibration protocol into a complete bench handoff packet."""

    frozen_at = profile_frozen_at or "set-after-fit-before-holdout"
    notebook = calibration_bench_notebook(
        config,
        protocol,
        fit_dataset_id=fit_dataset_id,
        holdout_dataset_id=holdout_dataset_id,
        profile_id=profile_id,
        profile_frozen_at=frozen_at,
    )
    fit_template = calibration_experiment_results_template(
        protocol,
        dataset_id=fit_dataset_id,
        dataset_role="fit",
        source_file_id=f"{fit_dataset_id}.json",
        profile_id=None,
        profile_frozen_at=None,
        provenance_notes="Fit file: collect before generating and freezing a profile.",
    ).to_dict()
    holdout_template = calibration_experiment_results_template(
        protocol,
        dataset_id=holdout_dataset_id,
        dataset_role="holdout",
        source_file_id=f"{holdout_dataset_id}.json",
        profile_id=profile_id,
        profile_frozen_at=frozen_at,
        provenance_notes=(
            "Holdout file: collect only after the fit-derived profile is frozen."
        ),
    ).to_dict()
    filenames = {
        "packet": "calibration_bench_packet.json",
        "notebook": "calibration_bench_notebook.json",
        "notebook_csv": "calibration_bench_notebook.csv",
        "protocol": "calibration_experiment_protocol.json",
        "fit_template": "calibration_fit_results_template.json",
        "holdout_template": "calibration_holdout_results_template.json",
        "readme": "README.md",
    }
    return {
        "schema": "rad-sim.calibration-bench-packet.v1",
        "method": (
            "bench handoff packet bundling calibration protocol, fit and "
            "holdout measurement templates, notebook scenario table, and "
            "validation instructions"
        ),
        "schemas": {
            "notebook": "rad-sim.calibration-bench-notebook.v1",
            "protocol": protocol.schema,
            "resultsTemplate": "rad-sim.calibration-experiment-results.v1",
            "datasetProvenance": "rad-sim.calibration-dataset-provenance.v1",
            "modelProfile": "rad-sim.calibration-model-profile.v1",
            "holdoutValidation": "rad-sim.calibration-model-profile-holdout-validation.v1",
            "splitMetadata": "rad-sim.calibration-train-holdout-split.v1",
        },
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "cellSize": config.cell_size,
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "pinHoleClearance": config.pin_hole_clearance,
        },
        "datasetPlan": notebook["datasetPlan"],
        "filenames": filenames,
        "benchNotebook": notebook,
        "benchNotebookCsv": export_calibration_bench_notebook_csv(notebook),
        "experimentProtocol": protocol.to_dict(),
        "fitResultsTemplate": fit_template,
        "holdoutResultsTemplate": holdout_template,
        "artifactManifest": [
            {
                "id": "packet",
                "filename": filenames["packet"],
                "schema": "rad-sim.calibration-bench-packet.v1",
                "purpose": "single JSON bundle for review and archival",
            },
            {
                "id": "notebook",
                "filename": filenames["notebook"],
                "schema": "rad-sim.calibration-bench-notebook.v1",
                "purpose": "human-readable protocol plan with phases and outputs",
            },
            {
                "id": "notebook-csv",
                "filename": filenames["notebook_csv"],
                "schema": "rad-sim.calibration-bench-notebook.v1",
                "purpose": "one-row-per-scenario table for lab notebook or spreadsheet",
            },
            {
                "id": "protocol",
                "filename": filenames["protocol"],
                "schema": protocol.schema,
                "purpose": "machine-readable protocol used by simulator comparisons",
            },
            {
                "id": "fit-template",
                "filename": filenames["fit_template"],
                "schema": "rad-sim.calibration-experiment-results.v1",
                "purpose": "blank fit measurement rows",
            },
            {
                "id": "holdout-template",
                "filename": filenames["holdout_template"],
                "schema": "rad-sim.calibration-experiment-results.v1",
                "purpose": "blank holdout measurement rows tied to frozen profile ID",
            },
        ],
        "validationInstructions": {
            "fitStep": "Fill the fit template first and use it to generate a model profile.",
            "freezeStep": "Assign profileId and profileFrozenAt before collecting holdout measurements.",
            "holdoutStep": "Fill the holdout template in a separate raw file after freezing the profile.",
            "compareFunction": "calibration_model_profile_holdout_validation",
            "csvFunction": "export_calibration_model_profile_holdout_validation_csv",
            "acceptanceRule": (
                "Treat independentValidationPass as meaningful only when residual "
                "validation passes and splitMetadata.missingEvidence is empty."
            ),
        },
        "formalization": {
            "targetId": "calibration_bench_packet_completeness",
            "leanStructure": "Mechanics.CalibrationBenchPacketCompletenessNat",
            "leanPredicate": "calibrationBenchPacketCompleteNat",
            "schema": "rad-sim.calibration-bench-packet.v1",
        },
        "claimLabels": {
            "packetAssembly": "bench protocol artifact",
            "fitTemplate": "experimentally unvalidated physical procedure until filled",
            "holdoutTemplate": "experimentally unvalidated physical procedure until independently filled",
            "formalCompleteness": "Lean-proven finite packet-completeness predicate, not physical accuracy",
        },
        "limitations": [
            "The packet separates fit and holdout files but cannot prove laboratory independence by itself.",
            "The templates are blank until real bench measurements replace null values.",
            "Passing simulator residual checks is not a proof of contact, friction, stiffness, gravity, or material physics.",
        ],
    }


def export_calibration_bench_packet_json(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    **kwargs: object,
) -> str:
    return json.dumps(
        calibration_bench_packet(config, protocol, **kwargs),
        indent=2,
    )


def _optional_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _optional_str(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


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
    provenance_raw = raw.get("provenance", {})
    provenance = provenance_raw if isinstance(provenance_raw, dict) else {}
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
        dataset_id=_optional_str(provenance.get("datasetId", raw.get("datasetId"))),
        dataset_role=str(
            provenance.get("datasetRole", raw.get("datasetRole", "unassigned"))
            or "unassigned"
        ),
        source_file_id=_optional_str(
            provenance.get("sourceFileId", raw.get("sourceFileId"))
        ),
        collected_at=_optional_str(
            provenance.get("collectedAt", raw.get("collectedAt"))
        ),
        operator=_optional_str(provenance.get("operator", raw.get("operator"))),
        profile_id=_optional_str(provenance.get("profileId", raw.get("profileId"))),
        profile_frozen_at=_optional_str(
            provenance.get("profileFrozenAt", raw.get("profileFrozenAt"))
        ),
        provenance_notes=str(provenance.get("notes", "")),
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


def _calibration_parameter_estimates(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    measurements: CalibrationExperimentMeasurements,
    comparisons: tuple[CalibrationExperimentComparison, ...],
    fit: dict[str, object],
    summary: dict[str, object],
) -> dict[str, object]:
    protocol_steps = {step.id: step for step in protocol.steps}
    z_ratios: list[float] = []
    slip_values: list[float] = []
    force_values: list[float] = []
    force_height_pairs: list[tuple[float, float]] = []

    for measured_step in measurements.steps:
        step = protocol_steps.get(measured_step.step_id)
        if step is None:
            continue
        command_cells = {command.cell for command in step.commands}
        direct_heights: list[float] = []
        neighbor_heights: list[float] = []
        for measurement in measured_step.cells:
            if measurement.pin_hole_slip_mm is not None:
                slip_values.append(float(measurement.pin_hole_slip_mm))
            if measurement.actuator_force_n is not None:
                force_values.append(float(measurement.actuator_force_n))
            if measurement.height_delta is not None:
                height_magnitude = abs(float(measurement.height_delta))
                if measurement.actuator_force_n is not None:
                    force_height_pairs.append(
                        (height_magnitude, abs(float(measurement.actuator_force_n)))
                    )
                if measurement.cell in command_cells:
                    direct_heights.append(height_magnitude)
                else:
                    neighbor_heights.append(height_magnitude)
        direct = _mean_optional(direct_heights)
        if direct is not None and direct > 1e-12:
            for neighbor in neighbor_heights:
                z_ratios.append(neighbor / direct)

    force_height_denominator = sum(height**2 for height, _ in force_height_pairs)
    force_height_slope = (
        sum(height * force for height, force in force_height_pairs)
        / force_height_denominator
        if force_height_denominator > 1e-12
        else None
    )
    z_estimate = _mean_optional(z_ratios)
    alpha_fit = fit.get("alpha", {}) if isinstance(fit, dict) else {}
    height_fit = fit.get("height", {}) if isinstance(fit, dict) else {}
    return {
        "schema": "rad-sim.calibration-parameter-estimates.v1",
        "method": (
            "conservative parameter bookkeeping from completed calibration "
            "experiment measurements; estimates are not treated as physical "
            "laws without separate validation"
        ),
        "claimLabels": {
            "alphaResponseGain": "simulator-derived empirical fit",
            "heightResponseGain": "simulator-derived empirical fit",
            "zCouplingGain": "experimentally unvalidated empirical estimate",
            "verticalFreePlay": "experimentally unvalidated physical assumption",
            "backlash": "not identifiable from current measurement set",
            "contactProxy": "experimentally unvalidated physical assumption",
        },
        "estimates": {
            "alphaResponse": {
                "suggestedGain": alpha_fit.get("suggestedGain"),
                "suggestedBias": alpha_fit.get("suggestedBias"),
                "sampleCount": alpha_fit.get("sampleCount", 0),
                "gainIdentifiable": alpha_fit.get("gainIdentifiable", False),
                "rmsRawError": alpha_fit.get("rmsRawError"),
                "rmsResidual": alpha_fit.get("rmsResidual"),
                "status": (
                    "estimated-from-results"
                    if int(alpha_fit.get("sampleCount") or 0) > 0
                    else "missing-measurements"
                ),
            },
            "heightResponse": {
                "suggestedGain": height_fit.get("suggestedGain"),
                "suggestedBias": height_fit.get("suggestedBias"),
                "sampleCount": height_fit.get("sampleCount", 0),
                "gainIdentifiable": height_fit.get("gainIdentifiable", False),
                "rmsRawError": height_fit.get("rmsRawError"),
                "rmsResidual": height_fit.get("rmsResidual"),
                "status": (
                    "estimated-from-results"
                    if int(height_fit.get("sampleCount") or 0) > 0
                    else "missing-measurements"
                ),
            },
            "zCouplingGain": {
                "estimate": z_estimate,
                "configured": config.z_coupling_gain,
                "sampleCount": len(z_ratios),
                "method": "mean neighbor/direct measured height ratio on z-actuated protocol steps",
                "status": (
                    "estimated-from-results" if z_ratios else "missing-pair-z-measurements"
                ),
            },
            "verticalFreePlay": {
                "configuredModelUnits": config.pin_hole_clearance,
                "meanMeasuredSlipMm": _mean_optional(slip_values),
                "sampleCount": len(slip_values),
                "status": (
                    "measured-slip-summary" if slip_values else "configured-clearance-only"
                ),
            },
            "backlash": {
                "configured": config.backlash,
                "estimate": None,
                "sampleCount": 0,
                "status": "not-identifiable-from-current-results",
                "requiredExperiment": "one-cell alpha command sweep across the dead-zone with reversal/hysteresis measurements",
            },
            "contactProxy": {
                "meanActuatorForceN": _mean_optional(force_values),
                "forceHeightSlopeNPerModelUnit": force_height_slope,
                "sampleCount": len(force_height_pairs),
                "status": (
                    "force-displacement-proxy" if force_height_pairs else "missing-force-measurements"
                ),
            },
        },
        "quality": {
            "comparisonCount": len(comparisons),
            "measuredCellCount": summary.get("measuredCellCount"),
            "missingObservationCount": summary.get("missingObservationCount"),
            "fitResidualMaxCombinedError": summary.get("fitResidualMaxCombinedError"),
            "maxCombinedError": summary.get("maxCombinedError"),
        },
        "suggestedModelUpdates": {
            "zCouplingGain": z_estimate,
            "alphaResponseGain": alpha_fit.get("suggestedGain"),
            "heightResponseGain": height_fit.get("suggestedGain"),
        },
        "limitations": [
            "The z-coupling estimate is a measured response ratio, not a derivation of the simulator z-coupling law.",
            "Backlash is not identifiable without a command sweep across the dead-zone and reversal path.",
            "The contact proxy uses force/displacement summaries only and does not calibrate friction or rigid-body contact.",
        ],
    }


def _estimate_section(
    parameter_estimates: dict[str, object], name: str
) -> dict[str, object]:
    estimates = parameter_estimates.get("estimates", {})
    if not isinstance(estimates, dict):
        return {}
    section = estimates.get(name, {})
    return section if isinstance(section, dict) else {}


def _profile_update(
    *,
    name: str,
    config_field: str | None,
    current: float | None,
    proposed: float | None,
    sample_count: int,
    source_status: str,
    claim_label: str,
    safe_to_apply: bool,
    reason: str,
) -> dict[str, object]:
    delta = (
        float(proposed - current)
        if proposed is not None and current is not None
        else None
    )
    return {
        "name": name,
        "configField": config_field,
        "current": current,
        "proposed": proposed,
        "delta": delta,
        "sampleCount": int(sample_count),
        "sourceStatus": source_status,
        "claimLabel": claim_label,
        "safeToApply": bool(safe_to_apply),
        "reason": reason,
    }


def _calibration_model_profile(
    config: LatticeConfig,
    parameter_estimates: dict[str, object],
    *,
    source_report_schema: str,
    min_samples: int = 1,
) -> dict[str, object]:
    """Build a conservative model-profile artifact from calibration estimates."""

    claims = parameter_estimates.get("claimLabels", {})
    if not isinstance(claims, dict):
        claims = {}
    quality = parameter_estimates.get("quality", {})
    if not isinstance(quality, dict):
        quality = {}

    z = _estimate_section(parameter_estimates, "zCouplingGain")
    z_estimate = _finite_or_none(z.get("estimate"))
    z_sample_count = int(z.get("sampleCount") or 0)
    z_status = str(z.get("status", "unknown"))
    z_safe = (
        z_estimate is not None
        and 0.0 <= z_estimate <= 1.0
        and z_sample_count >= min_samples
        and z_status == "estimated-from-results"
    )
    updates = [
        _profile_update(
            name="zCouplingGain",
            config_field="z_coupling_gain",
            current=float(config.z_coupling_gain),
            proposed=z_estimate,
            sample_count=z_sample_count,
            source_status=z_status,
            claim_label=str(
                claims.get(
                    "zCouplingGain",
                    "experimentally unvalidated empirical estimate",
                )
            ),
            safe_to_apply=z_safe,
            reason=(
                "bounded measured neighbor/direct z-response ratio"
                if z_safe
                else "requires a finite z-coupling estimate in [0, 1] with enough pair measurements"
            ),
        )
    ]

    diagnostic_names = (
        ("alphaResponseGain", "alphaResponse", "alphaResponseGain"),
        ("heightResponseGain", "heightResponse", "heightResponseGain"),
        ("verticalFreePlay", "verticalFreePlay", "verticalFreePlay"),
        ("backlash", "backlash", "backlash"),
        ("contactProxy", "contactProxy", "contactProxy"),
    )
    for name, section_name, claim_name in diagnostic_names:
        section = _estimate_section(parameter_estimates, section_name)
        proposed = _finite_or_none(
            section.get("suggestedGain")
            if section_name in {"alphaResponse", "heightResponse"}
            else section.get("meanMeasuredSlipMm")
            if section_name == "verticalFreePlay"
            else section.get("forceHeightSlopeNPerModelUnit")
            if section_name == "contactProxy"
            else section.get("estimate")
        )
        updates.append(
            _profile_update(
                name=name,
                config_field=None,
                current=None,
                proposed=proposed,
                sample_count=int(section.get("sampleCount") or 0),
                source_status=str(section.get("status", "unknown")),
                claim_label=str(claims.get(claim_name, "diagnostic only")),
                safe_to_apply=False,
                reason="diagnostic in v1; no corresponding LatticeConfig field is mutated",
            )
        )

    safe_updates = [update for update in updates if update["safeToApply"]]
    blockers = [
        {
            "name": update["name"],
            "reason": update["reason"],
            "sourceStatus": update["sourceStatus"],
        }
        for update in updates
        if not update["safeToApply"]
    ]
    return {
        "schema": "rad-sim.calibration-model-profile.v1",
        "sourceReportSchema": source_report_schema,
        "method": (
            "claim-labeled fitted-profile artifact; safe updates may adjust "
            "bounded simulator parameters, while non-configurable fitted "
            "response gains remain diagnostic evidence"
        ),
        "configBefore": {
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "pinHoleClearance": config.pin_hole_clearance,
        },
        "recommendedUpdates": updates,
        "safeUpdateCount": len(safe_updates),
        "diagnosticOnlyCount": len(updates) - len(safe_updates),
        "blockers": blockers,
        "quality": {
            "measuredCellCount": quality.get("measuredCellCount"),
            "missingObservationCount": quality.get("missingObservationCount"),
            "fitResidualMaxCombinedError": quality.get("fitResidualMaxCombinedError"),
            "maxCombinedError": quality.get("maxCombinedError"),
            "minSamples": int(min_samples),
        },
        "claimLabels": {
            "profile": "simulator-derived empirical law",
            "appliedConfigUpdates": "experimentally unvalidated empirical estimate",
            "diagnosticOnlyUpdates": "not a simulator mutation in v1",
        },
        "limitations": [
            "Applying this profile updates only bounded simulator parameters with direct config fields.",
            "Alpha and height response gains are recorded for model review but do not change the kinematic law in v1.",
            "Backlash, pin-hole clearance, friction, and rigid-body contact remain separate calibration targets.",
        ],
    }


def calibration_model_profile_from_report(
    config: LatticeConfig,
    report: dict[str, object],
    *,
    min_samples: int = 1,
) -> dict[str, object]:
    parameter_estimates = report.get("parameterEstimates", report)
    if not isinstance(parameter_estimates, dict):
        raise ValueError("calibration report must include parameterEstimates")
    source_schema = str(report.get("schema", parameter_estimates.get("schema", "")))
    return _calibration_model_profile(
        config,
        parameter_estimates,
        source_report_schema=source_schema,
        min_samples=min_samples,
    )


def apply_calibration_model_profile(
    config: LatticeConfig,
    profile: dict[str, object],
) -> tuple[LatticeConfig, dict[str, object]]:
    if profile.get("schema") != "rad-sim.calibration-model-profile.v1":
        raise ValueError("unsupported calibration model profile schema")
    updates = profile.get("recommendedUpdates", [])
    if not isinstance(updates, list):
        raise ValueError("calibration model profile requires recommendedUpdates")
    next_config = config
    applied: list[dict[str, object]] = []
    skipped: list[dict[str, object]] = []
    for update in updates:
        if not isinstance(update, dict):
            continue
        field = update.get("configField")
        proposed = _finite_or_none(update.get("proposed"))
        if update.get("safeToApply") is not True or proposed is None:
            skipped.append(
                {
                    "name": update.get("name"),
                    "configField": field,
                    "reason": update.get("reason", "not marked safe to apply"),
                }
            )
            continue
        if field == "z_coupling_gain":
            if not 0.0 <= proposed <= 1.0:
                skipped.append(
                    {
                        "name": update.get("name"),
                        "configField": field,
                        "reason": "z_coupling_gain must stay in [0, 1]",
                    }
                )
                continue
            next_config = replace(next_config, z_coupling_gain=proposed)
            applied.append(
                {
                    "name": update.get("name"),
                    "configField": field,
                    "value": proposed,
                }
            )
        else:
            skipped.append(
                {
                    "name": update.get("name"),
                    "configField": field,
                    "reason": "no safe v1 config mutation is defined for this field",
                }
            )
    audit = {
        "schema": "rad-sim.calibration-model-profile-application.v1",
        "sourceProfileSchema": profile.get("schema"),
        "appliedUpdates": applied,
        "skippedUpdates": skipped,
        "configBefore": {
            "zCouplingGain": config.z_coupling_gain,
        },
        "configAfter": {
            "zCouplingGain": next_config.z_coupling_gain,
        },
        "claimLabel": "simulator configuration update from empirical calibration profile",
    }
    return next_config, audit


def _comparison_summary_metrics(report: dict[str, object]) -> dict[str, object]:
    summary = report.get("summary", {})
    if not isinstance(summary, dict):
        summary = {}
    metrics = {
        "maxCombinedError": _finite_or_none(summary.get("maxCombinedError")),
        "fitResidualMaxCombinedError": _finite_or_none(
            summary.get("fitResidualMaxCombinedError")
        ),
        "heightRmseMean": _finite_or_none(summary.get("heightRmseMean")),
        "alphaRmseMean": _finite_or_none(summary.get("alphaRmseMean")),
        "missingObservationCount": int(summary.get("missingObservationCount") or 0),
        "measuredCellCount": int(summary.get("measuredCellCount") or 0),
    }
    finite_residuals = [
        value
        for key in (
            "maxCombinedError",
            "fitResidualMaxCombinedError",
            "heightRmseMean",
            "alphaRmseMean",
        )
        if (value := metrics[key]) is not None
    ]
    metrics["residualScore"] = (
        float(sum(finite_residuals)) if finite_residuals else None
    )
    return metrics


def _residual_delta(
    before: dict[str, object], after: dict[str, object]
) -> dict[str, object]:
    delta: dict[str, object] = {}
    for key in (
        "maxCombinedError",
        "fitResidualMaxCombinedError",
        "heightRmseMean",
        "alphaRmseMean",
        "residualScore",
    ):
        before_value = before.get(key)
        after_value = after.get(key)
        delta[key] = (
            float(after_value - before_value)
            if before_value is not None and after_value is not None
            else None
        )
    delta["missingObservationCount"] = int(after["missingObservationCount"]) - int(
        before["missingObservationCount"]
    )
    return delta


def _measurement_provenance(measurements: CalibrationExperimentMeasurements) -> dict[str, object]:
    return {
        "schema": "rad-sim.calibration-dataset-provenance.v1",
        "datasetId": measurements.dataset_id,
        "datasetRole": measurements.dataset_role,
        "sourceFileId": measurements.source_file_id,
        "collectedAt": measurements.collected_at,
        "operator": measurements.operator,
        "profileId": measurements.profile_id,
        "profileFrozenAt": measurements.profile_frozen_at,
        "notes": measurements.provenance_notes,
    }


def _profile_identifier(profile: dict[str, object]) -> str | None:
    candidate = _optional_str(profile.get("profileId"))
    if candidate is not None:
        return candidate
    source_schema = _optional_str(profile.get("sourceReportSchema"))
    safe_count = profile.get("safeUpdateCount")
    return f"{source_schema or 'profile'}:safe={safe_count}" if safe_count is not None else None


def _split_independence_metadata(
    fit_measurements: CalibrationExperimentMeasurements,
    holdout_measurements: CalibrationExperimentMeasurements,
    profile: dict[str, object],
    *,
    fit_measured_cells: object,
    holdout_measured_cells: object,
) -> dict[str, object]:
    fit_id = fit_measurements.dataset_id
    holdout_id = holdout_measurements.dataset_id
    source_file_overlap = (
        fit_measurements.source_file_id is not None
        and fit_measurements.source_file_id == holdout_measurements.source_file_id
    )
    dataset_overlap = fit_id is not None and fit_id == holdout_id
    known_overlap_count = 1 if dataset_overlap or source_file_overlap else 0
    profile_id = _profile_identifier(profile)
    profile_frozen = holdout_measurements.profile_frozen_at is not None
    profile_matches = (
        profile_id is not None
        and holdout_measurements.profile_id is not None
        and holdout_measurements.profile_id == profile_id
    )
    role_fit = fit_measurements.dataset_role == "fit"
    role_holdout = holdout_measurements.dataset_role == "holdout"
    required = {
        "fitDatasetId": fit_id is not None,
        "holdoutDatasetId": holdout_id is not None,
        "distinctDatasetIds": fit_id is not None and holdout_id is not None and fit_id != holdout_id,
        "distinctSourceFileIds": (
            fit_measurements.source_file_id is not None
            and holdout_measurements.source_file_id is not None
            and fit_measurements.source_file_id != holdout_measurements.source_file_id
        ),
        "fitRoleMarked": role_fit,
        "holdoutRoleMarked": role_holdout,
        "profileFrozenBeforeHoldout": profile_frozen,
        "holdoutProfileMatchesFrozenProfile": profile_matches,
    }
    independence_ready = all(required.values()) and known_overlap_count == 0
    missing = [key for key, value in required.items() if not value]
    return {
        "schema": "rad-sim.calibration-train-holdout-split.v1",
        "fitStepCount": len(fit_measurements.steps),
        "holdoutStepCount": len(holdout_measurements.steps),
        "fitMeasuredCellCount": fit_measured_cells,
        "holdoutMeasuredCellCount": holdout_measured_cells,
        "fitDataset": _measurement_provenance(fit_measurements),
        "holdoutDataset": _measurement_provenance(holdout_measurements),
        "profileId": profile_id,
        "knownOverlapCount": known_overlap_count,
        "profileFrozenBeforeHoldout": profile_frozen,
        "holdoutProfileMatchesFrozenProfile": profile_matches,
        "independenceReady": bool(independence_ready),
        "independenceStatus": (
            "documented-independent-holdout"
            if independence_ready
            else "requires-external-bench-protocol"
        ),
        "missingEvidence": missing,
        "requiredEvidence": required,
        "claimLabel": "experimentally unvalidated split metadata",
    }


def calibration_model_profile_residual_comparison(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    measurements: CalibrationExperimentMeasurements,
    profile: dict[str, object] | None = None,
    *,
    tolerance: float = 1e-9,
    min_improvement: float = 0.0,
) -> dict[str, object]:
    """Compare calibration residuals before and after applying a model profile."""

    before_report = calibration_experiment_comparison_report(
        config,
        protocol,
        measurements,
        tolerance=tolerance,
    )
    resolved_profile = profile or before_report["modelProfile"]
    next_config, application = apply_calibration_model_profile(config, resolved_profile)
    after_report = calibration_experiment_comparison_report(
        next_config,
        protocol,
        measurements,
        tolerance=tolerance,
    )
    before_metrics = _comparison_summary_metrics(before_report)
    after_metrics = _comparison_summary_metrics(after_report)
    delta = _residual_delta(before_metrics, after_metrics)
    score_delta = delta.get("residualScore")
    missing_delta = int(delta.get("missingObservationCount") or 0)
    improved = (
        score_delta is not None
        and float(score_delta) <= -float(min_improvement)
        and missing_delta <= 0
    )
    return {
        "schema": "rad-sim.calibration-model-profile-residual-comparison.v1",
        "profileSchema": resolved_profile.get("schema"),
        "applicationSchema": application["schema"],
        "application": application,
        "before": {
            "grid": before_report.get("grid"),
            "summary": before_report.get("summary"),
            "metrics": before_metrics,
        },
        "after": {
            "grid": after_report.get("grid"),
            "summary": after_report.get("summary"),
            "metrics": after_metrics,
        },
        "delta": delta,
        "improved": bool(improved),
        "selectionScore": after_metrics.get("residualScore"),
        "claimLabels": {
            "comparison": "simulator-derived empirical law",
            "profileApplication": application["claimLabel"],
        },
        "limitations": [
            "Residual improvement is evaluated against the same measurement file used to fit the profile.",
            "A lower simulator residual is not independent physical validation.",
            "Profile comparison does not calibrate backlash, friction, or rigid-body contact.",
        ],
    }


def calibration_model_profile_holdout_validation(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    fit_measurements: CalibrationExperimentMeasurements,
    holdout_measurements: CalibrationExperimentMeasurements,
    profile: dict[str, object] | None = None,
    *,
    tolerance: float = 1e-9,
    min_improvement: float = 0.0,
) -> dict[str, object]:
    """Validate a calibration profile against separate holdout measurements."""

    fit_report = calibration_experiment_comparison_report(
        config,
        protocol,
        fit_measurements,
        tolerance=tolerance,
    )
    resolved_profile = profile or fit_report["modelProfile"]
    fit_comparison = calibration_model_profile_residual_comparison(
        config,
        protocol,
        fit_measurements,
        resolved_profile,
        tolerance=tolerance,
        min_improvement=min_improvement,
    )
    holdout_comparison = calibration_model_profile_residual_comparison(
        config,
        protocol,
        holdout_measurements,
        resolved_profile,
        tolerance=tolerance,
        min_improvement=min_improvement,
    )

    def _passes(comparison: dict[str, object]) -> bool:
        delta = comparison.get("delta", {})
        application = comparison.get("application", {})
        applied = application.get("appliedUpdates", []) if isinstance(application, dict) else []
        score_delta = _finite_or_none(
            delta.get("residualScore") if isinstance(delta, dict) else None
        )
        missing_delta = (
            int(delta.get("missingObservationCount") or 0)
            if isinstance(delta, dict)
            else 0
        )
        return (
            isinstance(applied, list)
            and len(applied) > 0
            and score_delta is not None
            and score_delta <= -float(min_improvement)
            and missing_delta <= 0
    )

    fit_pass = _passes(fit_comparison)
    holdout_pass = _passes(holdout_comparison)
    fit_metrics = fit_comparison["after"]["metrics"]
    holdout_metrics = holdout_comparison["after"]["metrics"]
    split_metadata = _split_independence_metadata(
        fit_measurements,
        holdout_measurements,
        resolved_profile,
        fit_measured_cells=fit_metrics.get("measuredCellCount"),
        holdout_measured_cells=holdout_metrics.get("measuredCellCount"),
    )
    residual_pass = bool(fit_pass and holdout_pass)
    independent_pass = bool(residual_pass and split_metadata["independenceReady"])
    return {
        "schema": "rad-sim.calibration-model-profile-holdout-validation.v1",
        "profileSchema": resolved_profile.get("schema"),
        "modelProfile": resolved_profile,
        "fitComparisonSchema": fit_comparison["schema"],
        "holdoutComparisonSchema": holdout_comparison["schema"],
        "fit": fit_comparison,
        "holdout": holdout_comparison,
        "holdoutPass": bool(holdout_pass),
        "fitPass": bool(fit_pass),
        "residualValidationPass": residual_pass,
        "independentValidationPass": independent_pass,
        "splitMetadata": split_metadata,
        "provenanceWarnings": [
            f"missing or invalid split evidence: {field}"
            for field in split_metadata["missingEvidence"]
        ],
        "rule": {
            "name": "fit-profile-holdout-residual-nonincrease",
            "minImprovement": float(min_improvement),
            "requiresFitImprovement": True,
            "requiresHoldoutImprovement": True,
            "requiresAppliedUpdate": True,
            "requiresMissingObservationNonincrease": True,
        },
        "claimLabels": {
            "fit": "simulator-derived empirical law",
            "holdout": "holdout simulator validation",
            "independentValidationPass": "experimentally unvalidated until holdout measurements are independently collected",
        },
        "limitations": [
            "Holdout validation is meaningful only if the holdout file was not used to fit the profile.",
            "Passing holdout residual checks supports simulator calibration but does not prove contact, friction, or material laws.",
            "The current v1 profile mutates only bounded simulator fields such as z-coupling.",
        ],
    }


def _calibration_csv_scalar(value: object) -> str:
    if value is None:
        text = ""
    elif isinstance(value, bool):
        text = "true" if value else "false"
    else:
        text = str(value)
    if any(char in text for char in [",", '"', "\n", "\r"]):
        return '"' + text.replace('"', '""') + '"'
    return text


def _calibration_metric(
    comparison: dict[str, object],
    stage: str,
    key: str,
) -> object:
    stage_payload = comparison.get(stage, {})
    metrics = stage_payload.get("metrics", {}) if isinstance(stage_payload, dict) else {}
    return metrics.get(key) if isinstance(metrics, dict) else None


def export_calibration_model_profile_holdout_validation_csv(
    validation: dict[str, object],
) -> str:
    """Export a train/holdout profile validation as a bench summary CSV."""

    header = [
        "schema",
        "dataset",
        "pass",
        "independent_validation_pass",
        "residual_validation_pass",
        "profile_schema",
        "fit_dataset_id",
        "holdout_dataset_id",
        "applied_updates",
        "z_coupling_before",
        "z_coupling_after",
        "before_residual_score",
        "after_residual_score",
        "residual_score_delta",
        "before_missing_observations",
        "after_missing_observations",
        "missing_observation_delta",
        "before_measured_cells",
        "after_measured_cells",
        "profile_frozen_before_holdout",
        "known_overlap_count",
        "independence_status",
        "independence_ready",
        "missing_evidence",
        "claim_label",
    ]
    split = validation.get("splitMetadata", {})
    if not isinstance(split, dict):
        split = {}
    rows: list[list[object]] = [header]
    for dataset, key, pass_key in (
        ("fit", "fit", "fitPass"),
        ("holdout", "holdout", "holdoutPass"),
    ):
        comparison = validation.get(key, {})
        if not isinstance(comparison, dict):
            comparison = {}
        application = comparison.get("application", {})
        applied = application.get("appliedUpdates", []) if isinstance(application, dict) else []
        config_before = application.get("configBefore", {}) if isinstance(application, dict) else {}
        config_after = application.get("configAfter", {}) if isinstance(application, dict) else {}
        delta = comparison.get("delta", {})
        if not isinstance(delta, dict):
            delta = {}
        rows.append(
            [
                validation.get("schema", ""),
                dataset,
                validation.get(pass_key, ""),
                validation.get("independentValidationPass", ""),
                validation.get("residualValidationPass", ""),
                validation.get("profileSchema", ""),
                (split.get("fitDataset") or {}).get("datasetId")
                if isinstance(split.get("fitDataset"), dict)
                else "",
                (split.get("holdoutDataset") or {}).get("datasetId")
                if isinstance(split.get("holdoutDataset"), dict)
                else "",
                len(applied) if isinstance(applied, list) else 0,
                config_before.get("zCouplingGain") if isinstance(config_before, dict) else "",
                config_after.get("zCouplingGain") if isinstance(config_after, dict) else "",
                _calibration_metric(comparison, "before", "residualScore"),
                _calibration_metric(comparison, "after", "residualScore"),
                delta.get("residualScore"),
                _calibration_metric(comparison, "before", "missingObservationCount"),
                _calibration_metric(comparison, "after", "missingObservationCount"),
                delta.get("missingObservationCount"),
                _calibration_metric(comparison, "before", "measuredCellCount"),
                _calibration_metric(comparison, "after", "measuredCellCount"),
                split.get("profileFrozenBeforeHoldout"),
                split.get("knownOverlapCount"),
                split.get("independenceStatus"),
                split.get("independenceReady"),
                ";".join(str(item) for item in split.get("missingEvidence", []))
                if isinstance(split.get("missingEvidence"), list)
                else "",
                split.get("claimLabel"),
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row) for row in rows
    )


def calibration_bench_execution_validation(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    fit_measurements: CalibrationExperimentMeasurements,
    holdout_measurements: CalibrationExperimentMeasurements,
    profile: dict[str, object] | None = None,
    *,
    tolerance: float = 1e-9,
    min_improvement: float = 0.0,
) -> dict[str, object]:
    """Validate filled fit/holdout calibration packet measurements."""

    fit_report = calibration_experiment_comparison_report(
        config,
        protocol,
        fit_measurements,
        tolerance=tolerance,
    )
    resolved_profile = profile or fit_report["modelProfile"]
    holdout_validation = calibration_model_profile_holdout_validation(
        config,
        protocol,
        fit_measurements,
        holdout_measurements,
        resolved_profile,
        tolerance=tolerance,
        min_improvement=min_improvement,
    )
    split = holdout_validation.get("splitMetadata", {})
    split_metadata = split if isinstance(split, dict) else {}
    application = holdout_validation.get("fit", {})
    fit_application = (
        application.get("application", {})
        if isinstance(application, dict)
        else {}
    )
    applied_updates = (
        fit_application.get("appliedUpdates", [])
        if isinstance(fit_application, dict)
        else []
    )
    missing_evidence = split_metadata.get("missingEvidence", [])
    missing_count = len(missing_evidence) if isinstance(missing_evidence, list) else 0
    execution_validation_pass = bool(
        holdout_validation.get("independentValidationPass") and missing_count == 0
    )
    status = (
        "documented-independent-validation"
        if execution_validation_pass
        else "residual-pass-needs-provenance"
        if holdout_validation.get("residualValidationPass")
        else "residual-validation-failed"
    )
    return {
        "schema": "rad-sim.calibration-bench-execution-validation.v1",
        "method": (
            "filled fit/holdout calibration packet ingestion using bounded "
            "model-profile residual replay and split-provenance checks"
        ),
        "inputSchemas": {
            "protocol": protocol.schema,
            "fitResults": "rad-sim.calibration-experiment-results.v1",
            "holdoutResults": "rad-sim.calibration-experiment-results.v1",
            "datasetProvenance": "rad-sim.calibration-dataset-provenance.v1",
            "holdoutValidation": "rad-sim.calibration-model-profile-holdout-validation.v1",
        },
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "cellSize": config.cell_size,
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
            "pinHoleClearance": config.pin_hole_clearance,
        },
        "protocol": {
            "schema": protocol.schema,
            "hardwareProfile": protocol.hardware_profile_name,
            "stepCount": len(protocol.steps),
            "repeatCountTotal": sum(step.repeat_count for step in protocol.steps),
            "centerCell": _cell_dict(protocol.center_cell),
        },
        "fitDataset": _measurement_provenance(fit_measurements),
        "holdoutDataset": _measurement_provenance(holdout_measurements),
        "fitComparisonReport": fit_report,
        "modelProfile": resolved_profile,
        "holdoutValidation": holdout_validation,
        "holdoutValidationCsv": export_calibration_model_profile_holdout_validation_csv(
            holdout_validation
        ),
        "summary": {
            "status": status,
            "fitStepCount": len(fit_measurements.steps),
            "holdoutStepCount": len(holdout_measurements.steps),
            "appliedUpdateCount": len(applied_updates)
            if isinstance(applied_updates, list)
            else 0,
            "residualValidationPass": bool(
                holdout_validation.get("residualValidationPass")
            ),
            "independentValidationPass": bool(
                holdout_validation.get("independentValidationPass")
            ),
            "executionValidationPass": execution_validation_pass,
            "missingEvidenceCount": missing_count,
            "independenceStatus": split_metadata.get("independenceStatus"),
            "knownOverlapCount": split_metadata.get("knownOverlapCount"),
        },
        "formalization": {
            "targetId": "calibration_bench_executed_validation_gate",
            "leanStructure": "Mechanics.CalibrationBenchExecutedValidationNat",
            "leanPredicate": "calibrationBenchExecutedValidationReadyNat",
            "schema": "rad-sim.calibration-bench-execution-validation.v1",
        },
        "claimLabels": {
            "executionValidation": "bench-data residual validation bookkeeping",
            "modelProfile": "simulator-derived empirical law",
            "independentValidationPass": "experimentally unvalidated unless provenance is externally true",
            "physicalAccuracy": "experimentally unvalidated physical assumption",
        },
        "limitations": [
            "The report can detect missing provenance fields but cannot prove the lab followed the protocol.",
            "Residual agreement supports simulator-profile bookkeeping only; it does not validate contact, friction, stiffness, gravity, or material laws.",
            "The current bounded profile application mutates only implemented simulator fields.",
        ],
    }


def export_calibration_bench_execution_validation_json(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    fit_measurements: CalibrationExperimentMeasurements,
    holdout_measurements: CalibrationExperimentMeasurements,
    profile: dict[str, object] | None = None,
    *,
    tolerance: float = 1e-9,
    min_improvement: float = 0.0,
) -> str:
    return json.dumps(
        calibration_bench_execution_validation(
            config,
            protocol,
            fit_measurements,
            holdout_measurements,
            profile,
            tolerance=tolerance,
            min_improvement=min_improvement,
        ),
        indent=2,
    )


def export_calibration_bench_execution_validation_csv(
    validation: dict[str, object],
) -> str:
    summary = validation.get("summary", {})
    if not isinstance(summary, dict):
        summary = {}
    fit_dataset = validation.get("fitDataset", {})
    holdout_dataset = validation.get("holdoutDataset", {})
    if not isinstance(fit_dataset, dict):
        fit_dataset = {}
    if not isinstance(holdout_dataset, dict):
        holdout_dataset = {}
    holdout_validation = validation.get("holdoutValidation", {})
    split = (
        holdout_validation.get("splitMetadata", {})
        if isinstance(holdout_validation, dict)
        else {}
    )
    if not isinstance(split, dict):
        split = {}
    header = [
        "schema",
        "status",
        "execution_validation_pass",
        "residual_validation_pass",
        "independent_validation_pass",
        "fit_dataset_id",
        "holdout_dataset_id",
        "applied_updates",
        "fit_step_count",
        "holdout_step_count",
        "missing_evidence_count",
        "missing_evidence",
        "independence_status",
        "known_overlap_count",
        "claim_label",
    ]
    row = [
        validation.get("schema", ""),
        summary.get("status", ""),
        summary.get("executionValidationPass", ""),
        summary.get("residualValidationPass", ""),
        summary.get("independentValidationPass", ""),
        fit_dataset.get("datasetId", ""),
        holdout_dataset.get("datasetId", ""),
        summary.get("appliedUpdateCount", ""),
        summary.get("fitStepCount", ""),
        summary.get("holdoutStepCount", ""),
        summary.get("missingEvidenceCount", ""),
        ";".join(str(item) for item in split.get("missingEvidence", []))
        if isinstance(split.get("missingEvidence"), list)
        else "",
        summary.get("independenceStatus", ""),
        summary.get("knownOverlapCount", ""),
        (validation.get("claimLabels", {}) or {}).get("executionValidation")
        if isinstance(validation.get("claimLabels", {}), dict)
        else "",
    ]
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in values)
        for values in (header, row)
    )


def _physical_validation_vertical_proxy_counts(
    vertical_load_comparison_report: dict[str, object] | None,
) -> tuple[int, int]:
    if not isinstance(vertical_load_comparison_report, dict):
        return (0, 0)
    load_terms = 0
    contact_terms = 0
    for scenario in vertical_load_comparison_report.get("scenarios", []):
        if not isinstance(scenario, dict):
            continue
        validation = scenario.get("validation", {})
        if not isinstance(validation, dict):
            continue
        for side in ("intact", "removed"):
            side_payload = validation.get(side, {})
            if not isinstance(side_payload, dict):
                continue
            errors = side_payload.get("errors", {})
            measured = side_payload.get("measured", {})
            simulated = side_payload.get("simulatedObserved", {})
            if not (
                isinstance(errors, dict)
                and isinstance(measured, dict)
                and isinstance(simulated, dict)
            ):
                continue
            if (
                "signedLoadWork" in errors
                and "signedLoadWork" in measured
                and "signedLoadWork" in simulated
                and "loadWorkMagnitude" in errors
                and "loadWorkMagnitude" in measured
                and "loadWorkMagnitude" in simulated
            ):
                load_terms += 1
            if (
                "heightContactPenalty" in errors
                and "heightContactPenalty" in measured
                and "heightContactPenalty" in simulated
            ):
                contact_terms += 1
    return load_terms, contact_terms


def physical_validation_readiness_report(
    config: LatticeConfig,
    calibration_execution_validation: dict[str, object] | None = None,
    vertical_load_comparison_report: dict[str, object] | None = None,
) -> dict[str, object]:
    """Audit whether bench evidence is ready for physical validation review."""

    calibration_summary = (
        calibration_execution_validation.get("summary", {})
        if isinstance(calibration_execution_validation, dict)
        else {}
    )
    if not isinstance(calibration_summary, dict):
        calibration_summary = {}
    vertical_summary = (
        vertical_load_comparison_report.get("summary", {})
        if isinstance(vertical_load_comparison_report, dict)
        else {}
    )
    if not isinstance(vertical_summary, dict):
        vertical_summary = {}
    load_proxy_terms, contact_proxy_terms = _physical_validation_vertical_proxy_counts(
        vertical_load_comparison_report
    )
    calibration_pass = bool(calibration_summary.get("executionValidationPass"))
    vertical_scenarios = int(vertical_summary.get("scenarioCount") or 0)
    vertical_missing = int(vertical_summary.get("missingMeasurementCount") or 0)
    vertical_pass = bool(vertical_summary.get("allScenariosPassTolerance"))
    clearance_configured = bool(
        config.pin_radius > 0
        and config.hole_radius >= config.pin_radius
        and config.pin_hole_clearance >= 0
    )
    missing_evidence: list[str] = []
    if not isinstance(calibration_execution_validation, dict):
        missing_evidence.append("calibrationBenchExecutionValidation")
    elif not calibration_pass:
        missing_evidence.append("calibrationExecutionValidationPass")
    if not isinstance(vertical_load_comparison_report, dict):
        missing_evidence.append("verticalLoadComparisonReport")
    else:
        if vertical_scenarios <= 0:
            missing_evidence.append("verticalLoadScenarioCount")
        if vertical_missing > 0:
            missing_evidence.append("verticalLoadMissingMeasurements")
        if not vertical_pass:
            missing_evidence.append("verticalLoadTolerancePass")
    if load_proxy_terms <= 0:
        missing_evidence.append("signedLoadWorkAndMagnitudeTerms")
    if contact_proxy_terms <= 0:
        missing_evidence.append("heightContactPenaltyTerms")
    if not clearance_configured:
        missing_evidence.append("pinHoleClearanceConfiguration")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.physical-validation-readiness.v1",
        "method": (
            "conservative evidence gate combining executed calibration validation "
            "with filled vertical-load work/contact comparison reports"
        ),
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "pinHoleClearance": config.pin_hole_clearance,
            "zCouplingGain": config.z_coupling_gain,
        },
        "inputSchemas": {
            "calibrationExecution": "rad-sim.calibration-bench-execution-validation.v1",
            "verticalLoadComparison": "rad-sim.vertical-load-energy-comparison-report.v1",
        },
        "evidence": {
            "calibrationExecutionValidationPass": calibration_pass,
            "fitStepCount": calibration_summary.get("fitStepCount", 0),
            "holdoutStepCount": calibration_summary.get("holdoutStepCount", 0),
            "verticalLoadScenarioCount": vertical_scenarios,
            "verticalLoadMissingMeasurementCount": vertical_missing,
            "verticalLoadAllScenariosPassTolerance": vertical_pass,
            "loadProxyTermCount": load_proxy_terms,
            "contactProxyTermCount": contact_proxy_terms,
            "pinHoleClearanceConfigured": clearance_configured,
        },
        "summary": {
            "status": (
                "ready-for-physical-claim-review"
                if ready
                else "needs-bench-evidence"
            ),
            "physicalValidationReady": ready,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "physical_validation_readiness_gate",
            "leanStructure": "Mechanics.PhysicalValidationReadinessNat",
            "leanPredicate": "physicalValidationReadyNat",
            "schema": "rad-sim.physical-validation-readiness.v1",
        },
        "claimLabels": {
            "readinessGate": "bench-data evidence bookkeeping",
            "loadProxy": "Lean-proven finite scaffold plus simulator-derived empirical law",
            "contactProxy": "experimentally unvalidated physical assumption",
            "physicalAccuracy": "experimentally unvalidated physical assumption",
        },
        "limitations": [
            "Readiness means the required evidence artifacts are present and internally passing.",
            "It does not prove calibrated rigid-body contact, friction, stiffness, gravity, or material law.",
            "Physical claims still require lab audit of fixtures, sensors, raw files, and uncertainty.",
        ],
    }


def export_physical_validation_readiness_json(
    config: LatticeConfig,
    calibration_execution_validation: dict[str, object] | None = None,
    vertical_load_comparison_report: dict[str, object] | None = None,
) -> str:
    return json.dumps(
        physical_validation_readiness_report(
            config,
            calibration_execution_validation,
            vertical_load_comparison_report,
        ),
        indent=2,
    )


def export_physical_validation_readiness_csv(report: dict[str, object]) -> str:
    summary = report.get("summary", {})
    evidence = report.get("evidence", {})
    if not isinstance(summary, dict):
        summary = {}
    if not isinstance(evidence, dict):
        evidence = {}
    header = [
        "schema",
        "status",
        "physical_validation_ready",
        "missing_evidence_count",
        "missing_evidence",
        "calibration_execution_pass",
        "fit_step_count",
        "holdout_step_count",
        "vertical_load_scenarios",
        "vertical_load_missing_measurements",
        "vertical_load_pass",
        "load_proxy_terms",
        "contact_proxy_terms",
        "pin_hole_clearance_configured",
        "claim_label",
    ]
    row = [
        report.get("schema", ""),
        summary.get("status", ""),
        summary.get("physicalValidationReady", ""),
        summary.get("missingEvidenceCount", ""),
        ";".join(str(item) for item in summary.get("missingEvidence", []))
        if isinstance(summary.get("missingEvidence"), list)
        else "",
        evidence.get("calibrationExecutionValidationPass", ""),
        evidence.get("fitStepCount", ""),
        evidence.get("holdoutStepCount", ""),
        evidence.get("verticalLoadScenarioCount", ""),
        evidence.get("verticalLoadMissingMeasurementCount", ""),
        evidence.get("verticalLoadAllScenariosPassTolerance", ""),
        evidence.get("loadProxyTermCount", ""),
        evidence.get("contactProxyTermCount", ""),
        evidence.get("pinHoleClearanceConfigured", ""),
        (report.get("claimLabels", {}) or {}).get("readinessGate")
        if isinstance(report.get("claimLabels", {}), dict)
        else "",
    ]
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in values)
        for values in (header, row)
    )


def contact_state_abstraction_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    contact_stiffness: float = 1.0,
    vertical_load_comparison_report: dict[str, object] | None = None,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Build a normalized pin-hole/contact-state abstraction for the lattice."""

    if contact_stiffness < 0:
        raise ValueError("contact_stiffness must be non-negative")
    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    clearance = config.pin_hole_clearance
    cells: list[dict[str, object]] = []
    active_count = 0
    removed_count = 0
    engaged_count = 0
    total_penalty = 0.0
    for row in range(config.rows):
        for col in range(config.cols):
            removed = bool(base.removed_mask[row, col])
            z_command = float(base.z_actuator_grid[row, col])
            penetration = 0.0 if removed else max(0.0, abs(z_command) - clearance)
            penalty = 0.5 * float(contact_stiffness) * penetration**2
            mode = (
                "removed"
                if removed
                else "engaged"
                if penetration > tolerance
                else "free-clearance"
            )
            if removed:
                removed_count += 1
            else:
                active_count += 1
            if mode == "engaged":
                engaged_count += 1
            total_penalty += penalty
            cells.append(
                {
                    "row": row,
                    "col": col,
                    "removed": removed,
                    "zCommand": z_command,
                    "pinRadius": config.pin_radius,
                    "holeRadius": config.hole_radius,
                    "clearance": clearance,
                    "penetration": penetration,
                    "contactPenalty": penalty,
                    "mode": mode,
                    "contactState": {
                        "pinPresent": not removed,
                        "holePresent": not removed,
                        "unilateralContactActive": mode == "engaged",
                        "clearanceExceeded": penetration > tolerance,
                    },
                }
            )
    bench_load_terms, bench_contact_terms = _physical_validation_vertical_proxy_counts(
        vertical_load_comparison_report
    )
    clearance_configured = bool(
        config.pin_radius > 0
        and config.hole_radius >= config.pin_radius
        and clearance >= 0
    )
    missing_evidence: list[str] = []
    if active_count <= 0:
        missing_evidence.append("activeBodies")
    if config.pin_radius <= 0:
        missing_evidence.append("pinRadius")
    if config.hole_radius < config.pin_radius:
        missing_evidence.append("holeRadiusAtLeastPinRadius")
    if not clearance_configured:
        missing_evidence.append("clearanceConfiguration")
    if contact_stiffness <= 0:
        missing_evidence.append("positiveContactStiffness")
    if not cells:
        missing_evidence.append("contactStateRecords")
    abstraction_ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.contact-state-abstraction.v1",
        "method": (
            "normalized per-cell pin-hole clearance and unilateral contact-state "
            "bookkeeping derived from z commands, removed topology, and optional "
            "vertical-load comparison evidence"
        ),
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "pinHoleClearance": clearance,
            "contactStiffness": float(contact_stiffness),
        },
        "topology": {
            "activeBodyCount": active_count,
            "removedBodyCount": removed_count,
            "pinCount": active_count,
            "holeCount": active_count,
            "clearancePairCount": active_count,
            "contactStateRecordCount": len(cells),
            "penaltyTermCount": active_count,
        },
        "summary": {
            "status": (
                "contact-state-abstraction-ready"
                if abstraction_ready
                else "needs-contact-abstraction-inputs"
            ),
            "contactStateAbstractionReady": abstraction_ready,
            "engagedContactCount": engaged_count,
            "totalContactPenalty": float(total_penalty),
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "benchEvidence": {
            "verticalLoadComparisonAttached": isinstance(
                vertical_load_comparison_report, dict
            ),
            "loadProxyTermCount": bench_load_terms,
            "contactProxyTermCount": bench_contact_terms,
        },
        "cells": cells,
        "formalization": {
            "targetId": "contact_state_abstraction_gate",
            "leanStructure": "Mechanics.ContactStateAbstractionNat",
            "leanPredicate": "contactStateAbstractionReadyNat",
            "schema": "rad-sim.contact-state-abstraction.v1",
        },
        "claimLabels": {
            "contactState": "graph/contact bookkeeping abstraction",
            "contactPenalty": "simulator-derived empirical law",
            "physicalContact": "experimentally unvalidated physical assumption",
        },
        "limitations": [
            "The contact mode is derived from normalized z command and clearance, not collision detection.",
            "The penalty is a unilateral contact proxy until contact stiffness and friction are measured.",
            "Removed cells are deleted topology records, not simulated detached rigid bodies.",
        ],
    }


def export_contact_state_abstraction_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    contact_stiffness: float = 1.0,
    vertical_load_comparison_report: dict[str, object] | None = None,
    tolerance: float = 1e-9,
) -> str:
    return json.dumps(
        contact_state_abstraction_report(
            config,
            state,
            contact_stiffness=contact_stiffness,
            vertical_load_comparison_report=vertical_load_comparison_report,
            tolerance=tolerance,
        ),
        indent=2,
    )


def export_contact_state_abstraction_csv(report: dict[str, object]) -> str:
    header = [
        "row",
        "col",
        "removed",
        "mode",
        "z_command",
        "pin_radius",
        "hole_radius",
        "clearance",
        "penetration",
        "contact_penalty",
    ]
    rows = [header]
    for cell in report.get("cells", []):
        if not isinstance(cell, dict):
            continue
        rows.append(
            [
                cell.get("row", ""),
                cell.get("col", ""),
                cell.get("removed", ""),
                cell.get("mode", ""),
                cell.get("zCommand", ""),
                cell.get("pinRadius", ""),
                cell.get("holeRadius", ""),
                cell.get("clearance", ""),
                cell.get("penetration", ""),
                cell.get("contactPenalty", ""),
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row) for row in rows
    )


def _contact_cell_key(cell: dict[str, object]) -> tuple[int, int] | None:
    try:
        return (int(cell.get("row", -1)), int(cell.get("col", -1)))
    except (TypeError, ValueError):
        return None


def _graph_edges_for_state(
    config: LatticeConfig,
    state: LatticeState,
) -> tuple[list[dict[str, object]], int, int, int]:
    edges: list[dict[str, object]] = []
    active_edges = 0
    deleted_edges = 0
    removed_incident_active_edges = 0
    for row in range(config.rows):
        for col in range(config.cols):
            for nr, nc in ((row, col + 1), (row + 1, col)):
                if nr >= config.rows or nc >= config.cols:
                    continue
                touches_removed = bool(
                    state.removed_mask[row, col] or state.removed_mask[nr, nc]
                )
                active = not touches_removed
                if active:
                    active_edges += 1
                else:
                    deleted_edges += 1
                if active and touches_removed:
                    removed_incident_active_edges += 1
                edges.append(
                    {
                        "from": _cell_dict((row, col)),
                        "to": _cell_dict((nr, nc)),
                        "active": active,
                        "touchesRemoved": touches_removed,
                        "status": "active" if active else "deleted-by-removal",
                    }
                )
    return edges, active_edges, deleted_edges, removed_incident_active_edges


def _support_from_state(state: LatticeState, tolerance: float) -> tuple[tuple[int, int], ...]:
    support: list[tuple[int, int]] = []
    rows, cols = state.alpha_grid.shape
    for row in range(rows):
        for col in range(cols):
            if (
                abs(float(state.actuator_grid[row, col])) > tolerance
                or abs(float(state.z_actuator_grid[row, col])) > tolerance
                or bool(state.locked_mask[row, col])
            ):
                support.append((row, col))
    return tuple(support)


def contact_graph_consistency_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    contact_report: dict[str, object] | None = None,
    group_support: Iterable[tuple[int, int]] | None = None,
    contact_stiffness: float = 1.0,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Check graph-deletion consistency of contact-state records and supports."""

    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    resolved_contact = contact_report or contact_state_abstraction_report(
        config,
        base,
        contact_stiffness=contact_stiffness,
        tolerance=tolerance,
    )
    contact_cells = (
        resolved_contact.get("cells", []) if isinstance(resolved_contact, dict) else []
    )
    contact_by_cell: dict[tuple[int, int], dict[str, object]] = {}
    for cell in contact_cells:
        if not isinstance(cell, dict):
            continue
        key = _contact_cell_key(cell)
        if key is not None:
            contact_by_cell[key] = cell
    edges, active_edges, deleted_edges, removed_incident_active_edges = (
        _graph_edges_for_state(config, base)
    )
    removed_active_contacts = 0
    active_contact_records = 0
    missing_active_contact_records = 0
    for row in range(config.rows):
        for col in range(config.cols):
            key = (row, col)
            record = contact_by_cell.get(key)
            removed = bool(base.removed_mask[row, col])
            if record is None:
                if not removed:
                    missing_active_contact_records += 1
                continue
            mode = str(record.get("mode", ""))
            contact_state = record.get("contactState", {})
            active_contact = bool(
                isinstance(contact_state, dict)
                and contact_state.get("unilateralContactActive")
            )
            if not removed:
                active_contact_records += 1
            if removed and (active_contact or mode == "engaged"):
                removed_active_contacts += 1
    support_cells = tuple(group_support) if group_support is not None else _support_from_state(base, tolerance)
    valid_support: list[tuple[int, int]] = []
    skipped_support: list[tuple[int, int]] = []
    support_records = 0
    active_support = 0
    removed_support = 0
    for row, col in support_cells:
        if not (0 <= row < config.rows and 0 <= col < config.cols):
            skipped_support.append((row, col))
            continue
        valid_support.append((row, col))
        if (row, col) in contact_by_cell:
            support_records += 1
        if base.removed_mask[row, col]:
            removed_support += 1
        else:
            active_support += 1
    missing_evidence: list[str] = []
    if not isinstance(resolved_contact, dict):
        missing_evidence.append("contactStateAbstraction")
    elif not resolved_contact.get("summary", {}).get("contactStateAbstractionReady"):
        missing_evidence.append("contactStateAbstractionReady")
    if int(np.count_nonzero(~base.removed_mask)) <= 0:
        missing_evidence.append("activeBodies")
    if missing_active_contact_records > 0:
        missing_evidence.append("activeContactRecords")
    if removed_incident_active_edges > 0:
        missing_evidence.append("removedIncidentActiveEdges")
    if removed_active_contacts > 0:
        missing_evidence.append("removedActiveContacts")
    if support_records < len(valid_support):
        missing_evidence.append("supportContactRecords")
    consistency_ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.contact-graph-consistency.v1",
        "method": (
            "checks that graph edges touching removed cells are deleted, "
            "contact records are inactive on removed cells, and support cells "
            "are represented in the contact-state abstraction"
        ),
        "inputSchemas": {
            "contactState": "rad-sim.contact-state-abstraction.v1",
            "state": "rad-sim lattice state with removed/contact support masks",
        },
        "graph": {
            "activeBodyCount": int(np.count_nonzero(~base.removed_mask)),
            "removedBodyCount": int(np.count_nonzero(base.removed_mask)),
            "totalEdgeCount": len(edges),
            "activeEdgeCount": active_edges,
            "deletedEdgeCount": deleted_edges,
            "removedIncidentActiveEdgeCount": removed_incident_active_edges,
            "edges": edges,
        },
        "contact": {
            "contactRecordCount": len(contact_by_cell),
            "activeContactRecordCount": active_contact_records,
            "missingActiveContactRecordCount": missing_active_contact_records,
            "removedActiveContactCount": removed_active_contacts,
            "engagedContactCount": resolved_contact.get("summary", {}).get("engagedContactCount", 0)
            if isinstance(resolved_contact, dict)
            else 0,
        },
        "support": {
            "supportCells": [_cell_dict(cell) for cell in valid_support],
            "skippedSupportCells": [_cell_dict(cell) for cell in skipped_support],
            "supportCellCount": len(valid_support),
            "activeSupportCellCount": active_support,
            "removedSupportCellCount": removed_support,
            "supportContactRecordCount": support_records,
        },
        "summary": {
            "status": (
                "contact-graph-consistent"
                if consistency_ready
                else "needs-contact-graph-evidence"
            ),
            "contactGraphConsistent": consistency_ready,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "contact_graph_consistency_gate",
            "leanStructure": "Mechanics.ContactGraphConsistencyNat",
            "leanPredicate": "contactGraphConsistentNat",
            "schema": "rad-sim.contact-graph-consistency.v1",
        },
        "claimLabels": {
            "graphDeletion": "Lean-proven graph deletion theorem scaffold",
            "contactRecords": "graph/contact bookkeeping abstraction",
            "support": "simulator-derived event support diagnostic",
            "physicalContact": "experimentally unvalidated physical assumption",
        },
        "limitations": [
            "Consistency only checks graph/contact bookkeeping against removed topology.",
            "It does not prove that inferred contact modes match hardware contact.",
            "Removed support is reported as lost support, not treated as an error by itself.",
        ],
    }


def export_contact_graph_consistency_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    contact_report: dict[str, object] | None = None,
    group_support: Iterable[tuple[int, int]] | None = None,
    contact_stiffness: float = 1.0,
    tolerance: float = 1e-9,
) -> str:
    return json.dumps(
        contact_graph_consistency_report(
            config,
            state,
            contact_report=contact_report,
            group_support=group_support,
            contact_stiffness=contact_stiffness,
            tolerance=tolerance,
        ),
        indent=2,
    )


def export_contact_graph_consistency_csv(report: dict[str, object]) -> str:
    graph = report.get("graph", {})
    contact = report.get("contact", {})
    support = report.get("support", {})
    summary = report.get("summary", {})
    if not isinstance(graph, dict):
        graph = {}
    if not isinstance(contact, dict):
        contact = {}
    if not isinstance(support, dict):
        support = {}
    if not isinstance(summary, dict):
        summary = {}
    header = [
        "schema",
        "status",
        "contact_graph_consistent",
        "active_bodies",
        "removed_bodies",
        "active_edges",
        "deleted_edges",
        "removed_incident_active_edges",
        "contact_records",
        "removed_active_contacts",
        "support_cells",
        "support_contact_records",
        "removed_support_cells",
        "missing_evidence",
    ]
    row = [
        report.get("schema", ""),
        summary.get("status", ""),
        summary.get("contactGraphConsistent", ""),
        graph.get("activeBodyCount", ""),
        graph.get("removedBodyCount", ""),
        graph.get("activeEdgeCount", ""),
        graph.get("deletedEdgeCount", ""),
        graph.get("removedIncidentActiveEdgeCount", ""),
        contact.get("contactRecordCount", ""),
        contact.get("removedActiveContactCount", ""),
        support.get("supportCellCount", ""),
        support.get("supportContactRecordCount", ""),
        support.get("removedSupportCellCount", ""),
        ";".join(str(item) for item in summary.get("missingEvidence", []))
        if isinstance(summary.get("missingEvidence"), list)
        else "",
    ]
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in values)
        for values in (header, row)
    )


def _state_command_support(
    state: LatticeState, tolerance: float
) -> tuple[tuple[int, int], ...]:
    support: list[tuple[int, int]] = []
    rows, cols = state.alpha_grid.shape
    for row in range(rows):
        for col in range(cols):
            if (
                abs(float(state.actuator_grid[row, col])) > tolerance
                or abs(float(state.z_actuator_grid[row, col])) > tolerance
            ):
                support.append((row, col))
    return tuple(support)


def _infer_realization_events(
    config: LatticeConfig, state: LatticeState, tolerance: float
) -> tuple[ProgrammableDiscontinuityEvent, ...]:
    events: list[ProgrammableDiscontinuityEvent] = []
    for row in range(config.rows):
        for col in range(config.cols):
            alpha = float(state.actuator_grid[row, col])
            z = float(state.z_actuator_grid[row, col])
            if abs(alpha) > tolerance or abs(z) > tolerance:
                events.append(
                    ProgrammableDiscontinuityEvent(
                        kind="actuate",
                        cell=(row, col),
                        alpha=alpha,
                        z=z,
                    )
                )
            if bool(state.locked_mask[row, col]):
                events.append(ProgrammableDiscontinuityEvent(kind="lock", cell=(row, col)))
            if bool(state.removed_mask[row, col]):
                events.append(
                    ProgrammableDiscontinuityEvent(kind="remove_cell", cell=(row, col))
                )
    return tuple(events)


def _event_cell_candidates(
    config: LatticeConfig,
    state: LatticeState,
    event: ProgrammableDiscontinuityEvent,
    tolerance: float,
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    candidates: list[tuple[int, int]] = []
    if event.cell is not None:
        candidates.append((int(event.cell[0]), int(event.cell[1])))
    candidates.extend((int(row), int(col)) for row, col in event.cells)
    if event.kind == "clear_actuation" and event.cell is None and not candidates:
        candidates.extend(_state_command_support(state, tolerance))
    valid: list[tuple[int, int]] = []
    invalid: list[tuple[int, int]] = []
    for row, col in candidates:
        if 0 <= row < config.rows and 0 <= col < config.cols:
            valid.append((row, col))
        else:
            invalid.append((row, col))
    return valid, invalid


def _event_kind_to_operator_class(kind: str) -> str:
    return {
        "actuate": "local actuation operator",
        "group_actuate": "group actuation operator",
        "lock": "constraint activation operator",
        "release": "constraint release operator",
        "clear_actuation": "command reset operator",
        "remove_cell": "graph deletion operator",
        "restore_cell": "graph restoration operator",
    }.get(kind, "unknown programmable operator")


def _event_kind_to_hardware_channel(kind: str) -> str:
    return {
        "actuate": "linear alpha/z actuator channel",
        "group_actuate": "simultaneous multi-cell actuator channel",
        "lock": "cell lock or constraint latch",
        "release": "lock release channel",
        "clear_actuation": "actuator command reset",
        "remove_cell": "removable cell / deleted topology",
        "restore_cell": "restored cell / topology insertion",
    }.get(kind, "unmapped mechanism channel")


def _mask_change_count(first: np.ndarray, second: np.ndarray) -> int:
    return int(
        np.count_nonzero(np.asarray(first, dtype=bool) != np.asarray(second, dtype=bool))
    )


def _event_state_effect(
    before: LatticeState, after: LatticeState, tolerance: float
) -> dict[str, object]:
    command_alpha_delta = float(np.max(np.abs(after.actuator_grid - before.actuator_grid)))
    command_z_delta = float(np.max(np.abs(after.z_actuator_grid - before.z_actuator_grid)))
    lock_delta = _mask_change_count(before.locked_mask, after.locked_mask)
    removed_delta = _mask_change_count(before.removed_mask, after.removed_mask)
    active = bool(
        command_alpha_delta > tolerance
        or command_z_delta > tolerance
        or lock_delta > 0
        or removed_delta > 0
    )
    return {
        "active": active,
        "commandAlphaDelta": command_alpha_delta,
        "commandZDelta": command_z_delta,
        "lockDeltaCount": lock_delta,
        "removedDeltaCount": removed_delta,
    }


def _event_record(
    event: ProgrammableDiscontinuityEvent,
    index: int,
    before: LatticeState,
    after: LatticeState,
    support_cells: list[tuple[int, int]],
    invalid_cells: list[tuple[int, int]],
    error: str | None,
    tolerance: float,
) -> dict[str, object]:
    effect = _event_state_effect(before, after, tolerance)
    removed_support = [
        cell for cell in support_cells if bool(before.removed_mask[cell[0], cell[1]])
    ]
    evidence_gaps: list[str] = []
    if not support_cells:
        evidence_gaps.append("support")
    if invalid_cells:
        evidence_gaps.append("validSupport")
    if error:
        evidence_gaps.append("eventApplication")
    if not effect["active"]:
        evidence_gaps.append("stateEffect")
    if event.kind not in {
        "actuate",
        "group_actuate",
        "lock",
        "release",
        "clear_actuation",
        "remove_cell",
        "restore_cell",
    }:
        evidence_gaps.append("operatorKind")
    return {
        "index": index,
        "kind": event.kind,
        "operatorClass": _event_kind_to_operator_class(event.kind),
        "hardwareChannel": _event_kind_to_hardware_channel(event.kind),
        "cell": None if event.cell is None else _cell_dict(event.cell),
        "supportCells": [_cell_dict(cell) for cell in support_cells],
        "invalidSupportCells": [_cell_dict(cell) for cell in invalid_cells],
        "removedSupportCells": [_cell_dict(cell) for cell in removed_support],
        "alpha": float(event.alpha),
        "z": float(event.z),
        "stateEffect": effect,
        "realized": bool(effect["active"] and not error and not invalid_cells),
        "evidenceGaps": evidence_gaps,
        "claimLabel": (
            "simulator-derived realization map; hardware channel remains "
            "experimentally unvalidated"
        ),
        "error": error,
    }


def physical_realization_map_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    contact_graph_report: dict[str, object] | None = None,
    contact_report: dict[str, object] | None = None,
    contact_stiffness: float = 1.0,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Map abstract discontinuity operators to simulator state-change evidence."""

    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    sequence = (
        tuple(event_sequence)
        if event_sequence is not None
        else _infer_realization_events(config, base, tolerance)
    )
    current = base
    records: list[dict[str, object]] = []
    aggregate_support: list[tuple[int, int]] = []
    for index, event in enumerate(sequence):
        support, invalid = _event_cell_candidates(config, current, event, tolerance)
        aggregate_support.extend(support)
        error: str | None = None
        try:
            after = apply_programmable_event(config, current, event)
        except Exception as exc:  # pragma: no cover - defensive report path
            error = str(exc)
            after = current
        records.append(
            _event_record(
                event,
                index,
                current,
                after,
                support,
                invalid,
                error,
                tolerance,
            )
        )
        current = after

    unique_support = tuple(dict.fromkeys(aggregate_support))
    resolved_contact_graph = contact_graph_report or contact_graph_consistency_report(
        config,
        current,
        contact_report=contact_report,
        group_support=unique_support,
        contact_stiffness=contact_stiffness,
        tolerance=tolerance,
    )
    contact_graph_summary = (
        resolved_contact_graph.get("summary", {})
        if isinstance(resolved_contact_graph, dict)
        else {}
    )
    realized_count = sum(1 for record in records if record["realized"])
    support_record_count = sum(1 for record in records if record["supportCells"])
    state_effect_count = sum(
        1 for record in records if record["stateEffect"].get("active")
    )
    claim_label_count = sum(1 for record in records if record.get("claimLabel"))
    missing_evidence: list[str] = []
    if not records:
        missing_evidence.append("abstractOperators")
    if realized_count < len(records):
        missing_evidence.append("realizationRecords")
    if support_record_count < len(records):
        missing_evidence.append("supportRecords")
    if state_effect_count < len(records):
        missing_evidence.append("stateEffectRecords")
    if claim_label_count < len(records):
        missing_evidence.append("claimLabels")
    if not contact_graph_summary.get("contactGraphConsistent"):
        missing_evidence.append("contactGraphConsistency")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.physical-realization-map.v1",
        "method": (
            "finite map from abstract programmable discontinuity events to "
            "simulator state effects, support records, mechanism channels, "
            "and contact-graph consistency evidence"
        ),
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "backlash": config.backlash,
            "pinHoleClearance": config.pin_hole_clearance,
        },
        "operators": records,
        "support": {
            "supportCells": [_cell_dict(cell) for cell in unique_support],
            "supportCellCount": len(unique_support),
        },
        "contactGraph": {
            "attached": isinstance(resolved_contact_graph, dict),
            "schema": resolved_contact_graph.get("schema", "")
            if isinstance(resolved_contact_graph, dict)
            else "",
            "status": contact_graph_summary.get("status", ""),
            "contactGraphConsistent": bool(
                contact_graph_summary.get("contactGraphConsistent")
            ),
        },
        "summary": {
            "status": (
                "physical-realization-map-ready"
                if ready
                else "needs-realization-map-evidence"
            ),
            "physicalRealizationMapReady": ready,
            "abstractOperatorCount": len(records),
            "realizedOperatorCount": realized_count,
            "supportRecordCount": support_record_count,
            "stateEffectRecordCount": state_effect_count,
            "claimLabelCount": claim_label_count,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "physical_realization_map_gate",
            "leanStructure": "Mechanics.PhysicalRealizationMapNat",
            "leanPredicate": "physicalRealizationMapReadyNat",
            "schema": "rad-sim.physical-realization-map.v1",
        },
        "claimLabels": {
            "map": "simulator-derived physical realization map",
            "operatorSemantics": "Lean-proven finite operator/event scaffold",
            "contactGraph": "Lean-proven graph/contact bookkeeping scaffold",
            "hardwareRealization": "experimentally unvalidated physical assumption",
        },
        "limitations": [
            "The map proves finite bookkeeping, not that a real actuator produces the same state change.",
            "Mechanism channels are simulator labels until hardware experiments identify them.",
            "Contact-graph consistency is necessary evidence, not physical contact validation.",
        ],
    }


def export_physical_realization_map_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    contact_graph_report: dict[str, object] | None = None,
    contact_report: dict[str, object] | None = None,
    contact_stiffness: float = 1.0,
    tolerance: float = 1e-9,
) -> str:
    return json.dumps(
        physical_realization_map_report(
            config,
            state,
            event_sequence=event_sequence,
            contact_graph_report=contact_graph_report,
            contact_report=contact_report,
            contact_stiffness=contact_stiffness,
            tolerance=tolerance,
        ),
        indent=2,
    )


def export_physical_realization_map_csv(report: dict[str, object]) -> str:
    rows = [[
        "index",
        "kind",
        "operator_class",
        "hardware_channel",
        "realized",
        "support_cells",
        "removed_support_cells",
        "command_alpha_delta",
        "command_z_delta",
        "lock_delta_count",
        "removed_delta_count",
        "evidence_gaps",
    ]]
    for record in report.get("operators", []):
        if not isinstance(record, dict):
            continue
        effect = record.get("stateEffect", {})
        if not isinstance(effect, dict):
            effect = {}
        rows.append(
            [
                record.get("index", ""),
                record.get("kind", ""),
                record.get("operatorClass", ""),
                record.get("hardwareChannel", ""),
                record.get("realized", ""),
                ";".join(
                    f"{cell.get('row')}:{cell.get('col')}"
                    for cell in record.get("supportCells", [])
                    if isinstance(cell, dict)
                ),
                ";".join(
                    f"{cell.get('row')}:{cell.get('col')}"
                    for cell in record.get("removedSupportCells", [])
                    if isinstance(cell, dict)
                ),
                effect.get("commandAlphaDelta", ""),
                effect.get("commandZDelta", ""),
                effect.get("lockDeltaCount", ""),
                effect.get("removedDeltaCount", ""),
                ";".join(str(item) for item in record.get("evidenceGaps", []))
                if isinstance(record.get("evidenceGaps"), list)
                else "",
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row) for row in rows
    )


_EXTERNAL_ENGINE_CANDIDATES: tuple[dict[str, object], ...] = (
    {
        "engine": "MuJoCo",
        "pythonPackage": "mujoco",
        "role": "rigid-body contact dynamics with constraints and gravity",
        "supportsRigidBody": True,
        "supportsContact": True,
        "supportsGravity": True,
        "supportsJointConstraints": True,
        "supportsHeadless": True,
        "independentFromRadSolver": True,
        "notes": "Good candidate for reduced-order pin/contact validation once geometry is exported.",
    },
    {
        "engine": "PyBullet",
        "pythonPackage": "pybullet",
        "role": "open rigid-body/contact validation and quick gravity fixtures",
        "supportsRigidBody": True,
        "supportsContact": True,
        "supportsGravity": True,
        "supportsJointConstraints": True,
        "supportsHeadless": True,
        "independentFromRadSolver": True,
        "notes": "Useful for accessible fixed-cell and gravity sanity checks.",
    },
    {
        "engine": "Project Chrono",
        "pythonPackage": "pychrono",
        "role": "multibody dynamics and contact/friction validation",
        "supportsRigidBody": True,
        "supportsContact": True,
        "supportsGravity": True,
        "supportsJointConstraints": True,
        "supportsHeadless": True,
        "independentFromRadSolver": True,
        "notes": "Candidate for contact/friction-heavy studies if Python bindings are available.",
    },
)


def _engine_available(package_name: str, override: dict[str, bool] | None) -> bool:
    if override is not None and package_name in override:
        return bool(override[package_name])
    return importlib.util.find_spec(package_name) is not None


def external_physics_engine_audit_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    engine_availability: dict[str, bool] | None = None,
    required_features: Iterable[str] = (
        "supportsRigidBody",
        "supportsContact",
        "supportsGravity",
        "supportsJointConstraints",
        "supportsHeadless",
    ),
    contact_report: dict[str, object] | None = None,
    physical_realization_report: dict[str, object] | None = None,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Audit independent physics engines for future RAD contact validation."""

    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    resolved_contact = contact_report or contact_state_abstraction_report(
        config,
        base,
        tolerance=tolerance,
    )
    resolved_realization = physical_realization_report or physical_realization_map_report(
        config,
        base,
        contact_report=resolved_contact,
        event_sequence=event_sequence,
        tolerance=tolerance,
    )
    contact_summary = (
        resolved_contact.get("summary", {})
        if isinstance(resolved_contact, dict)
        else {}
    )
    contact_topology = (
        resolved_contact.get("topology", {})
        if isinstance(resolved_contact, dict)
        else {}
    )
    realization_summary = (
        resolved_realization.get("summary", {})
        if isinstance(resolved_realization, dict)
        else {}
    )
    features = tuple(dict.fromkeys(str(feature) for feature in required_features))
    engines: list[dict[str, object]] = []
    feasible_count = 0
    available_count = 0
    for candidate in _EXTERNAL_ENGINE_CANDIDATES:
        package_name = str(candidate["pythonPackage"])
        available = _engine_available(package_name, engine_availability)
        feature_pass = all(bool(candidate.get(feature)) for feature in features)
        feasible = bool(
            available
            and feature_pass
            and candidate.get("independentFromRadSolver")
        )
        if available:
            available_count += 1
        if feasible:
            feasible_count += 1
        engines.append(
            {
                **candidate,
                "availableInPython": available,
                "requiredFeaturePass": feature_pass,
                "feasibleForAudit": feasible,
                "missingFeatures": [
                    feature for feature in features if not bool(candidate.get(feature))
                ],
            }
        )

    active_bodies = int(contact_topology.get("activeBodyCount", 0) or 0)
    contact_records = int(contact_topology.get("contactStateRecordCount", 0) or 0)
    penalty_terms = int(contact_topology.get("penaltyTermCount", 0) or 0)
    abstract_operators = int(
        realization_summary.get("abstractOperatorCount", 0) or 0
    )
    scenario_records = int(active_bodies > 0) + int(abstract_operators > 0)
    missing_evidence: list[str] = []
    if not engines:
        missing_evidence.append("engineCandidates")
    if available_count <= 0:
        missing_evidence.append("availableExternalEngine")
    if not features:
        missing_evidence.append("requiredFeatureRecords")
    if scenario_records <= 0:
        missing_evidence.append("scenarioRecords")
    if not contact_summary.get("contactStateAbstractionReady"):
        missing_evidence.append("contactModelRecords")
    if contact_records <= 0 or penalty_terms <= 0:
        missing_evidence.append("contactPenaltyRecords")
    if feasible_count <= 0:
        missing_evidence.append("independentToolRecords")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.external-physics-engine-audit.v1",
        "method": (
            "finite feasibility audit for independent rigid-body/contact engines "
            "before treating RAD physical previews as externally validated"
        ),
        "requiredFeatures": list(features),
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "pinHoleClearance": config.pin_hole_clearance,
            "removedCells": int(np.count_nonzero(base.removed_mask)),
        },
        "engines": engines,
        "contactEvidence": {
            "schema": resolved_contact.get("schema", "")
            if isinstance(resolved_contact, dict)
            else "",
            "contactStateAbstractionReady": bool(
                contact_summary.get("contactStateAbstractionReady")
            ),
            "activeBodyCount": active_bodies,
            "contactStateRecordCount": contact_records,
            "penaltyTermCount": penalty_terms,
        },
        "realizationEvidence": {
            "schema": resolved_realization.get("schema", "")
            if isinstance(resolved_realization, dict)
            else "",
            "physicalRealizationMapReady": bool(
                realization_summary.get("physicalRealizationMapReady")
            ),
            "abstractOperatorCount": abstract_operators,
        },
        "summary": {
            "status": (
                "external-physics-engine-audit-ready"
                if ready
                else "needs-external-physics-engine"
            ),
            "externalPhysicsEngineAuditReady": ready,
            "engineCandidateCount": len(engines),
            "availableEngineCount": available_count,
            "feasibleEngineCount": feasible_count,
            "requiredFeatureCount": len(features),
            "scenarioRecordCount": scenario_records,
            "contactModelRecordCount": contact_records + penalty_terms,
            "independentToolRecordCount": feasible_count,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "external_physics_engine_audit_gate",
            "leanStructure": "Mechanics.ExternalPhysicsEngineAuditNat",
            "leanPredicate": "externalPhysicsEngineAuditReadyNat",
            "schema": "rad-sim.external-physics-engine-audit.v1",
        },
        "claimLabels": {
            "audit": "finite external-tool feasibility predicate",
            "engineAvailability": "local environment diagnostic",
            "contactScenario": "simulator-derived contact abstraction",
            "physicalAccuracy": "experimentally unvalidated until an external engine run and bench data agree",
        },
        "limitations": [
            "This audit does not run an external solver; it only checks feasibility and evidence coverage.",
            "Package availability is environment-specific unless explicitly supplied by the caller.",
            "A passing audit must be followed by exported geometry, engine-specific model construction, and bench comparison.",
        ],
    }


def export_external_physics_engine_audit_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    **kwargs: object,
) -> str:
    return json.dumps(
        external_physics_engine_audit_report(config, state, **kwargs),
        indent=2,
    )


def export_external_physics_engine_audit_csv(report: dict[str, object]) -> str:
    rows = [[
        "engine",
        "package",
        "available",
        "feasible",
        "rigid_body",
        "contact",
        "gravity",
        "joint_constraints",
        "headless",
        "missing_features",
    ]]
    for engine in report.get("engines", []):
        if not isinstance(engine, dict):
            continue
        rows.append(
            [
                engine.get("engine", ""),
                engine.get("pythonPackage", ""),
                engine.get("availableInPython", ""),
                engine.get("feasibleForAudit", ""),
                engine.get("supportsRigidBody", ""),
                engine.get("supportsContact", ""),
                engine.get("supportsGravity", ""),
                engine.get("supportsJointConstraints", ""),
                engine.get("supportsHeadless", ""),
                ";".join(str(item) for item in engine.get("missingFeatures", []))
                if isinstance(engine.get("missingFeatures"), list)
                else "",
            ]
        )
    summary = report.get("summary", {})
    if isinstance(summary, dict):
        rows.append(
            [
                "summary",
                report.get("schema", ""),
                summary.get("externalPhysicsEngineAuditReady", ""),
                summary.get("feasibleEngineCount", ""),
                summary.get("availableEngineCount", ""),
                summary.get("requiredFeatureCount", ""),
                summary.get("scenarioRecordCount", ""),
                "",
                "",
                ";".join(str(item) for item in summary.get("missingEvidence", []))
                if isinstance(summary.get("missingEvidence"), list)
                else "",
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row)
        for row in rows
    )


def _mujoco_cell_name(row: int, col: int) -> str:
    return f"rad_cell_{int(row)}_{int(col)}"


def _vec3(values: Iterable[float]) -> tuple[float, float, float]:
    items = tuple(float(value) for value in values)
    if len(items) < 3:
        return tuple(items + (0.0,) * (3 - len(items)))  # type: ignore[return-value]
    return (items[0], items[1], items[2])


def _mujoco_vec(values: Iterable[float]) -> str:
    return " ".join(f"{float(value):.9g}" for value in values)


def _mujoco_load_case(load_case: LoadCase | None) -> LoadCase:
    return load_case if load_case is not None else LoadCase()


def _tuple_with_default(
    values: Iterable[float] | None,
    default: tuple[float, ...],
) -> tuple[float, ...]:
    if values is None:
        return default
    items = tuple(float(value) for value in values)
    if len(items) >= len(default):
        return items[: len(default)]
    return items + default[len(items) :]


def _mujoco_contact_parameter_defaults(
    *,
    contact_stiffness: float = 1000.0,
    contact_damping: float = 2.0,
    friction: Iterable[float] | None = None,
    solref: Iterable[float] | None = None,
    solimp: Iterable[float] | None = None,
    contact_margin: float = 0.0,
    contact_gap: float = 0.0,
    condim: int = 3,
    calibrated: bool = False,
) -> dict[str, object]:
    return {
        "contactStiffness": max(0.0, float(contact_stiffness)),
        "contactDamping": max(0.0, float(contact_damping)),
        "friction": list(_tuple_with_default(friction, (0.4, 0.02, 0.001))),
        "solref": list(_tuple_with_default(solref, (0.02, 1.0))),
        "solimp": list(_tuple_with_default(solimp, (0.9, 0.95, 0.001, 0.5, 2.0))),
        "margin": max(0.0, float(contact_margin)),
        "gap": max(0.0, float(contact_gap)),
        "condim": max(1, int(condim)),
        "calibrated": bool(calibrated),
    }


def _mujoco_contact_parameter_xml_attributes(defaults: dict[str, object]) -> str:
    friction = _mujoco_vec(defaults.get("friction", (0.4, 0.02, 0.001)))  # type: ignore[arg-type]
    solref = _mujoco_vec(defaults.get("solref", (0.02, 1.0)))  # type: ignore[arg-type]
    solimp = _mujoco_vec(defaults.get("solimp", (0.9, 0.95, 0.001, 0.5, 2.0)))  # type: ignore[arg-type]
    return (
        f' condim="{int(defaults["condim"])}"'
        f' friction="{friction}"'
        f' solref="{solref}"'
        f' solimp="{solimp}"'
        f' margin="{float(defaults["margin"]):.9g}"'
        f' gap="{float(defaults["gap"]):.9g}"'
    )


def _mujoco_contact_parameter_profile_from_pairs(
    contact_pairs: Iterable[dict[str, object]],
    *,
    defaults: dict[str, object],
) -> dict[str, object]:
    pairs = [pair for pair in contact_pairs if isinstance(pair, dict)]
    xml_attributes = _mujoco_contact_parameter_xml_attributes(defaults)
    parameters: list[dict[str, object]] = []
    for pair in pairs:
        parameters.append(
            {
                "row": int(pair.get("row", 0)),
                "col": int(pair.get("col", 0)),
                "pivot": str(pair.get("pivot", "")),
                "pin": str(pair.get("pin", "")),
                "hole": str(pair.get("hole", "")),
                "contactStiffness": float(defaults["contactStiffness"]),
                "contactDamping": float(defaults["contactDamping"]),
                "friction": list(defaults["friction"]),  # type: ignore[arg-type]
                "solref": list(defaults["solref"]),  # type: ignore[arg-type]
                "solimp": list(defaults["solimp"]),  # type: ignore[arg-type]
                "margin": float(defaults["margin"]),
                "gap": float(defaults["gap"]),
                "condim": int(defaults["condim"]),
                "calibrated": bool(defaults["calibrated"]),
                "xmlAttributes": xml_attributes,
            }
        )
    contact_pair_count = len(pairs)
    parameter_count = len(parameters)
    friction_records = sum(1 for item in parameters if item["friction"])
    solver_records = sum(1 for item in parameters if item["solref"] and item["solimp"])
    stiffness_records = sum(1 for item in parameters if item["contactStiffness"] > 0)
    damping_records = sum(1 for item in parameters if item["contactDamping"] > 0)
    calibrated_records = sum(1 for item in parameters if item["calibrated"])
    missing_evidence: list[str] = []
    if contact_pair_count <= 0:
        missing_evidence.append("contactPairRecords")
    if parameter_count < contact_pair_count:
        missing_evidence.append("parameterRecords")
    if friction_records < contact_pair_count:
        missing_evidence.append("frictionRecords")
    if solver_records < contact_pair_count:
        missing_evidence.append("solverParameterRecords")
    if stiffness_records < contact_pair_count:
        missing_evidence.append("stiffnessRecords")
    if damping_records < contact_pair_count:
        missing_evidence.append("dampingRecords")
    if not xml_attributes:
        missing_evidence.append("xmlContactAttributes")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.mujoco-contact-parameter-profile.v1",
        "method": (
            "finite contact-parameter profile for MuJoCo proxy contact records; "
            "parameters are exported but remain uncalibrated unless calibrated=true"
        ),
        "defaults": defaults,
        "parameters": parameters,
        "xmlAttributes": xml_attributes,
        "summary": {
            "status": (
                "mujoco-contact-parameter-profile-ready"
                if ready
                else "needs-mujoco-contact-parameters"
            ),
            "mujocoContactParameterProfileReady": ready,
            "contactPairRecordCount": contact_pair_count,
            "parameterRecordCount": parameter_count,
            "frictionRecordCount": friction_records,
            "solverParameterRecordCount": solver_records,
            "stiffnessRecordCount": stiffness_records,
            "dampingRecordCount": damping_records,
            "calibratedRecordCount": calibrated_records,
            "xmlAttributeByteCount": len(xml_attributes.encode("utf-8")),
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "mujoco_contact_parameter_profile_gate",
            "leanStructure": "Mechanics.ExternalContactParameterProfileNat",
            "leanPredicate": "externalContactParameterProfileReadyNat",
            "schema": "rad-sim.mujoco-contact-parameter-profile.v1",
        },
        "claimLabels": {
            "parameters": "simulator-derived external contact parameter profile",
            "calibration": (
                "calibrated contact parameter"
                if bool(defaults["calibrated"])
                else "experimentally unvalidated contact parameter assumption"
            ),
            "externalEngine": "MuJoCo XML contact attributes for proxy geometry",
        },
        "limitations": [
            "Contact stiffness and damping are proxy metadata, not derived from MuJoCo constitutive laws.",
            "Friction, solref, and solimp values require bench calibration before physical claims.",
            "A complete parameter profile does not prove external-run or hardware agreement.",
        ],
    }


def _mujoco_contact_geometry_from_bodies(
    config: LatticeConfig,
    bodies: Iterable[dict[str, object]],
    *,
    plate_thickness: float,
    contact_stiffness: float = 1000.0,
    contact_damping: float = 2.0,
    friction: Iterable[float] | None = None,
    solref: Iterable[float] | None = None,
    solimp: Iterable[float] | None = None,
    contact_margin: float = 0.0,
    contact_gap: float = 0.0,
    condim: int = 3,
    calibrated_contact: bool = False,
) -> dict[str, object]:
    contact_defaults = _mujoco_contact_parameter_defaults(
        contact_stiffness=contact_stiffness,
        contact_damping=contact_damping,
        friction=friction,
        solref=solref,
        solimp=solimp,
        contact_margin=contact_margin,
        contact_gap=contact_gap,
        condim=condim,
        calibrated=calibrated_contact,
    )
    contact_attributes = _mujoco_contact_parameter_xml_attributes(contact_defaults)
    pin_radius = max(float(config.pin_radius), 1e-9)
    hole_radius = max(float(config.hole_radius), pin_radius)
    clearance = hole_radius - pin_radius
    half_height = max(float(plate_thickness) * 0.55, 1e-9)
    pin_records: list[dict[str, object]] = []
    hole_records: list[dict[str, object]] = []
    clearance_records: list[dict[str, object]] = []
    contact_pairs: list[dict[str, object]] = []
    body_geom_xml: dict[str, list[str]] = {}
    for body in bodies:
        name = str(body.get("name", ""))
        if not name:
            continue
        row = int(body.get("row", 0))
        col = int(body.get("col", 0))
        center = _vec3(body.get("position", (0.0, 0.0, 0.0)))  # type: ignore[arg-type]
        alpha = max(float(body.get("alpha", config.initial_alpha)), 0.1)
        pivot_span = max(config.cell_size * 0.18, 0.32 * config.cell_size * alpha)
        pivots = (
            ("nw", (-pivot_span, -pivot_span, 0.0)),
            ("ne", (pivot_span, -pivot_span, 0.0)),
            ("se", (pivot_span, pivot_span, 0.0)),
            ("sw", (-pivot_span, pivot_span, 0.0)),
        )
        body_geom_xml[name] = []
        for pivot_name, local in pivots:
            pin_name = f"{name}_{pivot_name}_pin"
            hole_name = f"{name}_{pivot_name}_hole"
            world = (
                center[0] + local[0],
                center[1] + local[1],
                center[2] + local[2],
            )
            pin_records.append(
                {
                    "row": row,
                    "col": col,
                    "pivot": pivot_name,
                    "name": pin_name,
                    "radius": pin_radius,
                    "halfHeight": half_height,
                    "position": list(world),
                }
            )
            hole_records.append(
                {
                    "row": row,
                    "col": col,
                    "pivot": pivot_name,
                    "name": hole_name,
                    "radius": hole_radius,
                    "halfHeight": half_height,
                    "position": list(world),
                }
            )
            clearance_records.append(
                {
                    "row": row,
                    "col": col,
                    "pivot": pivot_name,
                    "pinRadius": pin_radius,
                    "holeRadius": hole_radius,
                    "clearance": clearance,
                    "clearanceRatio": clearance / hole_radius if hole_radius > 0 else 0.0,
                }
            )
            contact_pairs.append(
                {
                    "row": row,
                    "col": col,
                    "pivot": pivot_name,
                    "pin": pin_name,
                    "hole": hole_name,
                    "clearance": clearance,
                    "active": clearance >= 0.0,
                }
            )
            body_geom_xml[name].append(
                f'      <geom name="{pin_name}" type="cylinder" pos="{_mujoco_vec(local)}" size="{pin_radius:.9g} {half_height:.9g}"{contact_attributes} rgba="0.05 0.05 0.05 1"/>'
            )
            body_geom_xml[name].append(
                f'      <geom name="{hole_name}_clearance" type="cylinder" pos="{_mujoco_vec(local)}" size="{hole_radius:.9g} {half_height * 1.05:.9g}" contype="0" conaffinity="0" rgba="0.1 0.7 0.9 0.18"/>'
            )
    xml_fragment = "\n".join(
        line
        for lines in body_geom_xml.values()
        for line in lines
    )
    missing_evidence: list[str] = []
    if not pin_records:
        missing_evidence.append("pinRecords")
    if not hole_records:
        missing_evidence.append("holeRecords")
    if not clearance_records:
        missing_evidence.append("clearanceRecords")
    if not contact_pairs:
        missing_evidence.append("contactPairRecords")
    if clearance < 0:
        missing_evidence.append("nonnegativeClearance")
    if not xml_fragment:
        missing_evidence.append("mjcfContactGeometry")
    contact_parameter_profile = _mujoco_contact_parameter_profile_from_pairs(
        contact_pairs,
        defaults=contact_defaults,
    )
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.mujoco-pin-hole-contact-geometry.v1",
        "method": (
            "finite pin-hole contact-geometry inventory for the MuJoCo RAD "
            "handoff; records geometry proxies but does not calibrate contact"
        ),
        "pins": pin_records,
        "holes": hole_records,
        "clearanceRecords": clearance_records,
        "contactPairs": contact_pairs,
        "contactParameterProfile": contact_parameter_profile,
        "bodyGeomXml": body_geom_xml,
        "xmlFragment": xml_fragment,
        "summary": {
            "status": (
                "mujoco-pin-hole-contact-geometry-ready"
                if ready
                else "needs-mujoco-pin-hole-contact-geometry"
            ),
            "mujocoPinHoleContactGeometryReady": ready,
            "pinRecordCount": len(pin_records),
            "holeRecordCount": len(hole_records),
            "clearanceRecordCount": len(clearance_records),
            "contactPairRecordCount": len(contact_pairs),
            "activeContactPairCount": sum(1 for pair in contact_pairs if pair["active"]),
            "contactParameterRecordCount": contact_parameter_profile["summary"]["parameterRecordCount"],
            "frictionRecordCount": contact_parameter_profile["summary"]["frictionRecordCount"],
            "solverParameterRecordCount": contact_parameter_profile["summary"]["solverParameterRecordCount"],
            "minClearance": min((record["clearance"] for record in clearance_records), default=0.0),
            "maxClearance": max((record["clearance"] for record in clearance_records), default=0.0),
            "xmlFragmentByteCount": len(xml_fragment.encode("utf-8")),
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "mujoco_pin_hole_contact_geometry_gate",
            "leanStructure": "Mechanics.ExternalContactGeometryNat",
            "leanPredicate": "externalContactGeometryReadyNat",
            "schema": "rad-sim.mujoco-pin-hole-contact-geometry.v1",
        },
        "claimLabels": {
            "geometry": "simulator-derived pin-hole proxy geometry",
            "clearance": "configured geometric clearance record",
            "contact": "external-engine contact candidate; stiffness/friction remain uncalibrated",
        },
        "limitations": [
            "Pin and clearance cylinders are proxy geometry, not measured CAD surfaces.",
            "The hole clearance shell is non-colliding in the exported MJCF until contact parameters are calibrated.",
            "Friction, compliance, assembly offsets, and wear are not represented.",
        ],
    }


def mujoco_model_export_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    load_case: LoadCase | None = None,
    gravity: tuple[float, float, float] = (0.0, 0.0, -9.81),
    timestep: float = 0.002,
    plate_thickness: float = 0.04,
    body_density: float = 900.0,
    contact_stiffness: float = 1000.0,
    contact_damping: float = 2.0,
    friction: Iterable[float] | None = None,
    solref: Iterable[float] | None = None,
    solimp: Iterable[float] | None = None,
    contact_margin: float = 0.0,
    contact_gap: float = 0.0,
    condim: int = 3,
    calibrated_contact: bool = False,
) -> dict[str, object]:
    """Export a coarse RAD sheet to MJCF for independent MuJoCo validation."""

    if timestep <= 0:
        raise ValueError("timestep must be positive")
    if plate_thickness <= 0:
        raise ValueError("plate_thickness must be positive")
    if body_density <= 0:
        raise ValueError("body_density must be positive")
    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    sequence = tuple(event_sequence or ())
    export_state = _apply_realization_event_sequence(config, base, sequence)
    sim = simulate_kinematic(config, export_state)
    load = _mujoco_load_case(load_case)
    fixed_cells = {tuple(cell) for cell in load.fixed_cells}
    force_cells = {tuple(cell) for cell in load.external_forces}
    bodies: list[dict[str, object]] = []
    removed_records: list[dict[str, int]] = []
    half_thickness = plate_thickness / 2.0
    min_half = max(config.cell_size * 0.05, 1e-6)
    for row in range(config.rows):
        for col in range(config.cols):
            if bool(export_state.removed_mask[row, col]):
                removed_records.append({"row": row, "col": col})
                continue
            alpha = float(sim.alpha[row, col])
            center = tuple(float(value) for value in sim.deformed_centers_3d[row, col])
            half_side = max(min_half, 0.24 * config.cell_size * max(alpha, 0.1))
            cell = (row, col)
            name = _mujoco_cell_name(row, col)
            fixed = cell in fixed_cells
            body = {
                "row": row,
                "col": col,
                "name": name,
                "position": list(center),
                "halfExtents": [half_side, half_side, half_thickness],
                "alpha": alpha,
                "thetaDegrees": float(sim.theta_degrees[row, col]),
                "commandAlpha": float(export_state.actuator_grid[row, col]),
                "commandZ": float(export_state.z_actuator_grid[row, col]),
                "locked": bool(export_state.locked_mask[row, col]),
                "fixed": fixed,
                "hasExternalForce": cell in force_cells,
            }
            bodies.append(body)

    contact_geometry = _mujoco_contact_geometry_from_bodies(
        config,
        bodies,
        plate_thickness=plate_thickness,
        contact_stiffness=contact_stiffness,
        contact_damping=contact_damping,
        friction=friction,
        solref=solref,
        solimp=solimp,
        contact_margin=contact_margin,
        contact_gap=contact_gap,
        condim=condim,
        calibrated_contact=calibrated_contact,
    )
    body_geom_xml = contact_geometry.get("bodyGeomXml", {})
    if not isinstance(body_geom_xml, dict):
        body_geom_xml = {}
    size = max(config.rows, config.cols, 1) * config.cell_size * 2.0
    xml_lines = [
        '<mujoco model="rad_lattice_external_validation">',
        '  <compiler angle="radian"/>',
        f'  <option timestep="{timestep:.9g}" gravity="{_mujoco_vec(gravity)}"/>',
        "  <worldbody>",
        f'    <geom name="floor" type="plane" size="{size:.9g} {size:.9g} 0.05" rgba="0.65 0.65 0.65 0.25"/>',
    ]
    for body in bodies:
        name = str(body["name"])
        position = _mujoco_vec(body["position"])  # type: ignore[arg-type]
        half_extents = _mujoco_vec(body["halfExtents"])  # type: ignore[arg-type]
        rgba = "0.20 0.45 0.85 1" if body["fixed"] else "0.85 0.55 0.18 1"
        xml_lines.append(f'    <body name="{name}" pos="{position}">')
        if not body["fixed"]:
            xml_lines.append(f'      <joint name="{name}_free" type="free" damping="0.05"/>')
        xml_lines.append(
            f'      <geom name="{name}_plate" type="box" size="{half_extents}" density="{body_density:.9g}" rgba="{rgba}"/>'
        )
        for line in body_geom_xml.get(name, []):
            xml_lines.append(str(line))
        xml_lines.append("    </body>")
    xml_lines.extend(["  </worldbody>", "</mujoco>"])
    xml = "\n".join(xml_lines)
    load_records = [
        {
            "row": int(cell[0]),
            "col": int(cell[1]),
            "force": list(_vec3(force)),
        }
        for cell, force in load.external_forces.items()
        if 0 <= int(cell[0]) < config.rows
        and 0 <= int(cell[1]) < config.cols
        and not bool(export_state.removed_mask[int(cell[0]), int(cell[1])])
    ]
    fixed_body_count = sum(1 for body in bodies if bool(body["fixed"]))
    gravity_enabled = any(abs(float(value)) > 0 for value in gravity)
    missing_evidence: list[str] = []
    if not bodies:
        missing_evidence.append("bodyRecords")
    if fixed_body_count <= 0:
        missing_evidence.append("fixedBodyRecords")
    if not gravity_enabled:
        missing_evidence.append("gravityRecord")
    if not xml.strip():
        missing_evidence.append("mjcfXml")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.mujoco-model-export.v1",
        "method": "coarse MJCF export for independent rigid-body/contact validation; not fabrication-accurate CAD",
        "engine": "MuJoCo",
        "xml": xml,
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "cellSize": config.cell_size,
            "backlash": config.backlash,
            "pinHoleClearance": config.pin_hole_clearance,
        },
        "settings": {
            "gravity": list(gravity),
            "timestep": timestep,
            "plateThickness": plate_thickness,
            "bodyDensity": body_density,
            "contactParameters": contact_geometry["contactParameterProfile"]["defaults"],
        },
        "bodies": bodies,
        "removedCells": removed_records,
        "contactGeometry": contact_geometry,
        "loads": load_records,
        "summary": {
            "status": "mujoco-model-export-ready" if ready else "needs-mujoco-export-evidence",
            "mujocoModelExportReady": ready,
            "bodyRecordCount": len(bodies),
            "fixedBodyCount": fixed_body_count,
            "removedBodyCount": len(removed_records),
            "loadRecordCount": len(load_records),
            "pinRecordCount": contact_geometry["summary"]["pinRecordCount"],
            "holeRecordCount": contact_geometry["summary"]["holeRecordCount"],
            "clearanceRecordCount": contact_geometry["summary"]["clearanceRecordCount"],
            "contactPairRecordCount": contact_geometry["summary"]["contactPairRecordCount"],
            "contactParameterRecordCount": contact_geometry["summary"]["contactParameterRecordCount"],
            "frictionRecordCount": contact_geometry["summary"]["frictionRecordCount"],
            "solverParameterRecordCount": contact_geometry["summary"]["solverParameterRecordCount"],
            "gravityRecordCount": int(gravity_enabled),
            "xmlByteCount": len(xml.encode("utf-8")),
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "mujoco_model_export_gate",
            "leanStructure": "Mechanics.ExternalPhysicsModelExportNat",
            "leanPredicate": "externalPhysicsModelExportReadyNat",
            "schema": "rad-sim.mujoco-model-export.v1",
        },
        "claimLabels": {
            "export": "simulator-derived coarse external-engine model export",
            "geometry": "normalized proxy geometry, not fabrication-accurate CAD",
            "physics": "experimentally unvalidated until MuJoCo run and bench comparison pass",
        },
        "limitations": [
            "Each cell is represented as one coarse box body rather than the full RAD pin/plate assembly.",
            "Fixed cells are implemented by omitting free joints; this is an external validation fixture, not a hardware mount design.",
            "Pin-hole clearance is recorded in metadata but not yet represented as exact contact geometry in MJCF.",
        ],
    }


def export_mujoco_model_xml(
    config: LatticeConfig,
    state: LatticeState | None = None,
    **kwargs: object,
) -> str:
    return str(mujoco_model_export_report(config, state, **kwargs)["xml"])


def export_mujoco_model_report_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    **kwargs: object,
) -> str:
    return json.dumps(mujoco_model_export_report(config, state, **kwargs), indent=2)


def mujoco_pin_hole_contact_geometry_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    export_report: dict[str, object] | None = None,
    plate_thickness: float = 0.04,
) -> dict[str, object]:
    """Return the pin, hole, clearance, and contact-pair records for MuJoCo export."""

    if export_report is not None:
        contact_geometry = export_report.get("contactGeometry")
        if isinstance(contact_geometry, dict):
            return contact_geometry
        bodies = export_report.get("bodies", [])
        if isinstance(bodies, list):
            return _mujoco_contact_geometry_from_bodies(
                config,
                (body for body in bodies if isinstance(body, dict)),
                plate_thickness=plate_thickness,
            )
    report = mujoco_model_export_report(
        config,
        state,
        event_sequence=event_sequence,
        plate_thickness=plate_thickness,
    )
    contact_geometry = report.get("contactGeometry", {})
    return contact_geometry if isinstance(contact_geometry, dict) else {}


def export_mujoco_pin_hole_contact_geometry_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    **kwargs: object,
) -> str:
    return json.dumps(
        mujoco_pin_hole_contact_geometry_report(config, state, **kwargs),
        indent=2,
    )


def export_mujoco_pin_hole_contact_geometry_csv(report: dict[str, object]) -> str:
    rows = [[
        "row",
        "col",
        "pivot",
        "pin",
        "hole",
        "pin_radius",
        "hole_radius",
        "clearance",
        "active",
    ]]
    pair_by_key = {
        (
            int(pair.get("row", -1)),
            int(pair.get("col", -1)),
            str(pair.get("pivot", "")),
        ): pair
        for pair in report.get("contactPairs", [])
        if isinstance(pair, dict)
    }
    for record in report.get("clearanceRecords", []):
        if not isinstance(record, dict):
            continue
        key = (
            int(record.get("row", -1)),
            int(record.get("col", -1)),
            str(record.get("pivot", "")),
        )
        pair = pair_by_key.get(key, {})
        rows.append(
            [
                record.get("row", ""),
                record.get("col", ""),
                record.get("pivot", ""),
                pair.get("pin", ""),
                pair.get("hole", ""),
                record.get("pinRadius", ""),
                record.get("holeRadius", ""),
                record.get("clearance", ""),
                pair.get("active", ""),
            ]
        )
    summary = report.get("summary", {})
    if isinstance(summary, dict):
        rows.append(
            [
                "summary",
                "",
                "",
                summary.get("pinRecordCount", ""),
                summary.get("holeRecordCount", ""),
                "",
                "",
                summary.get("minClearance", ""),
                summary.get("mujocoPinHoleContactGeometryReady", ""),
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row)
        for row in rows
    )


def mujoco_contact_parameter_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    export_report: dict[str, object] | None = None,
    contact_geometry_report: dict[str, object] | None = None,
    contact_stiffness: float = 1000.0,
    contact_damping: float = 2.0,
    friction: Iterable[float] | None = None,
    solref: Iterable[float] | None = None,
    solimp: Iterable[float] | None = None,
    contact_margin: float = 0.0,
    contact_gap: float = 0.0,
    condim: int = 3,
    calibrated_contact: bool = False,
) -> dict[str, object]:
    """Return the MuJoCo contact-parameter profile for pin-hole proxy geometry."""

    if contact_geometry_report is None:
        contact_geometry_report = mujoco_pin_hole_contact_geometry_report(
            config,
            state,
            event_sequence=event_sequence,
            export_report=export_report,
            plate_thickness=float(
                (export_report or {}).get("settings", {}).get("plateThickness", 0.04)
            )
            if isinstance((export_report or {}).get("settings", {}), dict)
            else 0.04,
        )
    if (
        friction is None
        and solref is None
        and solimp is None
        and contact_stiffness == 1000.0
        and contact_damping == 2.0
        and contact_margin == 0.0
        and contact_gap == 0.0
        and condim == 3
        and not calibrated_contact
    ):
        existing = contact_geometry_report.get("contactParameterProfile")
        if isinstance(existing, dict):
            return existing
    defaults = _mujoco_contact_parameter_defaults(
        contact_stiffness=contact_stiffness,
        contact_damping=contact_damping,
        friction=friction,
        solref=solref,
        solimp=solimp,
        contact_margin=contact_margin,
        contact_gap=contact_gap,
        condim=condim,
        calibrated=calibrated_contact,
    )
    return _mujoco_contact_parameter_profile_from_pairs(
        contact_geometry_report.get("contactPairs", []),
        defaults=defaults,
    )


def export_mujoco_contact_parameter_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    **kwargs: object,
) -> str:
    return json.dumps(
        mujoco_contact_parameter_report(config, state, **kwargs),
        indent=2,
    )


def export_mujoco_contact_parameter_csv(report: dict[str, object]) -> str:
    rows = [[
        "row",
        "col",
        "pivot",
        "pin",
        "hole",
        "contact_stiffness",
        "contact_damping",
        "friction",
        "solref",
        "solimp",
        "condim",
        "calibrated",
    ]]
    for parameter in report.get("parameters", []):
        if not isinstance(parameter, dict):
            continue
        rows.append(
            [
                parameter.get("row", ""),
                parameter.get("col", ""),
                parameter.get("pivot", ""),
                parameter.get("pin", ""),
                parameter.get("hole", ""),
                parameter.get("contactStiffness", ""),
                parameter.get("contactDamping", ""),
                ";".join(str(value) for value in parameter.get("friction", []))
                if isinstance(parameter.get("friction"), list)
                else "",
                ";".join(str(value) for value in parameter.get("solref", []))
                if isinstance(parameter.get("solref"), list)
                else "",
                ";".join(str(value) for value in parameter.get("solimp", []))
                if isinstance(parameter.get("solimp"), list)
                else "",
                parameter.get("condim", ""),
                parameter.get("calibrated", ""),
            ]
        )
    summary = report.get("summary", {})
    if isinstance(summary, dict):
        rows.append(
            [
                "summary",
                "",
                "",
                "",
                "",
                summary.get("stiffnessRecordCount", ""),
                summary.get("dampingRecordCount", ""),
                summary.get("frictionRecordCount", ""),
                summary.get("solverParameterRecordCount", ""),
                "",
                "",
                summary.get("calibratedRecordCount", ""),
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row)
        for row in rows
    )


def _contact_parameter_measurement_columns() -> list[str]:
    return [
        "dataset_id",
        "dataset_role",
        "scenario_id",
        "repeat_index",
        "row",
        "col",
        "pivot",
        "pin_radius_mm",
        "hole_radius_mm",
        "clearance_mm",
        "normal_load_n",
        "tangential_load_n",
        "imposed_z_mm",
        "measured_pin_hole_slip_mm",
        "measured_normal_force_n",
        "measured_tangent_force_n",
        "measured_rebound_ratio",
        "measured_contact_duration_s",
        "measured_static_friction_coeff",
        "measured_dynamic_friction_coeff",
        "fitted_contact_stiffness",
        "fitted_contact_damping",
        "fitted_solref_timeconst",
        "fitted_solref_damping_ratio",
        "fitted_solimp_width",
        "fixture_notes",
    ]


def _contact_parameter_template_rows(
    parameter_report: dict[str, object],
    *,
    dataset_id: str,
    dataset_role: str,
    repeat_count: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    parameters = parameter_report.get("parameters", [])
    if not isinstance(parameters, list):
        parameters = []
    for parameter in parameters:
        if not isinstance(parameter, dict):
            continue
        row = int(parameter.get("row", 0))
        col = int(parameter.get("col", 0))
        pivot = str(parameter.get("pivot", ""))
        for repeat_index in range(1, repeat_count + 1):
            scenario_id = (
                f"{dataset_role}-r{row}c{col}-{pivot}-contact-repeat-{repeat_index}"
            )
            rows.append(
                {
                    "datasetId": dataset_id,
                    "datasetRole": dataset_role,
                    "scenarioId": scenario_id,
                    "repeatIndex": repeat_index,
                    "row": row,
                    "col": col,
                    "pivot": pivot,
                    "pin": parameter.get("pin", ""),
                    "hole": parameter.get("hole", ""),
                    "assumedContactStiffness": parameter.get("contactStiffness", None),
                    "assumedContactDamping": parameter.get("contactDamping", None),
                    "assumedFriction": parameter.get("friction", []),
                    "assumedSolref": parameter.get("solref", []),
                    "assumedSolimp": parameter.get("solimp", []),
                    "assumedCondim": parameter.get("condim", None),
                    "requiredMeasurements": {
                        "pinRadiusMm": None,
                        "holeRadiusMm": None,
                        "clearanceMm": None,
                        "normalLoadN": None,
                        "tangentialLoadN": None,
                        "imposedZMm": None,
                        "measuredPinHoleSlipMm": None,
                        "measuredNormalForceN": None,
                        "measuredTangentForceN": None,
                        "measuredReboundRatio": None,
                        "measuredContactDurationS": None,
                        "measuredStaticFrictionCoeff": None,
                        "measuredDynamicFrictionCoeff": None,
                        "fittedContactStiffness": None,
                        "fittedContactDamping": None,
                        "fittedSolrefTimeconst": None,
                        "fittedSolrefDampingRatio": None,
                        "fittedSolimpWidth": None,
                        "fixtureNotes": "",
                    },
                }
            )
    return rows


def contact_parameter_calibration_packet(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    export_report: dict[str, object] | None = None,
    contact_geometry_report: dict[str, object] | None = None,
    contact_parameter_report: dict[str, object] | None = None,
    fit_dataset_id: str = "contact-fit-run-001",
    holdout_dataset_id: str = "contact-holdout-run-001",
    profile_id: str = "contact-parameter-profile-v1",
    profile_frozen_at: str | None = None,
    repeat_count: int = 3,
) -> dict[str, object]:
    """Bundle contact-parameter assumptions into a calibration measurement packet."""

    if repeat_count <= 0:
        raise ValueError("repeat_count must be positive")
    frozen_at = profile_frozen_at or "set-after-contact-fit-before-holdout"
    geometry = contact_geometry_report or mujoco_pin_hole_contact_geometry_report(
        config,
        state,
        event_sequence=event_sequence,
        export_report=export_report,
    )
    parameters = contact_parameter_report or mujoco_contact_parameter_report(
        config,
        state,
        event_sequence=event_sequence,
        export_report=export_report,
        contact_geometry_report=geometry,
    )
    summary = parameters.get("summary", {})
    if not isinstance(summary, dict):
        summary = {}
    contact_pairs = int(summary.get("contactPairRecordCount", 0) or 0)
    parameter_records = int(summary.get("parameterRecordCount", 0) or 0)
    profile_ready = bool(summary.get("mujocoContactParameterProfileReady", False))
    measurement_columns = _contact_parameter_measurement_columns()
    fit_rows = _contact_parameter_template_rows(
        parameters,
        dataset_id=fit_dataset_id,
        dataset_role="fit",
        repeat_count=repeat_count,
    )
    holdout_rows = _contact_parameter_template_rows(
        parameters,
        dataset_id=holdout_dataset_id,
        dataset_role="holdout",
        repeat_count=repeat_count,
    )
    filenames = {
        "packet": "contact_parameter_calibration_packet.json",
        "csv": "contact_parameter_calibration_template.csv",
        "contact_geometry": "mujoco_pin_hole_contact_geometry.json",
        "contact_parameters": "mujoco_contact_parameter_profile.json",
        "readme": "README.md",
    }
    artifact_manifest = [
        {
            "id": "packet",
            "filename": filenames["packet"],
            "schema": "rad-sim.contact-parameter-calibration-packet.v1",
            "purpose": "single JSON bundle for contact-parameter calibration review",
        },
        {
            "id": "template-csv",
            "filename": filenames["csv"],
            "schema": "rad-sim.contact-parameter-calibration-template.v1",
            "purpose": "blank fit and holdout measurement rows for contact pairs",
        },
        {
            "id": "contact-geometry",
            "filename": filenames["contact_geometry"],
            "schema": "rad-sim.mujoco-pin-hole-contact-geometry.v1",
            "purpose": "proxy pin-hole geometry that defines the contact pair inventory",
        },
        {
            "id": "contact-parameters",
            "filename": filenames["contact_parameters"],
            "schema": "rad-sim.mujoco-contact-parameter-profile.v1",
            "purpose": "explicit proxy stiffness, damping, friction, and solver settings",
        },
        {
            "id": "readme",
            "filename": filenames["readme"],
            "schema": "text/markdown",
            "purpose": "human instructions and claim limitations",
        },
    ]
    missing_evidence: list[str] = []
    if contact_pairs <= 0:
        missing_evidence.append("contactPairRecords")
    if parameter_records < contact_pairs:
        missing_evidence.append("parameterRecords")
    if not measurement_columns:
        missing_evidence.append("measurementColumns")
    if len(fit_rows) < contact_pairs:
        missing_evidence.append("fitTemplateRows")
    if len(holdout_rows) < contact_pairs:
        missing_evidence.append("holdoutTemplateRows")
    if len(artifact_manifest) < 5:
        missing_evidence.append("artifactManifest")
    if not profile_ready:
        missing_evidence.append("contactParameterProfile")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.contact-parameter-calibration-packet.v1",
        "method": (
            "bench handoff packet linking proxy pin-hole geometry, MuJoCo contact "
            "parameter assumptions, and blank fit/holdout measurement rows"
        ),
        "schemas": {
            "template": "rad-sim.contact-parameter-calibration-template.v1",
            "contactGeometry": "rad-sim.mujoco-pin-hole-contact-geometry.v1",
            "contactParameters": "rad-sim.mujoco-contact-parameter-profile.v1",
        },
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "cellSize": config.cell_size,
            "backlash": config.backlash,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "pinHoleClearance": config.pin_hole_clearance,
        },
        "profilePlan": {
            "profileId": profile_id,
            "profileFrozenAt": frozen_at,
            "fitDatasetId": fit_dataset_id,
            "holdoutDatasetId": holdout_dataset_id,
            "repeatCount": repeat_count,
            "freezeRule": (
                "Fit contact parameters from fit rows, freeze the profile ID, "
                "then collect holdout rows without modifying the fitted parameters."
            ),
        },
        "measurementColumns": measurement_columns,
        "fitTemplateRows": fit_rows,
        "holdoutTemplateRows": holdout_rows,
        "contactGeometry": geometry,
        "contactParameterProfile": parameters,
        "filenames": filenames,
        "artifactManifest": artifact_manifest,
        "validationInstructions": {
            "fitStep": "Measure pin radius, hole radius, slip, force, friction, and rebound fields in fit rows.",
            "freezeStep": "Estimate stiffness, damping, friction, solref, and solimp once, then freeze the profile.",
            "holdoutStep": "Repeat the same contact-pair measurements in a separate holdout dataset.",
            "acceptanceRule": (
                "Treat calibratedContact=true as meaningful only after fit and holdout rows "
                "are filled, residuals are compared, and missing measurement fields are zero."
            ),
        },
        "summary": {
            "status": (
                "contact-parameter-calibration-packet-ready"
                if ready
                else "needs-contact-parameter-calibration-evidence"
            ),
            "contactParameterCalibrationPacketReady": ready,
            "contactPairRecordCount": contact_pairs,
            "parameterRecordCount": parameter_records,
            "measurementColumnCount": len(measurement_columns),
            "fitTemplateRowCount": len(fit_rows),
            "holdoutTemplateRowCount": len(holdout_rows),
            "artifactManifestEntryCount": len(artifact_manifest),
            "profileEvidenceReady": profile_ready,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "contact_parameter_calibration_packet_completeness",
            "leanStructure": "Mechanics.ContactParameterCalibrationPacketNat",
            "leanPredicate": "contactParameterCalibrationPacketCompleteNat",
            "schema": "rad-sim.contact-parameter-calibration-packet.v1",
        },
        "claimLabels": {
            "packetAssembly": "bench protocol artifact",
            "parameterValues": "experimentally unvalidated physical assumption until fit/holdout measurements pass",
            "formalCompleteness": "Lean-proven finite packet-completeness predicate, not calibrated contact mechanics",
        },
        "limitations": [
            "The packet creates measurement rows but cannot prove that the measurements were collected.",
            "Fitted stiffness, damping, friction, and solver parameters remain empirical until independently validated.",
            "MuJoCo agreement is not hardware validation without bench measurements and provenance.",
        ],
    }


def export_contact_parameter_calibration_packet_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    **kwargs: object,
) -> str:
    return json.dumps(
        contact_parameter_calibration_packet(config, state, **kwargs),
        indent=2,
    )


def export_contact_parameter_calibration_packet_csv(packet: dict[str, object]) -> str:
    header = [
        "dataset_id",
        "dataset_role",
        "scenario_id",
        "repeat_index",
        "row",
        "col",
        "pivot",
        "pin",
        "hole",
        "assumed_contact_stiffness",
        "assumed_contact_damping",
        "assumed_friction",
        "assumed_solref",
        "assumed_solimp",
        "assumed_condim",
        "pin_radius_mm",
        "hole_radius_mm",
        "clearance_mm",
        "normal_load_n",
        "tangential_load_n",
        "imposed_z_mm",
        "measured_pin_hole_slip_mm",
        "measured_normal_force_n",
        "measured_tangent_force_n",
        "measured_static_friction_coeff",
        "measured_dynamic_friction_coeff",
        "fitted_contact_stiffness",
        "fitted_contact_damping",
        "fitted_solref_timeconst",
        "fitted_solref_damping_ratio",
        "fitted_solimp_width",
        "fixture_notes",
    ]
    rows = [header]
    for row in [
        *(packet.get("fitTemplateRows", []) if isinstance(packet.get("fitTemplateRows"), list) else []),
        *(packet.get("holdoutTemplateRows", []) if isinstance(packet.get("holdoutTemplateRows"), list) else []),
    ]:
        if not isinstance(row, dict):
            continue
        measured = row.get("requiredMeasurements", {})
        if not isinstance(measured, dict):
            measured = {}
        rows.append(
            [
                row.get("datasetId", ""),
                row.get("datasetRole", ""),
                row.get("scenarioId", ""),
                row.get("repeatIndex", ""),
                row.get("row", ""),
                row.get("col", ""),
                row.get("pivot", ""),
                row.get("pin", ""),
                row.get("hole", ""),
                row.get("assumedContactStiffness", ""),
                row.get("assumedContactDamping", ""),
                ";".join(str(value) for value in row.get("assumedFriction", []))
                if isinstance(row.get("assumedFriction"), list)
                else "",
                ";".join(str(value) for value in row.get("assumedSolref", []))
                if isinstance(row.get("assumedSolref"), list)
                else "",
                ";".join(str(value) for value in row.get("assumedSolimp", []))
                if isinstance(row.get("assumedSolimp"), list)
                else "",
                row.get("assumedCondim", ""),
                measured.get("pinRadiusMm", ""),
                measured.get("holeRadiusMm", ""),
                measured.get("clearanceMm", ""),
                measured.get("normalLoadN", ""),
                measured.get("tangentialLoadN", ""),
                measured.get("imposedZMm", ""),
                measured.get("measuredPinHoleSlipMm", ""),
                measured.get("measuredNormalForceN", ""),
                measured.get("measuredTangentForceN", ""),
                measured.get("measuredStaticFrictionCoeff", ""),
                measured.get("measuredDynamicFrictionCoeff", ""),
                measured.get("fittedContactStiffness", ""),
                measured.get("fittedContactDamping", ""),
                measured.get("fittedSolrefTimeconst", ""),
                measured.get("fittedSolrefDampingRatio", ""),
                measured.get("fittedSolimpWidth", ""),
                measured.get("fixtureNotes", ""),
            ]
        )
    summary = packet.get("summary", {})
    if isinstance(summary, dict):
        rows.append(
            [
                "summary",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                summary.get("parameterRecordCount", ""),
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                summary.get("contactParameterCalibrationPacketReady", ""),
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row)
        for row in rows
    )


CONTACT_PARAMETER_RESULT_FIELDS = (
    "pinRadiusMm",
    "holeRadiusMm",
    "clearanceMm",
    "normalLoadN",
    "tangentialLoadN",
    "imposedZMm",
    "measuredPinHoleSlipMm",
    "measuredNormalForceN",
    "measuredTangentForceN",
    "measuredReboundRatio",
    "measuredContactDurationS",
    "measuredStaticFrictionCoeff",
    "measuredDynamicFrictionCoeff",
    "fittedContactStiffness",
    "fittedContactDamping",
    "fittedSolrefTimeconst",
    "fittedSolrefDampingRatio",
    "fittedSolimpWidth",
)


def contact_parameter_calibration_results_template(
    packet: dict[str, object],
    *,
    dataset_id: str = "contact-parameter-results-001",
    dataset_role: str = "fit-holdout-results",
    source_file_id: str | None = None,
    operator: str = "",
    collected_at: str = "",
) -> dict[str, object]:
    """Create a fillable result artifact from a contact-parameter packet."""

    if packet.get("schema") != "rad-sim.contact-parameter-calibration-packet.v1":
        raise ValueError("unsupported contact-parameter calibration packet schema")
    fit_rows = packet.get("fitTemplateRows", [])
    holdout_rows = packet.get("holdoutTemplateRows", [])
    if not isinstance(fit_rows, list):
        fit_rows = []
    if not isinstance(holdout_rows, list):
        holdout_rows = []
    measurements = [
        deepcopy(row)
        for row in (*fit_rows, *holdout_rows)
        if isinstance(row, dict)
    ]
    return {
        "schema": "rad-sim.contact-parameter-calibration-results.v1",
        "sourcePacketSchema": packet.get("schema", ""),
        "sourcePacketTargetId": (
            packet.get("formalization", {}).get("targetId", "")
            if isinstance(packet.get("formalization"), dict)
            else ""
        ),
        "datasetId": dataset_id,
        "datasetRole": dataset_role,
        "sourceFileId": source_file_id or f"{dataset_id}.json",
        "operator": operator,
        "collectedAt": collected_at,
        "profilePlan": deepcopy(packet.get("profilePlan", {})),
        "measurementColumns": list(packet.get("measurementColumns", []))
        if isinstance(packet.get("measurementColumns"), list)
        else [],
        "provenanceNotes": (
            "Fill requiredMeasurements from bench or external contact fitting; "
            "do not edit assumedContact*, assumedFriction, assumedSolref, or assumedSolimp."
        ),
        "measurements": measurements,
    }


def export_contact_parameter_calibration_results_template_json(
    packet: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        contact_parameter_calibration_results_template(packet, **kwargs),
        indent=2,
    )


def contact_parameter_calibration_results_from_json(text: str) -> dict[str, object]:
    raw = json.loads(text)
    if not isinstance(raw, dict):
        raise ValueError("contact-parameter calibration results must be a JSON object")
    if raw.get("schema") != "rad-sim.contact-parameter-calibration-results.v1":
        raise ValueError("unsupported contact-parameter calibration results schema")
    if not isinstance(raw.get("measurements"), list):
        raise ValueError("contact-parameter calibration results require measurements")
    return raw


def _contact_parameter_result_float(row: dict[str, object], field: str) -> float | None:
    measurements = row.get("requiredMeasurements", {})
    if not isinstance(measurements, dict):
        return None
    try:
        value = float(measurements.get(field))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if not np.isfinite(value):
        return None
    return value


def _contact_parameter_row_float(row: dict[str, object], field: str) -> float | None:
    try:
        value = float(row.get(field))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if not np.isfinite(value):
        return None
    return value


def _contact_parameter_row_list_float(
    row: dict[str, object],
    field: str,
    index: int,
) -> float | None:
    values = row.get(field, [])
    if not isinstance(values, (list, tuple)) or index >= len(values):
        return None
    try:
        value = float(values[index])
    except (TypeError, ValueError):
        return None
    if not np.isfinite(value):
        return None
    return value


def _contact_parameter_expected_value(
    packet: dict[str, object],
    row: dict[str, object],
    field: str,
) -> float | None:
    grid = packet.get("grid", {})
    if not isinstance(grid, dict):
        grid = {}
    if field == "pinRadiusMm":
        return _contact_parameter_row_float({"value": grid.get("pinRadius")}, "value")
    if field == "holeRadiusMm":
        return _contact_parameter_row_float({"value": grid.get("holeRadius")}, "value")
    if field == "clearanceMm":
        return _contact_parameter_row_float(
            {"value": grid.get("pinHoleClearance")},
            "value",
        )
    if field == "fittedContactStiffness":
        return _contact_parameter_row_float(row, "assumedContactStiffness")
    if field == "fittedContactDamping":
        return _contact_parameter_row_float(row, "assumedContactDamping")
    if field == "measuredStaticFrictionCoeff":
        return _contact_parameter_row_list_float(row, "assumedFriction", 0)
    if field == "measuredDynamicFrictionCoeff":
        return _contact_parameter_row_list_float(row, "assumedFriction", 1)
    if field == "fittedSolrefTimeconst":
        return _contact_parameter_row_list_float(row, "assumedSolref", 0)
    if field == "fittedSolrefDampingRatio":
        return _contact_parameter_row_list_float(row, "assumedSolref", 1)
    if field == "fittedSolimpWidth":
        return _contact_parameter_row_list_float(row, "assumedSolimp", 0)
    return None


def _contact_parameter_key(row: dict[str, object]) -> tuple[int, int, str, str, str]:
    return (
        int(row.get("row", 0) or 0),
        int(row.get("col", 0) or 0),
        str(row.get("pivot", "")),
        str(row.get("pin", "")),
        str(row.get("hole", "")),
    )


def _mean(values: Iterable[float]) -> float | None:
    numeric = list(values)
    if not numeric:
        return None
    return float(sum(numeric) / len(numeric))


def compare_contact_parameter_calibration_results(
    packet: dict[str, object],
    results: dict[str, object] | None = None,
    *,
    tolerance: float = 1e-9,
    holdout_tolerance: float = 1e-6,
) -> dict[str, object]:
    """Compare filled contact-parameter rows against frozen profile assumptions."""

    if packet.get("schema") != "rad-sim.contact-parameter-calibration-packet.v1":
        raise ValueError("unsupported contact-parameter calibration packet schema")
    results = results or contact_parameter_calibration_results_template(packet)
    if results.get("schema") != "rad-sim.contact-parameter-calibration-results.v1":
        raise ValueError("unsupported contact-parameter calibration results schema")

    packet_summary = packet.get("summary", {})
    if not isinstance(packet_summary, dict):
        packet_summary = {}
    packet_ready = bool(packet_summary.get("contactParameterCalibrationPacketReady", False))
    expected_fit_rows = len(packet.get("fitTemplateRows", [])) if isinstance(packet.get("fitTemplateRows"), list) else 0
    expected_holdout_rows = len(packet.get("holdoutTemplateRows", [])) if isinstance(packet.get("holdoutTemplateRows"), list) else 0
    result_rows = [
        row for row in results.get("measurements", []) if isinstance(row, dict)
    ]
    fit_rows = [row for row in result_rows if row.get("datasetRole") == "fit"]
    holdout_rows = [row for row in result_rows if row.get("datasetRole") == "holdout"]
    compare_fields = (
        "pinRadiusMm",
        "holeRadiusMm",
        "clearanceMm",
        "measuredStaticFrictionCoeff",
        "measuredDynamicFrictionCoeff",
        "fittedContactStiffness",
        "fittedContactDamping",
        "fittedSolrefTimeconst",
        "fittedSolrefDampingRatio",
        "fittedSolimpWidth",
    )

    row_reports: list[dict[str, object]] = []
    role_values: dict[tuple[str, tuple[int, int, str, str, str], str], list[float]] = {}
    fit_residuals: list[float] = []
    holdout_residuals: list[float] = []
    missing_measurements = 0
    completed_measurement_rows = 0
    parameter_residual_count = 0
    for row in result_rows:
        missing_fields = [
            field
            for field in CONTACT_PARAMETER_RESULT_FIELDS
            if _contact_parameter_result_float(row, field) is None
        ]
        row_missing = len(missing_fields)
        missing_measurements += row_missing
        if row_missing == 0:
            completed_measurement_rows += 1
        residuals: list[dict[str, object]] = []
        max_row_residual = 0.0
        role = str(row.get("datasetRole", ""))
        key = _contact_parameter_key(row)
        for field in compare_fields:
            observed = _contact_parameter_result_float(row, field)
            expected = _contact_parameter_expected_value(packet, row, field)
            if observed is None or expected is None:
                continue
            residual = abs(observed - expected)
            parameter_residual_count += 1
            max_row_residual = max(max_row_residual, residual)
            residuals.append(
                {
                    "field": field,
                    "observed": observed,
                    "expected": expected,
                    "absResidual": residual,
                    "pass": residual <= (
                        holdout_tolerance if role == "holdout" else tolerance
                    ),
                }
            )
            role_values.setdefault((role, key, field), []).append(observed)
            if role == "holdout":
                holdout_residuals.append(residual)
            else:
                fit_residuals.append(residual)
        row_reports.append(
            {
                "datasetId": row.get("datasetId", ""),
                "datasetRole": role,
                "scenarioId": row.get("scenarioId", ""),
                "repeatIndex": row.get("repeatIndex", ""),
                "row": row.get("row", ""),
                "col": row.get("col", ""),
                "pivot": row.get("pivot", ""),
                "pin": row.get("pin", ""),
                "hole": row.get("hole", ""),
                "missingFields": missing_fields,
                "missingMeasurementCount": row_missing,
                "residuals": residuals,
                "maxParameterResidual": max_row_residual,
                "pass": row_missing == 0
                and all(bool(item.get("pass", False)) for item in residuals),
            }
        )

    holdout_pair_deltas: list[float] = []
    compared_pair_fields = 0
    for role, key, field in list(role_values):
        if role != "fit":
            continue
        fit_value = _mean(role_values.get(("fit", key, field), []))
        holdout_value = _mean(role_values.get(("holdout", key, field), []))
        if fit_value is None or holdout_value is None:
            continue
        compared_pair_fields += 1
        holdout_pair_deltas.append(abs(fit_value - holdout_value))

    max_fit_residual = max(fit_residuals) if fit_residuals else 0.0
    max_holdout_residual = max(holdout_residuals) if holdout_residuals else 0.0
    max_fit_holdout_delta = max(holdout_pair_deltas) if holdout_pair_deltas else 0.0
    fit_parameter_pass = (
        len(fit_rows) >= expected_fit_rows
        and len(fit_rows) > 0
        and all(report["missingMeasurementCount"] == 0 for report in row_reports if report["datasetRole"] == "fit")
        and bool(fit_residuals)
        and max_fit_residual <= tolerance
    )
    holdout_parameter_pass = (
        len(holdout_rows) >= expected_holdout_rows
        and len(holdout_rows) > 0
        and all(report["missingMeasurementCount"] == 0 for report in row_reports if report["datasetRole"] == "holdout")
        and bool(holdout_residuals)
        and max_holdout_residual <= holdout_tolerance
    )
    independent_holdout_pass = (
        compared_pair_fields > 0 and max_fit_holdout_delta <= holdout_tolerance
    )
    missing_evidence: list[str] = []
    if not packet_ready:
        missing_evidence.append("contactParameterCalibrationPacket")
    if expected_fit_rows <= 0:
        missing_evidence.append("fitTemplateRows")
    if expected_holdout_rows <= 0:
        missing_evidence.append("holdoutTemplateRows")
    if not result_rows:
        missing_evidence.append("resultRows")
    if len(fit_rows) < expected_fit_rows:
        missing_evidence.append("completedFitRows")
    if len(holdout_rows) < expected_holdout_rows:
        missing_evidence.append("completedHoldoutRows")
    if missing_measurements > 0:
        missing_evidence.append("completedMeasurements")
    if not fit_residuals:
        missing_evidence.append("fitParameterResiduals")
    if not holdout_residuals:
        missing_evidence.append("holdoutParameterResiduals")
    if compared_pair_fields <= 0:
        missing_evidence.append("independentHoldoutPairs")
    pass_flag = (
        len(missing_evidence) == 0
        and fit_parameter_pass
        and holdout_parameter_pass
        and independent_holdout_pass
    )
    return {
        "schema": "rad-sim.contact-parameter-bench-validation.v1",
        "method": (
            "comparison of filled fit and holdout contact-parameter measurements "
            "against the frozen proxy MuJoCo contact profile"
        ),
        "sourcePacket": {
            "schema": packet.get("schema", ""),
            "ready": packet_ready,
            "profilePlan": deepcopy(packet.get("profilePlan", {})),
        },
        "sourceResults": {
            "schema": results.get("schema", ""),
            "datasetId": results.get("datasetId", ""),
            "datasetRole": results.get("datasetRole", ""),
            "sourceFileId": results.get("sourceFileId", ""),
        },
        "metrics": {
            "resultRowCount": len(result_rows),
            "expectedFitRowCount": expected_fit_rows,
            "expectedHoldoutRowCount": expected_holdout_rows,
            "fitRowCount": len(fit_rows),
            "holdoutRowCount": len(holdout_rows),
            "completedMeasurementRowCount": completed_measurement_rows,
            "missingMeasurementCount": missing_measurements,
            "comparedPairFieldCount": compared_pair_fields,
            "parameterResidualCount": parameter_residual_count,
            "fitParameterResidualCount": len(fit_residuals),
            "holdoutParameterResidualCount": len(holdout_residuals),
            "maxFitParameterResidual": max_fit_residual,
            "maxHoldoutParameterResidual": max_holdout_residual,
            "maxFitHoldoutParameterDelta": max_fit_holdout_delta,
            "tolerance": float(tolerance),
            "holdoutTolerance": float(holdout_tolerance),
        },
        "summary": {
            "status": (
                "contact-parameter-bench-validation-pass"
                if pass_flag
                else "needs-contact-parameter-bench-review"
            ),
            "contactParameterBenchValidationPass": pass_flag,
            "fitParameterPass": fit_parameter_pass,
            "holdoutParameterPass": holdout_parameter_pass,
            "independentHoldoutPass": independent_holdout_pass,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "rows": row_reports,
        "formalization": {
            "targetId": "contact_parameter_bench_validation_gate",
            "leanStructure": "Mechanics.ContactParameterBenchValidationNat",
            "leanPredicate": "contactParameterBenchValidationReadyNat",
            "schema": "rad-sim.contact-parameter-bench-validation.v1",
        },
        "claimLabels": {
            "comparison": "bench comparison artifact",
            "parameterFit": "empirical comparison only; not a first-principles contact law",
            "holdout": "independent holdout check if collected after profile freeze",
            "formalCompleteness": "Lean-proven finite pass/fail predicate over result counts and flags",
        },
        "limitations": [
            "Passing this gate requires filled rows but does not prove friction, damping, or contact constitutive laws.",
            "Fit and holdout rows must come from independent collection after the profile is frozen.",
            "Radius and clearance fields are normalized proxy comparisons until the hardware unit scale is calibrated.",
        ],
    }


def export_contact_parameter_bench_validation_json(
    packet: dict[str, object],
    results: dict[str, object] | None = None,
    **kwargs: object,
) -> str:
    return json.dumps(
        compare_contact_parameter_calibration_results(packet, results, **kwargs),
        indent=2,
    )


def export_contact_parameter_bench_validation_csv(report: dict[str, object]) -> str:
    metrics = report.get("metrics", {})
    summary = report.get("summary", {})
    if not isinstance(metrics, dict):
        metrics = {}
    if not isinstance(summary, dict):
        summary = {}
    rows: list[list[object]] = [
        [
            "dataset_id",
            "dataset_role",
            "scenario_id",
            "repeat_index",
            "row",
            "col",
            "pivot",
            "pin",
            "hole",
            "missing_measurement_count",
            "max_parameter_residual",
            "pass",
        ]
    ]
    for row in report.get("rows", []):
        if not isinstance(row, dict):
            continue
        rows.append(
            [
                row.get("datasetId", ""),
                row.get("datasetRole", ""),
                row.get("scenarioId", ""),
                row.get("repeatIndex", ""),
                row.get("row", ""),
                row.get("col", ""),
                row.get("pivot", ""),
                row.get("pin", ""),
                row.get("hole", ""),
                row.get("missingMeasurementCount", ""),
                row.get("maxParameterResidual", ""),
                row.get("pass", ""),
            ]
        )
    rows.append(
        [
            "summary",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            metrics.get("missingMeasurementCount", ""),
            metrics.get("maxFitHoldoutParameterDelta", ""),
            summary.get("contactParameterBenchValidationPass", ""),
        ]
    )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row)
        for row in rows
    )


CONTACT_PARAMETER_INTERVAL_FIELDS = (
    ("pinRadius", "pinRadiusMm", "normalized radius"),
    ("holeRadius", "holeRadiusMm", "normalized radius"),
    ("pinHoleClearance", "clearanceMm", "normalized radius"),
    ("staticFrictionCoeff", "measuredStaticFrictionCoeff", "coefficient"),
    ("dynamicFrictionCoeff", "measuredDynamicFrictionCoeff", "coefficient"),
    ("contactStiffness", "fittedContactStiffness", "proxy stiffness"),
    ("contactDamping", "fittedContactDamping", "proxy damping"),
    ("solrefTimeconst", "fittedSolrefTimeconst", "MuJoCo solref"),
    ("solrefDampingRatio", "fittedSolrefDampingRatio", "MuJoCo solref"),
    ("solimpWidth", "fittedSolimpWidth", "MuJoCo solimp"),
)


def contact_parameter_interval_calibration_report(
    packet: dict[str, object],
    results: dict[str, object] | None = None,
    *,
    bench_validation: dict[str, object] | None = None,
    confidence_margin: float = 0.0,
    tolerance: float = 1e-9,
    holdout_tolerance: float = 1e-6,
) -> dict[str, object]:
    """Convert filled contact-parameter results into conservative parameter intervals."""

    if packet.get("schema") != "rad-sim.contact-parameter-calibration-packet.v1":
        raise ValueError("unsupported contact-parameter calibration packet schema")
    results = results or contact_parameter_calibration_results_template(packet)
    if results.get("schema") != "rad-sim.contact-parameter-calibration-results.v1":
        raise ValueError("unsupported contact-parameter calibration results schema")
    validation = bench_validation or compare_contact_parameter_calibration_results(
        packet,
        results,
        tolerance=tolerance,
        holdout_tolerance=holdout_tolerance,
    )
    validation_summary = validation.get("summary", {})
    if not isinstance(validation_summary, dict):
        validation_summary = {}
    validation_pass = bool(
        validation_summary.get("contactParameterBenchValidationPass", False)
    )
    result_rows = [
        row for row in results.get("measurements", []) if isinstance(row, dict)
    ]
    margin = max(0.0, float(confidence_margin))
    intervals: list[dict[str, object]] = []
    accepted_count = 0
    uncertainty_count = 0
    holdout_agreement_count = 0
    for parameter_name, measurement_field, unit in CONTACT_PARAMETER_INTERVAL_FIELDS:
        observed_values: list[float] = []
        fit_values: list[float] = []
        holdout_values: list[float] = []
        residuals: list[float] = []
        assumed_values: list[float] = []
        for row in result_rows:
            observed = _contact_parameter_result_float(row, measurement_field)
            expected = _contact_parameter_expected_value(packet, row, measurement_field)
            if expected is not None:
                assumed_values.append(expected)
            if observed is None:
                continue
            observed_values.append(observed)
            if row.get("datasetRole") == "holdout":
                holdout_values.append(observed)
            else:
                fit_values.append(observed)
            if expected is not None:
                residuals.append(abs(observed - expected))
        assumed = _mean(assumed_values)
        mean_value = _mean(observed_values)
        fit_mean = _mean(fit_values)
        holdout_mean = _mean(holdout_values)
        sample_min = min(observed_values) if observed_values else None
        sample_max = max(observed_values) if observed_values else None
        max_residual = max(residuals) if residuals else None
        holdout_delta = (
            abs(fit_mean - holdout_mean)
            if fit_mean is not None and holdout_mean is not None
            else None
        )
        if mean_value is None or sample_min is None or sample_max is None:
            half_width = None
            lower = None
            upper = None
        else:
            spread_half_width = max(abs(mean_value - sample_min), abs(sample_max - mean_value))
            evidence_half_width = max(
                spread_half_width,
                max_residual if max_residual is not None else 0.0,
                holdout_delta if holdout_delta is not None else 0.0,
                margin,
            )
            half_width = evidence_half_width
            lower = mean_value - evidence_half_width
            upper = mean_value + evidence_half_width
        accepted = (
            assumed is not None
            and lower is not None
            and upper is not None
            and lower <= assumed <= upper
        )
        if accepted:
            accepted_count += 1
        if half_width is not None:
            uncertainty_count += 1
        if holdout_delta is not None:
            holdout_agreement_count += 1
        intervals.append(
            {
                "parameter": parameter_name,
                "measurementField": measurement_field,
                "unit": unit,
                "assumedValue": assumed,
                "meanObserved": mean_value,
                "fitMean": fit_mean,
                "holdoutMean": holdout_mean,
                "sampleMin": sample_min,
                "sampleMax": sample_max,
                "halfWidth": half_width,
                "lowerBound": lower,
                "upperBound": upper,
                "maxResidualToAssumption": max_residual,
                "fitHoldoutDelta": holdout_delta,
                "sampleCount": len(observed_values),
                "fitSampleCount": len(fit_values),
                "holdoutSampleCount": len(holdout_values),
                "acceptedAssumptionInsideBounds": accepted,
                "claimLabel": "empirical interval estimate from filled fit/holdout rows",
            }
        )

    missing_evidence: list[str] = []
    if not validation_pass:
        missing_evidence.append("contactParameterBenchValidationPass")
    if not result_rows:
        missing_evidence.append("resultRows")
    if len(intervals) < len(CONTACT_PARAMETER_INTERVAL_FIELDS):
        missing_evidence.append("parameterIntervals")
    if accepted_count < len(intervals):
        missing_evidence.append("simulatorParametersInsideBounds")
    if uncertainty_count < len(intervals):
        missing_evidence.append("uncertaintyBounds")
    if holdout_agreement_count < len(intervals):
        missing_evidence.append("holdoutAgreement")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.contact-parameter-interval-calibration.v1",
        "method": (
            "conservative interval calibration over filled fit and holdout "
            "pin-hole/contact-parameter measurements"
        ),
        "sourcePacket": {
            "schema": packet.get("schema", ""),
            "profilePlan": deepcopy(packet.get("profilePlan", {})),
        },
        "sourceResults": {
            "schema": results.get("schema", ""),
            "datasetId": results.get("datasetId", ""),
            "datasetRole": results.get("datasetRole", ""),
            "sourceFileId": results.get("sourceFileId", ""),
        },
        "sourceBenchValidation": {
            "schema": validation.get("schema", ""),
            "passed": validation_pass,
            "targetId": (
                validation.get("formalization", {}).get("targetId", "")
                if isinstance(validation.get("formalization"), dict)
                else ""
            ),
        },
        "parameters": intervals,
        "summary": {
            "status": (
                "contact-parameter-interval-calibration-ready"
                if ready
                else "needs-contact-parameter-interval-review"
            ),
            "contactParameterIntervalCalibrationReady": ready,
            "benchValidationPass": validation_pass,
            "parameterIntervalCount": len(intervals),
            "acceptedParameterIntervalCount": accepted_count,
            "uncertaintyRecordCount": uncertainty_count,
            "holdoutAgreementRecordCount": holdout_agreement_count,
            "simulatorParametersInsideBounds": accepted_count == len(intervals),
            "confidenceMargin": margin,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "contact_parameter_interval_calibration_gate",
            "leanStructure": "Mechanics.ContactParameterIntervalCalibrationNat",
            "leanPredicate": "contactParameterIntervalCalibrationReadyNat",
            "schema": "rad-sim.contact-parameter-interval-calibration.v1",
        },
        "claimLabels": {
            "intervals": "empirical bounded-parameter artifact",
            "acceptance": "frozen simulator parameter lies inside measured interval",
            "formalCompleteness": "Lean-proven finite readiness predicate over interval evidence",
            "physics": "not a first-principles contact mechanics proof",
        },
        "limitations": [
            "Intervals are only as valid as the bench/result rows and measurement provenance.",
            "A parameter inside an empirical interval can still fail outside the tested load, velocity, wear, or assembly regime.",
            "This artifact bounds proxy MuJoCo/simulator parameters; it does not derive friction or damping from continuum mechanics.",
        ],
    }


def export_contact_parameter_interval_calibration_json(
    packet: dict[str, object],
    results: dict[str, object] | None = None,
    **kwargs: object,
) -> str:
    return json.dumps(
        contact_parameter_interval_calibration_report(packet, results, **kwargs),
        indent=2,
    )


def export_contact_parameter_interval_calibration_csv(report: dict[str, object]) -> str:
    summary = report.get("summary", {})
    if not isinstance(summary, dict):
        summary = {}
    rows: list[list[object]] = [
        [
            "parameter",
            "measurement_field",
            "unit",
            "assumed_value",
            "mean_observed",
            "fit_mean",
            "holdout_mean",
            "lower_bound",
            "upper_bound",
            "half_width",
            "max_residual_to_assumption",
            "fit_holdout_delta",
            "sample_count",
            "accepted_assumption_inside_bounds",
        ]
    ]
    for parameter in report.get("parameters", []):
        if not isinstance(parameter, dict):
            continue
        rows.append(
            [
                parameter.get("parameter", ""),
                parameter.get("measurementField", ""),
                parameter.get("unit", ""),
                parameter.get("assumedValue", ""),
                parameter.get("meanObserved", ""),
                parameter.get("fitMean", ""),
                parameter.get("holdoutMean", ""),
                parameter.get("lowerBound", ""),
                parameter.get("upperBound", ""),
                parameter.get("halfWidth", ""),
                parameter.get("maxResidualToAssumption", ""),
                parameter.get("fitHoldoutDelta", ""),
                parameter.get("sampleCount", ""),
                parameter.get("acceptedAssumptionInsideBounds", ""),
            ]
        )
    rows.append(
        [
            "summary",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            summary.get("confidenceMargin", ""),
            "",
            "",
            summary.get("parameterIntervalCount", ""),
            summary.get("contactParameterIntervalCalibrationReady", ""),
        ]
    )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row)
        for row in rows
    )


def mujoco_external_run_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    load_case: LoadCase | None = None,
    export_report: dict[str, object] | None = None,
    steps: int = 120,
    engine_availability: dict[str, bool] | None = None,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Run the exported MuJoCo model when the independent package is available."""

    if steps <= 0:
        raise ValueError("steps must be positive")
    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    load = _mujoco_load_case(load_case)
    resolved_export = export_report or mujoco_model_export_report(
        config,
        base,
        event_sequence=event_sequence,
        load_case=load,
    )
    export_summary = (
        resolved_export.get("summary", {}) if isinstance(resolved_export, dict) else {}
    )
    available = _engine_available("mujoco", engine_availability)
    missing_evidence: list[str] = []
    if not export_summary.get("mujocoModelExportReady"):
        missing_evidence.append("mujocoModelExport")
    if not available:
        missing_evidence.append("mujocoPackage")
    if int(export_summary.get("bodyRecordCount", 0) or 0) <= 0:
        missing_evidence.append("bodyRecords")
    if missing_evidence:
        return {
            "schema": "rad-sim.mujoco-external-run.v1",
            "method": "optional independent MuJoCo execution of the exported RAD proxy model",
            "engine": "MuJoCo",
            "modelExport": resolved_export,
            "solver": {
                "engineAvailable": available,
                "ran": False,
                "stepsRequested": steps,
                "stepsCompleted": 0,
                "error": "",
            },
            "results": {"bodies": []},
            "summary": {
                "status": "needs-mujoco-external-run-evidence",
                "mujocoExternalRunComplete": False,
                "engineAvailable": available,
                "bodyResultCount": 0,
                "expectedBodyResultCount": int(export_summary.get("bodyRecordCount", 0) or 0),
                "missingEvidenceCount": len(missing_evidence),
                "missingEvidence": missing_evidence,
            },
            "formalization": {
                "targetId": "mujoco_external_run_gate",
                "leanStructure": "Mechanics.ExternalPhysicsRunNat",
                "leanPredicate": "externalPhysicsRunReadyNat",
                "schema": "rad-sim.mujoco-external-run.v1",
            },
            "claimLabels": {
                "run": "external-engine run not completed",
                "physicalAccuracy": "experimentally unvalidated",
            },
        }

    bodies = [
        body for body in resolved_export.get("bodies", []) if isinstance(body, dict)
    ]
    results: list[dict[str, object]] = []
    error = ""
    steps_completed = 0
    try:
        mujoco = importlib.import_module("mujoco")
        model = mujoco.MjModel.from_xml_string(str(resolved_export["xml"]))
        data = mujoco.MjData(model)
        force_by_name: dict[str, tuple[float, float, float]] = {}
        for cell, force in load.external_forces.items():
            force_by_name[_mujoco_cell_name(int(cell[0]), int(cell[1]))] = _vec3(force)
        body_ids: dict[str, int] = {}
        for body in bodies:
            name = str(body["name"])
            body_id = int(mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, name))
            if body_id >= 0:
                body_ids[name] = body_id
        for _ in range(steps):
            data.xfrc_applied[:] = 0
            for name, force in force_by_name.items():
                body_id = body_ids.get(name)
                if body_id is None:
                    continue
                data.xfrc_applied[body_id, 0:3] = force
            mujoco.mj_step(model, data)
            steps_completed += 1
        for body in bodies:
            name = str(body["name"])
            body_id = body_ids.get(name)
            if body_id is None:
                continue
            initial = tuple(float(value) for value in body.get("position", (0, 0, 0)))
            final = tuple(float(value) for value in data.xpos[body_id, 0:3])
            displacement = tuple(final[index] - initial[index] for index in range(3))
            results.append(
                {
                    "row": int(body["row"]),
                    "col": int(body["col"]),
                    "name": name,
                    "initialPosition": list(initial),
                    "finalPosition": list(final),
                    "displacement": list(displacement),
                    "fixed": bool(body.get("fixed", False)),
                }
            )
    except Exception as exc:  # pragma: no cover - depends on optional external engine
        error = str(exc)
        missing_evidence.append("engineExecution")

    expected = int(export_summary.get("bodyRecordCount", len(bodies)) or 0)
    if len(results) < expected:
        missing_evidence.append("bodyResultCoverage")
    run_complete = (
        available
        and not error
        and steps_completed == steps
        and len(results) >= expected
        and len(missing_evidence) == 0
    )
    return {
        "schema": "rad-sim.mujoco-external-run.v1",
        "method": "optional independent MuJoCo execution of the exported RAD proxy model",
        "engine": "MuJoCo",
        "modelExport": resolved_export,
        "solver": {
            "engineAvailable": available,
            "ran": run_complete,
            "stepsRequested": steps,
            "stepsCompleted": steps_completed,
            "error": error,
            "tolerance": tolerance,
        },
        "results": {"bodies": results},
        "summary": {
            "status": "mujoco-external-run-complete" if run_complete else "needs-mujoco-external-run-evidence",
            "mujocoExternalRunComplete": run_complete,
            "engineAvailable": available,
            "bodyResultCount": len(results),
            "expectedBodyResultCount": expected,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "mujoco_external_run_gate",
            "leanStructure": "Mechanics.ExternalPhysicsRunNat",
            "leanPredicate": "externalPhysicsRunReadyNat",
            "schema": "rad-sim.mujoco-external-run.v1",
        },
        "claimLabels": {
            "run": "external-engine-derived result" if run_complete else "external-engine run incomplete",
            "physicalAccuracy": "experimentally unvalidated until bench comparison also passes",
        },
    }


def export_mujoco_external_run_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    **kwargs: object,
) -> str:
    return json.dumps(mujoco_external_run_report(config, state, **kwargs), indent=2)


def mujoco_external_comparison_report(
    config: LatticeConfig,
    state: LatticeState | None,
    run_report: dict[str, object],
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    tolerance: float = 1e-6,
) -> dict[str, object]:
    """Compare simulator centers against a MuJoCo run result or imported result."""

    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    sequence = tuple(event_sequence or ())
    comparison_state = _apply_realization_event_sequence(config, base, sequence)
    sim = simulate_kinematic(config, comparison_state)
    run_summary = run_report.get("summary", {}) if isinstance(run_report, dict) else {}
    body_results = (
        run_report.get("results", {}).get("bodies", [])
        if isinstance(run_report.get("results", {}), dict)
        else []
    )
    result_map: dict[tuple[int, int], dict[str, object]] = {}
    for result in body_results:
        if not isinstance(result, dict):
            continue
        result_map[(int(result.get("row", -1)), int(result.get("col", -1)))] = result

    records: list[dict[str, object]] = []
    squared_errors: list[float] = []
    missing_cells = 0
    for row in range(config.rows):
        for col in range(config.cols):
            if bool(comparison_state.removed_mask[row, col]):
                continue
            expected = tuple(float(value) for value in sim.deformed_centers_3d[row, col])
            result = result_map.get((row, col))
            if result is None:
                missing_cells += 1
                records.append(
                    {
                        "row": row,
                        "col": col,
                        "matched": False,
                        "simulatorPosition": list(expected),
                        "externalPosition": None,
                        "positionError": None,
                    }
                )
                continue
            external = _vec3(result.get("finalPosition", (0.0, 0.0, 0.0)))  # type: ignore[arg-type]
            error = float(
                np.linalg.norm(np.asarray(external, dtype=float) - np.asarray(expected, dtype=float))
            )
            squared_errors.append(error * error)
            records.append(
                {
                    "row": row,
                    "col": col,
                    "matched": True,
                    "simulatorPosition": list(expected),
                    "externalPosition": list(external),
                    "positionError": error,
                    "withinTolerance": error <= tolerance,
                }
            )
    rms_error = float(np.sqrt(_mean(squared_errors))) if squared_errors else 0.0
    max_error = max((float(record.get("positionError") or 0.0) for record in records), default=0.0)
    missing_evidence: list[str] = []
    if not run_summary.get("mujocoExternalRunComplete"):
        missing_evidence.append("mujocoExternalRun")
    if not records:
        missing_evidence.append("comparisonRecords")
    if missing_cells:
        missing_evidence.append("matchedBodyRecords")
    if max_error > tolerance:
        missing_evidence.append("positionTolerance")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.mujoco-external-comparison.v1",
        "method": "finite simulator-vs-MuJoCo body-center comparison for external validation review",
        "engine": "MuJoCo",
        "summary": {
            "status": "mujoco-external-comparison-ready" if ready else "needs-mujoco-external-comparison-evidence",
            "mujocoExternalComparisonReady": ready,
            "comparisonRecordCount": len(records),
            "matchedBodyRecordCount": len(records) - missing_cells,
            "missingBodyRecordCount": missing_cells,
            "rmsPositionError": rms_error,
            "maxPositionError": max_error,
            "tolerance": tolerance,
            "withinTolerance": max_error <= tolerance,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "records": records,
        "runSummary": run_summary,
        "formalization": {
            "targetId": "mujoco_external_comparison_gate",
            "leanStructure": "Mechanics.ExternalPhysicsComparisonNat",
            "leanPredicate": "externalPhysicsComparisonReadyNat",
            "schema": "rad-sim.mujoco-external-comparison.v1",
        },
        "claimLabels": {
            "comparison": "external-engine-derived comparison if the run report is complete",
            "simulatorAgreement": "numerical agreement diagnostic, not physical bench validation",
            "physicalAccuracy": "experimentally unvalidated until bench measurements agree",
        },
        "limitations": [
            "This compares body-center positions only.",
            "A passing external comparison does not validate contact, friction, or actuator mechanics without bench data.",
            "The comparison inherits the coarse MJCF export assumptions.",
        ],
    }


def export_mujoco_external_comparison_json(
    config: LatticeConfig,
    state: LatticeState | None,
    run_report: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        mujoco_external_comparison_report(config, state, run_report, **kwargs),
        indent=2,
    )


def export_mujoco_external_comparison_csv(report: dict[str, object]) -> str:
    rows = [[
        "row",
        "col",
        "matched",
        "simulator_position",
        "external_position",
        "position_error",
        "within_tolerance",
    ]]
    for record in report.get("records", []):
        if not isinstance(record, dict):
            continue
        rows.append(
            [
                record.get("row", ""),
                record.get("col", ""),
                record.get("matched", ""),
                ";".join(str(value) for value in record.get("simulatorPosition", []))
                if isinstance(record.get("simulatorPosition"), list)
                else "",
                ";".join(str(value) for value in record.get("externalPosition", []))
                if isinstance(record.get("externalPosition"), list)
                else "",
                record.get("positionError", ""),
                record.get("withinTolerance", ""),
            ]
        )
    summary = report.get("summary", {})
    if isinstance(summary, dict):
        rows.append(
            [
                "summary",
                "",
                summary.get("mujocoExternalComparisonReady", ""),
                "",
                "",
                summary.get("maxPositionError", ""),
                summary.get("withinTolerance", ""),
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row)
        for row in rows
    )


def _apply_realization_event_sequence(
    config: LatticeConfig,
    state: LatticeState,
    events: Iterable[ProgrammableDiscontinuityEvent] | None,
) -> LatticeState:
    current = state.normalized(config)
    for event in tuple(events or ()):
        current = apply_programmable_event(config, current, event)
    return current


def _metadata_energy_terms(metadata: dict[str, object]) -> dict[str, float]:
    return {
        "objectiveEnergy": float(metadata.get("objective_energy", metadata.get("energy", 0.0))),
        "storedEnergy": float(metadata.get("stored_energy", 0.0)),
        "springEnergy": float(metadata.get("spring_energy", 0.0)),
        "hingeEnergy": float(metadata.get("hinge_energy", 0.0)),
        "lockPenaltyEnergy": float(metadata.get("lock_penalty_energy", 0.0)),
        "externalPotentialEnergy": float(metadata.get("external_potential_energy", 0.0)),
    }


def equilibrium_relation_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    realization_report: dict[str, object] | None = None,
    contact_report: dict[str, object] | None = None,
    contact_graph_report: dict[str, object] | None = None,
    load_case: LoadCase | None = None,
    contact_stiffness: float = 1.0,
    tolerance: float = 1e-6,
    residual_tolerance: float | None = None,
) -> dict[str, object]:
    """Attach a finite equilibrium relation to realization, energy, and contact evidence."""

    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    sequence = tuple(event_sequence or ())
    equilibrium_state = _apply_realization_event_sequence(config, base, sequence)
    resolved_realization = realization_report or physical_realization_map_report(
        config,
        base,
        event_sequence=sequence if event_sequence is not None else None,
        contact_graph_report=contact_graph_report,
        contact_report=contact_report,
        contact_stiffness=contact_stiffness,
        tolerance=tolerance,
    )
    resolved_contact = contact_report or contact_state_abstraction_report(
        config,
        equilibrium_state,
        contact_stiffness=contact_stiffness,
        tolerance=tolerance,
    )
    solver = solve_spring_hinge_3d(config, equilibrium_state, load_case or LoadCase())
    metadata = dict(solver.metadata)
    residual_limit = float(tolerance if residual_tolerance is None else residual_tolerance)
    energies = _metadata_energy_terms(metadata)
    contact_penalty = 0.0
    if isinstance(resolved_contact, dict):
        contact_penalty = float(
            resolved_contact.get("summary", {}).get("totalContactPenalty", 0.0)
        )
    energy_terms = {
        **energies,
        "contactPenaltyEnergy": contact_penalty,
        "loadWorkMagnitudeProxy": abs(energies["externalPotentialEnergy"]),
    }
    nonnegative_keys = (
        "springEnergy",
        "hingeEnergy",
        "lockPenaltyEnergy",
        "contactPenaltyEnergy",
        "loadWorkMagnitudeProxy",
    )
    nonnegative_count = sum(
        1 for key in nonnegative_keys if energy_terms[key] >= -tolerance
    )
    objective_balance_error = abs(
        energies["objectiveEnergy"]
        - (energies["storedEnergy"] + energies["externalPotentialEnergy"])
    )
    residual = float(metadata.get("target_rms_error", 0.0))
    realization_summary = (
        resolved_realization.get("summary", {})
        if isinstance(resolved_realization, dict)
        else {}
    )
    contact_ready = bool(
        isinstance(resolved_contact, dict)
        and resolved_contact.get("summary", {}).get("contactStateAbstractionReady")
    )
    missing_evidence: list[str] = []
    if not realization_summary.get("physicalRealizationMapReady"):
        missing_evidence.append("physicalRealizationMap")
    if not bool(metadata.get("success")):
        missing_evidence.append("solverSuccess")
    if residual > residual_limit:
        missing_evidence.append("equilibriumResidual")
    if nonnegative_count < len(nonnegative_keys):
        missing_evidence.append("nonnegativeEnergyTerms")
    if objective_balance_error > tolerance:
        missing_evidence.append("energyBalance")
    if not contact_ready:
        missing_evidence.append("contactStateAbstraction")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.equilibrium-relation.v1",
        "method": (
            "finite equilibrium evidence gate tying a realized operator state "
            "to spring-hinge solver status, nonnegative stored/proxy energy "
            "terms, contact-state evidence, and residual tolerance"
        ),
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "backlash": config.backlash,
            "pinHoleClearance": config.pin_hole_clearance,
        },
        "solver": {
            "model": metadata.get("model", ""),
            "success": bool(metadata.get("success")),
            "message": str(metadata.get("message", "")),
            "iterations": int(metadata.get("iterations", 0)),
            "springEdges": int(metadata.get("spring_edges", 0)),
            "hingeTriples": int(metadata.get("hinge_triples", 0)),
            "removedCells": int(metadata.get("removed_cells", 0)),
        },
        "residual": {
            "targetRmsError": residual,
            "tolerance": residual_limit,
            "passesTolerance": residual <= residual_limit,
            "objectiveBalanceError": objective_balance_error,
        },
        "energy": {
            **energy_terms,
            "nonnegativeTermCount": nonnegative_count,
            "requiredNonnegativeTermCount": len(nonnegative_keys),
        },
        "realization": {
            "attached": isinstance(resolved_realization, dict),
            "schema": resolved_realization.get("schema", "")
            if isinstance(resolved_realization, dict)
            else "",
            "physicalRealizationMapReady": bool(
                realization_summary.get("physicalRealizationMapReady")
            ),
            "abstractOperatorCount": int(
                realization_summary.get("abstractOperatorCount", 0)
            ),
        },
        "contact": {
            "attached": isinstance(resolved_contact, dict),
            "schema": resolved_contact.get("schema", "")
            if isinstance(resolved_contact, dict)
            else "",
            "contactStateAbstractionReady": contact_ready,
            "contactPenaltyEnergy": contact_penalty,
        },
        "summary": {
            "status": "equilibrium-relation-ready" if ready else "needs-equilibrium-evidence",
            "equilibriumRelationReady": ready,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "equilibrium_relation_gate",
            "leanStructure": "Mechanics.EquilibriumRelationNat",
            "leanPredicate": "equilibriumRelationReadyNat",
            "schema": "rad-sim.equilibrium-relation.v1",
        },
        "claimLabels": {
            "equilibriumRelation": "Lean-proven finite evidence gate",
            "energyTerms": "simulator-derived spring/contact/load proxy",
            "solver": "numerical spring-hinge quasistatic preview",
            "physicalAccuracy": "experimentally unvalidated physical assumption",
        },
        "limitations": [
            "The equilibrium relation is a finite evidence gate over the current solver output.",
            "It does not prove continuous minimization, material behavior, friction, or contact mechanics.",
            "Residual tolerance is a user-chosen numerical threshold, not a hardware-calibrated bound.",
        ],
    }


def export_equilibrium_relation_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    realization_report: dict[str, object] | None = None,
    contact_report: dict[str, object] | None = None,
    contact_graph_report: dict[str, object] | None = None,
    load_case: LoadCase | None = None,
    contact_stiffness: float = 1.0,
    tolerance: float = 1e-6,
    residual_tolerance: float | None = None,
) -> str:
    return json.dumps(
        equilibrium_relation_report(
            config,
            state,
            event_sequence=event_sequence,
            realization_report=realization_report,
            contact_report=contact_report,
            contact_graph_report=contact_graph_report,
            load_case=load_case,
            contact_stiffness=contact_stiffness,
            tolerance=tolerance,
            residual_tolerance=residual_tolerance,
        ),
        indent=2,
    )


def export_equilibrium_relation_csv(report: dict[str, object]) -> str:
    solver = report.get("solver", {})
    residual = report.get("residual", {})
    energy = report.get("energy", {})
    summary = report.get("summary", {})
    if not isinstance(solver, dict):
        solver = {}
    if not isinstance(residual, dict):
        residual = {}
    if not isinstance(energy, dict):
        energy = {}
    if not isinstance(summary, dict):
        summary = {}
    header = [
        "schema",
        "status",
        "equilibrium_relation_ready",
        "solver_success",
        "target_rms_error",
        "residual_tolerance",
        "objective_balance_error",
        "stored_energy",
        "spring_energy",
        "hinge_energy",
        "lock_penalty_energy",
        "contact_penalty_energy",
        "load_work_magnitude_proxy",
        "missing_evidence",
    ]
    row = [
        report.get("schema", ""),
        summary.get("status", ""),
        summary.get("equilibriumRelationReady", ""),
        solver.get("success", ""),
        residual.get("targetRmsError", ""),
        residual.get("tolerance", ""),
        residual.get("objectiveBalanceError", ""),
        energy.get("storedEnergy", ""),
        energy.get("springEnergy", ""),
        energy.get("hingeEnergy", ""),
        energy.get("lockPenaltyEnergy", ""),
        energy.get("contactPenaltyEnergy", ""),
        energy.get("loadWorkMagnitudeProxy", ""),
        ";".join(str(item) for item in summary.get("missingEvidence", []))
        if isinstance(summary.get("missingEvidence"), list)
        else "",
    ]
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in values)
        for values in (header, row)
    )


def _event_actuator_cells(
    event: ProgrammableDiscontinuityEvent,
) -> tuple[tuple[int, int], ...]:
    if event.kind == "actuate" and event.cell is not None:
        return ((int(event.cell[0]), int(event.cell[1])),)
    if event.kind == "group_actuate":
        return tuple((int(row), int(col)) for row, col in event.cells)
    return ()


def _actuator_basis_for_reachability(
    config: LatticeConfig,
    initial: LatticeState,
    final: LatticeState,
    events: Iterable[ProgrammableDiscontinuityEvent],
    explicit: Iterable[tuple[int, int]] | None,
    tolerance: float,
) -> tuple[tuple[int, int], ...]:
    cells: list[tuple[int, int]] = []
    if explicit is not None:
        cells.extend((int(row), int(col)) for row, col in explicit)
    else:
        cells.extend(_state_command_support(initial, tolerance))
        cells.extend(_state_command_support(final, tolerance))
        for event in events:
            cells.extend(_event_actuator_cells(event))
    valid = (
        (row, col)
        for row, col in cells
        if 0 <= row < config.rows
        and 0 <= col < config.cols
        and not bool(final.removed_mask[row, col])
    )
    return _unique_cells(valid)


def _target_cells_for_reachability(
    config: LatticeConfig,
    state: LatticeState,
    target_cells: Iterable[tuple[int, int]] | None,
) -> tuple[tuple[int, int], ...]:
    if target_cells is None:
        return tuple(
            (row, col)
            for row in range(config.rows)
            for col in range(config.cols)
            if not bool(state.removed_mask[row, col])
        )
    return _unique_cells(
        (int(row), int(col))
        for row, col in target_cells
        if 0 <= int(row) < config.rows and 0 <= int(col) < config.cols
    )


def _reachable_map(matrix: np.ndarray, config: LatticeConfig, tolerance: float) -> np.ndarray:
    if matrix.size == 0 or matrix.shape[1] == 0:
        return np.zeros((config.rows, config.cols), dtype=bool)
    return np.any(np.abs(matrix) > tolerance, axis=1).reshape((config.rows, config.cols))


def _topology_component_reachable_map(
    topology: dict[str, object],
    actuator_cells: tuple[tuple[int, int], ...],
) -> np.ndarray:
    labels = np.asarray(topology["component_labels"], dtype=int)
    active_components = {
        int(labels[row, col])
        for row, col in actuator_cells
        if labels[row, col] >= 0
    }
    return np.isin(labels, list(active_components)) & (labels >= 0)


def reachable_equilibrium_controllability_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    actuator_cells: Iterable[tuple[int, int]] | None = None,
    target_cells: Iterable[tuple[int, int]] | None = None,
    equilibrium_report: dict[str, object] | None = None,
    load_case: LoadCase | None = None,
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    include_alpha: bool = True,
    include_z: bool = True,
    require_full_target_reachability: bool = False,
    tolerance: float = 1e-9,
    residual_tolerance: float | None = None,
) -> dict[str, object]:
    """Report finite reachable equilibrium evidence for a chosen actuator basis."""

    base = LatticeState.uniform(config) if state is None else state.normalized(config)
    sequence = tuple(event_sequence or ())
    final = _apply_realization_event_sequence(config, base, sequence)
    basis = _actuator_basis_for_reachability(
        config, base, final, sequence, actuator_cells, tolerance
    )
    targets = _target_cells_for_reachability(config, final, target_cells)
    matrix = build_response_matrix(
        config,
        basis,
        alpha_step=alpha_step,
        z_step=z_step,
        include_alpha=include_alpha,
        include_z=include_z,
        locked_cells=tuple(
            (row, col)
            for row in range(config.rows)
            for col in range(config.cols)
            if bool(final.locked_mask[row, col])
        ),
        removed_cells=tuple(
            (row, col)
            for row in range(config.rows)
            for col in range(config.cols)
            if bool(final.removed_mask[row, col])
        ),
        tolerance=tolerance,
    )
    topology = lattice_topology_diagnostic(config, final)
    alpha_map = _reachable_map(matrix.alpha, config, tolerance)
    height_map = _reachable_map(matrix.height, config, tolerance)
    topology_map = _topology_component_reachable_map(topology, basis)
    target_count = len(targets)
    alpha_target_reached = sum(1 for row, col in targets if bool(alpha_map[row, col]))
    height_target_reached = sum(1 for row, col in targets if bool(height_map[row, col]))
    topology_target_reached = sum(1 for row, col in targets if bool(topology_map[row, col]))
    alpha_under = max(0, target_count - alpha_target_reached)
    height_under = max(0, target_count - height_target_reached)
    topology_blocked = max(0, target_count - topology_target_reached)
    resolved_equilibrium = equilibrium_report or equilibrium_relation_report(
        config,
        base,
        event_sequence=sequence,
        load_case=load_case,
        tolerance=tolerance,
        residual_tolerance=residual_tolerance,
    )
    equilibrium_summary = (
        resolved_equilibrium.get("summary", {})
        if isinstance(resolved_equilibrium, dict)
        else {}
    )
    missing_evidence: list[str] = []
    if not equilibrium_summary.get("equilibriumRelationReady"):
        missing_evidence.append("equilibriumRelation")
    if not basis:
        missing_evidence.append("actuatorBasis")
    if len(matrix.commands) == 0:
        missing_evidence.append("responseColumns")
    if target_count <= 0:
        missing_evidence.append("targetCells")
    if matrix.reachable_alpha_cells(tolerance) + matrix.reachable_height_cells(tolerance) <= 0:
        missing_evidence.append("reachableResponse")
    if require_full_target_reachability and (alpha_under > 0 or height_under > 0 or topology_blocked > 0):
        missing_evidence.append("targetReachability")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.reachable-equilibrium-controllability.v1",
        "method": (
            "finite actuator-basis reachability report over an equilibrium "
            "state, response matrix, and removed-cell topology graph"
        ),
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "totalCells": config.rows * config.cols,
        },
        "actuatorBasis": {
            "cells": [_cell_dict(cell) for cell in basis],
            "cellCount": len(basis),
            "commandColumnCount": len(matrix.commands),
            "alphaStep": float(alpha_step),
            "zStep": float(z_step),
        },
        "targets": {
            "cells": [_cell_dict(cell) for cell in targets],
            "targetCellCount": target_count,
            "requireFullTargetReachability": bool(require_full_target_reachability),
        },
        "response": {
            "alphaRank": matrix.alpha_rank,
            "heightRank": matrix.height_rank,
            "reachableAlphaCells": matrix.reachable_alpha_cells(tolerance),
            "reachableHeightCells": matrix.reachable_height_cells(tolerance),
            "alphaUnderactuatedCells": matrix.alpha_underactuated_cells(tolerance),
            "heightUnderactuatedCells": matrix.height_underactuated_cells(tolerance),
            "targetAlphaReachableCells": alpha_target_reached,
            "targetHeightReachableCells": height_target_reached,
            "targetAlphaUnderactuatedCells": alpha_under,
            "targetHeightUnderactuatedCells": height_under,
            "tolerance": float(tolerance),
            "alphaReachableMap": alpha_map.tolist(),
            "heightReachableMap": height_map.tolist(),
        },
        "topology": {
            "componentCount": int(topology["component_count"]),
            "deletedEdgeCount": int(topology["deleted_edge_count"]),
            "removedCellCount": int(topology["removed_cell_count"]),
            "topologyReachableTargetCells": topology_target_reached,
            "topologyBlockedTargetCells": topology_blocked,
            "componentLabels": np.asarray(topology["component_labels"], dtype=int).tolist(),
            "topologyReachableMap": topology_map.tolist(),
        },
        "equilibrium": {
            "attached": isinstance(resolved_equilibrium, dict),
            "schema": resolved_equilibrium.get("schema", "")
            if isinstance(resolved_equilibrium, dict)
            else "",
            "equilibriumRelationReady": bool(
                equilibrium_summary.get("equilibriumRelationReady")
            ),
            "status": equilibrium_summary.get("status", ""),
        },
        "summary": {
            "status": (
                "reachable-equilibrium-controllability-ready"
                if ready
                else "needs-reachability-controllability-evidence"
            ),
            "reachableEquilibriumControllabilityReady": ready,
            "targetFullyReachable": bool(
                target_count > 0
                and alpha_under == 0
                and height_under == 0
                and topology_blocked == 0
            ),
            "targetUnderactuatedCells": max(alpha_under, height_under, topology_blocked),
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_controllability_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumControllabilityNat",
            "leanPredicate": "reachableEquilibriumControllabilityReadyNat",
            "schema": "rad-sim.reachable-equilibrium-controllability.v1",
        },
        "claimLabels": {
            "reachability": "finite response-matrix simulator diagnostic",
            "equilibrium": "Lean-proven finite equilibrium evidence gate",
            "topology": "graph deletion and component bookkeeping",
            "controllability": "linearized finite actuator-basis proxy",
            "physicalAccuracy": "experimentally unvalidated physical assumption",
        },
        "limitations": [
            "Reachability is measured through finite response columns, not full nonlinear controllability.",
            "Topology reachability is graph-component bookkeeping, not a proof of mechanical motion feasibility.",
            "Full target reachability is optional because underactuation is itself useful design information.",
        ],
    }


def export_reachable_equilibrium_controllability_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    actuator_cells: Iterable[tuple[int, int]] | None = None,
    target_cells: Iterable[tuple[int, int]] | None = None,
    equilibrium_report: dict[str, object] | None = None,
    load_case: LoadCase | None = None,
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    include_alpha: bool = True,
    include_z: bool = True,
    require_full_target_reachability: bool = False,
    tolerance: float = 1e-9,
    residual_tolerance: float | None = None,
) -> str:
    return json.dumps(
        reachable_equilibrium_controllability_report(
            config,
            state,
            event_sequence=event_sequence,
            actuator_cells=actuator_cells,
            target_cells=target_cells,
            equilibrium_report=equilibrium_report,
            load_case=load_case,
            alpha_step=alpha_step,
            z_step=z_step,
            include_alpha=include_alpha,
            include_z=include_z,
            require_full_target_reachability=require_full_target_reachability,
            tolerance=tolerance,
            residual_tolerance=residual_tolerance,
        ),
        indent=2,
    )


def export_reachable_equilibrium_controllability_csv(report: dict[str, object]) -> str:
    response = report.get("response", {})
    topology = report.get("topology", {})
    basis = report.get("actuatorBasis", {})
    targets = report.get("targets", {})
    summary = report.get("summary", {})
    if not isinstance(response, dict):
        response = {}
    if not isinstance(topology, dict):
        topology = {}
    if not isinstance(basis, dict):
        basis = {}
    if not isinstance(targets, dict):
        targets = {}
    if not isinstance(summary, dict):
        summary = {}
    header = [
        "schema",
        "status",
        "reachable_equilibrium_controllability_ready",
        "target_fully_reachable",
        "actuator_basis_cells",
        "command_columns",
        "target_cells",
        "alpha_rank",
        "height_rank",
        "reachable_alpha_cells",
        "reachable_height_cells",
        "target_alpha_underactuated_cells",
        "target_height_underactuated_cells",
        "topology_blocked_target_cells",
        "component_count",
        "deleted_edges",
        "missing_evidence",
    ]
    row = [
        report.get("schema", ""),
        summary.get("status", ""),
        summary.get("reachableEquilibriumControllabilityReady", ""),
        summary.get("targetFullyReachable", ""),
        basis.get("cellCount", ""),
        basis.get("commandColumnCount", ""),
        targets.get("targetCellCount", ""),
        response.get("alphaRank", ""),
        response.get("heightRank", ""),
        response.get("reachableAlphaCells", ""),
        response.get("reachableHeightCells", ""),
        response.get("targetAlphaUnderactuatedCells", ""),
        response.get("targetHeightUnderactuatedCells", ""),
        topology.get("topologyBlockedTargetCells", ""),
        topology.get("componentCount", ""),
        topology.get("deletedEdgeCount", ""),
        ";".join(str(item) for item in summary.get("missingEvidence", []))
        if isinstance(summary.get("missingEvidence"), list)
        else "",
    ]
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in values)
        for values in (header, row)
    )


def _report_cell_tuple(raw: object) -> tuple[int, int] | None:
    if isinstance(raw, dict):
        if "row" not in raw or "col" not in raw:
            return None
        return (int(raw["row"]), int(raw["col"]))
    if isinstance(raw, (tuple, list)) and len(raw) >= 2:
        return (int(raw[0]), int(raw[1]))
    return None


def _report_cells(raw: object) -> tuple[tuple[int, int], ...]:
    if not isinstance(raw, list):
        return ()
    cells: list[tuple[int, int]] = []
    for item in raw:
        cell = _report_cell_tuple(item)
        if cell is not None:
            cells.append(cell)
    return _unique_cells(cells)


def _bool_map_get(raw: object, cell: tuple[int, int]) -> bool:
    row, col = cell
    try:
        return bool(raw[row][col])  # type: ignore[index]
    except (IndexError, TypeError):
        return False


def _int_map_get(raw: object, cell: tuple[int, int], default: int = -1) -> int:
    row, col = cell
    try:
        return int(raw[row][col])  # type: ignore[index]
    except (IndexError, TypeError, ValueError):
        return default


def _reachable_equilibrium_target_summaries(
    report: dict[str, object],
) -> list[dict[str, object]]:
    targets_raw = report.get("targets", {})
    response = report.get("response", {})
    topology = report.get("topology", {})
    targets = _report_cells(
        targets_raw.get("cells", []) if isinstance(targets_raw, dict) else []
    )
    alpha_map = response.get("alphaReachableMap", []) if isinstance(response, dict) else []
    height_map = response.get("heightReachableMap", []) if isinstance(response, dict) else []
    topology_map = (
        topology.get("topologyReachableMap", []) if isinstance(topology, dict) else []
    )
    component_labels = (
        topology.get("componentLabels", []) if isinstance(topology, dict) else []
    )
    summaries: list[dict[str, object]] = []
    for cell in targets:
        alpha_reachable = _bool_map_get(alpha_map, cell)
        height_reachable = _bool_map_get(height_map, cell)
        topology_reachable = _bool_map_get(topology_map, cell)
        summaries.append(
            {
                **_cell_dict(cell),
                "alphaReachable": alpha_reachable,
                "heightReachable": height_reachable,
                "topologyReachable": topology_reachable,
                "componentLabel": _int_map_get(component_labels, cell),
                "blockedByTopology": not topology_reachable,
            }
        )
    return summaries


def _reachable_equilibrium_bench_step(
    *,
    step_id: str,
    scope: str,
    commands: Iterable[SourceCommand],
    observation_cells: Iterable[tuple[int, int]],
    target_cells: Iterable[tuple[int, int]],
    purpose: str,
    expected_response: str,
    pass_fail: str,
    repeat_count: int,
    claim_label: str,
) -> dict[str, object]:
    return {
        "id": step_id,
        "scope": scope,
        "commands": [
            {"row": row, "col": col, "alpha": command.alpha, "z": command.z}
            for command in commands
            for row, col in (command.cell,)
        ],
        "observationCells": [_cell_dict(cell) for cell in _unique_cells(observation_cells)],
        "targetCells": [_cell_dict(cell) for cell in _unique_cells(target_cells)],
        "requiredMeasurements": {
            "fields": list(REACHABLE_EQUILIBRIUM_BENCH_FIELDS),
            "coordinateFrame": "registered cell-center grid with pixel/mm or probe/mm scale",
        },
        "purpose": purpose,
        "expectedResponse": expected_response,
        "passFailCriterion": pass_fail,
        "repeatCount": int(repeat_count),
        "claimLabel": claim_label,
    }


def reachable_equilibrium_bench_protocol(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    controllability_report: dict[str, object] | None = None,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    actuator_cells: Iterable[tuple[int, int]] | None = None,
    target_cells: Iterable[tuple[int, int]] | None = None,
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    repeat_count: int = 3,
    residual_tolerance: float | None = None,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Build bench trials from a reachable-equilibrium controllability report."""

    report = controllability_report or reachable_equilibrium_controllability_report(
        config,
        state,
        event_sequence=event_sequence,
        actuator_cells=actuator_cells,
        target_cells=target_cells,
        alpha_step=alpha_step,
        z_step=z_step,
        residual_tolerance=residual_tolerance,
        tolerance=tolerance,
    )
    basis_raw = report.get("actuatorBasis", {})
    targets_raw = report.get("targets", {})
    summary_raw = report.get("summary", {})
    topology_raw = report.get("topology", {})
    basis = _report_cells(
        basis_raw.get("cells", []) if isinstance(basis_raw, dict) else []
    )
    targets = _report_cells(
        targets_raw.get("cells", []) if isinstance(targets_raw, dict) else []
    )
    target_summaries = _reachable_equilibrium_target_summaries(report)
    blocked_targets = tuple(
        (int(item["row"]), int(item["col"]))
        for item in target_summaries
        if item.get("blockedByTopology")
    )
    observation_cells = _unique_cells((*basis, *targets))
    steps: list[dict[str, object]] = [
        _reachable_equilibrium_bench_step(
            step_id="baseline_equilibrium_capture",
            scope="baseline",
            commands=(),
            observation_cells=observation_cells,
            target_cells=targets,
            purpose=(
                "Record the undeformed or post-event equilibrium state before "
                "any actuator-column perturbation."
            ),
            expected_response=(
                "Measured alpha, height, lock, contact, and component labels "
                "establish the zero command reference."
            ),
            pass_fail=(
                "baseline equilibrium residual is at or below the configured "
                "measurement tolerance before response columns are collected"
            ),
            repeat_count=repeat_count,
            claim_label="required boundary-condition and equilibrium evidence",
        )
    ]
    for index, cell in enumerate(basis):
        steps.append(
            _reachable_equilibrium_bench_step(
                step_id=f"basis_{index}_alpha_response_column",
                scope="actuator-column",
                commands=(SourceCommand(cell, alpha=float(alpha_step)),),
                observation_cells=observation_cells,
                target_cells=targets,
                purpose=(
                    "Measure the physical alpha response column for one "
                    "candidate actuator-basis cell."
                ),
                expected_response=(
                    "Target alpha deltas should match the simulator-predicted "
                    "reachable map within calibrated tolerance after backlash."
                ),
                pass_fail=(
                    "alpha target deltas, topology labels, and lock states are "
                    "recorded for every target cell"
                ),
                repeat_count=repeat_count,
                claim_label="simulator-derived empirical law until bench data exists",
            )
        )
        steps.append(
            _reachable_equilibrium_bench_step(
                step_id=f"basis_{index}_z_response_column",
                scope="actuator-column",
                commands=(SourceCommand(cell, z=float(z_step)),),
                observation_cells=observation_cells,
                target_cells=targets,
                purpose=(
                    "Measure the physical vertical response column and residual "
                    "neighbor movement for one actuator-basis cell."
                ),
                expected_response=(
                    "Height deltas should show direct actuation plus die-off "
                    "through pin-hole clearance and graph connectivity."
                ),
                pass_fail=(
                    "height target deltas and pin-hole slip measurements are "
                    "recorded for every target cell"
                ),
                repeat_count=repeat_count,
                claim_label="experimentally unvalidated physical assumption",
            )
        )
    if basis:
        group_commands = tuple(
            SourceCommand(cell, alpha=float(alpha_step), z=float(z_step))
            for cell in basis
        )
        steps.append(
            _reachable_equilibrium_bench_step(
                step_id="group_target_reachability_probe",
                scope="group",
                commands=group_commands,
                observation_cells=observation_cells,
                target_cells=targets,
                purpose=(
                    "Probe whether simultaneous basis actuation reaches the "
                    "target set predicted by finite response-column evidence."
                ),
                expected_response=(
                    "The measured group field should be compared against the "
                    "linear response-matrix span and a sequenced application."
                ),
                pass_fail=(
                    "group response is saved with simultaneous/sequenced order "
                    "metadata and target residuals"
                ),
                repeat_count=repeat_count,
                claim_label="finite response-matrix controllability proxy",
            )
        )
    if blocked_targets:
        steps.append(
            _reachable_equilibrium_bench_step(
                step_id="topology_blocked_target_control",
                scope="topology-control",
                commands=tuple(SourceCommand(cell, z=float(z_step)) for cell in basis),
                observation_cells=blocked_targets,
                target_cells=blocked_targets,
                purpose=(
                    "Check whether removed-cell graph deletion really blocks "
                    "the target component under physical actuation."
                ),
                expected_response=(
                    "Topology-blocked targets should remain below calibrated "
                    "noise unless the real sheet has bypass coupling."
                ),
                pass_fail=(
                    "blocked targets are explicitly measured and any nonzero "
                    "response is labeled as model-mismatch evidence"
                ),
                repeat_count=repeat_count,
                claim_label="topology/connectivity physical validation trial",
            )
        )
    pass_fail_criteria = {
        "baseline": "equilibrium residual and fixture registration are recorded before perturbation",
        "responseColumns": "each actuator-basis alpha and z response column has target measurements",
        "topologyControls": "topology-blocked targets are measured whenever the report predicts a block",
        "claimLimit": "passing the protocol validates a bench observation table, not nonlinear controllability",
    }
    topology_policy_satisfied = len(blocked_targets) == 0 or any(
        step["id"] == "topology_blocked_target_control" for step in steps
    )
    missing_evidence: list[str] = []
    if not (
        isinstance(summary_raw, dict)
        and summary_raw.get("reachableEquilibriumControllabilityReady")
    ):
        missing_evidence.append("reachableEquilibriumControllability")
    if not basis:
        missing_evidence.append("actuatorBasis")
    if not targets:
        missing_evidence.append("targetCells")
    if not steps:
        missing_evidence.append("protocolSteps")
    if not REACHABLE_EQUILIBRIUM_BENCH_FIELDS:
        missing_evidence.append("measurementFields")
    if not pass_fail_criteria:
        missing_evidence.append("passFailCriteria")
    if not topology_policy_satisfied:
        missing_evidence.append("topologyBlockedControl")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.reachable-equilibrium-bench-protocol.v1",
        "method": (
            "bench protocol derived from finite reachable-equilibrium "
            "controllability evidence"
        ),
        "grid": {"rows": config.rows, "cols": config.cols},
        "sourceReport": {
            "schema": report.get("schema", ""),
            "status": summary_raw.get("status", "") if isinstance(summary_raw, dict) else "",
            "targetFullyReachable": bool(
                summary_raw.get("targetFullyReachable")
                if isinstance(summary_raw, dict)
                else False
            ),
            "targetUnderactuatedCells": int(
                summary_raw.get("targetUnderactuatedCells", 0)
                if isinstance(summary_raw, dict)
                else 0
            ),
        },
        "actuatorBasis": [_cell_dict(cell) for cell in basis],
        "targetSummaries": target_summaries,
        "topology": {
            "componentCount": int(topology_raw.get("componentCount", 0))
            if isinstance(topology_raw, dict)
            else 0,
            "topologyBlockedTargetCells": len(blocked_targets),
            "topologyPolicySatisfied": topology_policy_satisfied,
        },
        "measurementFields": list(REACHABLE_EQUILIBRIUM_BENCH_FIELDS),
        "repeatCount": int(repeat_count),
        "steps": steps,
        "passFailCriteria": pass_fail_criteria,
        "summary": {
            "status": (
                "reachable-equilibrium-bench-protocol-ready"
                if ready
                else "needs-reachable-equilibrium-bench-evidence"
            ),
            "benchProtocolReady": ready,
            "stepCount": len(steps),
            "actuatorColumnTrialCount": max(0, 2 * len(basis)),
            "targetObservationCellCount": len(targets),
            "topologyBlockedControlCount": 1 if blocked_targets else 0,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_bench_protocol_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumBenchProtocolNat",
            "leanPredicate": "reachableEquilibriumBenchProtocolReadyNat",
            "schema": "rad-sim.reachable-equilibrium-bench-protocol.v1",
        },
        "claimLabels": {
            "protocol": "bench protocol artifact",
            "reachability": "finite simulator evidence to be tested physically",
            "blockedTargets": "graph-topology hypothesis requiring measurement",
            "physicalAccuracy": "experimentally unvalidated physical assumption",
        },
        "limitations": [
            "The protocol defines measurements but does not certify they were collected.",
            "No physical tolerance is calibrated until hardware scale, backlash, force, and clearance are measured.",
            "A passing bench protocol still does not prove continuous nonlinear controllability.",
        ],
    }


def export_reachable_equilibrium_bench_protocol_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    *,
    controllability_report: dict[str, object] | None = None,
    event_sequence: Iterable[ProgrammableDiscontinuityEvent] | None = None,
    actuator_cells: Iterable[tuple[int, int]] | None = None,
    target_cells: Iterable[tuple[int, int]] | None = None,
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    repeat_count: int = 3,
    residual_tolerance: float | None = None,
    tolerance: float = 1e-9,
) -> str:
    return json.dumps(
        reachable_equilibrium_bench_protocol(
            config,
            state,
            controllability_report=controllability_report,
            event_sequence=event_sequence,
            actuator_cells=actuator_cells,
            target_cells=target_cells,
            alpha_step=alpha_step,
            z_step=z_step,
            repeat_count=repeat_count,
            residual_tolerance=residual_tolerance,
            tolerance=tolerance,
        ),
        indent=2,
    )


def export_reachable_equilibrium_bench_protocol_csv(
    protocol: dict[str, object],
) -> str:
    summary = protocol.get("summary", {})
    topology = protocol.get("topology", {})
    if not isinstance(summary, dict):
        summary = {}
    if not isinstance(topology, dict):
        topology = {}
    header = [
        "schema",
        "protocol_ready",
        "step_id",
        "scope",
        "command_count",
        "observation_cells",
        "target_cells",
        "topology_blocked_target_cells",
        "repeat_count",
        "expected_response",
        "claim_label",
        "missing_evidence",
    ]
    rows: list[list[object]] = [header]
    for step in protocol.get("steps", []):
        if not isinstance(step, dict):
            continue
        rows.append(
            [
                protocol.get("schema", ""),
                summary.get("benchProtocolReady", ""),
                step.get("id", ""),
                step.get("scope", ""),
                len(step.get("commands", []))
                if isinstance(step.get("commands"), list)
                else "",
                json.dumps(step.get("observationCells", []), separators=(",", ":")),
                json.dumps(step.get("targetCells", []), separators=(",", ":")),
                topology.get("topologyBlockedTargetCells", ""),
                step.get("repeatCount", ""),
                step.get("expectedResponse", ""),
                step.get("claimLabel", ""),
                ";".join(str(item) for item in summary.get("missingEvidence", []))
                if isinstance(summary.get("missingEvidence"), list)
                else "",
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row) for row in rows
    )


def _bench_protocol_target_summary_map(
    protocol: dict[str, object],
) -> dict[tuple[int, int], dict[str, object]]:
    out: dict[tuple[int, int], dict[str, object]] = {}
    for raw in protocol.get("targetSummaries", []):
        if not isinstance(raw, dict):
            continue
        cell = _report_cell_tuple(raw)
        if cell is not None:
            out[cell] = raw
    return out


def _bench_step_measurement_mode(step: dict[str, object]) -> tuple[str, ...]:
    scope = str(step.get("scope", ""))
    step_id = str(step.get("id", ""))
    if scope == "baseline":
        return ("baseline",)
    if scope == "topology-control":
        return ("topology",)
    if scope == "group":
        return ("alpha", "height", "topology")
    if step_id.endswith("_alpha_response_column"):
        return ("alpha",)
    if step_id.endswith("_z_response_column"):
        return ("height",)
    return ("alpha", "height")


def reachable_equilibrium_bench_results_template(
    protocol: dict[str, object],
    *,
    dataset_id: str = "reachable-equilibrium-run-001",
    dataset_role: str = "bench",
    source_file_id: str | None = None,
    operator: str = "",
    collected_at: str = "",
) -> dict[str, object]:
    """Create a fillable results table for a reachable-equilibrium protocol."""

    target_summary = _bench_protocol_target_summary_map(protocol)
    rows: list[dict[str, object]] = []
    for step in protocol.get("steps", []):
        if not isinstance(step, dict):
            continue
        repeat_count = max(1, int(step.get("repeatCount", protocol.get("repeatCount", 1)) or 1))
        application_modes = (
            ("simultaneous", "sequenced")
            if str(step.get("scope", "")) == "group"
            else (str(step.get("scope", "single")),)
        )
        target_cells = _report_cells(step.get("targetCells", []))
        for repeat_index in range(repeat_count):
            for application_mode in application_modes:
                rows.append(
                    {
                        "stepId": step.get("id", ""),
                        "scope": step.get("scope", ""),
                        "repeatIndex": repeat_index,
                        "applicationMode": application_mode,
                        "commands": step.get("commands", []),
                        "expectedResponse": step.get("expectedResponse", ""),
                        "measurementMode": list(_bench_step_measurement_mode(step)),
                        "targetMeasurements": [
                            {
                                "row": row,
                                "col": col,
                                "predictedAlphaReachable": bool(
                                    target_summary.get((row, col), {}).get(
                                        "alphaReachable", False
                                    )
                                ),
                                "predictedHeightReachable": bool(
                                    target_summary.get((row, col), {}).get(
                                        "heightReachable", False
                                    )
                                ),
                                "predictedTopologyBlocked": bool(
                                    target_summary.get((row, col), {}).get(
                                        "blockedByTopology", False
                                    )
                                ),
                                "expectedTopologyComponentLabel": int(
                                    target_summary.get((row, col), {}).get(
                                        "componentLabel", -1
                                    )
                                ),
                                "measuredAlphaDelta": None,
                                "measuredHeightDelta": None,
                                "measuredTopologyComponentLabel": None,
                                "measuredPinHoleSlipMm": None,
                                "measuredActuatorForceN": None,
                                "notes": "",
                            }
                            for row, col in target_cells
                        ],
                    }
                )
    return {
        "schema": "rad-sim.reachable-equilibrium-bench-results.v1",
        "protocolSchema": protocol.get("schema", ""),
        "datasetId": dataset_id,
        "datasetRole": dataset_role,
        "sourceFileId": source_file_id or f"{dataset_id}.json",
        "operator": operator,
        "collectedAt": collected_at,
        "measurementFields": list(REACHABLE_EQUILIBRIUM_BENCH_FIELDS),
        "provenanceNotes": (
            "Fill measuredAlphaDelta and measuredHeightDelta from the registered "
            "cell-center frame; leave predictions unchanged."
        ),
        "measurements": rows,
    }


def export_reachable_equilibrium_bench_results_template_json(
    protocol: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        reachable_equilibrium_bench_results_template(protocol, **kwargs),
        indent=2,
    )


def reachable_equilibrium_bench_results_from_json(text: str) -> dict[str, object]:
    raw = json.loads(text)
    if not isinstance(raw, dict):
        raise ValueError("reachable equilibrium bench results must be a JSON object")
    if raw.get("schema") != "rad-sim.reachable-equilibrium-bench-results.v1":
        raise ValueError("unsupported reachable equilibrium bench results schema")
    if not isinstance(raw.get("measurements"), list):
        raise ValueError("reachable equilibrium bench results require measurements")
    return raw


def _bench_float(raw: object) -> float | None:
    try:
        value = float(raw)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if not np.isfinite(value):
        return None
    return value


def _bench_target_key(raw: object) -> tuple[int, int] | None:
    return _report_cell_tuple(raw)


def _bench_results_by_step_mode(
    results: dict[str, object],
) -> dict[tuple[str, int, str], dict[str, object]]:
    rows: dict[tuple[str, int, str], dict[str, object]] = {}
    for raw in results.get("measurements", []):
        if not isinstance(raw, dict):
            continue
        key = (
            str(raw.get("stepId", "")),
            int(raw.get("repeatIndex", 0) or 0),
            str(raw.get("applicationMode", "")),
        )
        rows[key] = raw
    return rows


def _bench_target_measurement_map(
    measurement_row: dict[str, object],
) -> dict[tuple[int, int], dict[str, object]]:
    out: dict[tuple[int, int], dict[str, object]] = {}
    for raw in measurement_row.get("targetMeasurements", []):
        if not isinstance(raw, dict):
            continue
        key = _bench_target_key(raw)
        if key is not None:
            out[key] = raw
    return out


def compare_reachable_equilibrium_bench_results(
    protocol: dict[str, object],
    results: dict[str, object],
    *,
    response_tolerance: float = 1e-6,
    blocked_tolerance: float = 1e-6,
    group_sequence_tolerance: float = 1e-6,
) -> dict[str, object]:
    """Compare filled reachable-equilibrium bench rows against protocol flags."""

    if results.get("schema") != "rad-sim.reachable-equilibrium-bench-results.v1":
        raise ValueError("unsupported reachable equilibrium bench results schema")
    step_by_id = {
        str(step.get("id", "")): step
        for step in protocol.get("steps", [])
        if isinstance(step, dict)
    }
    result_rows = [
        row for row in results.get("measurements", []) if isinstance(row, dict)
    ]
    missing_measurements = 0
    target_measurements = 0
    alpha_checks = 0
    height_checks = 0
    alpha_mismatches = 0
    height_mismatches = 0
    topology_leakage = 0
    component_mismatches = 0
    max_abs_alpha = 0.0
    max_abs_height = 0.0
    row_reports: list[dict[str, object]] = []
    for row in result_rows:
        step = step_by_id.get(str(row.get("stepId", "")), {})
        modes = set(
            str(item)
            for item in (
                row.get("measurementMode")
                if isinstance(row.get("measurementMode"), list)
                else _bench_step_measurement_mode(step)
            )
        )
        row_missing = 0
        row_alpha_mismatch = 0
        row_height_mismatch = 0
        row_topology_leakage = 0
        row_component_mismatch = 0
        for target in row.get("targetMeasurements", []):
            if not isinstance(target, dict):
                continue
            target_measurements += 1
            alpha = _bench_float(target.get("measuredAlphaDelta"))
            height = _bench_float(target.get("measuredHeightDelta"))
            if alpha is None:
                row_missing += 1
            else:
                max_abs_alpha = max(max_abs_alpha, abs(alpha))
            if height is None:
                row_missing += 1
            else:
                max_abs_height = max(max_abs_height, abs(height))
            if "baseline" in modes:
                if alpha is not None:
                    alpha_checks += 1
                    if abs(alpha) > response_tolerance:
                        row_alpha_mismatch += 1
                if height is not None:
                    height_checks += 1
                    if abs(height) > response_tolerance:
                        row_height_mismatch += 1
            if "alpha" in modes and alpha is not None:
                alpha_checks += 1
                expected = bool(target.get("predictedAlphaReachable"))
                observed = abs(alpha) > response_tolerance
                if expected != observed:
                    row_alpha_mismatch += 1
            if "height" in modes and height is not None:
                height_checks += 1
                expected = bool(target.get("predictedHeightReachable"))
                observed = abs(height) > response_tolerance
                if expected != observed:
                    row_height_mismatch += 1
            if "topology" in modes and bool(target.get("predictedTopologyBlocked")):
                leaked = (
                    (alpha is not None and abs(alpha) > blocked_tolerance)
                    or (height is not None and abs(height) > blocked_tolerance)
                )
                if leaked:
                    row_topology_leakage += 1
            measured_component = target.get("measuredTopologyComponentLabel")
            if measured_component is not None:
                try:
                    measured_label = int(measured_component)
                    expected_label = int(
                        target.get("expectedTopologyComponentLabel", measured_label)
                    )
                    if measured_label != expected_label:
                        row_component_mismatch += 1
                except (TypeError, ValueError):
                    row_component_mismatch += 1
        missing_measurements += row_missing
        alpha_mismatches += row_alpha_mismatch
        height_mismatches += row_height_mismatch
        topology_leakage += row_topology_leakage
        component_mismatches += row_component_mismatch
        row_reports.append(
            {
                "stepId": row.get("stepId", ""),
                "repeatIndex": row.get("repeatIndex", 0),
                "applicationMode": row.get("applicationMode", ""),
                "missingMeasurementCount": row_missing,
                "alphaReachabilityMismatchCount": row_alpha_mismatch,
                "heightReachabilityMismatchCount": row_height_mismatch,
                "topologyLeakageCount": row_topology_leakage,
                "componentMismatchCount": row_component_mismatch,
            }
        )
    rows_by_mode = _bench_results_by_step_mode(results)
    group_sequence_pairs = 0
    group_sequence_max_error = 0.0
    for step in protocol.get("steps", []):
        if not isinstance(step, dict) or str(step.get("scope", "")) != "group":
            continue
        repeat_count = max(1, int(step.get("repeatCount", protocol.get("repeatCount", 1)) or 1))
        for repeat_index in range(repeat_count):
            simultaneous = rows_by_mode.get(
                (str(step.get("id", "")), repeat_index, "simultaneous")
            )
            sequenced = rows_by_mode.get(
                (str(step.get("id", "")), repeat_index, "sequenced")
            )
            if simultaneous is None or sequenced is None:
                continue
            sim_targets = _bench_target_measurement_map(simultaneous)
            seq_targets = _bench_target_measurement_map(sequenced)
            for cell, sim_target in sim_targets.items():
                seq_target = seq_targets.get(cell)
                if seq_target is None:
                    continue
                group_sequence_pairs += 1
                for key in ("measuredAlphaDelta", "measuredHeightDelta"):
                    sim_value = _bench_float(sim_target.get(key))
                    seq_value = _bench_float(seq_target.get(key))
                    if sim_value is None or seq_value is None:
                        continue
                    group_sequence_max_error = max(
                        group_sequence_max_error,
                        abs(sim_value - seq_value),
                    )
    group_sequence_pass = group_sequence_max_error <= group_sequence_tolerance
    protocol_summary = protocol.get("summary", {})
    protocol_ready = bool(
        protocol_summary.get("benchProtocolReady")
        if isinstance(protocol_summary, dict)
        else False
    )
    missing_evidence: list[str] = []
    if not protocol_ready:
        missing_evidence.append("reachableEquilibriumBenchProtocol")
    if not result_rows:
        missing_evidence.append("resultRows")
    if target_measurements <= 0:
        missing_evidence.append("targetMeasurements")
    if missing_measurements > 0:
        missing_evidence.append("completedMeasurements")
    if alpha_checks + height_checks <= 0:
        missing_evidence.append("reachabilityChecks")
    if group_sequence_pairs <= 0:
        missing_evidence.append("groupSequenceChecks")
    pass_flag = (
        len(missing_evidence) == 0
        and alpha_mismatches == 0
        and height_mismatches == 0
        and topology_leakage == 0
        and component_mismatches == 0
        and group_sequence_pass
    )
    return {
        "schema": "rad-sim.reachable-equilibrium-bench-comparison.v1",
        "method": (
            "comparison of filled reachable-equilibrium bench measurements "
            "against simulator reachability flags and topology-blocked controls"
        ),
        "sourceProtocol": {
            "schema": protocol.get("schema", ""),
            "ready": protocol_ready,
        },
        "sourceResults": {
            "schema": results.get("schema", ""),
            "datasetId": results.get("datasetId", ""),
            "datasetRole": results.get("datasetRole", ""),
            "sourceFileId": results.get("sourceFileId", ""),
        },
        "metrics": {
            "resultRowCount": len(result_rows),
            "targetMeasurementCount": target_measurements,
            "missingMeasurementCount": missing_measurements,
            "alphaCheckCount": alpha_checks,
            "heightCheckCount": height_checks,
            "alphaReachabilityMismatchCount": alpha_mismatches,
            "heightReachabilityMismatchCount": height_mismatches,
            "topologyLeakageCount": topology_leakage,
            "componentMismatchCount": component_mismatches,
            "groupSequencePairCount": group_sequence_pairs,
            "groupSequenceMaxError": group_sequence_max_error,
            "maxAbsMeasuredAlphaDelta": max_abs_alpha,
            "maxAbsMeasuredHeightDelta": max_abs_height,
            "responseTolerance": float(response_tolerance),
            "blockedTolerance": float(blocked_tolerance),
            "groupSequenceTolerance": float(group_sequence_tolerance),
        },
        "summary": {
            "status": (
                "reachable-equilibrium-bench-comparison-pass"
                if pass_flag
                else "needs-reachable-equilibrium-bench-review"
            ),
            "benchComparisonPass": pass_flag,
            "topologyBlockedLeakagePass": topology_leakage == 0,
            "groupSequencePass": group_sequence_pass and group_sequence_pairs > 0,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "rows": row_reports,
        "formalization": {
            "targetId": "reachable_equilibrium_bench_validation_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumBenchValidationNat",
            "leanPredicate": "reachableEquilibriumBenchValidationReadyNat",
            "schema": "rad-sim.reachable-equilibrium-bench-comparison.v1",
        },
        "claimLabels": {
            "comparison": "bench result comparison artifact",
            "reachability": "simulator-derived empirical law compared to measurement",
            "topologyLeakage": "experimentally measured model-mismatch indicator",
            "physicalAccuracy": "requires calibrated hardware and repeated trials",
        },
        "limitations": [
            "Reachable flags are qualitative nonzero-response checks, not fitted amplitude laws.",
            "Topology leakage can indicate unmodeled bypass coupling, fixture compliance, or measurement noise.",
            "Passing this comparison does not prove nonlinear controllability or contact mechanics.",
        ],
    }


def export_reachable_equilibrium_bench_comparison_json(
    protocol: dict[str, object],
    results: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        compare_reachable_equilibrium_bench_results(protocol, results, **kwargs),
        indent=2,
    )


def export_reachable_equilibrium_bench_comparison_csv(
    report: dict[str, object],
) -> str:
    metrics = report.get("metrics", {})
    summary = report.get("summary", {})
    if not isinstance(metrics, dict):
        metrics = {}
    if not isinstance(summary, dict):
        summary = {}
    header = [
        "schema",
        "comparison_pass",
        "step_id",
        "repeat_index",
        "application_mode",
        "missing_measurements",
        "alpha_mismatches",
        "height_mismatches",
        "topology_leakage",
        "component_mismatches",
        "group_sequence_max_error",
        "missing_evidence",
    ]
    rows: list[list[object]] = [header]
    for row in report.get("rows", []):
        if not isinstance(row, dict):
            continue
        rows.append(
            [
                report.get("schema", ""),
                summary.get("benchComparisonPass", ""),
                row.get("stepId", ""),
                row.get("repeatIndex", ""),
                row.get("applicationMode", ""),
                row.get("missingMeasurementCount", ""),
                row.get("alphaReachabilityMismatchCount", ""),
                row.get("heightReachabilityMismatchCount", ""),
                row.get("topologyLeakageCount", ""),
                row.get("componentMismatchCount", ""),
                metrics.get("groupSequenceMaxError", ""),
                ";".join(str(item) for item in summary.get("missingEvidence", []))
                if isinstance(summary.get("missingEvidence"), list)
                else "",
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row) for row in rows
    )


def _bench_numeric_stats(
    values: Iterable[float],
    confidence_sigma: float,
) -> dict[str, object]:
    finite = [float(value) for value in values if np.isfinite(float(value))]
    if not finite:
        return {
            "sampleCount": 0,
            "mean": None,
            "sampleStd": None,
            "meanAbs": None,
            "maxAbs": None,
            "uncertaintyHalfWidth": None,
        }
    arr = np.asarray(finite, dtype=float)
    sample_std = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
    uncertainty = float(confidence_sigma) * sample_std / float(np.sqrt(arr.size))
    return {
        "sampleCount": int(arr.size),
        "mean": float(np.mean(arr)),
        "sampleStd": sample_std,
        "meanAbs": float(np.mean(np.abs(arr))),
        "maxAbs": float(np.max(np.abs(arr))),
        "uncertaintyHalfWidth": uncertainty,
    }


def _bench_qualitative_residual(
    mean_value: object,
    predicted_reachable: bool,
    tolerance: float,
) -> float:
    value = _bench_float(mean_value)
    if value is None:
        return float(tolerance)
    observed = abs(value) > tolerance
    if observed == predicted_reachable:
        return 0.0
    return abs(value) if observed else float(tolerance)


def _bench_estimate_key(raw: dict[str, object]) -> tuple[str, str, int, int]:
    return (
        str(raw.get("stepId", "")),
        str(raw.get("applicationMode", "")),
        int(raw.get("row", 0) or 0),
        int(raw.get("col", 0) or 0),
    )


def reachable_equilibrium_amplitude_calibration_report(
    protocol: dict[str, object],
    results: dict[str, object],
    *,
    comparison_report: dict[str, object] | None = None,
    response_tolerance: float = 1e-6,
    blocked_tolerance: float = 1e-6,
    group_sequence_tolerance: float = 1e-6,
    confidence_sigma: float = 2.0,
) -> dict[str, object]:
    """Estimate quantitative response amplitudes from repeated bench rows."""

    comparison = comparison_report or compare_reachable_equilibrium_bench_results(
        protocol,
        results,
        response_tolerance=response_tolerance,
        blocked_tolerance=blocked_tolerance,
        group_sequence_tolerance=group_sequence_tolerance,
    )
    grouped: dict[tuple[str, str, int, int], dict[str, object]] = {}
    step_by_id = {
        str(step.get("id", "")): step
        for step in protocol.get("steps", [])
        if isinstance(step, dict)
    }
    for row in results.get("measurements", []):
        if not isinstance(row, dict):
            continue
        modes = (
            tuple(str(item) for item in row.get("measurementMode", []))
            if isinstance(row.get("measurementMode"), list)
            else _bench_step_measurement_mode(step_by_id.get(str(row.get("stepId", "")), {}))
        )
        for target in row.get("targetMeasurements", []):
            if not isinstance(target, dict):
                continue
            cell = _bench_target_key(target)
            if cell is None:
                continue
            key = (
                str(row.get("stepId", "")),
                str(row.get("applicationMode", "")),
                cell[0],
                cell[1],
            )
            entry = grouped.setdefault(
                key,
                {
                    "stepId": str(row.get("stepId", "")),
                    "scope": str(row.get("scope", "")),
                    "applicationMode": str(row.get("applicationMode", "")),
                    "row": cell[0],
                    "col": cell[1],
                    "measurementMode": list(modes),
                    "predictedAlphaReachable": bool(
                        target.get("predictedAlphaReachable")
                    ),
                    "predictedHeightReachable": bool(
                        target.get("predictedHeightReachable")
                    ),
                    "predictedTopologyBlocked": bool(
                        target.get("predictedTopologyBlocked")
                    ),
                    "expectedTopologyComponentLabel": int(
                        target.get("expectedTopologyComponentLabel", -1) or -1
                    ),
                    "alphaValues": [],
                    "heightValues": [],
                    "slipValues": [],
                    "forceValues": [],
                },
            )
            alpha = _bench_float(target.get("measuredAlphaDelta"))
            height = _bench_float(target.get("measuredHeightDelta"))
            slip = _bench_float(target.get("measuredPinHoleSlipMm"))
            force = _bench_float(target.get("measuredActuatorForceN"))
            if alpha is not None:
                entry["alphaValues"].append(alpha)  # type: ignore[index, union-attr]
            if height is not None:
                entry["heightValues"].append(height)  # type: ignore[index, union-attr]
            if slip is not None:
                entry["slipValues"].append(slip)  # type: ignore[index, union-attr]
            if force is not None:
                entry["forceValues"].append(force)  # type: ignore[index, union-attr]
    estimates: list[dict[str, object]] = []
    residual_field: list[dict[str, object]] = []
    topology_bands: list[dict[str, object]] = []
    for raw in grouped.values():
        alpha_stats = _bench_numeric_stats(
            raw.get("alphaValues", []), confidence_sigma
        )
        height_stats = _bench_numeric_stats(
            raw.get("heightValues", []), confidence_sigma
        )
        slip_stats = _bench_numeric_stats(raw.get("slipValues", []), confidence_sigma)
        force_stats = _bench_numeric_stats(
            raw.get("forceValues", []), confidence_sigma
        )
        estimate = {
            "stepId": raw["stepId"],
            "scope": raw["scope"],
            "applicationMode": raw["applicationMode"],
            "row": raw["row"],
            "col": raw["col"],
            "measurementMode": raw["measurementMode"],
            "predictedAlphaReachable": raw["predictedAlphaReachable"],
            "predictedHeightReachable": raw["predictedHeightReachable"],
            "predictedTopologyBlocked": raw["predictedTopologyBlocked"],
            "expectedTopologyComponentLabel": raw["expectedTopologyComponentLabel"],
            "alpha": alpha_stats,
            "height": height_stats,
            "pinHoleSlipMm": slip_stats,
            "actuatorForceN": force_stats,
        }
        estimates.append(estimate)
        alpha_residual = _bench_qualitative_residual(
            alpha_stats["mean"],
            bool(raw["predictedAlphaReachable"]),
            response_tolerance,
        )
        height_residual = _bench_qualitative_residual(
            height_stats["mean"],
            bool(raw["predictedHeightReachable"]),
            response_tolerance,
        )
        residual_field.append(
            {
                "stepId": raw["stepId"],
                "applicationMode": raw["applicationMode"],
                "row": raw["row"],
                "col": raw["col"],
                "alphaQualitativeResidual": alpha_residual,
                "heightQualitativeResidual": height_residual,
                "combinedResidual": max(alpha_residual, height_residual),
            }
        )
        alpha_band = (
            None
            if alpha_stats["meanAbs"] is None
            else float(alpha_stats["meanAbs"])
            + float(alpha_stats["uncertaintyHalfWidth"] or 0.0)
        )
        height_band = (
            None
            if height_stats["meanAbs"] is None
            else float(height_stats["meanAbs"])
            + float(height_stats["uncertaintyHalfWidth"] or 0.0)
        )
        topology_bands.append(
            {
                "stepId": raw["stepId"],
                "applicationMode": raw["applicationMode"],
                "row": raw["row"],
                "col": raw["col"],
                "predictedTopologyBlocked": raw["predictedTopologyBlocked"],
                "alphaResponseBand": alpha_band,
                "heightResponseBand": height_band,
                "blockedTolerance": float(blocked_tolerance),
                "blockedLeakageBandPass": (
                    True
                    if not bool(raw["predictedTopologyBlocked"])
                    else (alpha_band or 0.0) <= blocked_tolerance
                    and (height_band or 0.0) <= blocked_tolerance
                ),
            }
        )
    rows_by_mode = _bench_results_by_step_mode(results)
    group_residuals: list[dict[str, object]] = []
    for step in protocol.get("steps", []):
        if not isinstance(step, dict) or str(step.get("scope", "")) != "group":
            continue
        repeat_count = max(1, int(step.get("repeatCount", protocol.get("repeatCount", 1)) or 1))
        for repeat_index in range(repeat_count):
            simultaneous = rows_by_mode.get(
                (str(step.get("id", "")), repeat_index, "simultaneous")
            )
            sequenced = rows_by_mode.get(
                (str(step.get("id", "")), repeat_index, "sequenced")
            )
            if simultaneous is None or sequenced is None:
                continue
            sim_targets = _bench_target_measurement_map(simultaneous)
            seq_targets = _bench_target_measurement_map(sequenced)
            for cell, sim_target in sim_targets.items():
                seq_target = seq_targets.get(cell)
                if seq_target is None:
                    continue
                alpha_error = 0.0
                height_error = 0.0
                sim_alpha = _bench_float(sim_target.get("measuredAlphaDelta"))
                seq_alpha = _bench_float(seq_target.get("measuredAlphaDelta"))
                sim_height = _bench_float(sim_target.get("measuredHeightDelta"))
                seq_height = _bench_float(seq_target.get("measuredHeightDelta"))
                if sim_alpha is not None and seq_alpha is not None:
                    alpha_error = abs(sim_alpha - seq_alpha)
                if sim_height is not None and seq_height is not None:
                    height_error = abs(sim_height - seq_height)
                group_residuals.append(
                    {
                        "stepId": step.get("id", ""),
                        "repeatIndex": repeat_index,
                        "row": cell[0],
                        "col": cell[1],
                        "alphaError": alpha_error,
                        "heightError": height_error,
                        "combinedError": max(alpha_error, height_error),
                    }
                )
    repeated_groups = sum(
        1
        for estimate in estimates
        if max(
            int(estimate["alpha"]["sampleCount"]),  # type: ignore[index]
            int(estimate["height"]["sampleCount"]),  # type: ignore[index]
        )
        >= 2
    )
    failed_topology_bands = sum(
        1 for band in topology_bands if not bool(band["blockedLeakageBandPass"])
    )
    max_residual = max(
        (float(item["combinedResidual"]) for item in residual_field),
        default=0.0,
    )
    max_group_error = max(
        (float(item["combinedError"]) for item in group_residuals),
        default=0.0,
    )
    comparison_summary = comparison.get("summary", {})
    comparison_pass = bool(
        comparison_summary.get("benchComparisonPass")
        if isinstance(comparison_summary, dict)
        else False
    )
    missing_evidence: list[str] = []
    if not comparison_pass:
        missing_evidence.append("benchComparison")
    if not estimates:
        missing_evidence.append("amplitudeEstimates")
    if repeated_groups <= 0:
        missing_evidence.append("repeatedTrials")
    if not residual_field:
        missing_evidence.append("residualField")
    if not topology_bands:
        missing_evidence.append("topologyBands")
    if failed_topology_bands > 0:
        missing_evidence.append("topologyLeakageBand")
    if not group_residuals:
        missing_evidence.append("groupSequenceResiduals")
    if max_group_error > group_sequence_tolerance:
        missing_evidence.append("groupSequenceTolerance")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
        "method": (
            "repeated-trial amplitude calibration report over filled "
            "reachable-equilibrium bench measurements"
        ),
        "sourceProtocol": {
            "schema": protocol.get("schema", ""),
            "stepCount": len(protocol.get("steps", []))
            if isinstance(protocol.get("steps"), list)
            else 0,
        },
        "sourceResults": {
            "schema": results.get("schema", ""),
            "datasetId": results.get("datasetId", ""),
            "datasetRole": results.get("datasetRole", ""),
            "sourceFileId": results.get("sourceFileId", ""),
        },
        "sourceComparison": {
            "schema": comparison.get("schema", ""),
            "benchComparisonPass": comparison_pass,
        },
        "parameters": {
            "responseTolerance": float(response_tolerance),
            "blockedTolerance": float(blocked_tolerance),
            "groupSequenceTolerance": float(group_sequence_tolerance),
            "confidenceSigma": float(confidence_sigma),
        },
        "amplitudeEstimates": estimates,
        "residualField": residual_field,
        "topologyResponseBands": topology_bands,
        "groupSequenceResiduals": group_residuals,
        "metrics": {
            "amplitudeEstimateCount": len(estimates),
            "repeatedTrialGroupCount": repeated_groups,
            "residualFieldCellCount": len(residual_field),
            "topologyBandCount": len(topology_bands),
            "failedTopologyBandCount": failed_topology_bands,
            "groupSequenceResidualCount": len(group_residuals),
            "maxQualitativeResidual": max_residual,
            "maxGroupSequenceResidual": max_group_error,
        },
        "summary": {
            "status": (
                "reachable-equilibrium-amplitude-calibration-ready"
                if ready
                else "needs-reachable-equilibrium-amplitude-review"
            ),
            "amplitudeCalibrationReady": ready,
            "topologyLeakageBandsPass": failed_topology_bands == 0,
            "groupSequenceResidualPass": max_group_error <= group_sequence_tolerance
            and len(group_residuals) > 0,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_amplitude_calibration_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumAmplitudeCalibrationNat",
            "leanPredicate": "reachableEquilibriumAmplitudeCalibrationReadyNat",
            "schema": "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
        },
        "claimLabels": {
            "amplitudes": "bench-measured empirical law",
            "uncertainty": "finite repeated-trial statistic",
            "topologyBands": "experimentally measured leakage bound",
            "physicalAccuracy": "requires external hardware calibration and uncertainty analysis",
        },
        "limitations": [
            "The report estimates measured amplitudes but does not derive a constitutive mechanics law.",
            "Uncertainty bands use a simple sigma multiplier and are not a full statistical model.",
            "Group residuals compare two command orderings but do not prove operator commutativity.",
        ],
    }


def export_reachable_equilibrium_amplitude_calibration_json(
    protocol: dict[str, object],
    results: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        reachable_equilibrium_amplitude_calibration_report(
            protocol,
            results,
            **kwargs,
        ),
        indent=2,
    )


def export_reachable_equilibrium_amplitude_calibration_csv(
    report: dict[str, object],
) -> str:
    summary = report.get("summary", {})
    metrics = report.get("metrics", {})
    if not isinstance(summary, dict):
        summary = {}
    if not isinstance(metrics, dict):
        metrics = {}
    header = [
        "schema",
        "calibration_ready",
        "step_id",
        "application_mode",
        "row",
        "col",
        "alpha_samples",
        "alpha_mean",
        "alpha_uncertainty",
        "height_samples",
        "height_mean",
        "height_uncertainty",
        "predicted_topology_blocked",
        "max_group_sequence_residual",
        "missing_evidence",
    ]
    rows: list[list[object]] = [header]
    for estimate in report.get("amplitudeEstimates", []):
        if not isinstance(estimate, dict):
            continue
        alpha = estimate.get("alpha", {})
        height = estimate.get("height", {})
        if not isinstance(alpha, dict):
            alpha = {}
        if not isinstance(height, dict):
            height = {}
        rows.append(
            [
                report.get("schema", ""),
                summary.get("amplitudeCalibrationReady", ""),
                estimate.get("stepId", ""),
                estimate.get("applicationMode", ""),
                estimate.get("row", ""),
                estimate.get("col", ""),
                alpha.get("sampleCount", ""),
                alpha.get("mean", ""),
                alpha.get("uncertaintyHalfWidth", ""),
                height.get("sampleCount", ""),
                height.get("mean", ""),
                height.get("uncertaintyHalfWidth", ""),
                estimate.get("predictedTopologyBlocked", ""),
                metrics.get("maxGroupSequenceResidual", ""),
                ";".join(str(item) for item in summary.get("missingEvidence", []))
                if isinstance(summary.get("missingEvidence"), list)
                else "",
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row) for row in rows
    )


def _amplitude_section_values(
    report: dict[str, object],
    channel: Literal["alpha", "height"],
) -> list[float]:
    values: list[float] = []
    reachable_key = (
        "predictedAlphaReachable" if channel == "alpha" else "predictedHeightReachable"
    )
    for estimate in report.get("amplitudeEstimates", []):
        if not isinstance(estimate, dict) or not bool(estimate.get(reachable_key)):
            continue
        section = estimate.get(channel, {})
        if not isinstance(section, dict):
            continue
        value = _finite_or_none(section.get("meanAbs"))
        if value is not None:
            values.append(value)
    return values


def _amplitude_uncertainty_values(report: dict[str, object]) -> list[float]:
    values: list[float] = []
    for estimate in report.get("amplitudeEstimates", []):
        if not isinstance(estimate, dict):
            continue
        for channel in ("alpha", "height"):
            section = estimate.get(channel, {})
            if not isinstance(section, dict):
                continue
            value = _finite_or_none(section.get("uncertaintyHalfWidth"))
            if value is not None:
                values.append(value)
    return values


def _amplitude_topology_band_values(report: dict[str, object]) -> list[float]:
    values: list[float] = []
    for band in report.get("topologyResponseBands", []):
        if not isinstance(band, dict) or not bool(band.get("predictedTopologyBlocked")):
            continue
        alpha = _finite_or_none(band.get("alphaResponseBand")) or 0.0
        height = _finite_or_none(band.get("heightResponseBand")) or 0.0
        values.append(max(alpha, height))
    return values


def _bounded_profile_proposal(
    name: str,
    proposed: float | None,
    lower: float,
    upper: float,
    sample_count: int,
    *,
    min_samples: int,
    claim_label: str,
    reason: str,
) -> dict[str, object]:
    finite = proposed is not None and np.isfinite(proposed)
    bounded = finite and lower <= float(proposed) <= upper
    safe = bool(bounded and sample_count >= min_samples)
    return {
        "name": name,
        "configField": None,
        "current": None,
        "proposed": proposed,
        "lowerBound": float(lower),
        "upperBound": float(upper),
        "sampleCount": int(sample_count),
        "bounded": bool(bounded),
        "safeToApply": safe,
        "claimLabel": claim_label,
        "reason": reason if safe else f"{reason}; requires finite bounded value and at least {min_samples} samples",
    }


def reachable_equilibrium_empirical_profile_from_amplitude(
    amplitude_report: dict[str, object],
    *,
    holdout_amplitude_report: dict[str, object] | None = None,
    safety_factor: float = 1.5,
    min_samples: int = 1,
    require_holdout: bool = False,
) -> dict[str, object]:
    """Fit bounded empirical profile proposals from amplitude calibration."""

    factor = max(1.0, float(safety_factor))
    alpha_values = _amplitude_section_values(amplitude_report, "alpha")
    height_values = _amplitude_section_values(amplitude_report, "height")
    uncertainty_values = _amplitude_uncertainty_values(amplitude_report)
    topology_values = _amplitude_topology_band_values(amplitude_report)
    metrics = amplitude_report.get("metrics", {})
    if not isinstance(metrics, dict):
        metrics = {}
    summary = amplitude_report.get("summary", {})
    if not isinstance(summary, dict):
        summary = {}
    alpha_scale = float(np.mean(alpha_values)) if alpha_values else None
    height_scale = float(np.mean(height_values)) if height_values else None
    topology_tolerance = (
        factor * max(topology_values) if topology_values else 0.0
    )
    group_sequence_residual = _finite_or_none(metrics.get("maxGroupSequenceResidual"))
    group_sequence_tolerance = (
        factor * group_sequence_residual
        if group_sequence_residual is not None
        else None
    )
    uncertainty_budget = (
        factor * max(uncertainty_values) if uncertainty_values else None
    )
    group_sequence_count = int(metrics.get("groupSequenceResidualCount") or 0)
    proposals = [
        _bounded_profile_proposal(
            "alphaResponseScale",
            alpha_scale,
            0.0,
            10.0,
            len(alpha_values),
            min_samples=min_samples,
            claim_label="bench-measured empirical alpha response scale",
            reason="mean absolute measured alpha response over reachable targets",
        ),
        _bounded_profile_proposal(
            "heightResponseScale",
            height_scale,
            0.0,
            10.0,
            len(height_values),
            min_samples=min_samples,
            claim_label="bench-measured empirical height response scale",
            reason="mean absolute measured height response over reachable targets",
        ),
        _bounded_profile_proposal(
            "topologyLeakageTolerance",
            topology_tolerance,
            0.0,
            10.0,
            len(topology_values),
            min_samples=0,
            claim_label="bench-measured empirical topology leakage bound",
            reason="safety-factor-scaled maximum blocked-target response band",
        ),
        _bounded_profile_proposal(
            "groupSequenceTolerance",
            group_sequence_tolerance,
            0.0,
            10.0,
            group_sequence_count,
            min_samples=min_samples,
            claim_label="bench-measured empirical group ordering tolerance",
            reason="safety-factor-scaled maximum simultaneous/sequenced residual",
        ),
        _bounded_profile_proposal(
            "uncertaintyBudget",
            uncertainty_budget,
            0.0,
            10.0,
            len(uncertainty_values),
            min_samples=min_samples,
            claim_label="finite repeated-trial uncertainty metadata",
            reason="safety-factor-scaled maximum alpha/height uncertainty half-width",
        ),
    ]
    holdout_summary: dict[str, object] = {
        "attached": holdout_amplitude_report is not None,
        "schema": holdout_amplitude_report.get("schema", "")
        if isinstance(holdout_amplitude_report, dict)
        else "",
        "holdoutReady": False,
        "holdoutPass": not require_holdout and holdout_amplitude_report is None,
        "reason": "holdout not required and not attached",
    }
    if isinstance(holdout_amplitude_report, dict):
        holdout_metrics = holdout_amplitude_report.get("metrics", {})
        holdout_summary_raw = holdout_amplitude_report.get("summary", {})
        if not isinstance(holdout_metrics, dict):
            holdout_metrics = {}
        if not isinstance(holdout_summary_raw, dict):
            holdout_summary_raw = {}
        holdout_ready = bool(holdout_summary_raw.get("amplitudeCalibrationReady"))
        holdout_topology = bool(holdout_summary_raw.get("topologyLeakageBandsPass"))
        holdout_group = _finite_or_none(
            holdout_metrics.get("maxGroupSequenceResidual")
        )
        group_limit = group_sequence_tolerance if group_sequence_tolerance is not None else 0.0
        holdout_pass = bool(
            holdout_ready
            and holdout_topology
            and holdout_group is not None
            and holdout_group <= group_limit
        )
        holdout_summary = {
            "attached": True,
            "schema": holdout_amplitude_report.get("schema", ""),
            "holdoutReady": holdout_ready,
            "holdoutPass": holdout_pass,
            "maxGroupSequenceResidual": holdout_group,
            "groupSequenceTolerance": group_limit,
            "topologyLeakageBandsPass": holdout_topology,
            "reason": (
                "holdout amplitude report fits within empirical profile tolerances"
                if holdout_pass
                else "holdout amplitude report is missing readiness or exceeds empirical tolerances"
            ),
        }
    safe_count = sum(1 for proposal in proposals if proposal["safeToApply"])
    bounded_count = sum(1 for proposal in proposals if proposal["bounded"])
    missing_evidence: list[str] = []
    if not summary.get("amplitudeCalibrationReady"):
        missing_evidence.append("amplitudeCalibration")
    if not alpha_values:
        missing_evidence.append("alphaResponseScale")
    if not height_values:
        missing_evidence.append("heightResponseScale")
    if group_sequence_residual is None:
        missing_evidence.append("groupSequenceResidual")
    if not uncertainty_values:
        missing_evidence.append("uncertaintyBudget")
    if safe_count <= 0:
        missing_evidence.append("safeBoundedProposals")
    if require_holdout and not holdout_summary.get("holdoutPass"):
        missing_evidence.append("holdoutValidation")
    ready = len(missing_evidence) == 0
    return {
        "schema": "rad-sim.reachable-equilibrium-empirical-profile.v1",
        "sourceReportSchema": amplitude_report.get("schema", ""),
        "method": (
            "bounded empirical profile proposals derived from repeated-trial "
            "reachable-equilibrium amplitude calibration"
        ),
        "parameters": {
            "safetyFactor": factor,
            "minSamples": int(min_samples),
            "requireHoldout": bool(require_holdout),
        },
        "recommendedUpdates": proposals,
        "holdoutValidation": holdout_summary,
        "metrics": {
            "proposalCount": len(proposals),
            "boundedProposalCount": bounded_count,
            "safeProposalCount": safe_count,
            "alphaSampleCount": len(alpha_values),
            "heightSampleCount": len(height_values),
            "topologyLeakageBandCount": len(topology_values),
            "uncertaintySampleCount": len(uncertainty_values),
            "groupSequenceResidualCount": group_sequence_count,
        },
        "summary": {
            "status": (
                "reachable-equilibrium-empirical-profile-ready"
                if ready
                else "needs-reachable-equilibrium-profile-review"
            ),
            "empiricalProfileReady": ready,
            "safeBoundedProposalCount": safe_count,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "formalization": {
            "targetId": "reachable_equilibrium_empirical_profile_gate",
            "leanStructure": "Mechanics.ReachableEquilibriumEmpiricalProfileNat",
            "leanPredicate": "reachableEquilibriumEmpiricalProfileReadyNat",
            "schema": "rad-sim.reachable-equilibrium-empirical-profile.v1",
        },
        "claimLabels": {
            "profile": "bench-measured empirical law",
            "safeProposals": "bounded finite-sample update proposals",
            "holdout": "optional empirical validation hook",
            "physicalAccuracy": "not a constitutive mechanics proof",
        },
        "limitations": [
            "Profile proposals are calibration metadata; no LatticeConfig field is mutated in v1.",
            "Bounds prevent nonsensical profile values but do not prove physical correctness.",
            "Holdout validation checks tolerance consistency, not statistical independence unless the lab supplied independent data.",
        ],
    }


def export_reachable_equilibrium_empirical_profile_json(
    amplitude_report: dict[str, object],
    **kwargs: object,
) -> str:
    return json.dumps(
        reachable_equilibrium_empirical_profile_from_amplitude(
            amplitude_report,
            **kwargs,
        ),
        indent=2,
    )


def export_reachable_equilibrium_empirical_profile_csv(
    profile: dict[str, object],
) -> str:
    summary = profile.get("summary", {})
    if not isinstance(summary, dict):
        summary = {}
    header = [
        "schema",
        "profile_ready",
        "name",
        "proposed",
        "lower_bound",
        "upper_bound",
        "sample_count",
        "bounded",
        "safe_to_apply",
        "claim_label",
        "missing_evidence",
    ]
    rows: list[list[object]] = [header]
    for update in profile.get("recommendedUpdates", []):
        if not isinstance(update, dict):
            continue
        rows.append(
            [
                profile.get("schema", ""),
                summary.get("empiricalProfileReady", ""),
                update.get("name", ""),
                update.get("proposed", ""),
                update.get("lowerBound", ""),
                update.get("upperBound", ""),
                update.get("sampleCount", ""),
                update.get("bounded", ""),
                update.get("safeToApply", ""),
                update.get("claimLabel", ""),
                ";".join(str(item) for item in summary.get("missingEvidence", []))
                if isinstance(summary.get("missingEvidence"), list)
                else "",
            ]
        )
    return "\n".join(
        ",".join(_calibration_csv_scalar(value) for value in row) for row in rows
    )


def select_calibration_model_profile(
    comparisons: Iterable[dict[str, object]],
    *,
    min_improvement: float = 0.0,
) -> dict[str, object]:
    """Choose the lowest-residual improving profile comparison conservatively."""

    candidates = list(comparisons)
    scored: list[tuple[float, int, int, dict[str, object]]] = []
    for index, candidate in enumerate(candidates):
        after = candidate.get("after", {})
        after_metrics = after.get("metrics", {}) if isinstance(after, dict) else {}
        delta = candidate.get("delta", {})
        score = _finite_or_none(
            after_metrics.get("residualScore")
            if isinstance(after_metrics, dict)
            else None
        )
        score_delta = _finite_or_none(
            delta.get("residualScore") if isinstance(delta, dict) else None
        )
        missing_delta = (
            int(delta.get("missingObservationCount") or 0)
            if isinstance(delta, dict)
            else 0
        )
        applied_count = 0
        application = candidate.get("application", {})
        if isinstance(application, dict):
            applied = application.get("appliedUpdates", [])
            applied_count = len(applied) if isinstance(applied, list) else 0
        conservative_pass = (
            score is not None
            and score_delta is not None
            and score_delta <= -float(min_improvement)
            and missing_delta <= 0
            and applied_count > 0
        )
        candidate["selectionEligible"] = bool(conservative_pass)
        candidate["selectionRejectionReason"] = (
            None
            if conservative_pass
            else "requires applied updates, non-increased missing count, and residual-score improvement"
        )
        if conservative_pass:
            scored.append((float(score), missing_delta, index, candidate))
    scored.sort(key=lambda item: (item[0], item[1], item[2]))
    selected = scored[0][3] if scored else None
    return {
        "schema": "rad-sim.calibration-model-profile-selection.v1",
        "candidateCount": len(candidates),
        "eligibleCount": len(scored),
        "selectedIndex": None if selected is None else candidates.index(selected),
        "selected": selected,
        "candidates": candidates,
        "rule": {
            "name": "lowest-after-residual-score-with-improvement",
            "minImprovement": float(min_improvement),
            "requiresAppliedUpdate": True,
            "requiresMissingObservationNonincrease": True,
        },
        "claimLabel": "simulator-derived empirical profile selection rule",
        "limitations": [
            "Selection is based on simulator residuals against available measurements.",
            "Independent validation data is still required before treating a profile as physically calibrated.",
        ],
    }


def export_calibration_model_profile_json(
    config: LatticeConfig,
    report: dict[str, object],
    *,
    min_samples: int = 1,
) -> str:
    return json.dumps(
        calibration_model_profile_from_report(
            config,
            report,
            min_samples=min_samples,
        ),
        indent=2,
    )


def export_calibration_model_profile_residual_comparison_json(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    measurements: CalibrationExperimentMeasurements,
    profile: dict[str, object] | None = None,
    *,
    tolerance: float = 1e-9,
    min_improvement: float = 0.0,
) -> str:
    return json.dumps(
        calibration_model_profile_residual_comparison(
            config,
            protocol,
            measurements,
            profile,
            tolerance=tolerance,
            min_improvement=min_improvement,
        ),
        indent=2,
    )


def export_calibration_model_profile_holdout_validation_json(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol,
    fit_measurements: CalibrationExperimentMeasurements,
    holdout_measurements: CalibrationExperimentMeasurements,
    profile: dict[str, object] | None = None,
    *,
    tolerance: float = 1e-9,
    min_improvement: float = 0.0,
) -> str:
    return json.dumps(
        calibration_model_profile_holdout_validation(
            config,
            protocol,
            fit_measurements,
            holdout_measurements,
            profile,
            tolerance=tolerance,
            min_improvement=min_improvement,
        ),
        indent=2,
    )


def export_calibration_model_profile_selection_json(
    comparisons: Iterable[dict[str, object]],
    *,
    min_improvement: float = 0.0,
) -> str:
    return json.dumps(
        select_calibration_model_profile(
            comparisons,
            min_improvement=min_improvement,
        ),
        indent=2,
    )


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
    parameter_estimates = _calibration_parameter_estimates(
        config,
        protocol,
        measurements,
        comparison_tuple,
        fit,
        summary,
    )
    model_profile = _calibration_model_profile(
        config,
        parameter_estimates,
        source_report_schema="rad-sim.calibration-comparison-report.v1",
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
        "parameterEstimates": parameter_estimates,
        "modelProfile": model_profile,
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


def build_response_atlas(
    config: LatticeConfig,
    protocol: CalibrationExperimentProtocol | None = None,
    *,
    center_cell: tuple[int, int] | None = None,
    alpha_step: float = -0.25,
    z_step: float = 0.30,
    physical: bool = False,
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
) -> ResponseAtlas:
    """Build a compact single/pair/cluster response atlas.

    The atlas is a simulator-side research artifact. It reuses the calibration
    protocol steps, records kinematic response fields at observation cells, and
    optionally attaches spring-hinge model-disagreement metrics for each step.
    """

    if protocol is None:
        protocol = build_calibration_experiment_protocol(
            config,
            center_cell=center_cell,
            alpha_step=alpha_step,
            z_step=z_step,
        )
    entries: list[ResponseAtlasEntry] = []
    for step in protocol.steps:
        locked_cells = tuple(step.locked_cells)
        commands = tuple(step.commands)
        response = characterize_response(
            config,
            commands,
            locked_cells=locked_cells,
            tolerance=tolerance,
        )
        alpha_error, height_error = _atlas_superposition_errors(
            config,
            commands,
            response,
            locked_cells,
            tolerance,
        )
        physical_height_rms_error: float | None = None
        physical_center_rms_error: float | None = None
        physical_success: bool | None = None
        physical_energy: float | None = None
        if physical:
            comparison = compare_physical_response(
                config,
                commands,
                locked_cells=locked_cells,
                load_case=load_case,
                tolerance=tolerance,
            )
            physical_height_rms_error = comparison.height_rms_error
            physical_center_rms_error = comparison.center_rms_error
            physical_success = comparison.physical_success
            physical_energy = comparison.physical_energy
        entries.append(
            ResponseAtlasEntry(
                step_id=step.id,
                scope=step.scope,
                purpose=step.purpose,
                expected_response=step.expected_response,
                commands=commands,
                locked_cells=locked_cells,
                observation_cells=_atlas_observation_cells(
                    response,
                    step.observation_cells,
                ),
                alpha_reach=response.alpha_reach,
                z_reach=response.z_reach,
                effective_alpha_die_off=response.effective_alpha_die_off,
                effective_z_die_off=response.effective_z_die_off,
                max_abs_alpha_delta=response.max_abs_alpha_delta,
                max_abs_height_delta=response.max_abs_height_delta,
                mean_abs_alpha_delta=float(np.mean(np.abs(response.alpha_delta))),
                mean_abs_height_delta=float(np.mean(np.abs(response.height_delta))),
                alpha_superposition_error=alpha_error,
                height_superposition_error=height_error,
                physical_height_rms_error=physical_height_rms_error,
                physical_center_rms_error=physical_center_rms_error,
                physical_success=physical_success,
                physical_energy=physical_energy,
            )
        )
    return ResponseAtlas(
        config_shape=(config.rows, config.cols),
        center_cell=protocol.center_cell,
        entries=tuple(entries),
        physical=physical,
    )


def sweep_response_atlas_parameters(
    config: LatticeConfig,
    *,
    backlash_values: Iterable[float] | None = None,
    clearance_values: Iterable[float] | None = None,
    protocol: CalibrationExperimentProtocol | None = None,
    center_cell: tuple[int, int] | None = None,
    alpha_step: float = -0.25,
    z_step: float = 0.30,
    physical: bool = False,
    load_case: LoadCase | None = None,
    tolerance: float = 1e-9,
) -> ResponseAtlasSweep:
    """Sweep backlash and pin-hole clearance through the atlas protocol.

    This is a numerical experiment for programmable-discontinuity parameter
    studies. The only varied fields are the alpha dead-zone (`backlash`) and
    vertical dead-zone (`hole_radius - pin_radius`); all other configuration
    values are inherited from the supplied base config.
    """

    backlash_tuple = tuple(
        float(value)
        for value in (
            (config.backlash,) if backlash_values is None else backlash_values
        )
    )
    clearance_tuple = tuple(
        float(value)
        for value in (
            (config.pin_hole_clearance,)
            if clearance_values is None
            else clearance_values
        )
    )
    if not backlash_tuple:
        raise ValueError("backlash_values must contain at least one value")
    if not clearance_tuple:
        raise ValueError("clearance_values must contain at least one value")
    if any(value < 0 for value in clearance_tuple):
        raise ValueError("clearance_values must be non-negative")
    if protocol is not None and protocol.config_shape != (config.rows, config.cols):
        raise ValueError("protocol config_shape must match config shape")

    center = (
        protocol.center_cell
        if protocol is not None
        else _clamp_cell(
            config,
            center_cell if center_cell is not None else (config.rows // 2, config.cols // 2),
        )
    )
    samples: list[ResponseAtlasSweepSample] = []
    for backlash in backlash_tuple:
        for clearance in clearance_tuple:
            sample_config = replace(
                config,
                backlash=backlash,
                hole_radius=config.pin_radius + clearance,
            )
            atlas = build_response_atlas(
                sample_config,
                protocol=protocol,
                center_cell=center,
                alpha_step=alpha_step,
                z_step=z_step,
                physical=physical,
                load_case=load_case,
                tolerance=tolerance,
            )
            samples.append(
                ResponseAtlasSweepSample(
                    backlash=sample_config.backlash,
                    pin_hole_clearance=sample_config.pin_hole_clearance,
                    pin_radius=sample_config.pin_radius,
                    hole_radius=sample_config.hole_radius,
                    atlas=atlas,
                )
            )

    return ResponseAtlasSweep(
        config_shape=(config.rows, config.cols),
        center_cell=center,
        backlash_values=backlash_tuple,
        clearance_values=clearance_tuple,
        samples=tuple(samples),
        physical=physical,
    )


def export_response_atlas_json(atlas: ResponseAtlas) -> str:
    return json.dumps(atlas.to_dict(), indent=2)


def export_response_atlas_sweep_json(sweep: ResponseAtlasSweep) -> str:
    return json.dumps(sweep.to_dict(), indent=2)


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
    removed_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
) -> ResponseMatrix:
    if not include_alpha and not include_z:
        raise ValueError("at least one command family must be included")
    removed = {(int(r), int(c)) for r, c in removed_cells}
    raw_cells = (
        actuator_cells
        if actuator_cells is not None
        else ((r, c) for r in range(config.rows) for c in range(config.cols))
    )
    cells = tuple(
        (int(r), int(c)) for r, c in raw_cells if (int(r), int(c)) not in removed
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
        response = characterize_response(
            config,
            (command,),
            locked_cells,
            tolerance,
            removed_cells=removed,
        )
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


def compare_removed_topology_reachability(
    config: LatticeConfig,
    removed_cells: Iterable[tuple[int, int]],
    *,
    actuator_cells: Iterable[tuple[int, int]] | None = None,
    locked_cells: Iterable[tuple[int, int]] = (),
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    include_alpha: bool = True,
    include_z: bool = True,
    tolerance: float = 1e-9,
) -> RemovedTopologyReachabilityComparison:
    """Compare finite-response reachability before and after deleting cells.

    The same candidate actuator cells are used in both response matrices, with
    removed cells excluded from that candidate set. The resulting losses isolate
    topology and propagation changes from the unrelated fact that a removed cell
    cannot itself host an actuator.
    """
    removed = tuple((int(r), int(c)) for r, c in removed_cells)
    removed_set = set(removed)
    locked = tuple((int(r), int(c)) for r, c in locked_cells)
    raw_cells = (
        actuator_cells
        if actuator_cells is not None
        else ((r, c) for r in range(config.rows) for c in range(config.cols))
    )
    candidates = tuple(
        (int(r), int(c))
        for r, c in raw_cells
        if (int(r), int(c)) not in removed_set
    )
    intact_state = _state_with_commands(config, (), locked)
    removed_state = _state_with_commands(config, (), locked, removed)
    intact_matrix = build_response_matrix(
        config,
        actuator_cells=candidates,
        alpha_step=alpha_step,
        z_step=z_step,
        include_alpha=include_alpha,
        include_z=include_z,
        locked_cells=locked,
        tolerance=tolerance,
    )
    removed_matrix = build_response_matrix(
        config,
        actuator_cells=candidates,
        alpha_step=alpha_step,
        z_step=z_step,
        include_alpha=include_alpha,
        include_z=include_z,
        locked_cells=locked,
        removed_cells=removed,
        tolerance=tolerance,
    )
    return RemovedTopologyReachabilityComparison(
        removed_cells=removed,
        actuator_cells=candidates,
        locked_cells=locked,
        intact_matrix=intact_matrix,
        removed_matrix=removed_matrix,
        intact_topology=lattice_topology_diagnostic(config, intact_state),
        removed_topology=lattice_topology_diagnostic(config, removed_state),
        tolerance=tolerance,
    )


def _neighbor_cells(config: LatticeConfig, cell: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    r, c = cell
    cells: list[tuple[int, int]] = []
    for rr, cc in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
        if 0 <= rr < config.rows and 0 <= cc < config.cols:
            cells.append((rr, cc))
    return tuple(cells)


def default_topology_experiment_scenarios(
    config: LatticeConfig,
    *,
    center_cell: tuple[int, int] | None = None,
    alpha: float = -0.25,
    z: float = 0.2,
) -> tuple[TopologyExperimentScenario, ...]:
    center = (
        (config.rows // 2, config.cols // 2)
        if center_cell is None
        else (int(center_cell[0]), int(center_cell[1]))
    )
    group_cells = (center, *_neighbor_cells(config, center))
    group_commands = tuple(SourceCommand(cell=cell, alpha=alpha, z=z) for cell in group_cells)
    return (
        TopologyExperimentScenario(
            name="intact-group-actuation",
            commands=group_commands,
            actuator_cells=group_cells,
        ),
        TopologyExperimentScenario(
            name="locked-center-group-actuation",
            commands=group_commands,
            locked_cells=(center,),
            actuator_cells=group_cells,
        ),
        TopologyExperimentScenario(
            name="removed-center-group-actuation",
            commands=group_commands,
            removed_cells=(center,),
            actuator_cells=group_cells,
        ),
    )


def _scenario_candidate_cells(
    config: LatticeConfig,
    scenario: TopologyExperimentScenario,
) -> tuple[tuple[int, int], ...]:
    removed = set(scenario.removed_cells)
    locked = set(scenario.locked_cells)
    raw_cells = (
        scenario.actuator_cells
        if scenario.actuator_cells is not None
        else ((r, c) for r in range(config.rows) for c in range(config.cols))
    )
    return tuple(
        (int(r), int(c))
        for r, c in raw_cells
        if (int(r), int(c)) not in removed and (int(r), int(c)) not in locked
    )


def build_topology_experiment_report(
    config: LatticeConfig,
    scenarios: Iterable[TopologyExperimentScenario] | None = None,
    *,
    target_alpha: np.ndarray | Iterable[Iterable[float]] | None = None,
    target_height: np.ndarray | Iterable[Iterable[float]] | None = None,
    alpha_step: float = 0.12,
    z_step: float = 0.12,
    include_alpha: bool = True,
    include_z: bool = True,
    tolerance: float = 1e-9,
) -> TopologyExperimentReport:
    """Run matched topology/reachability diagnostics across discontinuity states."""
    shape = (config.rows, config.cols)
    target_alpha_matrix = _target_matrix(target_alpha, shape, "target_alpha")
    target_height_matrix = _target_matrix(target_height, shape, "target_height")
    scenario_tuple = tuple(
        scenario.normalized()
        for scenario in (
            default_topology_experiment_scenarios(config)
            if scenarios is None
            else scenarios
        )
    )
    results: list[TopologyExperimentScenarioResult] = []
    for scenario in scenario_tuple:
        candidates = _scenario_candidate_cells(config, scenario)
        state = _state_with_commands(
            config,
            scenario.commands,
            scenario.locked_cells,
            scenario.removed_cells,
        )
        sim = simulate_kinematic(config, state)
        topology = lattice_topology_diagnostic(config, state)
        matrix = build_response_matrix(
            config,
            actuator_cells=candidates,
            alpha_step=alpha_step,
            z_step=z_step,
            include_alpha=include_alpha,
            include_z=include_z,
            locked_cells=scenario.locked_cells,
            removed_cells=scenario.removed_cells,
            tolerance=tolerance,
        )
        response = (
            None
            if not scenario.commands
            else characterize_response(
                config,
                scenario.commands,
                scenario.locked_cells,
                tolerance,
                removed_cells=scenario.removed_cells,
            )
        )
        component_reachable, component_blocked = _component_reach_counts(
            topology,
            candidates,
        )
        component_summaries = _component_response_summaries(
            topology,
            matrix,
            tolerance,
        )
        positive_height, negative_height = (
            _signed_height_reach_counts(
                config,
                candidates,
                scenario.locked_cells,
                scenario.removed_cells,
                z_step,
                tolerance,
            )
            if include_z
            else (0, 0)
        )
        alpha_residual = (
            None if target_alpha_matrix is None else target_alpha_matrix - sim.alpha
        )
        height_residual = (
            None
            if target_height_matrix is None
            else target_height_matrix - sim.metadata["height"]
        )
        results.append(
            TopologyExperimentScenarioResult(
                scenario=scenario,
                topology=topology,
                response_matrix=matrix,
                response=response,
                component_reachable_cells=component_reachable,
                component_blocked_cells=component_blocked,
                positive_height_reachable_cells=positive_height,
                negative_height_reachable_cells=negative_height,
                component_summaries=component_summaries,
                target_alpha_rms=_array_rms(alpha_residual),
                target_height_rms=_array_rms(height_residual),
                max_abs_target_alpha_residual=_array_max_abs(alpha_residual),
                max_abs_target_height_residual=_array_max_abs(height_residual),
                tolerance=tolerance,
            )
        )
    return TopologyExperimentReport(
        config_shape=shape,
        scenarios=tuple(results),
        target_alpha=target_alpha_matrix,
        target_height=target_height_matrix,
    )


def export_response_matrix_json(
    matrix: ResponseMatrix,
    *,
    tolerance: float = 1e-9,
) -> str:
    return json.dumps(matrix.to_dict(tolerance=tolerance), indent=2)


def export_removed_topology_reachability_comparison_json(
    comparison: RemovedTopologyReachabilityComparison,
    *,
    tolerance: float | None = None,
) -> str:
    return json.dumps(comparison.to_dict(tolerance=tolerance), indent=2)


def export_topology_experiment_report_json(
    report: TopologyExperimentReport,
    *,
    tolerance: float | None = None,
) -> str:
    return json.dumps(report.to_dict(tolerance=tolerance), indent=2)
