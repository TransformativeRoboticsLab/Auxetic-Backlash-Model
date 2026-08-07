from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

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


def characterize_cluster(
    config: LatticeConfig,
    commands: Iterable[SourceCommand],
    locked_cells: Iterable[tuple[int, int]] = (),
    tolerance: float = 1e-9,
) -> ResponseCharacterization:
    return characterize_response(config, tuple(commands), locked_cells, tolerance)


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
