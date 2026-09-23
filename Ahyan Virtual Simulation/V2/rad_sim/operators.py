from __future__ import annotations

import json
from collections import deque
from dataclasses import dataclass
from typing import Any, Iterable, Mapping

import numpy as np

from .coupling import backlash_activation, normalized_backlash_to_alpha_dead_zone
from .geometry import neighbor_indices
from .models import LatticeConfig, LatticeState, LoadCase


@dataclass(frozen=True)
class PropagationResult:
    field: np.ndarray
    die_off: np.ndarray


@dataclass(frozen=True)
class ProgrammableDiscontinuityEvent:
    kind: str
    cell: tuple[int, int] | None = None
    cells: tuple[tuple[int, int], ...] = ()
    alpha: float = 0.0
    z: float = 0.0


@dataclass(frozen=True)
class EventCommutativityDiagnostic:
    first_then_second: LatticeState
    second_then_first: LatticeState
    mode_commutes: bool
    command_commutes: bool
    alpha_grid_commutes: bool
    final_alpha_error: float
    final_height_error: float


@dataclass(frozen=True)
class StateDistanceDiagnostic:
    mode_changes: int
    command_alpha_error: float
    command_z_error: float
    command_error: float
    alpha_grid_error: float
    final_alpha_error: float
    final_height_error: float
    final_error: float


@dataclass(frozen=True)
class AdjacentSwapDiagnostic:
    index: int
    first_kind: str
    second_kind: str
    sensitive: bool
    distance: StateDistanceDiagnostic


@dataclass(frozen=True)
class SequenceOrderDiagnostic:
    event_count: int
    adjacent_pair_count: int
    noncommuting_adjacent_pairs: int
    order_sensitive: bool
    reverse: StateDistanceDiagnostic
    max_adjacent_alpha_error: float
    max_adjacent_height_error: float
    max_adjacent_command_error: float
    max_adjacent_alpha_grid_error: float
    max_order_error: float
    adjacent: tuple[AdjacentSwapDiagnostic, ...]


@dataclass(frozen=True)
class GroupActuationDecompositionDiagnostic:
    cells: tuple[tuple[int, int], ...]
    applied_cells: tuple[tuple[int, int], ...]
    skipped_removed_cells: tuple[tuple[int, int], ...]
    group_event: ProgrammableDiscontinuityEvent
    local_events: tuple[ProgrammableDiscontinuityEvent, ...]
    simultaneous_state: LatticeState
    sequenced_state: LatticeState
    distance: StateDistanceDiagnostic
    decomposes: bool


@dataclass(frozen=True)
class GroupRemovalInteractionDiagnostic:
    cells: tuple[tuple[int, int], ...]
    removed_cells: tuple[tuple[int, int], ...]
    removal_events: tuple[ProgrammableDiscontinuityEvent, ...]
    intact: GroupActuationDecompositionDiagnostic
    removed: GroupActuationDecompositionDiagnostic
    lost_applied_cells: tuple[tuple[int, int], ...]
    newly_skipped_cells: tuple[tuple[int, int], ...]
    intact_topology: dict[str, Any]
    removed_topology: dict[str, Any]
    distance: StateDistanceDiagnostic
    component_count_delta: int
    deleted_edge_delta: int


@dataclass(frozen=True)
class VerticalRemovalResidualDiagnostic:
    cells: tuple[tuple[int, int], ...]
    removed_cells: tuple[tuple[int, int], ...]
    fixed_cells: tuple[tuple[int, int], ...]
    removal_events: tuple[ProgrammableDiscontinuityEvent, ...]
    intact_state: LatticeState
    removed_state: LatticeState
    intact_topology: dict[str, Any]
    removed_topology: dict[str, Any]
    intact_z_residual: np.ndarray
    removed_z_residual: np.ndarray
    residual_delta: np.ndarray
    intact_z_die_off: np.ndarray
    removed_z_die_off: np.ndarray
    affected_neighbor_cells: tuple[tuple[int, int], ...]
    lost_residual_cells: tuple[tuple[int, int], ...]
    topology_blocked_cells: tuple[tuple[int, int], ...]
    removed_source_cells: tuple[tuple[int, int], ...]
    intact_reach_count: int
    removed_reach_count: int
    intact_neighbor_reach_count: int
    removed_neighbor_reach_count: int
    intact_die_off_radius: int
    removed_die_off_radius: int
    max_abs_neighbor_residual: float
    max_abs_residual_delta: float
    clearance: float
    contact_stiffness: float
    intact_contact_engaged_cells: tuple[tuple[int, int], ...]
    removed_contact_engaged_cells: tuple[tuple[int, int], ...]
    intact_contact_penalty_energy: float
    removed_contact_penalty_energy: float
    contact_penalty_delta: float
    external_z_load: float
    intact_load_active_cells: tuple[tuple[int, int], ...]
    removed_load_active_cells: tuple[tuple[int, int], ...]
    intact_load_work: float
    removed_load_work: float
    load_work_delta: float
    intact_load_work_magnitude: float
    removed_load_work_magnitude: float
    load_work_magnitude_delta: float
    intact_height_contact_engaged_cells: tuple[tuple[int, int], ...]
    removed_height_contact_engaged_cells: tuple[tuple[int, int], ...]
    intact_height_contact_penalty_energy: float
    removed_height_contact_penalty_energy: float
    height_contact_penalty_delta: float
    component_count_delta: int
    deleted_edge_delta: int


@dataclass(frozen=True)
class SolverEnergyBreakdown:
    objective_energy: float
    stored_energy: float
    spring_energy: float
    hinge_energy: float
    lock_penalty_energy: float
    external_potential_energy: float


@dataclass(frozen=True)
class VerticalRemovalPhysicalComparisonDiagnostic:
    vertical_diagnostic: VerticalRemovalResidualDiagnostic
    intact_load_case: LoadCase
    removed_load_case: LoadCase
    intact_solver_success: bool
    removed_solver_success: bool
    physical_success: bool
    intact_solver_energy: float
    removed_solver_energy: float
    intact_energy_breakdown: SolverEnergyBreakdown
    removed_energy_breakdown: SolverEnergyBreakdown
    intact_solver_iterations: int
    removed_solver_iterations: int
    intact_spring_edges: int
    removed_spring_edges: int
    spring_edge_delta: int
    intact_hinge_triples: int
    removed_hinge_triples: int
    hinge_triple_delta: int
    intact_solver_removed_cells: int
    removed_solver_removed_cells: int
    intact_solver_load_cells: tuple[tuple[int, int], ...]
    removed_solver_load_cells: tuple[tuple[int, int], ...]
    intact_unanchored_load_cells: tuple[tuple[int, int], ...]
    removed_unanchored_load_cells: tuple[tuple[int, int], ...]
    intact_physical_height: np.ndarray
    removed_physical_height: np.ndarray
    intact_kinematic_height: np.ndarray
    removed_kinematic_height: np.ndarray
    intact_model_error: np.ndarray
    removed_model_error: np.ndarray
    physical_height_delta: np.ndarray
    kinematic_height_delta: np.ndarray
    physical_vs_kinematic_delta_error: np.ndarray
    intact_height_rms_model_error: float
    removed_height_rms_model_error: float
    delta_rms_model_error: float
    max_abs_delta_model_error: float
    solver_topology_note: str


@dataclass(frozen=True)
class VerticalLoadPhysicalScenario:
    name: str
    cells: tuple[tuple[int, int], ...]
    removed_cells: tuple[tuple[int, int], ...] = ()
    fixed_cells: tuple[tuple[int, int], ...] = ()
    z: float = 0.0
    alpha: float = 0.0
    external_z_load: float = 0.0


@dataclass(frozen=True)
class DeadZonePropagationOperator:
    name: str
    dead_zone: float
    gain: float
    max_steps: int | None = None

    def transmit(self, signal: float) -> float:
        return float(backlash_activation(signal, self.dead_zone)) * self.gain

    def propagate(
        self,
        rows: int,
        cols: int,
        sources: Iterable[tuple[int, int, float]],
        active_mask: np.ndarray | None = None,
    ) -> PropagationResult:
        field = np.zeros((rows, cols), dtype=float)
        die_off = np.full((rows, cols), np.inf, dtype=float)
        steps = self.max_steps or (rows + cols)
        active = (
            np.ones((rows, cols), dtype=bool)
            if active_mask is None
            else np.asarray(active_mask, dtype=bool)
        )
        if active.shape != (rows, cols):
            raise ValueError(f"active_mask must have shape {(rows, cols)}")

        for r0, c0, value in sources:
            if not active[r0, c0]:
                continue
            if abs(value) < 1e-12:
                continue
            queue: deque[tuple[int, int, float, int]] = deque([(r0, c0, value, 0)])
            visited: dict[tuple[int, int], float] = {}
            while queue:
                r, c, signal, dist = queue.popleft()
                if not active[r, c]:
                    continue
                if dist > steps:
                    continue
                key = (r, c)
                if abs(signal) <= abs(visited.get(key, 0.0)):
                    continue
                visited[key] = signal
                field[r, c] += signal
                die_off[r, c] = min(die_off[r, c], dist)
                next_signal = self.transmit(signal)
                if abs(next_signal) < 1e-9:
                    continue
                for nr, nc in neighbor_indices(rows, cols, r, c):
                    if not active[nr, nc]:
                        continue
                    queue.append((nr, nc, next_signal, dist + 1))

        return PropagationResult(field=field, die_off=die_off)


def finite_die_off_radius(die_off: np.ndarray) -> int:
    finite = np.asarray(die_off)[np.isfinite(die_off)]
    if finite.size == 0:
        return 0
    return int(np.max(finite))


def lattice_topology_diagnostic(
    config: LatticeConfig, state: LatticeState
) -> dict[str, Any]:
    state = state.normalized(config)
    present = ~state.removed_mask
    labels = np.full((config.rows, config.cols), -1, dtype=int)
    component_sizes: list[int] = []
    component = 0

    for r0 in range(config.rows):
        for c0 in range(config.cols):
            if not present[r0, c0] or labels[r0, c0] >= 0:
                continue
            queue: deque[tuple[int, int]] = deque([(r0, c0)])
            labels[r0, c0] = component
            size = 0
            while queue:
                r, c = queue.popleft()
                size += 1
                for nr, nc in neighbor_indices(config.rows, config.cols, r, c):
                    if not present[nr, nc] or labels[nr, nc] >= 0:
                        continue
                    labels[nr, nc] = component
                    queue.append((nr, nc))
            component_sizes.append(size)
            component += 1

    total_edges = (
        config.rows * max(0, config.cols - 1)
        + max(0, config.rows - 1) * config.cols
    )
    active_edges = 0
    for r in range(config.rows):
        for c in range(config.cols):
            if c + 1 < config.cols and present[r, c] and present[r, c + 1]:
                active_edges += 1
            if r + 1 < config.rows and present[r, c] and present[r + 1, c]:
                active_edges += 1

    return {
        "present_mask": present,
        "component_labels": labels,
        "component_count": component,
        "component_sizes": tuple(component_sizes),
        "largest_component_size": max(component_sizes, default=0),
        "present_cell_count": int(np.count_nonzero(present)),
        "removed_cell_count": int(np.count_nonzero(state.removed_mask)),
        "total_edge_count": int(total_edges),
        "active_edge_count": int(active_edges),
        "deleted_edge_count": int(total_edges - active_edges),
    }


def nonzero_sources(values: np.ndarray) -> list[tuple[int, int, float]]:
    return [
        (int(r), int(c), float(values[r, c]))
        for r, c in zip(*np.nonzero(np.abs(values) > 1e-12), strict=False)
    ]


def alpha_backlash_operator(config: LatticeConfig) -> DeadZonePropagationOperator:
    return DeadZonePropagationOperator(
        name="alpha_backlash",
        dead_zone=float(normalized_backlash_to_alpha_dead_zone(config.backlash)),
        gain=config.coupling_gain,
        max_steps=config.max_coupling_steps,
    )


def vertical_clearance_operator(config: LatticeConfig) -> DeadZonePropagationOperator:
    return DeadZonePropagationOperator(
        name="vertical_clearance",
        dead_zone=config.pin_hole_clearance,
        gain=config.z_coupling_gain,
        max_steps=config.max_coupling_steps,
    )


def lock_projection(reference: np.ndarray, proposed: np.ndarray, locked: np.ndarray) -> np.ndarray:
    return np.where(locked, reference, proposed)


def local_actuation_event(
    cell: tuple[int, int], alpha: float = 0.0, z: float = 0.0
) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="actuate", cell=cell, alpha=alpha, z=z)


def group_actuation_event(
    cells: Iterable[tuple[int, int]], alpha: float = 0.0, z: float = 0.0
) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(
        kind="group_actuate",
        cells=tuple((int(r), int(c)) for r, c in cells),
        alpha=alpha,
        z=z,
    )


def lock_event(cell: tuple[int, int]) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="lock", cell=cell)


def release_event(cell: tuple[int, int]) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="release", cell=cell)


def clear_actuation_event(
    cell: tuple[int, int] | None = None,
) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="clear_actuation", cell=cell)


def remove_cell_event(cell: tuple[int, int]) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="remove_cell", cell=cell)


def restore_cell_event(cell: tuple[int, int]) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="restore_cell", cell=cell)


def _validate_cell(config: LatticeConfig, cell: tuple[int, int] | None) -> tuple[int, int]:
    if cell is None:
        raise ValueError("event requires a cell")
    r, c = int(cell[0]), int(cell[1])
    if not (0 <= r < config.rows and 0 <= c < config.cols):
        raise ValueError(f"cell {(r, c)} is outside lattice shape {(config.rows, config.cols)}")
    return r, c


def _validate_cells(
    config: LatticeConfig, cells: Iterable[tuple[int, int]]
) -> tuple[tuple[int, int], ...]:
    return tuple(_validate_cell(config, cell) for cell in cells)


def apply_programmable_event(
    config: LatticeConfig,
    state: LatticeState,
    event: ProgrammableDiscontinuityEvent,
) -> LatticeState:
    current = state.normalized(config)
    next_state = LatticeState(
        alpha_grid=current.alpha_grid.copy(),
        locked_mask=current.locked_mask.copy(),
        actuator_grid=current.actuator_grid.copy(),
        z_actuator_grid=current.z_actuator_grid.copy(),
        removed_mask=current.removed_mask.copy(),
        position_locked_mask=current.position_locked_mask.copy(),
        lock_z_grid=current.lock_z_grid.copy(),
    )

    if event.kind == "actuate":
        r, c = _validate_cell(config, event.cell)
        if next_state.removed_mask[r, c]:
            return next_state
        next_state.actuator_grid[r, c] += float(event.alpha)
        next_state.z_actuator_grid[r, c] += float(event.z)
    elif event.kind == "group_actuate":
        for r, c in _validate_cells(config, event.cells):
            if next_state.removed_mask[r, c]:
                continue
            next_state.actuator_grid[r, c] += float(event.alpha)
            next_state.z_actuator_grid[r, c] += float(event.z)
    elif event.kind == "lock":
        r, c = _validate_cell(config, event.cell)
        if next_state.removed_mask[r, c]:
            return next_state
        fields = evaluate_programmable_operators(config, current)
        next_state.alpha_grid[r, c] = fields["alpha"][r, c]
        next_state.lock_z_grid[r, c] = fields["height"][r, c]
        next_state.actuator_grid[r, c] = 0.0
        next_state.z_actuator_grid[r, c] = 0.0
        next_state.locked_mask[r, c] = True
    elif event.kind == "release":
        r, c = _validate_cell(config, event.cell)
        next_state.locked_mask[r, c] = False
    elif event.kind == "clear_actuation":
        if event.cell is None:
            next_state.actuator_grid[:, :] = 0.0
            next_state.z_actuator_grid[:, :] = 0.0
        else:
            r, c = _validate_cell(config, event.cell)
            next_state.actuator_grid[r, c] = 0.0
            next_state.z_actuator_grid[r, c] = 0.0
    elif event.kind == "remove_cell":
        r, c = _validate_cell(config, event.cell)
        next_state.alpha_grid[r, c] = config.initial_alpha
        next_state.locked_mask[r, c] = False
        next_state.actuator_grid[r, c] = 0.0
        next_state.z_actuator_grid[r, c] = 0.0
        next_state.lock_z_grid[r, c] = 0.0
        next_state.removed_mask[r, c] = True
        next_state.position_locked_mask[r, c] = False
    elif event.kind == "restore_cell":
        r, c = _validate_cell(config, event.cell)
        next_state.removed_mask[r, c] = False
    else:
        raise ValueError(f"unsupported programmable discontinuity event: {event.kind}")

    return next_state


def apply_event_sequence(
    config: LatticeConfig,
    state: LatticeState,
    events: Iterable[ProgrammableDiscontinuityEvent],
) -> LatticeState:
    current = state.normalized(config)
    for event in events:
        current = apply_programmable_event(config, current, event)
    return current


def compare_event_order(
    config: LatticeConfig,
    state: LatticeState,
    first: ProgrammableDiscontinuityEvent,
    second: ProgrammableDiscontinuityEvent,
    tolerance: float = 1e-9,
) -> EventCommutativityDiagnostic:
    first_then_second = apply_event_sequence(config, state, (first, second))
    second_then_first = apply_event_sequence(config, state, (second, first))
    first_fields = evaluate_programmable_operators(config, first_then_second)
    second_fields = evaluate_programmable_operators(config, second_then_first)
    final_alpha_error = float(np.max(np.abs(first_fields["alpha"] - second_fields["alpha"])))
    final_height_error = float(
        np.max(np.abs(first_fields["height"] - second_fields["height"]))
    )
    return EventCommutativityDiagnostic(
        first_then_second=first_then_second,
        second_then_first=second_then_first,
        mode_commutes=bool(
            np.array_equal(first_then_second.locked_mask, second_then_first.locked_mask)
            and np.array_equal(
                first_then_second.removed_mask, second_then_first.removed_mask
            )
        ),
        command_commutes=bool(
            np.allclose(
                first_then_second.actuator_grid,
                second_then_first.actuator_grid,
                atol=tolerance,
                rtol=0.0,
            )
            and np.allclose(
                first_then_second.z_actuator_grid,
                second_then_first.z_actuator_grid,
                atol=tolerance,
                rtol=0.0,
            )
        ),
        alpha_grid_commutes=bool(
            np.allclose(
                first_then_second.alpha_grid,
                second_then_first.alpha_grid,
                atol=tolerance,
                rtol=0.0,
            )
        ),
        final_alpha_error=final_alpha_error,
        final_height_error=final_height_error,
    )


def _mode_change_count(first: np.ndarray, second: np.ndarray) -> int:
    return int(
        np.count_nonzero(np.asarray(first, dtype=bool) != np.asarray(second, dtype=bool))
    )


def _state_distance(
    config: LatticeConfig,
    first: LatticeState,
    second: LatticeState,
    first_fields: dict[str, np.ndarray] | None = None,
    second_fields: dict[str, np.ndarray] | None = None,
) -> StateDistanceDiagnostic:
    first_state = first.normalized(config)
    second_state = second.normalized(config)
    if first_fields is None:
        first_fields = evaluate_programmable_operators(config, first_state)
    if second_fields is None:
        second_fields = evaluate_programmable_operators(config, second_state)
    command_alpha_error = float(
        np.max(np.abs(first_state.actuator_grid - second_state.actuator_grid))
    )
    command_z_error = float(
        np.max(np.abs(first_state.z_actuator_grid - second_state.z_actuator_grid))
    )
    alpha_grid_error = float(
        np.max(np.abs(first_state.alpha_grid - second_state.alpha_grid))
    )
    final_alpha_error = float(np.max(np.abs(first_fields["alpha"] - second_fields["alpha"])))
    final_height_error = float(
        np.max(np.abs(first_fields["height"] - second_fields["height"]))
    )
    command_error = max(command_alpha_error, command_z_error)
    return StateDistanceDiagnostic(
        mode_changes=(
            _mode_change_count(first_state.locked_mask, second_state.locked_mask)
            + _mode_change_count(first_state.removed_mask, second_state.removed_mask)
        ),
        command_alpha_error=command_alpha_error,
        command_z_error=command_z_error,
        command_error=command_error,
        alpha_grid_error=alpha_grid_error,
        final_alpha_error=final_alpha_error,
        final_height_error=final_height_error,
        final_error=max(final_alpha_error, final_height_error),
    )


def compare_sequence_order(
    config: LatticeConfig,
    state: LatticeState,
    events: Iterable[ProgrammableDiscontinuityEvent],
    tolerance: float = 1e-9,
) -> SequenceOrderDiagnostic:
    """Measure path dependence by reversing a sequence and swapping neighbors.

    This is a simulator diagnostic for programmable-discontinuity operator
    composition. It is not a paper-derived constitutive law.
    """

    sequence = tuple(events)
    base_final = apply_event_sequence(config, state, sequence)
    base_fields = evaluate_programmable_operators(config, base_final)
    reverse_final = apply_event_sequence(config, state, reversed(sequence))
    reverse = _state_distance(config, base_final, reverse_final, base_fields)

    adjacent: list[AdjacentSwapDiagnostic] = []
    max_adjacent_alpha_error = 0.0
    max_adjacent_height_error = 0.0
    max_adjacent_command_error = 0.0
    max_adjacent_alpha_grid_error = 0.0
    noncommuting_adjacent_pairs = 0
    for index in range(max(0, len(sequence) - 1)):
        swapped = list(sequence)
        swapped[index], swapped[index + 1] = swapped[index + 1], swapped[index]
        swapped_final = apply_event_sequence(config, state, swapped)
        distance = _state_distance(config, base_final, swapped_final, base_fields)
        sensitive = bool(
            distance.mode_changes > 0
            or distance.command_error > tolerance
            or distance.alpha_grid_error > tolerance
            or distance.final_alpha_error > tolerance
            or distance.final_height_error > tolerance
        )
        if sensitive:
            noncommuting_adjacent_pairs += 1
        max_adjacent_alpha_error = max(max_adjacent_alpha_error, distance.final_alpha_error)
        max_adjacent_height_error = max(max_adjacent_height_error, distance.final_height_error)
        max_adjacent_command_error = max(max_adjacent_command_error, distance.command_error)
        max_adjacent_alpha_grid_error = max(
            max_adjacent_alpha_grid_error, distance.alpha_grid_error
        )
        adjacent.append(
            AdjacentSwapDiagnostic(
                index=index,
                first_kind=sequence[index].kind,
                second_kind=sequence[index + 1].kind,
                sensitive=sensitive,
                distance=distance,
            )
        )

    max_order_error = max(
        reverse.final_error,
        reverse.command_error,
        reverse.alpha_grid_error,
        max_adjacent_alpha_error,
        max_adjacent_height_error,
        max_adjacent_command_error,
        max_adjacent_alpha_grid_error,
    )
    order_sensitive = bool(
        reverse.mode_changes > 0
        or max_order_error > tolerance
        or noncommuting_adjacent_pairs > 0
    )
    return SequenceOrderDiagnostic(
        event_count=len(sequence),
        adjacent_pair_count=max(0, len(sequence) - 1),
        noncommuting_adjacent_pairs=noncommuting_adjacent_pairs,
        order_sensitive=order_sensitive,
        reverse=reverse,
        max_adjacent_alpha_error=max_adjacent_alpha_error,
        max_adjacent_height_error=max_adjacent_height_error,
        max_adjacent_command_error=max_adjacent_command_error,
        max_adjacent_alpha_grid_error=max_adjacent_alpha_grid_error,
        max_order_error=max_order_error,
        adjacent=tuple(adjacent),
    )


def compare_group_actuation_decomposition(
    config: LatticeConfig,
    state: LatticeState,
    cells: Iterable[tuple[int, int]],
    alpha: float = 0.0,
    z: float = 0.0,
    tolerance: float = 1e-9,
) -> GroupActuationDecompositionDiagnostic:
    """Compare one group command with its local-event decomposition.

    This is an algebraic simulator diagnostic: it checks that the current
    group-actuation event semantics are equivalent to sequentially applying the
    corresponding local actuation events over the same declared support.
    Removed cells are tracked separately because both semantics intentionally
    skip them.
    """

    base = state.normalized(config)
    validated = _validate_cells(config, cells)
    applied = tuple((r, c) for r, c in validated if not base.removed_mask[r, c])
    skipped = tuple((r, c) for r, c in validated if base.removed_mask[r, c])
    group_event = group_actuation_event(validated, alpha=alpha, z=z)
    local_events = tuple(local_actuation_event(cell, alpha=alpha, z=z) for cell in validated)
    simultaneous_state = apply_event_sequence(config, base, (group_event,))
    sequenced_state = apply_event_sequence(config, base, local_events)
    distance = _state_distance(config, simultaneous_state, sequenced_state)
    decomposes = bool(
        distance.mode_changes == 0
        and distance.command_error <= tolerance
        and distance.alpha_grid_error <= tolerance
        and distance.final_alpha_error <= tolerance
        and distance.final_height_error <= tolerance
    )
    return GroupActuationDecompositionDiagnostic(
        cells=validated,
        applied_cells=applied,
        skipped_removed_cells=skipped,
        group_event=group_event,
        local_events=local_events,
        simultaneous_state=simultaneous_state,
        sequenced_state=sequenced_state,
        distance=distance,
        decomposes=decomposes,
    )


def compare_group_actuation_under_removal(
    config: LatticeConfig,
    state: LatticeState,
    cells: Iterable[tuple[int, int]],
    removed_cells: Iterable[tuple[int, int]],
    alpha: float = 0.0,
    z: float = 0.0,
    tolerance: float = 1e-9,
) -> GroupRemovalInteractionDiagnostic:
    """Compare a group command before and after deleting cells.

    This exposes the graph/topology effect of cell removal on a declared group
    support. It separates cells that are no longer actuatable from response
    changes caused by deleted couplings between still-present cells.
    """

    base = state.normalized(config)
    validated_cells = _validate_cells(config, cells)
    validated_removed = _validate_cells(config, removed_cells)
    removal_events = tuple(remove_cell_event(cell) for cell in validated_removed)
    intact = compare_group_actuation_decomposition(
        config, base, validated_cells, alpha=alpha, z=z, tolerance=tolerance
    )
    removed_base = apply_event_sequence(config, base, removal_events)
    removed = compare_group_actuation_decomposition(
        config, removed_base, validated_cells, alpha=alpha, z=z, tolerance=tolerance
    )
    intact_applied = set(intact.applied_cells)
    removed_applied = set(removed.applied_cells)
    intact_skipped = set(intact.skipped_removed_cells)
    removed_skipped = set(removed.skipped_removed_cells)
    intact_topology = lattice_topology_diagnostic(config, base)
    removed_topology = lattice_topology_diagnostic(config, removed_base)
    distance = _state_distance(
        config, intact.simultaneous_state, removed.simultaneous_state
    )
    return GroupRemovalInteractionDiagnostic(
        cells=validated_cells,
        removed_cells=validated_removed,
        removal_events=removal_events,
        intact=intact,
        removed=removed,
        lost_applied_cells=tuple(cell for cell in validated_cells if cell in intact_applied and cell not in removed_applied),
        newly_skipped_cells=tuple(cell for cell in validated_cells if cell in removed_skipped and cell not in intact_skipped),
        intact_topology=intact_topology,
        removed_topology=removed_topology,
        distance=distance,
        component_count_delta=int(removed_topology["component_count"])
        - int(intact_topology["component_count"]),
        deleted_edge_delta=int(removed_topology["deleted_edge_count"])
        - int(intact_topology["deleted_edge_count"]),
    )


def _cells_above_threshold(
    values: np.ndarray,
    threshold: float,
    exclude: set[tuple[int, int]] | None = None,
) -> tuple[tuple[int, int], ...]:
    excluded = exclude or set()
    rows, cols = values.shape
    return tuple(
        (r, c)
        for r in range(rows)
        for c in range(cols)
        if (r, c) not in excluded and abs(float(values[r, c])) > threshold
    )


def _vertical_contact_penalty(
    config: LatticeConfig,
    state: LatticeState,
    cells: tuple[tuple[int, int], ...],
    contact_stiffness: float,
    tolerance: float,
) -> tuple[float, tuple[tuple[int, int], ...]]:
    energy = 0.0
    engaged: list[tuple[int, int]] = []
    for r, c in cells:
        if state.removed_mask[r, c]:
            continue
        penetration = max(0.0, abs(float(state.z_actuator_grid[r, c])) - config.pin_hole_clearance)
        if penetration > tolerance:
            engaged.append((r, c))
        energy += 0.5 * contact_stiffness * penetration**2
    return float(energy), tuple(engaged)


def _load_active_cells(
    config: LatticeConfig,
    state: LatticeState,
    fixed_cells: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int], ...]:
    fixed = set(fixed_cells)
    return tuple(
        (r, c)
        for r in range(config.rows)
        for c in range(config.cols)
        if (r, c) not in fixed
        and not bool(state.removed_mask[r, c])
        and not bool(state.locked_mask[r, c])
    )


def _height_load_contact_metrics(
    config: LatticeConfig,
    state: LatticeState,
    height: np.ndarray,
    fixed_cells: tuple[tuple[int, int], ...],
    external_z_load: float,
    contact_stiffness: float,
    tolerance: float,
) -> tuple[
    tuple[tuple[int, int], ...],
    float,
    float,
    float,
    tuple[tuple[int, int], ...],
]:
    active = _load_active_cells(config, state, fixed_cells)
    work = 0.0
    work_magnitude = 0.0
    contact_energy = 0.0
    engaged: list[tuple[int, int]] = []
    for cell in active:
        displacement = float(height[cell])
        work += -external_z_load * displacement
        work_magnitude += abs(external_z_load) * abs(displacement)
        penetration = max(0.0, abs(displacement) - config.pin_hole_clearance)
        if penetration > tolerance:
            engaged.append(cell)
        contact_energy += 0.5 * contact_stiffness * penetration**2
    return (
        active,
        float(work),
        float(work_magnitude),
        float(contact_energy),
        tuple(engaged),
    )


def compare_vertical_residual_under_removal(
    config: LatticeConfig,
    state: LatticeState,
    cells: Iterable[tuple[int, int]],
    removed_cells: Iterable[tuple[int, int]],
    z: float,
    alpha: float = 0.0,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    fixed_cells: Iterable[tuple[int, int]] = (),
    external_z_load: float = 0.0,
) -> VerticalRemovalResidualDiagnostic:
    """Compare clearance-gated vertical residuals before and after cell removal.

    The contact penalty is an uncalibrated normalized proxy for contact onset:
    `0.5 * k * max(|z_command| - pin_hole_clearance, 0)^2` over the declared
    group support. Load work and height contact are also normalized diagnostics:
    a positive `external_z_load` means upward force and a negative value means
    downward gravity-like force. These are not validated RAD material laws.
    """

    if contact_stiffness < 0:
        raise ValueError("contact_stiffness must be non-negative")
    base = state.normalized(config)
    validated_cells = _validate_cells(config, cells)
    validated_removed = _validate_cells(config, removed_cells)
    validated_fixed = _validate_cells(config, fixed_cells)
    source_set = set(validated_cells)
    removal_events = tuple(remove_cell_event(cell) for cell in validated_removed)
    group_event = group_actuation_event(validated_cells, alpha=alpha, z=z)

    intact_state = apply_event_sequence(config, base, (group_event,))
    removed_base = apply_event_sequence(config, base, removal_events)
    removed_state = apply_event_sequence(config, removed_base, (group_event,))
    intact = evaluate_programmable_operators(config, intact_state)
    removed = evaluate_programmable_operators(config, removed_state)
    intact_z = np.asarray(intact["z_residual"], dtype=float)
    removed_z = np.asarray(removed["z_residual"], dtype=float)
    residual_delta = removed_z - intact_z
    intact_reach = _cells_above_threshold(intact_z, tolerance)
    removed_reach = _cells_above_threshold(removed_z, tolerance)
    intact_neighbors = _cells_above_threshold(intact_z, tolerance, exclude=source_set)
    removed_neighbors = _cells_above_threshold(removed_z, tolerance, exclude=source_set)
    removed_mask = removed_state.removed_mask
    lost_residual_cells = tuple(
        cell for cell in intact_reach if abs(float(removed_z[cell])) <= tolerance
    )
    topology_blocked_cells = tuple(
        cell for cell in lost_residual_cells if not bool(removed_mask[cell])
    )
    intact_energy, intact_engaged = _vertical_contact_penalty(
        config, intact_state, validated_cells, contact_stiffness, tolerance
    )
    removed_energy, removed_engaged = _vertical_contact_penalty(
        config, removed_state, validated_cells, contact_stiffness, tolerance
    )
    (
        intact_load_cells,
        intact_load_work,
        intact_load_work_magnitude,
        intact_height_contact_energy,
        intact_height_contact_engaged,
    ) = _height_load_contact_metrics(
        config,
        intact_state,
        np.asarray(intact["height"], dtype=float),
        validated_fixed,
        float(external_z_load),
        contact_stiffness,
        tolerance,
    )
    (
        removed_load_cells,
        removed_load_work,
        removed_load_work_magnitude,
        removed_height_contact_energy,
        removed_height_contact_engaged,
    ) = _height_load_contact_metrics(
        config,
        removed_state,
        np.asarray(removed["height"], dtype=float),
        validated_fixed,
        float(external_z_load),
        contact_stiffness,
        tolerance,
    )
    intact_topology = lattice_topology_diagnostic(config, base)
    removed_topology = lattice_topology_diagnostic(config, removed_base)
    neighbor_values = [
        abs(float(removed_z[cell]))
        for cell in removed_neighbors
    ]
    return VerticalRemovalResidualDiagnostic(
        cells=validated_cells,
        removed_cells=validated_removed,
        fixed_cells=validated_fixed,
        removal_events=removal_events,
        intact_state=intact_state,
        removed_state=removed_state,
        intact_topology=intact_topology,
        removed_topology=removed_topology,
        intact_z_residual=intact_z,
        removed_z_residual=removed_z,
        residual_delta=residual_delta,
        intact_z_die_off=np.asarray(intact["z_die_off"], dtype=float),
        removed_z_die_off=np.asarray(removed["z_die_off"], dtype=float),
        affected_neighbor_cells=removed_neighbors,
        lost_residual_cells=lost_residual_cells,
        topology_blocked_cells=topology_blocked_cells,
        removed_source_cells=tuple(cell for cell in validated_cells if removed_mask[cell]),
        intact_reach_count=len(intact_reach),
        removed_reach_count=len(removed_reach),
        intact_neighbor_reach_count=len(intact_neighbors),
        removed_neighbor_reach_count=len(removed_neighbors),
        intact_die_off_radius=finite_die_off_radius(intact["z_die_off"]),
        removed_die_off_radius=finite_die_off_radius(removed["z_die_off"]),
        max_abs_neighbor_residual=max(neighbor_values, default=0.0),
        max_abs_residual_delta=float(np.max(np.abs(residual_delta))),
        clearance=config.pin_hole_clearance,
        contact_stiffness=float(contact_stiffness),
        intact_contact_engaged_cells=intact_engaged,
        removed_contact_engaged_cells=removed_engaged,
        intact_contact_penalty_energy=intact_energy,
        removed_contact_penalty_energy=removed_energy,
        contact_penalty_delta=float(removed_energy - intact_energy),
        external_z_load=float(external_z_load),
        intact_load_active_cells=intact_load_cells,
        removed_load_active_cells=removed_load_cells,
        intact_load_work=intact_load_work,
        removed_load_work=removed_load_work,
        load_work_delta=float(removed_load_work - intact_load_work),
        intact_load_work_magnitude=intact_load_work_magnitude,
        removed_load_work_magnitude=removed_load_work_magnitude,
        load_work_magnitude_delta=float(
            removed_load_work_magnitude - intact_load_work_magnitude
        ),
        intact_height_contact_engaged_cells=intact_height_contact_engaged,
        removed_height_contact_engaged_cells=removed_height_contact_engaged,
        intact_height_contact_penalty_energy=intact_height_contact_energy,
        removed_height_contact_penalty_energy=removed_height_contact_energy,
        height_contact_penalty_delta=float(
            removed_height_contact_energy - intact_height_contact_energy
        ),
        component_count_delta=int(removed_topology["component_count"])
        - int(intact_topology["component_count"]),
        deleted_edge_delta=int(removed_topology["deleted_edge_count"])
        - int(intact_topology["deleted_edge_count"]),
    )


def _array_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(arr**2)))


def _anchored_load_cells(
    topology: dict[str, Any],
    load_cells: tuple[tuple[int, int], ...],
    fixed_cells: tuple[tuple[int, int], ...],
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    labels = np.asarray(topology["component_labels"], dtype=int)
    fixed_components = {
        int(labels[r, c])
        for r, c in fixed_cells
        if 0 <= int(labels[r, c])
    }
    anchored = tuple(
        cell for cell in load_cells if int(labels[cell]) in fixed_components
    )
    unanchored = tuple(cell for cell in load_cells if cell not in set(anchored))
    return anchored, unanchored


def _add_external_z_load(
    existing: tuple[float, ...] | np.ndarray | None, external_z_load: float
) -> tuple[float, float, float]:
    vector = np.zeros(3, dtype=float)
    if existing is not None:
        raw = np.asarray(existing, dtype=float)
        vector[: min(3, raw.size)] = raw[: min(3, raw.size)]
    vector[2] += float(external_z_load)
    return tuple(float(value) for value in vector)


def _load_case_with_solver_forces(
    base: LoadCase,
    fixed_cells: tuple[tuple[int, int], ...],
    load_cells: tuple[tuple[int, int], ...],
    external_z_load: float,
) -> LoadCase:
    external_forces = dict(base.external_forces)
    for cell in load_cells:
        external_forces[cell] = _add_external_z_load(
            external_forces.get(cell), external_z_load
        )
    return LoadCase(
        fixed_cells=fixed_cells,
        prescribed_displacements=dict(base.prescribed_displacements),
        external_forces=external_forces,
        axial_stiffness=base.axial_stiffness,
        hinge_stiffness=base.hinge_stiffness,
        lock_stiffness=base.lock_stiffness,
        maxiter=base.maxiter,
    )


def _solver_energy_breakdown_from_metadata(
    metadata: dict[str, Any],
) -> SolverEnergyBreakdown:
    return SolverEnergyBreakdown(
        objective_energy=float(metadata.get("objective_energy", metadata["energy"])),
        stored_energy=float(metadata.get("stored_energy", metadata["energy"])),
        spring_energy=float(metadata.get("spring_energy", 0.0)),
        hinge_energy=float(metadata.get("hinge_energy", 0.0)),
        lock_penalty_energy=float(metadata.get("lock_penalty_energy", 0.0)),
        external_potential_energy=float(metadata.get("external_potential_energy", 0.0)),
    )


def compare_vertical_residual_spring_hinge_3d(
    config: LatticeConfig,
    state: LatticeState,
    cells: Iterable[tuple[int, int]],
    removed_cells: Iterable[tuple[int, int]],
    z: float,
    alpha: float = 0.0,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    fixed_cells: Iterable[tuple[int, int]] | None = None,
    external_z_load: float = 0.0,
    load_case: LoadCase | None = None,
) -> VerticalRemovalPhysicalComparisonDiagnostic:
    """Compare kinematic vertical-removal diagnostics with 3D spring-hinge output.

    This is a physical-preview diagnostic. It uses the same removed-cell graph
    for the spring/hinge energy, applies z-forces only to components connected
    to fixed cells, and reports unanchored load cells separately.
    """

    base_load_case = load_case or LoadCase()
    selected_fixed = (
        base_load_case.fixed_cells if fixed_cells is None else tuple(fixed_cells)
    )
    diagnostic = compare_vertical_residual_under_removal(
        config,
        state,
        cells,
        removed_cells,
        z,
        alpha=alpha,
        tolerance=tolerance,
        contact_stiffness=contact_stiffness,
        fixed_cells=selected_fixed,
        external_z_load=external_z_load,
    )
    fixed = diagnostic.fixed_cells
    intact_load_cells, intact_unanchored = _anchored_load_cells(
        diagnostic.intact_topology,
        diagnostic.intact_load_active_cells,
        fixed,
    )
    removed_load_cells, removed_unanchored = _anchored_load_cells(
        diagnostic.removed_topology,
        diagnostic.removed_load_active_cells,
        fixed,
    )
    intact_load_case = _load_case_with_solver_forces(
        base_load_case,
        fixed,
        intact_load_cells,
        external_z_load,
    )
    removed_load_case = _load_case_with_solver_forces(
        base_load_case,
        fixed,
        removed_load_cells,
        external_z_load,
    )

    from .spring_hinge import solve_spring_hinge_3d

    intact_result = solve_spring_hinge_3d(
        config, diagnostic.intact_state, intact_load_case
    )
    removed_result = solve_spring_hinge_3d(
        config, diagnostic.removed_state, removed_load_case
    )
    intact_physical = np.asarray(intact_result.metadata["height"], dtype=float)
    removed_physical = np.asarray(removed_result.metadata["height"], dtype=float)
    intact_kinematic = np.asarray(
        intact_result.metadata["kinematic_height"], dtype=float
    )
    removed_kinematic = np.asarray(
        removed_result.metadata["kinematic_height"], dtype=float
    )
    intact_error = intact_physical - intact_kinematic
    removed_error = removed_physical - removed_kinematic
    physical_delta = removed_physical - intact_physical
    kinematic_delta = removed_kinematic - intact_kinematic
    delta_error = physical_delta - kinematic_delta
    return VerticalRemovalPhysicalComparisonDiagnostic(
        vertical_diagnostic=diagnostic,
        intact_load_case=intact_load_case,
        removed_load_case=removed_load_case,
        intact_solver_success=bool(intact_result.metadata["success"]),
        removed_solver_success=bool(removed_result.metadata["success"]),
        physical_success=bool(
            intact_result.metadata["success"] and removed_result.metadata["success"]
        ),
        intact_solver_energy=float(intact_result.metadata["energy"]),
        removed_solver_energy=float(removed_result.metadata["energy"]),
        intact_energy_breakdown=_solver_energy_breakdown_from_metadata(
            intact_result.metadata
        ),
        removed_energy_breakdown=_solver_energy_breakdown_from_metadata(
            removed_result.metadata
        ),
        intact_solver_iterations=int(intact_result.metadata["iterations"]),
        removed_solver_iterations=int(removed_result.metadata["iterations"]),
        intact_spring_edges=int(intact_result.metadata["spring_edges"]),
        removed_spring_edges=int(removed_result.metadata["spring_edges"]),
        spring_edge_delta=int(removed_result.metadata["spring_edges"])
        - int(intact_result.metadata["spring_edges"]),
        intact_hinge_triples=int(intact_result.metadata["hinge_triples"]),
        removed_hinge_triples=int(removed_result.metadata["hinge_triples"]),
        hinge_triple_delta=int(removed_result.metadata["hinge_triples"])
        - int(intact_result.metadata["hinge_triples"]),
        intact_solver_removed_cells=int(intact_result.metadata["removed_cells"]),
        removed_solver_removed_cells=int(removed_result.metadata["removed_cells"]),
        intact_solver_load_cells=intact_load_cells,
        removed_solver_load_cells=removed_load_cells,
        intact_unanchored_load_cells=intact_unanchored,
        removed_unanchored_load_cells=removed_unanchored,
        intact_physical_height=intact_physical,
        removed_physical_height=removed_physical,
        intact_kinematic_height=intact_kinematic,
        removed_kinematic_height=removed_kinematic,
        intact_model_error=intact_error,
        removed_model_error=removed_error,
        physical_height_delta=physical_delta,
        kinematic_height_delta=kinematic_delta,
        physical_vs_kinematic_delta_error=delta_error,
        intact_height_rms_model_error=_array_rms(intact_error),
        removed_height_rms_model_error=_array_rms(removed_error),
        delta_rms_model_error=_array_rms(delta_error),
        max_abs_delta_model_error=float(np.max(np.abs(delta_error))),
        solver_topology_note=(
            "spring_hinge_3d skips spring and hinge terms incident to removed "
            "cells; unanchored load components are reported and excluded from "
            "quasistatic external forcing"
        ),
    )


def _cell_dict(cell: tuple[int, int]) -> dict[str, int]:
    return {"row": int(cell[0]), "col": int(cell[1])}


def _cells_dict(cells: Iterable[tuple[int, int]]) -> list[dict[str, int]]:
    return [_cell_dict(cell) for cell in cells]


def _json_number(value: float) -> float | None:
    numeric = float(value)
    return numeric if np.isfinite(numeric) else None


def _matrix_json(values: np.ndarray) -> list[list[float | None]]:
    arr = np.asarray(values, dtype=float)
    return [[_json_number(value) for value in row] for row in arr]


def _load_case_dict(load_case: LoadCase) -> dict[str, object]:
    return {
        "fixedCells": _cells_dict(load_case.fixed_cells),
        "prescribedDisplacements": [
            {"cell": _cell_dict(cell), "displacement": [float(v) for v in values]}
            for cell, values in load_case.prescribed_displacements.items()
        ],
        "externalForces": [
            {"cell": _cell_dict(cell), "force": [float(v) for v in values]}
            for cell, values in load_case.external_forces.items()
        ],
        "axialStiffness": float(load_case.axial_stiffness),
        "hingeStiffness": float(load_case.hinge_stiffness),
        "lockStiffness": float(load_case.lock_stiffness),
        "maxiter": int(load_case.maxiter),
    }


def _energy_breakdown_dict(breakdown: SolverEnergyBreakdown) -> dict[str, float]:
    return {
        "objectiveEnergy": breakdown.objective_energy,
        "storedEnergy": breakdown.stored_energy,
        "springEnergy": breakdown.spring_energy,
        "hingeEnergy": breakdown.hinge_energy,
        "lockPenaltyEnergy": breakdown.lock_penalty_energy,
        "externalPotentialEnergy": breakdown.external_potential_energy,
    }


def mechanics_energy_certificate_to_dict(
    comparison: VerticalRemovalPhysicalComparisonDiagnostic,
    tolerance: float = 1e-8,
) -> dict[str, object]:
    """Return a claim-labeled energy certificate for a physical-preview result."""

    vertical = comparison.vertical_diagnostic

    def state_terms(
        breakdown: SolverEnergyBreakdown,
        clearance_contact: float,
        height_contact: float,
        load_work: float,
        load_work_magnitude: float,
    ) -> dict[str, object]:
        nonnegative_terms = {
            "springEnergy": breakdown.spring_energy,
            "hingeEnergy": breakdown.hinge_energy,
            "lockPenaltyEnergy": breakdown.lock_penalty_energy,
            "clearanceContactPenalty": clearance_contact,
            "heightContactPenalty": height_contact,
            "loadWorkMagnitude": load_work_magnitude,
        }
        signed_terms = {
            "solverObjectiveEnergy": breakdown.objective_energy,
            "externalPotentialEnergy": breakdown.external_potential_energy,
            "signedLoadWork": load_work,
        }
        max_violation = max(
            (max(0.0, -float(value)) for value in nonnegative_terms.values()),
            default=0.0,
        )
        return {
            "nonnegativeTerms": nonnegative_terms,
            "signedTerms": signed_terms,
            "nonnegativeProxyTotal": float(
                sum(float(value) for value in nonnegative_terms.values())
            ),
            "storedEnergy": breakdown.stored_energy,
            "maxNonnegativeViolation": float(max_violation),
            "passesNonnegativeCheck": bool(max_violation <= tolerance),
            "objectiveConsistencyError": abs(
                breakdown.objective_energy
                - (breakdown.stored_energy + breakdown.external_potential_energy)
            ),
        }

    intact = state_terms(
        comparison.intact_energy_breakdown,
        vertical.intact_contact_penalty_energy,
        vertical.intact_height_contact_penalty_energy,
        vertical.intact_load_work,
        vertical.intact_load_work_magnitude,
    )
    removed = state_terms(
        comparison.removed_energy_breakdown,
        vertical.removed_contact_penalty_energy,
        vertical.removed_height_contact_penalty_energy,
        vertical.removed_load_work,
        vertical.removed_load_work_magnitude,
    )
    return {
        "schema": "rad-sim.mechanics-energy-certificate.v1",
        "method": (
            "separates nonnegative stored/proxy terms from signed external "
            "work terms in the uncalibrated spring-hinge physical preview"
        ),
        "claimLabels": {
            "nonnegativeTerms": "Lean-proven theorem scaffold",
            "solverEnergyBreakdown": "simulator-derived empirical law",
            "signedExternalWork": "experimentally unvalidated physical assumption",
        },
        "leanTheorems": [
            "springEnergyNat_nonnegative",
            "hingeEnergyNat_nonnegative",
            "contactPenaltyNat_nonnegative",
            "loadWorkMagnitudeNat_nonnegative",
            "mechanicalStoredEnergyNat_nonnegative",
            "mechanicalStoredEnergyNat_zero_components",
        ],
        "tolerance": float(tolerance),
        "intact": intact,
        "removed": removed,
        "delta": {
            "nonnegativeProxyTotal": float(
                removed["nonnegativeProxyTotal"] - intact["nonnegativeProxyTotal"]
            ),
            "storedEnergy": (
                comparison.removed_energy_breakdown.stored_energy
                - comparison.intact_energy_breakdown.stored_energy
            ),
            "signedLoadWork": vertical.load_work_delta,
            "loadWorkMagnitude": vertical.load_work_magnitude_delta,
            "clearanceContactPenalty": vertical.contact_penalty_delta,
            "heightContactPenalty": vertical.height_contact_penalty_delta,
        },
        "passesNonnegativeCheck": bool(
            intact["passesNonnegativeCheck"] and removed["passesNonnegativeCheck"]
        ),
        "claimLimit": (
            "This certifies algebraic sign conventions for the normalized "
            "preview terms; it is not calibrated rigid-body contact validation."
        ),
    }


def _normalize_height_measurements(
    heights: Mapping[tuple[int, int], float],
) -> dict[tuple[int, int], float]:
    return {
        (int(cell[0]), int(cell[1])): float(value)
        for cell, value in heights.items()
    }


def _work_contact_metrics_from_heights(
    cells: Iterable[tuple[int, int]],
    heights: Mapping[tuple[int, int], float],
    external_z_load: float,
    clearance: float,
    contact_stiffness: float,
    tolerance: float,
) -> dict[str, object]:
    observed: list[tuple[int, int]] = []
    missing: list[tuple[int, int]] = []
    engaged: list[tuple[int, int]] = []
    signed_work = 0.0
    work_magnitude = 0.0
    contact_penalty = 0.0
    max_abs_height = 0.0
    for cell in cells:
        if cell not in heights:
            missing.append(cell)
            continue
        height = float(heights[cell])
        observed.append(cell)
        signed_work += -float(external_z_load) * height
        work_magnitude += abs(float(external_z_load)) * abs(height)
        penetration = max(0.0, abs(height) - float(clearance))
        if penetration > tolerance:
            engaged.append(cell)
        contact_penalty += 0.5 * float(contact_stiffness) * penetration**2
        max_abs_height = max(max_abs_height, abs(height))
    return {
        "observedCells": _cells_dict(observed),
        "missingCells": _cells_dict(missing),
        "observedCellCount": len(observed),
        "missingCellCount": len(missing),
        "signedLoadWork": float(signed_work),
        "loadWorkMagnitude": float(work_magnitude),
        "heightContactPenalty": float(contact_penalty),
        "contactEngagedCells": _cells_dict(engaged),
        "maxAbsHeight": float(max_abs_height),
    }


def _height_map_from_array(
    cells: Iterable[tuple[int, int]],
    height: np.ndarray,
) -> dict[tuple[int, int], float]:
    arr = np.asarray(height, dtype=float)
    return {cell: float(arr[cell]) for cell in cells}


def _energy_validation_state(
    load_cells: tuple[tuple[int, int], ...],
    simulated_height: np.ndarray,
    measured_height: Mapping[tuple[int, int], float],
    external_z_load: float,
    clearance: float,
    contact_stiffness: float,
    tolerance: float,
) -> dict[str, object]:
    measurements = _normalize_height_measurements(measured_height)
    observed_cells = tuple(cell for cell in load_cells if cell in measurements)
    simulated_observed = _work_contact_metrics_from_heights(
        observed_cells,
        _height_map_from_array(observed_cells, simulated_height),
        external_z_load,
        clearance,
        contact_stiffness,
        tolerance,
    )
    simulated_all = _work_contact_metrics_from_heights(
        load_cells,
        _height_map_from_array(load_cells, simulated_height),
        external_z_load,
        clearance,
        contact_stiffness,
        tolerance,
    )
    measured = _work_contact_metrics_from_heights(
        load_cells,
        measurements,
        external_z_load,
        clearance,
        contact_stiffness,
        tolerance,
    )
    errors = {
        "signedLoadWork": float(
            measured["signedLoadWork"] - simulated_observed["signedLoadWork"]
        ),
        "loadWorkMagnitude": float(
            measured["loadWorkMagnitude"] - simulated_observed["loadWorkMagnitude"]
        ),
        "heightContactPenalty": float(
            measured["heightContactPenalty"]
            - simulated_observed["heightContactPenalty"]
        ),
        "maxAbsHeight": float(
            measured["maxAbsHeight"] - simulated_observed["maxAbsHeight"]
        ),
    }
    return {
        "loadCells": _cells_dict(load_cells),
        "measured": measured,
        "simulatedObserved": simulated_observed,
        "simulatedAllLoaded": simulated_all,
        "errors": errors,
        "allLoadCellsMeasured": measured["missingCellCount"] == 0,
        "passesTolerance": bool(
            measured["missingCellCount"] == 0
            and all(abs(value) <= tolerance for value in errors.values())
        ),
    }


def validate_vertical_load_energy_measurements(
    comparison: VerticalRemovalPhysicalComparisonDiagnostic,
    intact_heights: Mapping[tuple[int, int], float],
    removed_heights: Mapping[tuple[int, int], float],
    tolerance: float = 1e-8,
    contact_stiffness: float | None = None,
) -> dict[str, object]:
    """Compare bench-measured load-cell heights against preview work/contact terms."""

    vertical = comparison.vertical_diagnostic
    stiffness = (
        vertical.contact_stiffness if contact_stiffness is None else contact_stiffness
    )
    intact = _energy_validation_state(
        comparison.intact_solver_load_cells,
        comparison.intact_physical_height,
        intact_heights,
        vertical.external_z_load,
        vertical.clearance,
        stiffness,
        tolerance,
    )
    removed = _energy_validation_state(
        comparison.removed_solver_load_cells,
        comparison.removed_physical_height,
        removed_heights,
        vertical.external_z_load,
        vertical.clearance,
        stiffness,
        tolerance,
    )
    error_values = [
        abs(float(state["errors"][key]))
        for state in (intact, removed)
        for key in ("signedLoadWork", "loadWorkMagnitude", "heightContactPenalty")
    ]
    return {
        "schema": "rad-sim.vertical-load-energy-validation.v1",
        "method": (
            "compares measured load-cell heights against spring-hinge predicted "
            "signed load work, load-work magnitude, and height-contact proxy terms"
        ),
        "claimLabels": {
            "zeroResidualIdentity": "Lean-proven theorem scaffold",
            "springHingePrediction": "simulator-derived empirical law",
            "benchMeasurement": "experimentally unvalidated physical assumption",
        },
        "leanTheorems": [
            "absoluteErrorNat_self",
            "measuredWorkResidualNat_zero_when_equal",
            "loadWorkMagnitudeNat_nonnegative",
        ],
        "tolerance": float(tolerance),
        "externalZLoad": vertical.external_z_load,
        "clearance": vertical.clearance,
        "contactStiffness": float(stiffness),
        "intact": intact,
        "removed": removed,
        "summary": {
            "requiredCellCount": (
                len(comparison.intact_solver_load_cells)
                + len(comparison.removed_solver_load_cells)
            ),
            "measuredCellCount": int(
                intact["measured"]["observedCellCount"]
                + removed["measured"]["observedCellCount"]
            ),
            "missingMeasurementCount": int(
                intact["measured"]["missingCellCount"]
                + removed["measured"]["missingCellCount"]
            ),
            "maxAbsWorkOrContactError": max(error_values, default=0.0),
            "allLoadCellsMeasured": bool(
                intact["allLoadCellsMeasured"] and removed["allLoadCellsMeasured"]
            ),
            "passesTolerance": bool(
                intact["passesTolerance"] and removed["passesTolerance"]
            ),
        },
        "claimLimit": (
            "A zero residual only proves agreement with the current normalized "
            "preview equations for the supplied measurements; it does not prove "
            "calibrated gravity, friction, or rigid-body contact."
        ),
    }


def export_vertical_load_energy_validation_json(
    validation: dict[str, object],
) -> str:
    return json.dumps(validation, indent=2)


def _cell_from_dict(payload: object) -> tuple[int, int]:
    if not isinstance(payload, dict):
        raise ValueError("cell payload must be an object with row and col")
    return (int(payload["row"]), int(payload["col"]))


def _cells_from_dicts(payload: object) -> tuple[tuple[int, int], ...]:
    if not isinstance(payload, list):
        return ()
    return tuple(_cell_from_dict(cell) for cell in payload)


def _scenario_input_dict(scenario: VerticalLoadPhysicalScenario) -> dict[str, object]:
    return {
        "cells": _cells_dict(scenario.cells),
        "removedCells": _cells_dict(scenario.removed_cells),
        "fixedCells": _cells_dict(scenario.fixed_cells),
        "alpha": scenario.alpha,
        "z": scenario.z,
        "externalZLoad": scenario.external_z_load,
    }


def _scenario_from_result_payload(
    scenario_payload: dict[str, object],
) -> VerticalLoadPhysicalScenario:
    input_payload = scenario_payload.get("input", {})
    if not isinstance(input_payload, dict):
        raise ValueError("vertical-load energy scenario is missing input")
    return VerticalLoadPhysicalScenario(
        name=str(scenario_payload.get("name", "unnamed-scenario")),
        cells=_cells_from_dicts(input_payload.get("cells", [])),
        removed_cells=_cells_from_dicts(input_payload.get("removedCells", [])),
        fixed_cells=_cells_from_dicts(input_payload.get("fixedCells", [])),
        alpha=float(input_payload.get("alpha", 0.0)),
        z=float(input_payload.get("z", 0.0)),
        external_z_load=float(input_payload.get("externalZLoad", 0.0)),
    )


def _measurement_template_cells(
    cells: Iterable[tuple[int, int]],
    height: np.ndarray,
) -> list[dict[str, object]]:
    arr = np.asarray(height, dtype=float)
    return [
        {
            "cell": _cell_dict(cell),
            "measuredHeight": None,
            "predictedHeight": float(arr[cell]),
            "notes": "",
        }
        for cell in cells
    ]


def _measurement_heights_from_cells(
    cells_payload: object,
) -> dict[tuple[int, int], float]:
    if not isinstance(cells_payload, list):
        return {}
    out: dict[tuple[int, int], float] = {}
    for item in cells_payload:
        if not isinstance(item, dict):
            continue
        cell_payload = item.get("cell", item)
        cell = _cell_from_dict(cell_payload)
        value = item.get("measuredHeight", item.get("height", None))
        if value is None:
            continue
        out[cell] = float(value)
    return out


def vertical_load_energy_measurement_template(
    config: LatticeConfig,
    state: LatticeState | None = None,
    scenarios: Iterable[VerticalLoadPhysicalScenario] | None = None,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    load_case: LoadCase | None = None,
) -> dict[str, object]:
    """Create a fillable bench-results template for vertical-load energy tests."""

    base_state = LatticeState.uniform(config) if state is None else state.normalized(config)
    scenario_tuple = (
        default_vertical_load_physical_scenarios(config)
        if scenarios is None
        else tuple(scenarios)
    )
    scenario_payloads: list[dict[str, object]] = []
    for scenario in scenario_tuple:
        comparison = compare_vertical_residual_spring_hinge_3d(
            config,
            base_state,
            scenario.cells,
            scenario.removed_cells,
            scenario.z,
            alpha=scenario.alpha,
            tolerance=tolerance,
            contact_stiffness=contact_stiffness,
            fixed_cells=scenario.fixed_cells,
            external_z_load=scenario.external_z_load,
            load_case=load_case,
        )
        scenario_payloads.append(
            {
                "name": scenario.name,
                "input": _scenario_input_dict(scenario),
                "requiredMeasurements": {
                    "intactLoadCells": _cells_dict(comparison.intact_solver_load_cells),
                    "removedLoadCells": _cells_dict(comparison.removed_solver_load_cells),
                    "unanchoredLoadCells": {
                        "intact": _cells_dict(comparison.intact_unanchored_load_cells),
                        "removed": _cells_dict(comparison.removed_unanchored_load_cells),
                    },
                },
                "measurements": {
                    "intact": _measurement_template_cells(
                        comparison.intact_solver_load_cells,
                        comparison.intact_physical_height,
                    ),
                    "removed": _measurement_template_cells(
                        comparison.removed_solver_load_cells,
                        comparison.removed_physical_height,
                    ),
                },
                "notes": "",
            }
        )
    return {
        "schema": "rad-sim.vertical-load-energy-measurement-results.v1",
        "templateSchema": "rad-sim.vertical-load-energy-measurement-template.v1",
        "method": (
            "fill measuredHeight for each listed load cell, then compare the "
            "bench measurements against vertical-load energy validation"
        ),
        "claimLabels": {
            "templateGeneration": "simulator-derived empirical law",
            "benchMeasurement": "experimentally unvalidated physical assumption",
        },
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "pinHoleClearance": config.pin_hole_clearance,
            "zCouplingGain": config.z_coupling_gain,
        },
        "contactStiffness": float(contact_stiffness),
        "tolerance": float(tolerance),
        "measurementUnits": "normalized simulator height units unless calibrated externally",
        "instructions": (
            "Do not copy predictedHeight into measuredHeight except for a simulator "
            "roundtrip smoke test. For bench use, measure each listed load-cell "
            "height after the fixed cells, removed cells, z command, and external "
            "load in the input block are applied."
        ),
        "scenarios": scenario_payloads,
    }


def export_vertical_load_energy_measurement_template_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    scenarios: Iterable[VerticalLoadPhysicalScenario] | None = None,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    load_case: LoadCase | None = None,
) -> str:
    return json.dumps(
        vertical_load_energy_measurement_template(
            config,
            state=state,
            scenarios=scenarios,
            tolerance=tolerance,
            contact_stiffness=contact_stiffness,
            load_case=load_case,
        ),
        indent=2,
    )


def vertical_load_energy_measurement_results_from_json(text: str) -> dict[str, object]:
    payload = json.loads(text)
    if not isinstance(payload, dict):
        raise ValueError("vertical-load energy measurement results must be a JSON object")
    if payload.get("schema") != "rad-sim.vertical-load-energy-measurement-results.v1":
        raise ValueError("unsupported vertical-load energy measurement results schema")
    return payload


def compare_vertical_load_energy_measurement_results(
    config: LatticeConfig,
    results: dict[str, object],
    state: LatticeState | None = None,
    tolerance: float | None = None,
    contact_stiffness: float | None = None,
    load_case: LoadCase | None = None,
) -> dict[str, object]:
    """Compare filled vertical-load bench-result rows against the simulator."""

    if results.get("schema") != "rad-sim.vertical-load-energy-measurement-results.v1":
        raise ValueError("unsupported vertical-load energy measurement results schema")
    base_state = LatticeState.uniform(config) if state is None else state.normalized(config)
    tol = float(tolerance if tolerance is not None else results.get("tolerance", 1e-8))
    stiffness = float(
        contact_stiffness
        if contact_stiffness is not None
        else results.get("contactStiffness", 1.0)
    )
    scenario_reports: list[dict[str, object]] = []
    for scenario_payload in results.get("scenarios", []):
        if not isinstance(scenario_payload, dict):
            continue
        scenario = _scenario_from_result_payload(scenario_payload)
        comparison = compare_vertical_residual_spring_hinge_3d(
            config,
            base_state,
            scenario.cells,
            scenario.removed_cells,
            scenario.z,
            alpha=scenario.alpha,
            tolerance=tol,
            contact_stiffness=stiffness,
            fixed_cells=scenario.fixed_cells,
            external_z_load=scenario.external_z_load,
            load_case=load_case,
        )
        measurements = scenario_payload.get("measurements", {})
        if not isinstance(measurements, dict):
            measurements = {}
        validation = validate_vertical_load_energy_measurements(
            comparison,
            _measurement_heights_from_cells(measurements.get("intact", [])),
            _measurement_heights_from_cells(measurements.get("removed", [])),
            tolerance=tol,
            contact_stiffness=stiffness,
        )
        scenario_reports.append(
            {
                "name": scenario.name,
                "input": _scenario_input_dict(scenario),
                "validation": validation,
            }
        )
    missing = sum(
        int(item["validation"]["summary"]["missingMeasurementCount"])
        for item in scenario_reports
    )
    max_error = max(
        (
            float(item["validation"]["summary"]["maxAbsWorkOrContactError"])
            for item in scenario_reports
        ),
        default=0.0,
    )
    return {
        "schema": "rad-sim.vertical-load-energy-comparison-report.v1",
        "measurementSchema": "rad-sim.vertical-load-energy-measurement-results.v1",
        "validationSchema": "rad-sim.vertical-load-energy-validation.v1",
        "method": (
            "runs each filled vertical-load measurement scenario through the "
            "same spring-hinge preview and measured energy validation helper"
        ),
        "claimLabels": {
            "zeroResidualIdentity": "Lean-proven theorem scaffold",
            "springHingePrediction": "simulator-derived empirical law",
            "benchMeasurement": "experimentally unvalidated physical assumption",
        },
        "summary": {
            "scenarioCount": len(scenario_reports),
            "missingMeasurementCount": int(missing),
            "maxAbsWorkOrContactError": float(max_error),
            "allScenariosPassTolerance": all(
                bool(item["validation"]["summary"]["passesTolerance"])
                for item in scenario_reports
            ),
        },
        "scenarios": scenario_reports,
    }


def export_vertical_load_energy_comparison_report_json(
    report: dict[str, object],
) -> str:
    return json.dumps(report, indent=2)


def vertical_load_energy_experiment_protocol(
    config: LatticeConfig,
    state: LatticeState | None = None,
    scenarios: Iterable[VerticalLoadPhysicalScenario] | None = None,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    load_case: LoadCase | None = None,
    repeat_count: int = 3,
) -> dict[str, object]:
    """Create a step-by-step protocol for vertical-load energy validation."""

    if repeat_count <= 0:
        raise ValueError("repeat_count must be positive")
    template = vertical_load_energy_measurement_template(
        config,
        state=state,
        scenarios=scenarios,
        tolerance=tolerance,
        contact_stiffness=contact_stiffness,
        load_case=load_case,
    )
    steps: list[dict[str, object]] = []
    for index, scenario in enumerate(template["scenarios"], start=1):
        if not isinstance(scenario, dict):
            continue
        input_payload = scenario.get("input", {})
        required = scenario.get("requiredMeasurements", {})
        if not isinstance(input_payload, dict) or not isinstance(required, dict):
            continue
        name = str(scenario.get("name", f"scenario-{index}"))
        steps.append(
            {
                "id": f"vertical-load-energy-{index:02d}-{name}",
                "name": name,
                "repeatCount": int(repeat_count),
                "fixture": {
                    "fixedCells": input_payload.get("fixedCells", []),
                    "removedCells": input_payload.get("removedCells", []),
                    "actuatedCells": input_payload.get("cells", []),
                    "alphaCommand": input_payload.get("alpha", 0.0),
                    "zCommand": input_payload.get("z", 0.0),
                    "externalZLoad": input_payload.get("externalZLoad", 0.0),
                },
                "setup": [
                    "Install the lattice in the same orientation as the simulator cell grid.",
                    "Clamp or otherwise constrain each fixed cell before applying the load.",
                    "Remove or disable each listed removed cell before applying actuation.",
                    "Apply the z command and external z load after the baseline height image is recorded.",
                ],
                "measurementOrder": [
                    "Record hardware profile, load fixture, and calibration notes.",
                    "Record baseline heights for every listed load cell.",
                    "Apply fixed-cell constraints and removed-cell changes.",
                    "Apply the listed z command and external z load.",
                    "Measure final height for each intact load cell.",
                    "Reset, then repeat for the removed-state load cells.",
                    "Repeat the entire scenario for the requested repeat count.",
                ],
                "requiredMeasurements": {
                    "intactLoadCells": required.get("intactLoadCells", []),
                    "removedLoadCells": required.get("removedLoadCells", []),
                    "unanchoredLoadCells": required.get("unanchoredLoadCells", {}),
                    "fields": [
                        "measuredHeight",
                        "appliedLoad",
                        "baselineHeight",
                        "repeatIndex",
                        "instrumentId",
                        "notes",
                    ],
                },
                "expectedResultFields": {
                    "scenarioName": name,
                    "measurements.intact[].measuredHeight": "required for each intact load cell",
                    "measurements.removed[].measuredHeight": "required for each removed-state load cell",
                    "input.fixedCells": "must match this protocol step",
                    "input.removedCells": "must match this protocol step",
                    "input.externalZLoad": "must match measured applied load sign convention",
                },
            }
        )
    return {
        "schema": "rad-sim.vertical-load-energy-experiment-protocol.v1",
        "measurementSchema": "rad-sim.vertical-load-energy-measurement-results.v1",
        "templateSchema": "rad-sim.vertical-load-energy-measurement-template.v1",
        "validationSchema": "rad-sim.vertical-load-energy-validation.v1",
        "comparisonReportSchema": "rad-sim.vertical-load-energy-comparison-report.v1",
        "method": (
            "physical protocol for filling vertical-load energy measurement rows "
            "and comparing them against the spring-hinge energy validation"
        ),
        "claimLabels": {
            "protocolDesign": "simulator-derived empirical law",
            "benchMeasurement": "experimentally unvalidated physical assumption",
            "zeroResidualIdentity": "Lean-proven theorem scaffold",
        },
        "grid": template["grid"],
        "contactStiffness": float(contact_stiffness),
        "tolerance": float(tolerance),
        "repeatCount": int(repeat_count),
        "requiredInstruments": [
            "height gauge or calibrated camera/motion-capture setup",
            "known z-load fixture or force sensor",
            "cell clamp or fixture for fixed cells",
            "actuator command readout for z displacement",
            "hardware profile record for pin, hole, plate, and stack dimensions",
        ],
        "safetyAndUncertaintyNotes": [
            "Keep loads below actuator, hinge, and printed-material limits.",
            "Record fixture compliance, backlash, friction, and visible slip as uncertainty sources.",
            "Treat contact and gravity terms as uncalibrated until repeated measurements agree.",
            "Do not compare signed work without preserving the simulator z-load sign convention.",
        ],
        "outputs": {
            "blankTemplate": "vertical_load_energy_measurement_template",
            "filledResultsSchema": "rad-sim.vertical-load-energy-measurement-results.v1",
            "comparisonFunction": "compare_vertical_load_energy_measurement_results",
            "comparisonReportSchema": "rad-sim.vertical-load-energy-comparison-report.v1",
        },
        "steps": steps,
    }


def export_vertical_load_energy_experiment_protocol_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    scenarios: Iterable[VerticalLoadPhysicalScenario] | None = None,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    load_case: LoadCase | None = None,
    repeat_count: int = 3,
) -> str:
    return json.dumps(
        vertical_load_energy_experiment_protocol(
            config,
            state=state,
            scenarios=scenarios,
            tolerance=tolerance,
            contact_stiffness=contact_stiffness,
            load_case=load_case,
            repeat_count=repeat_count,
        ),
        indent=2,
    )


def vertical_load_bench_packet(
    config: LatticeConfig,
    state: LatticeState | None = None,
    scenarios: Iterable[VerticalLoadPhysicalScenario] | None = None,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    load_case: LoadCase | None = None,
    repeat_count: int = 3,
    include_fields: bool = False,
    hardware_profile: object | None = None,
) -> dict[str, object]:
    """Bundle vertical-load preview, protocol, template, and comparison guidance."""

    from .cell_geometry import config_with_hardware_profile
    from .unit_scale import physical_unit_scale_metadata

    effective_config = (
        config_with_hardware_profile(config, hardware_profile)
        if hardware_profile is not None
        else config
    )
    scenario_tuple = (
        default_vertical_load_physical_scenarios(effective_config)
        if scenarios is None
        else tuple(scenarios)
    )
    base_state = (
        LatticeState.uniform(effective_config)
        if state is None
        else state.normalized(effective_config)
    )
    preview = build_vertical_load_physical_preview_report(
        effective_config,
        state=base_state,
        scenarios=scenario_tuple,
        tolerance=tolerance,
        contact_stiffness=contact_stiffness,
        load_case=load_case,
        include_fields=include_fields,
    )
    template = vertical_load_energy_measurement_template(
        effective_config,
        state=base_state,
        scenarios=scenario_tuple,
        tolerance=tolerance,
        contact_stiffness=contact_stiffness,
        load_case=load_case,
    )
    protocol = vertical_load_energy_experiment_protocol(
        effective_config,
        state=base_state,
        scenarios=scenario_tuple,
        tolerance=tolerance,
        contact_stiffness=contact_stiffness,
        load_case=load_case,
        repeat_count=repeat_count,
    )
    certificate_summaries: list[dict[str, object]] = []
    for scenario in preview.get("scenarios", []):
        if not isinstance(scenario, dict):
            continue
        comparison = scenario.get("comparison", {})
        if not isinstance(comparison, dict):
            continue
        certificate = comparison.get("mechanicsCertificate", {})
        if not isinstance(certificate, dict):
            continue
        certificate_summaries.append(
            {
                "scenario": scenario.get("name", ""),
                "schema": certificate.get("schema", ""),
                "passesNonnegativeCheck": certificate.get("passesNonnegativeCheck", False),
                "claimLabels": certificate.get("claimLabels", {}),
                "delta": certificate.get("delta", {}),
            }
        )
    return {
        "schema": "rad-sim.vertical-load-bench-packet.v1",
        "method": (
            "publication and bench packet bundling vertical-load physical preview, "
            "mechanics certificate summaries, fillable measurement template, "
            "experiment protocol, and comparison instructions"
        ),
        "claimLabels": {
            "zeroResidualIdentity": "Lean-proven theorem scaffold",
            "springHingePrediction": "simulator-derived empirical law",
            "benchMeasurement": "experimentally unvalidated physical assumption",
            "packetAssembly": "simulator-derived empirical law",
        },
        "schemas": {
            "preview": "rad-sim.vertical-load-physical-preview-report.v1",
            "comparison": "rad-sim.vertical-removal-physical-comparison.v1",
            "mechanicsCertificate": "rad-sim.mechanics-energy-certificate.v1",
            "protocol": "rad-sim.vertical-load-energy-experiment-protocol.v1",
            "measurementTemplate": "rad-sim.vertical-load-energy-measurement-results.v1",
            "energyValidation": "rad-sim.vertical-load-energy-validation.v1",
            "comparisonReport": "rad-sim.vertical-load-energy-comparison-report.v1",
            "csv": "vertical-load-physical-preview CSV",
            "hardwareProfile": "rad-sim.hardware-profile.v1",
        },
        "grid": {
            "rows": effective_config.rows,
            "cols": effective_config.cols,
            "pinHoleClearance": effective_config.pin_hole_clearance,
            "zCouplingGain": effective_config.z_coupling_gain,
        },
        "parameters": {
            "tolerance": float(tolerance),
            "contactStiffness": float(contact_stiffness),
            "repeatCount": int(repeat_count),
            "includeFields": bool(include_fields),
        },
        "unitScaleMetadata": physical_unit_scale_metadata(
            config,
            hardware_profile=hardware_profile,
        ),
        "previewReport": preview,
        "measurementTemplate": template,
        "experimentProtocol": protocol,
        "mechanicsCertificateSummaries": certificate_summaries,
        "comparisonInstructions": {
            "fillTemplate": "Fill measuredHeight in measurementTemplate.scenarios[].measurements.",
            "parseFunction": "vertical_load_energy_measurement_results_from_json",
            "compareFunction": "compare_vertical_load_energy_measurement_results",
            "exportFunction": "export_vertical_load_energy_comparison_report_json",
            "acceptanceRule": (
                "A scenario passes the simulator comparison only when all required "
                "load-cell heights are measured and work/contact errors are within tolerance."
            ),
        },
        "limitations": [
            "This packet is simulator-derived bench planning, not calibrated hardware validation.",
            "Contact, gravity, friction, and fixture compliance remain physical assumptions until measured.",
            "Predicted heights are included for planning and software roundtrips, not as bench measurements.",
        ],
    }


def export_vertical_load_bench_packet_json(
    config: LatticeConfig,
    state: LatticeState | None = None,
    scenarios: Iterable[VerticalLoadPhysicalScenario] | None = None,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    load_case: LoadCase | None = None,
    repeat_count: int = 3,
    include_fields: bool = False,
    hardware_profile: object | None = None,
) -> str:
    return json.dumps(
        vertical_load_bench_packet(
            config,
            state=state,
            scenarios=scenarios,
            tolerance=tolerance,
            contact_stiffness=contact_stiffness,
            load_case=load_case,
            repeat_count=repeat_count,
            include_fields=include_fields,
            hardware_profile=hardware_profile,
        ),
        indent=2,
    )


def _topology_dict(topology: dict[str, Any]) -> dict[str, object]:
    return {
        "componentLabels": np.asarray(topology["component_labels"], dtype=int).tolist(),
        "componentCount": int(topology["component_count"]),
        "componentSizes": list(topology["component_sizes"]),
        "largestComponentSize": int(topology["largest_component_size"]),
        "presentCellCount": int(topology["present_cell_count"]),
        "removedCellCount": int(topology["removed_cell_count"]),
        "totalEdgeCount": int(topology["total_edge_count"]),
        "activeEdgeCount": int(topology["active_edge_count"]),
        "deletedEdgeCount": int(topology["deleted_edge_count"]),
    }


def _vertical_residual_diagnostic_dict(
    diagnostic: VerticalRemovalResidualDiagnostic,
    include_fields: bool,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "cells": _cells_dict(diagnostic.cells),
        "removedCells": _cells_dict(diagnostic.removed_cells),
        "fixedCells": _cells_dict(diagnostic.fixed_cells),
        "removedSourceCells": _cells_dict(diagnostic.removed_source_cells),
        "affectedNeighborCells": _cells_dict(diagnostic.affected_neighbor_cells),
        "topologyBlockedCells": _cells_dict(diagnostic.topology_blocked_cells),
        "intactLoadActiveCells": _cells_dict(diagnostic.intact_load_active_cells),
        "removedLoadActiveCells": _cells_dict(diagnostic.removed_load_active_cells),
        "intactHeightContactEngagedCells": _cells_dict(
            diagnostic.intact_height_contact_engaged_cells
        ),
        "removedHeightContactEngagedCells": _cells_dict(
            diagnostic.removed_height_contact_engaged_cells
        ),
        "summary": {
            "intactReachCount": diagnostic.intact_reach_count,
            "removedReachCount": diagnostic.removed_reach_count,
            "intactNeighborReachCount": diagnostic.intact_neighbor_reach_count,
            "removedNeighborReachCount": diagnostic.removed_neighbor_reach_count,
            "intactDieOffRadius": diagnostic.intact_die_off_radius,
            "removedDieOffRadius": diagnostic.removed_die_off_radius,
            "maxAbsNeighborResidual": diagnostic.max_abs_neighbor_residual,
            "maxAbsResidualDelta": diagnostic.max_abs_residual_delta,
            "clearance": diagnostic.clearance,
            "contactStiffness": diagnostic.contact_stiffness,
            "contactPenaltyDelta": diagnostic.contact_penalty_delta,
            "externalZLoad": diagnostic.external_z_load,
            "loadWorkDelta": diagnostic.load_work_delta,
            "loadWorkMagnitudeDelta": diagnostic.load_work_magnitude_delta,
            "heightContactPenaltyDelta": diagnostic.height_contact_penalty_delta,
            "componentCountDelta": diagnostic.component_count_delta,
            "deletedEdgeDelta": diagnostic.deleted_edge_delta,
        },
        "topology": {
            "intact": _topology_dict(diagnostic.intact_topology),
            "removed": _topology_dict(diagnostic.removed_topology),
        },
    }
    if include_fields:
        payload["fields"] = {
            "intactZResidual": _matrix_json(diagnostic.intact_z_residual),
            "removedZResidual": _matrix_json(diagnostic.removed_z_residual),
            "residualDelta": _matrix_json(diagnostic.residual_delta),
            "intactZDieOff": _matrix_json(diagnostic.intact_z_die_off),
            "removedZDieOff": _matrix_json(diagnostic.removed_z_die_off),
        }
    return payload


def vertical_removal_physical_comparison_to_dict(
    comparison: VerticalRemovalPhysicalComparisonDiagnostic,
    include_fields: bool = False,
) -> dict[str, object]:
    """Serialize a vertical removal/spring-hinge comparison as a research artifact."""

    payload: dict[str, object] = {
        "schema": "rad-sim.vertical-removal-physical-comparison.v1",
        "method": (
            "matched intact-versus-removed vertical residual diagnostic compared "
            "against solve_spring_hinge_3d on the same removed-cell spring graph"
        ),
        "claimLabels": {
            "clearanceDeadZone": "Lean-proven theorem",
            "verticalResidual": "simulator-derived empirical law",
            "springHingePhysicalPreview": "simulator-derived empirical law",
            "gravityContact": "experimentally unvalidated physical assumption",
        },
        "verticalDiagnostic": _vertical_residual_diagnostic_dict(
            comparison.vertical_diagnostic, include_fields=include_fields
        ),
        "loadCases": {
            "intact": _load_case_dict(comparison.intact_load_case),
            "removed": _load_case_dict(comparison.removed_load_case),
        },
        "solver": {
            "model": "spring_hinge_3d",
            "physicalSuccess": comparison.physical_success,
            "intactSuccess": comparison.intact_solver_success,
            "removedSuccess": comparison.removed_solver_success,
            "intactEnergy": comparison.intact_solver_energy,
            "removedEnergy": comparison.removed_solver_energy,
            "intactEnergyBreakdown": _energy_breakdown_dict(
                comparison.intact_energy_breakdown
            ),
            "removedEnergyBreakdown": _energy_breakdown_dict(
                comparison.removed_energy_breakdown
            ),
            "intactIterations": comparison.intact_solver_iterations,
            "removedIterations": comparison.removed_solver_iterations,
            "intactSpringEdges": comparison.intact_spring_edges,
            "removedSpringEdges": comparison.removed_spring_edges,
            "springEdgeDelta": comparison.spring_edge_delta,
            "intactHingeTriples": comparison.intact_hinge_triples,
            "removedHingeTriples": comparison.removed_hinge_triples,
            "hingeTripleDelta": comparison.hinge_triple_delta,
            "intactRemovedCells": comparison.intact_solver_removed_cells,
            "removedRemovedCells": comparison.removed_solver_removed_cells,
            "intactSolverLoadCells": _cells_dict(comparison.intact_solver_load_cells),
            "removedSolverLoadCells": _cells_dict(comparison.removed_solver_load_cells),
            "intactUnanchoredLoadCells": _cells_dict(
                comparison.intact_unanchored_load_cells
            ),
            "removedUnanchoredLoadCells": _cells_dict(
                comparison.removed_unanchored_load_cells
            ),
            "topologyNote": comparison.solver_topology_note,
        },
        "metrics": {
            "intactHeightRmsModelError": comparison.intact_height_rms_model_error,
            "removedHeightRmsModelError": comparison.removed_height_rms_model_error,
            "deltaRmsModelError": comparison.delta_rms_model_error,
            "maxAbsDeltaModelError": comparison.max_abs_delta_model_error,
        },
        "mechanicsCertificate": mechanics_energy_certificate_to_dict(comparison),
    }
    if include_fields:
        payload["fields"] = {
            "intactPhysicalHeight": _matrix_json(comparison.intact_physical_height),
            "removedPhysicalHeight": _matrix_json(comparison.removed_physical_height),
            "intactKinematicHeight": _matrix_json(comparison.intact_kinematic_height),
            "removedKinematicHeight": _matrix_json(comparison.removed_kinematic_height),
            "intactModelError": _matrix_json(comparison.intact_model_error),
            "removedModelError": _matrix_json(comparison.removed_model_error),
            "physicalHeightDelta": _matrix_json(comparison.physical_height_delta),
            "kinematicHeightDelta": _matrix_json(comparison.kinematic_height_delta),
            "physicalVsKinematicDeltaError": _matrix_json(
                comparison.physical_vs_kinematic_delta_error
            ),
        }
    return payload


def export_vertical_removal_physical_comparison_json(
    comparison: VerticalRemovalPhysicalComparisonDiagnostic,
    include_fields: bool = False,
) -> str:
    return json.dumps(
        vertical_removal_physical_comparison_to_dict(
            comparison, include_fields=include_fields
        ),
        indent=2,
    )


def default_vertical_load_physical_scenarios(
    config: LatticeConfig,
    z: float = 0.25,
    external_z_load: float = -0.05,
) -> tuple[VerticalLoadPhysicalScenario, ...]:
    if config.cols < 3:
        raise ValueError("default vertical load physical scenarios require at least 3 columns")
    row = config.rows // 2
    fixed = ((row, 0),)
    return (
        VerticalLoadPhysicalScenario(
            name="one-cell-z-load",
            cells=((row, 1),),
            fixed_cells=fixed,
            z=z,
            external_z_load=external_z_load,
        ),
        VerticalLoadPhysicalScenario(
            name="two-cell-z-load",
            cells=((row, 1), (row, 2)),
            fixed_cells=fixed,
            z=z,
            external_z_load=external_z_load,
        ),
        VerticalLoadPhysicalScenario(
            name="removed-middle-z-load",
            cells=((row, 2),),
            removed_cells=((row, 1),),
            fixed_cells=fixed,
            z=z,
            external_z_load=external_z_load,
        ),
    )


def build_vertical_load_physical_preview_report(
    config: LatticeConfig,
    state: LatticeState | None = None,
    scenarios: Iterable[VerticalLoadPhysicalScenario] | None = None,
    tolerance: float = 1e-9,
    contact_stiffness: float = 1.0,
    load_case: LoadCase | None = None,
    include_fields: bool = False,
) -> dict[str, object]:
    base_state = LatticeState.uniform(config) if state is None else state.normalized(config)
    scenario_tuple = (
        default_vertical_load_physical_scenarios(config)
        if scenarios is None
        else tuple(scenarios)
    )
    comparisons: list[dict[str, object]] = []
    raw_comparisons: list[VerticalRemovalPhysicalComparisonDiagnostic] = []
    for scenario in scenario_tuple:
        comparison = compare_vertical_residual_spring_hinge_3d(
            config,
            base_state,
            scenario.cells,
            scenario.removed_cells,
            scenario.z,
            alpha=scenario.alpha,
            tolerance=tolerance,
            contact_stiffness=contact_stiffness,
            fixed_cells=scenario.fixed_cells,
            external_z_load=scenario.external_z_load,
            load_case=load_case,
        )
        raw_comparisons.append(comparison)
        comparisons.append(
            {
                "name": scenario.name,
                "input": {
                    "cells": _cells_dict(scenario.cells),
                    "removedCells": _cells_dict(scenario.removed_cells),
                    "fixedCells": _cells_dict(scenario.fixed_cells),
                    "alpha": scenario.alpha,
                    "z": scenario.z,
                    "externalZLoad": scenario.external_z_load,
                },
                "comparison": vertical_removal_physical_comparison_to_dict(
                    comparison, include_fields=include_fields
                ),
            }
        )
    return {
        "schema": "rad-sim.vertical-load-physical-preview-report.v1",
        "method": (
            "named z-load scenarios comparing vertical residual diagnostics "
            "with 3D spring-hinge physical previews"
        ),
        "claimLabels": {
            "springHingePhysicalPreview": "simulator-derived empirical law",
            "gravityContact": "experimentally unvalidated physical assumption",
        },
        "grid": {
            "rows": config.rows,
            "cols": config.cols,
            "pinHoleClearance": config.pin_hole_clearance,
            "zCouplingGain": config.z_coupling_gain,
        },
        "summary": {
            "scenarioCount": len(raw_comparisons),
            "allPhysicalSuccess": all(
                comparison.physical_success for comparison in raw_comparisons
            ),
            "maxDeltaRmsModelError": max(
                (comparison.delta_rms_model_error for comparison in raw_comparisons),
                default=0.0,
            ),
            "maxAbsDeltaModelError": max(
                (comparison.max_abs_delta_model_error for comparison in raw_comparisons),
                default=0.0,
            ),
            "totalUnanchoredLoadCells": sum(
                len(comparison.intact_unanchored_load_cells)
                + len(comparison.removed_unanchored_load_cells)
                for comparison in raw_comparisons
            ),
        },
        "scenarios": comparisons,
    }


def export_vertical_load_physical_preview_report_json(
    report: dict[str, object],
) -> str:
    return json.dumps(report, indent=2)


def _csv_scalar(value: object) -> str:
    if value is None:
        text = ""
    elif isinstance(value, bool):
        text = "true" if value else "false"
    else:
        text = str(value)
    if any(char in text for char in [",", '"', "\n", "\r"]):
        return '"' + text.replace('"', '""') + '"'
    return text


def _csv_cells(cells: object) -> str:
    if not isinstance(cells, list):
        return ""
    out: list[str] = []
    for cell in cells:
        if not isinstance(cell, dict):
            continue
        out.append(f"{int(cell.get('row', 0))}:{int(cell.get('col', 0))}")
    return ";".join(out)


def _claim_label_summary(claims: object) -> str:
    if not isinstance(claims, dict):
        return ""
    return ";".join(f"{key}={value}" for key, value in claims.items())


def export_vertical_load_physical_preview_report_csv(report: dict[str, object]) -> str:
    """Export a vertical-load physical-preview report as a bench-table CSV."""

    header = [
        "schema",
        "scenario",
        "command_cells",
        "fixed_cells",
        "removed_cells",
        "alpha_command",
        "z_command",
        "external_z_load",
        "physical_success",
        "intact_success",
        "removed_success",
        "intact_unanchored_load_cells",
        "removed_unanchored_load_cells",
        "spring_edge_delta",
        "hinge_triple_delta",
        "delta_rms_model_error",
        "max_abs_delta_model_error",
        "contact_penalty_delta",
        "load_work_delta",
        "load_work_magnitude_delta",
        "height_contact_penalty_delta",
        "claim_labels",
    ]
    rows = [header]
    report_claims = report.get("claimLabels", {})
    for scenario in report.get("scenarios", []):
        if not isinstance(scenario, dict):
            continue
        input_payload = scenario.get("input", {})
        comparison = scenario.get("comparison", {})
        if not isinstance(input_payload, dict) or not isinstance(comparison, dict):
            continue
        solver = comparison.get("solver", {})
        metrics = comparison.get("metrics", {})
        vertical = comparison.get("verticalDiagnostic", {})
        summary = vertical.get("summary", {}) if isinstance(vertical, dict) else {}
        row = [
            report.get("schema", ""),
            scenario.get("name", ""),
            _csv_cells(input_payload.get("cells")),
            _csv_cells(input_payload.get("fixedCells")),
            _csv_cells(input_payload.get("removedCells")),
            input_payload.get("alpha", ""),
            input_payload.get("z", ""),
            input_payload.get("externalZLoad", ""),
            solver.get("physicalSuccess", "") if isinstance(solver, dict) else "",
            solver.get("intactSuccess", "") if isinstance(solver, dict) else "",
            solver.get("removedSuccess", "") if isinstance(solver, dict) else "",
            _csv_cells(
                solver.get("intactUnanchoredLoadCells") if isinstance(solver, dict) else []
            ),
            _csv_cells(
                solver.get("removedUnanchoredLoadCells") if isinstance(solver, dict) else []
            ),
            solver.get("springEdgeDelta", "") if isinstance(solver, dict) else "",
            solver.get("hingeTripleDelta", "") if isinstance(solver, dict) else "",
            metrics.get("deltaRmsModelError", "") if isinstance(metrics, dict) else "",
            metrics.get("maxAbsDeltaModelError", "") if isinstance(metrics, dict) else "",
            summary.get("contactPenaltyDelta", "") if isinstance(summary, dict) else "",
            summary.get("loadWorkDelta", "") if isinstance(summary, dict) else "",
            summary.get("loadWorkMagnitudeDelta", "") if isinstance(summary, dict) else "",
            summary.get("heightContactPenaltyDelta", "") if isinstance(summary, dict) else "",
            _claim_label_summary(report_claims),
        ]
        rows.append(row)
    return "\n".join(",".join(_csv_scalar(value) for value in row) for row in rows)


def evaluate_programmable_operators(
    config: LatticeConfig, state: LatticeState
) -> dict[str, np.ndarray]:
    state = state.normalized(config)
    active_mask = ~state.removed_mask
    topology = lattice_topology_diagnostic(config, state)
    alpha = alpha_backlash_operator(config).propagate(
        config.rows,
        config.cols,
        nonzero_sources(state.actuator_grid),
        active_mask=active_mask,
    )
    vertical = vertical_clearance_operator(config).propagate(
        config.rows,
        config.cols,
        nonzero_sources(state.z_actuator_grid),
        active_mask=active_mask,
    )
    raw_alpha = state.alpha_grid + alpha.field
    locked_mask = state.locked_mask & active_mask
    projected_alpha = lock_projection(state.alpha_grid, raw_alpha, locked_mask)
    projected_alpha = np.where(active_mask, projected_alpha, config.initial_alpha)
    raw_height = -0.65 * config.cell_size * alpha.field + vertical.field
    height = lock_projection(state.lock_z_grid, raw_height, locked_mask)
    height = np.where(active_mask, height, 0.0)
    return {
        "actuator_influence": alpha.field,
        "die_off": alpha.die_off,
        "z_residual": vertical.field,
        "z_die_off": vertical.die_off,
        "raw_alpha": raw_alpha,
        "alpha": projected_alpha,
        "height": height,
        "removed_mask": state.removed_mask,
        "topology": topology,
    }
