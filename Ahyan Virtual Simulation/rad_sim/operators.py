from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .coupling import backlash_activation
from .geometry import neighbor_indices
from .models import LatticeConfig, LatticeState


@dataclass(frozen=True)
class PropagationResult:
    field: np.ndarray
    die_off: np.ndarray


@dataclass(frozen=True)
class ProgrammableDiscontinuityEvent:
    kind: str
    cell: tuple[int, int] | None = None
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
    ) -> PropagationResult:
        field = np.zeros((rows, cols), dtype=float)
        die_off = np.full((rows, cols), np.inf, dtype=float)
        steps = self.max_steps or (rows + cols)

        for r0, c0, value in sources:
            if abs(value) < 1e-12:
                continue
            queue: deque[tuple[int, int, float, int]] = deque([(r0, c0, value, 0)])
            visited: dict[tuple[int, int], float] = {}
            while queue:
                r, c, signal, dist = queue.popleft()
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
                    queue.append((nr, nc, next_signal, dist + 1))

        return PropagationResult(field=field, die_off=die_off)


def finite_die_off_radius(die_off: np.ndarray) -> int:
    finite = np.asarray(die_off)[np.isfinite(die_off)]
    if finite.size == 0:
        return 0
    return int(np.max(finite))


def nonzero_sources(values: np.ndarray) -> list[tuple[int, int, float]]:
    return [
        (int(r), int(c), float(values[r, c]))
        for r, c in zip(*np.nonzero(np.abs(values) > 1e-12), strict=False)
    ]


def alpha_backlash_operator(config: LatticeConfig) -> DeadZonePropagationOperator:
    return DeadZonePropagationOperator(
        name="alpha_backlash",
        dead_zone=config.backlash,
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


def lock_event(cell: tuple[int, int]) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="lock", cell=cell)


def release_event(cell: tuple[int, int]) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="release", cell=cell)


def clear_actuation_event(
    cell: tuple[int, int] | None = None,
) -> ProgrammableDiscontinuityEvent:
    return ProgrammableDiscontinuityEvent(kind="clear_actuation", cell=cell)


def _validate_cell(config: LatticeConfig, cell: tuple[int, int] | None) -> tuple[int, int]:
    if cell is None:
        raise ValueError("event requires a cell")
    r, c = int(cell[0]), int(cell[1])
    if not (0 <= r < config.rows and 0 <= c < config.cols):
        raise ValueError(f"cell {(r, c)} is outside lattice shape {(config.rows, config.cols)}")
    return r, c


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
    )

    if event.kind == "actuate":
        r, c = _validate_cell(config, event.cell)
        next_state.actuator_grid[r, c] += float(event.alpha)
        next_state.z_actuator_grid[r, c] += float(event.z)
    elif event.kind == "lock":
        r, c = _validate_cell(config, event.cell)
        fields = evaluate_programmable_operators(config, current)
        next_state.alpha_grid[r, c] = fields["alpha"][r, c]
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
        mode_changes=_mode_change_count(first_state.locked_mask, second_state.locked_mask),
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


def evaluate_programmable_operators(
    config: LatticeConfig, state: LatticeState
) -> dict[str, np.ndarray]:
    state = state.normalized(config)
    alpha = alpha_backlash_operator(config).propagate(
        config.rows, config.cols, nonzero_sources(state.actuator_grid)
    )
    vertical = vertical_clearance_operator(config).propagate(
        config.rows, config.cols, nonzero_sources(state.z_actuator_grid)
    )
    raw_alpha = state.alpha_grid + alpha.field
    projected_alpha = lock_projection(state.alpha_grid, raw_alpha, state.locked_mask)
    height = lock_projection(
        np.zeros_like(projected_alpha),
        -0.65 * config.cell_size * alpha.field + vertical.field,
        state.locked_mask,
    )
    return {
        "actuator_influence": alpha.field,
        "die_off": alpha.die_off,
        "z_residual": vertical.field,
        "z_die_off": vertical.die_off,
        "raw_alpha": raw_alpha,
        "alpha": projected_alpha,
        "height": height,
    }
