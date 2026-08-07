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
