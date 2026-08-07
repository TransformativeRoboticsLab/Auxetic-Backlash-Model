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
