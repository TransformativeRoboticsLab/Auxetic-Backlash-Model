from __future__ import annotations

from collections import deque

import numpy as np

from .coupling import alpha_to_theta, backlash_activation
from .geometry import cell_corners, grid_centers, neighbor_indices
from .models import LatticeConfig, LatticeState, SimulationResult


def _actuator_influence(config: LatticeConfig, actuator_grid: np.ndarray) -> np.ndarray:
    influence = np.zeros((config.rows, config.cols), dtype=float)
    steps = config.max_coupling_steps or (config.rows + config.cols)

    for source in zip(*np.nonzero(np.abs(actuator_grid) > 1e-12), strict=False):
        r0, c0 = int(source[0]), int(source[1])
        queue: deque[tuple[int, int, float, int]] = deque(
            [(r0, c0, float(actuator_grid[r0, c0]), 0)]
        )
        visited: dict[tuple[int, int], float] = {}
        while queue:
            r, c, signal, dist = queue.popleft()
            if dist > steps:
                continue
            key = (r, c)
            if abs(signal) <= abs(visited.get(key, 0.0)):
                continue
            visited[key] = signal
            influence[r, c] += signal
            next_signal = float(backlash_activation(signal, config.backlash)) * config.coupling_gain
            if abs(next_signal) < 1e-9:
                continue
            for nr, nc in neighbor_indices(config.rows, config.cols, r, c):
                queue.append((nr, nc, next_signal, dist + 1))

    return influence


def _vertical_residual(config: LatticeConfig, z_actuator_grid: np.ndarray) -> np.ndarray:
    residual = np.zeros((config.rows, config.cols), dtype=float)
    steps = config.max_coupling_steps or (config.rows + config.cols)
    dead_zone = min(0.18 * config.cell_size, config.backlash * 0.45)

    for source in zip(*np.nonzero(np.abs(z_actuator_grid) > 1e-12), strict=False):
        r0, c0 = int(source[0]), int(source[1])
        queue: deque[tuple[int, int, float, int]] = deque(
            [(r0, c0, float(z_actuator_grid[r0, c0]), 0)]
        )
        visited: dict[tuple[int, int], float] = {}
        while queue:
            r, c, signal, dist = queue.popleft()
            if dist > steps:
                continue
            key = (r, c)
            if abs(signal) <= abs(visited.get(key, 0.0)):
                continue
            visited[key] = signal
            residual[r, c] += signal
            next_signal = (
                float(backlash_activation(signal, dead_zone)) * config.z_coupling_gain
            )
            if abs(next_signal) < 1e-9:
                continue
            for nr, nc in neighbor_indices(config.rows, config.cols, r, c):
                queue.append((nr, nc, next_signal, dist + 1))

    return residual


def simulate_kinematic(config: LatticeConfig, state: LatticeState) -> SimulationResult:
    state = state.normalized(config)
    influence = _actuator_influence(config, state.actuator_grid)
    z_residual = _vertical_residual(config, state.z_actuator_grid)
    raw_alpha = state.alpha_grid + influence
    alpha = np.where(state.locked_mask, state.alpha_grid, raw_alpha)
    alpha = np.clip(alpha, config.alpha_min, config.alpha_max)
    theta = alpha_to_theta(alpha)

    original_centers = grid_centers(config)
    deformed_centers = grid_centers(config, alpha)
    original_corners = cell_corners(
        config, original_centers, np.full_like(alpha, config.initial_alpha)
    )
    deformed_corners = cell_corners(config, deformed_centers, alpha)
    height = np.where(state.locked_mask, 0.0, -0.65 * config.cell_size * influence + z_residual)
    original_centers_3d = np.dstack(
        [original_centers[..., 0], original_centers[..., 1], np.zeros_like(alpha)]
    )
    deformed_centers_3d = np.dstack(
        [deformed_centers[..., 0], deformed_centers[..., 1], height]
    )
    original_corners_3d = np.zeros((config.rows, config.cols, 5, 3), dtype=float)
    deformed_corners_3d = np.zeros((config.rows, config.cols, 5, 3), dtype=float)
    original_corners_3d[..., :2] = original_corners
    deformed_corners_3d[..., :2] = deformed_corners
    deformed_corners_3d[..., 2] = height[..., None]

    z0 = original_centers[..., 0] + 1j * original_centers[..., 1]
    z1 = deformed_centers[..., 0] + 1j * deformed_centers[..., 1]
    return SimulationResult(
        config=config,
        state=state,
        alpha=alpha,
        theta_degrees=theta,
        original_centers=original_centers,
        deformed_centers=deformed_centers,
        original_corners=original_corners,
        deformed_corners=deformed_corners,
        original_centers_3d=original_centers_3d,
        deformed_centers_3d=deformed_centers_3d,
        original_corners_3d=original_corners_3d,
        deformed_corners_3d=deformed_corners_3d,
        complex_original=z0,
        complex_deformed=z1,
        metadata={
            "model": "kinematic",
            "mean_alpha": float(np.mean(alpha)),
            "actuator_influence": influence,
            "z_residual": z_residual,
            "height": height,
        },
    )
