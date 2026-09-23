from __future__ import annotations

import numpy as np

from .coupling import alpha_to_theta
from .models import LatticeConfig


def alpha_pitch_scale(config: LatticeConfig, alpha: np.ndarray | float) -> np.ndarray | float:
    return np.maximum(np.asarray(alpha, dtype=float), 1e-9) / max(
        float(config.initial_alpha), 1e-9
    )


def shared_edge_pitch(
    config: LatticeConfig, first_alpha: np.ndarray | float, second_alpha: np.ndarray | float
) -> np.ndarray | float:
    return (
        0.5
        * float(config.cell_size)
        * (alpha_pitch_scale(config, first_alpha) + alpha_pitch_scale(config, second_alpha))
    )


def grid_centers(config: LatticeConfig, alpha: np.ndarray | None = None) -> np.ndarray:
    alpha_grid = (
        np.full((config.rows, config.cols), config.initial_alpha, dtype=float)
        if alpha is None
        else np.asarray(alpha, dtype=float)
    )
    if alpha_grid.shape != (config.rows, config.cols):
        raise ValueError(f"alpha must have shape {(config.rows, config.cols)}")

    x = np.zeros((config.rows, config.cols), dtype=float)
    y = np.zeros((config.rows, config.cols), dtype=float)
    if config.cols > 1:
        x[:, 1:] = np.cumsum(
            shared_edge_pitch(config, alpha_grid[:, :-1], alpha_grid[:, 1:]), axis=1
        )
    if config.rows > 1:
        y[1:, :] = np.cumsum(
            shared_edge_pitch(config, alpha_grid[:-1, :], alpha_grid[1:, :]), axis=0
        )

    centers = np.zeros((config.rows, config.cols, 2), dtype=float)
    centers[..., 0] = x
    centers[..., 1] = y
    centers[..., 0] -= centers[..., 0].mean()
    centers[..., 1] -= centers[..., 1].mean()
    return centers


def cell_corners(
    config: LatticeConfig, centers: np.ndarray, alpha: np.ndarray
) -> np.ndarray:
    theta = np.deg2rad(alpha_to_theta(alpha))
    side = 0.62 * config.cell_size * np.sqrt(np.clip(alpha, config.alpha_min, config.alpha_max))
    base = np.array(
        [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5], [-0.5, -0.5]],
        dtype=float,
    )
    corners = np.zeros((config.rows, config.cols, 5, 2), dtype=float)
    for r in range(config.rows):
        for c in range(config.cols):
            rot = np.array(
                [
                    [np.cos(theta[r, c]), -np.sin(theta[r, c])],
                    [np.sin(theta[r, c]), np.cos(theta[r, c])],
                ]
            )
            corners[r, c] = centers[r, c] + (base * side[r, c]) @ rot.T
    return corners


def neighbor_indices(rows: int, cols: int, r: int, c: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    if r > 0:
        out.append((r - 1, c))
    if r + 1 < rows:
        out.append((r + 1, c))
    if c > 0:
        out.append((r, c - 1))
    if c + 1 < cols:
        out.append((r, c + 1))
    return out
