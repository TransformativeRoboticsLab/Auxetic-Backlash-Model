from __future__ import annotations

import numpy as np

from .coupling import alpha_to_theta
from .models import LatticeConfig


def grid_centers(config: LatticeConfig, alpha: np.ndarray | None = None) -> np.ndarray:
    if alpha is None:
        row_spacing = np.full(config.rows, config.cell_size)
        col_spacing = np.full(config.cols, config.cell_size)
    else:
        row_spacing = np.mean(alpha, axis=1) * config.cell_size
        col_spacing = np.mean(alpha, axis=0) * config.cell_size

    y = np.zeros(config.rows)
    x = np.zeros(config.cols)
    if config.rows > 1:
        y[1:] = np.cumsum((row_spacing[:-1] + row_spacing[1:]) / 2.0)
    if config.cols > 1:
        x[1:] = np.cumsum((col_spacing[:-1] + col_spacing[1:]) / 2.0)

    centers = np.zeros((config.rows, config.cols, 2), dtype=float)
    centers[..., 0] = x[None, :]
    centers[..., 1] = y[:, None]
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

