from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

from .coupling import (
    alpha_to_theta,
    backlash_activation,
    normalized_backlash_to_theta_dead_zone,
)
from .geometry import cell_corners, grid_centers, neighbor_indices
from .models import LatticeConfig, LatticeState, SimulationResult
from .operators import evaluate_programmable_operators


@dataclass(frozen=True)
class ConstraintKinematicSettings:
    angle_alpha_weight: float = 20.0
    edge_distance_weight: float = 8.0
    backlash_weight: float = 2.0
    lock_weight: float = 80.0
    position_lock_weight: float = 100.0
    wall_weight: float = 70.0
    alpha_target_weight: float = 0.8
    actuator_alpha_weight: float = 35.0
    axis_alignment_weight: float = 6.0
    z_regularization_weight: float = 0.03
    xy_regularization_weight: float = 0.04
    continuation_weight: float = 0.02
    max_nfev: int = 800


def active_neighbor_edges(config: LatticeConfig, state: LatticeState) -> list[tuple[int, int, int, int]]:
    edges: list[tuple[int, int, int, int]] = []
    active = ~state.removed_mask
    for r in range(config.rows):
        for c in range(config.cols):
            if not active[r, c]:
                continue
            for nr, nc in neighbor_indices(config.rows, config.cols, r, c):
                if not active[nr, nc] or (nr, nc) < (r, c):
                    continue
                edges.append((r, c, nr, nc))
    return edges


def pair_distance_from_theta(config: LatticeConfig, theta_i_deg: float, theta_j_deg: float) -> float:
    theta_i = np.deg2rad(theta_i_deg)
    theta_j = np.deg2rad(theta_j_deg)
    return float(config.cell_size * (np.cos(theta_i) + np.cos(theta_j)))


def _boundary_value(boundary: dict[str, Any], key: str) -> float | None:
    value = boundary.get(key)
    if value is None or value == "":
        return None
    numeric = float(value)
    return numeric if np.isfinite(numeric) else None


def _initial_fields(config: LatticeConfig, state: LatticeState) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    operators = evaluate_programmable_operators(config, state)
    alpha = np.clip(operators["alpha"], config.alpha_min, config.alpha_max)
    theta = alpha_to_theta(alpha)
    centers_2d = grid_centers(config, alpha)
    height = operators["height"].copy()
    centers = np.dstack([centers_2d[..., 0], centers_2d[..., 1], height])
    return alpha, theta, centers


def _pack(alpha: np.ndarray, theta: np.ndarray, centers: np.ndarray) -> np.ndarray:
    return np.concatenate([alpha.ravel(), theta.ravel(), centers.reshape(-1)])


def _unpack(config: LatticeConfig, vector: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = config.rows * config.cols
    alpha = vector[:n].reshape((config.rows, config.cols))
    theta = vector[n : 2 * n].reshape((config.rows, config.cols))
    centers = vector[2 * n :].reshape((config.rows, config.cols, 3))
    return alpha, theta, centers


def _variable_indices(config: LatticeConfig, r: int, c: int) -> tuple[int, int, int, int, int]:
    idx = r * config.cols + c
    n = config.rows * config.cols
    center0 = 2 * n + 3 * idx
    return idx, n + idx, center0, center0 + 1, center0 + 2


def _wall_residuals(centers: np.ndarray, boundary: dict[str, Any]) -> list[float]:
    if not boundary or boundary.get("wallType") == "inactive":
        return []
    mode = boundary.get("mode", "free")
    residuals: list[float] = []
    if mode == "walls":
        x_min = _boundary_value(boundary, "xMin")
        x_max = _boundary_value(boundary, "xMax")
        y_min = _boundary_value(boundary, "yMin")
        y_max = _boundary_value(boundary, "yMax")
        for x, y in centers[..., :2].reshape((-1, 2)):
            if x_min is not None:
                residuals.append(max(0.0, x_min - float(x)))
            if x_max is not None:
                residuals.append(max(0.0, float(x) - x_max))
            if y_min is not None:
                residuals.append(max(0.0, y_min - float(y)))
            if y_max is not None:
                residuals.append(max(0.0, float(y) - y_max))
    elif mode == "channel":
        width = max(0.0, float(boundary.get("channelWidth", 0.0) or 0.0))
        if width <= 0:
            return residuals
        half = width / 2.0
        axis = "y" if boundary.get("channelAxis") == "y" else "x"
        blocked = 0 if axis == "y" else 1
        for point in centers.reshape((-1, 3)):
            residuals.append(max(0.0, abs(float(point[blocked])) - half))
    return residuals


def _wall_contact_count(centers: np.ndarray, boundary: dict[str, Any], tol: float = 1e-6) -> int:
    residuals = _wall_residuals(centers, boundary)
    return sum(1 for value in residuals if abs(value) > tol)


def _jacobian_sparsity(
    config: LatticeConfig,
    state: LatticeState,
    edges: list[tuple[int, int, int, int]],
    boundary: dict[str, Any],
) -> lil_matrix:
    n = config.rows * config.cols
    rows: list[set[int]] = []

    def add(*indices: int) -> None:
        rows.append({int(index) for index in indices})

    for r in range(config.rows):
        for c in range(config.cols):
            alpha_idx, theta_idx, x_idx, y_idx, z_idx = _variable_indices(config, r, c)
            if state.removed_mask[r, c]:
                add(alpha_idx)
                add(theta_idx)
                add(z_idx)
                continue
            add(theta_idx, alpha_idx)
            add(alpha_idx)
            if state.locked_mask[r, c]:
                add(alpha_idx)
                add(z_idx)
            if state.position_locked_mask[r, c]:
                add(x_idx)
                add(y_idx)
                add(z_idx)
            add(z_idx)
            add(x_idx)
            add(y_idx)
            add(x_idx)
            add(y_idx)
            add(z_idx)

    for r, c, nr, nc in edges:
        _, theta_i, x_i, y_i, z_i = _variable_indices(config, r, c)
        _, theta_j, x_j, y_j, z_j = _variable_indices(config, nr, nc)
        add(theta_i, theta_j, x_i, y_i, z_i, x_j, y_j, z_j)
        if r == nr:
            add(y_i, y_j)
        if c == nc:
            add(x_i, x_j)
        add(theta_i, theta_j)

    if boundary and boundary.get("wallType") != "inactive":
        mode = boundary.get("mode", "free")
        if mode == "walls":
            wall_axes = []
            if _boundary_value(boundary, "xMin") is not None:
                wall_axes.append(0)
            if _boundary_value(boundary, "xMax") is not None:
                wall_axes.append(0)
            if _boundary_value(boundary, "yMin") is not None:
                wall_axes.append(1)
            if _boundary_value(boundary, "yMax") is not None:
                wall_axes.append(1)
            for r in range(config.rows):
                for c in range(config.cols):
                    _, _, x_idx, y_idx, _ = _variable_indices(config, r, c)
                    for axis in wall_axes:
                        add(x_idx if axis == 0 else y_idx)
        elif mode == "channel":
            width = max(0.0, float(boundary.get("channelWidth", 0.0) or 0.0))
            if width > 0:
                blocked = 0 if boundary.get("channelAxis") == "y" else 1
                for r in range(config.rows):
                    for c in range(config.cols):
                        _, _, x_idx, y_idx, _ = _variable_indices(config, r, c)
                        add(x_idx if blocked == 0 else y_idx)

    matrix = lil_matrix((len(rows), 5 * n), dtype=int)
    for row_index, columns in enumerate(rows):
        for column in columns:
            matrix[row_index, column] = 1
    return matrix


def _seed_vertical_branch(
    config: LatticeConfig,
    state: LatticeState,
    centers: np.ndarray,
    boundary: dict[str, Any],
) -> np.ndarray:
    if not np.any(state.position_locked_mask) and not boundary:
        return centers
    seeded = centers.copy()
    rows = max(1, config.rows - 1)
    cols = max(1, config.cols - 1)
    amplitude = 0.18 * config.cell_size
    for r in range(config.rows):
        for c in range(config.cols):
            if state.removed_mask[r, c] or state.position_locked_mask[r, c]:
                continue
            row_shape = np.sin(np.pi * r / rows) if config.rows > 1 else 1.0
            col_shape = np.sin(np.pi * c / cols) if config.cols > 1 else 1.0
            seeded[r, c, 2] += amplitude * max(row_shape, col_shape)
    return seeded


def solve_constraint_kinematic(
    config: LatticeConfig,
    state: LatticeState,
    boundary: dict[str, Any] | None = None,
    previous: SimulationResult | None = None,
    settings: ConstraintKinematicSettings | None = None,
) -> SimulationResult:
    """Solve RAD cell states with paper-grounded pair constraints.

    This is a kinematic constraint solve, not a calibrated rigid-contact dynamics engine.
    """

    state = state.normalized(config)
    boundary = dict(boundary or {})
    settings = settings or ConstraintKinematicSettings()
    active_mask = ~state.removed_mask
    edges = active_neighbor_edges(config, state)
    initial_alpha, initial_theta, initial_centers = _initial_fields(config, state)
    initial_centers = _seed_vertical_branch(config, state, initial_centers, boundary)
    original_centers_2d = grid_centers(config)
    original_centers = np.dstack(
        [
            original_centers_2d[..., 0],
            original_centers_2d[..., 1],
            np.zeros((config.rows, config.cols), dtype=float),
        ]
    )
    theta_dead_zone = float(normalized_backlash_to_theta_dead_zone(config.backlash))
    x0 = _pack(initial_alpha, initial_theta, initial_centers)
    previous_centers = previous.deformed_centers_3d if previous is not None else initial_centers

    def residual(vector: np.ndarray) -> np.ndarray:
        alpha, theta, centers = _unpack(config, vector)
        out: list[float] = []
        for r in range(config.rows):
            for c in range(config.cols):
                if not active_mask[r, c]:
                    out.extend(
                        [
                            settings.lock_weight * (alpha[r, c] - config.initial_alpha),
                            settings.lock_weight * (theta[r, c] - alpha_to_theta(config.initial_alpha)),
                            settings.position_lock_weight * (centers[r, c, 2]),
                        ]
                    )
                    continue
                out.append(settings.angle_alpha_weight * (theta[r, c] - alpha_to_theta(alpha[r, c])))
                alpha_target_weight = (
                    settings.actuator_alpha_weight
                    if abs(float(state.actuator_grid[r, c])) > 1e-9
                    else settings.alpha_target_weight
                )
                out.append(alpha_target_weight * (alpha[r, c] - initial_alpha[r, c]))
                if state.locked_mask[r, c]:
                    out.append(settings.lock_weight * (alpha[r, c] - state.alpha_grid[r, c]))
                    out.append(settings.lock_weight * (centers[r, c, 2] - state.lock_z_grid[r, c]))
                if state.position_locked_mask[r, c]:
                    out.append(settings.position_lock_weight * (centers[r, c, 0] - original_centers[r, c, 0]))
                    out.append(settings.position_lock_weight * (centers[r, c, 1] - original_centers[r, c, 1]))
                    out.append(settings.position_lock_weight * (centers[r, c, 2] - state.lock_z_grid[r, c]))
                out.append(settings.z_regularization_weight * centers[r, c, 2])
                out.append(settings.xy_regularization_weight * (centers[r, c, 0] - initial_centers[r, c, 0]))
                out.append(settings.xy_regularization_weight * (centers[r, c, 1] - initial_centers[r, c, 1]))
                out.append(settings.continuation_weight * (centers[r, c, 0] - previous_centers[r, c, 0]))
                out.append(settings.continuation_weight * (centers[r, c, 1] - previous_centers[r, c, 1]))
                out.append(settings.continuation_weight * (centers[r, c, 2] - previous_centers[r, c, 2]))
        for r, c, nr, nc in edges:
            pi = centers[r, c]
            pj = centers[nr, nc]
            realized = float(np.linalg.norm(pi - pj))
            desired = pair_distance_from_theta(config, theta[r, c], theta[nr, nc])
            out.append(settings.edge_distance_weight * (realized - desired))
            if r == nr:
                out.append(settings.axis_alignment_weight * (pi[1] - pj[1]))
            if c == nc:
                out.append(settings.axis_alignment_weight * (pi[0] - pj[0]))
            backlash_residual = backlash_activation(theta[r, c] - theta[nr, nc], theta_dead_zone)
            out.append(settings.backlash_weight * float(backlash_residual))
        out.extend(settings.wall_weight * value for value in _wall_residuals(centers, boundary))
        return np.asarray(out, dtype=float)

    lower = np.full_like(x0, -np.inf, dtype=float)
    upper = np.full_like(x0, np.inf, dtype=float)
    n = config.rows * config.cols
    lower[:n] = config.alpha_min
    upper[:n] = config.alpha_max
    result = least_squares(
        residual,
        x0,
        bounds=(lower, upper),
        jac_sparsity=_jacobian_sparsity(config, state, edges, boundary),
        max_nfev=settings.max_nfev,
        xtol=1e-8,
        ftol=1e-8,
        gtol=1e-8,
    )
    alpha, theta, centers = _unpack(config, result.x)
    alpha = np.where(active_mask, alpha, config.initial_alpha)
    theta = np.where(active_mask, theta, alpha_to_theta(config.initial_alpha))
    centers = np.where(active_mask[..., None], centers, original_centers)
    centers_2d = centers[..., :2]
    height = centers[..., 2]
    corners = cell_corners(config, centers_2d, alpha)
    corners_3d = np.zeros((config.rows, config.cols, 5, 3), dtype=float)
    corners_3d[..., :2] = corners
    corners_3d[..., 2] = height[..., None]
    original_corners = cell_corners(
        config, original_centers_2d, np.full((config.rows, config.cols), config.initial_alpha)
    )
    original_corners_3d = np.zeros((config.rows, config.cols, 5, 3), dtype=float)
    original_corners_3d[..., :2] = original_corners

    edge_errors: list[float] = []
    backlash_errors: list[float] = []
    for r, c, nr, nc in edges:
        realized = float(np.linalg.norm(centers[r, c] - centers[nr, nc]))
        desired = pair_distance_from_theta(config, theta[r, c], theta[nr, nc])
        edge_errors.append(realized - desired)
        backlash_errors.append(float(backlash_activation(theta[r, c] - theta[nr, nc], theta_dead_zone)))
    max_edge_error = float(np.max(np.abs(edge_errors))) if edge_errors else 0.0
    mean_edge_error = float(np.mean(np.abs(edge_errors))) if edge_errors else 0.0
    max_backlash_residual = float(np.max(np.abs(backlash_errors))) if backlash_errors else 0.0
    mean_backlash_residual = float(np.mean(np.abs(backlash_errors))) if backlash_errors else 0.0
    residual_success = max(max_edge_error, max_backlash_residual) < 1e-3

    z0 = original_centers_2d[..., 0] + 1j * original_centers_2d[..., 1]
    z1 = centers[..., 0] + 1j * centers[..., 1]
    return SimulationResult(
        config=config,
        state=state,
        alpha=alpha,
        theta_degrees=theta,
        original_centers=original_centers_2d,
        deformed_centers=centers_2d,
        original_corners=original_corners,
        deformed_corners=corners,
        original_centers_3d=original_centers,
        deformed_centers_3d=centers,
        original_corners_3d=original_corners_3d,
        deformed_corners_3d=corners_3d,
        complex_original=z0,
        complex_deformed=z1,
        metadata={
            "model": "constraint_kinematic",
            "success": bool(result.success or residual_success),
            "optimizer_success": bool(result.success),
            "message": result.message,
            "iterations": int(result.nfev),
            "cost": float(result.cost),
            "theta_dead_zone_degrees": theta_dead_zone,
            "edge_count": len(edges),
            "max_edge_error": max_edge_error,
            "mean_edge_error": mean_edge_error,
            "max_backlash_residual": max_backlash_residual,
            "mean_backlash_residual": mean_backlash_residual,
            "boundary_contact_count": _wall_contact_count(centers, boundary),
            "height": height,
            "removed_mask": state.removed_mask,
            "boundary": boundary,
        },
    )
