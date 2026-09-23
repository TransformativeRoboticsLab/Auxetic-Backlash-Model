from __future__ import annotations

import numpy as np

from .coupling import alpha_to_theta
from .geometry import cell_corners, grid_centers
from .models import LatticeConfig, LatticeState, SimulationResult
from .operators import evaluate_programmable_operators


def simulate_kinematic(config: LatticeConfig, state: LatticeState) -> SimulationResult:
    state = state.normalized(config)
    operators = evaluate_programmable_operators(config, state)
    influence = operators["actuator_influence"]
    z_residual = operators["z_residual"]
    alpha = operators["alpha"]
    alpha = np.clip(alpha, config.alpha_min, config.alpha_max)
    theta = alpha_to_theta(alpha)

    original_centers = grid_centers(config)
    deformed_centers = grid_centers(config, alpha)
    height = operators["height"].copy()
    position_locked = state.position_locked_mask & ~state.removed_mask
    deformed_centers[position_locked] = original_centers[position_locked]
    height[position_locked] = 0.0
    original_corners = cell_corners(
        config, original_centers, np.full_like(alpha, config.initial_alpha)
    )
    deformed_corners = cell_corners(config, deformed_centers, alpha)
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
            "die_off": operators["die_off"],
            "z_residual": z_residual,
            "z_die_off": operators["z_die_off"],
            "z_dead_zone": config.pin_hole_clearance,
            "height": height,
            "removed_mask": operators["removed_mask"],
            "topology": operators["topology"],
        },
    )
