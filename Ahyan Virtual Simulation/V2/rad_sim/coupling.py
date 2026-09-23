from __future__ import annotations

import numpy as np


def normalized_backlash_to_theta_dead_zone(backlash: np.ndarray | float) -> np.ndarray | float:
    """Paper RAD relation: angular dead-zone Δφ = asin(b/L), in degrees."""
    normalized = np.clip(np.asarray(backlash, dtype=float), 0.0, 1.0)
    return np.degrees(np.arcsin(normalized))


def normalized_backlash_to_alpha_dead_zone(backlash: np.ndarray | float) -> np.ndarray | float:
    """Convert Δφ through theta[degrees] = 70 * alpha - 60."""
    return normalized_backlash_to_theta_dead_zone(backlash) / 70.0


def backlash_activation(x: np.ndarray | float, backlash: float) -> np.ndarray | float:
    """Dead-zone backlash model from the RAD preprint."""
    return np.maximum(0.0, np.asarray(x) - backlash) + np.minimum(
        np.asarray(x) + backlash, 0.0
    )


def alpha_to_theta(alpha: np.ndarray | float) -> np.ndarray | float:
    """Default RAD relation: theta[degrees] = 70 * alpha - 60."""
    return 70.0 * np.asarray(alpha) - 60.0


def theta_to_alpha(theta_degrees: np.ndarray | float) -> np.ndarray | float:
    return (np.asarray(theta_degrees) + 60.0) / 70.0
