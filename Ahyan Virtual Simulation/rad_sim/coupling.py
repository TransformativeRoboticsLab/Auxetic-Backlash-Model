from __future__ import annotations

import numpy as np


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

