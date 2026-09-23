from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class LatticeConfig:
    rows: int = 5
    cols: int = 5
    backlash: float = 0.1
    cell_size: float = 1.0
    initial_alpha: float = 1.0
    coupling_gain: float = 0.55
    z_coupling_gain: float = 0.32
    pin_radius: float = 0.18
    hole_radius: float = 0.225
    max_coupling_steps: int | None = None
    alpha_min: float = 0.25
    alpha_max: float = 1.75

    def __post_init__(self) -> None:
        if self.rows <= 0 or self.cols <= 0:
            raise ValueError("rows and cols must be positive")
        if self.cell_size <= 0:
            raise ValueError("cell_size must be positive")
        if self.backlash < 0:
            raise ValueError("backlash must be non-negative")
        if not (0 <= self.coupling_gain <= 1):
            raise ValueError("coupling_gain must be in [0, 1]")
        if not (0 <= self.z_coupling_gain <= 1):
            raise ValueError("z_coupling_gain must be in [0, 1]")
        if self.pin_radius < 0:
            raise ValueError("pin_radius must be non-negative")
        if self.hole_radius < self.pin_radius:
            raise ValueError("hole_radius must be greater than or equal to pin_radius")
        if self.alpha_min <= 0 or self.alpha_min >= self.alpha_max:
            raise ValueError("alpha bounds must satisfy 0 < alpha_min < alpha_max")

    @property
    def pin_hole_clearance(self) -> float:
        return self.hole_radius - self.pin_radius


@dataclass
class LatticeState:
    alpha_grid: np.ndarray
    locked_mask: np.ndarray | None = None
    actuator_grid: np.ndarray | None = None
    z_actuator_grid: np.ndarray | None = None

    @classmethod
    def uniform(cls, config: LatticeConfig, alpha: float | None = None) -> "LatticeState":
        value = config.initial_alpha if alpha is None else alpha
        shape = (config.rows, config.cols)
        return cls(
            alpha_grid=np.full(shape, value, dtype=float),
            locked_mask=np.zeros(shape, dtype=bool),
            actuator_grid=np.zeros(shape, dtype=float),
            z_actuator_grid=np.zeros(shape, dtype=float),
        )

    def normalized(self, config: LatticeConfig) -> "LatticeState":
        shape = (config.rows, config.cols)
        alpha = np.asarray(self.alpha_grid, dtype=float)
        if alpha.shape != shape:
            raise ValueError(f"alpha_grid must have shape {shape}")
        locked = (
            np.zeros(shape, dtype=bool)
            if self.locked_mask is None
            else np.asarray(self.locked_mask, dtype=bool)
        )
        actuators = (
            np.zeros(shape, dtype=float)
            if self.actuator_grid is None
            else np.asarray(self.actuator_grid, dtype=float)
        )
        z_actuators = (
            np.zeros(shape, dtype=float)
            if self.z_actuator_grid is None
            else np.asarray(self.z_actuator_grid, dtype=float)
        )
        if locked.shape != shape:
            raise ValueError(f"locked_mask must have shape {shape}")
        if actuators.shape != shape:
            raise ValueError(f"actuator_grid must have shape {shape}")
        if z_actuators.shape != shape:
            raise ValueError(f"z_actuator_grid must have shape {shape}")
        return LatticeState(alpha.copy(), locked.copy(), actuators.copy(), z_actuators.copy())


@dataclass(frozen=True)
class LoadCase:
    fixed_cells: tuple[tuple[int, int], ...] = ((0, 0),)
    prescribed_displacements: dict[tuple[int, int], tuple[float, ...]] = field(
        default_factory=dict
    )
    external_forces: dict[tuple[int, int], tuple[float, ...]] = field(default_factory=dict)
    axial_stiffness: float = 25.0
    hinge_stiffness: float = 1.0
    lock_stiffness: float = 100.0
    maxiter: int = 500


@dataclass
class SimulationResult:
    config: LatticeConfig
    state: LatticeState
    alpha: np.ndarray
    theta_degrees: np.ndarray
    original_centers: np.ndarray
    deformed_centers: np.ndarray
    original_corners: np.ndarray
    deformed_corners: np.ndarray
    original_centers_3d: np.ndarray
    deformed_centers_3d: np.ndarray
    original_corners_3d: np.ndarray
    deformed_corners_3d: np.ndarray
    complex_original: np.ndarray
    complex_deformed: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)
