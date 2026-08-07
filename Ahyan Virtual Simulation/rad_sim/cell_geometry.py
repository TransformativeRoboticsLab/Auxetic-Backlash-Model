from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .coupling import alpha_to_theta
from .models import LatticeConfig


@dataclass(frozen=True)
class PaperRADReference:
    side_length_mm: float = 35.0
    normalized_backlash: float = 0.1
    poisson_ratio: float = -0.4
    concentric_parts: int = 2
    joints_per_part: int = 4


@dataclass(frozen=True)
class RADJointGeometry:
    part: Literal["outer", "inner"]
    index: int
    position: np.ndarray
    pin_radius: float
    hole_radius: float

    @property
    def clearance(self) -> float:
        return float(self.hole_radius - self.pin_radius)


@dataclass(frozen=True)
class PaperRADCellGeometry:
    reference: PaperRADReference
    center: np.ndarray
    alpha: float
    theta_degrees: float
    outer_part: np.ndarray
    inner_part: np.ndarray
    outer_joints: tuple[RADJointGeometry, ...]
    inner_joints: tuple[RADJointGeometry, ...]
    lock_sites: np.ndarray
    alpha_actuator_axis: np.ndarray
    z_actuator_axis: np.ndarray
    backlash_gap: float
    vertical_free_play: float

    @property
    def joint_count(self) -> int:
        return len(self.outer_joints) + len(self.inner_joints)


PAPER_RAD_REFERENCE = PaperRADReference()


def _square_vertices(
    center: np.ndarray,
    side: float,
    theta_degrees: float,
    z: float,
) -> np.ndarray:
    half = side / 2.0
    base = np.array(
        [[-half, -half], [half, -half], [half, half], [-half, half]],
        dtype=float,
    )
    theta = np.deg2rad(theta_degrees)
    rot = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]],
        dtype=float,
    )
    xy = center[:2] + base @ rot.T
    return np.column_stack([xy, np.full(4, z, dtype=float)])


def build_paper_rad_cell_geometry(
    config: LatticeConfig,
    *,
    alpha: float | None = None,
    center: tuple[float, float] | np.ndarray = (0.0, 0.0),
    z: float = 0.0,
    reference: PaperRADReference = PAPER_RAD_REFERENCE,
) -> PaperRADCellGeometry:
    """Return the normalized paper-grounded RAD unit-cell geometry.

    The source paper specifies rotating-square behavior, two concentric parts,
    four joints per part, b = 0.1 reference backlash, nu = -0.4, and a 35 mm
    prototype side length. Plate thickness and exact CAD clearances are not in
    the extracted text, so this geometry keeps those quantities normalized and
    configurable through ``LatticeConfig``.
    """
    alpha_value = config.initial_alpha if alpha is None else float(alpha)
    alpha_value = float(np.clip(alpha_value, config.alpha_min, config.alpha_max))
    theta = float(alpha_to_theta(alpha_value))
    center_2d = np.asarray(center, dtype=float)
    if center_2d.shape != (2,):
        raise ValueError("center must be a 2D point")
    center_3d = np.array([center_2d[0], center_2d[1], z], dtype=float)

    active_side = 0.62 * config.cell_size * np.sqrt(alpha_value)
    outer_side = active_side + config.backlash * config.cell_size * 0.5
    inner_side = max(0.18 * config.cell_size, active_side * 0.58)
    inner_z = z + 0.08 * config.cell_size

    outer_part = _square_vertices(center_3d, outer_side, 0.0, z)
    inner_part = _square_vertices(center_3d, inner_side, theta, inner_z)
    outer_joints = tuple(
        RADJointGeometry(
            "outer", i, outer_part[i].copy(), config.pin_radius, config.hole_radius
        )
        for i in range(4)
    )
    inner_joints = tuple(
        RADJointGeometry(
            "inner", i, inner_part[i].copy(), config.pin_radius, config.hole_radius
        )
        for i in range(4)
    )
    half_outer = outer_side / 2.0
    alpha_axis_z = z + 0.12 * config.cell_size
    alpha_axis = np.array(
        [
            [center_2d[0] - half_outer, center_2d[1], alpha_axis_z],
            [center_2d[0] + half_outer, center_2d[1], alpha_axis_z],
        ],
        dtype=float,
    )
    z_axis = np.array(
        [center_3d, center_3d + np.array([0.0, 0.0, config.pin_hole_clearance])],
        dtype=float,
    )
    return PaperRADCellGeometry(
        reference=reference,
        center=center_3d,
        alpha=alpha_value,
        theta_degrees=theta,
        outer_part=outer_part,
        inner_part=inner_part,
        outer_joints=outer_joints,
        inner_joints=inner_joints,
        lock_sites=outer_part.copy(),
        alpha_actuator_axis=alpha_axis,
        z_actuator_axis=z_axis,
        backlash_gap=config.backlash * config.cell_size,
        vertical_free_play=config.pin_hole_clearance,
    )
