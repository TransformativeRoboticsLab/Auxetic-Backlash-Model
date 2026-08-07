from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal

import numpy as np

from .coupling import alpha_to_theta
from .kinematic import simulate_kinematic
from .models import LatticeConfig, LatticeState, SimulationResult


@dataclass(frozen=True)
class PaperRADReference:
    side_length_mm: float = 35.0
    normalized_backlash: float = 0.1
    poisson_ratio: float = -0.4
    concentric_parts: int = 2
    joints_per_part: int = 4
    fabrication_hole_tolerance_mm: float = 0.1

    @property
    def reference_backlash_mm(self) -> float:
        return self.normalized_backlash * self.side_length_mm

    def model_length_to_mm(self, config: LatticeConfig, value: float) -> float:
        return float(value) * self.side_length_mm / config.cell_size

    def mm_to_model_length(self, config: LatticeConfig, value_mm: float) -> float:
        return float(value_mm) * config.cell_size / self.side_length_mm


PAPER_RAD_REFERENCE = PaperRADReference()


@dataclass(frozen=True)
class PaperRADCalibration:
    reference: PaperRADReference
    model_cell_size: float
    mm_per_model_unit: float
    side_length_mm: float
    configured_backlash_mm: float
    reference_backlash_mm: float
    pin_radius_mm: float
    hole_radius_mm: float
    pin_hole_clearance_mm: float
    fabrication_hole_tolerance_mm: float
    fabrication_hole_tolerance_model: float
    poisson_ratio: float
    concentric_parts: int
    joints_per_part: int

    def model_length_to_mm(self, value: float) -> float:
        return float(value) * self.mm_per_model_unit

    def mm_to_model_length(self, value_mm: float) -> float:
        return float(value_mm) / self.mm_per_model_unit


HARDWARE_PROFILE_DIMENSIONS: tuple[str, ...] = (
    "pin_radius_mm",
    "hole_radius_mm",
    "plate_thickness_mm",
    "joint_stack_height_mm",
    "boss_radius_mm",
)
CALIBRATION_VISUAL_FIELDS: tuple[str, ...] = (
    "pin_radius_mm",
    "hole_radius_mm",
    "plate_thickness_mm",
    "joint_stack_height_mm",
)
CALIBRATION_MESH_FIELDS: tuple[str, ...] = HARDWARE_PROFILE_DIMENSIONS
CALIBRATION_SOLVER_GAPS: tuple[str, ...] = (
    "measured axial and hinge stiffness",
    "actuator force/stroke calibration",
    "friction and contact characterization",
    "measured single, pair, and cluster response data",
)
CalibrationReadinessLevel = Literal[
    "paper-scale",
    "partial-measured",
    "visual-calibrated",
    "mesh-calibrated",
]


@dataclass(frozen=True)
class RADHardwareProfile:
    """Measured hardware dimensions for moving from paper-scale to real-cell CAD.

    ``side_length_mm`` and ``fabrication_hole_tolerance_mm`` can come from the
    paper reference. The remaining fields are treated as measured hardware
    dimensions only when they are explicitly provided.
    """

    name: str = "paper-reference"
    source: str = "RAD preprint defaults"
    side_length_mm: float = PAPER_RAD_REFERENCE.side_length_mm
    fabrication_hole_tolerance_mm: float = (
        PAPER_RAD_REFERENCE.fabrication_hole_tolerance_mm
    )
    backlash_mm: float | None = None
    pin_radius_mm: float | None = None
    hole_radius_mm: float | None = None
    plate_thickness_mm: float | None = None
    joint_stack_height_mm: float | None = None
    boss_radius_mm: float | None = None
    notes: str = ""

    def __post_init__(self) -> None:
        if self.side_length_mm <= 0:
            raise ValueError("side_length_mm must be positive")
        if self.fabrication_hole_tolerance_mm < 0:
            raise ValueError("fabrication_hole_tolerance_mm must be non-negative")
        for field_name in ("backlash_mm", *HARDWARE_PROFILE_DIMENSIONS):
            value = getattr(self, field_name)
            if value is not None and value < 0:
                raise ValueError(f"{field_name} must be non-negative when provided")
        if (
            self.pin_radius_mm is not None
            and self.hole_radius_mm is not None
            and self.hole_radius_mm < self.pin_radius_mm
        ):
            raise ValueError("hole_radius_mm must be greater than or equal to pin_radius_mm")

    @property
    def measured_fields(self) -> tuple[str, ...]:
        return tuple(
            field_name
            for field_name in HARDWARE_PROFILE_DIMENSIONS
            if getattr(self, field_name) is not None
        )

    @property
    def missing_fields(self) -> tuple[str, ...]:
        return tuple(
            field_name
            for field_name in HARDWARE_PROFILE_DIMENSIONS
            if getattr(self, field_name) is None
        )

    @property
    def coverage_ratio(self) -> float:
        return len(self.measured_fields) / len(HARDWARE_PROFILE_DIMENSIONS)

    @property
    def pin_hole_clearance_mm(self) -> float | None:
        if self.pin_radius_mm is None or self.hole_radius_mm is None:
            return None
        return self.hole_radius_mm - self.pin_radius_mm

    def to_reference(self) -> PaperRADReference:
        return PaperRADReference(
            side_length_mm=self.side_length_mm,
            normalized_backlash=PAPER_RAD_REFERENCE.normalized_backlash,
            poisson_ratio=PAPER_RAD_REFERENCE.poisson_ratio,
            concentric_parts=PAPER_RAD_REFERENCE.concentric_parts,
            joints_per_part=PAPER_RAD_REFERENCE.joints_per_part,
            fabrication_hole_tolerance_mm=self.fabrication_hole_tolerance_mm,
        )

    def mm_to_model_length(self, config: LatticeConfig, value_mm: float) -> float:
        return float(value_mm) * config.cell_size / self.side_length_mm


@dataclass(frozen=True)
class RADCalibrationReadiness:
    profile_name: str
    level: CalibrationReadinessLevel
    measured_fields: tuple[str, ...]
    missing_fields: tuple[str, ...]
    visual_missing_fields: tuple[str, ...]
    mesh_missing_fields: tuple[str, ...]
    solver_gaps: tuple[str, ...]
    coverage_ratio: float
    visual_ready: bool
    mesh_ready: bool
    solver_ready: bool

    @property
    def summary(self) -> str:
        if self.mesh_ready:
            return "mesh-calibrated geometry; solver still needs physical response calibration"
        if self.visual_ready:
            return "visual-calibrated geometry; mesh export still has missing dimensions"
        if self.measured_fields:
            return "partial measured geometry; solver still needs geometry and response calibration"
        return "paper-scale defaults only; solver still lacks measured hardware geometry"


def calibration_readiness(profile: RADHardwareProfile) -> RADCalibrationReadiness:
    measured = profile.measured_fields
    visual_missing = tuple(
        field for field in CALIBRATION_VISUAL_FIELDS if getattr(profile, field) is None
    )
    mesh_missing = tuple(
        field for field in CALIBRATION_MESH_FIELDS if getattr(profile, field) is None
    )
    visual_ready = len(visual_missing) == 0
    mesh_ready = len(mesh_missing) == 0
    if mesh_ready:
        level: CalibrationReadinessLevel = "mesh-calibrated"
    elif visual_ready:
        level = "visual-calibrated"
    elif measured:
        level = "partial-measured"
    else:
        level = "paper-scale"
    return RADCalibrationReadiness(
        profile_name=profile.name,
        level=level,
        measured_fields=measured,
        missing_fields=profile.missing_fields,
        visual_missing_fields=visual_missing,
        mesh_missing_fields=mesh_missing,
        solver_gaps=CALIBRATION_SOLVER_GAPS,
        coverage_ratio=profile.coverage_ratio,
        visual_ready=visual_ready,
        mesh_ready=mesh_ready,
        solver_ready=False,
    )


def config_with_hardware_profile(
    config: LatticeConfig,
    profile: RADHardwareProfile,
) -> LatticeConfig:
    """Return a config whose available radii/backlash come from measurements."""

    updates: dict[str, float] = {}
    if profile.backlash_mm is not None:
        updates["backlash"] = profile.backlash_mm / profile.side_length_mm
    pin_radius = config.pin_radius
    hole_radius = config.hole_radius
    if profile.pin_radius_mm is not None:
        pin_radius = profile.mm_to_model_length(config, profile.pin_radius_mm)
        updates["pin_radius"] = pin_radius
    if profile.hole_radius_mm is not None:
        hole_radius = profile.mm_to_model_length(config, profile.hole_radius_mm)
    if hole_radius < pin_radius:
        hole_radius = pin_radius
    if profile.hole_radius_mm is not None or "pin_radius" in updates:
        updates["hole_radius"] = hole_radius
    return replace(config, **updates) if updates else config


def hardware_profile_from_config(
    config: LatticeConfig,
    reference: PaperRADReference = PAPER_RAD_REFERENCE,
) -> RADHardwareProfile:
    """Expose the current normalized radii as a configured hardware estimate."""

    calibration = calibrate_paper_rad_config(config, reference)
    return RADHardwareProfile(
        name="current-config-estimate",
        source="normalized simulator controls",
        side_length_mm=calibration.side_length_mm,
        fabrication_hole_tolerance_mm=calibration.fabrication_hole_tolerance_mm,
        backlash_mm=calibration.configured_backlash_mm,
        pin_radius_mm=calibration.pin_radius_mm,
        hole_radius_mm=calibration.hole_radius_mm,
        notes="Derived from current normalized grid values; verify against hardware.",
    )


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
    calibration: PaperRADCalibration
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

    @property
    def backlash_gap_mm(self) -> float:
        return self.calibration.model_length_to_mm(self.backlash_gap)

    @property
    def vertical_free_play_mm(self) -> float:
        return self.calibration.model_length_to_mm(self.vertical_free_play)


@dataclass(frozen=True)
class PaperRADCellRecord:
    index: tuple[int, int]
    geometry: PaperRADCellGeometry
    locked: bool
    command_alpha: float
    command_z: float

    @property
    def is_actuated(self) -> bool:
        return abs(self.command_alpha) > 1e-12 or abs(self.command_z) > 1e-12


@dataclass(frozen=True)
class PaperRADConnectorGeometry:
    first: tuple[int, int]
    second: tuple[int, int]
    axis: Literal["x", "y"]
    start: np.ndarray
    end: np.ndarray
    backlash_gap: float
    vertical_clearance: float
    alpha_delta: float
    height_delta: float

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.end - self.start))


@dataclass(frozen=True)
class PaperRADLatticeGeometry:
    cells: tuple[PaperRADCellRecord, ...]
    connectors: tuple[PaperRADConnectorGeometry, ...]
    centers: np.ndarray
    alpha: np.ndarray
    height: np.ndarray
    locked_mask: np.ndarray

    @property
    def shape(self) -> tuple[int, int]:
        return self.alpha.shape

    @property
    def cell_count(self) -> int:
        return len(self.cells)

    @property
    def connector_count(self) -> int:
        return len(self.connectors)

    @property
    def active_actuator_count(self) -> int:
        return sum(1 for cell in self.cells if cell.is_actuated)

    def cell_at(self, row: int, col: int) -> PaperRADCellRecord:
        rows, cols = self.shape
        if not (0 <= row < rows and 0 <= col < cols):
            raise IndexError("cell index out of range")
        return self.cells[row * cols + col]

def calibrate_paper_rad_config(
    config: LatticeConfig,
    reference: PaperRADReference = PAPER_RAD_REFERENCE,
) -> PaperRADCalibration:
    """Convert the normalized simulator parameters to paper prototype units.

    The paper gives a 35 mm prototype side length, normalized backlash b = 0.1,
    and 0.1 mm hole fabrication tolerance. Pin and hole radii remain current
    simulator parameters, so their millimeter values are configured estimates
    until measured hardware dimensions are added.
    """

    mm_per_model_unit = reference.side_length_mm / config.cell_size
    return PaperRADCalibration(
        reference=reference,
        model_cell_size=config.cell_size,
        mm_per_model_unit=mm_per_model_unit,
        side_length_mm=reference.side_length_mm,
        configured_backlash_mm=reference.model_length_to_mm(
            config, config.backlash * config.cell_size
        ),
        reference_backlash_mm=reference.reference_backlash_mm,
        pin_radius_mm=reference.model_length_to_mm(config, config.pin_radius),
        hole_radius_mm=reference.model_length_to_mm(config, config.hole_radius),
        pin_hole_clearance_mm=reference.model_length_to_mm(
            config, config.pin_hole_clearance
        ),
        fabrication_hole_tolerance_mm=reference.fabrication_hole_tolerance_mm,
        fabrication_hole_tolerance_model=reference.mm_to_model_length(
            config, reference.fabrication_hole_tolerance_mm
        ),
        poisson_ratio=reference.poisson_ratio,
        concentric_parts=reference.concentric_parts,
        joints_per_part=reference.joints_per_part,
    )


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
    calibration = calibrate_paper_rad_config(config, reference)
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
        calibration=calibration,
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


def _side_midpoint(
    cell: PaperRADCellGeometry,
    side: Literal["left", "right", "top", "bottom"],
) -> np.ndarray:
    index_pairs = {
        "bottom": (0, 1),
        "right": (1, 2),
        "top": (2, 3),
        "left": (3, 0),
    }
    i, j = index_pairs[side]
    return (cell.outer_part[i] + cell.outer_part[j]) / 2.0


def _simulation_result(
    config: LatticeConfig,
    state_or_result: LatticeState | SimulationResult | None,
) -> SimulationResult:
    if state_or_result is None:
        return simulate_kinematic(config, LatticeState.uniform(config))
    if isinstance(state_or_result, SimulationResult):
        return state_or_result
    return simulate_kinematic(config, state_or_result)


def build_paper_rad_lattice_geometry(
    config: LatticeConfig,
    state_or_result: LatticeState | SimulationResult | None = None,
    *,
    reference: PaperRADReference = PAPER_RAD_REFERENCE,
) -> PaperRADLatticeGeometry:
    result = _simulation_result(config, state_or_result)
    state = result.state.normalized(config)
    cells: list[PaperRADCellRecord] = []
    geometry_grid: list[list[PaperRADCellGeometry]] = []
    centers_3d = result.deformed_centers_3d
    height = centers_3d[..., 2]

    for r in range(config.rows):
        row: list[PaperRADCellGeometry] = []
        for c in range(config.cols):
            center = centers_3d[r, c]
            geometry = build_paper_rad_cell_geometry(
                config,
                alpha=float(result.alpha[r, c]),
                center=center[:2],
                z=float(center[2]),
                reference=reference,
            )
            row.append(geometry)
            cells.append(
                PaperRADCellRecord(
                    index=(r, c),
                    geometry=geometry,
                    locked=bool(state.locked_mask[r, c]),
                    command_alpha=float(state.actuator_grid[r, c]),
                    command_z=float(state.z_actuator_grid[r, c]),
                )
            )
        geometry_grid.append(row)

    connectors: list[PaperRADConnectorGeometry] = []
    for r in range(config.rows):
        for c in range(config.cols):
            first = geometry_grid[r][c]
            if c + 1 < config.cols:
                second = geometry_grid[r][c + 1]
                connectors.append(
                    PaperRADConnectorGeometry(
                        first=(r, c),
                        second=(r, c + 1),
                        axis="x",
                        start=_side_midpoint(first, "right"),
                        end=_side_midpoint(second, "left"),
                        backlash_gap=config.backlash * config.cell_size,
                        vertical_clearance=config.pin_hole_clearance,
                        alpha_delta=float(result.alpha[r, c + 1] - result.alpha[r, c]),
                        height_delta=float(height[r, c + 1] - height[r, c]),
                    )
                )
            if r + 1 < config.rows:
                second = geometry_grid[r + 1][c]
                connectors.append(
                    PaperRADConnectorGeometry(
                        first=(r, c),
                        second=(r + 1, c),
                        axis="y",
                        start=_side_midpoint(first, "top"),
                        end=_side_midpoint(second, "bottom"),
                        backlash_gap=config.backlash * config.cell_size,
                        vertical_clearance=config.pin_hole_clearance,
                        alpha_delta=float(result.alpha[r + 1, c] - result.alpha[r, c]),
                        height_delta=float(height[r + 1, c] - height[r, c]),
                    )
                )

    return PaperRADLatticeGeometry(
        cells=tuple(cells),
        connectors=tuple(connectors),
        centers=centers_3d.copy(),
        alpha=result.alpha.copy(),
        height=height.copy(),
        locked_mask=state.locked_mask.copy(),
    )
