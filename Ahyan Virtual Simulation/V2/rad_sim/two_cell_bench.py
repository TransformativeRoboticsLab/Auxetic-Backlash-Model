from __future__ import annotations

import csv
import hashlib
import io
import importlib
import importlib.util
import json
import math
import re
import zipfile
from collections import Counter
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Iterable
from xml.sax.saxutils import escape

from .cad_mesh_intake import SUPPORTED_CAD_BOUND_FORMATS, cad_mesh_audit
from .cell_geometry import RADHardwareProfile
from .coupling import alpha_to_theta, backlash_activation
from .geometry import shared_edge_pitch
from .models import LatticeConfig


TWO_CELL_BENCH_SCHEMA = "rad-sim.two-cell-physical-bench.v1"
TWO_CELL_SWEEP_SCHEMA = "rad-sim.two-cell-backlash-sweep.v1"
TWO_CELL_PACKET_SCHEMA = "rad-sim.two-cell-physical-test-packet.v1"
TWO_CELL_ACTUATION_SWEEP_SCHEMA = "rad-sim.two-cell-actuation-sweep.v1"
TWO_CELL_MJCF_SCHEMA = "rad-sim.two-cell-cad-mjcf-proxy.v1"
TWO_CELL_MUJOCO_RUN_SCHEMA = "rad-sim.two-cell-mujoco-proxy-run.v1"
TWO_CELL_QUASISTATIC_SCHEMA = "rad-sim.two-cell-quasistatic-physics.v1"
TWO_CELL_QUASISTATIC_SWEEP_SCHEMA = "rad-sim.two-cell-quasistatic-sweep.v1"
TWO_CELL_EXTERNAL_TEMPLATE_SCHEMA = "rad-sim.two-cell-external-results-template.v1"
TWO_CELL_EXTERNAL_COMPARISON_SCHEMA = "rad-sim.two-cell-external-comparison.v1"
TWO_CELL_PHYSICS_VALIDATION_SCHEMA = "rad-sim.two-cell-physics-validation-report.v1"
CAD_RAD_CELL_LAYOUT_SCHEMA = "rad-sim.cad-rad-cell-layout.v1"
TWO_CELL_CONNECTOR_CONTACT_SCHEMA = "rad-sim.two-cell-connector-contact-report.v1"
TWO_CELL_CONNECTOR_CONTACT_SWEEP_SCHEMA = "rad-sim.two-cell-connector-contact-sweep.v1"
TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA = "rad-sim.two-cell-physical-simulation-suite.v1"
TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA = "rad-sim.two-cell-physical-response-atlas.v1"
TWO_CELL_CONNECTOR_MEASUREMENT_TEMPLATE_SCHEMA = "rad-sim.two-cell-connector-measurement-template.v1"
TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA = "rad-sim.two-cell-physical-fidelity-matrix.v1"
TWO_CELL_CONTACT_PHASE_MAP_SCHEMA = "rad-sim.two-cell-contact-phase-map.v1"
TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA = "rad-sim.two-cell-radius-backlash-phase-diagram.v1"
TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA = "rad-sim.two-cell-radius-backlash-transition-report.v1"
TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA = "rad-sim.two-cell-fidelity-matrix-measurement-template.v1"
TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA = "rad-sim.two-cell-external-fidelity-matrix-manifest.v1"
TWO_CELL_EXTERNAL_FIDELITY_MJCF_RUN_SCHEMA = "rad-sim.two-cell-external-fidelity-matrix-mjcf-run.v1"
TWO_CELL_SEGMENTED_CAD_READINESS_SCHEMA = "rad-sim.two-cell-segmented-cad-readiness.v1"
TWO_CELL_SEGMENTED_CAD_INTAKE_TEMPLATES_SCHEMA = "rad-sim.two-cell-segmented-cad-intake-templates.v1"
TWO_CELL_SEGMENTED_CAD_INTAKE_VALIDATION_SCHEMA = "rad-sim.two-cell-segmented-cad-intake-validation.v1"
TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA = "rad-sim.two-cell-cad-contact-decomposition.v1"
TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA = "rad-sim.two-cell-exact-contact-handoff-plan.v1"
CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA = "rad-sim.cad-rad-cell-archive-audit.v1"
CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA = "rad-sim.cad-rad-cell-reference-profile.v1"
ONE_CELL_CAD_EXPORT_AUDIT_SCHEMA = "rad-sim.one-cell-cad-export-audit.v1"


@dataclass(frozen=True)
class CADRadCellReference:
    schema: str = "rad-sim.cad-rad-cell-reference.v1"
    source: str = "Autodesk A360 public share https://a360.co/4bMlzip"
    fusion_archive: str = "assets/cad/RADs_unit_cell.f3d"
    viewer_summary: str = "assets/cad/RADs_unit_cell_reference.json"
    preview_image: str = "assets/cad/RADs_unit_cell_preview.png"
    units: str = "mm"
    width_x_mm: float = 55.604331129396634
    width_y_mm: float = 55.60433117189373
    height_z_mm: float = 19.99999621152464
    nominal_hole_diameter_mm: float = 3.4
    archive_display_name: str = "RADs unit cell"
    linked_cell_display_name: str = "RADs free cell 4mm tall 3.4mm hole"
    viewer_title: str = "RADs unit cell - AUTODESK FUSION"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "source": self.source,
            "localFiles": {
                "fusionArchive": self.fusion_archive,
                "viewerSummary": self.viewer_summary,
                "preview": self.preview_image,
            },
            "units": self.units,
            "boundingBoxMm": {
                "widthX": self.width_x_mm,
                "widthY": self.width_y_mm,
                "heightZ": self.height_z_mm,
            },
            "nominalHoleDiameterMm": self.nominal_hole_diameter_mm,
            "fusionManifestLabels": {
                "archiveDisplayName": self.archive_display_name,
                "linkedCellDisplayName": self.linked_cell_display_name,
            },
            "viewerAccess": {
                "shareUrl": "https://a360.co/4bMlzip",
                "resolvedTitle": self.viewer_title,
                "downloadFormatsObserved": [
                    "Fusion Archive",
                    "Inventor 2025",
                    "IGES",
                    "SAT",
                    "SMT",
                    "STEP",
                    "DWG",
                    "DXF",
                    "STL",
                    "FBX",
                    "SketchUp",
                    "OBJ",
                ],
                "downloadRequiresEmailExport": True,
            },
            "visibleFeatures": [
                "central circular crown and screw stack",
                "eight radial arms",
                "outer circular pin-hole pads",
                "single-cell monolithic upper/lower body fragments",
            ],
            "fragments": [
                {
                    "id": 0,
                    "name": "RADs free cell 4mm tall 3.4mm hole body",
                    "zMin": 3.9999988079071045,
                    "zMax": 8.000000953674316,
                },
                {"id": 1, "name": "cell lower body", "zMin": 0.0, "zMax": 4.0},
                {
                    "id": 2,
                    "name": "18-8 stainless screw",
                    "zMin": -9.999998092651367,
                    "zMax": 9.999998092651367,
                },
            ],
            "physicalInterpretation": (
                "CAD-derived visual/dimensional reference. The reduced two-cell model "
                "does not replace segmented rigid-body contact simulation."
            ),
        }


@dataclass(frozen=True)
class TwoCellBenchControls:
    alpha_command: float = -0.35
    z_command: float = 0.35
    hole_sweep_max: float = 0.5
    hole_sweep_steps: int = 9
    left_position_locked: bool = True
    right_locked: bool = False
    right_position_locked: bool = False


CAD_RAD_CELL_REFERENCE = CADRadCellReference()
SEGMENTED_CAD_INTAKE_DIR = "assets/cad/segmented"
CAD_RAD_SITE_ORDER = (
    ("n", 90.0),
    ("ne", 45.0),
    ("e", 0.0),
    ("se", -45.0),
    ("s", -90.0),
    ("sw", -135.0),
    ("w", 180.0),
    ("nw", 135.0),
)
CAD_RAD_X_CONNECTOR_PAIRS = (
    ("upper", "ne", "nw"),
    ("middle", "e", "w"),
    ("lower", "se", "sw"),
)
SEGMENTED_CAD_EXPECTED_ASSETS = {
    "upper_free_cell_body": ("upper_free_cell_body.step", "upper_free_cell_body.stl", "upper_free_cell_body.obj"),
    "lower_cell_body": ("lower_cell_body.step", "lower_cell_body.stl", "lower_cell_body.obj"),
    "screw_pin_body": ("screw_pin_body.step", "screw_pin_body.stl", "screw_pin_body.obj"),
    "radial_pad_hole_surfaces": ("radial_pad_hole_surfaces.step", "radial_pad_hole_surfaces.json"),
    "central_alpha_rotation_axis": ("joint_axes.json", "central_alpha_rotation_axis.json"),
    "vertical_slide_axis": ("joint_axes.json", "vertical_slide_axis.json"),
    "eight_pin_hole_contact_axes": ("joint_axes.json", "pin_hole_contact_axes.json"),
    "two_cell_connector_pairs": ("two_cell_connector_pairs.json",),
    "mass_inertia": ("mass_inertia.json",),
    "friction_contact_parameters": ("contact_parameters.json",),
    "lock_crown_geometry": ("lock_crown_geometry.json", "lock_crown_geometry.csv"),
    "actuator_force_displacement": ("actuator_force_displacement.json", "actuator_force_displacement.csv"),
    "bench_coordinate_truth": ("bench_coordinate_truth.json", "bench_coordinate_truth.csv"),
    "MuJoCo": ("mujoco/rad_two_cell.xml", "rad_two_cell.xml"),
    "Gazebo/Ignition": ("gazebo/rad_two_cell.sdf", "rad_two_cell.sdf", "rad_two_cell.urdf"),
    "Isaac Sim": ("isaac/rad_two_cell.usd", "rad_two_cell.usd"),
}
EXACT_CONTACT_HANDOFF_REQUIRED_CAD_ASSETS = (
    "upper_free_cell_body",
    "lower_cell_body",
    "screw_pin_body",
    "radial_pad_hole_surfaces",
    "central_alpha_rotation_axis",
    "vertical_slide_axis",
    "eight_pin_hole_contact_axes",
    "two_cell_connector_pairs",
    "mass_inertia",
    "friction_contact_parameters",
)
EXACT_CONTACT_HANDOFF_REQUIRED_MEASURED_INPUTS = (
    "lock_crown_geometry",
    "actuator_force_displacement",
    "bench_coordinate_truth",
    "pin_hole_friction",
    "normal_contact_stiffness",
)


def _archive_text_tokens(data: bytes, *, limit: int = 80) -> list[str]:
    tokens: set[str] = set()
    for encoding in ("utf-8", "utf-16le"):
        text = data.decode(encoding, errors="ignore")
        for match in re.findall(r"[A-Za-z0-9][A-Za-z0-9_ .:/?=&+\\-]{3,96}", text):
            token = " ".join(match.strip().split())
            if token:
                tokens.add(token)
    return sorted(tokens)[:limit]


def cad_rad_cell_archive_audit(
    path: str | Path | None = None,
    *,
    root: Path | None = None,
) -> dict[str, Any]:
    """Audit the local Fusion archive evidence for the real one-cell CAD file."""

    project_root = root or Path(__file__).resolve().parents[1]
    archive_path = Path(path or CAD_RAD_CELL_REFERENCE.fusion_archive)
    resolved = archive_path if archive_path.is_absolute() else project_root / archive_path
    local_record = _asset_record(str(archive_path).replace("\\", "/"), root=project_root)
    if not resolved.exists():
        return {
            "schema": CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA,
            "source": CAD_RAD_CELL_REFERENCE.source,
            "archive": local_record,
            "summary": {
                "archivePresent": False,
                "zipReadable": False,
                "physicalAccuracyValidated": False,
                "status": "missing-fusion-archive",
                "missingEvidence": ["fusionArchive"],
            },
            "entries": [],
            "brepEntries": [],
            "previewEntries": [],
            "manifestEntries": [],
            "textEvidence": [],
            "limitations": ["No local Fusion archive was available to inspect."],
        }

    archive_bytes = resolved.read_bytes()
    sha256 = hashlib.sha256(archive_bytes).hexdigest()
    zip_readable = zipfile.is_zipfile(resolved)
    entries: list[dict[str, Any]] = []
    brep_entries: list[dict[str, Any]] = []
    preview_entries: list[dict[str, Any]] = []
    manifest_entries: list[dict[str, Any]] = []
    text_evidence: list[dict[str, Any]] = []
    linked_labels: set[str] = set()
    if zip_readable:
        with zipfile.ZipFile(resolved) as archive:
            for info in archive.infolist():
                entry = {
                    "name": info.filename,
                    "bytes": info.file_size,
                    "compressedBytes": info.compress_size,
                    "isDirectory": info.is_dir(),
                }
                entries.append(entry)
                lower = info.filename.lower()
                if "breps.blobparts" in lower and not info.is_dir():
                    brep_entries.append(entry)
                if "/previews/" in lower and not info.is_dir():
                    preview_entries.append(entry)
                if info.filename.endswith("Manifest.dat") or info.filename.endswith(".json"):
                    manifest_entries.append(entry)
                if not info.is_dir() and (
                    info.filename.endswith("Manifest.dat")
                    or info.filename.endswith(".json")
                    or info.filename.endswith("MetaStream.dat")
                    or info.filename.endswith("BulkStream.dat")
                ):
                    data = archive.read(info.filename)
                    tokens = _archive_text_tokens(data)
                    interesting = [
                        token
                        for token in tokens
                        if any(
                            needle in token.lower()
                            for needle in ("rad", "cell", "hole", "screw", "body", "fusion", "brep", "3.4", "4mm")
                        )
                    ][:20]
                    if interesting:
                        for token in interesting:
                            if "rad" in token.lower() or "hole" in token.lower():
                                linked_labels.add(token)
                        text_evidence.append(
                            {
                                "entry": info.filename,
                                "tokens": interesting,
                            }
                        )

    expected_entries = {
        "Manifest.dat",
        "FusionAssetName[Active]/Manifest.dat",
        "FusionAssetName[Active]/Previews/small.png",
    }
    entry_names = {entry["name"] for entry in entries}
    missing_expected_entries = sorted(expected_entries - entry_names)
    exact_contact_ready = False
    missing_evidence = []
    if not zip_readable:
        missing_evidence.append("zipReadableFusionArchive")
    if not brep_entries:
        missing_evidence.append("brepBlobEntries")
    if not preview_entries:
        missing_evidence.append("previewEntry")
    if missing_expected_entries:
        missing_evidence.append("expectedFusionArchiveEntries")
    missing_evidence.extend(
        [
            "segmentedBodyMeshExport",
            "exactJointAxes",
            "exactPinHoleContactSurfaces",
            "materialContactParameters",
            "benchCoordinateTruth",
        ]
    )
    return {
        "schema": CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA,
        "source": CAD_RAD_CELL_REFERENCE.source,
        "archive": {
            **local_record,
            "sha256": sha256,
            "zipReadable": zip_readable,
        },
        "summary": {
            "archivePresent": True,
            "zipReadable": zip_readable,
            "entryCount": len(entries),
            "fileEntryCount": sum(1 for entry in entries if not entry["isDirectory"]),
            "brepEntryCount": len(brep_entries),
            "previewEntryCount": len(preview_entries),
            "manifestEntryCount": len(manifest_entries),
            "linkedLabelCount": len(linked_labels),
            "expectedEntriesPresent": not missing_expected_entries,
            "canDeriveVisualReference": zip_readable and bool(preview_entries or brep_entries),
            "canAttemptExactRigidBodyContact": exact_contact_ready,
            "physicalAccuracyValidated": False,
            "status": "archive-audited-needs-segmentation",
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "entries": entries,
        "brepEntries": brep_entries,
        "previewEntries": preview_entries,
        "manifestEntries": manifest_entries,
        "missingExpectedEntries": missing_expected_entries,
        "textEvidence": text_evidence,
        "linkedLabels": sorted(linked_labels),
        "limitations": [
            "Fusion archive entries prove local CAD evidence exists, not exact simulation fidelity.",
            "BREP blobs are not segmented into named moving bodies or pin-hole contact surfaces here.",
            "Exact real-life mechanics still require segmented exports, external contact simulation, and bench comparison.",
        ],
    }


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, float(value)))


def _dead_zone(value: float, gap: float) -> float:
    return float(backlash_activation(float(value), float(gap)))


def _two_cell_rest_pitch(config: LatticeConfig, left_alpha: float, right_alpha: float) -> float:
    return float(shared_edge_pitch(config, left_alpha, right_alpha))


def _bounded_pattern_search(
    objective,
    x0: list[float],
    bounds: list[tuple[float, float]],
    *,
    tolerance: float,
    maxiter: int,
) -> dict[str, Any]:
    values = [_clamp(float(value), float(lower), float(upper)) for value, (lower, upper) in zip(x0, bounds)]
    best = float(objective(values))
    spans = [max(0.0, float(upper) - float(lower)) for lower, upper in bounds]
    steps = [max(span * 0.25, 0.05) if span > 0.0 else 0.0 for span in spans]
    threshold = max(math.sqrt(max(float(tolerance), 1e-15)), 1e-6)
    iteration = 0

    for iteration in range(1, max(1, int(maxiter)) + 1):
        improved = False
        for index, step in enumerate(steps):
            if step <= 0.0:
                continue
            lower, upper = bounds[index]
            for direction in (-1.0, 1.0):
                candidate = list(values)
                candidate[index] = _clamp(candidate[index] + direction * step, float(lower), float(upper))
                if candidate[index] == values[index]:
                    continue
                score = float(objective(candidate))
                if math.isfinite(score) and score + max(float(tolerance), 1e-12) < best:
                    values = candidate
                    best = score
                    improved = True
        if not improved:
            steps = [step * 0.5 for step in steps]
            if max(steps, default=0.0) <= threshold:
                break

    return {
        "x": values,
        "success": math.isfinite(best),
        "message": "bounded pattern-search fallback completed",
        "nit": iteration,
        "fun": best,
    }


def _controls(controls: TwoCellBenchControls | None = None, **overrides: Any) -> TwoCellBenchControls:
    base = asdict(controls or TwoCellBenchControls())
    base.update({key: value for key, value in overrides.items() if value is not None})
    return TwoCellBenchControls(
        alpha_command=_clamp(base["alpha_command"], -1.2, 1.2),
        z_command=_clamp(base["z_command"], -1.4, 1.4),
        hole_sweep_max=_clamp(base["hole_sweep_max"], 0.02, 1.25),
        hole_sweep_steps=max(3, min(41, int(round(base["hole_sweep_steps"])))),
        left_position_locked=bool(base["left_position_locked"]),
        right_locked=bool(base["right_locked"]),
        right_position_locked=bool(base["right_position_locked"]),
    )


def _radii(config: LatticeConfig, hole_radius: float | None = None, pin_radius: float | None = None) -> dict[str, float]:
    pin = max(0.0, float(config.pin_radius if pin_radius is None else pin_radius))
    hole = max(pin, float(config.hole_radius if hole_radius is None else hole_radius))
    return {"pinRadius": pin, "holeRadius": hole, "pinHoleClearance": hole - pin}


def _dead_zones(config: LatticeConfig, radii: dict[str, float]) -> dict[str, float]:
    clearance = radii["pinHoleClearance"]
    return {
        "backlash": config.backlash,
        "clearance": clearance,
        "alphaDeadZone": config.backlash + 0.5 * clearance,
        "verticalDeadZone": clearance,
        "effectiveBacklash": config.backlash + clearance,
    }


def _cad_rad_layout_dimensions(
    *,
    width_mm: float,
    height_mm: float,
    pin_radius_mm: float,
    hole_radius_mm: float,
) -> dict[str, float]:
    pad_radius = max(hole_radius_mm * 1.85, width_mm * 0.055)
    hub_radius = max(hole_radius_mm * 1.6, width_mm * 0.065)
    site_radius = max(0.0, 0.5 * width_mm - pad_radius)
    arm_width = max(0.85, width_mm * 0.018)
    return {
        "widthMm": width_mm,
        "heightMm": height_mm,
        "bodyThicknessMm": 4.0,
        "lowerBodyThicknessMm": 4.0,
        "stackHeightMm": height_mm,
        "nominalHoleRadiusMm": 0.5 * CAD_RAD_CELL_REFERENCE.nominal_hole_diameter_mm,
        "pinRadiusMm": pin_radius_mm,
        "holeRadiusMm": hole_radius_mm,
        "pinHoleClearanceMm": max(0.0, hole_radius_mm - pin_radius_mm),
        "padRadiusMm": pad_radius,
        "hubRadiusMm": hub_radius,
        "siteRadiusMm": site_radius,
        "connectorPitchMm": 2.0 * site_radius,
        "armStartRadiusMm": hub_radius * 0.78,
        "armWidthMm": arm_width,
        "screwRadiusMm": max(0.65 * pin_radius_mm, 0.8),
    }


def _cad_rad_site_records(dimensions: dict[str, float]) -> list[dict[str, float | str]]:
    sites: list[dict[str, float | str]] = []
    radius = dimensions["siteRadiusMm"]
    pad_radius = dimensions["padRadiusMm"]
    half_width = 0.5 * dimensions["widthMm"]
    for name, angle_degrees in CAD_RAD_SITE_ORDER:
        angle_radians = math.radians(angle_degrees)
        x = radius * math.cos(angle_radians)
        y = radius * math.sin(angle_radians)
        sites.append(
            {
                "name": name,
                "angleDegrees": angle_degrees,
                "xMm": x,
                "yMm": y,
                "padRadiusMm": pad_radius,
                "holeRadiusMm": dimensions["holeRadiusMm"],
                "pinRadiusMm": dimensions["pinRadiusMm"],
                "outerEnvelopeRadiusMm": min(half_width, math.hypot(x, y) + pad_radius),
            }
        )
    return sites


def cad_rad_cell_layout(
    config: LatticeConfig | None = None,
    *,
    hole_radius: float | None = None,
    pin_radius: float | None = None,
) -> dict[str, Any]:
    """Dimension-consistent one-cell CAD abstraction from the A360 reference."""

    config = config or LatticeConfig(rows=1, cols=1)
    cad = CAD_RAD_CELL_REFERENCE
    radii = _radii(config, hole_radius=hole_radius, pin_radius=pin_radius)
    reference_hole_radius = max(float(config.hole_radius), 1e-9)
    nominal_hole_radius_mm = 0.5 * cad.nominal_hole_diameter_mm
    width_mm = min(cad.width_x_mm, cad.width_y_mm)
    pin_radius_mm = nominal_hole_radius_mm * radii["pinRadius"] / reference_hole_radius
    hole_radius_mm = nominal_hole_radius_mm * radii["holeRadius"] / reference_hole_radius
    dimensions = _cad_rad_layout_dimensions(
        width_mm=width_mm,
        height_mm=cad.height_z_mm,
        pin_radius_mm=pin_radius_mm,
        hole_radius_mm=hole_radius_mm,
    )
    sites = _cad_rad_site_records(dimensions)
    half_width = 0.5 * dimensions["widthMm"]
    max_envelope = max((abs(float(site["xMm"])) + dimensions["padRadiusMm"] for site in sites), default=0.0)
    max_envelope = max(
        max_envelope,
        max((abs(float(site["yMm"])) + dimensions["padRadiusMm"] for site in sites), default=0.0),
    )
    width_scale = max(dimensions["widthMm"], 1e-9)
    pitch_scale = max(dimensions["connectorPitchMm"], 1e-9)
    model_scale = max(float(config.cell_size), 1e-9) / width_scale
    connector_pairs = [
        {
            "name": name,
            "leftSite": left_site,
            "rightSite": right_site,
            "axis": "x",
            "leftPinGeom": f"left_cell_{left_site}_pin",
            "rightHoleGeom": f"right_cell_{right_site}_hole_clearance",
        }
        for name, left_site, right_site in CAD_RAD_X_CONNECTOR_PAIRS
    ]
    return {
        "schema": CAD_RAD_CELL_LAYOUT_SCHEMA,
        "cadReference": cad.to_dict(),
        "sourceEvidence": [
            "A360 Fusion archive bounding box",
            "linked cell label: RADs free cell 4mm tall 3.4mm hole",
            "eight radial pad sites visible in extracted preview",
        ],
        "limitations": [
            "Site layout is dimension-consistent and topology-matched, but not exact BREP vertex extraction.",
            "Fusion BREP streams still need segmentation/export before exact rigid-body contact simulation.",
        ],
        "dimensionsMm": dimensions,
        "padSites": sites,
        "twoCellConnectorPairs": connector_pairs,
        "envelopeCheck": {
            "maxOuterEnvelopeMm": max_envelope,
            "halfWidthMm": half_width,
            "matchesBoundingBox": abs(max_envelope - half_width) <= 1e-9,
        },
        "visualModel": {
            "siteRadiusToWidth": dimensions["siteRadiusMm"] / width_scale,
            "padRadiusToWidth": dimensions["padRadiusMm"] / width_scale,
            "hubRadiusToWidth": dimensions["hubRadiusMm"] / width_scale,
            "holeRadiusToWidth": dimensions["holeRadiusMm"] / width_scale,
            "pinRadiusToWidth": dimensions["pinRadiusMm"] / width_scale,
            "armStartRadiusToWidth": dimensions["armStartRadiusMm"] / width_scale,
            "armWidthToWidth": dimensions["armWidthMm"] / width_scale,
            "bodyThicknessToWidth": dimensions["bodyThicknessMm"] / width_scale,
            "stackHeightToWidth": dimensions["stackHeightMm"] / width_scale,
            "siteRadiusToPitch": dimensions["siteRadiusMm"] / pitch_scale,
            "padRadiusToPitch": dimensions["padRadiusMm"] / pitch_scale,
            "hubRadiusToPitch": dimensions["hubRadiusMm"] / pitch_scale,
            "holeRadiusToPitch": dimensions["holeRadiusMm"] / pitch_scale,
            "pinRadiusToPitch": dimensions["pinRadiusMm"] / pitch_scale,
            "holeToPadRadiusRatio": dimensions["holeRadiusMm"] / max(dimensions["padRadiusMm"], 1e-9),
            "pinToHoleRadiusRatio": dimensions["pinRadiusMm"] / max(dimensions["holeRadiusMm"], 1e-9),
            "screwToHubRadiusRatio": dimensions["screwRadiusMm"] / max(dimensions["hubRadiusMm"], 1e-9),
            "hubToPadRadiusRatio": dimensions["hubRadiusMm"] / max(dimensions["padRadiusMm"], 1e-9),
            "armStartRadiusToPitch": dimensions["armStartRadiusMm"] / pitch_scale,
            "armWidthToPitch": dimensions["armWidthMm"] / pitch_scale,
            "bodyThicknessToPitch": dimensions["bodyThicknessMm"] / pitch_scale,
            "stackHeightToPitch": dimensions["stackHeightMm"] / pitch_scale,
            "modelUnitsPerMm": model_scale,
        },
    }


def cad_rad_cell_reference_profile(
    config: LatticeConfig | None = None,
    *,
    hole_radius: float | None = None,
    pin_radius: float | None = None,
) -> dict[str, Any]:
    """A360 one-cell dimensions packaged as a hardware-profile candidate.

    This profile is intentionally marked as CAD-derived. It is useful for visual
    scale and packet handoff, but it is not a substitute for segmented bodies,
    contact surfaces, or bench-calibrated physical parameters.
    """

    config = config or LatticeConfig(rows=1, cols=1)
    layout = cad_rad_cell_layout(config, hole_radius=hole_radius, pin_radius=pin_radius)
    dimensions = layout["dimensionsMm"]
    hardware_profile = RADHardwareProfile(
        name="a360-rads-unit-cell-visual-profile",
        source=CAD_RAD_CELL_REFERENCE.source,
        side_length_mm=float(dimensions["widthMm"]),
        fabrication_hole_tolerance_mm=0.1,
        backlash_mm=None,
        pin_radius_mm=float(dimensions["pinRadiusMm"]),
        hole_radius_mm=float(dimensions["holeRadiusMm"]),
        plate_thickness_mm=float(dimensions["bodyThicknessMm"]),
        joint_stack_height_mm=float(dimensions["stackHeightMm"]),
        boss_radius_mm=float(dimensions["hubRadiusMm"]),
        notes=(
            "A360 CAD-derived visual profile from the unit-cell bounding box, "
            "3.4 mm nominal hole label, archive z-fragments, and current normalized pin/hole ratio. "
            "Use for CAD-like display and external-engine handoff only until segmented contact and bench data exist."
        ),
    ).to_dict()
    return {
        "schema": CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA,
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "cadLayout": layout,
        "hardwareProfile": hardware_profile,
        "sourceEvidence": [
            "A360 public share metadata identifies the file as RADs unit cell.",
            "Fusion archive/preview provide bounding box and visible radial-pad topology.",
            "Linked model label reports a 3.4 mm hole.",
            "Archive fragments show 4 mm upper body, 4 mm lower body, and 20 mm screw stack envelope.",
        ],
        "assumptionLabels": {
            "pinRadiusMm": "derived from the current normalized pin/hole ratio and the 3.4 mm nominal hole label",
            "holeRadiusMm": "derived from the 3.4 mm nominal hole label",
            "plateThicknessMm": "from archive fragment z-ranges, not segmented body mass properties",
            "bossRadiusMm": "visual hub radius inferred from bounding-box-scaled layout",
            "backlashMm": "not claimed from CAD; keep configurable or measure directly",
        },
        "claimBoundary": {
            "allowedClaim": "A360-dimensioned visual/profile candidate for real-cell UI and external packet generation",
            "blockedClaim": "exact real-cell dynamics, contact, friction, or lock behavior without segmented CAD and bench validation",
        },
    }


def export_cad_rad_cell_reference_profile_json(
    config: LatticeConfig | None = None,
    *,
    hole_radius: float | None = None,
    pin_radius: float | None = None,
) -> str:
    return json.dumps(
        cad_rad_cell_reference_profile(config, hole_radius=hole_radius, pin_radius=pin_radius),
        indent=2,
    )


def export_one_cell_cad_export_audit_json(
    config: LatticeConfig | None = None,
    **audit_overrides: Any,
) -> str:
    return json.dumps(one_cell_cad_export_audit(config, **audit_overrides), indent=2)


def _asset_record(path: str, *, root: Path | None = None) -> dict[str, Any]:
    base = root or Path(__file__).resolve().parents[1]
    candidate = Path(path)
    resolved = candidate if candidate.is_absolute() else base / candidate
    return {
        "path": path,
        "exists": resolved.exists(),
        "bytes": resolved.stat().st_size if resolved.exists() else None,
    }


def one_cell_cad_export_audit(
    config: LatticeConfig | None = None,
    *,
    cad_dir: str | Path = "assets/cad",
    basename: str = "RADs_unit_cell",
) -> dict[str, Any]:
    """Find a one-cell STEP/STL/OBJ export and derive bounds/profile evidence."""

    config = config or LatticeConfig(rows=1, cols=2)
    project_root = Path(__file__).resolve().parents[1]
    source_dir = Path(cad_dir)
    if not source_dir.is_absolute():
        source_dir = project_root / source_dir
    candidates = [source_dir / f"{basename}{suffix}" for suffix in SUPPORTED_CAD_BOUND_FORMATS]

    def _display_path(path: Path) -> str:
        try:
            return str(path.relative_to(project_root)).replace("\\", "/")
        except ValueError:
            return str(path).replace("\\", "/")

    candidate_records = [_asset_record(_display_path(path), root=project_root) for path in candidates]
    detected = [path for path in candidates if path.exists()]
    preferred = sorted(
        detected,
        key=lambda path: {".step": 0, ".stp": 1, ".obj": 2, ".stl": 3}.get(path.suffix.lower(), 9),
    )
    mesh_audit = None
    error = ""
    if preferred:
        try:
            profile = cad_rad_cell_reference_profile(config)["hardwareProfile"]["dimensionsMm"]
            mesh_audit = cad_mesh_audit(
                preferred[0],
                name="one-cell-cad-export-derived",
                source=_display_path(preferred[0]),
                pin_radius_mm=profile.get("pinRadiusMm"),
                hole_radius_mm=profile.get("holeRadiusMm"),
                backlash_mm=profile.get("backlashMm"),
                boss_radius_mm=profile.get("bossRadiusMm"),
            )
        except Exception as exc:
            error = str(exc)
    export_detected = bool(mesh_audit)
    missing_evidence = []
    if not export_detected:
        missing_evidence.append("oneCellStepOrMeshExport")
    if error:
        missing_evidence.append("cadExportParseFailure")
    missing_evidence.extend(
        [
            "segmentedMovingBodies",
            "jointAxes",
            "pinHoleContactSurfaces",
            "benchCoordinateHoldout",
        ]
    )
    return {
        "schema": ONE_CELL_CAD_EXPORT_AUDIT_SCHEMA,
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "supportedFormats": list(SUPPORTED_CAD_BOUND_FORMATS),
        "search": {
            "directory": _display_path(source_dir),
            "basename": basename,
            "candidates": candidate_records,
        },
        "selectedExport": _display_path(preferred[0]) if preferred else "",
        "cadMeshAudit": mesh_audit,
        "summary": {
            "status": "one-cell-cad-export-detected-needs-segmentation"
            if export_detected
            else "needs-one-cell-step-or-mesh-export",
            "exportDetected": export_detected,
            "candidateCount": len(candidates),
            "detectedCandidateCount": len(detected),
            "physicalAccuracyValidated": False,
            "error": error,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimBoundary": {
            "allowedClaim": "one-cell CAD export bounds can calibrate visual scale and hardware profile fields",
            "blockedClaim": "exact real-cell contact dynamics until the CAD is segmented and validated against bench data",
            "physicalAccuracyValidated": False,
        },
    }


def _candidate_asset_records(
    task_id: str,
    *,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
    root: Path | None = None,
) -> list[dict[str, Any]]:
    candidates = SEGMENTED_CAD_EXPECTED_ASSETS.get(task_id, ())
    prefix = Path(intake_dir)
    return [_asset_record(str(prefix / candidate).replace("\\", "/"), root=root) for candidate in candidates]


def _detected_assets(task_id: str, *, root: Path) -> list[dict[str, Any]]:
    return [record for record in _candidate_asset_records(task_id, root=root) if record["exists"]]


def _with_intake_status(
    task: dict[str, Any],
    *,
    task_id: str | None = None,
    detected_status: str = "asset-detected-needs-validation",
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
    root: Path,
) -> dict[str, Any]:
    key = task_id or str(task.get("id", task.get("engine", "")))
    expected = _candidate_asset_records(key, intake_dir=intake_dir, root=root)
    detected = [record for record in expected if record["exists"]]
    updated = dict(task)
    updated["expectedAssets"] = expected
    updated["detectedAssets"] = detected
    if detected:
        updated["status"] = detected_status
    return updated


def two_cell_segmented_cad_readiness_report(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
) -> dict[str, Any]:
    """List the CAD segmentation evidence needed for physically exact two-cell contact."""

    config = replace(config or LatticeConfig(rows=1, cols=2), rows=1, cols=2)
    controls = _controls(controls)
    layout = cad_rad_cell_layout(config)
    cad = CAD_RAD_CELL_REFERENCE
    project_root = Path(__file__).resolve().parents[1]
    local_assets = {
        "fusionArchive": _asset_record(cad.fusion_archive, root=project_root),
        "viewerSummary": _asset_record(cad.viewer_summary, root=project_root),
        "previewImage": _asset_record(cad.preview_image, root=project_root),
    }
    one_cell_export = one_cell_cad_export_audit(config)
    fragments = cad.to_dict()["fragments"]
    body_tasks = [
        {
            "id": "upper_free_cell_body",
            "requiredFor": "moving auxetic cell body and alpha/z actuation",
            "currentEvidence": "Fusion manifest fragment label and z-range only",
            "neededAsset": "separated STEP/STL/OBJ mesh with local frame, mass, inertia, and collision surface",
            "status": "missing-segmented-body",
        },
        {
            "id": "lower_cell_body",
            "requiredFor": "stacked lower body / crown support and ground-plane contacts",
            "currentEvidence": "Fusion manifest fragment label and z-range only",
            "neededAsset": "separated body mesh plus crown/lock interface surfaces",
            "status": "missing-segmented-body",
        },
        {
            "id": "screw_pin_body",
            "requiredFor": "pin-hole backlash, vertical free play, screw-stack contact",
            "currentEvidence": "Fusion manifest screw fragment label and z-range only",
            "neededAsset": "cylindrical pin/screw collision mesh with axis, radius, length, material",
            "status": "missing-segmented-body",
        },
        {
            "id": "radial_pad_hole_surfaces",
            "requiredFor": "eight local pin/hole contact sites and two-cell connector pairing",
            "currentEvidence": "derived eight-pad layout from A360 bounding box and nominal 3.4 mm hole label",
            "neededAsset": "hole wall surfaces or analytic cylinders for each pad site",
            "status": "derived-not-extracted",
        },
    ]
    body_tasks = [
        _with_intake_status(task, intake_dir=intake_dir, root=project_root)
        for task in body_tasks
    ]
    joint_tasks = [
        {
            "id": "central_alpha_rotation_axis",
            "neededEvidence": "axis origin/direction for rotating-square auxetic degree of freedom",
            "status": "missing-measured-axis",
        },
        {
            "id": "vertical_slide_axis",
            "neededEvidence": "allowed screw/pin z travel and hard stops from CAD or bench measurement",
            "status": "missing-measured-axis",
        },
        {
            "id": "eight_pin_hole_contact_axes",
            "neededEvidence": "per-pad cylinder axes/radii and local wall normals",
            "count": len(layout["padSites"]),
            "status": "derived-not-extracted",
        },
        {
            "id": "two_cell_connector_pairs",
            "neededEvidence": "left/right mating site transforms for upper, middle, lower connector pairs",
            "count": len(layout["twoCellConnectorPairs"]),
            "status": "dimension-derived",
        },
    ]
    joint_tasks = [
        _with_intake_status(task, intake_dir=intake_dir, root=project_root)
        for task in joint_tasks
    ]
    physics_tasks = [
        {
            "id": "mass_inertia",
            "neededEvidence": "mass and inertia tensor for each separated body",
            "status": "missing",
        },
        {
            "id": "friction_contact_parameters",
            "neededEvidence": "pin/hole friction, restitution, normal stiffness, damping, solver impedance",
            "status": "missing",
        },
        {
            "id": "lock_crown_geometry",
            "neededEvidence": "actual crown angle, arm width correction, lock fall-out behavior, and fixture placement",
            "status": "missing",
        },
        {
            "id": "actuator_force_displacement",
            "neededEvidence": "alpha and z actuator command to displacement/force curves",
            "status": "missing",
        },
        {
            "id": "bench_coordinate_truth",
            "neededEvidence": "tracked 3D marker coordinates for one-cell and two-cell actuation/lock cases",
            "status": "missing",
        },
    ]
    physics_tasks = [
        _with_intake_status(task, intake_dir=intake_dir, root=project_root)
        for task in physics_tasks
    ]
    engine_tasks = [
        {
            "engine": "MuJoCo",
            "neededAsset": "MJCF with separated bodies, joints, analytic pins/holes, contact pairs, materials",
            "currentProxy": "two_cell_cad_proxy.xml",
            "status": "proxy-only",
        },
        {
            "engine": "Gazebo/Ignition",
            "neededAsset": "URDF/SDF with collision meshes, joints, inertials, and contact parameters",
            "currentProxy": None,
            "status": "not-exported",
        },
        {
            "engine": "Isaac Sim",
            "neededAsset": "USD articulation or rigid-body assembly with collision approximations",
            "currentProxy": None,
            "status": "not-exported",
        },
    ]
    engine_tasks = [
        _with_intake_status(
            task,
            task_id=str(task["engine"]),
            detected_status="engine-handoff-detected-needs-run",
            intake_dir=intake_dir,
            root=project_root,
        )
        for task in engine_tasks
    ]
    required_assets = [*body_tasks, *joint_tasks, *physics_tasks, *engine_tasks]
    local_asset_count = sum(1 for asset in local_assets.values() if asset["exists"])
    detected_intake_count = sum(1 for task in required_assets if task.get("detectedAssets"))
    satisfied_statuses = {
        "dimension-derived",
        "asset-detected-needs-validation",
        "engine-handoff-detected-needs-run",
    }
    missing_required_count = sum(1 for task in required_assets if task["status"] not in satisfied_statuses)
    segmented_ready = (
        local_assets["fusionArchive"]["exists"]
        and missing_required_count == 0
    )
    readiness_status = "segmented-cad-assets-detected-needs-validation" if segmented_ready else "needs-segmented-cad-export"
    return {
        "schema": TWO_CELL_SEGMENTED_CAD_READINESS_SCHEMA,
        "model": "a360-one-cell-to-two-cell-segmented-cad-readiness",
        "cadReference": cad.to_dict(),
        "cadLayout": layout,
        "controls": {
            "alphaCommand": controls.alpha_command,
            "zCommand": controls.z_command,
            "leftPositionLocked": controls.left_position_locked,
            "rightLocked": controls.right_locked,
            "rightPositionLocked": controls.right_position_locked,
        },
        "localAssets": local_assets,
        "oneCellCadExportAudit": one_cell_export,
        "segmentedCadIntake": {
            "directory": intake_dir,
            "expectedAssetMap": {
                task_id: [str(Path(intake_dir) / candidate).replace("\\", "/") for candidate in candidates]
                for task_id, candidates in SEGMENTED_CAD_EXPECTED_ASSETS.items()
            },
            "detectedAssetCount": detected_intake_count,
        },
        "fusionFragments": fragments,
        "bodySegmentationTasks": body_tasks,
        "jointAndContactTasks": joint_tasks,
        "physicsParameterTasks": physics_tasks,
        "engineHandoffTasks": engine_tasks,
        "summary": {
            "status": readiness_status,
            "localAssetCount": local_asset_count,
            "oneCellCadExportDetected": bool(one_cell_export["summary"]["exportDetected"]),
            "detectedSegmentedCadAssetCount": detected_intake_count,
            "requiredTaskCount": len(required_assets),
            "missingRequiredTaskCount": missing_required_count,
            "segmentedCadReady": segmented_ready,
            "physicalAccuracyValidated": False,
            "canRunReducedProxy": True,
            "canRunExactRigidBodyContact": False,
            "nextAction": (
                "Run external rigid-body contact and compare against bench coordinates."
                if segmented_ready
                else "Export separated STEP/STL/OBJ bodies from Fusion with pin/hole axes, "
                "then attach friction/contact and bench coordinate data."
            ),
        },
        "claimLabels": {
            "availableCad": "A360 one-cell archive, preview, bounding box, fragment labels, and dimension-derived layout",
            "missingCad": "segmented rigid bodies, exact contact surfaces, joint axes, material/contact parameters",
            "simulation": "current two-cell physics remains a reduced proxy until this manifest is satisfied",
        },
    }


def _json_template_text(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2)


def _csv_template_text(rows: list[dict[str, Any]], fieldnames: list[str]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def _resolved_asset_path(record: dict[str, Any], *, root: Path) -> Path:
    path = Path(str(record["path"]))
    return path if path.is_absolute() else root / path


def _field_is_filled(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (int, float, bool)):
        return True
    if isinstance(value, (list, tuple)):
        return bool(value) and all(_field_is_filled(item) for item in value)
    if isinstance(value, dict):
        return bool(value) and all(_field_is_filled(item) for item in value.values())
    return True


def _validation_entry(
    *,
    template_id: str,
    path: str,
    kind: str,
    description: str,
    exists: bool,
    payload: Any = None,
    rows: list[dict[str, Any]] | None = None,
    parse_error: str | None = None,
    errors: list[str] | None = None,
    filled_fields: int = 0,
    required_fields: int = 0,
) -> dict[str, Any]:
    issues = list(errors or [])
    if parse_error:
        status = "parse-error"
        issues.insert(0, parse_error)
    elif not exists:
        status = "missing"
        issues.insert(0, "file missing")
    elif issues or filled_fields < required_fields:
        status = "incomplete"
    else:
        status = "complete"
    return {
        "templateId": template_id,
        "path": path,
        "kind": kind,
        "description": description,
        "exists": exists,
        "status": status,
        "requiredFieldCount": required_fields,
        "filledFieldCount": filled_fields,
        "errorCount": len(issues),
        "errors": issues,
        "rowCount": len(rows or []),
        "topLevelKeys": sorted(payload.keys()) if isinstance(payload, dict) else [],
    }


def _parse_json_intake(path: Path) -> tuple[dict[str, Any] | None, str | None]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        return None, f"read failed: {exc}"
    except json.JSONDecodeError as exc:
        return None, f"json parse failed: {exc.msg}"
    if not isinstance(payload, dict):
        return None, "json root must be an object"
    return payload, None


def _parse_csv_intake(path: Path) -> tuple[list[dict[str, Any]] | None, str | None]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        return None, f"read failed: {exc}"
    try:
        rows = list(csv.DictReader(io.StringIO(text)))
    except csv.Error as exc:
        return None, f"csv parse failed: {exc}"
    return rows, None


def _validation_counter() -> dict[str, Any]:
    return {"required": 0, "filled": 0, "errors": []}


def _require(counter: dict[str, Any], label: str, value: Any) -> None:
    counter["required"] += 1
    if _field_is_filled(value):
        counter["filled"] += 1
    else:
        counter["errors"].append(f"{label} missing")


def _require_present(counter: dict[str, Any], label: str, value: Any) -> None:
    counter["required"] += 1
    present = value is not None and not (isinstance(value, (list, tuple, dict)) and len(value) == 0)
    if present:
        counter["filled"] += 1
    else:
        counter["errors"].append(f"{label} missing")


def _row_by_key(rows: Iterable[dict[str, Any]], key: str, value: str) -> dict[str, Any] | None:
    for row in rows:
        if str(row.get(key, "")) == value:
            return row
    return None


def _validate_joint_axes_payload(payload: dict[str, Any], layout: dict[str, Any]) -> dict[str, Any]:
    counter = _validation_counter()
    axes = payload.get("axes")
    pin_axes = payload.get("pinHoleAxes")
    if not isinstance(axes, list):
        counter["errors"].append("axes must be a list")
        axes = []
    if not isinstance(pin_axes, list):
        counter["errors"].append("pinHoleAxes must be a list")
        pin_axes = []
    for axis_id in ("central_alpha_rotation_axis", "vertical_slide_axis"):
        axis = _row_by_key(axes, "id", axis_id)
        _require_present(counter, f"axes.{axis_id}", axis)
        if isinstance(axis, dict):
            _require(counter, f"axes.{axis_id}.originMm", axis.get("originMm"))
            _require(counter, f"axes.{axis_id}.direction", axis.get("direction"))
            if axis_id == "vertical_slide_axis":
                _require(counter, "axes.vertical_slide_axis.travelMinMm", axis.get("travelMinMm"))
                _require(counter, "axes.vertical_slide_axis.travelMaxMm", axis.get("travelMaxMm"))
    expected_sites = [site["name"] for site in layout["padSites"]]
    if len(pin_axes) != len(expected_sites):
        counter["errors"].append(f"pinHoleAxes expected {len(expected_sites)} rows, found {len(pin_axes)}")
    for site in expected_sites:
        row = _row_by_key(pin_axes, "site", site)
        _require_present(counter, f"pinHoleAxes.{site}", row)
        if isinstance(row, dict):
            _require(counter, f"pinHoleAxes.{site}.originMm", row.get("originMm"))
            _require(counter, f"pinHoleAxes.{site}.direction", row.get("direction"))
            _require(counter, f"pinHoleAxes.{site}.pinRadiusMm", row.get("pinRadiusMm"))
            _require(counter, f"pinHoleAxes.{site}.holeRadiusMm", row.get("holeRadiusMm"))
    return counter


def _validate_radial_surfaces_payload(payload: dict[str, Any], layout: dict[str, Any]) -> dict[str, Any]:
    counter = _validation_counter()
    surfaces = payload.get("surfaces")
    if not isinstance(surfaces, list):
        counter["errors"].append("surfaces must be a list")
        surfaces = []
    expected_sites = [site["name"] for site in layout["padSites"]]
    if len(surfaces) != len(expected_sites):
        counter["errors"].append(f"surfaces expected {len(expected_sites)} rows, found {len(surfaces)}")
    for site in expected_sites:
        row = _row_by_key(surfaces, "site", site)
        _require_present(counter, f"surfaces.{site}", row)
        if isinstance(row, dict):
            _require(counter, f"surfaces.{site}.centerMm", row.get("centerMm"))
            _require(counter, f"surfaces.{site}.axis", row.get("axis"))
            _require(counter, f"surfaces.{site}.holeRadiusMm", row.get("holeRadiusMm"))
            _require(
                counter,
                f"surfaces.{site}.surfaceAssetOrFusionFace",
                row.get("surfaceAsset") or row.get("fusionBodyOrFaceId"),
            )
    return counter


def _validate_connector_pairs_payload(payload: dict[str, Any], layout: dict[str, Any]) -> dict[str, Any]:
    counter = _validation_counter()
    pairs = payload.get("pairs")
    if not isinstance(pairs, list):
        counter["errors"].append("pairs must be a list")
        pairs = []
    expected_pairs = list(layout["twoCellConnectorPairs"])
    if len(pairs) != len(expected_pairs):
        counter["errors"].append(f"pairs expected {len(expected_pairs)} rows, found {len(pairs)}")
    for expected in expected_pairs:
        name = expected["name"]
        row = _row_by_key(pairs, "name", name)
        _require_present(counter, f"pairs.{name}", row)
        if isinstance(row, dict):
            if row.get("leftSite") != expected["leftSite"]:
                counter["errors"].append(f"pairs.{name}.leftSite expected {expected['leftSite']}")
            if row.get("rightSite") != expected["rightSite"]:
                counter["errors"].append(f"pairs.{name}.rightSite expected {expected['rightSite']}")
            _require(counter, f"pairs.{name}.axis", row.get("axis"))
            _require(counter, f"pairs.{name}.measuredLeftOriginMm", row.get("measuredLeftOriginMm"))
            _require(counter, f"pairs.{name}.measuredRightOriginMm", row.get("measuredRightOriginMm"))
    return counter


def _validate_mass_inertia_payload(payload: dict[str, Any]) -> dict[str, Any]:
    counter = _validation_counter()
    bodies = payload.get("bodies")
    if not isinstance(bodies, list):
        counter["errors"].append("bodies must be a list")
        bodies = []
    for body_id in ("upper_free_cell_body", "lower_cell_body", "screw_pin_body"):
        row = _row_by_key(bodies, "bodyId", body_id)
        _require_present(counter, f"bodies.{body_id}", row)
        if isinstance(row, dict):
            _require(counter, f"bodies.{body_id}.massKg", row.get("massKg"))
            _require(counter, f"bodies.{body_id}.centerOfMassMm", row.get("centerOfMassMm"))
            _require(counter, f"bodies.{body_id}.inertiaTensorKgMm2", row.get("inertiaTensorKgMm2"))
    return counter


def _validate_contact_parameters_payload(payload: dict[str, Any]) -> dict[str, Any]:
    counter = _validation_counter()
    contacts = payload.get("contacts")
    if not isinstance(contacts, list):
        counter["errors"].append("contacts must be a list")
        contacts = []
    for contact_id in ("pin_hole_wall", "lock_crown_contact", "ground_plane"):
        row = _row_by_key(contacts, "contactId", contact_id)
        _require_present(counter, f"contacts.{contact_id}", row)
        if isinstance(row, dict):
            _require(counter, f"contacts.{contact_id}.normalStiffness", row.get("normalStiffness"))
            _require(counter, f"contacts.{contact_id}.normalDamping", row.get("normalDamping"))
            _require(counter, f"contacts.{contact_id}.staticFriction", row.get("staticFriction"))
            _require(counter, f"contacts.{contact_id}.dynamicFriction", row.get("dynamicFriction"))
    return counter


def _validate_lock_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    counter = _validation_counter()
    if not rows:
        counter["errors"].append("lock_crown_geometry.csv must include at least one lock row")
    for index, row in enumerate(rows, start=1):
        prefix = f"rows.{index}"
        _require(counter, f"{prefix}.lockId", row.get("lockId"))
        _require(counter, f"{prefix}.nominalCrownAngleDeg", row.get("nominalCrownAngleDeg"))
        _require(counter, f"{prefix}.effectiveLatticeThetaDeg", row.get("effectiveLatticeThetaDeg"))
        _require(counter, f"{prefix}.armWidthCorrectionDeg", row.get("armWidthCorrectionDeg"))
        _require(counter, f"{prefix}.falloutObserved", row.get("falloutObserved"))
    return counter


def _validate_actuator_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    counter = _validation_counter()
    actuator_ids = {row.get("actuatorId", "") for row in rows}
    for actuator_id in ("alpha_actuator", "z_actuator"):
        if actuator_id not in actuator_ids:
            counter["errors"].append(f"{actuator_id} rows missing")
    for index, row in enumerate(rows, start=1):
        prefix = f"rows.{index}"
        _require(counter, f"{prefix}.actuatorId", row.get("actuatorId"))
        _require(counter, f"{prefix}.command", row.get("command"))
        _require(counter, f"{prefix}.measuredDisplacementMm", row.get("measuredDisplacementMm"))
        _require(counter, f"{prefix}.measuredForceN", row.get("measuredForceN"))
    return counter


def _validate_bench_truth_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    counter = _validation_counter()
    if not rows:
        counter["errors"].append("bench_coordinate_truth.csv must include marker coordinate rows")
    for index, row in enumerate(rows, start=1):
        prefix = f"rows.{index}"
        _require(counter, f"{prefix}.caseId", row.get("caseId"))
        _require(counter, f"{prefix}.bodyId", row.get("bodyId"))
        _require(counter, f"{prefix}.markerId", row.get("markerId"))
        _require(counter, f"{prefix}.xMm", row.get("xMm"))
        _require(counter, f"{prefix}.yMm", row.get("yMm"))
        _require(counter, f"{prefix}.zMm", row.get("zMm"))
        _require(counter, f"{prefix}.measurementSource", row.get("measurementSource"))
    return counter


def two_cell_segmented_cad_intake_templates(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
) -> dict[str, Any]:
    """Generate starter files for segmented CAD/contact intake."""

    config = replace(config or LatticeConfig(rows=1, cols=2), rows=1, cols=2)
    controls = _controls(controls)
    layout = cad_rad_cell_layout(config)
    dimensions = layout["dimensionsMm"]
    site_rows = list(layout["padSites"])
    connector_pairs = list(layout["twoCellConnectorPairs"])
    joint_axes = {
        "schema": "rad-sim.segmented-cad-joint-axes-template.v1",
        "units": "mm",
        "cadSource": CAD_RAD_CELL_REFERENCE.source,
        "instructions": [
            "Replace derived origins/directions with values measured from segmented Fusion exports.",
            "Use a right-handed cell-local frame with z along the screw/pin stack.",
        ],
        "axes": [
            {
                "id": "central_alpha_rotation_axis",
                "originMm": [0.0, 0.0, 0.0],
                "direction": [0.0, 0.0, 1.0],
                "source": "derived-placeholder",
            },
            {
                "id": "vertical_slide_axis",
                "originMm": [0.0, 0.0, 0.0],
                "direction": [0.0, 0.0, 1.0],
                "travelMinMm": "",
                "travelMaxMm": "",
                "source": "measure-from-CAD-or-bench",
            },
        ],
        "pinHoleAxes": [
            {
                "site": site["name"],
                "originMm": [site["xMm"], site["yMm"], 0.0],
                "direction": [0.0, 0.0, 1.0],
                "pinRadiusMm": site["pinRadiusMm"],
                "holeRadiusMm": site["holeRadiusMm"],
                "source": "layout-derived-placeholder",
            }
            for site in site_rows
        ],
    }
    radial_holes = {
        "schema": "rad-sim.segmented-cad-radial-pad-hole-surfaces-template.v1",
        "units": "mm",
        "instructions": [
            "Replace analytic placeholders with exact hole-wall surface IDs or mesh filenames from Fusion.",
            "Keep site names n/ne/e/se/s/sw/w/nw so connector pairing remains stable.",
        ],
        "surfaces": [
            {
                "site": site["name"],
                "centerMm": [site["xMm"], site["yMm"], 0.0],
                "axis": [0.0, 0.0, 1.0],
                "holeRadiusMm": site["holeRadiusMm"],
                "surfaceAsset": "",
                "fusionBodyOrFaceId": "",
            }
            for site in site_rows
        ],
    }
    connector_template = {
        "schema": "rad-sim.segmented-cad-two-cell-connector-pairs-template.v1",
        "units": "mm",
        "pairs": [
            {
                "name": pair["name"],
                "leftSite": pair["leftSite"],
                "rightSite": pair["rightSite"],
                "axis": pair["axis"],
                "leftBody": "left_upper_free_cell_body",
                "rightBody": "right_upper_free_cell_body",
                "measuredLeftOriginMm": "",
                "measuredRightOriginMm": "",
                "source": "layout-derived-placeholder",
            }
            for pair in connector_pairs
        ],
    }
    mass_inertia = {
        "schema": "rad-sim.segmented-cad-mass-inertia-template.v1",
        "units": {"mass": "kg", "length": "mm", "inertia": "kg*mm^2"},
        "bodies": [
            {
                "bodyId": body_id,
                "massKg": "",
                "centerOfMassMm": ["", "", ""],
                "inertiaTensorKgMm2": [["", "", ""], ["", "", ""], ["", "", ""]],
                "source": "Fusion physical properties or measured mass",
            }
            for body_id in ("upper_free_cell_body", "lower_cell_body", "screw_pin_body")
        ],
    }
    contact_parameters = {
        "schema": "rad-sim.segmented-cad-contact-parameters-template.v1",
        "units": {"stiffness": "N/m", "damping": "N*s/m", "friction": "unitless"},
        "contacts": [
            {
                "contactId": "pin_hole_wall",
                "normalStiffness": "",
                "normalDamping": "",
                "staticFriction": "",
                "dynamicFriction": "",
                "restitution": "",
                "source": "bench fit or material test",
            },
            {
                "contactId": "lock_crown_contact",
                "normalStiffness": "",
                "normalDamping": "",
                "staticFriction": "",
                "dynamicFriction": "",
                "restitution": "",
                "source": "lock experiment fit",
            },
            {
                "contactId": "ground_plane",
                "normalStiffness": "",
                "normalDamping": "",
                "staticFriction": "",
                "dynamicFriction": "",
                "restitution": "",
                "source": "bench fixture",
            },
        ],
    }
    lock_rows = [
        {
            "lockId": "lock_30deg",
            "nominalCrownAngleDeg": 30,
            "effectiveLatticeThetaDeg": "",
            "armWidthCorrectionDeg": "",
            "falloutObserved": "",
            "notes": "",
        },
        {
            "lockId": "lock_40deg",
            "nominalCrownAngleDeg": 40,
            "effectiveLatticeThetaDeg": "",
            "armWidthCorrectionDeg": "",
            "falloutObserved": "",
            "notes": "",
        },
    ]
    actuator_rows = [
        {
            "actuatorId": "alpha_actuator",
            "command": command,
            "measuredDisplacementMm": "",
            "measuredForceN": "",
            "measuredAlpha": "",
            "measuredThetaDeg": "",
            "notes": "",
        }
        for command in (-0.5, -0.25, 0.0, 0.25, 0.5)
    ] + [
        {
            "actuatorId": "z_actuator",
            "command": command,
            "measuredDisplacementMm": "",
            "measuredForceN": "",
            "measuredAlpha": "",
            "measuredThetaDeg": "",
            "notes": "",
        }
        for command in (-0.45, 0.0, 0.45)
    ]
    bench_rows = []
    for case in two_cell_physical_simulation_suite(config, controls)["rows"]:
        for body_id in ("left", "right"):
            bench_rows.append(
                {
                    "caseId": case["caseId"],
                    "bodyId": body_id,
                    "markerId": f"{body_id}_center",
                    "xMm": "",
                    "yMm": "",
                    "zMm": "",
                    "alpha": "",
                    "thetaDeg": "",
                    "lockHeldObserved": "",
                    "measurementSource": "",
                    "notes": "",
                }
            )
    templates = [
        {
            "path": f"{intake_dir}/joint_axes.json",
            "kind": "json",
            "description": "cell-local alpha/z axes and eight pin-hole axes",
            "content": _json_template_text(joint_axes),
        },
        {
            "path": f"{intake_dir}/radial_pad_hole_surfaces.json",
            "kind": "json",
            "description": "exact or analytic hole-wall surfaces for all eight pad sites",
            "content": _json_template_text(radial_holes),
        },
        {
            "path": f"{intake_dir}/two_cell_connector_pairs.json",
            "kind": "json",
            "description": "upper/middle/lower mating connector site transforms",
            "content": _json_template_text(connector_template),
        },
        {
            "path": f"{intake_dir}/mass_inertia.json",
            "kind": "json",
            "description": "mass, center of mass, and inertia tensors for segmented bodies",
            "content": _json_template_text(mass_inertia),
        },
        {
            "path": f"{intake_dir}/contact_parameters.json",
            "kind": "json",
            "description": "normal/friction/contact parameters for external engines",
            "content": _json_template_text(contact_parameters),
        },
        {
            "path": f"{intake_dir}/lock_crown_geometry.csv",
            "kind": "csv",
            "description": "nominal/effective lock crown angles and fallout observations",
            "content": _csv_template_text(
                lock_rows,
                [
                    "lockId",
                    "nominalCrownAngleDeg",
                    "effectiveLatticeThetaDeg",
                    "armWidthCorrectionDeg",
                    "falloutObserved",
                    "notes",
                ],
            ),
        },
        {
            "path": f"{intake_dir}/actuator_force_displacement.csv",
            "kind": "csv",
            "description": "alpha/z actuator command-to-displacement/force measurements",
            "content": _csv_template_text(
                actuator_rows,
                [
                    "actuatorId",
                    "command",
                    "measuredDisplacementMm",
                    "measuredForceN",
                    "measuredAlpha",
                    "measuredThetaDeg",
                    "notes",
                ],
            ),
        },
        {
            "path": f"{intake_dir}/bench_coordinate_truth.csv",
            "kind": "csv",
            "description": "tracked center/marker coordinates for physical two-cell cases",
            "content": _csv_template_text(
                bench_rows,
                [
                    "caseId",
                    "bodyId",
                    "markerId",
                    "xMm",
                    "yMm",
                    "zMm",
                    "alpha",
                    "thetaDeg",
                    "lockHeldObserved",
                    "measurementSource",
                    "notes",
                ],
            ),
        },
    ]
    return {
        "schema": TWO_CELL_SEGMENTED_CAD_INTAKE_TEMPLATES_SCHEMA,
        "model": "two-cell-segmented-cad-intake-template-files",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "cadLayout": layout,
        "intakeDirectory": intake_dir,
        "templateCount": len(templates),
        "templates": templates,
        "summary": {
            "jsonTemplateCount": sum(1 for item in templates if item["kind"] == "json"),
            "csvTemplateCount": sum(1 for item in templates if item["kind"] == "csv"),
            "benchCoordinateRows": len(bench_rows),
            "actuatorRows": len(actuator_rows),
            "lockRows": len(lock_rows),
        },
    }


def two_cell_segmented_cad_intake_validation_report(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
) -> dict[str, Any]:
    """Parse filled segmented CAD/contact intake files and block placeholder physics claims."""

    config = replace(config or LatticeConfig(rows=1, cols=2), rows=1, cols=2)
    controls = _controls(controls)
    project_root = Path(__file__).resolve().parents[1]
    layout = cad_rad_cell_layout(config)
    readiness = two_cell_segmented_cad_readiness_report(config, controls, intake_dir=intake_dir)
    templates = two_cell_segmented_cad_intake_templates(config, controls, intake_dir=intake_dir)
    validators = {
        "joint_axes.json": lambda payload: _validate_joint_axes_payload(payload, layout),
        "radial_pad_hole_surfaces.json": lambda payload: _validate_radial_surfaces_payload(payload, layout),
        "two_cell_connector_pairs.json": lambda payload: _validate_connector_pairs_payload(payload, layout),
        "mass_inertia.json": _validate_mass_inertia_payload,
        "contact_parameters.json": _validate_contact_parameters_payload,
    }
    csv_validators = {
        "lock_crown_geometry.csv": _validate_lock_rows,
        "actuator_force_displacement.csv": _validate_actuator_rows,
        "bench_coordinate_truth.csv": _validate_bench_truth_rows,
    }
    template_rows = []
    for item in templates["templates"]:
        template_id = Path(str(item["path"])).name
        record = _asset_record(str(item["path"]), root=project_root)
        path = _resolved_asset_path(record, root=project_root)
        payload: dict[str, Any] | None = None
        rows: list[dict[str, Any]] | None = None
        parse_error = None
        counter = _validation_counter()
        if record["exists"]:
            if item["kind"] == "json":
                payload, parse_error = _parse_json_intake(path)
                if payload is not None:
                    counter = validators[template_id](payload)
            else:
                rows, parse_error = _parse_csv_intake(path)
                if rows is not None:
                    counter = csv_validators[template_id](rows)
        template_rows.append(
            _validation_entry(
                template_id=template_id,
                path=str(item["path"]),
                kind=item["kind"],
                description=item["description"],
                exists=bool(record["exists"]),
                payload=payload,
                rows=rows,
                parse_error=parse_error,
                errors=list(counter["errors"]),
                filled_fields=int(counter["filled"]),
                required_fields=int(counter["required"]),
            )
        )

    body_asset_ids = ("upper_free_cell_body", "lower_cell_body", "screw_pin_body")
    body_asset_rows = []
    for body_id in body_asset_ids:
        candidates = _candidate_asset_records(body_id, intake_dir=intake_dir, root=project_root)
        detected = [record for record in candidates if record["exists"] and (record["bytes"] or 0) > 0]
        body_asset_rows.append(
            {
                "assetId": body_id,
                "requiredFor": "segmented rigid body collision and inertia frame",
                "status": "complete" if detected else "missing",
                "expectedAssets": candidates,
                "detectedAssets": detected,
            }
        )
    radial_surface_assets = _candidate_asset_records("radial_pad_hole_surfaces", intake_dir=intake_dir, root=project_root)
    detected_radial_surface_assets = [
        record for record in radial_surface_assets if record["exists"] and (record["bytes"] or 0) > 0
    ]
    engine_asset_rows = []
    for engine_id in ("MuJoCo", "Gazebo/Ignition", "Isaac Sim"):
        candidates = _candidate_asset_records(engine_id, intake_dir=intake_dir, root=project_root)
        detected = [record for record in candidates if record["exists"] and (record["bytes"] or 0) > 0]
        engine_asset_rows.append(
            {
                "engine": engine_id,
                "status": "ready-for-external-run" if detected else "missing",
                "expectedAssets": candidates,
                "detectedAssets": detected,
            }
        )

    complete_template_count = sum(1 for row in template_rows if row["status"] == "complete")
    existing_template_count = sum(1 for row in template_rows if row["exists"])
    missing_template_count = sum(1 for row in template_rows if row["status"] == "missing")
    parse_error_count = sum(1 for row in template_rows if row["status"] == "parse-error")
    complete_body_asset_count = sum(1 for row in body_asset_rows if row["status"] == "complete")
    engine_handoff_asset_count = sum(len(row["detectedAssets"]) for row in engine_asset_rows)
    required_field_count = sum(int(row["requiredFieldCount"]) for row in template_rows)
    filled_field_count = sum(int(row["filledFieldCount"]) for row in template_rows)
    incomplete_paths = [row["path"] for row in template_rows if row["status"] != "complete"]
    missing_body_assets = [
        row["assetId"] for row in body_asset_rows if row["status"] != "complete"
    ]
    intake_validation_ready = (
        complete_template_count == len(template_rows)
        and complete_body_asset_count == len(body_asset_rows)
        and bool(detected_radial_surface_assets)
    )
    can_attempt_exact_contact = intake_validation_ready and engine_handoff_asset_count > 0
    missing_evidence = []
    missing_evidence.extend(incomplete_paths)
    missing_evidence.extend(missing_body_assets)
    if not detected_radial_surface_assets:
        missing_evidence.append("radial_pad_hole_surfaces.step-or-filled-json")
    if engine_handoff_asset_count == 0:
        missing_evidence.append("external-engine-handoff-file")
    missing_evidence.append("external-engine-run-and-bench-comparison")
    return {
        "schema": TWO_CELL_SEGMENTED_CAD_INTAKE_VALIDATION_SCHEMA,
        "model": "two-cell-segmented-cad-intake-validation",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "intakeDirectory": intake_dir,
        "readinessSummary": readiness["summary"],
        "templateValidationRows": template_rows,
        "bodyAssetValidationRows": body_asset_rows,
        "radialSurfaceAssetValidation": {
            "status": "complete" if detected_radial_surface_assets else "missing",
            "expectedAssets": radial_surface_assets,
            "detectedAssets": detected_radial_surface_assets,
        },
        "engineHandoffValidationRows": engine_asset_rows,
        "summary": {
            "templateCount": len(template_rows),
            "existingTemplateCount": existing_template_count,
            "completeTemplateCount": complete_template_count,
            "missingTemplateCount": missing_template_count,
            "parseErrorCount": parse_error_count,
            "requiredFieldCount": required_field_count,
            "filledFieldCount": filled_field_count,
            "completeBodyAssetCount": complete_body_asset_count,
            "requiredBodyAssetCount": len(body_asset_rows),
            "radialSurfaceAssetDetected": bool(detected_radial_surface_assets),
            "engineHandoffAssetCount": engine_handoff_asset_count,
            "intakeValidationReady": intake_validation_ready,
            "canAttemptExactRigidBodyContact": can_attempt_exact_contact,
            "physicalAccuracyValidated": False,
            "missingEvidence": missing_evidence,
            "nextAction": (
                "Run an external rigid-body contact simulation and compare it against bench coordinates."
                if can_attempt_exact_contact
                else "Fill segmented CAD/body/contact templates and add at least one external engine handoff."
            ),
        },
        "claimLabels": {
            "oneCellCadSource": CAD_RAD_CELL_REFERENCE.source,
            "validationBoundary": "complete intake only permits an external contact run; it does not validate physical accuracy",
            "blockedClaim": "exact real-life mechanics until contact simulation and measured bench comparison pass",
        },
    }


def _two_cell_cad_instance_origins(layout: dict[str, Any]) -> dict[str, dict[str, Any]]:
    pitch_mm = float(layout["dimensionsMm"]["connectorPitchMm"])
    return {
        "left": {
            "cell": "left",
            "originMm": [0.0, 0.0, 0.0],
            "upperBody": "left_upper_free_cell_body",
            "lowerBody": "left_lower_cell_body",
            "pinBody": "left_screw_pin_body",
            "fixtureRule": "position-lockable fixture; left cell is fixed in the default two-cell bench",
        },
        "right": {
            "cell": "right",
            "originMm": [pitch_mm, 0.0, 0.0],
            "upperBody": "right_upper_free_cell_body",
            "lowerBody": "right_lower_cell_body",
            "pinBody": "right_screw_pin_body",
            "fixtureRule": "free, state-locked, or position-locked depending on the selected bench case",
        },
    }


def _shifted_site_center(site: dict[str, Any], origin_mm: list[float]) -> list[float]:
    return [
        float(origin_mm[0]) + float(site["xMm"]),
        float(origin_mm[1]) + float(site["yMm"]),
        float(origin_mm[2]),
    ]


def _hole_ring_segments(
    *,
    cell: str,
    site: dict[str, Any],
    ring_segments: int,
    origin_mm: list[float],
    upper_body: str,
    pin_body: str,
) -> list[dict[str, Any]]:
    center = _shifted_site_center(site, origin_mm)
    step_degrees = 360.0 / float(ring_segments)
    rows = []
    for index in range(ring_segments):
        rows.append(
            {
                "segmentIndex": index,
                "angleStartDeg": index * step_degrees,
                "angleEndDeg": (index + 1) * step_degrees,
                "centerMm": center,
                "axis": [0.0, 0.0, 1.0],
                "bodyA": pin_body,
                "bodyB": upper_body,
                "primitive": "convex_ring_sector_or_capsule_wall",
                "source": "layout-derived-placeholder-needs-CAD-surface",
                "physicalAccuracyValidated": False,
            }
        )
    return rows


def two_cell_cad_contact_decomposition_spec(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    ring_segments: int = 12,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
) -> dict[str, Any]:
    """CAD-to-contact contract for converting the one-cell A360 model into exact solver inputs."""

    config = replace(config or LatticeConfig(rows=1, cols=2), rows=1, cols=2)
    controls = _controls(controls)
    ring_segments = max(4, int(ring_segments))
    layout = cad_rad_cell_layout(config)
    dimensions = layout["dimensionsMm"]
    validation = two_cell_segmented_cad_intake_validation_report(config, controls, intake_dir=intake_dir)
    instances = _two_cell_cad_instance_origins(layout)
    site_by_name = {str(site["name"]): site for site in layout["padSites"]}
    clearance_mm = float(dimensions["pinHoleClearanceMm"])
    body_roles = [
        {
            "role": "upper_free_cell_body",
            "instances": [instances["left"]["upperBody"], instances["right"]["upperBody"]],
            "expectedAssets": list(SEGMENTED_CAD_EXPECTED_ASSETS["upper_free_cell_body"]),
            "solverUse": "moving rigid body carrying radial arms, annular pads, pin-hole surfaces, alpha rotation, and vertical free-play",
        },
        {
            "role": "lower_cell_body",
            "instances": [instances["left"]["lowerBody"], instances["right"]["lowerBody"]],
            "expectedAssets": list(SEGMENTED_CAD_EXPECTED_ASSETS["lower_cell_body"]),
            "solverUse": "fixture or base body defining the ground z reference and position-lock frame",
        },
        {
            "role": "screw_pin_body",
            "instances": [instances["left"]["pinBody"], instances["right"]["pinBody"]],
            "expectedAssets": list(SEGMENTED_CAD_EXPECTED_ASSETS["screw_pin_body"]),
            "solverUse": "pin or screw stack participating in radial clearance and vertical slip contacts",
        },
    ]
    hole_decomposition = []
    for cell, instance in instances.items():
        for site in layout["padSites"]:
            hole_decomposition.append(
                {
                    "cell": cell,
                    "site": site["name"],
                    "centerMm": _shifted_site_center(site, instance["originMm"]),
                    "axis": [0.0, 0.0, 1.0],
                    "holeWallBody": instance["upperBody"],
                    "pinBody": instance["pinBody"],
                    "holeRadiusMm": float(site["holeRadiusMm"]),
                    "pinRadiusMm": float(site["pinRadiusMm"]),
                    "clearanceMm": clearance_mm,
                    "ringSegmentCount": ring_segments,
                    "recommendedCollisionPrimitive": "convex_ring_sector_or_capsule_wall",
                    "ringSegments": _hole_ring_segments(
                        cell=cell,
                        site=site,
                        ring_segments=ring_segments,
                        origin_mm=instance["originMm"],
                        upper_body=instance["upperBody"],
                        pin_body=instance["pinBody"],
                    ),
                    "source": "A360 layout-derived site; replace with segmented radial_pad_hole_surfaces evidence before exact claims",
                    "physicalAccuracyValidated": False,
                }
            )
    connector_pairs = []
    marker_z = {
        "lower": 0.5 * float(dimensions["lowerBodyThicknessMm"]),
        "middle": float(dimensions["lowerBodyThicknessMm"]),
        "upper": float(dimensions["bodyThicknessMm"]),
    }
    for name, left_site_name, right_site_name in CAD_RAD_X_CONNECTOR_PAIRS:
        left_site = site_by_name[left_site_name]
        right_site = site_by_name[right_site_name]
        connector_pairs.append(
            {
                "name": name,
                "leftCell": "left",
                "rightCell": "right",
                "leftSite": left_site_name,
                "rightSite": right_site_name,
                "leftCenterMm": _shifted_site_center(left_site, instances["left"]["originMm"]),
                "rightCenterMm": _shifted_site_center(right_site, instances["right"]["originMm"]),
                "axis": [1.0, 0.0, 0.0],
                "leftPinBody": f"left_{left_site_name}_connector_pin_or_screw",
                "rightHoleBody": instances["right"]["upperBody"],
                "holeRadiusMm": float(dimensions["holeRadiusMm"]),
                "pinRadiusMm": float(dimensions["pinRadiusMm"]),
                "clearanceMm": clearance_mm,
                "requiredSurfacePair": "two_cell_connector_pairs.json plus radial_pad_hole_surfaces.json",
                "observableMarkers": [
                    {
                        "marker": f"{name}_{marker}",
                        "zMm": z_mm,
                        "description": f"{marker} connector marker coordinate for bench and engine comparison",
                    }
                    for marker, z_mm in marker_z.items()
                ],
                "source": "layout-derived X-neighbor connector pair; replace with measured connector origins before exact claims",
                "physicalAccuracyValidated": False,
            }
        )
    missing_evidence = sorted(
        set(
            [
                *validation["summary"].get("missingEvidence", []),
                "convex-or-analytic-hole-wall-contact-primitives",
                "body-collision-exclusion-mask",
                "joint-limit-stop-surfaces",
                "external-engine-contact-run",
                "bench-coordinate-holdout-comparison",
            ]
        )
    )
    collision_primitive_count = len(hole_decomposition) * ring_segments + len(connector_pairs)
    can_build = bool(validation["summary"].get("canAttemptExactRigidBodyContact")) and not missing_evidence
    return {
        "schema": TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA,
        "model": "two-cell-cad-contact-decomposition-contract",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "layout": layout,
        "cellInstances": list(instances.values()),
        "bodyRoles": body_roles,
        "jointSpec": {
            "alphaRotation": {
                "axis": [0.0, 0.0, 1.0],
                "requiresEvidence": ["joint_axes.json.central_alpha_rotation_axis"],
                "stateVariable": "alpha",
                "thetaMapping": "theta_degrees = 70 * alpha - 60 in the reduced model; replace by measured crown/arm relation when available",
            },
            "verticalSlide": {
                "axis": [0.0, 0.0, 1.0],
                "requiresEvidence": ["joint_axes.json.vertical_slide_axis", "contact_parameters.json.verticalTravelMm"],
                "stateVariable": "z",
                "contactRole": "pin-hole radius difference permits vertical residual motion under actuation and gravity",
            },
            "lockStops": {
                "requiresEvidence": ["lock_crown_geometry.csv", "lock_crown_geometry.json"],
                "note": "Crown angle is not equal to lattice theta because the arm has nonzero width; use measured effective stop theta.",
            },
        },
        "holeDecomposition": hole_decomposition,
        "twoCellContactPairs": connector_pairs,
        "collisionPrimitiveContract": {
            "holeWallPrimitive": "split each annular hole wall into convex sectors or analytic capsule/cylinder contacts",
            "avoid": "single concave hole mesh as the active collision body",
            "ringSegmentCount": ring_segments,
            "minimumRecommendedRingSegments": 8,
            "collisionPrimitiveCount": collision_primitive_count,
        },
        "engineAdapters": [
            {
                "engine": "MuJoCo",
                "preferredFile": "assets/cad/segmented/mujoco/rad_two_cell.xml",
                "contactGuidance": "use named bodies, slide/hinge joints, convex contact geoms, and margin/gap to expose pin-hole free-play before contact",
                "requiredOutputs": ["qpos", "body_xpos", "contact pairs", "constraint forces", "final marker xyz"],
            },
            {
                "engine": "Gazebo/Ignition",
                "preferredFile": "assets/cad/segmented/gazebo/rad_two_cell.sdf",
                "contactGuidance": "use simplified convex collision meshes separate from visual meshes; log joint states and contact wrenches",
                "requiredOutputs": ["link poses", "joint states", "contact wrenches", "final marker xyz"],
            },
            {
                "engine": "Isaac Sim",
                "preferredFile": "assets/cad/segmented/isaac/rad_two_cell.usd",
                "contactGuidance": "use USD visual mesh plus convex decomposition collision meshes with PhysX contact reporting",
                "requiredOutputs": ["body poses", "joint states", "contact reports", "final marker xyz"],
            },
        ],
        "benchTruthRequirements": {
            "coordinateFrame": "left lower-cell frame, millimeters, z=0 at the lower body ground plane",
            "observables": [
                "left/right cell center xyz",
                "left/right alpha",
                "left/right theta",
                "upper/middle/lower connector marker xyz",
                "lock held/released/fell",
                "stable state index",
            ],
            "holdoutPolicy": "do not tune on every bench case; reserve some lock/actuation cases for validation",
        },
        "summary": {
            "status": "ready-for-contact-model-build" if can_build else "needs-cad-contact-decomposition-inputs",
            "cellCount": len(instances),
            "holeCount": len(hole_decomposition),
            "twoCellConnectorCount": len(connector_pairs),
            "ringSegments": ring_segments,
            "collisionPrimitiveCount": collision_primitive_count,
            "canBuildExactContactModel": can_build,
            "physicalAccuracyValidated": False,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimBoundary": {
            "allowedClaim": "solver-ready decomposition contract derived from the A360 one-cell reference and current two-cell topology",
            "blockedClaim": "fabrication-accurate contact mechanics until segmented CAD primitives and bench holdouts validate the model",
            "physicalAccuracyValidated": False,
        },
    }


def simulate_two_cell_bench(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    hole_radius: float | None = None,
    pin_radius: float | None = None,
    axial_stiffness: float = 8.0,
    vertical_stiffness: float = 5.0,
    hinge_stiffness: float = 2.5,
    lock_stiffness: float = 18.0,
    contact_stiffness: float = 24.0,
    **control_overrides: Any,
) -> dict[str, Any]:
    """Reduced two-cell contact proxy used before segmented CAD/contact simulation exists."""

    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls, **control_overrides)
    radii = _radii(config, hole_radius=hole_radius, pin_radius=pin_radius)
    zones = _dead_zones(config, radii)
    alpha_free = _dead_zone(controls.alpha_command, zones["alphaDeadZone"])
    z_free = _dead_zone(controls.z_command, zones["verticalDeadZone"])
    alpha_transmission = config.coupling_gain * math.exp(-2.25 * zones["alphaDeadZone"])
    z_transmission = config.z_coupling_gain * math.exp(-3.1 * zones["verticalDeadZone"])
    right_blocked = controls.right_locked or controls.right_position_locked
    right_alpha_delta = 0.0 if right_blocked else alpha_free
    right_z = 0.0 if right_blocked else z_free
    left_alpha_delta = 0.0 if controls.left_position_locked else right_alpha_delta * alpha_transmission
    left_z = 0.0 if controls.left_position_locked else right_z * z_transmission
    left_alpha = _clamp(config.initial_alpha + left_alpha_delta, config.alpha_min, config.alpha_max)
    right_alpha = _clamp(config.initial_alpha + right_alpha_delta, config.alpha_min, config.alpha_max)
    pitch = config.cell_size
    rest_pitch = _two_cell_rest_pitch(config, left_alpha, right_alpha)
    free_left_x = 0.0 if controls.left_position_locked else 0.5 * (pitch - rest_pitch)
    free_right_x = rest_pitch if controls.left_position_locked else 0.5 * (pitch + rest_pitch)
    left_x = 0.0 if controls.left_position_locked else pitch - rest_pitch if controls.right_position_locked else free_left_x
    right_x = pitch if controls.right_position_locked else free_right_x
    distance = math.hypot(right_x - left_x, right_z - left_z)
    axial_strain = distance / max(rest_pitch, 1e-9) - 1.0
    left_theta = float(alpha_to_theta(left_alpha))
    right_theta = float(alpha_to_theta(right_alpha))
    theta_delta_radians = math.radians(right_theta - left_theta)
    vertical_shear = right_z - left_z
    contact_penetration = max(0.0, abs(controls.z_command) - zones["verticalDeadZone"])
    alpha_unrealized = abs(alpha_free - right_alpha_delta)
    z_unrealized = abs(z_free - right_z)
    position_unrealized = math.hypot(free_right_x - pitch, z_free) if controls.right_position_locked else 0.0
    locked_unrealized = math.hypot(alpha_free, z_free) if controls.right_locked else 0.0
    axial_energy = 0.5 * axial_stiffness * axial_strain**2
    hinge_energy = 0.5 * hinge_stiffness * theta_delta_radians**2
    vertical_energy = 0.5 * vertical_stiffness * vertical_shear**2
    contact_energy = 0.5 * contact_stiffness * contact_penetration**2
    lock_penalty_energy = 0.5 * lock_stiffness * (
        alpha_unrealized**2 + z_unrealized**2 + position_unrealized**2 + locked_unrealized**2
    )
    if right_blocked:
        contact_mode = "blocked-by-lock"
    elif abs(controls.alpha_command) <= zones["alphaDeadZone"] and abs(controls.z_command) <= zones["verticalDeadZone"]:
        contact_mode = "inside-backlash"
    elif contact_penetration > 0.0:
        contact_mode = "clearance-taken-up"
    else:
        contact_mode = "alpha-contact-only"
    return {
        "schema": TWO_CELL_BENCH_SCHEMA,
        "model": "cad-derived-two-cell-rigid-contact-proxy",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "units": "normalized-cell-units",
        "controls": {
            "alphaCommand": controls.alpha_command,
            "zCommand": controls.z_command,
            "holeSweepMax": controls.hole_sweep_max,
            "holeSweepSteps": controls.hole_sweep_steps,
            "leftPositionLocked": controls.left_position_locked,
            "rightLocked": controls.right_locked,
            "rightPositionLocked": controls.right_position_locked,
        },
        "dimensions": {"cellSize": pitch, **radii},
        "deadZones": zones,
        "cells": [
            {
                "id": "left",
                "positionLocked": controls.left_position_locked,
                "locked": False,
                "alpha": left_alpha,
                "theta": left_theta,
                "center": {"x": left_x, "y": 0.0, "z": left_z},
                "residualAlpha": left_alpha_delta,
                "residualZ": left_z,
            },
            {
                "id": "right",
                "positionLocked": controls.right_position_locked,
                "locked": controls.right_locked,
                "alpha": right_alpha,
                "theta": right_theta,
                "center": {"x": right_x, "y": 0.0, "z": right_z},
                "residualAlpha": right_alpha_delta,
                "residualZ": right_z,
            },
        ],
        "connector": {
            "pitch": pitch,
            "restPitch": rest_pitch,
            "distance": distance,
            "axialStrain": axial_strain,
            "verticalShear": vertical_shear,
            "contactMode": contact_mode,
            "contactPenetration": contact_penetration,
            "alphaTransmission": alpha_transmission,
            "zTransmission": z_transmission,
        },
        "energy": {
            "axialEnergy": axial_energy,
            "hingeEnergy": hinge_energy,
            "verticalEnergy": vertical_energy,
            "contactEnergy": contact_energy,
            "lockPenaltyEnergy": lock_penalty_energy,
            "totalEnergy": axial_energy
            + hinge_energy
            + vertical_energy
            + contact_energy
            + lock_penalty_energy,
        },
        "lockChecks": {
            "leftFixtureError": math.hypot(left_x, left_z) if controls.left_position_locked else None,
            "rightFixtureError": math.hypot(right_x - pitch, right_z) if controls.right_position_locked else None,
            "rightStateMotionSuppressed": right_blocked and abs(right_alpha_delta) < 1e-12 and abs(right_z) < 1e-12,
        },
        "interpretation": {
            "calibratedAccuracy": "not-yet-calibrated",
            "limitation": (
                "Reduced proxy: no segmented moving CAD bodies, no solved collisions, "
                "no friction, and no measured joint stiffness."
            ),
            "nextEvidenceNeeded": [
                "two-cell CAD assembly or STEP with separated moving bodies",
                "pin-hole contact radius and friction measurements",
                "alpha/z actuation traces for free, state-locked, and position-locked pairs",
            ],
        },
    }


def sweep_two_cell_backlash(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    hole_sweep_max: float | None = None,
    hole_sweep_steps: int | None = None,
    **overrides: Any,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(
        controls,
        hole_sweep_max=hole_sweep_max,
        hole_sweep_steps=hole_sweep_steps,
        **overrides,
    )
    radii = _radii(config)
    end = max(radii["pinRadius"], controls.hole_sweep_max)
    rows = []
    for index in range(controls.hole_sweep_steps):
        t = index / max(1, controls.hole_sweep_steps - 1)
        hole_radius = radii["pinRadius"] + (end - radii["pinRadius"]) * t
        result = simulate_two_cell_bench(config, controls, hole_radius=hole_radius)
        rows.append(
            {
                "index": index,
                "holeRadius": hole_radius,
                "clearance": result["dimensions"]["pinHoleClearance"],
                "effectiveBacklash": result["deadZones"]["effectiveBacklash"],
                "alphaDeadZone": result["deadZones"]["alphaDeadZone"],
                "verticalDeadZone": result["deadZones"]["verticalDeadZone"],
                "leftResidualZ": result["cells"][0]["residualZ"],
                "rightZ": result["cells"][1]["center"]["z"],
                "axialStrain": result["connector"]["axialStrain"],
                "contactMode": result["connector"]["contactMode"],
                "totalEnergy": result["energy"]["totalEnergy"],
            }
        )
    first = rows[0]
    last = rows[-1]
    return {
        "schema": TWO_CELL_SWEEP_SCHEMA,
        "model": "cad-derived-two-cell-clearance-sweep",
        "controls": simulate_two_cell_bench(config, controls)["controls"],
        "pinRadius": radii["pinRadius"],
        "holeRadiusStart": first["holeRadius"],
        "holeRadiusEnd": last["holeRadius"],
        "rows": rows,
        "trend": {
            "neighborZStart": first["leftResidualZ"],
            "neighborZEnd": last["leftResidualZ"],
            "rightZStart": first["rightZ"],
            "rightZEnd": last["rightZ"],
            "energyStart": first["totalEnergy"],
            "energyEnd": last["totalEnergy"],
            "clearanceIncreasesBacklash": last["effectiveBacklash"] >= first["effectiveBacklash"],
            "neighborResponseDropsWithClearance": abs(last["leftResidualZ"]) <= abs(first["leftResidualZ"]) + 1e-12,
        },
    }


def _lock_mode_controls(controls: TwoCellBenchControls, mode: str) -> TwoCellBenchControls:
    if mode == "free":
        return _controls(controls, left_position_locked=True, right_locked=False, right_position_locked=False)
    if mode == "right_state_locked":
        return _controls(controls, right_locked=True, right_position_locked=False)
    if mode == "right_position_locked":
        return _controls(controls, right_locked=False, right_position_locked=True)
    if mode == "left_free":
        return _controls(controls, left_position_locked=False, right_locked=False, right_position_locked=False)
    raise ValueError(f"Unknown two-cell lock mode: {mode}")


def _lock_mode_from_controls(controls: TwoCellBenchControls) -> str:
    if controls.right_position_locked:
        return "right_position_locked"
    if controls.right_locked:
        return "right_state_locked"
    if not controls.left_position_locked:
        return "left_free"
    return "free"


def _named_case_controls(controls: TwoCellBenchControls) -> list[tuple[str, TwoCellBenchControls]]:
    return [
        ("free_contract_lift", controls),
        ("free_expand_pushdown", _controls(controls, alpha_command=0.35, z_command=-0.35)),
        ("right_state_locked", _controls(controls, right_locked=True, right_position_locked=False)),
        ("right_position_locked", _controls(controls, right_locked=False, right_position_locked=True)),
        ("left_free_neighbor_residual", _controls(controls, left_position_locked=False)),
    ]


def _unique_floats(values: list[float] | tuple[float, ...]) -> list[float]:
    ordered: list[float] = []
    for value in values:
        numeric = float(value)
        if not any(math.isclose(numeric, existing, rel_tol=0.0, abs_tol=1e-12) for existing in ordered):
            ordered.append(numeric)
    return ordered


def two_cell_actuation_sweep(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    alpha_commands: list[float] | tuple[float, ...] | None = None,
    z_commands: list[float] | tuple[float, ...] | None = None,
    hole_radii: list[float] | tuple[float, ...] | None = None,
    lock_modes: list[str] | tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """Sweep two-cell actuation and lock cases for physical calibration planning."""

    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    radii = _radii(config)
    alpha_values = _unique_floats(list(alpha_commands or (-0.5, -0.25, 0.0, 0.25, 0.5)))
    z_values = _unique_floats(list(z_commands or (-0.45, 0.0, 0.45)))
    hole_values = _unique_floats(
        list(hole_radii or (radii["pinRadius"], radii["holeRadius"], controls.hole_sweep_max))
    )
    modes = list(lock_modes or ("free", "right_state_locked", "right_position_locked", "left_free"))
    rows = []
    index = 0
    for mode in modes:
        mode_base = _lock_mode_controls(controls, mode)
        for hole_radius in hole_values:
            for alpha_command in alpha_values:
                for z_command in z_values:
                    case_controls = _controls(
                        mode_base,
                        alpha_command=alpha_command,
                        z_command=z_command,
                    )
                    result = simulate_two_cell_bench(config, case_controls, hole_radius=hole_radius)
                    left, right = result["cells"]
                    row = {
                        "index": index,
                        "lockMode": mode,
                        "alphaCommand": alpha_command,
                        "zCommand": z_command,
                        "pinRadius": result["dimensions"]["pinRadius"],
                        "holeRadius": result["dimensions"]["holeRadius"],
                        "clearance": result["dimensions"]["pinHoleClearance"],
                        "effectiveBacklash": result["deadZones"]["effectiveBacklash"],
                        "contactMode": result["connector"]["contactMode"],
                        "leftX": left["center"]["x"],
                        "leftY": left["center"]["y"],
                        "leftZ": left["center"]["z"],
                        "rightX": right["center"]["x"],
                        "rightY": right["center"]["y"],
                        "rightZ": right["center"]["z"],
                        "leftAlpha": left["alpha"],
                        "rightAlpha": right["alpha"],
                        "leftTheta": left["theta"],
                        "rightTheta": right["theta"],
                        "verticalShear": result["connector"]["verticalShear"],
                        "axialStrain": result["connector"]["axialStrain"],
                        "contactPenetration": result["connector"]["contactPenetration"],
                        "totalEnergy": result["energy"]["totalEnergy"],
                        "leftFixtureError": result["lockChecks"]["leftFixtureError"],
                        "rightFixtureError": result["lockChecks"]["rightFixtureError"],
                        "rightStateMotionSuppressed": result["lockChecks"]["rightStateMotionSuppressed"],
                    }
                    rows.append(row)
                    index += 1
    contact_counts: dict[str, int] = {}
    for row in rows:
        contact_counts[row["contactMode"]] = contact_counts.get(row["contactMode"], 0) + 1
    clearances = [row["clearance"] for row in rows]
    right_z_values = [row["rightZ"] for row in rows]
    return {
        "schema": TWO_CELL_ACTUATION_SWEEP_SCHEMA,
        "model": "cad-derived-two-cell-actuation-lock-clearance-sweep",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "controls": simulate_two_cell_bench(config, controls)["controls"],
        "grid": {
            "alphaCommands": alpha_values,
            "zCommands": z_values,
            "holeRadii": hole_values,
            "lockModes": modes,
        },
        "rows": rows,
        "summary": {
            "rowCount": len(rows),
            "lockModes": modes,
            "clearanceMin": min(clearances) if clearances else 0.0,
            "clearanceMax": max(clearances) if clearances else 0.0,
            "rightZMin": min(right_z_values) if right_z_values else 0.0,
            "rightZMax": max(right_z_values) if right_z_values else 0.0,
            "maxAbsAxialStrain": max((abs(row["axialStrain"]) for row in rows), default=0.0),
            "maxTotalEnergy": max((row["totalEnergy"] for row in rows), default=0.0),
            "contactModeCounts": contact_counts,
            "requiresMeasuredData": True,
            "requiresSegmentedPhysics": True,
        },
    }


def _xml_attr(value: object) -> str:
    return escape(str(value), {'"': "&quot;"})


def _fmt(value: float) -> str:
    return f"{float(value):.9g}"


def _fmt_vec(values: list[float] | tuple[float, ...]) -> str:
    return " ".join(_fmt(value) for value in values)


def _two_cell_scale_mm(config: LatticeConfig, bench: dict[str, Any]) -> dict[str, float]:
    cad = CAD_RAD_CELL_REFERENCE
    reference_hole_radius = max(float(config.hole_radius), 1e-9)
    nominal_hole_radius = 0.5 * cad.nominal_hole_diameter_mm
    hole_radius = float(bench["dimensions"]["holeRadius"])
    pin_radius = float(bench["dimensions"]["pinRadius"])
    pin_radius_mm = nominal_hole_radius * pin_radius / reference_hole_radius
    hole_radius_mm = nominal_hole_radius * hole_radius / reference_hole_radius
    layout_dimensions = _cad_rad_layout_dimensions(
        width_mm=min(cad.width_x_mm, cad.width_y_mm),
        height_mm=cad.height_z_mm,
        pin_radius_mm=pin_radius_mm,
        hole_radius_mm=hole_radius_mm,
    )
    connector_pitch_mm = layout_dimensions["connectorPitchMm"]
    return {
        "xyMmPerModelUnit": connector_pitch_mm / max(float(config.cell_size), 1e-9),
        "zMmPerModelUnit": cad.height_z_mm,
        "cadWidthMm": cad.width_x_mm,
        "cadHeightMm": cad.height_z_mm,
        "connectorPitchMm": connector_pitch_mm,
        "bodyThicknessMm": 4.0,
        "nominalHoleRadiusMm": nominal_hole_radius,
        "pinRadiusMm": pin_radius_mm,
        "holeRadiusMm": hole_radius_mm,
        "clearanceMm": max(0.0, hole_radius_mm - pin_radius_mm),
    }


def _scale_center_mm(center: dict[str, float], scale: dict[str, float]) -> list[float]:
    return [
        float(center["x"]) * scale["xyMmPerModelUnit"],
        float(center["y"]) * scale["xyMmPerModelUnit"],
        float(center["z"]) * scale["zMmPerModelUnit"],
    ]


def _mjcf_site_layout(width_mm: float, hole_radius_mm: float, pin_radius_mm: float) -> list[tuple[str, float, float]]:
    dimensions = _cad_rad_layout_dimensions(
        width_mm=width_mm,
        height_mm=CAD_RAD_CELL_REFERENCE.height_z_mm,
        hole_radius_mm=hole_radius_mm,
        pin_radius_mm=pin_radius_mm,
    )
    return [
        (str(site["name"]), float(site["xMm"]), float(site["yMm"]))
        for site in _cad_rad_site_records(dimensions)
    ]


def _mjcf_body_block(
    *,
    body: dict[str, Any],
    cell: dict[str, Any],
    scale: dict[str, float],
    layout: list[tuple[str, float, float]],
    fixed: bool,
    actuated: bool,
) -> tuple[list[str], list[dict[str, Any]], list[dict[str, Any]]]:
    name = str(body["name"])
    rgba = str(body["rgba"])
    center = body["positionMm"]
    half_height = 0.5 * scale["bodyThicknessMm"]
    width = scale["cadWidthMm"]
    hub_radius = max(scale["holeRadiusMm"] * 1.6, 0.065 * width)
    pad_radius = max(scale["holeRadiusMm"] * 1.85, 0.055 * width)
    arm_radius = max(0.85, 0.018 * width)
    hole_radius = scale["holeRadiusMm"]
    pin_radius = scale["pinRadiusMm"]
    screw_radius = max(0.65 * pin_radius, 0.8)
    screw_half_height = max(0.5, 0.25 * scale["cadHeightMm"])
    lines = [f'    <body name="{_xml_attr(name)}" pos="{_fmt_vec(center)}">']
    joint_records: list[dict[str, Any]] = []
    if not fixed:
        joint_specs = [
            ("slide_x", "slide", (1.0, 0.0, 0.0), True),
            ("slide_z", "slide", (0.0, 0.0, 1.0), True),
            ("yaw", "hinge", (0.0, 0.0, 1.0), actuated),
        ]
        for suffix, joint_type, axis, actuator_enabled in joint_specs:
            joint_name = f"{name}_{suffix}"
            lines.append(
                f'      <joint name="{_xml_attr(joint_name)}" type="{joint_type}" axis="{_fmt_vec(axis)}" damping="0.35"/>'
            )
            joint_records.append(
                {
                    "body": name,
                    "joint": joint_name,
                    "type": joint_type,
                    "axis": list(axis),
                    "actuatorEnabled": bool(actuator_enabled),
                }
            )
    lines.append(
        f'      <geom name="{_xml_attr(name)}_hub" type="cylinder" pos="0 0 0" size="{_fmt(hub_radius)} {_fmt(half_height)}" rgba="{rgba}"/>'
    )
    lines.append(
        f'      <geom name="{_xml_attr(name)}_screw" type="cylinder" pos="0 0 {_fmt(half_height + screw_half_height)}" size="{_fmt(screw_radius)} {_fmt(screw_half_height)}" rgba="0.16 0.16 0.15 1"/>'
    )
    geom_records: list[dict[str, Any]] = []
    for site, x, y in layout:
        arm = f"{name}_{site}_arm"
        pad = f"{name}_{site}_pad"
        pin = f"{name}_{site}_pin"
        hole = f"{name}_{site}_hole_clearance"
        marker = f"{name}_{site}_marker"
        lines.append(
            f'      <geom name="{_xml_attr(arm)}" type="capsule" fromto="0 0 0 {_fmt(x)} {_fmt(y)} 0" size="{_fmt(arm_radius)}" rgba="{rgba}"/>'
        )
        lines.append(
            f'      <geom name="{_xml_attr(pad)}" type="cylinder" pos="{_fmt(x)} {_fmt(y)} 0" size="{_fmt(pad_radius)} {_fmt(half_height)}" rgba="{rgba}"/>'
        )
        lines.append(
            f'      <geom name="{_xml_attr(pin)}" type="cylinder" pos="{_fmt(x)} {_fmt(y)} 0" size="{_fmt(pin_radius)} {_fmt(half_height * 1.18)}" rgba="0.05 0.05 0.05 1"/>'
        )
        lines.append(
            f'      <geom name="{_xml_attr(hole)}" type="cylinder" pos="{_fmt(x)} {_fmt(y)} 0" size="{_fmt(hole_radius)} {_fmt(half_height * 1.22)}" rgba="0.1 0.65 0.9 0.22"/>'
        )
        lines.append(
            f'      <site name="{_xml_attr(marker)}" pos="{_fmt(x)} {_fmt(y)} 0" size="{_fmt(max(0.35, pin_radius * 0.35))}" rgba="1 0.92 0.15 1"/>'
        )
        geom_records.append(
            {
                "body": name,
                "site": site,
                "arm": arm,
                "pad": pad,
                "pin": pin,
                "hole": hole,
                "marker": marker,
                "localPositionMm": [x, y, 0.0],
            }
        )
    lines.append("    </body>")
    return lines, joint_records, geom_records


def two_cell_mjcf_proxy_report(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    hole_radius: float | None = None,
    pin_radius: float | None = None,
    lock_mode: str | None = None,
    timestep: float = 0.001,
    gravity: tuple[float, float, float] = (0.0, 0.0, -9.81),
    **control_overrides: Any,
) -> dict[str, Any]:
    """CAD-dimensioned two-cell MuJoCo XML proxy for external physics handoff."""

    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls, **control_overrides)
    if lock_mode is not None:
        controls = _lock_mode_controls(controls, lock_mode)
    bench = simulate_two_cell_bench(
        config,
        controls,
        hole_radius=hole_radius,
        pin_radius=pin_radius,
    )
    scale = _two_cell_scale_mm(config, bench)
    left_cell, right_cell = bench["cells"]
    bodies = [
        {
            "id": "left",
            "name": "left_cell",
            "positionMm": _scale_center_mm(left_cell["center"], scale),
            "positionLocked": bool(left_cell["positionLocked"]),
            "stateLocked": bool(left_cell["locked"]),
            "alpha": left_cell["alpha"],
            "theta": left_cell["theta"],
            "rgba": "0.12 0.28 0.42 1",
        },
        {
            "id": "right",
            "name": "right_cell",
            "positionMm": _scale_center_mm(right_cell["center"], scale),
            "positionLocked": bool(right_cell["positionLocked"]),
            "stateLocked": bool(right_cell["locked"]),
            "alpha": right_cell["alpha"],
            "theta": right_cell["theta"],
            "rgba": "0.70 0.34 0.12 1",
        },
    ]
    for body in bodies:
        body["fixed"] = bool(body["positionLocked"] or body["stateLocked"])
    cad_layout = cad_rad_cell_layout(
        config,
        hole_radius=bench["dimensions"]["holeRadius"],
        pin_radius=bench["dimensions"]["pinRadius"],
    )
    layout = _mjcf_site_layout(scale["cadWidthMm"], scale["holeRadiusMm"], scale["pinRadiusMm"])
    body_lines: list[str] = []
    joints: list[dict[str, Any]] = []
    geoms: list[dict[str, Any]] = []
    for body, cell in zip(bodies, (left_cell, right_cell)):
        lines, body_joints, body_geoms = _mjcf_body_block(
            body=body,
            cell=cell,
            scale=scale,
            layout=layout,
            fixed=bool(body["fixed"]),
            actuated=body["id"] == "right",
        )
        body_lines.extend(lines)
        joints.extend(body_joints)
        geoms.extend(body_geoms)
    connector_pairs = [
        ("upper", "left_cell_ne_pin", "right_cell_nw_hole_clearance"),
        ("middle", "left_cell_e_pin", "right_cell_w_hole_clearance"),
        ("lower", "left_cell_se_pin", "right_cell_sw_hole_clearance"),
    ]
    contact_gap = max(0.0, 0.25 * scale["clearanceMm"])
    contact_margin = max(0.0, scale["clearanceMm"])
    contact_records = [
        {
            "connector": connector,
            "pinGeom": pin_geom,
            "holeGeom": hole_geom,
            "clearanceMm": scale["clearanceMm"],
            "margin": contact_margin,
            "gap": contact_gap,
        }
        for connector, pin_geom, hole_geom in connector_pairs
    ]
    contact_lines = [
        f'    <pair name="{_xml_attr(record["connector"])}_pin_hole" geom1="{_xml_attr(record["pinGeom"])}" geom2="{_xml_attr(record["holeGeom"])}" margin="{_fmt(record["margin"])}" gap="{_fmt(record["gap"])}" condim="3" friction="0.4 0.02 0.001"/>'
        for record in contact_records
    ]
    actuator_records = [
        {
            "name": "right_cell_alpha_slide",
            "joint": "right_cell_slide_x",
            "channel": "alpha",
            "command": controls.alpha_command,
            "targetMm": bench["cells"][1]["center"]["x"] * scale["xyMmPerModelUnit"],
        },
        {
            "name": "right_cell_vertical_slide",
            "joint": "right_cell_slide_z",
            "channel": "z",
            "command": controls.z_command,
            "targetMm": bench["cells"][1]["center"]["z"] * scale["zMmPerModelUnit"],
        },
        {
            "name": "right_cell_yaw_theta",
            "joint": "right_cell_yaw",
            "channel": "theta",
            "command": bench["cells"][1]["theta"],
            "targetMm": 0.0,
        },
    ]
    free_joint_names = {record["joint"] for record in joints}
    actuator_records = [
        record for record in actuator_records if record["joint"] in free_joint_names
    ]
    actuator_lines = [
        f'    <motor name="{_xml_attr(record["name"])}" joint="{_xml_attr(record["joint"])}" gear="1" ctrlrange="-1 1"/>'
        for record in actuator_records
    ]
    xml_lines = [
        '<mujoco model="rad_two_cell_cad_proxy">',
        '  <compiler angle="degree" coordinate="local"/>',
        f'  <option timestep="{_fmt(max(float(timestep), 1e-6))}" gravity="{_fmt_vec(gravity)}"/>',
        '  <worldbody>',
        '    <geom name="ground" type="plane" pos="0 0 -2" size="90 45 1" rgba="0.82 0.84 0.82 1"/>',
        *body_lines,
        '  </worldbody>',
        '  <contact>',
        *contact_lines,
        '  </contact>',
    ]
    if actuator_lines:
        xml_lines.extend(['  <actuator>', *actuator_lines, '  </actuator>'])
    xml_lines.append("</mujoco>")
    missing_evidence = [
        "segmentedStepOrMeshExport",
        "measuredJointAxes",
        "measuredFriction",
        "measuredContactStiffness",
        "externalMuJoCoRun",
    ]
    lock_records = [
        {
            "body": body["name"],
            "positionLocked": body["positionLocked"],
            "stateLocked": body["stateLocked"],
            "fixedInMjcf": body["fixed"],
        }
        for body in bodies
    ]
    return {
        "schema": TWO_CELL_MJCF_SCHEMA,
        "model": "cad-dimensioned-two-cell-mujoco-proxy",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "cadLayout": cad_layout,
        "bench": bench,
        "scale": scale,
        "bodies": bodies,
        "joints": joints,
        "actuators": actuator_records,
        "geoms": geoms,
        "contacts": contact_records,
        "locks": lock_records,
        "xml": "\n".join(xml_lines) + "\n",
        "summary": {
            "status": "ready-for-external-mujoco-run-with-proxy-caveats",
            "bodyCount": len(bodies),
            "fixedBodyCount": sum(1 for body in bodies if body["fixed"]),
            "jointCount": len(joints),
            "actuatorCount": len(actuator_records),
            "padSiteCount": len(geoms),
            "connectorContactPairCount": len(contact_records),
            "pinRadiusMm": scale["pinRadiusMm"],
            "holeRadiusMm": scale["holeRadiusMm"],
            "clearanceMm": scale["clearanceMm"],
            "rightStateMotionSuppressed": bench["lockChecks"]["rightStateMotionSuppressed"],
            "externalRunComplete": False,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimLabels": {
            "geometry": "CAD-dimensioned proxy using A360 one-cell bounding box and manifest hole label",
            "dynamics": "not yet externally solved in this environment",
            "physicalAccuracy": "unvalidated until STEP/mesh segmentation and bench data agree",
        },
        "limitations": [
            "The .f3d archive is not parsed into separate rigid bodies here.",
            "Arms, pads, pins, and hole shells are generated proxy geoms, not exact BREP surfaces.",
            "Locking is represented by omitting free joints on fixed bodies; real lock crown compliance is not modeled.",
            "MuJoCo is not installed in this environment, so this report exports a runnable candidate but not a completed external run.",
        ],
    }


def _engine_available(engine_name: str, engine_availability: dict[str, bool] | None = None) -> bool:
    if engine_availability is not None and engine_name in engine_availability:
        return bool(engine_availability[engine_name])
    return importlib.util.find_spec(engine_name) is not None


def two_cell_mujoco_proxy_run_report(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    steps: int = 120,
    timestep: float = 0.001,
    engine_availability: dict[str, bool] | None = None,
    **proxy_overrides: Any,
) -> dict[str, Any]:
    """Run the two-cell CAD proxy with MuJoCo when the optional package exists."""

    if steps <= 0:
        raise ValueError("steps must be positive")
    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    proxy = two_cell_mjcf_proxy_report(
        config,
        controls,
        timestep=timestep,
        **proxy_overrides,
    )
    proxy_summary = proxy["summary"]
    expected_body_count = int(proxy_summary.get("bodyCount", 0) or 0)
    available = _engine_available("mujoco", engine_availability)
    missing_evidence: list[str] = []
    if expected_body_count <= 0:
        missing_evidence.append("twoCellMjcfBodies")
    if "<mujoco" not in str(proxy.get("xml", "")):
        missing_evidence.append("twoCellMjcfXml")
    if not available:
        missing_evidence.extend(["mujocoPythonPackage", "externalMuJoCoRun"])

    solver = {
        "engineAvailable": available,
        "mujocoAvailable": available,
        "ran": False,
        "stepsRequested": int(steps),
        "stepsCompleted": 0,
        "timestep": float(timestep),
        "engineVersion": "",
        "error": "",
    }
    results: list[dict[str, Any]] = []
    qpos: list[float] = []
    qvel: list[float] = []
    ctrl: list[float] = []

    if not missing_evidence:
        try:
            mujoco = importlib.import_module("mujoco")
            solver["engineVersion"] = str(getattr(mujoco, "__version__", ""))
            model = mujoco.MjModel.from_xml_string(str(proxy["xml"]))
            data = mujoco.MjData(model)

            actuator_lookup: dict[str, int] = {}
            for actuator in proxy["actuators"]:
                actuator_id = int(
                    mujoco.mj_name2id(
                        model,
                        mujoco.mjtObj.mjOBJ_ACTUATOR,
                        str(actuator["name"]),
                    )
                )
                if actuator_id >= 0:
                    actuator_lookup[str(actuator["name"])] = actuator_id
                    data.ctrl[actuator_id] = _clamp(float(actuator["command"]), -1.0, 1.0)
            ctrl = [float(value) for value in data.ctrl.tolist()]

            body_ids: dict[str, int] = {}
            for body in proxy["bodies"]:
                body_id = int(
                    mujoco.mj_name2id(
                        model,
                        mujoco.mjtObj.mjOBJ_BODY,
                        str(body["name"]),
                    )
                )
                if body_id >= 0:
                    body_ids[str(body["name"])] = body_id
            mujoco.mj_forward(model, data)
            for _ in range(int(steps)):
                mujoco.mj_step(model, data)
                solver["stepsCompleted"] = int(solver["stepsCompleted"]) + 1

            for body in proxy["bodies"]:
                name = str(body["name"])
                body_id = body_ids.get(name)
                if body_id is None:
                    continue
                initial = [float(value) for value in body["positionMm"]]
                final = [float(value) for value in data.xpos[body_id, 0:3].tolist()]
                displacement = [final[index] - initial[index] for index in range(3)]
                results.append(
                    {
                        "bodyId": body["id"],
                        "bodyName": name,
                        "initialPositionMm": initial,
                        "finalPositionMm": final,
                        "displacementMm": displacement,
                        "displacementNormMm": float(math.sqrt(sum(value * value for value in displacement))),
                        "positionLocked": bool(body["positionLocked"]),
                        "stateLocked": bool(body["stateLocked"]),
                        "fixedInMjcf": bool(body["fixed"]),
                    }
                )
            qpos = [float(value) for value in data.qpos.tolist()]
            qvel = [float(value) for value in data.qvel.tolist()]
        except Exception as exc:  # pragma: no cover - depends on optional external engine
            solver["error"] = str(exc)
            missing_evidence.extend(["mujocoExecution", "externalMuJoCoRun"])

    if len(results) < expected_body_count:
        missing_evidence.append("mujocoBodyResults")
    missing_evidence = sorted(set(missing_evidence))
    run_complete = (
        available
        and not solver["error"]
        and int(solver["stepsCompleted"]) == int(steps)
        and len(results) >= expected_body_count
        and not missing_evidence
    )
    solver["ran"] = bool(run_complete)
    fixed_displacements = [
        float(result["displacementNormMm"])
        for result in results
        if bool(result["fixedInMjcf"])
    ]
    max_fixture_displacement = max(fixed_displacements, default=0.0)
    return {
        "schema": TWO_CELL_MUJOCO_RUN_SCHEMA,
        "method": "optional MuJoCo execution of the A360-dimensioned two-cell RAD proxy",
        "engine": "MuJoCo",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "mjcfProxy": {
            "schema": proxy["schema"],
            "model": proxy["model"],
            "summary": proxy_summary,
            "scale": proxy["scale"],
            "bodies": proxy["bodies"],
            "actuators": proxy["actuators"],
            "contacts": proxy["contacts"],
            "sourceXmlLength": len(str(proxy["xml"])),
            "sourceXmlArtifact": "outputs/two_cell_bench_packet/two_cell_cad_proxy.xml",
        },
        "solver": solver,
        "results": {
            "bodies": results,
            "qpos": qpos,
            "qvel": qvel,
            "ctrl": ctrl,
        },
        "summary": {
            "status": "two-cell-mujoco-proxy-run-complete" if run_complete else "needs-two-cell-mujoco-run",
            "twoCellMujocoProxyRunComplete": run_complete,
            "mujocoProxyRunComplete": run_complete,
            "mujocoAvailable": available,
            "bodyResultCount": len(results),
            "expectedBodyResultCount": expected_body_count,
            "fixedBodyResultCount": sum(1 for result in results if bool(result["fixedInMjcf"])),
            "maxFixtureDisplacementMm": max_fixture_displacement,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "claimLabels": {
            "cadGeometry": "A360 one-cell dimensions drive the proxy scale, not exact segmented BREP contact",
            "run": "external-engine-derived result" if run_complete else "external-engine run not completed",
            "physicalAccuracy": "not established until segmented CAD/contact setup and bench coordinates agree",
        },
        "limitations": [
            "The MJCF model is a reduced proxy, not a direct import of separated Fusion rigid bodies.",
            "When MuJoCo is unavailable this report is an explicit missing-evidence artifact.",
            "A completed proxy run still needs measured friction, contact stiffness, joint axes, and bench validation.",
        ],
    }


def _external_fidelity_mjcf_index_default_path() -> Path:
    return Path("outputs/two_cell_bench_packet/two_cell_external_fidelity_matrix_mjcf_index.json")


def _resolve_index_relative_path(index_path: Path, candidate: str | Path) -> Path:
    path = Path(candidate)
    if path.is_absolute():
        return path
    return index_path.parent / path


def _site_world_position(mujoco: Any, model: Any, data: Any, site_name: str) -> list[float] | None:
    site_id = int(mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, site_name))
    if site_id < 0:
        return None
    return [float(value) for value in data.site_xpos[site_id, 0:3].tolist()]


def _connector_mode_from_slip(lateral_slip: float, vertical_slip: float, clearance: float) -> str:
    lateral_excess = max(0.0, lateral_slip - clearance)
    vertical_excess = max(0.0, vertical_slip - clearance)
    if lateral_excess > 1e-9 and vertical_excess > 1e-9:
        return "lateral+vertical-contact"
    if lateral_excess > 1e-9:
        return "lateral-contact"
    if vertical_excess > 1e-9:
        return "vertical-contact"
    return "inside-clearance"


def _body_result_for(body_results: list[dict[str, Any]], body_id: str) -> dict[str, Any] | None:
    for result in body_results:
        if str(result.get("bodyId")) == body_id:
            return result
    return None


def _normalized_body_position(body_result: dict[str, Any] | None, case: dict[str, Any]) -> dict[str, float | str]:
    if body_result is None:
        return {"x": "", "y": "", "z": ""}
    position = list(body_result.get("finalPositionMm", []))
    if len(position) < 3:
        return {"x": "", "y": "", "z": ""}
    xy_scale = float(case.get("xyMmPerModelUnit", 0.0) or 0.0)
    z_scale = float(case.get("zMmPerModelUnit", 0.0) or 0.0)
    return {
        "x": position[0] / xy_scale if xy_scale > 0.0 else "",
        "y": position[1] / xy_scale if xy_scale > 0.0 else "",
        "z": position[2] / z_scale if z_scale > 0.0 else "",
    }


def _lock_held_observation(case: dict[str, Any], body_results: list[dict[str, Any]]) -> str:
    lock_mode = str(case.get("lockMode", "")).strip()
    if lock_mode not in {"right_state_locked", "right_position_locked"}:
        return ""
    right = _body_result_for(body_results, "right")
    if right is None:
        return ""
    displacement = float(right.get("displacementNormMm", math.inf))
    return "held" if displacement <= 1e-6 else "failed"


def _case_connector_measurement_rows(case_result: dict[str, Any]) -> list[dict[str, Any]]:
    case = dict(case_result.get("case", {}))
    body_results = list(case_result.get("bodies", []))
    left_cell = _normalized_body_position(_body_result_for(body_results, "left"), case)
    right_cell = _normalized_body_position(_body_result_for(body_results, "right"), case)
    marker_positions = dict(case_result.get("markerPositionsMm", {}))
    lock_observation = _lock_held_observation(case, body_results)
    clearance = float(case.get("pinHoleClearanceMm", 0.0) or 0.0)
    rows: list[dict[str, Any]] = []
    for connector, left_site, right_site in CAD_RAD_X_CONNECTOR_PAIRS:
        left_position = marker_positions.get(f"left_cell_{left_site}_marker")
        right_position = marker_positions.get(f"right_cell_{right_site}_marker")
        if left_position is None or right_position is None:
            continue
        slip_x = float(right_position[0]) - float(left_position[0])
        slip_y = float(right_position[1]) - float(left_position[1])
        slip_z = float(right_position[2]) - float(left_position[2])
        lateral_slip = math.hypot(slip_x, slip_y)
        vertical_slip = abs(slip_z)
        total_slip = math.sqrt(lateral_slip**2 + vertical_slip**2)
        rows.append(
            {
                "caseId": case_result.get("caseId", ""),
                "matrixIndex": case.get("matrixIndex", ""),
                "connector": connector,
                "leftSite": left_site,
                "rightSite": right_site,
                "actuationCase": case.get("actuationCase", ""),
                "lockMode": case.get("lockMode", ""),
                "backlash": case.get("backlash", ""),
                "alphaCommand": case.get("alphaCommand", ""),
                "zCommand": case.get("zCommand", ""),
                "gravityForce": case.get("gravityForce", ""),
                "pinRadius": case.get("pinRadius", ""),
                "holeRadius": case.get("holeRadius", ""),
                "pinHoleClearanceMm": case.get("pinHoleClearanceMm", ""),
                "observedLeftCellX": left_cell["x"],
                "observedLeftCellY": left_cell["y"],
                "observedLeftCellZ": left_cell["z"],
                "observedRightCellX": right_cell["x"],
                "observedRightCellY": right_cell["y"],
                "observedRightCellZ": right_cell["z"],
                "observedLeftAlpha": "",
                "observedRightAlpha": "",
                "observedLeftTheta": "",
                "observedRightTheta": "",
                "observedLeftXmm": float(left_position[0]),
                "observedLeftYmm": float(left_position[1]),
                "observedLeftZmm": float(left_position[2]),
                "observedRightXmm": float(right_position[0]),
                "observedRightYmm": float(right_position[1]),
                "observedRightZmm": float(right_position[2]),
                "observedLateralSlipMm": lateral_slip,
                "observedVerticalSlipMm": vertical_slip,
                "observedTotalSlipMm": total_slip,
                "observedContactMode": _connector_mode_from_slip(lateral_slip, vertical_slip, clearance),
                "lockHeldObserved": lock_observation,
                "measurementSource": "mujoco-mjcf-marker-sites",
                "notes": "alpha/theta intentionally blank until marker-to-angle calibration is defined",
            }
        )
    return rows


def _run_mujoco_xml_case(
    mujoco: Any,
    *,
    xml_text: str,
    case: dict[str, Any],
    steps: int,
    timestep: float,
) -> dict[str, Any]:
    model = mujoco.MjModel.from_xml_string(xml_text)
    model.opt.timestep = float(timestep)
    data = mujoco.MjData(model)

    for name, value in {
        "right_cell_alpha_slide": case.get("alphaCommand", 0.0),
        "right_cell_vertical_slide": case.get("zCommand", 0.0),
    }.items():
        actuator_id = int(
            mujoco.mj_name2id(
                model,
                mujoco.mjtObj.mjOBJ_ACTUATOR,
                name,
            )
        )
        if actuator_id >= 0:
            data.ctrl[actuator_id] = _clamp(float(value), -1.0, 1.0)

    body_ids: dict[str, int] = {}
    for name in ("left_cell", "right_cell"):
        body_id = int(
            mujoco.mj_name2id(
                model,
                mujoco.mjtObj.mjOBJ_BODY,
                name,
            )
        )
        if body_id >= 0:
            body_ids[name] = body_id

    mujoco.mj_forward(model, data)
    initial_positions = {
        name: [float(value) for value in data.xpos[body_id, 0:3].tolist()]
        for name, body_id in body_ids.items()
    }
    for _ in range(int(steps)):
        mujoco.mj_step(model, data)

    body_results = []
    for name, body_id in body_ids.items():
        initial = initial_positions[name]
        final = [float(value) for value in data.xpos[body_id, 0:3].tolist()]
        displacement = [final[index] - initial[index] for index in range(3)]
        body_results.append(
            {
                "bodyName": name,
                "bodyId": "right" if name.startswith("right") else "left",
                "initialPositionMm": initial,
                "finalPositionMm": final,
                "displacementMm": displacement,
                "displacementNormMm": float(math.sqrt(sum(value * value for value in displacement))),
            }
        )
    marker_positions = {}
    for body_name in ("left_cell", "right_cell"):
        for site, _, _ in _mjcf_site_layout(
            CAD_RAD_CELL_REFERENCE.width_x_mm,
            float(case.get("holeRadiusMm", CAD_RAD_CELL_REFERENCE.nominal_hole_diameter_mm * 0.5) or 0.0),
            float(case.get("pinRadiusMm", CAD_RAD_CELL_REFERENCE.nominal_hole_diameter_mm * 0.5) or 0.0),
        ):
            marker_name = f"{body_name}_{site}_marker"
            position = _site_world_position(mujoco, model, data, marker_name)
            if position is not None:
                marker_positions[marker_name] = position

    result = {
        "caseId": str(case.get("caseId", "")),
        "matrixIndex": case.get("matrixIndex"),
        "case": dict(case),
        "sourceEngine": "mujoco",
        "stepsCompleted": int(steps),
        "timestep": float(timestep),
        "lockMode": case.get("lockMode", ""),
        "alphaCommand": case.get("alphaCommand", 0.0),
        "zCommand": case.get("zCommand", 0.0),
        "holeRadius": case.get("holeRadius", 0.0),
        "backlash": case.get("backlash", 0.0),
        "bodies": body_results,
        "markerPositionsMm": marker_positions,
        "qpos": [float(value) for value in data.qpos.tolist()],
        "qvel": [float(value) for value in data.qvel.tolist()],
        "ctrl": [float(value) for value in data.ctrl.tolist()],
    }
    result["measurementRows"] = _case_connector_measurement_rows(result)
    return result


def two_cell_external_fidelity_mjcf_run_report(
    mjcf_index_path: str | Path | None = None,
    *,
    output_dir: str | Path | None = None,
    steps: int = 120,
    timestep: float = 0.001,
    engine_availability: dict[str, bool] | None = None,
    max_cases: int | None = None,
    case_ids: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Run the dense external-fidelity MJCF matrix with MuJoCo when available."""

    if steps <= 0:
        raise ValueError("steps must be positive")
    if max_cases is not None and max_cases < 0:
        raise ValueError("max_cases must be non-negative")

    index_path = Path(mjcf_index_path) if mjcf_index_path is not None else _external_fidelity_mjcf_index_default_path()
    index_payload: dict[str, Any] = {}
    case_files: list[dict[str, Any]] = []
    missing_evidence: list[str] = []
    errors: list[dict[str, Any]] = []

    if index_path.exists():
        try:
            index_payload = json.loads(index_path.read_text(encoding="utf-8"))
            case_files = list(index_payload.get("files", []))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append({"scope": "index", "path": str(index_path), "error": str(exc)})
            missing_evidence.append("externalFidelityMjcfIndexReadable")
    else:
        missing_evidence.append("externalFidelityMjcfIndex")

    case_count = int(index_payload.get("caseCount", len(case_files)) or len(case_files))
    requested_case_ids = [str(case_id).strip() for case_id in list(case_ids or []) if str(case_id).strip()]
    requested_case_id_set = set(requested_case_ids)
    if requested_case_id_set:
        selected_case_files = [case for case in case_files if str(case.get("caseId", "")) in requested_case_id_set]
        missing_case_ids = sorted(requested_case_id_set - {str(case.get("caseId", "")) for case in selected_case_files})
    else:
        selected_case_files = list(case_files)
        missing_case_ids = []
    selected_case_files = selected_case_files[:max_cases] if max_cases is not None else selected_case_files
    expected_body_result_count = sum(int(case.get("bodyCount", 2) or 2) for case in selected_case_files)
    expected_connector_measurement_row_count = 3 * len(selected_case_files)
    available = _engine_available("mujoco", engine_availability)
    if not available:
        missing_evidence.extend(["mujocoPythonPackage", "externalFidelityMjcfRun"])
    if case_count <= 0:
        missing_evidence.append("externalFidelityMjcfCases")
    if missing_case_ids:
        missing_evidence.append("requestedExternalFidelityCaseIds")
    if requested_case_id_set or (max_cases is not None and max_cases < case_count):
        missing_evidence.append("fullExternalFidelityMatrixRun")

    output_path = Path(output_dir) if output_dir is not None else None
    result_dir = output_path / "external_fidelity_matrix_results" if output_path is not None else None
    if result_dir is not None:
        result_dir.mkdir(parents=True, exist_ok=True)

    solver = {
        "engineAvailable": available,
        "mujocoAvailable": available,
        "ran": False,
        "stepsRequested": int(steps),
        "timestep": float(timestep),
        "engineVersion": "",
        "caseCountInIndex": case_count,
        "caseCountSelected": len(selected_case_files),
        "caseIdsRequested": requested_case_ids,
        "caseIdsMissing": missing_case_ids,
        "caseRunCount": 0,
        "caseErrorCount": 0,
        "maxCases": max_cases,
        "outputDir": str(output_path) if output_path is not None else "",
    }
    case_results: list[dict[str, Any]] = []

    if available and index_path.exists() and selected_case_files:
        try:
            mujoco = importlib.import_module("mujoco")
            solver["engineVersion"] = str(getattr(mujoco, "__version__", ""))
            for case in selected_case_files:
                case_id = str(case.get("caseId", ""))
                xml_path = _resolve_index_relative_path(index_path, str(case.get("path", "")))
                try:
                    xml_text = xml_path.read_text(encoding="utf-8")
                    result = _run_mujoco_xml_case(
                        mujoco,
                        xml_text=xml_text,
                        case=case,
                        steps=steps,
                        timestep=timestep,
                    )
                    case_results.append(result)
                    solver["caseRunCount"] = int(solver["caseRunCount"]) + 1
                    if result_dir is not None:
                        (result_dir / f"{case_id}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
                except Exception as exc:  # pragma: no cover - depends on optional external engine
                    errors.append({"scope": "case", "caseId": case_id, "path": str(xml_path), "error": str(exc)})
        except Exception as exc:  # pragma: no cover - depends on optional external engine
            errors.append({"scope": "engine", "engine": "mujoco", "error": str(exc)})
            missing_evidence.append("mujocoExecution")

    solver["caseErrorCount"] = len(errors)
    body_result_count = sum(len(result.get("bodies", [])) for result in case_results)
    measurement_rows = [
        dict(row)
        for result in case_results
        for row in result.get("measurementRows", [])
        if isinstance(row, dict)
    ]
    if len(case_results) < len(selected_case_files):
        missing_evidence.append("mujocoCaseResults")
    if body_result_count < expected_body_result_count:
        missing_evidence.append("mujocoBodyResults")
    if len(measurement_rows) < expected_connector_measurement_row_count:
        missing_evidence.append("mujocoConnectorMeasurements")

    missing_evidence = sorted(set(missing_evidence))
    run_complete = (
        available
        and not errors
        and max_cases is None
        and case_count > 0
        and len(case_results) == case_count
        and body_result_count >= expected_body_result_count
        and len(measurement_rows) >= expected_connector_measurement_row_count
        and not missing_evidence
    )
    partial_run = available and len(case_results) > 0 and not run_complete
    solver["ran"] = bool(run_complete or partial_run)
    status = (
        "external-fidelity-mjcf-run-complete"
        if run_complete
        else "external-fidelity-mjcf-partial-run"
        if partial_run
        else "needs-external-fidelity-mjcf-run"
    )

    return {
        "schema": TWO_CELL_EXTERNAL_FIDELITY_MJCF_RUN_SCHEMA,
        "method": "optional MuJoCo execution of the dense two-cell external-fidelity MJCF matrix",
        "engine": "MuJoCo",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "mjcfIndex": {
            "path": str(index_path),
            "exists": index_path.exists(),
            "schema": index_payload.get("schema", ""),
            "caseCount": case_count,
            "directory": index_payload.get("directory", ""),
        },
        "solver": solver,
        "results": {
            "cases": case_results,
            "measurementRows": measurement_rows,
            "resultDirectory": str(result_dir) if result_dir is not None else "",
        },
        "errors": errors,
        "summary": {
            "status": status,
            "externalFidelityMjcfRunComplete": run_complete,
            "externalFidelityMjcfPartialRun": partial_run,
            "mujocoAvailable": available,
            "caseCount": case_count,
            "caseSelectedCount": len(selected_case_files),
            "caseRunCount": len(case_results),
            "bodyResultCount": body_result_count,
            "expectedBodyResultCount": expected_body_result_count,
            "connectorMeasurementRowCount": len(measurement_rows),
            "expectedConnectorMeasurementRowCount": expected_connector_measurement_row_count,
            "physicalAccuracyValidated": False,
            "missingEvidenceCount": len(missing_evidence),
            "missingEvidence": missing_evidence,
        },
        "claimLabels": {
            "run": "external-engine-derived batch result" if run_complete else "external batch run not complete",
            "geometry": "dense MJCF proxies generated from A360-dimensional reduced two-cell cases",
            "physicalAccuracy": "false until segmented CAD/contact setup and bench coordinates agree",
        },
        "limitations": [
            "The dense batch runs generated MJCF proxy cases, not exact segmented Fusion rigid bodies.",
            "The runner records body-center results; connector-level benchmark rows still require marker extraction or bench measurement.",
            "A completed MuJoCo batch remains a model result until friction, stiffness, gravity scaling, and lock behavior are calibrated.",
        ],
    }


def solve_two_cell_quasistatic(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    hole_radius: float | None = None,
    pin_radius: float | None = None,
    lock_mode: str | None = None,
    axial_stiffness: float = 90.0,
    hinge_stiffness: float = 6.0,
    alpha_actuator_stiffness: float = 55.0,
    z_actuator_stiffness: float = 45.0,
    neighbor_residual_stiffness: float = 4.0,
    contact_stiffness: float = 600.0,
    ground_stiffness: float = 500.0,
    gravity_force: float = 0.025,
    tolerance: float = 1e-9,
    maxiter: int = 800,
    **control_overrides: Any,
) -> dict[str, Any]:
    """Solve a reduced two-cell quasistatic energy with hard lock constraints."""

    minimize = None
    solver_engine = "rad_sim.bounded_pattern_search"
    solver_method = "coordinate-pattern-search"
    try:
        from scipy.optimize import minimize
        solver_engine = "scipy.optimize.minimize"
        solver_method = "L-BFGS-B"
    except ImportError:
        pass

    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls, **control_overrides)
    if lock_mode is not None:
        controls = _lock_mode_controls(controls, lock_mode)
    radii = _radii(config, hole_radius=hole_radius, pin_radius=pin_radius)
    zones = _dead_zones(config, radii)
    pitch = float(config.cell_size)
    initial_alpha = float(config.initial_alpha)
    alpha_drive = _dead_zone(controls.alpha_command, zones["alphaDeadZone"])
    z_drive = _dead_zone(controls.z_command, zones["verticalDeadZone"])
    alpha_transmission = config.coupling_gain * math.exp(-2.25 * zones["alphaDeadZone"])
    z_transmission = config.z_coupling_gain * math.exp(-3.1 * zones["verticalDeadZone"])
    left_alpha_target = _clamp(
        initial_alpha + (0.0 if controls.left_position_locked else alpha_drive * alpha_transmission),
        config.alpha_min,
        config.alpha_max,
    )
    right_alpha_target = _clamp(
        initial_alpha if controls.right_locked else initial_alpha + alpha_drive,
        config.alpha_min,
        config.alpha_max,
    )
    left_z_target = 0.0 if controls.left_position_locked else z_drive * z_transmission
    right_z_target = z_drive
    rest_pitch_target = _two_cell_rest_pitch(config, left_alpha_target, right_alpha_target)
    initial = {
        "leftX": 0.0,
        "leftZ": left_z_target,
        "leftAlpha": left_alpha_target,
        "rightX": rest_pitch_target,
        "rightZ": right_z_target,
        "rightAlpha": right_alpha_target,
    }
    fixed = {}
    if controls.left_position_locked:
        fixed["leftX"] = 0.0
        fixed["leftZ"] = 0.0
    if controls.right_position_locked:
        fixed["rightX"] = pitch
        fixed["rightZ"] = 0.0
    if controls.right_locked:
        fixed["rightAlpha"] = initial_alpha
    names = ["leftX", "leftZ", "leftAlpha", "rightX", "rightZ", "rightAlpha"]
    bounds_by_name = {
        "leftX": (-2.0 * pitch, 2.0 * pitch),
        "rightX": (-1.0 * pitch, 3.0 * pitch),
        "leftZ": (-2.0, 2.0),
        "rightZ": (-2.0, 2.0),
        "leftAlpha": (config.alpha_min, config.alpha_max),
        "rightAlpha": (config.alpha_min, config.alpha_max),
    }
    free_names = [name for name in names if name not in fixed]
    x0 = [float(initial[name]) for name in free_names]
    bounds = [bounds_by_name[name] for name in free_names]

    def unpack(values: list[float] | tuple[float, ...]) -> dict[str, float]:
        out = {name: float(initial[name]) for name in names}
        out.update({name: float(value) for name, value in zip(free_names, values)})
        out.update(fixed)
        return out

    def energy_breakdown(q: dict[str, float]) -> dict[str, float]:
        dx = q["rightX"] - q["leftX"]
        dz = q["rightZ"] - q["leftZ"]
        distance = math.hypot(dx, dz)
        rest_pitch = _two_cell_rest_pitch(config, q["leftAlpha"], q["rightAlpha"])
        axial_gap = zones["effectiveBacklash"] * pitch
        axial_error = distance - rest_pitch
        axial_excess = _dead_zone(axial_error, axial_gap)
        vertical_clearance = zones["verticalDeadZone"]
        vertical_excess = max(0.0, abs(dz) - vertical_clearance)
        left_floor_penetration = max(0.0, -(q["leftZ"] + vertical_clearance))
        right_floor_penetration = max(0.0, -(q["rightZ"] + vertical_clearance))
        theta_delta = math.radians(float(alpha_to_theta(q["rightAlpha"])) - float(alpha_to_theta(q["leftAlpha"])))
        return {
            "axialEnergy": 0.5 * max(0.0, axial_stiffness) * axial_excess**2,
            "hingeEnergy": 0.5 * max(0.0, hinge_stiffness) * theta_delta**2,
            "rightAlphaActuatorEnergy": 0.5
            * max(0.0, alpha_actuator_stiffness)
            * (q["rightAlpha"] - right_alpha_target) ** 2,
            "rightZActuatorEnergy": 0.5
            * max(0.0, z_actuator_stiffness)
            * (q["rightZ"] - right_z_target) ** 2,
            "leftResidualEnergy": 0.5
            * max(0.0, neighbor_residual_stiffness)
            * ((q["leftAlpha"] - left_alpha_target) ** 2 + (q["leftZ"] - left_z_target) ** 2),
            "verticalContactEnergy": 0.5 * max(0.0, contact_stiffness) * vertical_excess**2,
            "groundContactEnergy": 0.5
            * max(0.0, ground_stiffness)
            * (left_floor_penetration**2 + right_floor_penetration**2),
            "gravityPotential": max(0.0, gravity_force) * (q["leftZ"] + q["rightZ"]),
            "axialError": axial_error,
            "axialExcess": axial_excess,
            "verticalShear": dz,
            "verticalExcess": vertical_excess,
            "leftFloorPenetration": left_floor_penetration,
            "rightFloorPenetration": right_floor_penetration,
            "distance": distance,
            "restPitch": rest_pitch,
        }

    def objective(values: list[float] | tuple[float, ...]) -> float:
        breakdown = energy_breakdown(unpack(values))
        return (
            breakdown["axialEnergy"]
            + breakdown["hingeEnergy"]
            + breakdown["rightAlphaActuatorEnergy"]
            + breakdown["rightZActuatorEnergy"]
            + breakdown["leftResidualEnergy"]
            + breakdown["verticalContactEnergy"]
            + breakdown["groundContactEnergy"]
            + breakdown["gravityPotential"]
        )

    if free_names and minimize is not None:
        result = minimize(
            objective,
            x0,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": max(1, int(maxiter)), "ftol": max(float(tolerance), 1e-15)},
        )
        q = unpack([float(value) for value in result.x])
        success = bool(result.success)
        message = str(result.message)
        iterations = int(getattr(result, "nit", 0))
    elif free_names:
        result = _bounded_pattern_search(objective, x0, bounds, tolerance=tolerance, maxiter=maxiter)
        q = unpack([float(value) for value in result["x"]])
        success = bool(result["success"])
        message = str(result["message"])
        iterations = int(result["nit"])
    else:
        q = unpack(())
        success = True
        message = "all degrees of freedom fixed"
        iterations = 0
    breakdown = energy_breakdown(q)
    stored_energy = (
        breakdown["axialEnergy"]
        + breakdown["hingeEnergy"]
        + breakdown["rightAlphaActuatorEnergy"]
        + breakdown["rightZActuatorEnergy"]
        + breakdown["leftResidualEnergy"]
        + breakdown["verticalContactEnergy"]
        + breakdown["groundContactEnergy"]
    )
    total_energy = stored_energy + breakdown["gravityPotential"]
    left_theta = float(alpha_to_theta(q["leftAlpha"]))
    right_theta = float(alpha_to_theta(q["rightAlpha"]))
    left_fixture_error = math.hypot(q["leftX"], q["leftZ"]) if controls.left_position_locked else None
    right_fixture_error = math.hypot(q["rightX"] - pitch, q["rightZ"]) if controls.right_position_locked else None
    right_state_alpha_error = abs(q["rightAlpha"] - initial_alpha) if controls.right_locked else None
    contact_parts = []
    if abs(breakdown["verticalExcess"]) <= max(tolerance, 1e-12):
        contact_parts.append("inside-vertical-clearance")
    else:
        contact_parts.append("pin-hole-vertical-contact")
    if abs(breakdown["axialExcess"]) > max(tolerance, 1e-12):
        contact_parts.append("axial-backlash-contact")
    if breakdown["leftFloorPenetration"] > max(tolerance, 1e-12) or breakdown["rightFloorPenetration"] > max(tolerance, 1e-12):
        contact_parts.append("ground-contact")
    if controls.right_position_locked:
        contact_parts.append("right-position-locked")
    if controls.right_locked:
        contact_parts.append("right-state-locked")
    return {
        "schema": TWO_CELL_QUASISTATIC_SCHEMA,
        "model": "two-cell-hard-lock-quasistatic-clearance-energy",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "controls": {
            "alphaCommand": controls.alpha_command,
            "zCommand": controls.z_command,
            "leftPositionLocked": controls.left_position_locked,
            "rightLocked": controls.right_locked,
            "rightPositionLocked": controls.right_position_locked,
        },
        "dimensions": {"cellSize": pitch, **radii},
        "deadZones": zones,
        "targets": {
            "leftAlpha": left_alpha_target,
            "rightAlpha": right_alpha_target,
            "leftZ": left_z_target,
            "rightZ": right_z_target,
            "restPitch": rest_pitch_target,
        },
        "cells": [
            {
                "id": "left",
                "center": {"x": q["leftX"], "y": 0.0, "z": q["leftZ"]},
                "alpha": q["leftAlpha"],
                "theta": left_theta,
                "positionLocked": controls.left_position_locked,
                "stateLocked": False,
            },
            {
                "id": "right",
                "center": {"x": q["rightX"], "y": 0.0, "z": q["rightZ"]},
                "alpha": q["rightAlpha"],
                "theta": right_theta,
                "positionLocked": controls.right_position_locked,
                "stateLocked": controls.right_locked,
            },
        ],
        "contact": {
            "mode": "+".join(contact_parts),
            "distance": breakdown["distance"],
            "restPitch": breakdown["restPitch"],
            "axialBacklashGap": zones["effectiveBacklash"] * pitch,
            "axialError": breakdown["axialError"],
            "axialExcess": breakdown["axialExcess"],
            "verticalClearance": zones["verticalDeadZone"],
            "verticalShear": breakdown["verticalShear"],
            "verticalExcess": breakdown["verticalExcess"],
            "leftFloorPenetration": breakdown["leftFloorPenetration"],
            "rightFloorPenetration": breakdown["rightFloorPenetration"],
        },
        "energy": {
            **{
                key: value
                for key, value in breakdown.items()
                if key.endswith("Energy") or key == "gravityPotential"
            },
            "storedEnergy": stored_energy,
            "totalEnergy": total_energy,
        },
        "constraints": {
            "hardFixedDofs": fixed,
            "freeDofs": free_names,
            "leftFixtureError": left_fixture_error,
            "rightFixtureError": right_fixture_error,
            "rightStateAlphaError": right_state_alpha_error,
            "positionLocksHold": (
                (left_fixture_error is None or left_fixture_error <= 1e-8)
                and (right_fixture_error is None or right_fixture_error <= 1e-8)
            ),
            "stateLocksHold": right_state_alpha_error is None or right_state_alpha_error <= 1e-8,
            "rightLockedVerticalMotionAllowed": bool(
                controls.right_locked and not controls.right_position_locked and abs(q["rightZ"]) > 1e-9
            ),
        },
        "solver": {
            "engine": solver_engine,
            "method": solver_method,
            "success": success,
            "message": message,
            "iterations": iterations,
            "tolerance": tolerance,
            "maxiter": maxiter,
            "scipyAvailable": minimize is not None,
        },
        "claimLabels": {
            "physics": "internal quasistatic energy solve with hard coordinate constraints",
            "contact": "clearance penalty proxy, not collision-resolved CAD contact",
            "accuracy": "requires external solver and bench measurement comparison before physical claims",
        },
    }


def _site_by_name(layout: dict[str, Any], site_name: str) -> dict[str, Any]:
    for site in layout.get("padSites", []):
        if isinstance(site, dict) and site.get("name") == site_name:
            return site
    raise KeyError(f"unknown CAD site {site_name!r}")


def _connector_site_position_mm(
    *,
    cell: dict[str, Any],
    site: dict[str, Any],
    scale: dict[str, float],
    config: LatticeConfig,
    theta_yaw_gain: float,
) -> dict[str, float]:
    center = _scale_center_mm(cell["center"], scale)
    alpha_scale = math.sqrt(max(float(cell["alpha"]), 1e-9) / max(float(config.initial_alpha), 1e-9))
    theta_delta_degrees = float(cell["theta"]) - float(alpha_to_theta(config.initial_alpha))
    yaw = math.radians(theta_delta_degrees * theta_yaw_gain)
    local_x = float(site["xMm"]) * alpha_scale
    local_y = float(site["yMm"]) * alpha_scale
    cos_yaw = math.cos(yaw)
    sin_yaw = math.sin(yaw)
    return {
        "x": center[0] + cos_yaw * local_x - sin_yaw * local_y,
        "y": center[1] + sin_yaw * local_x + cos_yaw * local_y,
        "z": center[2],
    }


def two_cell_connector_contact_report(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    hole_radius: float | None = None,
    pin_radius: float | None = None,
    lock_mode: str | None = None,
    lateral_contact_stiffness: float = 12.0,
    vertical_contact_stiffness: float = 8.0,
    theta_yaw_gain: float = 0.18,
    tolerance: float = 1e-9,
    **solve_overrides: Any,
) -> dict[str, Any]:
    """Evaluate CAD-derived connector-level pin-hole slip for a two-cell pair."""

    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    if lock_mode is not None:
        controls = _lock_mode_controls(controls, lock_mode)
    solve = solve_two_cell_quasistatic(
        config,
        controls,
        hole_radius=hole_radius,
        pin_radius=pin_radius,
        tolerance=tolerance,
        **solve_overrides,
    )
    layout = cad_rad_cell_layout(
        config,
        hole_radius=solve["dimensions"]["holeRadius"],
        pin_radius=solve["dimensions"]["pinRadius"],
    )
    scale = _two_cell_scale_mm(
        config,
        {
            **solve,
            "dimensions": {
                "pinRadius": solve["dimensions"]["pinRadius"],
                "holeRadius": solve["dimensions"]["holeRadius"],
            },
        },
    )
    left_cell, right_cell = solve["cells"]
    clearance = float(layout["dimensionsMm"]["pinHoleClearanceMm"])
    connector_rows = []
    total_penalty = 0.0
    active_contacts = 0
    for pair in layout["twoCellConnectorPairs"]:
        left_site = _site_by_name(layout, str(pair["leftSite"]))
        right_site = _site_by_name(layout, str(pair["rightSite"]))
        left_pos = _connector_site_position_mm(
            cell=left_cell,
            site=left_site,
            scale=scale,
            config=config,
            theta_yaw_gain=theta_yaw_gain,
        )
        right_pos = _connector_site_position_mm(
            cell=right_cell,
            site=right_site,
            scale=scale,
            config=config,
            theta_yaw_gain=theta_yaw_gain,
        )
        slip_x = right_pos["x"] - left_pos["x"]
        slip_y = right_pos["y"] - left_pos["y"]
        slip_z = right_pos["z"] - left_pos["z"]
        lateral_slip = math.hypot(slip_x, slip_y)
        vertical_slip = abs(slip_z)
        total_slip = math.sqrt(lateral_slip**2 + slip_z**2)
        lateral_excess = max(0.0, lateral_slip - clearance)
        vertical_excess = max(0.0, vertical_slip - clearance)
        lateral_active = lateral_excess > tolerance
        vertical_active = vertical_excess > tolerance
        if lateral_active or vertical_active:
            active_contacts += 1
        penalty = 0.5 * float(lateral_contact_stiffness) * lateral_excess**2 + 0.5 * float(
            vertical_contact_stiffness
        ) * vertical_excess**2
        total_penalty += penalty
        if lateral_active and vertical_active:
            mode = "lateral+vertical-contact"
        elif lateral_active:
            mode = "lateral-contact"
        elif vertical_active:
            mode = "vertical-contact"
        else:
            mode = "inside-clearance"
        connector_rows.append(
            {
                "connector": pair["name"],
                "leftSite": pair["leftSite"],
                "rightSite": pair["rightSite"],
                "leftPositionMm": left_pos,
                "rightPositionMm": right_pos,
                "slipMm": {"x": slip_x, "y": slip_y, "z": slip_z},
                "lateralSlipMm": lateral_slip,
                "verticalSlipMm": vertical_slip,
                "totalSlipMm": total_slip,
                "lateralClearanceMm": clearance,
                "verticalFreePlayMm": clearance,
                "lateralExcessMm": lateral_excess,
                "verticalExcessMm": vertical_excess,
                "contactPenalty": penalty,
                "contactMode": mode,
            }
        )
    max_lateral_slip = max((row["lateralSlipMm"] for row in connector_rows), default=0.0)
    max_vertical_slip = max((row["verticalSlipMm"] for row in connector_rows), default=0.0)
    max_lateral_excess = max((row["lateralExcessMm"] for row in connector_rows), default=0.0)
    max_vertical_excess = max((row["verticalExcessMm"] for row in connector_rows), default=0.0)
    constraints = solve["constraints"]
    ready = (
        bool(solve["solver"]["success"])
        and len(connector_rows) == len(CAD_RAD_X_CONNECTOR_PAIRS)
        and bool(constraints["positionLocksHold"])
        and bool(constraints["stateLocksHold"])
    )
    return {
        "schema": TWO_CELL_CONNECTOR_CONTACT_SCHEMA,
        "model": "cad-layout-two-cell-connector-clearance-contact",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "cadLayout": layout,
        "scale": scale,
        "controls": solve["controls"],
        "dimensions": {
            **solve["dimensions"],
            "connectorPitchMm": scale["connectorPitchMm"],
            "pinRadiusMm": scale["pinRadiusMm"],
            "holeRadiusMm": scale["holeRadiusMm"],
            "pinHoleClearanceMm": clearance,
        },
        "quasistatic": solve,
        "connectors": connector_rows,
        "summary": {
            "status": "connector-contact-internal-ready" if ready else "needs-connector-contact-inputs",
            "connectorContactReady": ready,
            "connectorCount": len(connector_rows),
            "activeContactCount": active_contacts,
            "insideClearanceCount": sum(1 for row in connector_rows if row["contactMode"] == "inside-clearance"),
            "maxLateralSlipMm": max_lateral_slip,
            "maxVerticalSlipMm": max_vertical_slip,
            "maxLateralExcessMm": max_lateral_excess,
            "maxVerticalExcessMm": max_vertical_excess,
            "totalContactPenalty": total_penalty,
            "solverSuccess": bool(solve["solver"]["success"]),
            "positionLocksHold": bool(constraints["positionLocksHold"]),
            "stateLocksHold": bool(constraints["stateLocksHold"]),
            "rightLockedVerticalMotionAllowed": bool(constraints["rightLockedVerticalMotionAllowed"]),
            "requiresExternalEngineRun": True,
            "requiresMeasuredData": True,
        },
        "claimLabels": {
            "connectorGeometry": "CAD-layout pin-hole site inventory with connector-pitch scaling",
            "contact": "clearance-excess penalty over three connector pairs",
            "physicalAccuracy": "internal diagnostic until segmented CAD contact or measured connector slips agree",
        },
        "limitations": [
            "Connector sites are derived from the A360 envelope and preview topology, not exact BREP vertices.",
            "Alpha changes scale connector radius with a normalized sqrt(alpha) law until CAD mechanism kinematics are measured.",
            "Friction, crown compliance, screw preload, and pin/hole wall contact are not externally solved here.",
        ],
    }


def sweep_two_cell_connector_contact(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    alpha_commands: list[float] | tuple[float, ...] | None = None,
    z_commands: list[float] | tuple[float, ...] | None = None,
    hole_radii: list[float] | tuple[float, ...] | None = None,
    lock_modes: list[str] | tuple[str, ...] | None = None,
    **report_overrides: Any,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    radii = _radii(config)
    alpha_values = _unique_floats(list(alpha_commands or (-0.5, 0.0, 0.5)))
    z_values = _unique_floats(list(z_commands or (-0.35, 0.0, 0.35)))
    hole_values = _unique_floats(
        list(hole_radii or (radii["pinRadius"], radii["holeRadius"], controls.hole_sweep_max))
    )
    modes = list(lock_modes or ("free", "right_state_locked", "right_position_locked", "left_free"))
    rows = []
    connector_rows = []
    index = 0
    for mode in modes:
        mode_base = _lock_mode_controls(controls, mode)
        for hole in hole_values:
            for alpha in alpha_values:
                for z in z_values:
                    case_controls = _controls(mode_base, alpha_command=alpha, z_command=z)
                    report = two_cell_connector_contact_report(
                        config,
                        case_controls,
                        hole_radius=hole,
                        **report_overrides,
                    )
                    summary = report["summary"]
                    row = {
                        "index": index,
                        "lockMode": mode,
                        "alphaCommand": alpha,
                        "zCommand": z,
                        "pinRadius": report["dimensions"]["pinRadius"],
                        "holeRadius": report["dimensions"]["holeRadius"],
                        "clearance": report["dimensions"]["pinHoleClearance"],
                        "pinHoleClearanceMm": report["dimensions"]["pinHoleClearanceMm"],
                        "connectorPitchMm": report["dimensions"]["connectorPitchMm"],
                        "connectorCount": summary["connectorCount"],
                        "activeContactCount": summary["activeContactCount"],
                        "maxLateralSlipMm": summary["maxLateralSlipMm"],
                        "maxVerticalSlipMm": summary["maxVerticalSlipMm"],
                        "maxLateralExcessMm": summary["maxLateralExcessMm"],
                        "maxVerticalExcessMm": summary["maxVerticalExcessMm"],
                        "totalContactPenalty": summary["totalContactPenalty"],
                        "solverSuccess": summary["solverSuccess"],
                        "positionLocksHold": summary["positionLocksHold"],
                        "stateLocksHold": summary["stateLocksHold"],
                    }
                    rows.append(row)
                    for connector in report["connectors"]:
                        connector_rows.append(
                            {
                                "index": index,
                                "lockMode": mode,
                                "alphaCommand": alpha,
                                "zCommand": z,
                                "holeRadius": report["dimensions"]["holeRadius"],
                                "pinHoleClearanceMm": report["dimensions"]["pinHoleClearanceMm"],
                                "connector": connector["connector"],
                                "leftSite": connector["leftSite"],
                                "rightSite": connector["rightSite"],
                                "pinHoleClearanceMm": report["dimensions"]["pinHoleClearanceMm"],
                                "leftPositionMm": connector["leftPositionMm"],
                                "rightPositionMm": connector["rightPositionMm"],
                                "slipMm": connector["slipMm"],
                                "lateralSlipMm": connector["lateralSlipMm"],
                                "verticalSlipMm": connector["verticalSlipMm"],
                                "totalSlipMm": connector["totalSlipMm"],
                                "lateralExcessMm": connector["lateralExcessMm"],
                                "verticalExcessMm": connector["verticalExcessMm"],
                                "contactPenalty": connector["contactPenalty"],
                                "contactMode": connector["contactMode"],
                            }
                        )
                    index += 1
    failures = [row for row in rows if not row["solverSuccess"]]
    return {
        "schema": TWO_CELL_CONNECTOR_CONTACT_SWEEP_SCHEMA,
        "model": "two-cell-connector-clearance-lock-actuation-sweep",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "grid": {
            "alphaCommands": alpha_values,
            "zCommands": z_values,
            "holeRadii": hole_values,
            "lockModes": modes,
        },
        "rows": rows,
        "connectorRows": connector_rows,
        "summary": {
            "rowCount": len(rows),
            "connectorRowCount": len(connector_rows),
            "solverSuccessCount": len(rows) - len(failures),
            "solverFailureCount": len(failures),
            "maxLateralSlipMm": max((row["maxLateralSlipMm"] for row in rows), default=0.0),
            "maxVerticalSlipMm": max((row["maxVerticalSlipMm"] for row in rows), default=0.0),
            "maxLateralExcessMm": max((row["maxLateralExcessMm"] for row in rows), default=0.0),
            "maxVerticalExcessMm": max((row["maxVerticalExcessMm"] for row in rows), default=0.0),
            "maxContactPenalty": max((row["totalContactPenalty"] for row in rows), default=0.0),
            "requiresExternalEngineRun": True,
            "requiresMeasuredData": True,
        },
    }


def _two_cell_case_config(
    config: LatticeConfig,
    *,
    backlash: float | None = None,
    hole_radius: float | None = None,
) -> LatticeConfig:
    hole = float(config.hole_radius if hole_radius is None else hole_radius)
    hole = max(float(config.pin_radius), hole)
    return replace(
        config,
        rows=1,
        cols=2,
        backlash=max(0.0, float(config.backlash if backlash is None else backlash)),
        hole_radius=hole,
    )


def two_cell_physical_simulation_suite(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    gravity_force: float = 0.025,
    tolerance: float = 1e-7,
) -> dict[str, Any]:
    """Run named two-cell reduced-physics tests for RAD calibration planning."""

    config = config or LatticeConfig(rows=1, cols=2)
    config = replace(config, rows=1, cols=2)
    controls = _controls(controls)
    tight_hole = float(config.pin_radius)
    nominal_hole = float(config.hole_radius)
    loose_hole = max(float(controls.hole_sweep_max), nominal_hole)
    zero_backlash = 0.0
    high_backlash = max(float(config.backlash) * 3.0, 0.25)
    alpha_drive = controls.alpha_command if abs(controls.alpha_command) > 0.05 else -0.35
    z_drive = controls.z_command if abs(controls.z_command) > 0.05 else 0.35
    case_specs = [
        {
            "caseId": "free_contract_lift_nominal",
            "purpose": "baseline two-cell contraction plus vertical lift",
            "lockMode": "free",
            "alpha": alpha_drive,
            "z": z_drive,
            "hole": nominal_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "tight_clearance_contract_lift",
            "purpose": "minimum pin-hole clearance response",
            "lockMode": "free",
            "alpha": alpha_drive,
            "z": z_drive,
            "hole": tight_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "loose_clearance_contract_lift",
            "purpose": "large pin-hole clearance / backlash response",
            "lockMode": "free",
            "alpha": alpha_drive,
            "z": z_drive,
            "hole": loose_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "free_expand_pushdown_nominal",
            "purpose": "opposite alpha and vertical actuation polarity",
            "lockMode": "free",
            "alpha": abs(alpha_drive),
            "z": -abs(z_drive),
            "hole": nominal_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "alpha_only_contract",
            "purpose": "horizontal/dilation actuation without vertical command",
            "lockMode": "free",
            "alpha": alpha_drive,
            "z": 0.0,
            "hole": nominal_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "z_only_lift",
            "purpose": "vertical actuation and residual contact without alpha command",
            "lockMode": "free",
            "alpha": 0.0,
            "z": abs(z_drive),
            "hole": nominal_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "right_state_locked_lift",
            "purpose": "state lock freezes alpha while vertical pin-hole motion remains free",
            "lockMode": "right_state_locked",
            "alpha": alpha_drive,
            "z": abs(z_drive),
            "hole": nominal_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "right_position_locked_lift",
            "purpose": "position lock fixes fixture coordinates under actuation",
            "lockMode": "right_position_locked",
            "alpha": alpha_drive,
            "z": abs(z_drive),
            "hole": nominal_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "left_free_neighbor_residual_lift",
            "purpose": "residual vertical motion transmitted into the neighboring free cell",
            "lockMode": "left_free",
            "alpha": alpha_drive,
            "z": abs(z_drive),
            "hole": nominal_hole,
            "backlash": config.backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "zero_backlash_left_free_lift",
            "purpose": "neighbor response with backlash removed",
            "lockMode": "left_free",
            "alpha": alpha_drive,
            "z": abs(z_drive),
            "hole": nominal_hole,
            "backlash": zero_backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "high_backlash_left_free_lift",
            "purpose": "neighbor response after increasing backlash",
            "lockMode": "left_free",
            "alpha": alpha_drive,
            "z": abs(z_drive),
            "hole": nominal_hole,
            "backlash": high_backlash,
            "gravity": gravity_force,
        },
        {
            "caseId": "gravity_sag_loose_clearance",
            "purpose": "zero-command gravity sag with large pin-hole clearance",
            "lockMode": "free",
            "alpha": 0.0,
            "z": 0.0,
            "hole": loose_hole,
            "backlash": config.backlash,
            "gravity": max(0.5, float(gravity_force)),
        },
    ]

    rows: list[dict[str, Any]] = []
    detailed_cases: list[dict[str, Any]] = []
    for index, spec in enumerate(case_specs):
        case_config = _two_cell_case_config(
            config,
            backlash=float(spec["backlash"]),
            hole_radius=float(spec["hole"]),
        )
        case_controls = _lock_mode_controls(
            _controls(controls, alpha_command=float(spec["alpha"]), z_command=float(spec["z"])),
            str(spec["lockMode"]),
        )
        bench = simulate_two_cell_bench(case_config, case_controls)
        contact = two_cell_connector_contact_report(
            case_config,
            case_controls,
            gravity_force=float(spec["gravity"]),
            tolerance=tolerance,
        )
        solve = contact["quasistatic"]
        left, right = solve["cells"]
        summary = contact["summary"]
        constraints = solve["constraints"]
        lock_mode = str(spec["lockMode"])
        if lock_mode == "right_position_locked":
            lock_pass = bool(constraints["positionLocksHold"]) and float(
                constraints["rightFixtureError"] or 0.0
            ) <= tolerance
        elif lock_mode == "right_state_locked":
            lock_pass = bool(constraints["stateLocksHold"]) and float(
                constraints["rightStateAlphaError"] or 0.0
            ) <= tolerance
        else:
            lock_pass = True
        row = {
            "index": index,
            "caseId": spec["caseId"],
            "purpose": spec["purpose"],
            "lockMode": lock_mode,
            "backlash": case_config.backlash,
            "alphaCommand": case_controls.alpha_command,
            "zCommand": case_controls.z_command,
            "gravityForce": float(spec["gravity"]),
            "pinRadius": solve["dimensions"]["pinRadius"],
            "holeRadius": solve["dimensions"]["holeRadius"],
            "clearance": solve["dimensions"]["pinHoleClearance"],
            "effectiveBacklash": solve["deadZones"]["effectiveBacklash"],
            "pinHoleClearanceMm": contact["dimensions"]["pinHoleClearanceMm"],
            "solverSuccess": bool(solve["solver"]["success"]),
            "lockPass": lock_pass,
            "positionLocksHold": bool(constraints["positionLocksHold"]),
            "stateLocksHold": bool(constraints["stateLocksHold"]),
            "rightLockedVerticalMotionAllowed": bool(constraints["rightLockedVerticalMotionAllowed"]),
            "leftX": left["center"]["x"],
            "leftZ": left["center"]["z"],
            "rightX": right["center"]["x"],
            "rightZ": right["center"]["z"],
            "leftAlpha": left["alpha"],
            "rightAlpha": right["alpha"],
            "leftAlphaDeltaFromInitial": left["alpha"] - case_config.initial_alpha,
            "rightAlphaDeltaFromInitial": right["alpha"] - case_config.initial_alpha,
            "leftTheta": left["theta"],
            "rightTheta": right["theta"],
            "benchLeftResidualZ": bench["cells"][0]["residualZ"],
            "verticalShear": solve["contact"]["verticalShear"],
            "verticalExcess": solve["contact"]["verticalExcess"],
            "axialExcess": solve["contact"]["axialExcess"],
            "contactMode": solve["contact"]["mode"],
            "connectorActiveContactCount": summary["activeContactCount"],
            "connectorMaxVerticalSlipMm": summary["maxVerticalSlipMm"],
            "connectorMaxVerticalExcessMm": summary["maxVerticalExcessMm"],
            "connectorTotalContactPenalty": summary["totalContactPenalty"],
            "totalEnergy": solve["energy"]["totalEnergy"],
            "externalMeasurementRequired": True,
        }
        rows.append(row)
        detailed_cases.append(
            {
                "caseId": spec["caseId"],
                "purpose": spec["purpose"],
                "config": {
                    "backlash": case_config.backlash,
                    "pinRadius": case_config.pin_radius,
                    "holeRadius": case_config.hole_radius,
                },
                "controls": solve["controls"],
                "quasistatic": solve,
                "connectorContact": contact,
            }
        )

    by_id = {str(row["caseId"]): row for row in rows}
    tight = by_id["tight_clearance_contract_lift"]
    loose = by_id["loose_clearance_contract_lift"]
    zero = by_id["zero_backlash_left_free_lift"]
    high = by_id["high_backlash_left_free_lift"]
    clearance_trend = {
        "tightClearanceMm": tight["pinHoleClearanceMm"],
        "looseClearanceMm": loose["pinHoleClearanceMm"],
        "looseHasLargerClearance": loose["pinHoleClearanceMm"] >= tight["pinHoleClearanceMm"],
        "tightContactPenalty": tight["connectorTotalContactPenalty"],
        "looseContactPenalty": loose["connectorTotalContactPenalty"],
        "loosePenaltyNotHigher": loose["connectorTotalContactPenalty"] <= tight["connectorTotalContactPenalty"] + tolerance,
    }
    backlash_trend = {
        "zeroBacklash": zero["backlash"],
        "highBacklash": high["backlash"],
        "zeroNeighborResidualZ": zero["benchLeftResidualZ"],
        "highNeighborResidualZ": high["benchLeftResidualZ"],
        "zeroLeftAlphaDelta": zero["leftAlphaDeltaFromInitial"],
        "highLeftAlphaDelta": high["leftAlphaDeltaFromInitial"],
        "highBacklashReducesNeighborResidual": abs(high["benchLeftResidualZ"])
        <= abs(zero["benchLeftResidualZ"]) + tolerance,
        "highBacklashReducesAlphaResponse": abs(high["leftAlphaDeltaFromInitial"])
        <= abs(zero["leftAlphaDeltaFromInitial"]) + tolerance,
    }
    solver_success_count = sum(1 for row in rows if row["solverSuccess"])
    lock_rows = [row for row in rows if row["lockMode"] in {"right_state_locked", "right_position_locked"}]
    missing_evidence = [
        "segmentedStepOrMeshExport",
        "externalMuJoCoRun",
        "benchCoordinates",
        "measuredPinHoleClearance",
        "measuredConnectorSlip",
        "measuredFriction",
        "measuredContactStiffness",
        "measuredJointAxes",
    ]
    internal_ready = (
        solver_success_count == len(rows)
        and all(bool(row["lockPass"]) for row in lock_rows)
        and bool(clearance_trend["looseHasLargerClearance"])
        and bool(backlash_trend["highBacklashReducesAlphaResponse"])
    )
    return {
        "schema": TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA,
        "model": "named-two-cell-reduced-physics-test-suite",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "cadLayout": cad_rad_cell_layout(config),
        "rows": rows,
        "cases": detailed_cases,
        "measurementColumns": [
            "caseId",
            "leftX",
            "leftY",
            "leftZ",
            "rightX",
            "rightY",
            "rightZ",
            "leftAlpha",
            "rightAlpha",
            "leftTheta",
            "rightTheta",
            "upperConnectorSlipMm",
            "middleConnectorSlipMm",
            "lowerConnectorSlipMm",
            "lockHeldObserved",
            "notes",
        ],
        "summary": {
            "status": "reduced-two-cell-suite-ready-needs-external-data"
            if internal_ready
            else "reduced-two-cell-suite-internal-review-needed",
            "internalSuiteReady": internal_ready,
            "physicalAccuracyValidated": False,
            "caseCount": len(rows),
            "solverSuccessCount": solver_success_count,
            "solverFailureCount": len(rows) - solver_success_count,
            "lockCaseCount": len(lock_rows),
            "lockPassCount": sum(1 for row in lock_rows if row["lockPass"]),
            "maxVerticalSlipMm": max((row["connectorMaxVerticalSlipMm"] for row in rows), default=0.0),
            "maxVerticalExcessMm": max((row["connectorMaxVerticalExcessMm"] for row in rows), default=0.0),
            "maxContactPenalty": max((row["connectorTotalContactPenalty"] for row in rows), default=0.0),
            "clearanceTrend": clearance_trend,
            "backlashTrend": backlash_trend,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimLabels": {
            "simulation": "reduced quasistatic/contact suite for two connected RAD cells",
            "cadGeometry": "A360-derived dimensional abstraction, not exact segmented CAD dynamics",
            "physicalAccuracy": "not validated until the listed external engine and bench measurements are attached",
        },
    }


def _default_fidelity_actuation_cases(controls: TwoCellBenchControls) -> list[dict[str, Any]]:
    alpha_drive = abs(controls.alpha_command) if abs(controls.alpha_command) > 0.05 else 0.35
    z_drive = abs(controls.z_command) if abs(controls.z_command) > 0.05 else 0.35
    return [
        {"actuationCase": "contract_lift", "alpha": -alpha_drive, "z": z_drive, "gravity": 0.025},
        {"actuationCase": "expand_pushdown", "alpha": alpha_drive, "z": -z_drive, "gravity": 0.025},
        {"actuationCase": "alpha_contract", "alpha": -alpha_drive, "z": 0.0, "gravity": 0.025},
        {"actuationCase": "alpha_expand", "alpha": alpha_drive, "z": 0.0, "gravity": 0.025},
        {"actuationCase": "z_lift", "alpha": 0.0, "z": z_drive, "gravity": 0.025},
        {"actuationCase": "z_pushdown", "alpha": 0.0, "z": -z_drive, "gravity": 0.025},
        {"actuationCase": "gravity_sag", "alpha": 0.0, "z": 0.0, "gravity": 0.5},
    ]


def _matrix_lock_pass(lock_mode: str, solve: dict[str, Any], tolerance: float) -> bool:
    constraints = solve["constraints"]
    if lock_mode == "right_position_locked":
        return bool(constraints["positionLocksHold"]) and float(constraints["rightFixtureError"] or 0.0) <= tolerance
    if lock_mode == "right_state_locked":
        return bool(constraints["stateLocksHold"]) and float(constraints["rightStateAlphaError"] or 0.0) <= tolerance
    return True


def _row_lookup(
    rows: list[dict[str, Any]],
    *,
    actuation_case: str,
    lock_mode: str,
    backlash: float,
    hole_radius: float,
    tolerance: float,
) -> dict[str, Any] | None:
    for row in rows:
        if (
            row["actuationCase"] == actuation_case
            and row["lockMode"] == lock_mode
            and abs(float(row["backlash"]) - float(backlash)) <= tolerance
            and abs(float(row["holeRadius"]) - float(hole_radius)) <= tolerance
        ):
            return row
    return None


def _free_play_utilization(slip_mm: float, clearance_mm: float) -> float | None:
    slip = abs(float(slip_mm))
    clearance = max(0.0, float(clearance_mm))
    if clearance <= 1e-12:
        return None
    return slip / clearance


def two_cell_physical_fidelity_matrix(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    backlash_values: list[float] | tuple[float, ...] | None = None,
    hole_radii: list[float] | tuple[float, ...] | None = None,
    lock_modes: list[str] | tuple[str, ...] | None = None,
    actuation_cases: list[dict[str, Any]] | tuple[dict[str, Any], ...] | None = None,
    tolerance: float = 1e-7,
) -> dict[str, Any]:
    """Dense reduced-physics matrix for two connected CAD-derived RAD cells."""

    config = replace(config or LatticeConfig(rows=1, cols=2), rows=1, cols=2)
    controls = _controls(controls)
    nominal_backlash = float(config.backlash)
    nominal_hole = float(config.hole_radius)
    loose_hole = max(nominal_hole, controls.hole_sweep_max)
    backlash_axis = _unique_floats(
        list(backlash_values or (0.0, nominal_backlash, max(0.25, 3.0 * nominal_backlash)))
    )
    hole_axis = _unique_floats(list(hole_radii or (config.pin_radius, nominal_hole, loose_hole)))
    modes = list(lock_modes or ("free", "right_state_locked", "right_position_locked", "left_free"))
    cases = list(actuation_cases or _default_fidelity_actuation_cases(controls))

    rows: list[dict[str, Any]] = []
    connector_rows: list[dict[str, Any]] = []
    index = 0
    for backlash in backlash_axis:
        for hole in hole_axis:
            case_config = _two_cell_case_config(config, backlash=backlash, hole_radius=hole)
            for lock_mode in modes:
                locked_base = _lock_mode_controls(controls, lock_mode)
                for actuation in cases:
                    case_controls = _lock_mode_controls(
                        _controls(
                            locked_base,
                            alpha_command=float(actuation["alpha"]),
                            z_command=float(actuation["z"]),
                        ),
                        lock_mode,
                    )
                    contact = two_cell_connector_contact_report(
                        case_config,
                        case_controls,
                        gravity_force=float(actuation.get("gravity", 0.025)),
                        tolerance=tolerance,
                    )
                    solve = contact["quasistatic"]
                    left, right = solve["cells"]
                    constraints = solve["constraints"]
                    summary = contact["summary"]
                    lock_pass = _matrix_lock_pass(lock_mode, solve, tolerance)
                    row = {
                        "index": index,
                        "actuationCase": str(actuation["actuationCase"]),
                        "lockMode": lock_mode,
                        "backlash": case_config.backlash,
                        "alphaCommand": case_controls.alpha_command,
                        "zCommand": case_controls.z_command,
                        "gravityForce": float(actuation.get("gravity", 0.025)),
                        "pinRadius": solve["dimensions"]["pinRadius"],
                        "holeRadius": solve["dimensions"]["holeRadius"],
                        "clearance": solve["dimensions"]["pinHoleClearance"],
                        "effectiveBacklash": solve["deadZones"]["effectiveBacklash"],
                        "pinHoleClearanceMm": contact["dimensions"]["pinHoleClearanceMm"],
                        "solverSuccess": bool(solve["solver"]["success"]),
                        "lockPass": lock_pass,
                        "positionLocksHold": bool(constraints["positionLocksHold"]),
                        "stateLocksHold": bool(constraints["stateLocksHold"]),
                        "rightLockedVerticalMotionAllowed": bool(
                            constraints["rightLockedVerticalMotionAllowed"]
                        ),
                        "leftX": left["center"]["x"],
                        "leftZ": left["center"]["z"],
                        "rightX": right["center"]["x"],
                        "rightZ": right["center"]["z"],
                        "leftAlpha": left["alpha"],
                        "rightAlpha": right["alpha"],
                        "leftAlphaDeltaFromInitial": left["alpha"] - case_config.initial_alpha,
                        "rightAlphaDeltaFromInitial": right["alpha"] - case_config.initial_alpha,
                        "leftTheta": left["theta"],
                        "rightTheta": right["theta"],
                        "verticalShear": solve["contact"]["verticalShear"],
                        "verticalExcess": solve["contact"]["verticalExcess"],
                        "axialExcess": solve["contact"]["axialExcess"],
                        "contactMode": solve["contact"]["mode"],
                        "connectorActiveContactCount": summary["activeContactCount"],
                        "connectorInsideClearanceCount": summary["insideClearanceCount"],
                        "connectorMaxLateralSlipMm": summary["maxLateralSlipMm"],
                        "connectorMaxVerticalSlipMm": summary["maxVerticalSlipMm"],
                        "connectorMaxLateralExcessMm": summary["maxLateralExcessMm"],
                        "connectorMaxVerticalExcessMm": summary["maxVerticalExcessMm"],
                        "connectorTotalContactPenalty": summary["totalContactPenalty"],
                        "totalEnergy": solve["energy"]["totalEnergy"],
                        "rightFixtureError": constraints["rightFixtureError"],
                        "rightStateAlphaError": constraints["rightStateAlphaError"],
                    }
                    rows.append(row)
                    for connector in contact["connectors"]:
                        connector_rows.append(
                            {
                                "matrixIndex": index,
                                "actuationCase": row["actuationCase"],
                                "lockMode": lock_mode,
                                "backlash": case_config.backlash,
                                "holeRadius": solve["dimensions"]["holeRadius"],
                                "connector": connector["connector"],
                                "leftSite": connector["leftSite"],
                                "rightSite": connector["rightSite"],
                                "pinHoleClearanceMm": contact["dimensions"]["pinHoleClearanceMm"],
                                "leftPositionMm": connector["leftPositionMm"],
                                "rightPositionMm": connector["rightPositionMm"],
                                "slipMm": connector["slipMm"],
                                "lateralSlipMm": connector["lateralSlipMm"],
                                "verticalSlipMm": connector["verticalSlipMm"],
                                "totalSlipMm": connector["totalSlipMm"],
                                "lateralExcessMm": connector["lateralExcessMm"],
                                "verticalExcessMm": connector["verticalExcessMm"],
                                "contactPenalty": connector["contactPenalty"],
                                "contactMode": connector["contactMode"],
                            }
                        )
                    index += 1

    nominal_backlash_for_lookup = min(backlash_axis, key=lambda value: abs(value - nominal_backlash))
    nominal_hole_for_lookup = min(hole_axis, key=lambda value: abs(value - nominal_hole))
    clearance_probe_rows = sorted(
        [
            row
            for row in rows
            if row["actuationCase"] == "contract_lift"
            and row["lockMode"] == "free"
            and abs(float(row["backlash"]) - nominal_backlash_for_lookup) <= tolerance
        ],
        key=lambda row: float(row["holeRadius"]),
    )
    backlash_probe_rows = sorted(
        [
            row
            for row in rows
            if row["actuationCase"] == "contract_lift"
            and row["lockMode"] == "left_free"
            and abs(float(row["holeRadius"]) - nominal_hole_for_lookup) <= tolerance
        ],
        key=lambda row: float(row["backlash"]),
    )
    backlash_left_z_abs = [abs(float(row["leftZ"])) for row in backlash_probe_rows]
    backlash_vertical_materially_flat = (
        max(backlash_left_z_abs, default=0.0) - min(backlash_left_z_abs, default=0.0) <= max(1e-5, 10.0 * tolerance)
    )
    clearance_reference_slip_mm = (
        float(clearance_probe_rows[0]["connectorMaxVerticalSlipMm"]) if clearance_probe_rows else 0.0
    )
    clearance_relief_penalty_mm2 = [
        0.5
        * max(
            0.0,
            clearance_reference_slip_mm - float(row["pinHoleClearanceMm"]),
        )
        ** 2
        for row in clearance_probe_rows
    ]
    clearance_trend = {
        "holeRadii": [row["holeRadius"] for row in clearance_probe_rows],
        "clearances": [row["clearance"] for row in clearance_probe_rows],
        "verticalSlipMm": [row["connectorMaxVerticalSlipMm"] for row in clearance_probe_rows],
        "verticalExcessMm": [row["connectorMaxVerticalExcessMm"] for row in clearance_probe_rows],
        "verticalFreePlayUtilization": [
            _free_play_utilization(row["connectorMaxVerticalSlipMm"], row["pinHoleClearanceMm"])
            for row in clearance_probe_rows
        ],
        "zeroClearanceContact": [
            float(row["pinHoleClearanceMm"]) <= 1e-12 and float(row["connectorMaxVerticalSlipMm"]) > tolerance
            for row in clearance_probe_rows
        ],
        "contactPenalty": [row["connectorTotalContactPenalty"] for row in clearance_probe_rows],
        "clearanceReliefReferenceSlipMm": clearance_reference_slip_mm,
        "clearanceReliefPenaltyMm2": clearance_relief_penalty_mm2,
        "clearanceNondecreasing": _nondecreasing(
            [float(row["clearance"]) for row in clearance_probe_rows],
            tolerance,
        ),
        "rawVerticalSlipNonincreasing": _nonincreasing(
            [float(row["connectorMaxVerticalSlipMm"]) for row in clearance_probe_rows],
            tolerance,
        ),
        "verticalExcessNonincreasing": _nonincreasing(
            [float(row["connectorMaxVerticalExcessMm"]) for row in clearance_probe_rows],
            tolerance,
        ),
        "penaltyNonincreasing": _nonincreasing(clearance_relief_penalty_mm2, tolerance),
    }
    backlash_trend = {
        "backlashes": [row["backlash"] for row in backlash_probe_rows],
        "leftAlphaDeltaAbs": [abs(float(row["leftAlphaDeltaFromInitial"])) for row in backlash_probe_rows],
        "leftZAbs": backlash_left_z_abs,
        "alphaResponseNonincreasing": _nonincreasing(
            [abs(float(row["leftAlphaDeltaFromInitial"])) for row in backlash_probe_rows],
            tolerance,
        ),
        "verticalResponseNonincreasing": False
        if backlash_vertical_materially_flat
        else _nonincreasing(backlash_left_z_abs, tolerance),
        "verticalResponseMateriallyFlat": backlash_vertical_materially_flat,
    }
    lock_rows = [row for row in rows if row["lockMode"] in {"right_state_locked", "right_position_locked"}]
    state_lock_rows = [row for row in rows if row["lockMode"] == "right_state_locked"]
    position_lock_rows = [row for row in rows if row["lockMode"] == "right_position_locked"]
    lock_checks = {
        "lockCaseCount": len(lock_rows),
        "lockPassCount": sum(1 for row in lock_rows if row["lockPass"]),
        "stateLockAlphaHold": all(
            bool(row["stateLocksHold"]) and abs(float(row["rightStateAlphaError"] or 0.0)) <= tolerance
            for row in state_lock_rows
        ),
        "stateLockVerticalMotionSamples": sum(
            1
            for row in state_lock_rows
            if abs(float(row["zCommand"])) > tolerance and abs(float(row["rightZ"])) > tolerance
        ),
        "positionLockFixtureHold": all(
            bool(row["positionLocksHold"]) and abs(float(row["rightFixtureError"] or 0.0)) <= tolerance
            for row in position_lock_rows
        ),
    }
    lift = _row_lookup(
        rows,
        actuation_case="z_lift",
        lock_mode="free",
        backlash=nominal_backlash_for_lookup,
        hole_radius=nominal_hole_for_lookup,
        tolerance=tolerance,
    )
    pushdown = _row_lookup(
        rows,
        actuation_case="z_pushdown",
        lock_mode="free",
        backlash=nominal_backlash_for_lookup,
        hole_radius=nominal_hole_for_lookup,
        tolerance=tolerance,
    )
    contract = _row_lookup(
        rows,
        actuation_case="alpha_contract",
        lock_mode="free",
        backlash=nominal_backlash_for_lookup,
        hole_radius=nominal_hole_for_lookup,
        tolerance=tolerance,
    )
    expand = _row_lookup(
        rows,
        actuation_case="alpha_expand",
        lock_mode="free",
        backlash=nominal_backlash_for_lookup,
        hole_radius=nominal_hole_for_lookup,
        tolerance=tolerance,
    )
    actuation_polarity = {
        "zLiftPositive": bool(lift) and float(lift["rightZ"]) > tolerance,
        "zPushdownNegative": bool(pushdown) and float(pushdown["rightZ"]) < -tolerance,
        "alphaContractBelowInitial": bool(contract) and float(contract["rightAlpha"]) < config.initial_alpha,
        "alphaExpandAboveInitial": bool(expand) and float(expand["rightAlpha"]) > config.initial_alpha,
    }
    solver_success_count = sum(1 for row in rows if row["solverSuccess"])
    missing_evidence = [
        "segmentedCADBodies",
        "externalRigidBodyContactRun",
        "filledConnectorMeasurements",
        "benchCoordinates",
        "measuredFriction",
        "measuredContactStiffness",
        "measuredJointAxes",
        "measuredActuatorForceDisplacement",
    ]
    internal_ready = (
        solver_success_count == len(rows)
        and all(bool(row["lockPass"]) for row in lock_rows)
        and bool(clearance_trend["clearanceNondecreasing"])
        and bool(clearance_trend["penaltyNonincreasing"])
        and bool(backlash_trend["alphaResponseNonincreasing"])
        and all(bool(value) for value in actuation_polarity.values())
    )
    return {
        "schema": TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA,
        "model": "cad-derived-two-cell-hole-backlash-lock-actuation-matrix",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "cadLayout": cad_rad_cell_layout(config),
        "axes": {
            "backlashValues": backlash_axis,
            "holeRadii": hole_axis,
            "lockModes": modes,
            "actuationCases": [case["actuationCase"] for case in cases],
        },
        "rows": rows,
        "connectorRows": connector_rows,
        "summary": {
            "status": "reduced-physics-fidelity-matrix-ready-needs-external-data"
            if internal_ready
            else "reduced-physics-fidelity-matrix-internal-review-needed",
            "internalMatrixReady": internal_ready,
            "physicalAccuracyValidated": False,
            "rowCount": len(rows),
            "connectorRowCount": len(connector_rows),
            "solverSuccessCount": solver_success_count,
            "solverFailureCount": len(rows) - solver_success_count,
            "maxConnectorVerticalSlipMm": max(
                (float(row["connectorMaxVerticalSlipMm"]) for row in rows),
                default=0.0,
            ),
            "maxConnectorVerticalExcessMm": max(
                (float(row["connectorMaxVerticalExcessMm"]) for row in rows),
                default=0.0,
            ),
            "maxContactPenalty": max((float(row["connectorTotalContactPenalty"]) for row in rows), default=0.0),
            "clearanceTrend": clearance_trend,
            "backlashTrend": backlash_trend,
            "lockChecks": lock_checks,
            "actuationPolarity": actuation_polarity,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimLabels": {
            "simulation": "dense internal reduced-physics sweep for two connected RAD cells",
            "holeRadius": "hole radius changes pin-hole clearance and therefore backlash/free-play",
            "locks": "state locks freeze alpha; position locks freeze fixture coordinates",
            "physicalAccuracy": "not validated until external contact or bench connector measurements agree",
        },
    }


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _group_key(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.9g}"
    return str(value)


def _summarize_row_group(rows: list[dict[str, Any]]) -> dict[str, Any]:
    right_z = [float(row["rightZ"]) for row in rows]
    right_alpha_delta = [float(row["rightAlphaDeltaFromInitial"]) for row in rows]
    left_z = [float(row["leftZ"]) for row in rows]
    vertical_slip = [float(row["connectorMaxVerticalSlipMm"]) for row in rows]
    contact_penalty = [float(row["connectorTotalContactPenalty"]) for row in rows]
    lock_rows = [row for row in rows if row["lockMode"] in {"right_state_locked", "right_position_locked"}]
    return {
        "caseCount": len(rows),
        "solverSuccessCount": sum(1 for row in rows if bool(row["solverSuccess"])),
        "maxAbsRightZ": max((abs(value) for value in right_z), default=0.0),
        "meanAbsRightZ": _mean([abs(value) for value in right_z]),
        "maxAbsLeftZ": max((abs(value) for value in left_z), default=0.0),
        "maxAbsRightAlphaDelta": max((abs(value) for value in right_alpha_delta), default=0.0),
        "meanAbsRightAlphaDelta": _mean([abs(value) for value in right_alpha_delta]),
        "maxConnectorVerticalSlipMm": max(vertical_slip, default=0.0),
        "meanConnectorVerticalSlipMm": _mean(vertical_slip),
        "maxContactPenalty": max(contact_penalty, default=0.0),
        "meanContactPenalty": _mean(contact_penalty),
        "lockCaseCount": len(lock_rows),
        "lockPassCount": sum(1 for row in lock_rows if bool(row["lockPass"])),
    }


def _group_summary(rows: list[dict[str, Any]], field: str) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(_group_key(row[field]), []).append(row)
    out = []
    for key, group in grouped.items():
        out.append({field: group[0][field], **_summarize_row_group(group)})
    return sorted(out, key=lambda row: str(row[field]))


def _case_digest(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "index",
        "actuationCase",
        "lockMode",
        "backlash",
        "holeRadius",
        "pinHoleClearanceMm",
        "alphaCommand",
        "zCommand",
        "rightZ",
        "leftZ",
        "rightAlpha",
        "rightAlphaDeltaFromInitial",
        "connectorMaxVerticalSlipMm",
        "connectorMaxVerticalExcessMm",
        "connectorTotalContactPenalty",
        "totalEnergy",
        "lockPass",
        "positionLocksHold",
        "stateLocksHold",
        "rightLockedVerticalMotionAllowed",
    ]
    digest = {key: row[key] for key in keys}
    digest["caseId"] = _fidelity_matrix_case_id(row)
    digest["matrixIndex"] = int(row["index"])
    return digest


def two_cell_physical_response_atlas(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    rank_limit: int = 12,
    **matrix_overrides: Any,
) -> dict[str, Any]:
    """Condense the two-cell matrix into bench-facing response and lock evidence."""

    matrix = two_cell_physical_fidelity_matrix(config, controls, **matrix_overrides)
    rows = list(matrix["rows"])
    connector_rows = list(matrix["connectorRows"])
    limit = max(1, int(rank_limit))
    lock_failure_rows = [
        row
        for row in rows
        if row["lockMode"] in {"right_state_locked", "right_position_locked"} and not bool(row["lockPass"])
    ]
    state_lock_vertical_rows = [
        row
        for row in rows
        if row["lockMode"] == "right_state_locked"
        and abs(float(row["zCommand"])) > 1e-12
        and bool(row["rightLockedVerticalMotionAllowed"])
    ]
    position_lock_rows = [row for row in rows if row["lockMode"] == "right_position_locked"]
    ranked = {
        "maxVerticalSlip": [
            _case_digest(row)
            for row in sorted(rows, key=lambda row: float(row["connectorMaxVerticalSlipMm"]), reverse=True)[:limit]
        ],
        "maxContactPenalty": [
            _case_digest(row)
            for row in sorted(rows, key=lambda row: float(row["connectorTotalContactPenalty"]), reverse=True)[:limit]
        ],
        "maxRightHeight": [
            _case_digest(row)
            for row in sorted(rows, key=lambda row: abs(float(row["rightZ"])), reverse=True)[:limit]
        ],
        "lockFailures": [_case_digest(row) for row in lock_failure_rows[:limit]],
        "stateLockVerticalMotion": [_case_digest(row) for row in state_lock_vertical_rows[:limit]],
        "positionLockFixtureCases": [_case_digest(row) for row in position_lock_rows[:limit]],
    }
    summary = matrix["summary"]
    invariants = {
        "solverAllConverged": summary["solverSuccessCount"] == summary["rowCount"],
        "lockAllPass": summary["lockChecks"]["lockPassCount"] == summary["lockChecks"]["lockCaseCount"],
        "stateLockAlphaHolds": bool(summary["lockChecks"]["stateLockAlphaHold"]),
        "stateLockAllowsVerticalSamples": int(summary["lockChecks"]["stateLockVerticalMotionSamples"]),
        "positionLockFixtureHolds": bool(summary["lockChecks"]["positionLockFixtureHold"]),
        "clearanceNondecreasing": bool(summary["clearanceTrend"]["clearanceNondecreasing"]),
        "clearancePenaltyNonincreasing": bool(summary["clearanceTrend"]["penaltyNonincreasing"]),
        "backlashAlphaResponseNonincreasing": bool(summary["backlashTrend"]["alphaResponseNonincreasing"]),
        "zLiftPositive": bool(summary["actuationPolarity"]["zLiftPositive"]),
        "zPushdownNegative": bool(summary["actuationPolarity"]["zPushdownNegative"]),
        "alphaContractBelowInitial": bool(summary["actuationPolarity"]["alphaContractBelowInitial"]),
        "alphaExpandAboveInitial": bool(summary["actuationPolarity"]["alphaExpandAboveInitial"]),
        "physicalAccuracyValidated": False,
    }
    return {
        "schema": TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA,
        "model": "two-cell-response-atlas-from-cad-proxy-quasistatic-contact-matrix",
        "cadReference": matrix["cadReference"],
        "cadLayout": matrix["cadLayout"],
        "axes": matrix["axes"],
        "matrixSummary": summary,
        "groupSummaries": {
            "byLockMode": _group_summary(rows, "lockMode"),
            "byHoleRadius": _group_summary(rows, "holeRadius"),
            "byBacklash": _group_summary(rows, "backlash"),
            "byActuationCase": _group_summary(rows, "actuationCase"),
        },
        "rankedCases": ranked,
        "invariants": invariants,
        "connectorSummary": {
            "connectorRowCount": len(connector_rows),
            "contactModeCounts": dict(
                Counter(str(row["contactMode"]) for row in connector_rows)
            ),
            "maxVerticalSlipMm": max((float(row["verticalSlipMm"]) for row in connector_rows), default=0.0),
            "maxTotalSlipMm": max((float(row["totalSlipMm"]) for row in connector_rows), default=0.0),
            "maxContactPenalty": max((float(row["contactPenalty"]) for row in connector_rows), default=0.0),
        },
        "benchPriority": {
            "firstCasesToMeasure": ranked["maxVerticalSlip"][: min(5, len(ranked["maxVerticalSlip"]))],
            "lockCasesToMeasure": ranked["stateLockVerticalMotion"][: min(5, len(ranked["stateLockVerticalMotion"]))],
            "clearanceCasesToMeasure": [
                _case_digest(row)
                for row in sorted(
                    [
                        row
                        for row in rows
                        if row["actuationCase"] == "contract_lift" and row["lockMode"] == "free"
                    ],
                    key=lambda row: (float(row["backlash"]), float(row["holeRadius"])),
                )[:limit]
            ],
        },
        "claimBoundary": {
            "physicalAccuracyValidated": False,
            "allowedClaim": "internal reduced two-cell response atlas for selecting physical tests",
            "blockedClaim": "exact real-cell mechanics until segmented CAD contact and bench-coordinate holdouts agree",
            "remainingEvidence": summary["missingEvidence"],
        },
    }


def _phase_map_value(row: dict[str, Any], key: str, default: float = 0.0) -> float:
    value = row.get(key, default)
    if value is None:
        return default
    return float(value)


def _two_cell_contact_phase(row: dict[str, Any], tolerance: float) -> str:
    lock_mode = str(row.get("lockMode", "free"))
    contact_mode = str(row.get("contactMode", "inside-clearance"))
    vertical_excess = max(
        abs(_phase_map_value(row, "verticalExcess")),
        abs(_phase_map_value(row, "connectorMaxVerticalExcessMm")),
    )
    axial_excess = max(
        abs(_phase_map_value(row, "axialExcess")),
        abs(_phase_map_value(row, "connectorMaxLateralExcessMm")),
    )
    active_contact = (
        int(row.get("connectorActiveContactCount", 0) or 0) > 0
        or _phase_map_value(row, "connectorTotalContactPenalty") > tolerance
    )
    vertical_contact = vertical_excess > tolerance or "pin-hole-vertical-contact" in contact_mode
    axial_contact = axial_excess > tolerance or "axial-backlash-contact" in contact_mode
    gravity_sag = (
        _phase_map_value(row, "gravityForce") > tolerance
        and abs(_phase_map_value(row, "zCommand")) <= tolerance
        and min(_phase_map_value(row, "leftZ"), _phase_map_value(row, "rightZ")) < -tolerance
    )
    if lock_mode == "right_position_locked":
        return "position-locked"
    if lock_mode == "right_state_locked":
        return (
            "state-locked-vertical-free"
            if bool(row.get("rightLockedVerticalMotionAllowed"))
            and abs(_phase_map_value(row, "rightZ")) > tolerance
            else "state-locked"
        )
    if axial_contact and vertical_contact:
        return "axial-and-vertical-contact"
    if vertical_contact:
        return "vertical-contact"
    if axial_contact or active_contact:
        return "axial-contact"
    if gravity_sag:
        return "gravity-sag"
    return "free-play"


def _phase_map_case_digest(row: dict[str, Any], category: str) -> dict[str, Any]:
    return {
        "caseId": _fidelity_matrix_case_id(row),
        "category": category,
        "matrixIndex": int(row["index"]),
        "phase": row["phase"],
        "actuationCase": row["actuationCase"],
        "lockMode": row["lockMode"],
        "backlash": row["backlash"],
        "holeRadius": row["holeRadius"],
        "pinHoleClearanceMm": row["pinHoleClearanceMm"],
        "effectiveBacklash": row.get("effectiveBacklash", row["backlash"]),
        "effectiveHoleRadius": row.get("effectiveHoleRadius", row["holeRadius"]),
        "effectivePinHoleClearanceMm": row.get(
            "effectivePinHoleClearanceMm",
            row["pinHoleClearanceMm"],
        ),
        "alphaCommand": row["alphaCommand"],
        "zCommand": row["zCommand"],
        "gravityForce": row["gravityForce"],
        "leftZ": row["leftZ"],
        "rightZ": row["rightZ"],
        "rightAlpha": row["rightAlpha"],
        "contactMode": row["contactMode"],
        "connectorActiveContactCount": row["connectorActiveContactCount"],
        "connectorMaxVerticalSlipMm": row["connectorMaxVerticalSlipMm"],
        "connectorMaxVerticalExcessMm": row["connectorMaxVerticalExcessMm"],
        "connectorTotalContactPenalty": row["connectorTotalContactPenalty"],
        "lockPass": row["lockPass"],
        "physicalAccuracyValidated": False,
    }


def _phase_priority_cases(rows: list[dict[str, Any]], limit: int) -> list[dict[str, Any]]:
    categories = [
        (
            "vertical-contact",
            lambda row: row["phase"] in {"vertical-contact", "axial-and-vertical-contact"},
            lambda row: (
                abs(float(row["connectorMaxVerticalExcessMm"])),
                abs(float(row["connectorMaxVerticalSlipMm"])),
                abs(float(row["rightZ"])),
            ),
        ),
        (
            "axial-contact",
            lambda row: row["phase"] in {"axial-contact", "axial-and-vertical-contact"},
            lambda row: (
                abs(float(row.get("connectorMaxLateralExcessMm", 0.0))),
                abs(float(row.get("axialExcess", 0.0))),
                abs(float(row["connectorTotalContactPenalty"])),
            ),
        ),
        (
            "gravity-sag",
            lambda row: row["phase"] == "gravity-sag",
            lambda row: (abs(float(row["leftZ"])), abs(float(row["rightZ"]))),
        ),
        (
            "state-lock-vertical-free",
            lambda row: row["phase"] == "state-locked-vertical-free",
            lambda row: (abs(float(row["rightZ"])), abs(float(row["zCommand"]))),
        ),
        (
            "position-lock",
            lambda row: row["phase"] == "position-locked",
            lambda row: (
                abs(float(row["connectorMaxVerticalSlipMm"])),
                abs(float(row["connectorTotalContactPenalty"])),
            ),
        ),
        (
            "free-play-boundary",
            lambda row: row["phase"] == "free-play",
            lambda row: (abs(float(row["connectorMaxVerticalSlipMm"])), abs(float(row["leftAlphaDeltaFromInitial"]))),
        ),
    ]
    out: list[dict[str, Any]] = []
    seen: set[int] = set()
    for category, predicate, sort_key in categories:
        matches = [row for row in rows if predicate(row)]
        matches.sort(key=sort_key, reverse=True)
        for row in matches[: max(1, min(3, limit))]:
            index = int(row["index"])
            if index in seen:
                continue
            seen.add(index)
            out.append(_phase_map_case_digest(row, category))
            if len(out) >= limit:
                return out
    return out


def two_cell_contact_phase_map(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    priority_limit: int = 12,
    tolerance: float = 1e-7,
    **matrix_overrides: Any,
) -> dict[str, Any]:
    """Classify two-cell matrix rows into free-play, contact, and lock phases."""

    matrix = two_cell_physical_fidelity_matrix(
        config,
        controls,
        tolerance=tolerance,
        **matrix_overrides,
    )
    rows: list[dict[str, Any]] = []
    for source in matrix["rows"]:
        phase = _two_cell_contact_phase(source, tolerance)
        row = {
            "caseId": _fidelity_matrix_case_id(source),
            **source,
            "phase": phase,
            "axialContact": phase in {"axial-contact", "axial-and-vertical-contact"},
            "verticalContact": phase in {"vertical-contact", "axial-and-vertical-contact"},
            "lockPhase": phase in {"position-locked", "state-locked", "state-locked-vertical-free"},
            "freePlayPhase": phase == "free-play",
            "gravitySagPhase": phase == "gravity-sag",
            "physicalAccuracyValidated": False,
        }
        rows.append(row)

    phase_counts = dict(Counter(str(row["phase"]) for row in rows))
    dominant_phase = max(sorted(phase_counts), key=lambda phase: phase_counts[phase]) if phase_counts else "none"
    active_contact_count = sum(1 for row in rows if int(row["connectorActiveContactCount"]) > 0)
    lock_rows = [row for row in rows if bool(row["lockPhase"])]
    priority_limit = max(1, int(priority_limit))
    return {
        "schema": TWO_CELL_CONTACT_PHASE_MAP_SCHEMA,
        "model": "two-cell-contact-lock-phase-map-from-fidelity-matrix",
        "cadReference": matrix["cadReference"],
        "cadLayout": matrix["cadLayout"],
        "axes": matrix["axes"],
        "matrixSummary": matrix["summary"],
        "rows": rows,
        "summary": {
            "status": "contact-phase-map-ready-needs-bench-or-external-validation"
            if matrix["summary"]["solverFailureCount"] == 0
            else "contact-phase-map-review-solver-failures",
            "rowCount": len(rows),
            "phaseCounts": phase_counts,
            "dominantPhase": dominant_phase,
            "solverFailureCount": matrix["summary"]["solverFailureCount"],
            "lockCaseCount": len(lock_rows),
            "lockFailureCount": sum(1 for row in lock_rows if not bool(row["lockPass"])),
            "activeContactCaseCount": active_contact_count,
            "freePlayCaseCount": phase_counts.get("free-play", 0),
            "gravitySagCaseCount": phase_counts.get("gravity-sag", 0),
            "maxVerticalSlipMm": max((abs(float(row["connectorMaxVerticalSlipMm"])) for row in rows), default=0.0),
            "maxVerticalExcessMm": max((abs(float(row["connectorMaxVerticalExcessMm"])) for row in rows), default=0.0),
            "maxAxialExcess": max((abs(float(row["axialExcess"])) for row in rows), default=0.0),
            "maxContactPenalty": max((abs(float(row["connectorTotalContactPenalty"])) for row in rows), default=0.0),
            "measurementPriority": _phase_priority_cases(rows, priority_limit),
            "physicalAccuracyValidated": False,
            "remainingEvidence": matrix["summary"]["missingEvidence"],
        },
        "claimBoundary": {
            "allowedClaim": "interpretable reduced-model phase classification for selecting two-cell tests",
            "blockedClaim": "exact real-cell phase boundaries without segmented CAD contact and bench validation",
            "physicalAccuracyValidated": False,
        },
    }


def _axis_values(start: float, stop: float, steps: int) -> list[float]:
    count = max(2, int(steps))
    if math.isclose(float(start), float(stop), rel_tol=0.0, abs_tol=1e-12):
        return [float(start)]
    return [float(start) + (float(stop) - float(start)) * index / (count - 1) for index in range(count)]


def _fidelity_actuation_case(controls: TwoCellBenchControls, actuation_case: str) -> dict[str, Any]:
    cases = _default_fidelity_actuation_cases(controls)
    for case in cases:
        if str(case["actuationCase"]) == str(actuation_case):
            return dict(case)
    known = ", ".join(str(case["actuationCase"]) for case in cases)
    raise ValueError(f"unknown actuation_case {actuation_case!r}; expected one of {known}")


def two_cell_radius_backlash_phase_diagram(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    actuation_case: str = "contract_lift",
    lock_mode: str = "free",
    hole_min: float | None = None,
    hole_max: float | None = None,
    hole_steps: int = 9,
    backlash_min: float = 0.0,
    backlash_max: float | None = None,
    backlash_steps: int = 9,
    hole_radius_bias: float = 0.0,
    backlash_bias: float = 0.0,
    tolerance: float = 1e-7,
) -> dict[str, Any]:
    """Dense two-cell phase diagram over physical hole radius and backlash."""

    config = replace(config or LatticeConfig(rows=1, cols=2), rows=1, cols=2)
    controls = _controls(controls)
    h_min = float(config.pin_radius if hole_min is None else hole_min)
    h_max = float(max(config.hole_radius, controls.hole_sweep_max) if hole_max is None else hole_max)
    b_min = float(backlash_min)
    b_max = float(max(0.25, 3.0 * config.backlash) if backlash_max is None else backlash_max)
    hole_axis = _unique_floats(_axis_values(h_min, h_max, hole_steps))
    backlash_axis = _unique_floats(_axis_values(b_min, b_max, backlash_steps))
    h_bias = float(hole_radius_bias)
    b_bias = float(backlash_bias)
    effective_hole_axis = [max(float(config.pin_radius), float(hole) + h_bias) for hole in hole_axis]
    effective_backlash_axis = [max(0.0, float(backlash) + b_bias) for backlash in backlash_axis]
    actuation = _fidelity_actuation_case(controls, actuation_case)
    matrix = two_cell_physical_fidelity_matrix(
        config,
        controls,
        backlash_values=effective_backlash_axis,
        hole_radii=effective_hole_axis,
        lock_modes=(lock_mode,),
        actuation_cases=(actuation,),
        tolerance=tolerance,
    )
    rows: list[dict[str, Any]] = []
    grid: list[list[str]] = []
    by_coord = {
        (float(row["backlash"]), float(row["holeRadius"])): row
        for row in matrix["rows"]
    }
    for row_index, backlash in enumerate(backlash_axis):
        grid_row: list[str] = []
        effective_backlash = effective_backlash_axis[row_index]
        for column_index, hole in enumerate(hole_axis):
            effective_hole = effective_hole_axis[column_index]
            source = by_coord[(float(effective_backlash), float(effective_hole))]
            phase = _two_cell_contact_phase(source, tolerance)
            row = {
                "caseId": f"PD_{len(rows):03d}",
                "matrixIndex": int(source["index"]),
                **source,
                "backlash": backlash,
                "holeRadius": hole,
                "effectiveBacklash": source["backlash"],
                "effectiveHoleRadius": source["holeRadius"],
                "effectivePinHoleClearanceMm": source["pinHoleClearanceMm"],
                "phase": phase,
                "axialContact": phase in {"axial-contact", "axial-and-vertical-contact"},
                "verticalContact": phase in {"vertical-contact", "axial-and-vertical-contact"},
                "lockPhase": phase in {"position-locked", "state-locked", "state-locked-vertical-free"},
                "freePlayPhase": phase == "free-play",
                "gravitySagPhase": phase == "gravity-sag",
                "physicalAccuracyValidated": False,
            }
            rows.append(row)
            grid_row.append(phase)
        grid.append(grid_row)

    phase_counts = dict(Counter(row["phase"] for row in rows))
    transition_rows = []
    for row_index, backlash in enumerate(backlash_axis):
        phases = grid[row_index]
        for column_index in range(1, len(phases)):
            if phases[column_index] != phases[column_index - 1]:
                transition_rows.append(
                    {
                        "axis": "holeRadius",
                        "backlash": backlash,
                        "leftHoleRadius": hole_axis[column_index - 1],
                        "rightHoleRadius": hole_axis[column_index],
                        "fromPhase": phases[column_index - 1],
                        "toPhase": phases[column_index],
                    }
                )
    for column_index, hole in enumerate(hole_axis):
        column_phases = [grid[row_index][column_index] for row_index in range(len(backlash_axis))]
        for row_index in range(1, len(column_phases)):
            if column_phases[row_index] != column_phases[row_index - 1]:
                transition_rows.append(
                    {
                        "axis": "backlash",
                        "holeRadius": hole,
                        "lowerBacklash": backlash_axis[row_index - 1],
                        "upperBacklash": backlash_axis[row_index],
                        "fromPhase": column_phases[row_index - 1],
                        "toPhase": column_phases[row_index],
                    }
                )
    dominant_phase = max(sorted(phase_counts), key=lambda phase: phase_counts[phase]) if phase_counts else "none"
    return {
        "schema": TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA,
        "model": "dense-two-cell-hole-radius-backlash-contact-phase-diagram",
        "cadReference": matrix["cadReference"],
        "cadLayout": matrix["cadLayout"],
        "axes": {
            "holeRadii": hole_axis,
            "backlashValues": backlash_axis,
            "effectiveHoleRadii": effective_hole_axis,
            "effectiveBacklashValues": effective_backlash_axis,
            "holeRadiusBias": h_bias,
            "backlashBias": b_bias,
            "actuationCase": str(actuation_case),
            "lockMode": str(lock_mode),
        },
        "actuation": {
            "actuationCase": actuation["actuationCase"],
            "alphaCommand": actuation["alpha"],
            "zCommand": actuation["z"],
            "gravityForce": actuation["gravity"],
        },
        "phaseGrid": grid,
        "rows": rows,
        "transitions": transition_rows,
        "matrixSummary": matrix["summary"],
        "summary": {
            "status": "radius-backlash-phase-diagram-ready-needs-bench-or-external-validation",
            "rowCount": len(rows),
            "holeStepCount": len(hole_axis),
            "backlashStepCount": len(backlash_axis),
            "phaseCounts": phase_counts,
            "dominantPhase": dominant_phase,
            "transitionCount": len(transition_rows),
            "activeContactCaseCount": sum(1 for row in rows if int(row["connectorActiveContactCount"]) > 0),
            "maxVerticalSlipMm": max((abs(float(row["connectorMaxVerticalSlipMm"])) for row in rows), default=0.0),
            "maxVerticalExcessMm": max((abs(float(row["connectorMaxVerticalExcessMm"])) for row in rows), default=0.0),
            "maxContactPenalty": max((abs(float(row["connectorTotalContactPenalty"])) for row in rows), default=0.0),
            "physicalAccuracyValidated": False,
            "remainingEvidence": matrix["summary"]["missingEvidence"],
        },
        "claimBoundary": {
            "allowedClaim": "dense reduced-model phase diagram for choosing hole-radius and backlash experiments",
            "blockedClaim": "real transition boundaries until segmented CAD contact and bench sweeps are compared",
            "physicalAccuracyValidated": False,
        },
    }


def _phase_diagram_coord_key(backlash: float, hole_radius: float) -> tuple[float, float]:
    return (round(float(backlash), 12), round(float(hole_radius), 12))


def _phase_diagram_case_for_transition(
    rows_by_coord: dict[tuple[float, float], dict[str, Any]],
    *,
    backlash: float,
    hole_radius: float,
) -> dict[str, Any] | None:
    return rows_by_coord.get(_phase_diagram_coord_key(backlash, hole_radius))


def _transition_case_digest(row: dict[str, Any], *, category: str, rank: int) -> dict[str, Any]:
    digest = _phase_map_case_digest(row, category)
    digest.update(
        {
            "rank": rank,
            "phaseDiagramCaseId": row.get("caseId", ""),
            "measurementTarget": "two-cell center pose plus upper/middle/lower connector marker coordinates",
            "observedPhase": "",
            "observedRightZ": "",
            "observedRightAlpha": "",
            "observedConnectorMaxVerticalSlipMm": "",
            "observedConnectorMaxVerticalExcessMm": "",
            "observedLockHeld": "",
            "notes": "",
        }
    )
    return digest


def two_cell_radius_backlash_transition_report(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    priority_limit: int = 16,
    **diagram_overrides: Any,
) -> dict[str, Any]:
    """Convert a dense phase diagram into bench-measurable transition brackets."""

    diagram = two_cell_radius_backlash_phase_diagram(config, controls, **diagram_overrides)
    rows_by_coord = {
        _phase_diagram_coord_key(row["backlash"], row["holeRadius"]): row
        for row in diagram["rows"]
    }
    brackets: list[dict[str, Any]] = []
    candidate_rows: dict[str, tuple[str, dict[str, Any]]] = {}
    for index, transition in enumerate(diagram["transitions"], start=1):
        if transition["axis"] == "holeRadius":
            left = _phase_diagram_case_for_transition(
                rows_by_coord,
                backlash=transition["backlash"],
                hole_radius=transition["leftHoleRadius"],
            )
            right = _phase_diagram_case_for_transition(
                rows_by_coord,
                backlash=transition["backlash"],
                hole_radius=transition["rightHoleRadius"],
            )
            span = abs(float(transition["rightHoleRadius"]) - float(transition["leftHoleRadius"]))
            midpoint = 0.5 * (float(transition["rightHoleRadius"]) + float(transition["leftHoleRadius"]))
            bracket = {
                "transitionId": f"TR_{index:03d}",
                "axis": "holeRadius",
                "fixedBacklash": transition["backlash"],
                "lowerHoleRadius": transition["leftHoleRadius"],
                "upperHoleRadius": transition["rightHoleRadius"],
                "midHoleRadius": midpoint,
                "bracketWidth": span,
                "fromPhase": transition["fromPhase"],
                "toPhase": transition["toPhase"],
                "leftCaseId": left.get("caseId", "") if left else "",
                "rightCaseId": right.get("caseId", "") if right else "",
                "leftMatrixCaseId": _fidelity_matrix_case_id(left) if left else "",
                "rightMatrixCaseId": _fidelity_matrix_case_id(right) if right else "",
                "physicalAccuracyValidated": False,
            }
        else:
            left = _phase_diagram_case_for_transition(
                rows_by_coord,
                backlash=transition["lowerBacklash"],
                hole_radius=transition["holeRadius"],
            )
            right = _phase_diagram_case_for_transition(
                rows_by_coord,
                backlash=transition["upperBacklash"],
                hole_radius=transition["holeRadius"],
            )
            span = abs(float(transition["upperBacklash"]) - float(transition["lowerBacklash"]))
            midpoint = 0.5 * (float(transition["upperBacklash"]) + float(transition["lowerBacklash"]))
            bracket = {
                "transitionId": f"TR_{index:03d}",
                "axis": "backlash",
                "fixedHoleRadius": transition["holeRadius"],
                "lowerBacklash": transition["lowerBacklash"],
                "upperBacklash": transition["upperBacklash"],
                "midBacklash": midpoint,
                "bracketWidth": span,
                "fromPhase": transition["fromPhase"],
                "toPhase": transition["toPhase"],
                "leftCaseId": left.get("caseId", "") if left else "",
                "rightCaseId": right.get("caseId", "") if right else "",
                "leftMatrixCaseId": _fidelity_matrix_case_id(left) if left else "",
                "rightMatrixCaseId": _fidelity_matrix_case_id(right) if right else "",
                "physicalAccuracyValidated": False,
            }
        if left:
            candidate_rows[str(left["caseId"])] = ("transition-lower-side", left)
        if right:
            candidate_rows[str(right["caseId"])] = ("transition-upper-side", right)
        brackets.append(bracket)

    ranked_by_slip = sorted(
        diagram["rows"],
        key=lambda row: abs(float(row.get("connectorMaxVerticalSlipMm", 0.0))),
        reverse=True,
    )
    ranked_by_penalty = sorted(
        diagram["rows"],
        key=lambda row: abs(float(row.get("connectorTotalContactPenalty", 0.0))),
        reverse=True,
    )
    if ranked_by_slip:
        candidate_rows.setdefault(str(ranked_by_slip[0]["caseId"]), ("max-vertical-slip", ranked_by_slip[0]))
    if ranked_by_penalty:
        candidate_rows.setdefault(str(ranked_by_penalty[0]["caseId"]), ("max-contact-penalty", ranked_by_penalty[0]))
    if not candidate_rows and diagram["rows"]:
        candidate_rows[str(diagram["rows"][0]["caseId"])] = ("no-transition-baseline", diagram["rows"][0])

    priority_limit = max(1, int(priority_limit))
    measurement_cases = [
        _transition_case_digest(row, category=category, rank=rank)
        for rank, (category, row) in enumerate(candidate_rows.values(), start=1)
    ][:priority_limit]
    bracket_widths = [float(item["bracketWidth"]) for item in brackets]
    phase_pairs = dict(Counter(f"{item['fromPhase']}->{item['toPhase']}" for item in brackets))
    return {
        "schema": TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA,
        "model": "two-cell-radius-backlash-transition-bracket-calibration-report",
        "diagram": {
            "schema": diagram["schema"],
            "axes": diagram["axes"],
            "summary": diagram["summary"],
        },
        "cadReference": diagram["cadReference"],
        "transitionBrackets": brackets,
        "measurementCases": measurement_cases,
        "measurementProtocol": {
            "goal": "Measure the two-cell cases adjacent to predicted phase changes and fit the real pin-hole backlash boundary.",
            "requiredObservables": [
                "left and right cell center xyz",
                "left and right cell alpha/theta",
                "upper, middle, and lower connector marker xyz",
                "observed contact/free-play phase",
                "whether state locks and position locks held",
            ],
            "notes": [
                "Use the same pin radius, hole radius, and backlash labels recorded in each case.",
                "Do not treat the reduced-model bracket midpoint as a real threshold until measured.",
                "Repeat transition-adjacent cases after any CAD export or hardware change.",
            ],
        },
        "summary": {
            "status": "transition-report-ready-needs-bench-sweep-validation",
            "transitionBracketCount": len(brackets),
            "measurementCaseCount": len(measurement_cases),
            "phasePairCounts": phase_pairs,
            "minBracketWidth": min(bracket_widths) if bracket_widths else 0.0,
            "maxBracketWidth": max(bracket_widths) if bracket_widths else 0.0,
            "physicalAccuracyValidated": False,
            "remainingEvidence": [
                "measured transition-adjacent two-cell coordinates",
                "segmented one-cell and two-cell CAD contact geometry",
                "pin-hole friction and normal stiffness calibration",
                "independent holdout transition sweeps",
            ],
        },
        "claimBoundary": {
            "allowedClaim": "reduced-model transition brackets and prioritized bench cases",
            "blockedClaim": "real RAD phase-transition law without measured transition sweeps and exact contact geometry",
            "physicalAccuracyValidated": False,
        },
    }


def sweep_two_cell_quasistatic(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    alpha_commands: list[float] | tuple[float, ...] | None = None,
    z_commands: list[float] | tuple[float, ...] | None = None,
    hole_radii: list[float] | tuple[float, ...] | None = None,
    lock_modes: list[str] | tuple[str, ...] | None = None,
    **solve_overrides: Any,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    radii = _radii(config)
    alpha_values = _unique_floats(list(alpha_commands or (-0.5, 0.0, 0.5)))
    z_values = _unique_floats(list(z_commands or (-0.35, 0.0, 0.35)))
    hole_values = _unique_floats(
        list(hole_radii or (radii["pinRadius"], radii["holeRadius"], controls.hole_sweep_max))
    )
    modes = list(lock_modes or ("free", "right_state_locked", "right_position_locked", "left_free"))
    rows = []
    index = 0
    for mode in modes:
        mode_base = _lock_mode_controls(controls, mode)
        for hole in hole_values:
            for alpha in alpha_values:
                for z in z_values:
                    case_controls = _controls(mode_base, alpha_command=alpha, z_command=z)
                    result = solve_two_cell_quasistatic(
                        config,
                        case_controls,
                        hole_radius=hole,
                        **solve_overrides,
                    )
                    left, right = result["cells"]
                    row = {
                        "index": index,
                        "lockMode": mode,
                        "alphaCommand": alpha,
                        "zCommand": z,
                        "pinRadius": result["dimensions"]["pinRadius"],
                        "holeRadius": result["dimensions"]["holeRadius"],
                        "clearance": result["dimensions"]["pinHoleClearance"],
                        "contactMode": result["contact"]["mode"],
                        "solverSuccess": result["solver"]["success"],
                        "leftX": left["center"]["x"],
                        "leftZ": left["center"]["z"],
                        "leftAlpha": left["alpha"],
                        "rightX": right["center"]["x"],
                        "rightZ": right["center"]["z"],
                        "rightAlpha": right["alpha"],
                        "verticalShear": result["contact"]["verticalShear"],
                        "verticalExcess": result["contact"]["verticalExcess"],
                        "axialExcess": result["contact"]["axialExcess"],
                        "storedEnergy": result["energy"]["storedEnergy"],
                        "totalEnergy": result["energy"]["totalEnergy"],
                        "leftFixtureError": result["constraints"]["leftFixtureError"],
                        "rightFixtureError": result["constraints"]["rightFixtureError"],
                        "rightStateAlphaError": result["constraints"]["rightStateAlphaError"],
                    }
                    rows.append(row)
                    index += 1
    clearances = [row["clearance"] for row in rows]
    right_z_values = [row["rightZ"] for row in rows]
    failures = [row for row in rows if not row["solverSuccess"]]
    return {
        "schema": TWO_CELL_QUASISTATIC_SWEEP_SCHEMA,
        "model": "two-cell-quasistatic-lock-actuation-clearance-sweep",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "grid": {
            "alphaCommands": alpha_values,
            "zCommands": z_values,
            "holeRadii": hole_values,
            "lockModes": modes,
        },
        "rows": rows,
        "summary": {
            "rowCount": len(rows),
            "solverSuccessCount": len(rows) - len(failures),
            "solverFailureCount": len(failures),
            "clearanceMin": min(clearances) if clearances else 0.0,
            "clearanceMax": max(clearances) if clearances else 0.0,
            "rightZMin": min(right_z_values) if right_z_values else 0.0,
            "rightZMax": max(right_z_values) if right_z_values else 0.0,
            "maxVerticalExcess": max((abs(row["verticalExcess"]) for row in rows), default=0.0),
            "maxAxialExcess": max((abs(row["axialExcess"]) for row in rows), default=0.0),
            "maxStoredEnergy": max((row["storedEnergy"] for row in rows), default=0.0),
            "requiresExternalEngineRun": True,
            "requiresMeasuredData": True,
        },
    }


def two_cell_physical_test_packet(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    cases = _named_case_controls(controls)
    results = [
        {"caseId": case_id, "bench": simulate_two_cell_bench(config, case_controls)}
        for case_id, case_controls in cases
    ]
    return {
        "schema": TWO_CELL_PACKET_SCHEMA,
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "cadReferenceProfile": cad_rad_cell_reference_profile(config),
        "cadArchiveAudit": cad_rad_cell_archive_audit(),
        "oneCellCadExportAudit": one_cell_cad_export_audit(config),
        "cadLayout": cad_rad_cell_layout(config),
        "config": {
            "cellSize": config.cell_size,
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "initialAlpha": config.initial_alpha,
        },
        "purpose": (
            "Bench packet for measuring whether the reduced RAD two-cell model "
            "matches real free, locked, position-locked, and clearance-swept pairs."
        ),
        "cases": results,
        "clearanceSweep": sweep_two_cell_backlash(config, controls),
        "actuationSweep": two_cell_actuation_sweep(config, controls),
        "quasistaticPhysics": solve_two_cell_quasistatic(config, controls),
        "quasistaticSweep": sweep_two_cell_quasistatic(config, controls),
        "physicalSimulationSuite": two_cell_physical_simulation_suite(config, controls),
        "physicalResponseAtlas": two_cell_physical_response_atlas(config, controls),
        "contactPhaseMap": two_cell_contact_phase_map(config, controls),
        "radiusBacklashPhaseDiagram": two_cell_radius_backlash_phase_diagram(config, controls),
        "radiusBacklashTransitionReport": two_cell_radius_backlash_transition_report(config, controls),
        "cadContactDecomposition": two_cell_cad_contact_decomposition_spec(config, controls),
        "exactContactHandoffPlan": two_cell_exact_contact_handoff_plan(config, controls),
        "connectorMeasurementTemplate": two_cell_connector_measurement_template(config, controls),
        "fidelityMatrixMeasurementTemplate": two_cell_fidelity_matrix_measurement_template(config, controls),
        "externalFidelityMatrixManifest": two_cell_external_fidelity_matrix_manifest(config, controls),
        "connectorContact": two_cell_connector_contact_report(config, controls),
        "connectorContactSweep": sweep_two_cell_connector_contact(config, controls),
        "segmentedCadReadiness": two_cell_segmented_cad_readiness_report(config, controls),
        "segmentedCadIntakeTemplates": two_cell_segmented_cad_intake_templates(config, controls),
        "segmentedCadIntakeValidation": two_cell_segmented_cad_intake_validation_report(config, controls),
        "mjcfProxy": two_cell_mjcf_proxy_report(config, controls),
        "mujocoProxyRun": two_cell_mujoco_proxy_run_report(config, controls),
        "externalResultsTemplate": two_cell_external_results_template(config, controls),
        "physicsValidationReport": two_cell_physics_validation_report(config, controls),
        "measurementTemplateColumns": [
            "caseId",
            "alphaCommand",
            "zCommand",
            "pinRadius",
            "holeRadius",
            "leftX",
            "leftY",
            "leftZ",
            "rightX",
            "rightY",
            "rightZ",
            "rightAlpha",
            "rightTheta",
            "contactModeObserved",
            "notes",
        ],
        "measurementTemplatePredictedColumns": [
            "predictedLeftX",
            "predictedLeftY",
            "predictedLeftZ",
            "predictedRightX",
            "predictedRightY",
            "predictedRightZ",
            "predictedRightAlpha",
            "predictedRightTheta",
            "contactModePredicted",
        ],
        "acceptance": {
            "requiresMeasuredData": True,
            "positionLockMustHoldFixtureCoordinates": True,
            "clearanceSweepExpectedTrend": "larger hole radius increases dead-zone and reduces transmitted neighbor response",
        },
    }


def two_cell_measurement_template_rows(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> list[dict[str, Any]]:
    packet = two_cell_physical_test_packet(config, controls)
    rows = []
    for case in packet["cases"]:
        bench = case["bench"]
        controls_dict = bench["controls"]
        left = bench["cells"][0]
        right = bench["cells"][1]
        dimensions = bench["dimensions"]
        rows.append(
            {
                "caseId": case["caseId"],
                "alphaCommand": controls_dict["alphaCommand"],
                "zCommand": controls_dict["zCommand"],
                "pinRadius": dimensions["pinRadius"],
                "holeRadius": dimensions["holeRadius"],
                "leftX": "",
                "leftY": "",
                "leftZ": "",
                "rightX": "",
                "rightY": "",
                "rightZ": "",
                "rightAlpha": "",
                "rightTheta": "",
                "contactModeObserved": "",
                "notes": "",
                "predictedLeftX": left["center"]["x"],
                "predictedLeftY": left["center"]["y"],
                "predictedLeftZ": left["center"]["z"],
                "predictedRightX": right["center"]["x"],
                "predictedRightY": right["center"]["y"],
                "predictedRightZ": right["center"]["z"],
                "predictedRightAlpha": right["alpha"],
                "predictedRightTheta": right["theta"],
                "contactModePredicted": bench["connector"]["contactMode"],
            }
        )
    return rows


def two_cell_connector_measurement_template_rows(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> list[dict[str, Any]]:
    suite = two_cell_physical_simulation_suite(config, controls)
    rows: list[dict[str, Any]] = []
    for case in suite["cases"]:
        case_id = str(case["caseId"])
        case_controls = case["controls"]
        case_config = case["config"]
        contact = case["connectorContact"]
        lock_mode = (
            "right_position_locked"
            if case_controls["rightPositionLocked"]
            else "right_state_locked"
            if case_controls["rightLocked"]
            else "left_free"
            if not case_controls["leftPositionLocked"]
            else "free"
        )
        for connector in contact["connectors"]:
            left_position = connector["leftPositionMm"]
            right_position = connector["rightPositionMm"]
            slip = connector["slipMm"]
            rows.append(
                {
                    "caseId": case_id,
                    "connector": connector["connector"],
                    "leftSite": connector["leftSite"],
                    "rightSite": connector["rightSite"],
                    "lockMode": lock_mode,
                    "alphaCommand": case_controls["alphaCommand"],
                    "zCommand": case_controls["zCommand"],
                    "backlash": case_config["backlash"],
                    "pinRadius": case_config["pinRadius"],
                    "holeRadius": case_config["holeRadius"],
                    "pinHoleClearanceMm": contact["dimensions"]["pinHoleClearanceMm"],
                    "predictedLeftXmm": left_position["x"],
                    "predictedLeftYmm": left_position["y"],
                    "predictedLeftZmm": left_position["z"],
                    "predictedRightXmm": right_position["x"],
                    "predictedRightYmm": right_position["y"],
                    "predictedRightZmm": right_position["z"],
                    "predictedSlipXmm": slip["x"],
                    "predictedSlipYmm": slip["y"],
                    "predictedSlipZmm": slip["z"],
                    "predictedLateralSlipMm": connector["lateralSlipMm"],
                    "predictedVerticalSlipMm": connector["verticalSlipMm"],
                    "predictedTotalSlipMm": connector["totalSlipMm"],
                    "predictedLateralExcessMm": connector["lateralExcessMm"],
                    "predictedVerticalExcessMm": connector["verticalExcessMm"],
                    "predictedContactMode": connector["contactMode"],
                    "observedLeftXmm": "",
                    "observedLeftYmm": "",
                    "observedLeftZmm": "",
                    "observedRightXmm": "",
                    "observedRightYmm": "",
                    "observedRightZmm": "",
                    "observedLateralSlipMm": "",
                    "observedVerticalSlipMm": "",
                    "observedTotalSlipMm": "",
                    "observedContactMode": "",
                    "lockHeldObserved": "",
                    "measurementSource": "",
                    "notes": "",
                }
            )
    return rows


def two_cell_connector_measurement_template(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> dict[str, Any]:
    rows = two_cell_connector_measurement_template_rows(config, controls)
    return {
        "schema": TWO_CELL_CONNECTOR_MEASUREMENT_TEMPLATE_SCHEMA,
        "model": "two-cell-suite-connector-slip-fillable-template",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "rows": rows,
        "summary": {
            "caseCount": len({row["caseId"] for row in rows}),
            "connectorRowCount": len(rows),
            "connectorNames": sorted({row["connector"] for row in rows}),
            "measurementRequired": True,
            "physicalAccuracyValidated": False,
        },
        "measurementColumns": [
            "observedLeftXmm",
            "observedLeftYmm",
            "observedLeftZmm",
            "observedRightXmm",
            "observedRightYmm",
            "observedRightZmm",
            "observedLateralSlipMm",
            "observedVerticalSlipMm",
            "observedTotalSlipMm",
            "observedContactMode",
            "lockHeldObserved",
            "measurementSource",
            "notes",
        ],
        "instructions": [
            "Measure connector marker positions in millimeters for each upper, middle, and lower pair.",
            "Keep the caseId and connector columns unchanged so comparisons can align rows.",
            "Use lockHeldObserved to record whether state or position lock behavior matched the fixture.",
        ],
    }


def _fidelity_matrix_case_id(row: dict[str, Any]) -> str:
    return f"FM_{int(row['index']):03d}"


def two_cell_fidelity_matrix_measurement_template_rows(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **matrix_overrides: Any,
) -> list[dict[str, Any]]:
    matrix = two_cell_physical_fidelity_matrix(config, controls, **matrix_overrides)
    matrix_rows = {int(row["index"]): row for row in matrix["rows"]}
    rows: list[dict[str, Any]] = []
    for connector in matrix["connectorRows"]:
        matrix_index = int(connector["matrixIndex"])
        parent = matrix_rows[matrix_index]
        left_position = connector["leftPositionMm"]
        right_position = connector["rightPositionMm"]
        slip = connector["slipMm"]
        rows.append(
            {
                "caseId": _fidelity_matrix_case_id(parent),
                "matrixIndex": matrix_index,
                "connector": connector["connector"],
                "leftSite": connector["leftSite"],
                "rightSite": connector["rightSite"],
                "actuationCase": parent["actuationCase"],
                "lockMode": parent["lockMode"],
                "backlash": parent["backlash"],
                "alphaCommand": parent["alphaCommand"],
                "zCommand": parent["zCommand"],
                "gravityForce": parent["gravityForce"],
                "pinRadius": parent["pinRadius"],
                "holeRadius": parent["holeRadius"],
                "pinHoleClearanceMm": parent["pinHoleClearanceMm"],
                "predictedLeftCellX": parent["leftX"],
                "predictedLeftCellY": 0.0,
                "predictedLeftCellZ": parent["leftZ"],
                "predictedRightCellX": parent["rightX"],
                "predictedRightCellY": 0.0,
                "predictedRightCellZ": parent["rightZ"],
                "predictedLeftAlpha": parent["leftAlpha"],
                "predictedRightAlpha": parent["rightAlpha"],
                "predictedLeftTheta": parent["leftTheta"],
                "predictedRightTheta": parent["rightTheta"],
                "predictedLeftXmm": left_position["x"],
                "predictedLeftYmm": left_position["y"],
                "predictedLeftZmm": left_position["z"],
                "predictedRightXmm": right_position["x"],
                "predictedRightYmm": right_position["y"],
                "predictedRightZmm": right_position["z"],
                "predictedSlipXmm": slip["x"],
                "predictedSlipYmm": slip["y"],
                "predictedSlipZmm": slip["z"],
                "predictedLateralSlipMm": connector["lateralSlipMm"],
                "predictedVerticalSlipMm": connector["verticalSlipMm"],
                "predictedTotalSlipMm": connector["totalSlipMm"],
                "predictedLateralExcessMm": connector["lateralExcessMm"],
                "predictedVerticalExcessMm": connector["verticalExcessMm"],
                "predictedContactMode": connector["contactMode"],
                "predictedLockPass": parent["lockPass"],
                "observedLeftCellX": "",
                "observedLeftCellY": "",
                "observedLeftCellZ": "",
                "observedRightCellX": "",
                "observedRightCellY": "",
                "observedRightCellZ": "",
                "observedLeftAlpha": "",
                "observedRightAlpha": "",
                "observedLeftTheta": "",
                "observedRightTheta": "",
                "observedLeftXmm": "",
                "observedLeftYmm": "",
                "observedLeftZmm": "",
                "observedRightXmm": "",
                "observedRightYmm": "",
                "observedRightZmm": "",
                "observedLateralSlipMm": "",
                "observedVerticalSlipMm": "",
                "observedTotalSlipMm": "",
                "observedContactMode": "",
                "lockHeldObserved": "",
                "measurementSource": "",
                "notes": "",
            }
        )
    return rows


def two_cell_fidelity_matrix_measurement_template(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **matrix_overrides: Any,
) -> dict[str, Any]:
    matrix = two_cell_physical_fidelity_matrix(config, controls, **matrix_overrides)
    rows = two_cell_fidelity_matrix_measurement_template_rows(config, controls, **matrix_overrides)
    return {
        "schema": TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA,
        "model": "two-cell-fidelity-matrix-connector-and-cell-fillable-template",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "matrixSummary": matrix["summary"],
        "rows": rows,
        "summary": {
            "matrixRowCount": matrix["summary"]["rowCount"],
            "connectorRowCount": len(rows),
            "connectorNames": sorted({row["connector"] for row in rows}),
            "measurementRequired": True,
            "physicalAccuracyValidated": False,
        },
        "measurementColumns": [
            "observedLeftCellX",
            "observedLeftCellY",
            "observedLeftCellZ",
            "observedRightCellX",
            "observedRightCellY",
            "observedRightCellZ",
            "observedLeftAlpha",
            "observedRightAlpha",
            "observedLeftTheta",
            "observedRightTheta",
            "observedLeftXmm",
            "observedLeftYmm",
            "observedLeftZmm",
            "observedRightXmm",
            "observedRightYmm",
            "observedRightZmm",
            "observedLateralSlipMm",
            "observedVerticalSlipMm",
            "observedTotalSlipMm",
            "observedContactMode",
            "lockHeldObserved",
            "measurementSource",
            "notes",
        ],
        "instructions": [
            "Keep caseId, matrixIndex, and connector unchanged so rows align with the prediction matrix.",
            "Record normalized cell-center motion if tracked in simulator units, and connector marker positions in millimeters.",
            "Use lockHeldObserved for right_state_locked and right_position_locked cases.",
        ],
    }


def _manifest_lock_flags(lock_mode: str) -> dict[str, bool]:
    return {
        "leftPositionLocked": lock_mode != "left_free",
        "rightStateLocked": lock_mode == "right_state_locked",
        "rightPositionLocked": lock_mode == "right_position_locked",
    }


def two_cell_external_fidelity_matrix_manifest_rows(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **matrix_overrides: Any,
) -> list[dict[str, Any]]:
    matrix = two_cell_physical_fidelity_matrix(config, controls, **matrix_overrides)
    rows = []
    for row in matrix["rows"]:
        case_id = _fidelity_matrix_case_id(row)
        rows.append(
            {
                "caseId": case_id,
                "matrixIndex": row["index"],
                "actuationCase": row["actuationCase"],
                "lockMode": row["lockMode"],
                "backlash": row["backlash"],
                "alphaCommand": row["alphaCommand"],
                "zCommand": row["zCommand"],
                "gravityForce": row["gravityForce"],
                "pinRadius": row["pinRadius"],
                "holeRadius": row["holeRadius"],
                "pinHoleClearanceMm": row["pinHoleClearanceMm"],
                **_manifest_lock_flags(str(row["lockMode"])),
                "expectedLeftX": row["leftX"],
                "expectedLeftY": 0.0,
                "expectedLeftZ": row["leftZ"],
                "expectedRightX": row["rightX"],
                "expectedRightY": 0.0,
                "expectedRightZ": row["rightZ"],
                "expectedLeftAlpha": row["leftAlpha"],
                "expectedRightAlpha": row["rightAlpha"],
                "expectedLeftTheta": row["leftTheta"],
                "expectedRightTheta": row["rightTheta"],
                "expectedContactMode": row["contactMode"],
                "expectedConnectorMaxLateralSlipMm": row["connectorMaxLateralSlipMm"],
                "expectedConnectorMaxVerticalSlipMm": row["connectorMaxVerticalSlipMm"],
                "expectedConnectorMaxVerticalExcessMm": row["connectorMaxVerticalExcessMm"],
                "expectedConnectorTotalContactPenalty": row["connectorTotalContactPenalty"],
                "mjcfProxyPath": f"external_fidelity_matrix_mjcf/{case_id}.xml",
                "engineResultPath": f"external_fidelity_matrix_results/{case_id}.json",
                "measurementRows": 3,
                "runStatus": "pending-external-engine",
            }
        )
    return rows


def two_cell_external_fidelity_matrix_manifest(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    engine: str = "mujoco",
    intake_dir: str | Path | None = None,
    **matrix_overrides: Any,
) -> dict[str, Any]:
    config = replace(config or LatticeConfig(rows=1, cols=2), rows=1, cols=2)
    controls = _controls(controls)
    matrix = two_cell_physical_fidelity_matrix(config, controls, **matrix_overrides)
    measurement_template = two_cell_fidelity_matrix_measurement_template(config, controls, **matrix_overrides)
    readiness_kwargs = {"intake_dir": intake_dir} if intake_dir is not None else {}
    readiness = two_cell_segmented_cad_readiness_report(config, controls, **readiness_kwargs)
    intake_validation = two_cell_segmented_cad_intake_validation_report(config, controls, **readiness_kwargs)
    case_rows = two_cell_external_fidelity_matrix_manifest_rows(config, controls, **matrix_overrides)
    connector_rows = measurement_template["rows"]
    engine_key = str(engine or "mujoco").strip().lower()
    proxy_engine_ready = engine_key in {
        "mujoco",
        "gazebo",
        "isaac",
        "isaac-sim",
        "pybullet",
        "chrono",
        "project-chrono",
    }
    exact_cad_ready = bool(intake_validation["summary"].get("intakeValidationReady"))
    missing_evidence = sorted(
        set(
            [
                *readiness["summary"].get("missingEvidence", []),
                *intake_validation["summary"].get("missingEvidence", []),
                "externalEngineRunResults",
                "filledFidelityMatrixMeasurements",
                "benchCoordinateHoldout",
            ]
        )
    )
    return {
        "schema": TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA,
        "model": "dense-two-cell-external-engine-fidelity-matrix-handoff",
        "engine": {
            "requested": engine_key,
            "proxyEngineReady": proxy_engine_ready,
            "exactSegmentedCadReady": exact_cad_ready,
            "runScope": "252 two-cell cases with three connector measurement rows per case",
        },
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "segmentedCadReadiness": readiness["summary"],
        "segmentedCadIntakeValidation": intake_validation["summary"],
        "axes": matrix["axes"],
        "caseRows": case_rows,
        "connectorMeasurementRows": connector_rows,
        "resultSchema": {
            "caseResultJson": {
                "caseId": "FM_000",
                "sourceEngine": engine_key,
                "finalBodies": [
                    {
                        "bodyId": "left",
                        "x": "float model units",
                        "y": "float model units",
                        "z": "float model units",
                        "alpha": "float",
                        "theta": "float degrees",
                    },
                    {
                        "bodyId": "right",
                        "x": "float model units",
                        "y": "float model units",
                        "z": "float model units",
                        "alpha": "float",
                        "theta": "float degrees",
                    },
                ],
                "connectors": [
                    {
                        "connector": "upper|middle|lower",
                        "observedLeftXmm": "float",
                        "observedLeftYmm": "float",
                        "observedLeftZmm": "float",
                        "observedRightXmm": "float",
                        "observedRightYmm": "float",
                        "observedRightZmm": "float",
                        "observedLateralSlipMm": "float",
                        "observedVerticalSlipMm": "float",
                        "observedTotalSlipMm": "float",
                        "observedContactMode": "inside_clearance|clearance_edge|contacting",
                        "lockHeldObserved": "held|released|failed",
                    }
                ],
            },
            "fillableCsv": "two_cell_fidelity_matrix_measurement_template.csv",
        },
        "runInstructions": [
            "Generate or export one MJCF/SDF/USD case per caseId using the listed controls, locks, backlash, and hole radius.",
            "Run each case to static equilibrium under the listed gravityForce and actuator commands.",
            "Record final left/right body centers, alpha, theta, connector marker positions, slip, contact mode, and lock-held status.",
            "Fill the connectorMeasurementRows-compatible CSV, then compare with compare_two_cell_fidelity_matrix_measurements.",
            "Do not claim exact physical fidelity unless segmented CAD intake, external engine results, and bench holdout measurements pass.",
        ],
        "summary": {
            "status": "external-fidelity-manifest-ready-needs-engine-results"
            if proxy_engine_ready
            else "external-fidelity-manifest-engine-review-needed",
            "caseCount": len(case_rows),
            "connectorMeasurementRowCount": len(connector_rows),
            "matrixRowCount": matrix["summary"]["rowCount"],
            "proxyEngineReady": proxy_engine_ready,
            "exactSegmentedCadReady": exact_cad_ready,
            "physicalAccuracyValidated": False,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimLabels": {
            "handoff": "external engine run manifest aligned to the dense two-cell fidelity matrix",
            "geometry": "proxy-ready now; exact geometry requires filled segmented CAD intake",
            "physicalAccuracy": "false until external solver and bench holdout measurements agree",
        },
    }


def _handoff_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value] if value else []
    if isinstance(value, Iterable):
        return [str(item) for item in value if str(item)]
    return [str(value)]


def _handoff_case_key(case: dict[str, Any], source: str, rank: int) -> str:
    return str(
        case.get("caseId")
        or case.get("phaseDiagramCaseId")
        or case.get("matrixCaseId")
        or case.get("index")
        or f"{source}-{rank}"
    )


def _handoff_lock_flags(case: dict[str, Any]) -> dict[str, bool]:
    lock_mode = str(case.get("lockMode", "free"))
    return {
        "leftPositionLocked": bool(case.get("leftPositionLocked", lock_mode != "left_free")),
        "rightStateLocked": bool(case.get("rightStateLocked", lock_mode == "right_state_locked")),
        "rightPositionLocked": bool(case.get("rightPositionLocked", lock_mode == "right_position_locked")),
    }


def _exact_contact_handoff_case_row(
    case: dict[str, Any],
    *,
    source: str,
    category: str,
    rank: int,
    can_run_exact_contact: bool,
) -> dict[str, Any]:
    lock_flags = _handoff_lock_flags(case)
    return {
        "rank": rank,
        "caseId": _handoff_case_key(case, source, rank),
        "source": source,
        "category": category,
        "phase": str(case.get("phase", "")),
        "actuationCase": str(case.get("actuationCase", "")),
        "lockMode": str(case.get("lockMode", "free")),
        "backlash": case.get("backlash", ""),
        "holeRadius": case.get("holeRadius", ""),
        "pinHoleClearanceMm": case.get("pinHoleClearanceMm", ""),
        "effectiveBacklash": case.get("effectiveBacklash", case.get("backlash", "")),
        "effectiveHoleRadius": case.get("effectiveHoleRadius", case.get("holeRadius", "")),
        "effectivePinHoleClearanceMm": case.get(
            "effectivePinHoleClearanceMm",
            case.get("pinHoleClearanceMm", ""),
        ),
        "alphaCommand": case.get("alphaCommand", ""),
        "zCommand": case.get("zCommand", ""),
        "gravityForce": case.get("gravityForce", ""),
        "leftZ": case.get("leftZ", ""),
        "rightZ": case.get("rightZ", ""),
        "rightAlpha": case.get("rightAlpha", ""),
        "connectorActiveContactCount": case.get("connectorActiveContactCount", ""),
        "connectorMaxVerticalSlipMm": case.get("connectorMaxVerticalSlipMm", ""),
        "connectorMaxVerticalExcessMm": case.get("connectorMaxVerticalExcessMm", ""),
        "connectorTotalContactPenalty": case.get("connectorTotalContactPenalty", ""),
        **lock_flags,
        "expectedConnectorCount": 3,
        "requiredCadAssets": list(EXACT_CONTACT_HANDOFF_REQUIRED_CAD_ASSETS),
        "requiredMeasuredInputs": list(EXACT_CONTACT_HANDOFF_REQUIRED_MEASURED_INPUTS),
        "engineHandoffReady": can_run_exact_contact,
        "canRunExactContact": can_run_exact_contact,
        "physicalAccuracyValidated": False,
        "notes": (
            "Run with segmented rigid bodies, explicit pin-hole contacts, gravity, locks, and actuator commands; "
            "compare final cell centers and upper/middle/lower connector markers against bench data."
        ),
    }


def two_cell_exact_contact_handoff_plan(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    priority_limit: int = 24,
    intake_dir: str | Path | None = None,
    engine: str = "mujoco",
) -> dict[str, Any]:
    """Rank the first exact-contact two-cell cases once segmented CAD is available."""

    config = replace(config or LatticeConfig(rows=1, cols=2), rows=1, cols=2)
    controls = _controls(controls)
    readiness_kwargs = {"intake_dir": intake_dir} if intake_dir is not None else {}
    readiness = two_cell_segmented_cad_readiness_report(config, controls, **readiness_kwargs)
    intake_validation = two_cell_segmented_cad_intake_validation_report(config, controls, **readiness_kwargs)
    atlas = two_cell_physical_response_atlas(config, controls)
    transition_report = two_cell_radius_backlash_transition_report(config, controls)
    manifest = two_cell_external_fidelity_matrix_manifest(
        config,
        controls,
        engine=engine,
        intake_dir=intake_dir,
    )
    contact_decomposition = two_cell_cad_contact_decomposition_spec(config, controls, **readiness_kwargs)
    engine_targets = [
        {
            "engine": str(row.get("engine", "")),
            "status": str(row.get("status", "")),
            "expectedAssets": [asset["path"] for asset in row.get("expectedAssets", [])],
            "detectedAssets": [asset["path"] for asset in row.get("detectedAssets", [])],
        }
        for row in intake_validation.get("engineHandoffValidationRows", [])
    ]
    engine_handoff_ready = any(target["status"] == "ready-for-external-run" for target in engine_targets)
    can_run_exact_contact = bool(intake_validation["summary"].get("canAttemptExactRigidBodyContact")) and engine_handoff_ready
    missing_evidence = sorted(
        set(
            [
                *readiness["summary"].get("missingEvidence", []),
                *intake_validation["summary"].get("missingEvidence", []),
                *manifest["summary"].get("missingEvidence", []),
                *contact_decomposition["summary"].get("missingEvidence", []),
                "segmentedCADBodies",
                "externalRigidBodyContactResults",
                "benchCoordinateHoldout",
            ]
        )
    )

    rows_by_id: dict[str, dict[str, Any]] = {}
    ordered_rows: list[dict[str, Any]] = []

    def add_cases(source: str, category: str, cases: list[dict[str, Any]]) -> None:
        for case in cases:
            key = _handoff_case_key(case, source, len(ordered_rows) + 1)
            row = _exact_contact_handoff_case_row(
                case,
                source=source,
                category=str(case.get("category") or category),
                rank=len(ordered_rows) + 1,
                can_run_exact_contact=can_run_exact_contact,
            )
            existing = rows_by_id.get(key)
            if existing is None:
                rows_by_id[key] = row
                ordered_rows.append(row)
            else:
                existing_sources = set(str(existing["source"]).split(";"))
                existing_sources.add(source)
                existing["source"] = ";".join(sorted(existing_sources))
                existing_categories = set(str(existing["category"]).split(";"))
                existing_categories.add(str(case.get("category") or category))
                existing["category"] = ";".join(sorted(existing_categories))

    add_cases("transitionReport", "transition-boundary", list(transition_report["measurementCases"]))
    bench_priority = atlas.get("benchPriority", {})
    add_cases("responseAtlas.first", "high-response", list(bench_priority.get("firstCasesToMeasure", [])))
    add_cases("responseAtlas.lock", "lock-diagnostic", list(bench_priority.get("lockCasesToMeasure", [])))
    add_cases("responseAtlas.clearance", "clearance-dieoff", list(bench_priority.get("clearanceCasesToMeasure", [])))

    limit = max(1, int(priority_limit))
    case_rows = ordered_rows[:limit]
    for rank, row in enumerate(case_rows, start=1):
        row["rank"] = rank

    return {
        "schema": TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA,
        "model": "two-cell-exact-rigid-body-contact-handoff-plan",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "segmentedCadReadiness": readiness["summary"],
        "segmentedCadIntakeValidation": intake_validation["summary"],
        "externalFidelityManifest": manifest["summary"],
        "cadContactDecomposition": contact_decomposition["summary"],
        "engineTargets": engine_targets,
        "caseRows": case_rows,
        "protocol": {
            "startFrom": "one-cell A360 CAD reference plus two connected-cell X-neighbor assembly",
            "requiredRigidBodies": [
                "left upper moving cell body",
                "left lower fixture body",
                "right upper moving cell body",
                "right lower fixture body",
                "connector pins or screw bodies",
            ],
            "requiredContacts": [
                "pin-hole radial clearance contacts",
                "vertical slide/free-play limits",
                "lock crown contact/stop angles",
                "ground/support contacts under gravity",
            ],
            "collisionPrimitiveContract": contact_decomposition["collisionPrimitiveContract"],
            "observables": [
                "left/right cell center xyz",
                "left/right alpha and theta",
                "upper/middle/lower connector marker xyz",
                "contact mode",
                "lock held/released/failed",
            ],
        },
        "summary": {
            "status": "ready-for-exact-contact-run" if can_run_exact_contact else "needs-segmented-cad-contact-inputs",
            "caseCount": len(case_rows),
            "transitionCaseCount": sum(1 for row in case_rows if "transitionReport" in str(row["source"])),
            "atlasCaseCount": sum(1 for row in case_rows if "responseAtlas" in str(row["source"])),
            "engineHandoffReady": engine_handoff_ready,
            "intakeValidationReady": bool(intake_validation["summary"].get("intakeValidationReady")),
            "canAttemptExactRigidBodyContact": bool(
                intake_validation["summary"].get("canAttemptExactRigidBodyContact")
            ),
            "canRunExactContact": can_run_exact_contact,
            "physicalAccuracyValidated": False,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimBoundary": {
            "allowedClaim": "ranked exact-contact and bench handoff plan derived from current reduced two-cell diagnostics",
            "blockedClaim": "exact real RAD mechanics until segmented CAD, external contact results, and bench holdouts pass",
            "physicalAccuracyValidated": False,
        },
    }


def two_cell_external_results_template_rows(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> list[dict[str, Any]]:
    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    rows = []
    for case_id, case_controls in _named_case_controls(controls):
        solved = solve_two_cell_quasistatic(config, case_controls)
        for cell in solved["cells"]:
            center = cell["center"]
            rows.append(
                {
                    "caseId": case_id,
                    "bodyId": cell["id"],
                    "bodyName": f"{cell['id']}_cell",
                    "lockMode": _lock_mode_from_controls(case_controls),
                    "alphaCommand": case_controls.alpha_command,
                    "zCommand": case_controls.z_command,
                    "pinRadius": solved["dimensions"]["pinRadius"],
                    "holeRadius": solved["dimensions"]["holeRadius"],
                    "expectedX": center["x"],
                    "expectedY": center["y"],
                    "expectedZ": center["z"],
                    "expectedAlpha": cell["alpha"],
                    "expectedTheta": cell["theta"],
                    "externalX": "",
                    "externalY": "",
                    "externalZ": "",
                    "externalAlpha": "",
                    "externalTheta": "",
                    "sourceEngine": "",
                    "notes": "",
                }
            )
    return rows


def two_cell_external_results_template(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> dict[str, Any]:
    rows = two_cell_external_results_template_rows(config, controls)
    return {
        "schema": TWO_CELL_EXTERNAL_TEMPLATE_SCHEMA,
        "model": "two-cell-external-results-fillable-template",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "rows": rows,
        "measurementColumns": [
            "externalX",
            "externalY",
            "externalZ",
            "externalAlpha",
            "externalTheta",
            "sourceEngine",
            "notes",
        ],
        "summary": {
            "rowCount": len(rows),
            "bodyRowsPerCase": 2,
            "caseCount": len({row["caseId"] for row in rows}),
            "requiresExternalEngineOrBenchData": True,
        },
    }


def _optional_number(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _row_value(row: dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
    return None


def two_cell_external_results_from_rows(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    normalized = []
    for row in rows:
        case_id = str(_row_value(row, "caseId", "case_id") or "").strip()
        body_id = str(_row_value(row, "bodyId", "body_id", "body", "bodyName", "body_name") or "").strip()
        if not case_id or not body_id:
            continue
        body_key = "right" if "right" in body_id.lower() else "left"
        normalized.append(
            {
                "caseId": case_id,
                "bodyId": body_key,
                "bodyName": str(_row_value(row, "bodyName", "body_name") or f"{body_key}_cell"),
                "lockMode": str(_row_value(row, "lockMode", "lock_mode") or "").strip(),
                "alphaCommand": _optional_number(_row_value(row, "alphaCommand", "alpha_command")),
                "zCommand": _optional_number(_row_value(row, "zCommand", "z_command")),
                "pinRadius": _optional_number(_row_value(row, "pinRadius", "pin_radius")),
                "holeRadius": _optional_number(_row_value(row, "holeRadius", "hole_radius")),
                "externalX": _optional_number(_row_value(row, "externalX", "x", "finalX", "final_x")),
                "externalY": _optional_number(_row_value(row, "externalY", "y", "finalY", "final_y")),
                "externalZ": _optional_number(_row_value(row, "externalZ", "z", "finalZ", "final_z")),
                "externalAlpha": _optional_number(_row_value(row, "externalAlpha", "alpha", "finalAlpha", "final_alpha")),
                "externalTheta": _optional_number(_row_value(row, "externalTheta", "theta", "finalTheta", "final_theta")),
                "sourceEngine": str(_row_value(row, "sourceEngine", "source_engine", "engine") or ""),
                "notes": str(_row_value(row, "notes") or ""),
            }
        )
    return normalized


def two_cell_external_results_from_csv(text: str) -> list[dict[str, Any]]:
    return two_cell_external_results_from_rows(csv.DictReader(io.StringIO(text)))


def two_cell_external_results_from_json(text: str) -> list[dict[str, Any]]:
    payload = json.loads(text)
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict) and isinstance(payload.get("rows"), list):
        rows = payload["rows"]
    elif isinstance(payload, dict) and isinstance(payload.get("results"), list):
        rows = payload["results"]
    else:
        rows = []
    return two_cell_external_results_from_rows(rows)


def _controls_for_external_row(
    row: dict[str, Any],
    controls_by_case: dict[str, TwoCellBenchControls],
    base_controls: TwoCellBenchControls,
) -> TwoCellBenchControls:
    controls = controls_by_case.get(str(row["caseId"]), base_controls)
    lock_mode = str(row.get("lockMode") or _lock_mode_from_controls(controls))
    controls = _lock_mode_controls(controls, lock_mode)
    return _controls(
        controls,
        alpha_command=row["alphaCommand"] if row.get("alphaCommand") is not None else controls.alpha_command,
        z_command=row["zCommand"] if row.get("zCommand") is not None else controls.z_command,
    )


def compare_two_cell_external_results(
    external_results: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-6,
) -> dict[str, Any]:
    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    controls_by_case = {case_id: case_controls for case_id, case_controls in _named_case_controls(controls)}
    rows = []
    scalar_squared = 0.0
    scalar_count = 0
    position_squared = 0.0
    position_count = 0
    max_abs_scalar = 0.0
    max_position = 0.0
    missing_results = 0
    lock_violations = 0
    normalized = two_cell_external_results_from_rows(external_results)
    for row in normalized:
        case_controls = _controls_for_external_row(row, controls_by_case, controls)
        solved = solve_two_cell_quasistatic(
            config,
            case_controls,
            pin_radius=row["pinRadius"] if row.get("pinRadius") is not None else None,
            hole_radius=row["holeRadius"] if row.get("holeRadius") is not None else None,
        )
        cell = solved["cells"][1 if row["bodyId"] == "right" else 0]
        expected = {
            "externalX": cell["center"]["x"],
            "externalY": cell["center"]["y"],
            "externalZ": cell["center"]["z"],
            "externalAlpha": cell["alpha"],
            "externalTheta": cell["theta"],
        }
        residual = {}
        observed_fields = {}
        for key, expected_value in expected.items():
            observed = row.get(key)
            if observed is None:
                continue
            observed_fields[key] = observed
            delta = float(observed) - float(expected_value)
            residual[key] = delta
            scalar_squared += delta * delta
            scalar_count += 1
            max_abs_scalar = max(max_abs_scalar, abs(delta))
        if not observed_fields:
            missing_results += 1
        if all(key in observed_fields for key in ("externalX", "externalY", "externalZ")):
            position_error = math.sqrt(
                residual["externalX"] ** 2 + residual["externalY"] ** 2 + residual["externalZ"] ** 2
            )
            position_squared += position_error * position_error
            position_count += 1
            max_position = max(max_position, position_error)
        else:
            position_error = None
        expected_fixed = bool(cell["positionLocked"])
        expected_state_locked = bool(cell["stateLocked"])
        lock_violation = False
        if expected_fixed and position_error is not None and position_error > tolerance:
            lock_violation = True
        if expected_state_locked and abs(residual.get("externalAlpha", 0.0)) > tolerance:
            lock_violation = True
        if lock_violation:
            lock_violations += 1
        rows.append(
            {
                "caseId": row["caseId"],
                "bodyId": row["bodyId"],
                "bodyName": row["bodyName"],
                "lockMode": _lock_mode_from_controls(case_controls),
                "sourceEngine": row["sourceEngine"],
                "expected": expected,
                "observed": observed_fields,
                "residual": residual,
                "positionError": position_error,
                "withinTolerance": (
                    bool(observed_fields)
                    and max((abs(value) for value in residual.values()), default=0.0) <= tolerance
                ),
                "lockViolation": lock_violation,
            }
        )
    rms_scalar = math.sqrt(scalar_squared / scalar_count) if scalar_count else None
    rms_position = math.sqrt(position_squared / position_count) if position_count else None
    missing_evidence = []
    if missing_results:
        missing_evidence.append("externalResultRows")
    if scalar_count == 0:
        missing_evidence.append("observedExternalScalars")
    if lock_violations:
        missing_evidence.append("lockConstraintAgreement")
    if max_abs_scalar > tolerance:
        missing_evidence.append("externalTolerance")
    ready = len(missing_evidence) == 0
    return {
        "schema": TWO_CELL_EXTERNAL_COMPARISON_SCHEMA,
        "model": "two-cell-external-vs-quasistatic-comparison",
        "tolerance": tolerance,
        "rows": rows,
        "summary": {
            "status": "two-cell-external-comparison-ready" if ready else "needs-two-cell-external-comparison-evidence",
            "externalComparisonReady": ready,
            "rowCount": len(rows),
            "observedScalarCount": scalar_count,
            "positionResultCount": position_count,
            "missingExternalResultCount": missing_results,
            "lockViolationCount": lock_violations,
            "rmsScalarError": rms_scalar,
            "rmsPositionError": rms_position,
            "maxAbsScalarError": max_abs_scalar if scalar_count else None,
            "maxPositionError": max_position if position_count else None,
            "withinTolerance": ready,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimLabels": {
            "comparison": "external engine or bench coordinates compared to internal quasistatic solve",
            "physicalAccuracy": "agreement diagnostic only; physical validity also needs real measured data",
        },
    }


def _validation_record(
    test_id: str,
    *,
    passed: bool,
    status: str | None = None,
    metrics: dict[str, Any] | None = None,
    evidence: list[str] | None = None,
    limitations: list[str] | None = None,
    missing_evidence: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "testId": test_id,
        "status": status or ("pass" if passed else "fail"),
        "pass": bool(passed),
        "metrics": metrics or {},
        "evidence": evidence or [],
        "limitations": limitations or [],
        "missingEvidence": missing_evidence or [],
    }


def _nondecreasing(values: list[float], tolerance: float) -> bool:
    return all(b + tolerance >= a for a, b in zip(values, values[1:]))


def _nonincreasing(values: list[float], tolerance: float) -> bool:
    return all(b <= a + tolerance for a, b in zip(values, values[1:]))


def two_cell_physics_validation_report(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    tolerance: float = 1e-7,
) -> dict[str, Any]:
    """Consolidated internal/external gate report for the two-cell CAD proxy."""

    config = config or LatticeConfig(rows=1, cols=2)
    controls = _controls(controls)
    probe_controls = _controls(
        controls,
        alpha_command=controls.alpha_command if abs(controls.alpha_command) > 0.1 else -0.35,
        z_command=controls.z_command if abs(controls.z_command) > 0.1 else 0.35,
    )
    tests: list[dict[str, Any]] = []

    cad_layout = cad_rad_cell_layout(config)
    envelope = cad_layout["envelopeCheck"]
    cad_layout_pass = (
        cad_layout["schema"] == CAD_RAD_CELL_LAYOUT_SCHEMA
        and len(cad_layout["padSites"]) == 8
        and bool(envelope["matchesBoundingBox"])
    )
    tests.append(
        _validation_record(
            "cad_layout_envelope",
            passed=cad_layout_pass,
            metrics={
                "padSiteCount": len(cad_layout["padSites"]),
                "maxOuterEnvelopeMm": envelope["maxOuterEnvelopeMm"],
                "halfWidthMm": envelope["halfWidthMm"],
                "matchesBoundingBox": envelope["matchesBoundingBox"],
            },
            evidence=["A360 bounding box and nominal 3.4 mm hole label"],
            limitations=["dimension-consistent layout, not direct BREP vertex extraction"],
        )
    )

    left_free_sweep = sweep_two_cell_backlash(config, _lock_mode_controls(probe_controls, "left_free"))
    sweep_rows = left_free_sweep["rows"]
    clearances = [float(row["clearance"]) for row in sweep_rows]
    neighbor_z = [abs(float(row["leftResidualZ"])) for row in sweep_rows]
    clearance_pass = (
        bool(sweep_rows)
        and _nondecreasing(clearances, tolerance)
        and _nonincreasing(neighbor_z, tolerance)
        and bool(left_free_sweep["trend"]["clearanceIncreasesBacklash"])
        and bool(left_free_sweep["trend"]["neighborResponseDropsWithClearance"])
    )
    tests.append(
        _validation_record(
            "hole_clearance_monotonicity",
            passed=clearance_pass,
            metrics={
                "rowCount": len(sweep_rows),
                "clearanceStart": clearances[0] if clearances else None,
                "clearanceEnd": clearances[-1] if clearances else None,
                "neighborResidualZStart": neighbor_z[0] if neighbor_z else None,
                "neighborResidualZEnd": neighbor_z[-1] if neighbor_z else None,
            },
            evidence=[
                "sweep_two_cell_backlash with left cell free",
                "larger pin-hole clearance increases dead-zone and reduces neighbor residual z",
            ],
        )
    )

    quasistatic_sweep = sweep_two_cell_quasistatic(
        config,
        probe_controls,
        alpha_commands=(-abs(probe_controls.alpha_command), 0.0, abs(probe_controls.alpha_command)),
        z_commands=(-abs(probe_controls.z_command), 0.0, abs(probe_controls.z_command)),
        hole_radii=(config.pin_radius, config.hole_radius, max(config.hole_radius, controls.hole_sweep_max)),
        lock_modes=("free", "right_state_locked", "right_position_locked", "left_free"),
    )
    quasistatic_summary = quasistatic_sweep["summary"]
    quasistatic_pass = quasistatic_summary["solverFailureCount"] == 0
    tests.append(
        _validation_record(
            "quasistatic_sweep_success",
            passed=quasistatic_pass,
            metrics={
                "rowCount": quasistatic_summary["rowCount"],
                "solverSuccessCount": quasistatic_summary["solverSuccessCount"],
                "solverFailureCount": quasistatic_summary["solverFailureCount"],
                "maxStoredEnergy": quasistatic_summary["maxStoredEnergy"],
            },
            evidence=["targeted alpha/z/hole/lock-mode quasistatic sweep"],
        )
    )

    connector_tight = two_cell_connector_contact_report(
        config,
        probe_controls,
        hole_radius=config.pin_radius,
    )
    connector_loose = two_cell_connector_contact_report(
        config,
        probe_controls,
        hole_radius=max(config.hole_radius, controls.hole_sweep_max),
    )
    connector_pass = (
        connector_tight["summary"]["connectorContactReady"]
        and connector_loose["summary"]["connectorContactReady"]
        and connector_loose["dimensions"]["pinHoleClearanceMm"] >= connector_tight["dimensions"]["pinHoleClearanceMm"]
        and connector_loose["summary"]["maxVerticalExcessMm"] <= connector_tight["summary"]["maxVerticalExcessMm"] + tolerance
    )
    tests.append(
        _validation_record(
            "connector_clearance_contact_dieoff",
            passed=connector_pass,
            metrics={
                "tightClearanceMm": connector_tight["dimensions"]["pinHoleClearanceMm"],
                "looseClearanceMm": connector_loose["dimensions"]["pinHoleClearanceMm"],
                "tightMaxVerticalExcessMm": connector_tight["summary"]["maxVerticalExcessMm"],
                "looseMaxVerticalExcessMm": connector_loose["summary"]["maxVerticalExcessMm"],
                "connectorCount": connector_tight["summary"]["connectorCount"],
            },
            evidence=["three CAD-derived two-cell connector pairs"],
            limitations=["connector contact is clearance-excess bookkeeping, not solved BREP collision"],
        )
    )

    state_locked = solve_two_cell_quasistatic(config, _lock_mode_controls(probe_controls, "right_state_locked"))
    state_right = state_locked["cells"][1]
    state_constraints = state_locked["constraints"]
    state_alpha_error = abs(float(state_constraints["rightStateAlphaError"] or 0.0))
    state_vertical_motion = abs(float(state_right["center"]["z"]))
    state_lock_pass = (
        bool(state_locked["solver"]["success"])
        and bool(state_constraints["stateLocksHold"])
        and state_alpha_error <= tolerance
        and state_vertical_motion > tolerance
        and bool(state_constraints["rightLockedVerticalMotionAllowed"])
    )
    tests.append(
        _validation_record(
            "state_lock_alpha_hold_vertical_free",
            passed=state_lock_pass,
            metrics={
                "rightAlpha": state_right["alpha"],
                "rightZ": state_right["center"]["z"],
                "rightStateAlphaError": state_alpha_error,
                "verticalMotionAllowed": state_constraints["rightLockedVerticalMotionAllowed"],
            },
            evidence=["state lock constrains alpha but does not force the right cell back to ground z"],
        )
    )

    position_locked = solve_two_cell_quasistatic(config, _lock_mode_controls(probe_controls, "right_position_locked"))
    position_right = position_locked["cells"][1]
    position_constraints = position_locked["constraints"]
    right_fixture_error = float(position_constraints["rightFixtureError"] or 0.0)
    position_pass = (
        bool(position_locked["solver"]["success"])
        and bool(position_constraints["positionLocksHold"])
        and right_fixture_error <= tolerance
        and abs(float(position_right["center"]["x"]) - config.cell_size) <= tolerance
        and abs(float(position_right["center"]["z"])) <= tolerance
    )
    tests.append(
        _validation_record(
            "position_lock_fixture_hold",
            passed=position_pass,
            metrics={
                "rightX": position_right["center"]["x"],
                "rightZ": position_right["center"]["z"],
                "rightFixtureError": right_fixture_error,
            },
            evidence=["position lock hard-fixes right-cell x/z fixture coordinates"],
        )
    )

    gravity_controls = _controls(probe_controls, alpha_command=0.0, z_command=0.0)
    low_clearance = solve_two_cell_quasistatic(
        config,
        gravity_controls,
        hole_radius=config.pin_radius,
    )
    high_hole = max(config.hole_radius, config.pin_radius + 0.1)
    high_clearance = solve_two_cell_quasistatic(
        config,
        gravity_controls,
        hole_radius=high_hole,
    )
    low_sag = max(0.0, -float(low_clearance["cells"][1]["center"]["z"]))
    high_sag = max(0.0, -float(high_clearance["cells"][1]["center"]["z"]))
    gravity_pass = (
        bool(low_clearance["solver"]["success"])
        and bool(high_clearance["solver"]["success"])
        and high_sag + tolerance >= low_sag
    )
    tests.append(
        _validation_record(
            "gravity_clearance_sag",
            passed=gravity_pass,
            metrics={
                "lowClearanceHoleRadius": config.pin_radius,
                "highClearanceHoleRadius": high_hole,
                "lowClearanceSag": low_sag,
                "highClearanceSag": high_sag,
            },
            evidence=["zero-command quasistatic gravity probe"],
            limitations=["gravity/contact parameters are uncalibrated reduced-model values"],
        )
    )

    z_positive = solve_two_cell_quasistatic(config, _controls(probe_controls, alpha_command=0.0, z_command=abs(probe_controls.z_command)))
    z_negative = solve_two_cell_quasistatic(config, _controls(probe_controls, alpha_command=0.0, z_command=-abs(probe_controls.z_command)))
    alpha_negative = solve_two_cell_quasistatic(
        config,
        _controls(probe_controls, alpha_command=-abs(probe_controls.alpha_command), z_command=0.0),
    )
    alpha_positive = solve_two_cell_quasistatic(
        config,
        _controls(probe_controls, alpha_command=abs(probe_controls.alpha_command), z_command=0.0),
    )
    polarity_pass = (
        bool(z_positive["solver"]["success"])
        and bool(z_negative["solver"]["success"])
        and bool(alpha_negative["solver"]["success"])
        and bool(alpha_positive["solver"]["success"])
        and float(z_positive["cells"][1]["center"]["z"]) > tolerance
        and float(z_negative["cells"][1]["center"]["z"]) < -tolerance
        and float(alpha_negative["cells"][1]["alpha"]) < config.initial_alpha
        and float(alpha_positive["cells"][1]["alpha"]) > config.initial_alpha
    )
    tests.append(
        _validation_record(
            "actuation_polarity",
            passed=polarity_pass,
            metrics={
                "positiveZ": z_positive["cells"][1]["center"]["z"],
                "negativeZ": z_negative["cells"][1]["center"]["z"],
                "negativeAlpha": alpha_negative["cells"][1]["alpha"],
                "positiveAlpha": alpha_positive["cells"][1]["alpha"],
            },
            evidence=["positive/negative z commands and contraction/expansion alpha commands"],
        )
    )

    mjcf = two_cell_mjcf_proxy_report(config, probe_controls)
    mjcf_summary = mjcf["summary"]
    mjcf_pass = (
        mjcf_summary["bodyCount"] == 2
        and mjcf_summary["connectorContactPairCount"] >= 3
        and mjcf_summary["actuatorCount"] >= 1
        and "<mujoco" in mjcf["xml"]
    )
    tests.append(
        _validation_record(
            "mjcf_proxy_ready",
            passed=mjcf_pass,
            metrics={
                "bodyCount": mjcf_summary["bodyCount"],
                "jointCount": mjcf_summary["jointCount"],
                "actuatorCount": mjcf_summary["actuatorCount"],
                "connectorContactPairCount": mjcf_summary["connectorContactPairCount"],
                "externalRunComplete": mjcf_summary["externalRunComplete"],
            },
            evidence=["CAD-dimensioned MuJoCo XML proxy exported from A360 one-cell metadata"],
            limitations=mjcf["limitations"],
            missing_evidence=list(mjcf_summary["missingEvidence"]),
        )
    )

    mujoco_run = two_cell_mujoco_proxy_run_report(config, probe_controls, steps=30)
    mujoco_run_summary = mujoco_run["summary"]
    tests.append(
        _validation_record(
            "mujoco_proxy_run_gate",
            passed=bool(mujoco_run_summary["twoCellMujocoProxyRunComplete"]),
            status="pass" if mujoco_run_summary["twoCellMujocoProxyRunComplete"] else "incomplete",
            metrics={
                "mujocoAvailable": mujoco_run_summary["mujocoAvailable"],
                "bodyResultCount": mujoco_run_summary["bodyResultCount"],
                "expectedBodyResultCount": mujoco_run_summary["expectedBodyResultCount"],
                "stepsCompleted": mujoco_run["solver"]["stepsCompleted"],
                "maxFixtureDisplacementMm": mujoco_run_summary["maxFixtureDisplacementMm"],
            },
            evidence=["optional MuJoCo execution path for the two-cell CAD proxy"],
            limitations=mujoco_run["limitations"],
            missing_evidence=list(mujoco_run_summary["missingEvidence"]),
        )
    )

    external_template = two_cell_external_results_template(config, probe_controls)
    external_comparison = compare_two_cell_external_results(
        external_template["rows"],
        config,
        probe_controls,
        tolerance=tolerance,
    )
    external_summary = external_comparison["summary"]
    tests.append(
        _validation_record(
            "external_result_gate",
            passed=bool(external_summary["externalComparisonReady"]),
            status="pass" if external_summary["externalComparisonReady"] else "incomplete",
            metrics={
                "rowCount": external_summary["rowCount"],
                "observedScalarCount": external_summary["observedScalarCount"],
                "missingExternalResultCount": external_summary["missingExternalResultCount"],
                "lockViolationCount": external_summary["lockViolationCount"],
            },
            evidence=["external results template and comparison gate"],
            limitations=[
                "No external MuJoCo/Gazebo/Isaac run or bench coordinate data has been loaded into this report."
            ],
            missing_evidence=list(external_summary["missingEvidence"]),
        )
    )

    external_test_ids = {"external_result_gate", "mujoco_proxy_run_gate"}
    internal_tests = [test for test in tests if test["testId"] not in external_test_ids]
    pass_count = sum(1 for test in tests if test["pass"])
    fail_count = sum(1 for test in tests if test["status"] == "fail")
    incomplete_count = sum(1 for test in tests if test["status"] == "incomplete")
    missing_evidence = sorted(
        {
            item
            for test in tests
            for item in test["missingEvidence"]
        }
    )
    internal_validation_pass = all(test["pass"] for test in internal_tests)
    external_validation_pass = bool(external_summary["externalComparisonReady"]) and bool(
        mujoco_run_summary["twoCellMujocoProxyRunComplete"]
    )
    if internal_validation_pass and external_validation_pass:
        overall_status = "validated-against-external-data"
    elif internal_validation_pass:
        overall_status = "needs-external-physical-validation"
    else:
        overall_status = "internal-validation-failed"
    return {
        "schema": TWO_CELL_PHYSICS_VALIDATION_SCHEMA,
        "model": "two-cell-cad-proxy-validation-gates",
        "cadReference": CAD_RAD_CELL_REFERENCE.to_dict(),
        "config": {
            "cellSize": config.cell_size,
            "backlash": config.backlash,
            "couplingGain": config.coupling_gain,
            "zCouplingGain": config.z_coupling_gain,
            "pinRadius": config.pin_radius,
            "holeRadius": config.hole_radius,
            "initialAlpha": config.initial_alpha,
        },
        "controls": {
            "alphaCommand": probe_controls.alpha_command,
            "zCommand": probe_controls.z_command,
            "holeSweepMax": probe_controls.hole_sweep_max,
            "holeSweepSteps": probe_controls.hole_sweep_steps,
        },
        "tests": tests,
        "summary": {
            "status": overall_status,
            "testCount": len(tests),
            "passCount": pass_count,
            "failCount": fail_count,
            "incompleteCount": incomplete_count,
            "internalValidationPass": internal_validation_pass,
            "externalValidationPass": external_validation_pass,
            "missingEvidence": missing_evidence,
            "missingEvidenceCount": len(missing_evidence),
        },
        "claimLabels": {
            "internalModel": "reduced two-cell proxy is numerically self-consistent for the listed gates",
            "cadReference": "A360 one-cell archive supplies dimensions and visual topology only",
            "physicalAccuracy": "not established until segmented CAD/contact simulation and bench data pass external_result_gate",
        },
    }


def export_two_cell_bench_json(config: LatticeConfig | None = None, controls: TwoCellBenchControls | None = None) -> str:
    return json.dumps(simulate_two_cell_bench(config, controls), indent=2)


def export_cad_rad_cell_layout_json(config: LatticeConfig | None = None, **layout_overrides: Any) -> str:
    return json.dumps(cad_rad_cell_layout(config, **layout_overrides), indent=2)


def export_cad_rad_cell_archive_audit_json(path: str | Path | None = None) -> str:
    return json.dumps(cad_rad_cell_archive_audit(path), indent=2)


def export_two_cell_segmented_cad_readiness_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    return json.dumps(two_cell_segmented_cad_readiness_report(config, controls), indent=2)


def export_two_cell_segmented_cad_readiness_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    report = two_cell_segmented_cad_readiness_report(config, controls)
    rows = []
    for section, key in (
        ("bodySegmentation", "bodySegmentationTasks"),
        ("jointAndContact", "jointAndContactTasks"),
        ("physicsParameter", "physicsParameterTasks"),
        ("engineHandoff", "engineHandoffTasks"),
    ):
        for task in report[key]:
            rows.append(
                {
                    "section": section,
                    "id": task.get("id", task.get("engine", "")),
                    "status": task.get("status", ""),
                    "neededEvidence": task.get("neededEvidence", task.get("neededAsset", "")),
                    "currentEvidence": task.get("currentEvidence", task.get("currentProxy", "")),
                    "expectedAssets": ";".join(asset["path"] for asset in task.get("expectedAssets", [])),
                    "detectedAssets": ";".join(asset["path"] for asset in task.get("detectedAssets", [])),
                }
            )
    fieldnames = [
        "section",
        "id",
        "status",
        "neededEvidence",
        "currentEvidence",
        "expectedAssets",
        "detectedAssets",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_segmented_cad_intake_templates_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
) -> str:
    return json.dumps(two_cell_segmented_cad_intake_templates(config, controls, intake_dir=intake_dir), indent=2)


def export_two_cell_segmented_cad_intake_templates_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
) -> str:
    templates = two_cell_segmented_cad_intake_templates(config, controls, intake_dir=intake_dir)
    rows = [
        {
            "path": item["path"],
            "kind": item["kind"],
            "description": item["description"],
            "bytes": len(str(item["content"]).encode("utf-8")),
        }
        for item in templates["templates"]
    ]
    fieldnames = ["path", "kind", "description", "bytes"]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_segmented_cad_intake_validation_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
) -> str:
    return json.dumps(
        two_cell_segmented_cad_intake_validation_report(config, controls, intake_dir=intake_dir),
        indent=2,
    )


def export_two_cell_segmented_cad_intake_validation_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    *,
    intake_dir: str = SEGMENTED_CAD_INTAKE_DIR,
) -> str:
    report = two_cell_segmented_cad_intake_validation_report(config, controls, intake_dir=intake_dir)
    rows = []
    for row in report["templateValidationRows"]:
        rows.append(
            {
                "section": "template",
                "id": row["templateId"],
                "path": row["path"],
                "status": row["status"],
                "requiredFieldCount": row["requiredFieldCount"],
                "filledFieldCount": row["filledFieldCount"],
                "errorCount": row["errorCount"],
                "errors": ";".join(row["errors"]),
            }
        )
    for row in report["bodyAssetValidationRows"]:
        rows.append(
            {
                "section": "bodyAsset",
                "id": row["assetId"],
                "path": ";".join(asset["path"] for asset in row["expectedAssets"]),
                "status": row["status"],
                "requiredFieldCount": "",
                "filledFieldCount": "",
                "errorCount": 0 if row["status"] == "complete" else 1,
                "errors": "" if row["status"] == "complete" else "segmented body asset missing",
            }
        )
    rows.append(
        {
            "section": "radialSurfaceAsset",
            "id": "radial_pad_hole_surfaces",
            "path": ";".join(asset["path"] for asset in report["radialSurfaceAssetValidation"]["expectedAssets"]),
            "status": report["radialSurfaceAssetValidation"]["status"],
            "requiredFieldCount": "",
            "filledFieldCount": "",
            "errorCount": 0 if report["radialSurfaceAssetValidation"]["status"] == "complete" else 1,
            "errors": ""
            if report["radialSurfaceAssetValidation"]["status"] == "complete"
            else "hole surface asset or filled surface json missing",
        }
    )
    for row in report["engineHandoffValidationRows"]:
        rows.append(
            {
                "section": "engineHandoff",
                "id": row["engine"],
                "path": ";".join(asset["path"] for asset in row["expectedAssets"]),
                "status": row["status"],
                "requiredFieldCount": "",
                "filledFieldCount": "",
                "errorCount": 0 if row["status"] == "ready-for-external-run" else 1,
                "errors": "" if row["status"] == "ready-for-external-run" else "engine handoff asset missing",
            }
        )
    fieldnames = [
        "section",
        "id",
        "path",
        "status",
        "requiredFieldCount",
        "filledFieldCount",
        "errorCount",
        "errors",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_backlash_sweep_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    sweep = sweep_two_cell_backlash(config, controls)
    fieldnames = [
        "index",
        "holeRadius",
        "clearance",
        "effectiveBacklash",
        "alphaDeadZone",
        "verticalDeadZone",
        "leftResidualZ",
        "rightZ",
        "axialStrain",
        "contactMode",
        "totalEnergy",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(sweep["rows"])
    return buffer.getvalue()


def export_two_cell_actuation_sweep_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **sweep_overrides: Any,
) -> str:
    return json.dumps(two_cell_actuation_sweep(config, controls, **sweep_overrides), indent=2)


def export_two_cell_actuation_sweep_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **sweep_overrides: Any,
) -> str:
    sweep = two_cell_actuation_sweep(config, controls, **sweep_overrides)
    fieldnames = [
        "index",
        "lockMode",
        "alphaCommand",
        "zCommand",
        "pinRadius",
        "holeRadius",
        "clearance",
        "effectiveBacklash",
        "contactMode",
        "leftX",
        "leftY",
        "leftZ",
        "rightX",
        "rightY",
        "rightZ",
        "leftAlpha",
        "rightAlpha",
        "leftAlphaDeltaFromInitial",
        "rightAlphaDeltaFromInitial",
        "leftTheta",
        "rightTheta",
        "verticalShear",
        "axialStrain",
        "contactPenetration",
        "totalEnergy",
        "leftFixtureError",
        "rightFixtureError",
        "rightStateMotionSuppressed",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(sweep["rows"])
    return buffer.getvalue()


def export_two_cell_mjcf_proxy_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **report_overrides: Any,
) -> str:
    return json.dumps(two_cell_mjcf_proxy_report(config, controls, **report_overrides), indent=2)


def export_two_cell_mjcf_xml(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **report_overrides: Any,
) -> str:
    return str(two_cell_mjcf_proxy_report(config, controls, **report_overrides)["xml"])


def export_two_cell_mujoco_proxy_run_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **run_overrides: Any,
) -> str:
    return json.dumps(two_cell_mujoco_proxy_run_report(config, controls, **run_overrides), indent=2)


def export_two_cell_external_fidelity_mjcf_run_json(
    mjcf_index_path: str | Path | None = None,
    **run_overrides: Any,
) -> str:
    return json.dumps(two_cell_external_fidelity_mjcf_run_report(mjcf_index_path, **run_overrides), indent=2)


def two_cell_external_fidelity_mjcf_measurement_rows(report: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in report.get("results", {}).get("measurementRows", [])
        if isinstance(row, dict)
    ]


def two_cell_external_fidelity_mjcf_measurements_csv(report: dict[str, Any]) -> str:
    rows = two_cell_external_fidelity_mjcf_measurement_rows(report)
    fieldnames = [
        "caseId",
        "matrixIndex",
        "connector",
        "leftSite",
        "rightSite",
        "actuationCase",
        "lockMode",
        "backlash",
        "alphaCommand",
        "zCommand",
        "gravityForce",
        "pinRadius",
        "holeRadius",
        "pinHoleClearanceMm",
        "observedLeftCellX",
        "observedLeftCellY",
        "observedLeftCellZ",
        "observedRightCellX",
        "observedRightCellY",
        "observedRightCellZ",
        "observedLeftAlpha",
        "observedRightAlpha",
        "observedLeftTheta",
        "observedRightTheta",
        "observedLeftXmm",
        "observedLeftYmm",
        "observedLeftZmm",
        "observedRightXmm",
        "observedRightYmm",
        "observedRightZmm",
        "observedLateralSlipMm",
        "observedVerticalSlipMm",
        "observedTotalSlipMm",
        "observedContactMode",
        "lockHeldObserved",
        "measurementSource",
        "notes",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row.get(field, "") for field in fieldnames} for row in rows)
    return buffer.getvalue()


def export_two_cell_external_fidelity_mjcf_measurements_csv(
    mjcf_index_path: str | Path | None = None,
    **run_overrides: Any,
) -> str:
    report = two_cell_external_fidelity_mjcf_run_report(mjcf_index_path, **run_overrides)
    return two_cell_external_fidelity_mjcf_measurements_csv(report)


def export_two_cell_quasistatic_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **solve_overrides: Any,
) -> str:
    return json.dumps(solve_two_cell_quasistatic(config, controls, **solve_overrides), indent=2)


def export_two_cell_quasistatic_sweep_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **sweep_overrides: Any,
) -> str:
    return json.dumps(sweep_two_cell_quasistatic(config, controls, **sweep_overrides), indent=2)


def export_two_cell_quasistatic_sweep_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **sweep_overrides: Any,
) -> str:
    sweep = sweep_two_cell_quasistatic(config, controls, **sweep_overrides)
    fieldnames = [
        "index",
        "lockMode",
        "alphaCommand",
        "zCommand",
        "pinRadius",
        "holeRadius",
        "clearance",
        "contactMode",
        "solverSuccess",
        "leftX",
        "leftZ",
        "leftAlpha",
        "rightX",
        "rightZ",
        "rightAlpha",
        "verticalShear",
        "verticalExcess",
        "axialExcess",
        "storedEnergy",
        "totalEnergy",
        "leftFixtureError",
        "rightFixtureError",
        "rightStateAlphaError",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(sweep["rows"])
    return buffer.getvalue()


def export_two_cell_connector_contact_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **report_overrides: Any,
) -> str:
    return json.dumps(two_cell_connector_contact_report(config, controls, **report_overrides), indent=2)


def export_two_cell_connector_contact_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **report_overrides: Any,
) -> str:
    report = two_cell_connector_contact_report(config, controls, **report_overrides)
    fieldnames = [
        "connector",
        "leftSite",
        "rightSite",
        "lateralClearanceMm",
        "verticalFreePlayMm",
        "lateralSlipMm",
        "verticalSlipMm",
        "totalSlipMm",
        "lateralExcessMm",
        "verticalExcessMm",
        "contactPenalty",
        "contactMode",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(
        {
            "connector": row["connector"],
            "leftSite": row["leftSite"],
            "rightSite": row["rightSite"],
            "lateralClearanceMm": row["lateralClearanceMm"],
            "verticalFreePlayMm": row["verticalFreePlayMm"],
            "lateralSlipMm": row["lateralSlipMm"],
            "verticalSlipMm": row["verticalSlipMm"],
            "totalSlipMm": row["totalSlipMm"],
            "lateralExcessMm": row["lateralExcessMm"],
            "verticalExcessMm": row["verticalExcessMm"],
            "contactPenalty": row["contactPenalty"],
            "contactMode": row["contactMode"],
        }
        for row in report["connectors"]
    )
    return buffer.getvalue()


def export_two_cell_connector_contact_sweep_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **sweep_overrides: Any,
) -> str:
    return json.dumps(sweep_two_cell_connector_contact(config, controls, **sweep_overrides), indent=2)


def export_two_cell_connector_contact_sweep_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **sweep_overrides: Any,
) -> str:
    sweep = sweep_two_cell_connector_contact(config, controls, **sweep_overrides)
    fieldnames = [
        "index",
        "lockMode",
        "alphaCommand",
        "zCommand",
        "pinRadius",
        "holeRadius",
        "clearance",
        "pinHoleClearanceMm",
        "connectorPitchMm",
        "connectorCount",
        "activeContactCount",
        "maxLateralSlipMm",
        "maxVerticalSlipMm",
        "maxLateralExcessMm",
        "maxVerticalExcessMm",
        "totalContactPenalty",
        "solverSuccess",
        "positionLocksHold",
        "stateLocksHold",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(sweep["rows"])
    return buffer.getvalue()


def export_two_cell_physical_simulation_suite_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **suite_overrides: Any,
) -> str:
    return json.dumps(two_cell_physical_simulation_suite(config, controls, **suite_overrides), indent=2)


def export_two_cell_physical_simulation_suite_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **suite_overrides: Any,
) -> str:
    suite = two_cell_physical_simulation_suite(config, controls, **suite_overrides)
    fieldnames = [
        "index",
        "caseId",
        "lockMode",
        "backlash",
        "alphaCommand",
        "zCommand",
        "gravityForce",
        "pinRadius",
        "holeRadius",
        "clearance",
        "effectiveBacklash",
        "pinHoleClearanceMm",
        "solverSuccess",
        "lockPass",
        "positionLocksHold",
        "stateLocksHold",
        "rightLockedVerticalMotionAllowed",
        "leftX",
        "leftZ",
        "rightX",
        "rightZ",
        "leftAlpha",
        "rightAlpha",
        "leftAlphaDeltaFromInitial",
        "rightAlphaDeltaFromInitial",
        "leftTheta",
        "rightTheta",
        "benchLeftResidualZ",
        "verticalShear",
        "verticalExcess",
        "axialExcess",
        "contactMode",
        "connectorActiveContactCount",
        "connectorMaxVerticalSlipMm",
        "connectorMaxVerticalExcessMm",
        "connectorTotalContactPenalty",
        "totalEnergy",
        "externalMeasurementRequired",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row[field] for field in fieldnames} for row in suite["rows"])
    return buffer.getvalue()


def export_two_cell_physical_fidelity_matrix_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **matrix_overrides: Any,
) -> str:
    return json.dumps(two_cell_physical_fidelity_matrix(config, controls, **matrix_overrides), indent=2)


def export_two_cell_physical_fidelity_matrix_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **matrix_overrides: Any,
) -> str:
    matrix = two_cell_physical_fidelity_matrix(config, controls, **matrix_overrides)
    fieldnames = [
        "index",
        "actuationCase",
        "lockMode",
        "backlash",
        "alphaCommand",
        "zCommand",
        "gravityForce",
        "pinRadius",
        "holeRadius",
        "clearance",
        "effectiveBacklash",
        "pinHoleClearanceMm",
        "solverSuccess",
        "lockPass",
        "positionLocksHold",
        "stateLocksHold",
        "rightLockedVerticalMotionAllowed",
        "leftX",
        "leftZ",
        "rightX",
        "rightZ",
        "leftAlpha",
        "rightAlpha",
        "leftAlphaDeltaFromInitial",
        "rightAlphaDeltaFromInitial",
        "leftTheta",
        "rightTheta",
        "verticalShear",
        "verticalExcess",
        "axialExcess",
        "contactMode",
        "connectorActiveContactCount",
        "connectorInsideClearanceCount",
        "connectorMaxLateralSlipMm",
        "connectorMaxVerticalSlipMm",
        "connectorMaxLateralExcessMm",
        "connectorMaxVerticalExcessMm",
        "connectorTotalContactPenalty",
        "totalEnergy",
        "rightFixtureError",
        "rightStateAlphaError",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row[field] for field in fieldnames} for row in matrix["rows"])
    return buffer.getvalue()


def export_two_cell_physical_response_atlas_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **atlas_overrides: Any,
) -> str:
    return json.dumps(two_cell_physical_response_atlas(config, controls, **atlas_overrides), indent=2)


def export_two_cell_physical_response_atlas_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **atlas_overrides: Any,
) -> str:
    atlas = two_cell_physical_response_atlas(config, controls, **atlas_overrides)
    fieldnames = [
        "category",
        "rank",
        "caseId",
        "matrixIndex",
        "index",
        "actuationCase",
        "lockMode",
        "backlash",
        "holeRadius",
        "pinHoleClearanceMm",
        "alphaCommand",
        "zCommand",
        "rightZ",
        "leftZ",
        "rightAlpha",
        "rightAlphaDeltaFromInitial",
        "connectorMaxVerticalSlipMm",
        "connectorMaxVerticalExcessMm",
        "connectorTotalContactPenalty",
        "totalEnergy",
        "lockPass",
        "positionLocksHold",
        "stateLocksHold",
        "rightLockedVerticalMotionAllowed",
    ]
    rows = []
    for category, cases in atlas["benchPriority"].items():
        for rank, case in enumerate(cases, start=1):
            rows.append({"category": category, "rank": rank, **case})
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row.get(field, "") for field in fieldnames} for row in rows)
    return buffer.getvalue()


def export_two_cell_contact_phase_map_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **phase_map_overrides: Any,
) -> str:
    return json.dumps(two_cell_contact_phase_map(config, controls, **phase_map_overrides), indent=2)


def export_two_cell_contact_phase_map_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **phase_map_overrides: Any,
) -> str:
    phase_map = two_cell_contact_phase_map(config, controls, **phase_map_overrides)
    fieldnames = [
        "caseId",
        "index",
        "phase",
        "actuationCase",
        "lockMode",
        "backlash",
        "alphaCommand",
        "zCommand",
        "gravityForce",
        "pinRadius",
        "holeRadius",
        "pinHoleClearanceMm",
        "leftZ",
        "rightZ",
        "leftAlphaDeltaFromInitial",
        "rightAlphaDeltaFromInitial",
        "verticalShear",
        "verticalExcess",
        "axialExcess",
        "contactMode",
        "connectorActiveContactCount",
        "connectorMaxLateralExcessMm",
        "connectorMaxVerticalSlipMm",
        "connectorMaxVerticalExcessMm",
        "connectorTotalContactPenalty",
        "totalEnergy",
        "lockPass",
        "positionLocksHold",
        "stateLocksHold",
        "rightLockedVerticalMotionAllowed",
        "axialContact",
        "verticalContact",
        "lockPhase",
        "freePlayPhase",
        "gravitySagPhase",
        "physicalAccuracyValidated",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row.get(field, "") for field in fieldnames} for row in phase_map["rows"])
    return buffer.getvalue()


def export_two_cell_radius_backlash_phase_diagram_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **diagram_overrides: Any,
) -> str:
    return json.dumps(two_cell_radius_backlash_phase_diagram(config, controls, **diagram_overrides), indent=2)


def export_two_cell_radius_backlash_phase_diagram_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **diagram_overrides: Any,
) -> str:
    diagram = two_cell_radius_backlash_phase_diagram(config, controls, **diagram_overrides)
    fieldnames = [
        "caseId",
        "matrixIndex",
        "phase",
        "actuationCase",
        "lockMode",
        "backlash",
        "holeRadius",
        "pinHoleClearanceMm",
        "effectiveBacklash",
        "effectiveHoleRadius",
        "effectivePinHoleClearanceMm",
        "alphaCommand",
        "zCommand",
        "gravityForce",
        "leftZ",
        "rightZ",
        "rightAlpha",
        "contactMode",
        "connectorActiveContactCount",
        "connectorMaxLateralExcessMm",
        "connectorMaxVerticalSlipMm",
        "connectorMaxVerticalExcessMm",
        "connectorTotalContactPenalty",
        "totalEnergy",
        "axialContact",
        "verticalContact",
        "lockPhase",
        "freePlayPhase",
        "gravitySagPhase",
        "physicalAccuracyValidated",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row.get(field, "") for field in fieldnames} for row in diagram["rows"])
    return buffer.getvalue()


def export_two_cell_radius_backlash_transition_report_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **report_overrides: Any,
) -> str:
    return json.dumps(two_cell_radius_backlash_transition_report(config, controls, **report_overrides), indent=2)


def export_two_cell_radius_backlash_transition_report_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **report_overrides: Any,
) -> str:
    report = two_cell_radius_backlash_transition_report(config, controls, **report_overrides)
    fieldnames = [
        "rank",
        "category",
        "phaseDiagramCaseId",
        "caseId",
        "matrixIndex",
        "phase",
        "actuationCase",
        "lockMode",
        "backlash",
        "holeRadius",
        "pinHoleClearanceMm",
        "effectiveBacklash",
        "effectiveHoleRadius",
        "effectivePinHoleClearanceMm",
        "alphaCommand",
        "zCommand",
        "gravityForce",
        "leftZ",
        "rightZ",
        "rightAlpha",
        "contactMode",
        "connectorActiveContactCount",
        "connectorMaxVerticalSlipMm",
        "connectorMaxVerticalExcessMm",
        "connectorTotalContactPenalty",
        "lockPass",
        "measurementTarget",
        "observedPhase",
        "observedRightZ",
        "observedRightAlpha",
        "observedConnectorMaxVerticalSlipMm",
        "observedConnectorMaxVerticalExcessMm",
        "observedLockHeld",
        "notes",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row.get(field, "") for field in fieldnames} for row in report["measurementCases"])
    return buffer.getvalue()


def export_two_cell_fidelity_matrix_measurement_template_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **matrix_overrides: Any,
) -> str:
    return json.dumps(
        two_cell_fidelity_matrix_measurement_template(config, controls, **matrix_overrides),
        indent=2,
    )


def export_two_cell_fidelity_matrix_measurement_template_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **matrix_overrides: Any,
) -> str:
    rows = two_cell_fidelity_matrix_measurement_template_rows(config, controls, **matrix_overrides)
    fieldnames = list(rows[0]) if rows else []
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_external_fidelity_matrix_manifest_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **manifest_overrides: Any,
) -> str:
    return json.dumps(
        two_cell_external_fidelity_matrix_manifest(config, controls, **manifest_overrides),
        indent=2,
    )


def export_two_cell_external_fidelity_matrix_manifest_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **matrix_overrides: Any,
) -> str:
    rows = two_cell_external_fidelity_matrix_manifest_rows(config, controls, **matrix_overrides)
    fieldnames = [
        "caseId",
        "matrixIndex",
        "actuationCase",
        "lockMode",
        "backlash",
        "alphaCommand",
        "zCommand",
        "gravityForce",
        "pinRadius",
        "holeRadius",
        "pinHoleClearanceMm",
        "leftPositionLocked",
        "rightStateLocked",
        "rightPositionLocked",
        "expectedLeftX",
        "expectedLeftY",
        "expectedLeftZ",
        "expectedRightX",
        "expectedRightY",
        "expectedRightZ",
        "expectedLeftAlpha",
        "expectedRightAlpha",
        "expectedLeftTheta",
        "expectedRightTheta",
        "expectedContactMode",
        "expectedConnectorMaxLateralSlipMm",
        "expectedConnectorMaxVerticalSlipMm",
        "expectedConnectorMaxVerticalExcessMm",
        "expectedConnectorTotalContactPenalty",
        "mjcfProxyPath",
        "engineResultPath",
        "measurementRows",
        "runStatus",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_cad_contact_decomposition_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **spec_overrides: Any,
) -> str:
    return json.dumps(two_cell_cad_contact_decomposition_spec(config, controls, **spec_overrides), indent=2)


def export_two_cell_cad_contact_decomposition_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **spec_overrides: Any,
) -> str:
    spec = two_cell_cad_contact_decomposition_spec(config, controls, **spec_overrides)
    rows = []
    for role in spec["bodyRoles"]:
        rows.append(
            {
                "section": "bodyRole",
                "cell": "",
                "connector": "",
                "site": "",
                "segmentIndex": "",
                "angleStartDeg": "",
                "angleEndDeg": "",
                "bodyA": ";".join(role["instances"]),
                "bodyB": "",
                "primitive": "",
                "centerXmm": "",
                "centerYmm": "",
                "centerZmm": "",
                "axisX": "",
                "axisY": "",
                "axisZ": "",
                "pinRadiusMm": "",
                "holeRadiusMm": "",
                "clearanceMm": "",
                "requiredAsset": ";".join(role["expectedAssets"]),
                "source": role["solverUse"],
                "physicalAccuracyValidated": False,
            }
        )
    for hole in spec["holeDecomposition"]:
        for segment in hole["ringSegments"]:
            center = segment["centerMm"]
            axis = segment["axis"]
            rows.append(
                {
                    "section": "holeSegment",
                    "cell": hole["cell"],
                    "connector": "",
                    "site": hole["site"],
                    "segmentIndex": segment["segmentIndex"],
                    "angleStartDeg": segment["angleStartDeg"],
                    "angleEndDeg": segment["angleEndDeg"],
                    "bodyA": segment["bodyA"],
                    "bodyB": segment["bodyB"],
                    "primitive": segment["primitive"],
                    "centerXmm": center[0],
                    "centerYmm": center[1],
                    "centerZmm": center[2],
                    "axisX": axis[0],
                    "axisY": axis[1],
                    "axisZ": axis[2],
                    "pinRadiusMm": hole["pinRadiusMm"],
                    "holeRadiusMm": hole["holeRadiusMm"],
                    "clearanceMm": hole["clearanceMm"],
                    "requiredAsset": "radial_pad_hole_surfaces.json",
                    "source": segment["source"],
                    "physicalAccuracyValidated": segment["physicalAccuracyValidated"],
                }
            )
    for pair in spec["twoCellContactPairs"]:
        left_center = pair["leftCenterMm"]
        axis = pair["axis"]
        rows.append(
            {
                "section": "connectorPair",
                "cell": f"{pair['leftCell']}-{pair['rightCell']}",
                "connector": pair["name"],
                "site": f"{pair['leftSite']}-{pair['rightSite']}",
                "segmentIndex": "",
                "angleStartDeg": "",
                "angleEndDeg": "",
                "bodyA": pair["leftPinBody"],
                "bodyB": pair["rightHoleBody"],
                "primitive": "explicit_pin_hole_connector_contact",
                "centerXmm": left_center[0],
                "centerYmm": left_center[1],
                "centerZmm": left_center[2],
                "axisX": axis[0],
                "axisY": axis[1],
                "axisZ": axis[2],
                "pinRadiusMm": pair["pinRadiusMm"],
                "holeRadiusMm": pair["holeRadiusMm"],
                "clearanceMm": pair["clearanceMm"],
                "requiredAsset": pair["requiredSurfacePair"],
                "source": pair["source"],
                "physicalAccuracyValidated": pair["physicalAccuracyValidated"],
            }
        )
    fieldnames = [
        "section",
        "cell",
        "connector",
        "site",
        "segmentIndex",
        "angleStartDeg",
        "angleEndDeg",
        "bodyA",
        "bodyB",
        "primitive",
        "centerXmm",
        "centerYmm",
        "centerZmm",
        "axisX",
        "axisY",
        "axisZ",
        "pinRadiusMm",
        "holeRadiusMm",
        "clearanceMm",
        "requiredAsset",
        "source",
        "physicalAccuracyValidated",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_exact_contact_handoff_plan_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **plan_overrides: Any,
) -> str:
    return json.dumps(two_cell_exact_contact_handoff_plan(config, controls, **plan_overrides), indent=2)


def export_two_cell_exact_contact_handoff_plan_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **plan_overrides: Any,
) -> str:
    plan = two_cell_exact_contact_handoff_plan(config, controls, **plan_overrides)
    fieldnames = [
        "rank",
        "caseId",
        "source",
        "category",
        "phase",
        "actuationCase",
        "lockMode",
        "backlash",
        "holeRadius",
        "pinHoleClearanceMm",
        "effectiveBacklash",
        "effectiveHoleRadius",
        "effectivePinHoleClearanceMm",
        "alphaCommand",
        "zCommand",
        "gravityForce",
        "leftZ",
        "rightZ",
        "rightAlpha",
        "connectorActiveContactCount",
        "connectorMaxVerticalSlipMm",
        "connectorMaxVerticalExcessMm",
        "connectorTotalContactPenalty",
        "leftPositionLocked",
        "rightStateLocked",
        "rightPositionLocked",
        "expectedConnectorCount",
        "requiredCadAssets",
        "requiredMeasuredInputs",
        "engineHandoffReady",
        "canRunExactContact",
        "physicalAccuracyValidated",
        "notes",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    for row in plan["caseRows"]:
        writer.writerow(
            {
                field: ";".join(row.get(field, []))
                if field in {"requiredCadAssets", "requiredMeasuredInputs"}
                else row.get(field, "")
                for field in fieldnames
            }
        )
    return buffer.getvalue()


def export_two_cell_physical_test_packet_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    return json.dumps(two_cell_physical_test_packet(config, controls), indent=2)


def export_two_cell_physics_validation_report_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **report_overrides: Any,
) -> str:
    return json.dumps(two_cell_physics_validation_report(config, controls, **report_overrides), indent=2)


def export_two_cell_physics_validation_report_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **report_overrides: Any,
) -> str:
    report = two_cell_physics_validation_report(config, controls, **report_overrides)
    rows = [
        {
            "testId": test["testId"],
            "status": test["status"],
            "pass": test["pass"],
            "metrics": json.dumps(test["metrics"], sort_keys=True, separators=(",", ":")),
            "missingEvidence": ";".join(test["missingEvidence"]),
            "limitations": ";".join(test["limitations"]),
        }
        for test in report["tests"]
    ]
    fieldnames = ["testId", "status", "pass", "metrics", "missingEvidence", "limitations"]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_measurement_template_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    rows = two_cell_measurement_template_rows(config, controls)
    fieldnames = [
        "caseId",
        "alphaCommand",
        "zCommand",
        "pinRadius",
        "holeRadius",
        "leftX",
        "leftY",
        "leftZ",
        "rightX",
        "rightY",
        "rightZ",
        "rightAlpha",
        "rightTheta",
        "contactModeObserved",
        "notes",
        "predictedLeftX",
        "predictedLeftY",
        "predictedLeftZ",
        "predictedRightX",
        "predictedRightY",
        "predictedRightZ",
        "predictedRightAlpha",
        "predictedRightTheta",
        "contactModePredicted",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_connector_measurement_template_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    return json.dumps(two_cell_connector_measurement_template(config, controls), indent=2)


def export_two_cell_connector_measurement_template_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    rows = two_cell_connector_measurement_template_rows(config, controls)
    fieldnames = [
        "caseId",
        "connector",
        "leftSite",
        "rightSite",
        "lockMode",
        "alphaCommand",
        "zCommand",
        "backlash",
        "pinRadius",
        "holeRadius",
        "pinHoleClearanceMm",
        "predictedLeftXmm",
        "predictedLeftYmm",
        "predictedLeftZmm",
        "predictedRightXmm",
        "predictedRightYmm",
        "predictedRightZmm",
        "predictedSlipXmm",
        "predictedSlipYmm",
        "predictedSlipZmm",
        "predictedLateralSlipMm",
        "predictedVerticalSlipMm",
        "predictedTotalSlipMm",
        "predictedLateralExcessMm",
        "predictedVerticalExcessMm",
        "predictedContactMode",
        "observedLeftXmm",
        "observedLeftYmm",
        "observedLeftZmm",
        "observedRightXmm",
        "observedRightYmm",
        "observedRightZmm",
        "observedLateralSlipMm",
        "observedVerticalSlipMm",
        "observedTotalSlipMm",
        "observedContactMode",
        "lockHeldObserved",
        "measurementSource",
        "notes",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_external_results_template_json(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    return json.dumps(two_cell_external_results_template(config, controls), indent=2)


def export_two_cell_external_results_template_csv(
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
) -> str:
    rows = two_cell_external_results_template_rows(config, controls)
    fieldnames = [
        "caseId",
        "bodyId",
        "bodyName",
        "lockMode",
        "alphaCommand",
        "zCommand",
        "pinRadius",
        "holeRadius",
        "expectedX",
        "expectedY",
        "expectedZ",
        "expectedAlpha",
        "expectedTheta",
        "externalX",
        "externalY",
        "externalZ",
        "externalAlpha",
        "externalTheta",
        "sourceEngine",
        "notes",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def export_two_cell_external_comparison_json(
    external_results: Iterable[dict[str, Any]],
    config: LatticeConfig | None = None,
    controls: TwoCellBenchControls | None = None,
    **comparison_overrides: Any,
) -> str:
    return json.dumps(
        compare_two_cell_external_results(
            external_results,
            config,
            controls,
            **comparison_overrides,
        ),
        indent=2,
    )
