from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


ProvenanceStatus = Literal[
    "paper-supported",
    "implementation-assumption",
    "simulator-diagnostic",
    "calibration-gap",
]


@dataclass(frozen=True)
class ModelProvenanceItem:
    id: str
    label: str
    status: ProvenanceStatus
    source: str
    page_refs: tuple[str, ...]
    evidence: str
    implementation: str
    limitation: str


MODEL_PROVENANCE: tuple[ModelProvenanceItem, ...] = (
    ModelProvenanceItem(
        id="rotating_square_dilation",
        label="Rotating-square dilation",
        status="paper-supported",
        source="RLABS proceedings and RAD preprint",
        page_refs=("3002882.pdf pp. 3-4", "RAD preprint pp. 2, 5"),
        evidence="Cells are modeled as rotating-square auxetic mechanisms with a spatial alpha(x, y) dilation field.",
        implementation="alpha_to_theta, build_paper_rad_cell_geometry, browser Paper RAD cell view",
        limitation="The exact CAD travel curve still needs hardware measurement.",
    ),
    ModelProvenanceItem(
        id="alpha_theta_law",
        label="Theta-alpha law",
        status="paper-supported",
        source="RAD preprint",
        page_refs=("RAD preprint p. 5 Eq. 3",),
        evidence="theta[degrees] = 70 * alpha - 60 in the simplified coupling model.",
        implementation="alpha_to_theta and theta_to_alpha",
        limitation="Use as the v1 simplified model, not as a universal calibrated CAD relation.",
    ),
    ModelProvenanceItem(
        id="backlash_dead_zone",
        label="Backlash dead zone",
        status="paper-supported",
        source="RAD preprint and RLABS proceedings",
        page_refs=("RAD preprint p. 5 Eq. 2", "3002882.pdf pp. 6-7"),
        evidence="f(x) = max(0, x - b) + min(x + b, 0) with zero response inside the backlash gap.",
        implementation="backlash_activation, alpha_backlash_operator, browser RAD.backlashActivation",
        limitation="Friction, wear, and asymmetric clearances are not included yet.",
    ),
    ModelProvenanceItem(
        id="paper_rad_cell_topology",
        label="Paper RAD cell topology",
        status="paper-supported",
        source="RAD preprint",
        page_refs=("RAD preprint p. 6",),
        evidence="The CAD unit cell has two concentric parts with four joints each, b = 0.1, 35 mm side length, and target nu = -0.4.",
        implementation="PaperRADReference, build_paper_rad_cell_geometry, browser paperRad visual mode",
        limitation="Plate thickness, bosses, screw hardware, and exact hole dimensions remain configurable.",
    ),
    ModelProvenanceItem(
        id="vertical_residual_coupling",
        label="Vertical residual coupling",
        status="implementation-assumption",
        source="User hardware observation plus RAD pin-hole clearance context",
        page_refs=("RAD preprint p. 6",),
        evidence="The papers support backlash and hole tolerance, but do not give a calibrated z die-off law.",
        implementation="z_coupling_gain, vertical_clearance_operator, z residual overlay",
        limitation="Requires physical measurement of pin/hole clearance, load, and neighbor-state dependence.",
    ),
    ModelProvenanceItem(
        id="spring_hinge_preview",
        label="Spring-hinge physical preview",
        status="implementation-assumption",
        source="Reduced-order mechanics model",
        page_refs=("programmable_mechanics_monograph p. 52",),
        evidence="The validation ladder calls for spring, hinge, frame, and auxetic checks before higher-fidelity mechanics.",
        implementation="solve_spring_hinge_3d and browser springPreview mode",
        limitation="Not a finite-element or calibrated rigid-body simulation.",
    ),
    ModelProvenanceItem(
        id="operator_diagnostics",
        label="Operator interaction diagnostics",
        status="simulator-diagnostic",
        source="Programmable-discontinuity framework",
        page_refs=("programmable_mechanics_monograph pp. 23, 42, 55",),
        evidence="Event operators should be studied for composition, commutation, locality, and reachable behavior.",
        implementation="pairwise interaction residuals, hotspot maps, degree maps, event-order checks",
        limitation="Diagnostic values reveal non-additivity; they are not material constitutive laws.",
    ),
    ModelProvenanceItem(
        id="inverse_design_linearization",
        label="Linearized inverse design",
        status="simulator-diagnostic",
        source="Current numerical architecture",
        page_refs=("programmable_mechanics_monograph p. 23", "programmable_mechanics_monograph p. 52"),
        evidence="Compiling event sequences or commands from a target shape is an open research question.",
        implementation="build_response_matrix, solve_inverse_design, browser Solve Linear Fit",
        limitation="Command plans remain proposals until checked against physical preview or measured response.",
    ),
    ModelProvenanceItem(
        id="calibrated_cad_dimensions",
        label="Calibrated CAD dimensions",
        status="calibration-gap",
        source="Physical RAD hardware",
        page_refs=("RAD preprint p. 6",),
        evidence="The paper gives prototype scale and fabrication tolerance, but not all fabrication part dimensions.",
        implementation="Future calibrated CAD-like visual and solver parameters",
        limitation="Measure pin radius, hole radius, plate thickness, boss dimensions, joint stack height, and friction.",
    ),
)


def model_provenance() -> tuple[ModelProvenanceItem, ...]:
    return MODEL_PROVENANCE


def provenance_by_status(status: ProvenanceStatus) -> tuple[ModelProvenanceItem, ...]:
    return tuple(item for item in MODEL_PROVENANCE if item.status == status)


def provenance_summary() -> dict[str, int]:
    summary: dict[str, int] = {
        "paper-supported": 0,
        "implementation-assumption": 0,
        "simulator-diagnostic": 0,
        "calibration-gap": 0,
    }
    for item in MODEL_PROVENANCE:
        summary[item.status] += 1
    return summary
