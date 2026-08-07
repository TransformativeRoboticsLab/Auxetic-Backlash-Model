(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  const MODEL_PROVENANCE = Object.freeze([
    {
      id: "rotating_square_dilation",
      label: "Rotating-square dilation",
      status: "paper-supported",
      source: "RLABS proceedings and RAD preprint",
      pageRefs: ["3002882.pdf pp. 3-4", "RAD preprint pp. 2, 5"],
      evidence: "Cells are rotating-square auxetic mechanisms with an alpha(x, y) dilation field.",
      implementation: "alpha/theta math, Paper RAD cell view",
      limitation: "Exact CAD travel curve still needs hardware measurement.",
    },
    {
      id: "alpha_theta_law",
      label: "Theta-alpha law",
      status: "paper-supported",
      source: "RAD preprint",
      pageRefs: ["RAD preprint p. 5 Eq. 3"],
      evidence: "theta[degrees] = 70 * alpha - 60 in the simplified coupling model.",
      implementation: "RAD.alphaToTheta and RAD.thetaToAlpha",
      limitation: "Use as v1 simplified model, not a universal calibrated CAD relation.",
    },
    {
      id: "backlash_dead_zone",
      label: "Backlash dead zone",
      status: "paper-supported",
      source: "RAD preprint and RLABS proceedings",
      pageRefs: ["RAD preprint p. 5 Eq. 2", "3002882.pdf pp. 6-7"],
      evidence: "f(x) = max(0, x - b) + min(x + b, 0) with zero response inside backlash.",
      implementation: "RAD.backlashActivation and operator propagation",
      limitation: "Friction, wear, and asymmetric clearances are not included yet.",
    },
    {
      id: "paper_rad_cell_topology",
      label: "Paper RAD cell topology",
      status: "paper-supported",
      source: "RAD preprint",
      pageRefs: ["RAD preprint p. 6"],
      evidence: "Two concentric parts with four joints each, b = 0.1, 35 mm side length, target nu = -0.4.",
      implementation: "Paper RAD cell visual mode and mesh export",
      limitation: "Exact plate, boss, screw, and hole dimensions remain configurable.",
    },
    {
      id: "vertical_residual_coupling",
      label: "Vertical residual coupling",
      status: "implementation-assumption",
      source: "User hardware observation plus pin-hole clearance context",
      pageRefs: ["RAD preprint p. 6"],
      evidence: "Papers support backlash and hole tolerance, but not a calibrated z die-off law.",
      implementation: "z coupling, vertical residual overlay",
      limitation: "Requires measurement of clearance, load, and neighbor-state dependence.",
    },
    {
      id: "spring_hinge_preview",
      label: "Spring-hinge physical preview",
      status: "implementation-assumption",
      source: "Reduced-order mechanics model",
      pageRefs: ["monograph p. 52"],
      evidence: "Validation ladder calls for spring, hinge, frame, and auxetic checks.",
      implementation: "browser springPreview mode and Python spring-hinge solver",
      limitation: "Not a finite-element or calibrated rigid-body simulation.",
    },
    {
      id: "operator_diagnostics",
      label: "Operator diagnostics",
      status: "simulator-diagnostic",
      source: "Programmable-discontinuity framework",
      pageRefs: ["monograph pp. 23, 42, 55"],
      evidence: "Event operators should be studied for composition, commutation, locality, and reachability.",
      implementation: "pair residuals, hotspot/degree maps, event-order checks, calibration experiment protocols and result comparisons",
      limitation: "Reveals disagreement and non-additivity; not a material constitutive law.",
    },
    {
      id: "inverse_design_linearization",
      label: "Linearized inverse design",
      status: "simulator-diagnostic",
      source: "Current numerical architecture",
      pageRefs: ["monograph pp. 23, 52"],
      evidence: "Compiling commands from target shapes is an open research problem.",
      implementation: "response matrix, inverse plan, linear fit",
      limitation: "Plans remain proposals until physical validation or measurement.",
    },
    {
      id: "calibrated_cad_dimensions",
      label: "Calibrated CAD dimensions",
      status: "calibration-gap",
      source: "Physical RAD hardware",
      pageRefs: ["RAD preprint p. 6"],
      evidence: "Prototype scale and tolerance are known, but full part dimensions are not.",
      implementation: "hardwareProfile state, calibrationReadiness, calibrationMeasurementPlan, calibrationExperimentProtocol/results, calibratedRad visual mode, profile-aware OBJ export, future calibrated solver",
      limitation: "Visual/export dimensions can use measurements; contact, friction, and dynamics still need hardware tests.",
    },
  ]);

  function cloneItem(item) {
    return { ...item, pageRefs: [...item.pageRefs] };
  }

  function modelProvenance() {
    return MODEL_PROVENANCE.map(cloneItem);
  }

  function provenanceByStatus(status) {
    return MODEL_PROVENANCE.filter((item) => item.status === status).map(cloneItem);
  }

  function provenanceSummary() {
    const summary = {
      "paper-supported": 0,
      "implementation-assumption": 0,
      "simulator-diagnostic": 0,
      "calibration-gap": 0,
    };
    for (const item of MODEL_PROVENANCE) summary[item.status] = (summary[item.status] || 0) + 1;
    return summary;
  }

  RAD.MODEL_PROVENANCE = MODEL_PROVENANCE;
  RAD.modelProvenance = modelProvenance;
  RAD.provenanceByStatus = provenanceByStatus;
  RAD.provenanceSummary = provenanceSummary;
})();
