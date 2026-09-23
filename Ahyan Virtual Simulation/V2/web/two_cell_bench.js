(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  const TWO_CELL_BENCH_SCHEMA = "rad-sim.two-cell-physical-bench.v1";
  const TWO_CELL_SWEEP_SCHEMA = "rad-sim.two-cell-backlash-sweep.v1";
  const CAD_RAD_CELL_LAYOUT_SCHEMA = "rad-sim.cad-rad-cell-layout.v1";
  const CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA = "rad-sim.cad-rad-cell-reference-profile.v1";
  const TWO_CELL_CONNECTOR_CONTACT_SCHEMA = "rad-sim.two-cell-connector-contact-report.v1";
  const TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA = "rad-sim.two-cell-physical-simulation-suite.v1";
  const TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA = "rad-sim.two-cell-physical-fidelity-matrix.v1";
  const TWO_CELL_CONTACT_PHASE_MAP_SCHEMA = "rad-sim.two-cell-contact-phase-map.v1";
  const TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA = "rad-sim.two-cell-radius-backlash-phase-diagram.v1";
  const TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA = "rad-sim.two-cell-radius-backlash-transition-report.v1";
  const TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA = "rad-sim.two-cell-radius-backlash-transition-comparison.v1";
  const TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA = "rad-sim.two-cell-radius-backlash-transition-rerun.v1";
  const TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA = "rad-sim.two-cell-physical-response-atlas.v1";
  const TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA = "rad-sim.two-cell-cad-contact-decomposition.v1";
  const TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA = "rad-sim.two-cell-exact-contact-handoff-plan.v1";
  const TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA = "rad-sim.two-cell-fidelity-matrix-measurement-template.v1";
  const TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA = "rad-sim.browser-two-cell-fidelity-matrix-measurement-comparison.v1";
  const TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA = "rad-sim.browser-two-cell-fidelity-matrix-parameter-calibration.v1";
  const TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA = "rad-sim.browser-two-cell-external-fidelity-matrix-manifest.v1";
  const CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA = "rad-sim.cad-rad-cell-archive-audit.v1";
  const TWO_CELL_PHYSICAL_FIDELITY_STATUS_SCHEMA = "rad-sim.browser-two-cell-physical-fidelity-status.v1";
  const TWO_CELL_EXTERNAL_FIDELITY_WEB_SUMMARY_SCHEMA = "rad-sim.web-two-cell-external-fidelity-summary.v1";

  const CAD_RAD_CELL_REFERENCE = Object.freeze({
    schema: "rad-sim.cad-rad-cell-reference.v1",
    source: "Autodesk A360 public share https://a360.co/4bMlzip",
    localFiles: {
      fusionArchive: "assets/cad/RADs_unit_cell.f3d",
      viewerSummary: "assets/cad/RADs_unit_cell_reference.json",
      preview: "assets/cad/RADs_unit_cell_preview.png",
    },
    units: "mm",
    boundingBoxMm: {
      widthX: 55.604331129396634,
      widthY: 55.60433117189373,
      heightZ: 19.99999621152464,
    },
    nominalHoleDiameterMm: 3.4,
    viewerAccess: {
      shareUrl: "https://a360.co/4bMlzip",
      resolvedTitle: "RADs unit cell - AUTODESK FUSION",
      downloadFormatsObserved: ["Fusion Archive", "Inventor 2025", "IGES", "SAT", "SMT", "STEP", "DWG", "DXF", "STL", "FBX", "SketchUp", "OBJ"],
      downloadRequiresEmailExport: true,
    },
    visibleFeatures: [
      "central circular crown and screw stack",
      "eight radial arms",
      "outer circular pin-hole pads",
      "single-cell monolithic upper/lower body fragments",
    ],
    fragments: [
      { id: 0, name: "RADs free cell 4mm tall 3.4mm hole body", zMin: 3.9999988079071045, zMax: 8.000000953674316 },
      { id: 1, name: "cell lower body", zMin: 0, zMax: 4 },
      { id: 2, name: "18-8 stainless screw", zMin: -9.999998092651367, zMax: 9.999998092651367 },
    ],
    physicalInterpretation:
      "CAD-derived visual/dimensional reference. The simulator still uses a reduced rigid/contact proxy until the CAD is segmented into moving bodies with joint constraints.",
  });

  const CAD_RAD_CELL_ARCHIVE_AUDIT = Object.freeze({
    schema: CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA,
    source: CAD_RAD_CELL_REFERENCE.source,
    archive: {
      path: CAD_RAD_CELL_REFERENCE.localFiles.fusionArchive,
      exists: true,
      bytes: 106245,
      sha256: "d103feef5c1b437919a4bea87e3bcca97d333948fd7cd9445525a57a01ab1967",
      zipReadable: true,
    },
    summary: {
      archivePresent: true,
      zipReadable: true,
      entryCount: 21,
      fileEntryCount: 13,
      brepEntryCount: 2,
      previewEntryCount: 1,
      manifestEntryCount: 3,
      linkedLabelCount: 2,
      expectedEntriesPresent: true,
      canDeriveVisualReference: true,
      canAttemptExactRigidBodyContact: false,
      physicalAccuracyValidated: false,
      status: "archive-audited-needs-segmentation",
      missingEvidence: [
        "segmentedBodyMeshExport",
        "exactJointAxes",
        "exactPinHoleContactSurfaces",
        "materialContactParameters",
        "benchCoordinateTruth",
      ],
      missingEvidenceCount: 5,
    },
    brepEntries: [
      { name: "FusionAssetName[Active]/Breps.BlobParts/BREP.7d2e0b8b-fbfb-4c1a-a973-6c500ff05860.smb", bytes: 181604, compressedBytes: 21489 },
      { name: "FusionAssetName[Active]/Breps.BlobParts/BREP.f5f9711e-f41d-4b1f-bef9-939b224e2b81.smbh", bytes: 113571, compressedBytes: 14626 },
    ],
    previewEntries: [
      { name: "FusionAssetName[Active]/Previews/small.png", bytes: 12116, compressedBytes: 11366 },
    ],
    linkedLabels: [
      "RADs free cell 4mm tall 3.4mm hole",
      "n/F/_RADs unit cell.04389982-da15-4623-9731-c443d9e44290.f3d",
    ],
    limitations: [
      "Fusion archive entries prove local CAD evidence exists, not exact simulation fidelity.",
      "BREP blobs are not segmented into named moving bodies or pin-hole contact surfaces here.",
      "Exact real-life mechanics still require segmented exports, external contact simulation, and bench comparison.",
    ],
  });

  const CAD_RAD_SITE_ORDER = Object.freeze([
    ["n", 90],
    ["ne", 45],
    ["e", 0],
    ["se", -45],
    ["s", -90],
    ["sw", -135],
    ["w", 180],
    ["nw", 135],
  ]);

  const EXACT_CONTACT_HANDOFF_REQUIRED_CAD_ASSETS = Object.freeze([
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
  ]);
  const EXACT_CONTACT_HANDOFF_REQUIRED_MEASURED_INPUTS = Object.freeze([
    "lock_crown_geometry",
    "actuator_force_displacement",
    "bench_coordinate_truth",
    "pin_hole_friction",
    "normal_contact_stiffness",
  ]);

  const CAD_RAD_X_CONNECTOR_PAIRS = Object.freeze([
    ["upper", "ne", "nw"],
    ["middle", "e", "w"],
    ["lower", "se", "sw"],
  ]);

  const DEFAULT_TWO_CELL_CONTROLS = Object.freeze({
    alphaCommand: -0.35,
    zCommand: 0.35,
    holeSweepMax: 0.5,
    holeSweepSteps: 9,
    leftPositionLocked: true,
    rightLocked: false,
    rightPositionLocked: false,
  });

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function finiteNumber(value, fallback = 0) {
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : fallback;
  }

  function uniqueFloats(values) {
    const out = [];
    for (const value of values || []) {
      const numeric = finiteNumber(value, 0);
      if (!out.some((existing) => Math.abs(existing - numeric) <= 1e-12)) out.push(numeric);
    }
    return out;
  }

  function csvValue(value) {
    if (value === null || value === undefined) return "";
    if (Array.isArray(value)) return `"${value.map((item) => String(item)).join(";").replace(/"/g, '""')}"`;
    if (typeof value === "string") return `"${value.replace(/"/g, '""')}"`;
    return String(value);
  }

  function parseCsv(text) {
    const rows = [];
    let row = [];
    let cell = "";
    let quoted = false;
    const input = String(text || "");
    for (let i = 0; i < input.length; i += 1) {
      const char = input[i];
      if (quoted) {
        if (char === '"' && input[i + 1] === '"') {
          cell += '"';
          i += 1;
        } else if (char === '"') {
          quoted = false;
        } else {
          cell += char;
        }
      } else if (char === '"') {
        quoted = true;
      } else if (char === ",") {
        row.push(cell);
        cell = "";
      } else if (char === "\n") {
        row.push(cell);
        rows.push(row);
        row = [];
        cell = "";
      } else if (char !== "\r") {
        cell += char;
      }
    }
    if (cell !== "" || row.length) {
      row.push(cell);
      rows.push(row);
    }
    const header = rows.shift() || [];
    return rows
      .filter((items) => items.some((item) => String(item || "").trim() !== ""))
      .map((items) => Object.fromEntries(header.map((name, index) => [name, items[index] ?? ""])));
  }

  function numericOrNull(value) {
    if (value === null || value === undefined || value === "") return null;
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : null;
  }

  function nondecreasing(values, tolerance = 1e-9) {
    for (let i = 1; i < values.length; i += 1) {
      if (finiteNumber(values[i], 0) + tolerance < finiteNumber(values[i - 1], 0)) return false;
    }
    return true;
  }

  function nonincreasing(values, tolerance = 1e-9) {
    for (let i = 1; i < values.length; i += 1) {
      if (finiteNumber(values[i], 0) > finiteNumber(values[i - 1], 0) + tolerance) return false;
    }
    return true;
  }

  function twoCellBenchControls(state, overrides = {}) {
    const raw = {
      ...DEFAULT_TWO_CELL_CONTROLS,
      ...(state.experiment?.twoCellBenchControls || {}),
      ...(overrides || {}),
    };
    return {
      alphaCommand: clamp(finiteNumber(raw.alphaCommand, DEFAULT_TWO_CELL_CONTROLS.alphaCommand), -1.2, 1.2),
      zCommand: clamp(finiteNumber(raw.zCommand, DEFAULT_TWO_CELL_CONTROLS.zCommand), -1.4, 1.4),
      holeSweepMax: clamp(finiteNumber(raw.holeSweepMax, DEFAULT_TWO_CELL_CONTROLS.holeSweepMax), 0.02, 1.25),
      holeSweepSteps: Math.max(3, Math.min(41, Math.round(finiteNumber(raw.holeSweepSteps, DEFAULT_TWO_CELL_CONTROLS.holeSweepSteps)))),
      leftPositionLocked: raw.leftPositionLocked !== false,
      rightLocked: raw.rightLocked === true,
      rightPositionLocked: raw.rightPositionLocked === true,
    };
  }

  function stateRadii(state, overrides = {}) {
    const pinRadius = Math.max(0, finiteNumber(overrides.pinRadius ?? state.grid.pinRadius, 0.18));
    const holeRadius = Math.max(pinRadius, finiteNumber(overrides.holeRadius ?? state.grid.holeRadius, 0.225));
    return { pinRadius, holeRadius, clearance: Math.max(0, holeRadius - pinRadius) };
  }

  function cadLayoutDimensions(widthMm, heightMm, pinRadiusMm, holeRadiusMm) {
    const padRadiusMm = Math.max(holeRadiusMm * 1.85, widthMm * 0.055);
    const hubRadiusMm = Math.max(holeRadiusMm * 1.6, widthMm * 0.065);
    return {
      widthMm,
      heightMm,
      bodyThicknessMm: 4,
      lowerBodyThicknessMm: 4,
      stackHeightMm: heightMm,
      nominalHoleRadiusMm: CAD_RAD_CELL_REFERENCE.nominalHoleDiameterMm * 0.5,
      pinRadiusMm,
      holeRadiusMm,
      pinHoleClearanceMm: Math.max(0, holeRadiusMm - pinRadiusMm),
      padRadiusMm,
      hubRadiusMm,
      siteRadiusMm: Math.max(0, 0.5 * widthMm - padRadiusMm),
      connectorPitchMm: 2 * Math.max(0, 0.5 * widthMm - padRadiusMm),
      armStartRadiusMm: hubRadiusMm * 0.78,
      armWidthMm: Math.max(0.85, widthMm * 0.018),
      screwRadiusMm: Math.max(0.65 * pinRadiusMm, 0.8),
    };
  }

  function cadRadCellLayout(state, options = {}) {
    const radii = stateRadii(state, options);
    const widthMm = Math.min(CAD_RAD_CELL_REFERENCE.boundingBoxMm.widthX, CAD_RAD_CELL_REFERENCE.boundingBoxMm.widthY);
    const referenceHoleRadius = Math.max(finiteNumber(state.grid.holeRadius, 0.225), 1e-9);
    const nominalHoleRadiusMm = CAD_RAD_CELL_REFERENCE.nominalHoleDiameterMm * 0.5;
    const pinRadiusMm = (nominalHoleRadiusMm * radii.pinRadius) / referenceHoleRadius;
    const holeRadiusMm = (nominalHoleRadiusMm * radii.holeRadius) / referenceHoleRadius;
    const dimensions = cadLayoutDimensions(widthMm, CAD_RAD_CELL_REFERENCE.boundingBoxMm.heightZ, pinRadiusMm, holeRadiusMm);
    const halfWidth = 0.5 * dimensions.widthMm;
    const padSites = CAD_RAD_SITE_ORDER.map(([name, angleDegrees]) => {
      const radians = (angleDegrees * Math.PI) / 180;
      const xMm = dimensions.siteRadiusMm * Math.cos(radians);
      const yMm = dimensions.siteRadiusMm * Math.sin(radians);
      return {
        name,
        angleDegrees,
        xMm,
        yMm,
        padRadiusMm: dimensions.padRadiusMm,
        holeRadiusMm: dimensions.holeRadiusMm,
        pinRadiusMm: dimensions.pinRadiusMm,
        outerEnvelopeRadiusMm: Math.min(halfWidth, Math.hypot(xMm, yMm) + dimensions.padRadiusMm),
      };
    });
    const maxOuterEnvelopeMm = padSites.reduce(
      (max, site) =>
        Math.max(max, Math.abs(site.xMm) + dimensions.padRadiusMm, Math.abs(site.yMm) + dimensions.padRadiusMm),
      0
    );
    const widthScale = Math.max(dimensions.widthMm, 1e-9);
    return {
      schema: CAD_RAD_CELL_LAYOUT_SCHEMA,
      cadReference: CAD_RAD_CELL_REFERENCE,
      sourceEvidence: [
        "A360 Fusion archive bounding box",
        "linked cell label: RADs free cell 4mm tall 3.4mm hole",
        "eight radial pad sites visible in extracted preview",
      ],
      limitations: [
        "Site layout is dimension-consistent and topology-matched, but not exact BREP vertex extraction.",
        "Fusion BREP streams still need segmentation/export before exact rigid-body contact simulation.",
      ],
      dimensionsMm: dimensions,
      padSites,
      twoCellConnectorPairs: CAD_RAD_X_CONNECTOR_PAIRS.map(([name, leftSite, rightSite]) => ({
        name,
        leftSite,
        rightSite,
        axis: "x",
        leftPinGeom: `left_cell_${leftSite}_pin`,
        rightHoleGeom: `right_cell_${rightSite}_hole_clearance`,
      })),
      envelopeCheck: {
        maxOuterEnvelopeMm,
        halfWidthMm: halfWidth,
        matchesBoundingBox: Math.abs(maxOuterEnvelopeMm - halfWidth) <= 1e-9,
      },
      visualModel: {
        siteRadiusToWidth: dimensions.siteRadiusMm / widthScale,
        padRadiusToWidth: dimensions.padRadiusMm / widthScale,
        hubRadiusToWidth: dimensions.hubRadiusMm / widthScale,
        holeRadiusToWidth: dimensions.holeRadiusMm / widthScale,
        pinRadiusToWidth: dimensions.pinRadiusMm / widthScale,
        armStartRadiusToWidth: dimensions.armStartRadiusMm / widthScale,
        armWidthToWidth: dimensions.armWidthMm / widthScale,
        bodyThicknessToWidth: dimensions.bodyThicknessMm / widthScale,
        stackHeightToWidth: dimensions.stackHeightMm / widthScale,
        siteRadiusToPitch: dimensions.siteRadiusMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        padRadiusToPitch: dimensions.padRadiusMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        hubRadiusToPitch: dimensions.hubRadiusMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        holeRadiusToPitch: dimensions.holeRadiusMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        pinRadiusToPitch: dimensions.pinRadiusMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        holeToPadRadiusRatio: dimensions.holeRadiusMm / Math.max(dimensions.padRadiusMm, 1e-9),
        pinToHoleRadiusRatio: dimensions.pinRadiusMm / Math.max(dimensions.holeRadiusMm, 1e-9),
        screwToHubRadiusRatio: dimensions.screwRadiusMm / Math.max(dimensions.hubRadiusMm, 1e-9),
        hubToPadRadiusRatio: dimensions.hubRadiusMm / Math.max(dimensions.padRadiusMm, 1e-9),
        armStartRadiusToPitch: dimensions.armStartRadiusMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        armWidthToPitch: dimensions.armWidthMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        bodyThicknessToPitch: dimensions.bodyThicknessMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        stackHeightToPitch: dimensions.stackHeightMm / Math.max(dimensions.connectorPitchMm, 1e-9),
        modelUnitsPerMm: Math.max(finiteNumber(state.grid.cellSize, 1), 1e-9) / widthScale,
      },
    };
  }

  function cadRadCellReferenceProfile(state, options = {}) {
    const layout = cadRadCellLayout(state, options);
    const dimensions = layout.dimensionsMm;
    const hardwareProfile = {
      schema: "rad-sim.hardware-profile.v1",
      name: "a360-rads-unit-cell-visual-profile",
      source: CAD_RAD_CELL_REFERENCE.source,
      units: "mm",
      dimensionsMm: {
        sideLengthMm: dimensions.widthMm,
        fabricationHoleToleranceMm: 0.1,
        backlashMm: null,
        pinRadiusMm: dimensions.pinRadiusMm,
        holeRadiusMm: dimensions.holeRadiusMm,
        plateThicknessMm: dimensions.bodyThicknessMm,
        jointStackHeightMm: dimensions.stackHeightMm,
        bossRadiusMm: dimensions.hubRadiusMm,
      },
      derived: {
        pinHoleClearanceMm: dimensions.pinHoleClearanceMm,
      },
      measuredFields: ["pinRadiusMm", "holeRadiusMm", "plateThicknessMm", "jointStackHeightMm", "bossRadiusMm"],
      missingFields: [],
      coverageRatio: 1,
      notes:
        "A360 CAD-derived visual profile from the unit-cell bounding box, 3.4 mm nominal hole label, archive z-fragments, and current normalized pin/hole ratio. Use for CAD-like display and external-engine handoff only until segmented contact and bench data exist.",
    };
    return {
      schema: CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA,
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadLayout: layout,
      hardwareProfile,
      sourceEvidence: [
        "A360 public share metadata identifies the file as RADs unit cell.",
        "Fusion archive/preview provide bounding box and visible radial-pad topology.",
        "Linked model label reports a 3.4 mm hole.",
        "Archive fragments show 4 mm upper body, 4 mm lower body, and 20 mm screw stack envelope.",
      ],
      assumptionLabels: {
        pinRadiusMm: "derived from the current normalized pin/hole ratio and the 3.4 mm nominal hole label",
        holeRadiusMm: "derived from the 3.4 mm nominal hole label",
        plateThicknessMm: "from archive fragment z-ranges, not segmented body mass properties",
        bossRadiusMm: "visual hub radius inferred from bounding-box-scaled layout",
        backlashMm: "not claimed from CAD; keep configurable or measure directly",
      },
      claimBoundary: {
        allowedClaim: "A360-dimensioned visual/profile candidate for real-cell UI and external packet generation",
        blockedClaim: "exact real-cell dynamics, contact, friction, or lock behavior without segmented CAD and bench validation",
      },
    };
  }

  function exportCadRadCellReferenceProfileJson(state, options = {}) {
    return JSON.stringify(cadRadCellReferenceProfile(state, options), null, 2);
  }

  function cadRadCellArchiveAudit() {
    return JSON.parse(JSON.stringify(CAD_RAD_CELL_ARCHIVE_AUDIT));
  }

  function twoCellExternalFidelityWebSummary() {
    const summary = window.RAD_EXTERNAL_FIDELITY_SUMMARY || RAD.EXTERNAL_FIDELITY_SUMMARY || null;
    if (summary && summary.schema === TWO_CELL_EXTERNAL_FIDELITY_WEB_SUMMARY_SCHEMA) {
      return JSON.parse(JSON.stringify(summary));
    }
    return {
      schema: TWO_CELL_EXTERNAL_FIDELITY_WEB_SUMMARY_SCHEMA,
      model: "browser-compact-two-cell-external-fidelity-evidence",
      summary: {
        available: false,
        externalRunComplete: false,
        coversFullFidelityMatrix: false,
        readyForProxyPreview: false,
        physicalAccuracyValidated: false,
        remainingEvidence: ["webExternalFidelitySummaryExport"],
        remainingEvidenceCount: 1,
      },
      claimLabels: {
        physicalAccuracy: "false until segmented CAD contact and independent bench holdout pass",
      },
    };
  }

  function numericValues(rows, key) {
    return rows.map((row) => Number(row?.[key])).filter((value) => Number.isFinite(value));
  }

  function meanField(rows, key) {
    const values = numericValues(rows, key);
    return values.length ? values.reduce((total, value) => total + value, 0) / values.length : null;
  }

  function maxField(rows, key) {
    const values = numericValues(rows, key);
    return values.length ? Math.max(...values) : null;
  }

  function lockModeFlags(lockMode) {
    const mode = String(lockMode || "free");
    return {
      leftPositionLocked: mode !== "left_free",
      rightLocked: mode === "right_state_locked",
      rightPositionLocked: mode === "right_position_locked",
    };
  }

  function twoCellExternalFidelityCaseOptions(summary = twoCellExternalFidelityWebSummary()) {
    const rows = Array.isArray(summary?.representativeRows) ? summary.representativeRows : [];
    const byCase = new Map();
    for (const row of rows) {
      const caseId = String(row?.caseId || "");
      if (!caseId) continue;
      if (!byCase.has(caseId)) byCase.set(caseId, []);
      byCase.get(caseId).push(row);
    }
    return [...byCase.entries()]
      .map(([caseId, caseRows]) => {
        const first = caseRows[0] || {};
        const correctedRightCellZ = meanField(caseRows, "correctedRightCellZ");
        const observedRightCellZ = meanField(caseRows, "observedRightCellZ");
        const correctedVerticalSlipMm = maxField(caseRows, "correctedVerticalSlipMm");
        const observedVerticalSlipMm = maxField(caseRows, "observedVerticalSlipMm");
        const correctedTotalSlipMm = maxField(caseRows, "correctedTotalSlipMm");
        const observedTotalSlipMm = maxField(caseRows, "observedTotalSlipMm");
        return {
          caseId,
          matrixIndex: Number(first.matrixIndex ?? 0),
          actuationCase: first.actuationCase || "",
          lockMode: first.lockMode || "free",
          backlash: finiteNumber(first.backlash, 0),
          pinRadius: finiteNumber(first.pinRadius, 0.18),
          holeRadius: finiteNumber(first.holeRadius, 0.225),
          alphaCommand: finiteNumber(first.alphaCommand, 0),
          zCommand: finiteNumber(first.zCommand, 0),
          gravityForce: finiteNumber(first.gravityForce, 0),
          pinHoleClearanceMm: finiteNumber(first.pinHoleClearanceMm, 0),
          connectorCount: caseRows.length,
          contactModes: [...new Set(caseRows.map((row) => row.observedContactMode || row.predictedContactMode).filter(Boolean))],
          correctedRightCellZ,
          observedRightCellZ,
          correctedVerticalSlipMm,
          observedVerticalSlipMm,
          correctedTotalSlipMm,
          observedTotalSlipMm,
          rows: JSON.parse(JSON.stringify(caseRows)),
          label: `${caseId} ${first.actuationCase || "case"} ${first.lockMode || "free"} b${Number(first.backlash || 0).toFixed(2)} h${Number(first.holeRadius || 0).toFixed(3)}`,
        };
      })
      .sort((a, b) => a.matrixIndex - b.matrixIndex || a.caseId.localeCompare(b.caseId));
  }

  function twoCellExternalFidelityCase(caseId, summary = twoCellExternalFidelityWebSummary()) {
    const id = String(caseId || "");
    return twoCellExternalFidelityCaseOptions(summary).find((item) => item.caseId === id) || null;
  }

  function twoCellExternalFidelityCaseControls(caseId, summary = twoCellExternalFidelityWebSummary(), baseControls = {}) {
    const selected = typeof caseId === "object" && caseId !== null ? caseId : twoCellExternalFidelityCase(caseId, summary);
    if (!selected) return null;
    const flags = lockModeFlags(selected.lockMode);
    const controls = twoCellBenchControls(nullStateForControls(), {
      ...baseControls,
      alphaCommand: selected.alphaCommand,
      zCommand: selected.zCommand,
      holeSweepMax: Math.max(selected.holeRadius, finiteNumber(baseControls.holeSweepMax, selected.holeRadius)),
      ...flags,
    });
    return {
      schema: "rad-sim.browser-two-cell-external-fidelity-case-preview.v1",
      caseId: selected.caseId,
      sourceSchema: summary?.schema || TWO_CELL_EXTERNAL_FIDELITY_WEB_SUMMARY_SCHEMA,
      controls,
      grid: {
        backlash: selected.backlash,
        pinRadius: selected.pinRadius,
        holeRadius: Math.max(selected.pinRadius, selected.holeRadius),
      },
      correctedPreview: {
        rightCellZ: selected.correctedRightCellZ,
        verticalSlipMm: selected.correctedVerticalSlipMm,
        totalSlipMm: selected.correctedTotalSlipMm,
      },
      observedProxy: {
        rightCellZ: selected.observedRightCellZ,
        verticalSlipMm: selected.observedVerticalSlipMm,
        totalSlipMm: selected.observedTotalSlipMm,
      },
      contactModes: selected.contactModes,
      connectorCount: selected.connectorCount,
      connectorRows: JSON.parse(JSON.stringify(selected.rows || [])),
      physicalAccuracyValidated: false,
      remainingEvidence: summary?.summary?.remainingEvidence || ["segmentedCadContactModel", "benchCoordinateHoldout"],
    };
  }

  function exportTwoCellExternalFidelityWebSummaryJson() {
    return JSON.stringify(twoCellExternalFidelityWebSummary(), null, 2);
  }

  function twoCellPhysicalFidelityStatus(state, options = {}, artifacts = {}) {
    const audit = cadRadCellArchiveAudit();
    const suite = artifacts.suite || state.experiment?.twoCellPhysicalSuite || state.experiment?.twoCellBench?.suite || null;
    const matrix = artifacts.matrix || state.experiment?.twoCellPhysicalFidelityMatrix || state.experiment?.twoCellBench?.fidelityMatrix || null;
    const comparison = artifacts.comparison || state.experiment?.twoCellFidelityMatrixMeasurementComparison || null;
    const calibration = artifacts.calibration || state.experiment?.twoCellFidelityMatrixParameterCalibration || null;
    const missing = new Set([
      ...(audit.summary?.missingEvidence || []),
      ...(suite?.summary?.missingEvidence || []),
      ...(matrix?.summary?.missingEvidence || []),
    ]);
    if (!comparison?.summary?.comparisonReady) missing.add("filledFidelityMatrixMeasurements");
    if (comparison?.summary?.comparisonReady && !calibration?.summary?.readyForReducedProxyCalibration) {
      missing.add("fidelityMatrixParameterCalibration");
    }
    const reducedProxyReady = Boolean(
      suite?.summary?.internalSuiteReady || matrix?.summary?.internalMatrixReady || state.experiment?.twoCellBench?.bench
    );
    const matrixMeasurementRowsExpected = matrix?.summary?.connectorRowCount || 756;
    return {
      schema: TWO_CELL_PHYSICAL_FIDELITY_STATUS_SCHEMA,
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadArchiveAudit: audit,
      suiteSummary: suite?.summary || null,
      matrixSummary: matrix?.summary || null,
      measurementComparisonSummary: comparison?.summary || null,
      parameterCalibrationSummary: calibration?.summary || null,
      proposedReducedProxyUpdates: calibration?.proposedReducedProxyUpdates || null,
      summary: {
        status: audit.summary.canAttemptExactRigidBodyContact
          ? "exact-cad-ready-needs-bench-validation"
          : "reduced-proxy-active-exact-cad-missing-segmentation",
        reducedProxyReady,
        archiveBrepsPresent: audit.summary.brepEntryCount > 0,
        exactGeometryReady: audit.summary.canAttemptExactRigidBodyContact,
        physicalAccuracyValidated: false,
        matrixRunReady: Boolean(matrix?.summary?.internalMatrixReady),
        suiteRunReady: Boolean(suite?.summary?.internalSuiteReady),
        measurementComparisonReady: Boolean(comparison?.summary?.comparisonReady),
        parameterCalibrationReady: Boolean(calibration?.summary?.readyForReducedProxyCalibration),
        reducedProxyMeasurementValidated: Boolean(comparison?.summary?.passesTolerance),
        matrixMeasurementRowsExpected,
        missingEvidence: [...missing].sort(),
        missingEvidenceCount: missing.size,
      },
      nextActions: [
        "Export segmented moving body meshes and exact pin-hole surfaces from Fusion.",
        "Fill the 756-row fidelity matrix measurement template from bench or external contact runs.",
        "Run MuJoCo/Gazebo/Isaac contact with the segmented bodies and compare connector residuals.",
        "Only mark physical accuracy validated after measured coordinates and contact modes match tolerance.",
      ],
    };
  }

  function exportCadRadCellArchiveAuditJson() {
    return JSON.stringify(cadRadCellArchiveAudit(), null, 2);
  }

  function exportTwoCellPhysicalFidelityStatusJson(state, options = {}, artifacts = {}) {
    return JSON.stringify(twoCellPhysicalFidelityStatus(state, options, artifacts), null, 2);
  }

  function siteByName(layout, name) {
    return (layout.padSites || []).find((site) => site.name === name);
  }

  function connectorSitePositionMm(cell, site, layout, cellSize, initialAlpha, thetaYawGain = 0.18) {
    const pitchMm = layout.dimensionsMm.connectorPitchMm;
    const modelToMm = pitchMm / Math.max(cellSize, 1e-9);
    const centerX = finiteNumber(cell.center?.x, 0) * modelToMm;
    const centerY = finiteNumber(cell.center?.y, 0) * modelToMm;
    const centerZ = finiteNumber(cell.center?.z, 0) * layout.dimensionsMm.heightMm;
    const alphaScale = Math.sqrt(Math.max(finiteNumber(cell.alpha, 1), 1e-9) / Math.max(initialAlpha, 1e-9));
    const thetaDelta = finiteNumber(cell.theta, alphaToTheta(initialAlpha)) - alphaToTheta(initialAlpha);
    const yaw = (thetaDelta * thetaYawGain * Math.PI) / 180;
    const localX = finiteNumber(site?.xMm, 0) * alphaScale;
    const localY = finiteNumber(site?.yMm, 0) * alphaScale;
    return {
      x: centerX + Math.cos(yaw) * localX - Math.sin(yaw) * localY,
      y: centerY + Math.sin(yaw) * localX + Math.cos(yaw) * localY,
      z: centerZ,
    };
  }

  function twoCellConnectorContactReport(state, options = {}) {
    const bench = simulateTwoCellBench(state, options);
    const layout = bench.cadLayout || cadRadCellLayout(state, options);
    const clearance = finiteNumber(layout.dimensionsMm.pinHoleClearanceMm, 0);
    const initialAlpha = finiteNumber(state.grid.initialAlpha, 1);
    const cellSize = Math.max(1e-9, finiteNumber(state.grid.cellSize, 1));
    const [left, right] = bench.cells;
    const connectors = (layout.twoCellConnectorPairs || []).map((pair) => {
      const leftSite = siteByName(layout, pair.leftSite);
      const rightSite = siteByName(layout, pair.rightSite);
      const leftPositionMm = connectorSitePositionMm(left, leftSite, layout, cellSize, initialAlpha);
      const rightPositionMm = connectorSitePositionMm(right, rightSite, layout, cellSize, initialAlpha);
      const slip = {
        x: rightPositionMm.x - leftPositionMm.x,
        y: rightPositionMm.y - leftPositionMm.y,
        z: rightPositionMm.z - leftPositionMm.z,
      };
      const lateralSlipMm = Math.hypot(slip.x, slip.y);
      const verticalSlipMm = Math.abs(slip.z);
      const totalSlipMm = Math.hypot(lateralSlipMm, slip.z);
      const lateralExcessMm = Math.max(0, lateralSlipMm - clearance);
      const verticalExcessMm = Math.max(0, verticalSlipMm - clearance);
      const contactPenalty = 0.5 * 12 * lateralExcessMm * lateralExcessMm + 0.5 * 8 * verticalExcessMm * verticalExcessMm;
      const contactMode =
        lateralExcessMm > 1e-9 && verticalExcessMm > 1e-9
          ? "lateral+vertical-contact"
          : lateralExcessMm > 1e-9
            ? "lateral-contact"
            : verticalExcessMm > 1e-9
              ? "vertical-contact"
              : "inside-clearance";
      return {
        connector: pair.name,
        leftSite: pair.leftSite,
        rightSite: pair.rightSite,
        leftPositionMm,
        rightPositionMm,
        slipMm: slip,
        lateralSlipMm,
        verticalSlipMm,
        totalSlipMm,
        lateralClearanceMm: clearance,
        verticalFreePlayMm: clearance,
        lateralExcessMm,
        verticalExcessMm,
        contactPenalty,
        contactMode,
      };
    });
    return {
      schema: TWO_CELL_CONNECTOR_CONTACT_SCHEMA,
      model: "cad-layout-two-cell-connector-clearance-contact",
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadLayout: layout,
      controls: bench.controls,
      dimensions: {
        ...bench.dimensions,
        connectorPitchMm: layout.dimensionsMm.connectorPitchMm,
        pinRadiusMm: layout.dimensionsMm.pinRadiusMm,
        holeRadiusMm: layout.dimensionsMm.holeRadiusMm,
        pinHoleClearanceMm: clearance,
      },
      bench,
      connectors,
      summary: {
        connectorContactReady: connectors.length === 3,
        connectorCount: connectors.length,
        activeContactCount: connectors.filter((row) => row.contactMode !== "inside-clearance").length,
        maxLateralSlipMm: Math.max(0, ...connectors.map((row) => row.lateralSlipMm)),
        maxVerticalSlipMm: Math.max(0, ...connectors.map((row) => row.verticalSlipMm)),
        maxLateralExcessMm: Math.max(0, ...connectors.map((row) => row.lateralExcessMm)),
        maxVerticalExcessMm: Math.max(0, ...connectors.map((row) => row.verticalExcessMm)),
        totalContactPenalty: connectors.reduce((sum, row) => sum + row.contactPenalty, 0),
        requiresExternalEngineRun: true,
        requiresMeasuredData: true,
      },
    };
  }

  function benchDeadZones(state, radii) {
    const backlash = Math.max(0, finiteNumber(state.grid.backlash, 0.1));
    const clearance = Math.max(0, radii.clearance);
    const verticalGap = typeof RAD.verticalDeadZone === "function" ? RAD.verticalDeadZone({ ...state, grid: { ...state.grid, pinRadius: radii.pinRadius, holeRadius: radii.holeRadius } }) : clearance;
    return {
      backlash,
      clearance,
      alphaDeadZone: backlash + 0.5 * clearance,
      verticalDeadZone: Math.max(0, verticalGap),
      effectiveBacklash: backlash + clearance,
    };
  }

  function alphaToTheta(alpha) {
    return typeof RAD.alphaToTheta === "function" ? RAD.alphaToTheta(alpha) : 70 * alpha - 60;
  }

  function deadZone(value, gap) {
    return typeof RAD.backlashActivation === "function"
      ? RAD.backlashActivation(value, gap)
      : Math.max(0, value - gap) + Math.min(value + gap, 0);
  }

  function twoCellRestPitch(state, cellSize, initialAlpha, leftAlpha, rightAlpha) {
    if (typeof RAD.sharedEdgePitch === "function") return RAD.sharedEdgePitch(state, leftAlpha, rightAlpha);
    const denominator = Math.max(1e-9, initialAlpha);
    return 0.5 * cellSize * (leftAlpha / denominator + rightAlpha / denominator);
  }

  function simulateTwoCellBench(state, options = {}) {
    const controls = twoCellBenchControls(state, options);
    const radii = stateRadii(state, options);
    const zones = benchDeadZones(state, radii);
    const cellSize = Math.max(1e-9, finiteNumber(state.grid.cellSize, 1));
    const initialAlpha = finiteNumber(state.grid.initialAlpha, 1);
    const couplingGain = clamp(finiteNumber(state.grid.couplingGain, 0.55), 0, 1);
    const zCouplingGain = clamp(finiteNumber(state.grid.zCouplingGain, 0.32), 0, 1);
    const axialStiffness = Math.max(1e-9, finiteNumber(options.axialStiffness, 8));
    const verticalStiffness = Math.max(1e-9, finiteNumber(options.verticalStiffness, 5));
    const hingeStiffness = Math.max(1e-9, finiteNumber(options.hingeStiffness, 2.5));
    const lockStiffness = Math.max(1e-9, finiteNumber(options.lockStiffness, 18));
    const contactStiffness = Math.max(1e-9, finiteNumber(options.contactStiffness, 24));

    const alphaFree = deadZone(controls.alphaCommand, zones.alphaDeadZone);
    const zFree = deadZone(controls.zCommand, zones.verticalDeadZone);
    const alphaTransmission = couplingGain * Math.exp(-2.25 * zones.alphaDeadZone);
    const zTransmission = zCouplingGain * Math.exp(-3.1 * zones.verticalDeadZone);

    const leftFixed = controls.leftPositionLocked;
    const rightAlphaBlocked = controls.rightLocked || controls.rightPositionLocked;
    const rightZBlocked = controls.rightPositionLocked;
    const rightAlphaDelta = rightAlphaBlocked ? 0 : alphaFree;
    const rightZ = rightZBlocked ? 0 : zFree;
    const leftAlphaDelta = leftFixed ? 0 : rightAlphaDelta * alphaTransmission;
    const leftZ = leftFixed ? 0 : rightZ * zTransmission;
    const leftAlpha = clamp(initialAlpha + leftAlphaDelta, finiteNumber(state.grid.alphaMin, 0.25), finiteNumber(state.grid.alphaMax, 1.75));
    const rightAlpha = clamp(initialAlpha + rightAlphaDelta, finiteNumber(state.grid.alphaMin, 0.25), finiteNumber(state.grid.alphaMax, 1.75));
    const pitch = cellSize;
    const restPitch = twoCellRestPitch(state, cellSize, initialAlpha, leftAlpha, rightAlpha);
    const freeLeftX = leftFixed ? 0 : 0.5 * (pitch - restPitch);
    const freeRightX = leftFixed ? restPitch : 0.5 * (pitch + restPitch);
    const leftX = leftFixed ? 0 : controls.rightPositionLocked ? pitch - restPitch : freeLeftX;
    const rightX = controls.rightPositionLocked ? pitch : freeRightX;
    const distance = Math.hypot(rightX - leftX, rightZ - leftZ);
    const axialStrain = distance / Math.max(1e-9, restPitch) - 1;
    const leftTheta = alphaToTheta(leftAlpha);
    const rightTheta = alphaToTheta(rightAlpha);
    const thetaDeltaRadians = ((rightTheta - leftTheta) * Math.PI) / 180;
    const verticalShear = rightZ - leftZ;
    const contactPenetration = Math.max(0, Math.abs(controls.zCommand) - zones.verticalDeadZone);
    const alphaUnrealized = Math.abs(alphaFree - rightAlphaDelta);
    const zUnrealized = Math.abs(zFree - rightZ);
    const positionUnrealized = controls.rightPositionLocked ? Math.hypot(freeRightX - pitch, zFree) : 0;
    const lockedUnrealized = controls.rightLocked ? Math.abs(alphaFree) : 0;
    const axialEnergy = 0.5 * axialStiffness * axialStrain * axialStrain;
    const hingeEnergy = 0.5 * hingeStiffness * thetaDeltaRadians * thetaDeltaRadians;
    const verticalEnergy = 0.5 * verticalStiffness * verticalShear * verticalShear;
    const contactEnergy = 0.5 * contactStiffness * contactPenetration * contactPenetration;
    const lockPenaltyEnergy = 0.5 * lockStiffness * (alphaUnrealized * alphaUnrealized + zUnrealized * zUnrealized + positionUnrealized * positionUnrealized + lockedUnrealized * lockedUnrealized);
    const totalEnergy = axialEnergy + hingeEnergy + verticalEnergy + contactEnergy + lockPenaltyEnergy;
    const contactMode = controls.rightPositionLocked
      ? "blocked-by-position-lock"
      : controls.rightLocked
        ? "state-lock-alpha-only"
      : Math.abs(controls.alphaCommand) <= zones.alphaDeadZone && Math.abs(controls.zCommand) <= zones.verticalDeadZone
        ? "inside-backlash"
        : contactPenetration > 0
          ? "clearance-taken-up"
          : "alpha-contact-only";

    return {
      schema: TWO_CELL_BENCH_SCHEMA,
      model: "cad-derived-two-cell-rigid-contact-proxy",
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadLayout: cadRadCellLayout(state, options),
      units: "normalized-cell-units",
      controls,
      dimensions: {
        cellSize,
        pinRadius: radii.pinRadius,
        holeRadius: radii.holeRadius,
        pinHoleClearance: radii.clearance,
        pinRadiusMm: typeof RAD.modelLengthToMm === "function" ? RAD.modelLengthToMm(state, radii.pinRadius) : null,
        holeRadiusMm: typeof RAD.modelLengthToMm === "function" ? RAD.modelLengthToMm(state, radii.holeRadius) : null,
        clearanceMm: typeof RAD.modelLengthToMm === "function" ? RAD.modelLengthToMm(state, radii.clearance) : null,
      },
      deadZones: zones,
      cells: [
        { id: "left", positionLocked: leftFixed, locked: false, alpha: leftAlpha, theta: leftTheta, center: { x: leftX, y: 0, z: leftZ }, residualAlpha: leftAlphaDelta, residualZ: leftZ },
        { id: "right", positionLocked: controls.rightPositionLocked, locked: controls.rightLocked, alpha: rightAlpha, theta: rightTheta, center: { x: rightX, y: 0, z: rightZ }, residualAlpha: rightAlphaDelta, residualZ: rightZ },
      ],
      connector: {
        pitch,
        restPitch,
        distance,
        axialStrain,
        verticalShear,
        contactMode,
        contactPenetration,
        alphaTransmission,
        zTransmission,
      },
      energy: {
        axialEnergy,
        hingeEnergy,
        verticalEnergy,
        contactEnergy,
        lockPenaltyEnergy,
        totalEnergy,
      },
      interpretation: {
        calibratedAccuracy: "not-yet-calibrated",
        limitation: "This is a reduced two-body proxy driven by CAD dimensions and measured-state controls, not a solved rigid-body collision model.",
        nextEvidenceNeeded: ["segmented STEP assembly", "pin-hole contact friction", "joint stiffness", "two-cell measured actuation traces"],
      },
    };
  }

  function sweepTwoCellBacklash(state, options = {}) {
    const controls = twoCellBenchControls(state, options);
    const radii = stateRadii(state, options);
    const maxHole = Math.max(radii.pinRadius, finiteNumber(options.holeSweepMax ?? controls.holeSweepMax, controls.holeSweepMax));
    const steps = Math.max(3, Math.min(41, Math.round(finiteNumber(options.holeSweepSteps ?? controls.holeSweepSteps, controls.holeSweepSteps))));
    const rows = [];
    for (let i = 0; i < steps; i += 1) {
      const t = steps === 1 ? 0 : i / (steps - 1);
      const holeRadius = radii.pinRadius + (maxHole - radii.pinRadius) * t;
      const result = simulateTwoCellBench(state, { ...controls, ...options, holeRadius });
      rows.push({
        index: i,
        holeRadius,
        clearance: result.dimensions.pinHoleClearance,
        effectiveBacklash: result.deadZones.effectiveBacklash,
        alphaDeadZone: result.deadZones.alphaDeadZone,
        verticalDeadZone: result.deadZones.verticalDeadZone,
        leftResidualZ: result.cells[0].residualZ,
        rightZ: result.cells[1].center.z,
        axialStrain: result.connector.axialStrain,
        contactMode: result.connector.contactMode,
        totalEnergy: result.energy.totalEnergy,
      });
    }
    const first = rows[0] || {};
    const last = rows[rows.length - 1] || {};
    return {
      schema: TWO_CELL_SWEEP_SCHEMA,
      model: "cad-derived-two-cell-clearance-sweep",
      controls,
      pinRadius: radii.pinRadius,
      holeRadiusStart: rows[0]?.holeRadius ?? radii.holeRadius,
      holeRadiusEnd: rows.at(-1)?.holeRadius ?? radii.holeRadius,
      rows,
      trend: {
        neighborZStart: first.leftResidualZ ?? 0,
        neighborZEnd: last.leftResidualZ ?? 0,
        rightZStart: first.rightZ ?? 0,
        rightZEnd: last.rightZ ?? 0,
        energyStart: first.totalEnergy ?? 0,
        energyEnd: last.totalEnergy ?? 0,
        clearanceIncreasesBacklash: (last.effectiveBacklash ?? 0) >= (first.effectiveBacklash ?? 0),
        neighborResponseDropsWithClearance: Math.abs(last.leftResidualZ ?? 0) <= Math.abs(first.leftResidualZ ?? 0) + 1e-9,
      },
    };
  }

  function controlsForLockMode(controls, mode) {
    const base = { ...controls };
    if (mode === "right_state_locked") {
      return twoCellBenchControls(nullStateForControls(), { ...base, leftPositionLocked: true, rightLocked: true, rightPositionLocked: false });
    }
    if (mode === "right_position_locked") {
      return twoCellBenchControls(nullStateForControls(), { ...base, leftPositionLocked: true, rightLocked: false, rightPositionLocked: true });
    }
    if (mode === "left_free") {
      return twoCellBenchControls(nullStateForControls(), { ...base, leftPositionLocked: false, rightLocked: false, rightPositionLocked: false });
    }
    return twoCellBenchControls(nullStateForControls(), { ...base, leftPositionLocked: true, rightLocked: false, rightPositionLocked: false });
  }

  function nullStateForControls() {
    return { grid: {}, experiment: { twoCellBenchControls: {} } };
  }

  function twoCellPhysicalSuiteCaseSpecs(state, options = {}) {
    const controls = twoCellBenchControls(state, options);
    const radii = stateRadii(state, options);
    const nominalHole = radii.holeRadius;
    const tightHole = radii.pinRadius;
    const looseHole = Math.max(controls.holeSweepMax, nominalHole);
    const backlash = finiteNumber(state.grid.backlash, 0.1);
    const highBacklash = Math.max(backlash * 3, 0.25);
    const alphaDrive = Math.abs(controls.alphaCommand) > 0.05 ? controls.alphaCommand : -0.35;
    const zDrive = Math.abs(controls.zCommand) > 0.05 ? controls.zCommand : 0.35;
    return [
      ["free_contract_lift_nominal", "baseline contraction plus vertical lift", "free", alphaDrive, zDrive, nominalHole, backlash, 0.025],
      ["tight_clearance_contract_lift", "minimum pin-hole clearance response", "free", alphaDrive, zDrive, tightHole, backlash, 0.025],
      ["loose_clearance_contract_lift", "large pin-hole clearance response", "free", alphaDrive, zDrive, looseHole, backlash, 0.025],
      ["free_expand_pushdown_nominal", "opposite actuation polarity", "free", Math.abs(alphaDrive), -Math.abs(zDrive), nominalHole, backlash, 0.025],
      ["alpha_only_contract", "horizontal/dilation actuation without z command", "free", alphaDrive, 0, nominalHole, backlash, 0.025],
      ["z_only_lift", "vertical actuation without alpha command", "free", 0, Math.abs(zDrive), nominalHole, backlash, 0.025],
      ["right_state_locked_lift", "state lock freezes alpha while z remains free", "right_state_locked", alphaDrive, Math.abs(zDrive), nominalHole, backlash, 0.025],
      ["right_position_locked_lift", "position lock fixes fixture coordinates", "right_position_locked", alphaDrive, Math.abs(zDrive), nominalHole, backlash, 0.025],
      ["left_free_neighbor_residual_lift", "residual motion transmitted into free neighbor", "left_free", alphaDrive, Math.abs(zDrive), nominalHole, backlash, 0.025],
      ["zero_backlash_left_free_lift", "neighbor response with backlash removed", "left_free", alphaDrive, Math.abs(zDrive), nominalHole, 0, 0.025],
      ["high_backlash_left_free_lift", "neighbor response with high backlash", "left_free", alphaDrive, Math.abs(zDrive), nominalHole, highBacklash, 0.025],
      ["gravity_sag_loose_clearance", "zero-command loose-clearance gravity proxy", "free", 0, 0, looseHole, backlash, 0.5],
    ].map(([caseId, purpose, lockMode, alphaCommand, zCommand, holeRadius, caseBacklash, gravityForce], index) => ({
      index,
      caseId,
      purpose,
      lockMode,
      alphaCommand,
      zCommand,
      holeRadius,
      backlash: caseBacklash,
      gravityForce,
    }));
  }

  function stateForTwoCellCase(state, spec) {
    return {
      ...state,
      grid: {
        ...state.grid,
        rows: 1,
        cols: 2,
        backlash: Math.max(0, finiteNumber(spec.backlash, state.grid.backlash)),
        holeRadius: Math.max(finiteNumber(state.grid.pinRadius, 0.18), finiteNumber(spec.holeRadius, state.grid.holeRadius)),
      },
    };
  }

  function twoCellPhysicalSuiteCaseOptions(state, caseId, options = {}) {
    const spec = twoCellPhysicalSuiteCaseSpecs(state, options).find((item) => item.caseId === caseId) || twoCellPhysicalSuiteCaseSpecs(state, options)[0];
    const controls = controlsForLockMode(
      twoCellBenchControls(state, {
        ...options,
        alphaCommand: spec.alphaCommand,
        zCommand: spec.zCommand,
      }),
      spec.lockMode
    );
    return {
      caseId: spec.caseId,
      purpose: spec.purpose,
      backlash: spec.backlash,
      holeRadius: spec.holeRadius,
      controls,
    };
  }

  function twoCellPhysicalSimulationSuite(state, options = {}) {
    const initialAlpha = finiteNumber(state.grid.initialAlpha, 1);
    const specs = twoCellPhysicalSuiteCaseSpecs(state, options);
    const rows = specs.map((spec) => {
      const caseState = stateForTwoCellCase(state, spec);
      const caseControls = controlsForLockMode(
        twoCellBenchControls(state, {
          ...options,
          alphaCommand: spec.alphaCommand,
          zCommand: spec.zCommand,
        }),
        spec.lockMode
      );
      const bench = simulateTwoCellBench(caseState, caseControls);
      const contact = twoCellConnectorContactReport(caseState, caseControls);
      const left = bench.cells[0];
      const right = bench.cells[1];
      const lockPass =
        spec.lockMode === "right_position_locked"
          ? Math.abs((right.center?.x || 0) - finiteNumber(caseState.grid.cellSize, 1)) <= 1e-9 && Math.abs(right.center?.z || 0) <= 1e-9
          : spec.lockMode === "right_state_locked"
            ? Math.abs((right.alpha || initialAlpha) - initialAlpha) <= 1e-9
            : true;
      return {
        index: spec.index,
        caseId: spec.caseId,
        purpose: spec.purpose,
        lockMode: spec.lockMode,
        backlash: caseState.grid.backlash,
        alphaCommand: caseControls.alphaCommand,
        zCommand: caseControls.zCommand,
        gravityForce: spec.gravityForce,
        pinRadius: bench.dimensions.pinRadius,
        holeRadius: bench.dimensions.holeRadius,
        clearance: bench.dimensions.pinHoleClearance,
        effectiveBacklash: bench.deadZones.effectiveBacklash,
        pinHoleClearanceMm: contact.dimensions.pinHoleClearanceMm,
        solverSuccess: true,
        lockPass,
        positionLocksHold: spec.lockMode !== "right_position_locked" || lockPass,
        stateLocksHold: spec.lockMode !== "right_state_locked" || lockPass,
        rightLockedVerticalMotionAllowed: spec.lockMode === "right_state_locked" && Math.abs(right.center?.z || 0) > 1e-9,
        leftX: left.center.x,
        leftZ: left.center.z,
        rightX: right.center.x,
        rightZ: right.center.z,
        leftAlpha: left.alpha,
        rightAlpha: right.alpha,
        leftAlphaDeltaFromInitial: left.alpha - initialAlpha,
        rightAlphaDeltaFromInitial: right.alpha - initialAlpha,
        leftTheta: left.theta,
        rightTheta: right.theta,
        benchLeftResidualZ: left.residualZ,
        verticalShear: bench.connector.verticalShear,
        axialStrain: bench.connector.axialStrain,
        contactMode: bench.connector.contactMode,
        connectorActiveContactCount: contact.summary.activeContactCount,
        connectorMaxVerticalSlipMm: contact.summary.maxVerticalSlipMm,
        connectorMaxVerticalExcessMm: contact.summary.maxVerticalExcessMm,
        connectorTotalContactPenalty: contact.summary.totalContactPenalty,
        totalEnergy: bench.energy.totalEnergy,
        externalMeasurementRequired: true,
      };
    });
    const byId = Object.fromEntries(rows.map((row) => [row.caseId, row]));
    const tight = byId.tight_clearance_contract_lift || {};
    const loose = byId.loose_clearance_contract_lift || {};
    const zero = byId.zero_backlash_left_free_lift || {};
    const high = byId.high_backlash_left_free_lift || {};
    const clearanceTrend = {
      tightClearanceMm: finiteNumber(tight.pinHoleClearanceMm, 0),
      looseClearanceMm: finiteNumber(loose.pinHoleClearanceMm, 0),
      looseHasLargerClearance: finiteNumber(loose.pinHoleClearanceMm, 0) >= finiteNumber(tight.pinHoleClearanceMm, 0),
      tightContactPenalty: finiteNumber(tight.connectorTotalContactPenalty, 0),
      looseContactPenalty: finiteNumber(loose.connectorTotalContactPenalty, 0),
      loosePenaltyNotHigher: finiteNumber(loose.connectorTotalContactPenalty, 0) <= finiteNumber(tight.connectorTotalContactPenalty, 0) + 1e-9,
    };
    const backlashTrend = {
      zeroBacklash: finiteNumber(zero.backlash, 0),
      highBacklash: finiteNumber(high.backlash, 0),
      zeroNeighborResidualZ: finiteNumber(zero.benchLeftResidualZ, 0),
      highNeighborResidualZ: finiteNumber(high.benchLeftResidualZ, 0),
      zeroLeftAlphaDelta: finiteNumber(zero.leftAlphaDeltaFromInitial, 0),
      highLeftAlphaDelta: finiteNumber(high.leftAlphaDeltaFromInitial, 0),
      highBacklashReducesAlphaResponse: Math.abs(finiteNumber(high.leftAlphaDeltaFromInitial, 0)) <= Math.abs(finiteNumber(zero.leftAlphaDeltaFromInitial, 0)) + 1e-9,
    };
    const lockRows = rows.filter((row) => row.lockMode === "right_state_locked" || row.lockMode === "right_position_locked");
    const internalSuiteReady =
      rows.every((row) => row.solverSuccess) &&
      lockRows.every((row) => row.lockPass) &&
      clearanceTrend.looseHasLargerClearance &&
      backlashTrend.highBacklashReducesAlphaResponse;
    return {
      schema: TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA,
      model: "browser-named-two-cell-reduced-physics-test-suite",
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadLayout: cadRadCellLayout(state, options),
      rows,
      summary: {
        status: internalSuiteReady ? "reduced-two-cell-suite-ready-needs-external-data" : "reduced-two-cell-suite-internal-review-needed",
        internalSuiteReady,
        physicalAccuracyValidated: false,
        caseCount: rows.length,
        solverSuccessCount: rows.filter((row) => row.solverSuccess).length,
        solverFailureCount: rows.filter((row) => !row.solverSuccess).length,
        lockCaseCount: lockRows.length,
        lockPassCount: lockRows.filter((row) => row.lockPass).length,
        maxVerticalSlipMm: Math.max(0, ...rows.map((row) => row.connectorMaxVerticalSlipMm)),
        maxVerticalExcessMm: Math.max(0, ...rows.map((row) => row.connectorMaxVerticalExcessMm)),
        maxContactPenalty: Math.max(0, ...rows.map((row) => row.connectorTotalContactPenalty)),
        clearanceTrend,
        backlashTrend,
        missingEvidence: [
          "segmentedStepOrMeshExport",
          "externalMuJoCoRun",
          "benchCoordinates",
          "measuredPinHoleClearance",
          "measuredConnectorSlip",
          "measuredFriction",
          "measuredContactStiffness",
          "measuredJointAxes",
        ],
      },
      claimLabels: {
        simulation: "browser reduced two-cell suite",
        physicalAccuracy: "not validated until external engine and bench measurements are attached",
      },
    };
  }

  function exportTwoCellBenchJson(state, options = {}) {
    return JSON.stringify(
      {
        bench: simulateTwoCellBench(state, options),
        sweep: sweepTwoCellBacklash(state, options),
        physicalSuite: twoCellPhysicalSimulationSuite(state, options),
      },
      null,
      2
    );
  }

  function exportTwoCellBacklashSweepCsv(state, options = {}) {
    const sweep = sweepTwoCellBacklash(state, options);
    const columns = [
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
    ];
    return [
      columns.join(","),
      ...sweep.rows.map((row) =>
        columns
          .map((column) => {
            const value = row[column];
            return typeof value === "string" ? `"${value.replace(/"/g, '""')}"` : String(value);
          })
          .join(",")
      ),
    ].join("\n");
  }

  function exportTwoCellPhysicalSuiteJson(state, options = {}) {
    return JSON.stringify(twoCellPhysicalSimulationSuite(state, options), null, 2);
  }

  function exportTwoCellPhysicalSuiteCsv(state, options = {}) {
    const suite = twoCellPhysicalSimulationSuite(state, options);
    const columns = [
      "index",
      "caseId",
      "lockMode",
      "backlash",
      "alphaCommand",
      "zCommand",
      "pinRadius",
      "holeRadius",
      "clearance",
      "pinHoleClearanceMm",
      "lockPass",
      "rightLockedVerticalMotionAllowed",
      "leftZ",
      "rightZ",
      "leftAlphaDeltaFromInitial",
      "rightAlphaDeltaFromInitial",
      "contactMode",
      "connectorMaxVerticalSlipMm",
      "connectorMaxVerticalExcessMm",
      "connectorTotalContactPenalty",
      "externalMeasurementRequired",
    ];
    return [
      columns.join(","),
      ...suite.rows.map((row) =>
        columns
          .map((column) => {
            const value = row[column];
            return typeof value === "string" ? `"${value.replace(/"/g, '""')}"` : String(value);
          })
          .join(",")
      ),
    ].join("\n");
  }

  function defaultFidelityActuationCases(controls) {
    const alphaDrive = Math.abs(controls.alphaCommand) > 0.05 ? Math.abs(controls.alphaCommand) : 0.35;
    const zDrive = Math.abs(controls.zCommand) > 0.05 ? Math.abs(controls.zCommand) : 0.35;
    return [
      { actuationCase: "contract_lift", alpha: -alphaDrive, z: zDrive, gravity: 0.025 },
      { actuationCase: "expand_pushdown", alpha: alphaDrive, z: -zDrive, gravity: 0.025 },
      { actuationCase: "alpha_contract", alpha: -alphaDrive, z: 0, gravity: 0.025 },
      { actuationCase: "alpha_expand", alpha: alphaDrive, z: 0, gravity: 0.025 },
      { actuationCase: "z_lift", alpha: 0, z: zDrive, gravity: 0.025 },
      { actuationCase: "z_pushdown", alpha: 0, z: -zDrive, gravity: 0.025 },
      { actuationCase: "gravity_sag", alpha: 0, z: 0, gravity: 0.5 },
    ];
  }

  function matrixLockPass(lockMode, bench, tolerance = 1e-9) {
    const initialAlpha = finiteNumber(bench?.dimensions?.initialAlpha, finiteNumber(bench?.cadLayout?.grid?.initialAlpha, 1));
    const right = bench?.cells?.[1] || {};
    if (lockMode === "right_position_locked") {
      const cellSize = finiteNumber(bench?.dimensions?.cellSize, 1);
      return Math.abs(finiteNumber(right.center?.x, 0) - cellSize) <= tolerance && Math.abs(finiteNumber(right.center?.z, 0)) <= tolerance;
    }
    if (lockMode === "right_state_locked") {
      return Math.abs(finiteNumber(right.alpha, initialAlpha) - initialAlpha) <= tolerance;
    }
    return true;
  }

  function matrixRowLookup(rows, query, tolerance = 1e-9) {
    return rows.find(
      (row) =>
        row.actuationCase === query.actuationCase &&
        row.lockMode === query.lockMode &&
        Math.abs(finiteNumber(row.backlash, 0) - finiteNumber(query.backlash, 0)) <= tolerance &&
        Math.abs(finiteNumber(row.holeRadius, 0) - finiteNumber(query.holeRadius, 0)) <= tolerance
    );
  }

  function freePlayUtilization(slipMm, clearanceMm) {
    const clearance = Math.max(0, finiteNumber(clearanceMm, 0));
    if (clearance <= 1e-12) return null;
    return Math.abs(finiteNumber(slipMm, 0)) / clearance;
  }

  function twoCellPhysicalFidelityMatrix(state, options = {}) {
    const controls = twoCellBenchControls(state, options);
    const radii = stateRadii(state, options);
    const nominalBacklash = finiteNumber(state.grid.backlash, 0.1);
    const nominalHole = radii.holeRadius;
    const looseHole = Math.max(nominalHole, controls.holeSweepMax);
    const backlashAxis = uniqueFloats(options.backlashValues || [0, nominalBacklash, Math.max(0.25, 3 * nominalBacklash)]);
    const holeAxis = uniqueFloats(options.holeRadii || [radii.pinRadius, nominalHole, looseHole]);
    const lockModes = Array.from(options.lockModes || ["free", "right_state_locked", "right_position_locked", "left_free"]);
    const actuationCases = Array.from(options.actuationCases || defaultFidelityActuationCases(controls));
    const initialAlpha = finiteNumber(state.grid.initialAlpha, 1);
    const rows = [];
    const connectorRows = [];
    let index = 0;
    for (const backlash of backlashAxis) {
      for (const holeRadius of holeAxis) {
        const caseBaseState = stateForTwoCellCase(state, { backlash, holeRadius });
        for (const lockMode of lockModes) {
          const lockBase = controlsForLockMode(controls, lockMode);
          for (const actuation of actuationCases) {
            const caseControls = controlsForLockMode(
              twoCellBenchControls(state, {
                ...lockBase,
                alphaCommand: finiteNumber(actuation.alpha, 0),
                zCommand: finiteNumber(actuation.z, 0),
              }),
              lockMode
            );
            const contact = twoCellConnectorContactReport(caseBaseState, caseControls);
            const bench = contact.bench;
            const left = bench.cells[0];
            const right = bench.cells[1];
            const lockPass = matrixLockPass(lockMode, bench);
            const row = {
              index,
              actuationCase: String(actuation.actuationCase),
              lockMode,
              backlash: caseBaseState.grid.backlash,
              alphaCommand: caseControls.alphaCommand,
              zCommand: caseControls.zCommand,
              gravityForce: finiteNumber(actuation.gravity, 0.025),
              pinRadius: bench.dimensions.pinRadius,
              holeRadius: bench.dimensions.holeRadius,
              clearance: bench.dimensions.pinHoleClearance,
              effectiveBacklash: bench.deadZones.effectiveBacklash,
              pinHoleClearanceMm: contact.dimensions.pinHoleClearanceMm,
              solverSuccess: true,
              lockPass,
              positionLocksHold: lockMode !== "right_position_locked" || lockPass,
              stateLocksHold: lockMode !== "right_state_locked" || lockPass,
              rightLockedVerticalMotionAllowed: lockMode === "right_state_locked" && Math.abs(finiteNumber(right.center?.z, 0)) > 1e-9,
              leftX: left.center.x,
              leftZ: left.center.z,
              rightX: right.center.x,
              rightZ: right.center.z,
              leftAlpha: left.alpha,
              rightAlpha: right.alpha,
              leftAlphaDeltaFromInitial: left.alpha - initialAlpha,
              rightAlphaDeltaFromInitial: right.alpha - initialAlpha,
              leftTheta: left.theta,
              rightTheta: right.theta,
              verticalShear: bench.connector.verticalShear,
              axialStrain: bench.connector.axialStrain,
              contactMode: bench.connector.contactMode,
              connectorActiveContactCount: contact.summary.activeContactCount,
              connectorMaxLateralSlipMm: contact.summary.maxLateralSlipMm,
              connectorMaxVerticalSlipMm: contact.summary.maxVerticalSlipMm,
              connectorMaxLateralExcessMm: contact.summary.maxLateralExcessMm,
              connectorMaxVerticalExcessMm: contact.summary.maxVerticalExcessMm,
              connectorTotalContactPenalty: contact.summary.totalContactPenalty,
              totalEnergy: bench.energy.totalEnergy,
              externalMeasurementRequired: true,
            };
            rows.push(row);
            for (const connector of contact.connectors || []) {
              connectorRows.push({
                index,
                actuationCase: row.actuationCase,
                lockMode,
                backlash: row.backlash,
                holeRadius: row.holeRadius,
                pinHoleClearanceMm: row.pinHoleClearanceMm,
                connector: connector.connector,
                leftSite: connector.leftSite,
                rightSite: connector.rightSite,
                leftPositionMm: connector.leftPositionMm,
                rightPositionMm: connector.rightPositionMm,
                slipMm: connector.slipMm,
                lateralSlipMm: connector.lateralSlipMm,
                verticalSlipMm: connector.verticalSlipMm,
                totalSlipMm: connector.totalSlipMm,
                lateralExcessMm: connector.lateralExcessMm,
                verticalExcessMm: connector.verticalExcessMm,
                contactPenalty: connector.contactPenalty,
                contactMode: connector.contactMode,
              });
            }
            index += 1;
          }
        }
      }
    }

    const nominalQuery = {
      actuationCase: "contract_lift",
      lockMode: "left_free",
      backlash: nominalBacklash,
    };
    const clearanceRows = holeAxis
      .map((holeRadius) => matrixRowLookup(rows, { ...nominalQuery, holeRadius }))
      .filter(Boolean);
    const backlashRows = backlashAxis
      .map((backlash) => matrixRowLookup(rows, { actuationCase: "contract_lift", lockMode: "left_free", backlash, holeRadius: nominalHole }))
      .filter(Boolean);
    const stateLockRows = rows.filter((row) => row.lockMode === "right_state_locked");
    const positionLockRows = rows.filter((row) => row.lockMode === "right_position_locked");
    const contract = matrixRowLookup(rows, { actuationCase: "alpha_contract", lockMode: "free", backlash: nominalBacklash, holeRadius: nominalHole });
    const expand = matrixRowLookup(rows, { actuationCase: "alpha_expand", lockMode: "free", backlash: nominalBacklash, holeRadius: nominalHole });
    const lift = matrixRowLookup(rows, { actuationCase: "z_lift", lockMode: "free", backlash: nominalBacklash, holeRadius: nominalHole });
    const push = matrixRowLookup(rows, { actuationCase: "z_pushdown", lockMode: "free", backlash: nominalBacklash, holeRadius: nominalHole });
    const clearanceReliefReferenceSlipMm = clearanceRows.length ? finiteNumber(clearanceRows[0].connectorMaxVerticalSlipMm, 0) : 0;
    const clearanceReliefPenaltyMm2 = clearanceRows.map((row) => {
      const residual = Math.max(0, clearanceReliefReferenceSlipMm - finiteNumber(row.pinHoleClearanceMm, 0));
      return 0.5 * residual * residual;
    });
    const clearanceTrend = {
      holeRadii: clearanceRows.map((row) => row.holeRadius),
      clearancesMm: clearanceRows.map((row) => row.pinHoleClearanceMm),
      maxVerticalSlipMm: clearanceRows.map((row) => row.connectorMaxVerticalSlipMm),
      maxVerticalExcessMm: clearanceRows.map((row) => row.connectorMaxVerticalExcessMm),
      contactPenalty: clearanceRows.map((row) => row.connectorTotalContactPenalty),
      clearanceReliefReferenceSlipMm,
      clearanceReliefPenaltyMm2,
      verticalFreePlayUtilization: clearanceRows.map((row) =>
        freePlayUtilization(row.connectorMaxVerticalSlipMm, row.pinHoleClearanceMm)
      ),
      clearanceNondecreasing: nondecreasing(clearanceRows.map((row) => row.pinHoleClearanceMm)),
      verticalExcessNonincreasing: nonincreasing(clearanceRows.map((row) => row.connectorMaxVerticalExcessMm)),
      penaltyNonincreasing: nonincreasing(clearanceReliefPenaltyMm2),
    };
    const backlashTrend = {
      backlashValues: backlashRows.map((row) => row.backlash),
      leftAlphaDelta: backlashRows.map((row) => row.leftAlphaDeltaFromInitial),
      leftZ: backlashRows.map((row) => row.leftZ),
      alphaResponseNonincreasing: nonincreasing(backlashRows.map((row) => Math.abs(row.leftAlphaDeltaFromInitial))),
      verticalResponseMateriallyFlat:
        Math.max(0, ...backlashRows.map((row) => Math.abs(row.leftZ))) -
          Math.min(...backlashRows.map((row) => Math.abs(row.leftZ))) <=
        1e-6,
    };
    const lockChecks = {
      stateLockAlphaHold: stateLockRows.every((row) => Math.abs(row.rightAlpha - initialAlpha) <= 1e-9),
      stateLockVerticalMotionSamples: stateLockRows.filter((row) => Math.abs(row.rightZ) > 1e-9).length,
      positionLockFixtureHold: positionLockRows.every(
        (row) => Math.abs(row.rightX - finiteNumber(state.grid.cellSize, 1)) <= 1e-9 && Math.abs(row.rightZ) <= 1e-9
      ),
    };
    const actuationPolarity = {
      contractDecreasesAlpha: contract ? contract.rightAlpha < initialAlpha : false,
      expandIncreasesAlpha: expand ? expand.rightAlpha > initialAlpha : false,
      positiveZLifts: lift ? lift.rightZ > 0 : false,
      negativeZPushesDown: push ? push.rightZ < 0 : false,
    };
    const solverSuccessCount = rows.filter((row) => row.solverSuccess).length;
    const internalMatrixReady =
      solverSuccessCount === rows.length &&
      clearanceTrend.clearanceNondecreasing &&
      clearanceTrend.penaltyNonincreasing &&
      backlashTrend.alphaResponseNonincreasing &&
      lockChecks.stateLockAlphaHold &&
      lockChecks.positionLockFixtureHold &&
      Object.values(actuationPolarity).every(Boolean);
    return {
      schema: TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA,
      model: "browser-two-cell-radius-backlash-lock-actuation-fidelity-matrix",
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadLayout: cadRadCellLayout(state, options),
      axes: {
        backlashValues: backlashAxis,
        holeRadii: holeAxis,
        lockModes,
        actuationCases: actuationCases.map((item) => item.actuationCase),
      },
      rows,
      connectorRows,
      summary: {
        status: internalMatrixReady ? "reduced-physics-fidelity-matrix-ready-needs-external-data" : "reduced-physics-fidelity-matrix-review-needed",
        internalMatrixReady,
        physicalAccuracyValidated: false,
        rowCount: rows.length,
        connectorRowCount: connectorRows.length,
        solverSuccessCount,
        solverFailureCount: rows.length - solverSuccessCount,
        maxVerticalSlipMm: Math.max(0, ...rows.map((row) => row.connectorMaxVerticalSlipMm)),
        maxVerticalExcessMm: Math.max(0, ...rows.map((row) => row.connectorMaxVerticalExcessMm)),
        maxContactPenalty: Math.max(0, ...rows.map((row) => row.connectorTotalContactPenalty)),
        clearanceTrend,
        backlashTrend,
        lockChecks,
        actuationPolarity,
        missingEvidence: [
          "filledSegmentedCadIntake",
          "externalRigidBodyContactRun",
          "filledConnectorMeasurements",
          "benchCoordinates",
          "measuredPinHoleClearance",
          "measuredFriction",
          "measuredContactStiffness",
        ],
      },
      claimLabels: {
        simulation: "browser reduced two-cell matrix over radius, backlash, locks, and actuation",
        physicalAccuracy: "proxy only until segmented CAD/contact and bench measurements are attached",
      },
    };
  }

  function exportTwoCellPhysicalFidelityMatrixJson(state, options = {}) {
    return JSON.stringify(twoCellPhysicalFidelityMatrix(state, options), null, 2);
  }

  function exportTwoCellPhysicalFidelityMatrixCsv(state, options = {}) {
    const matrix = twoCellPhysicalFidelityMatrix(state, options);
    const columns = [
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
      "pinHoleClearanceMm",
      "lockPass",
      "rightLockedVerticalMotionAllowed",
      "leftZ",
      "rightZ",
      "leftAlphaDeltaFromInitial",
      "rightAlphaDeltaFromInitial",
      "contactMode",
      "connectorActiveContactCount",
      "connectorMaxVerticalSlipMm",
      "connectorMaxVerticalExcessMm",
      "connectorTotalContactPenalty",
      "totalEnergy",
      "externalMeasurementRequired",
    ];
    return [columns.join(","), ...matrix.rows.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function contactPhaseForRow(row, tolerance = 1e-7) {
    const lockMode = String(row.lockMode || "free");
    const contactMode = String(row.contactMode || "inside-clearance");
    const verticalExcess = Math.max(Math.abs(finiteNumber(row.verticalExcess, 0)), Math.abs(finiteNumber(row.connectorMaxVerticalExcessMm, 0)));
    const axialExcess = Math.max(
      Math.abs(finiteNumber(row.axialExcess, 0)),
      Math.abs(finiteNumber(row.connectorMaxLateralExcessMm, 0)),
      Math.abs(finiteNumber(row.axialStrain, 0))
    );
    const activeContact =
      Math.round(finiteNumber(row.connectorActiveContactCount, 0)) > 0 || finiteNumber(row.connectorTotalContactPenalty, 0) > tolerance;
    const verticalContact = verticalExcess > tolerance || contactMode.includes("pin-hole-vertical-contact");
    const axialContact = axialExcess > tolerance || contactMode.includes("axial-backlash-contact");
    const gravitySag =
      finiteNumber(row.gravityForce, 0) > tolerance &&
      Math.abs(finiteNumber(row.zCommand, 0)) <= tolerance &&
      Math.min(finiteNumber(row.leftZ, 0), finiteNumber(row.rightZ, 0)) < -tolerance;
    if (lockMode === "right_position_locked") return "position-locked";
    if (lockMode === "right_state_locked") {
      return row.rightLockedVerticalMotionAllowed === true && Math.abs(finiteNumber(row.rightZ, 0)) > tolerance
        ? "state-locked-vertical-free"
        : "state-locked";
    }
    if (axialContact && verticalContact) return "axial-and-vertical-contact";
    if (verticalContact) return "vertical-contact";
    if (axialContact || activeContact) return "axial-contact";
    if (gravitySag) return "gravity-sag";
    return "free-play";
  }

  function phaseMapDigest(row, category) {
    return {
      caseId: fidelityMatrixCaseId(row),
      category,
      matrixIndex: row.index,
      phase: row.phase,
      actuationCase: row.actuationCase,
      lockMode: row.lockMode,
      backlash: row.backlash,
      holeRadius: row.holeRadius,
      pinHoleClearanceMm: row.pinHoleClearanceMm,
      effectiveBacklash: row.effectiveBacklash ?? row.backlash,
      effectiveHoleRadius: row.effectiveHoleRadius ?? row.holeRadius,
      effectivePinHoleClearanceMm: row.effectivePinHoleClearanceMm ?? row.pinHoleClearanceMm,
      alphaCommand: row.alphaCommand,
      zCommand: row.zCommand,
      gravityForce: row.gravityForce,
      leftZ: row.leftZ,
      rightZ: row.rightZ,
      rightAlpha: row.rightAlpha,
      contactMode: row.contactMode,
      connectorActiveContactCount: row.connectorActiveContactCount,
      connectorMaxVerticalSlipMm: row.connectorMaxVerticalSlipMm,
      connectorMaxVerticalExcessMm: row.connectorMaxVerticalExcessMm,
      connectorTotalContactPenalty: row.connectorTotalContactPenalty,
      lockPass: row.lockPass,
      physicalAccuracyValidated: false,
    };
  }

  function phasePriorityCases(rows, limit = 12) {
    const categories = [
      [
        "vertical-contact",
        (row) => row.phase === "vertical-contact" || row.phase === "axial-and-vertical-contact",
        (row) => [
          Math.abs(finiteNumber(row.connectorMaxVerticalExcessMm, 0)),
          Math.abs(finiteNumber(row.connectorMaxVerticalSlipMm, 0)),
          Math.abs(finiteNumber(row.rightZ, 0)),
        ],
      ],
      [
        "axial-contact",
        (row) => row.phase === "axial-contact" || row.phase === "axial-and-vertical-contact",
        (row) => [
          Math.abs(finiteNumber(row.connectorMaxLateralExcessMm, 0)),
          Math.abs(finiteNumber(row.axialExcess, row.axialStrain || 0)),
          Math.abs(finiteNumber(row.connectorTotalContactPenalty, 0)),
        ],
      ],
      ["gravity-sag", (row) => row.phase === "gravity-sag", (row) => [Math.abs(finiteNumber(row.leftZ, 0)), Math.abs(finiteNumber(row.rightZ, 0))]],
      [
        "state-lock-vertical-free",
        (row) => row.phase === "state-locked-vertical-free",
        (row) => [Math.abs(finiteNumber(row.rightZ, 0)), Math.abs(finiteNumber(row.zCommand, 0))],
      ],
      [
        "position-lock",
        (row) => row.phase === "position-locked",
        (row) => [Math.abs(finiteNumber(row.connectorMaxVerticalSlipMm, 0)), Math.abs(finiteNumber(row.connectorTotalContactPenalty, 0))],
      ],
      [
        "free-play-boundary",
        (row) => row.phase === "free-play",
        (row) => [Math.abs(finiteNumber(row.connectorMaxVerticalSlipMm, 0)), Math.abs(finiteNumber(row.leftAlphaDeltaFromInitial, 0))],
      ],
    ];
    const seen = new Set();
    const out = [];
    const maxCount = Math.max(1, Math.round(finiteNumber(limit, 12)));
    for (const [category, predicate, sortKey] of categories) {
      const matches = rows
        .filter(predicate)
        .slice()
        .sort((left, right) => {
          const leftKey = sortKey(left);
          const rightKey = sortKey(right);
          for (let i = 0; i < Math.max(leftKey.length, rightKey.length); i += 1) {
            const delta = finiteNumber(rightKey[i], 0) - finiteNumber(leftKey[i], 0);
            if (Math.abs(delta) > 1e-12) return delta;
          }
          return finiteNumber(left.index, 0) - finiteNumber(right.index, 0);
        })
        .slice(0, Math.max(1, Math.min(3, maxCount)));
      for (const row of matches) {
        const index = Math.round(finiteNumber(row.index, 0));
        if (seen.has(index)) continue;
        seen.add(index);
        out.push(phaseMapDigest(row, category));
        if (out.length >= maxCount) return out;
      }
    }
    return out;
  }

  function twoCellContactPhaseMap(state, options = {}) {
    const tolerance = finiteNumber(options.tolerance, 1e-7);
    const matrix = twoCellPhysicalFidelityMatrix(state, options);
    const rows = (matrix.rows || []).map((source) => {
      const phase = contactPhaseForRow(source, tolerance);
      return {
        caseId: fidelityMatrixCaseId(source),
        ...source,
        phase,
        axialContact: phase === "axial-contact" || phase === "axial-and-vertical-contact",
        verticalContact: phase === "vertical-contact" || phase === "axial-and-vertical-contact",
        lockPhase: phase === "position-locked" || phase === "state-locked" || phase === "state-locked-vertical-free",
        freePlayPhase: phase === "free-play",
        gravitySagPhase: phase === "gravity-sag",
        physicalAccuracyValidated: false,
      };
    });
    const phaseCounts = rows.reduce((counts, row) => {
      counts[row.phase] = (counts[row.phase] || 0) + 1;
      return counts;
    }, {});
    const dominantPhase =
      Object.keys(phaseCounts)
        .sort()
        .sort((left, right) => phaseCounts[right] - phaseCounts[left])[0] || "none";
    const lockRows = rows.filter((row) => row.lockPhase);
    const activeContactCaseCount = rows.filter((row) => Math.round(finiteNumber(row.connectorActiveContactCount, 0)) > 0).length;
    return {
      schema: TWO_CELL_CONTACT_PHASE_MAP_SCHEMA,
      model: "browser-two-cell-contact-lock-phase-map-from-fidelity-matrix",
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadLayout: matrix.cadLayout,
      axes: matrix.axes,
      matrixSummary: matrix.summary,
      rows,
      summary: {
        status:
          matrix.summary?.solverFailureCount === 0
            ? "contact-phase-map-ready-needs-bench-or-external-validation"
            : "contact-phase-map-review-solver-failures",
        rowCount: rows.length,
        phaseCounts,
        dominantPhase,
        solverFailureCount: matrix.summary?.solverFailureCount || 0,
        lockCaseCount: lockRows.length,
        lockFailureCount: lockRows.filter((row) => row.lockPass !== true).length,
        activeContactCaseCount,
        freePlayCaseCount: phaseCounts["free-play"] || 0,
        gravitySagCaseCount: phaseCounts["gravity-sag"] || 0,
        maxVerticalSlipMm: Math.max(0, ...rows.map((row) => Math.abs(finiteNumber(row.connectorMaxVerticalSlipMm, 0)))),
        maxVerticalExcessMm: Math.max(0, ...rows.map((row) => Math.abs(finiteNumber(row.connectorMaxVerticalExcessMm, 0)))),
        maxAxialExcess: Math.max(
          0,
          ...rows.map((row) => Math.max(Math.abs(finiteNumber(row.axialExcess, 0)), Math.abs(finiteNumber(row.axialStrain, 0))))
        ),
        maxContactPenalty: Math.max(0, ...rows.map((row) => Math.abs(finiteNumber(row.connectorTotalContactPenalty, 0)))),
        measurementPriority: phasePriorityCases(rows, options.priorityLimit || 12),
        physicalAccuracyValidated: false,
        remainingEvidence: matrix.summary?.missingEvidence || [],
      },
      claimBoundary: {
        allowedClaim: "interpretable reduced-model phase classification for selecting two-cell tests",
        blockedClaim: "exact real-cell phase boundaries without segmented CAD contact and bench validation",
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellContactPhaseMapJson(state, options = {}) {
    return JSON.stringify(twoCellContactPhaseMap(state, options), null, 2);
  }

  function exportTwoCellContactPhaseMapCsv(state, options = {}) {
    const phaseMap = twoCellContactPhaseMap(state, options);
    const columns = [
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
    ];
    return [columns.join(","), ...phaseMap.rows.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function numericAxis(start, stop, steps) {
    const count = Math.max(2, Math.round(finiteNumber(steps, 9)));
    const a = finiteNumber(start, 0);
    const b = finiteNumber(stop, a);
    if (Math.abs(a - b) <= 1e-12) return [a];
    return Array.from({ length: count }, (_, index) => a + ((b - a) * index) / (count - 1));
  }

  function fidelityActuationCase(controls, actuationCase) {
    const cases = defaultFidelityActuationCases(controls);
    const selected = cases.find((item) => String(item.actuationCase) === String(actuationCase));
    if (!selected) throw new Error(`unknown actuation case ${actuationCase}`);
    return { ...selected };
  }

  function twoCellRadiusBacklashPhaseDiagram(state, options = {}) {
    const controls = twoCellBenchControls(state, options);
    const radii = stateRadii(state, options);
    const actuationCase = options.actuationCase || "contract_lift";
    const lockMode = options.lockMode || "free";
    const actuation = fidelityActuationCase(controls, actuationCase);
    const holeMin = finiteNumber(options.holeMin, radii.pinRadius);
    const holeMax = finiteNumber(options.holeMax, Math.max(radii.holeRadius, controls.holeSweepMax));
    const backlashMin = finiteNumber(options.backlashMin, 0);
    const backlashMax = finiteNumber(options.backlashMax, Math.max(0.25, 3 * finiteNumber(state.grid.backlash, 0.1)));
    const holeAxis = uniqueFloats(numericAxis(holeMin, holeMax, options.holeSteps || 9));
    const backlashAxis = uniqueFloats(numericAxis(backlashMin, backlashMax, options.backlashSteps || 9));
    const holeRadiusBias = finiteNumber(options.holeRadiusBias ?? options.hole_radius_bias, 0);
    const backlashBias = finiteNumber(options.backlashBias ?? options.backlash_bias, 0);
    const effectiveHoleAxis = holeAxis.map((holeRadius) => Math.max(radii.pinRadius, finiteNumber(holeRadius, radii.pinRadius) + holeRadiusBias));
    const effectiveBacklashAxis = backlashAxis.map((backlash) => Math.max(0, finiteNumber(backlash, 0) + backlashBias));
    const matrix = twoCellPhysicalFidelityMatrix(state, {
      ...options,
      backlashValues: effectiveBacklashAxis,
      holeRadii: effectiveHoleAxis,
      lockModes: [lockMode],
      actuationCases: [actuation],
    });
    const byCoord = new Map();
    for (const row of matrix.rows || []) byCoord.set(`${finiteNumber(row.backlash, 0)}:${finiteNumber(row.holeRadius, 0)}`, row);
    const rows = [];
    const phaseGrid = [];
    for (let rowIndex = 0; rowIndex < backlashAxis.length; rowIndex += 1) {
      const backlash = backlashAxis[rowIndex];
      const effectiveBacklash = effectiveBacklashAxis[rowIndex];
      const gridRow = [];
      for (let columnIndex = 0; columnIndex < holeAxis.length; columnIndex += 1) {
        const holeRadius = holeAxis[columnIndex];
        const effectiveHoleRadius = effectiveHoleAxis[columnIndex];
        const source = byCoord.get(`${finiteNumber(effectiveBacklash, 0)}:${finiteNumber(effectiveHoleRadius, 0)}`);
        if (!source) continue;
        const phase = contactPhaseForRow(source, finiteNumber(options.tolerance, 1e-7));
        const row = {
          caseId: `PD_${String(rows.length).padStart(3, "0")}`,
          matrixIndex: source.index,
          ...source,
          backlash,
          holeRadius,
          effectiveBacklash: source.backlash,
          effectiveHoleRadius: source.holeRadius,
          effectivePinHoleClearanceMm: source.pinHoleClearanceMm,
          phase,
          axialContact: phase === "axial-contact" || phase === "axial-and-vertical-contact",
          verticalContact: phase === "vertical-contact" || phase === "axial-and-vertical-contact",
          lockPhase: phase === "position-locked" || phase === "state-locked" || phase === "state-locked-vertical-free",
          freePlayPhase: phase === "free-play",
          gravitySagPhase: phase === "gravity-sag",
          physicalAccuracyValidated: false,
        };
        rows.push(row);
        gridRow.push(phase);
      }
      phaseGrid.push(gridRow);
    }
    const phaseCounts = rows.reduce((counts, row) => {
      counts[row.phase] = (counts[row.phase] || 0) + 1;
      return counts;
    }, {});
    const transitions = [];
    for (let rowIndex = 0; rowIndex < phaseGrid.length; rowIndex += 1) {
      const phases = phaseGrid[rowIndex] || [];
      for (let columnIndex = 1; columnIndex < phases.length; columnIndex += 1) {
        if (phases[columnIndex] !== phases[columnIndex - 1]) {
          transitions.push({
            axis: "holeRadius",
            backlash: backlashAxis[rowIndex],
            leftHoleRadius: holeAxis[columnIndex - 1],
            rightHoleRadius: holeAxis[columnIndex],
            fromPhase: phases[columnIndex - 1],
            toPhase: phases[columnIndex],
          });
        }
      }
    }
    for (let columnIndex = 0; columnIndex < holeAxis.length; columnIndex += 1) {
      for (let rowIndex = 1; rowIndex < phaseGrid.length; rowIndex += 1) {
        const previous = phaseGrid[rowIndex - 1]?.[columnIndex];
        const current = phaseGrid[rowIndex]?.[columnIndex];
        if (current !== previous) {
          transitions.push({
            axis: "backlash",
            holeRadius: holeAxis[columnIndex],
            lowerBacklash: backlashAxis[rowIndex - 1],
            upperBacklash: backlashAxis[rowIndex],
            fromPhase: previous,
            toPhase: current,
          });
        }
      }
    }
    const dominantPhase =
      Object.keys(phaseCounts)
        .sort()
        .sort((left, right) => phaseCounts[right] - phaseCounts[left])[0] || "none";
    return {
      schema: TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA,
      model: "browser-dense-two-cell-hole-radius-backlash-contact-phase-diagram",
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadLayout: matrix.cadLayout,
      axes: {
        holeRadii: holeAxis,
        backlashValues: backlashAxis,
        effectiveHoleRadii: effectiveHoleAxis,
        effectiveBacklashValues: effectiveBacklashAxis,
        holeRadiusBias,
        backlashBias,
        actuationCase,
        lockMode,
      },
      actuation: {
        actuationCase: actuation.actuationCase,
        alphaCommand: actuation.alpha,
        zCommand: actuation.z,
        gravityForce: actuation.gravity,
      },
      phaseGrid,
      rows,
      transitions,
      matrixSummary: matrix.summary,
      summary: {
        status: "radius-backlash-phase-diagram-ready-needs-bench-or-external-validation",
        rowCount: rows.length,
        holeStepCount: holeAxis.length,
        backlashStepCount: backlashAxis.length,
        phaseCounts,
        dominantPhase,
        transitionCount: transitions.length,
        activeContactCaseCount: rows.filter((row) => Math.round(finiteNumber(row.connectorActiveContactCount, 0)) > 0).length,
        maxVerticalSlipMm: Math.max(0, ...rows.map((row) => Math.abs(finiteNumber(row.connectorMaxVerticalSlipMm, 0)))),
        maxVerticalExcessMm: Math.max(0, ...rows.map((row) => Math.abs(finiteNumber(row.connectorMaxVerticalExcessMm, 0)))),
        maxContactPenalty: Math.max(0, ...rows.map((row) => Math.abs(finiteNumber(row.connectorTotalContactPenalty, 0)))),
        physicalAccuracyValidated: false,
        remainingEvidence: matrix.summary?.missingEvidence || [],
      },
      claimBoundary: {
        allowedClaim: "dense reduced-model phase diagram for choosing hole-radius and backlash experiments",
        blockedClaim: "real transition boundaries until segmented CAD contact and bench sweeps are compared",
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellRadiusBacklashPhaseDiagramJson(state, options = {}) {
    return JSON.stringify(twoCellRadiusBacklashPhaseDiagram(state, options), null, 2);
  }

  function exportTwoCellRadiusBacklashPhaseDiagramCsv(state, options = {}) {
    const diagram = twoCellRadiusBacklashPhaseDiagram(state, options);
    const columns = [
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
    ];
    return [columns.join(","), ...diagram.rows.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function phaseDiagramCoordKey(backlash, holeRadius) {
    return `${finiteNumber(backlash, 0).toFixed(12)}:${finiteNumber(holeRadius, 0).toFixed(12)}`;
  }

  function transitionCaseDigest(row, category, rank) {
    return {
      ...phaseMapDigest(row, category),
      rank,
      phaseDiagramCaseId: row.caseId || "",
      measurementTarget: "two-cell center pose plus upper/middle/lower connector marker coordinates",
      observedPhase: "",
      observedRightZ: "",
      observedRightAlpha: "",
      observedConnectorMaxVerticalSlipMm: "",
      observedConnectorMaxVerticalExcessMm: "",
      observedLockHeld: "",
      notes: "",
    };
  }

  function twoCellRadiusBacklashTransitionReport(state, options = {}) {
    const diagram = twoCellRadiusBacklashPhaseDiagram(state, options);
    const rowsByCoord = new Map();
    for (const row of diagram.rows || []) rowsByCoord.set(phaseDiagramCoordKey(row.backlash, row.holeRadius), row);
    const transitionBrackets = [];
    const candidateRows = new Map();
    for (let index = 0; index < (diagram.transitions || []).length; index += 1) {
      const transition = diagram.transitions[index];
      let left = null;
      let right = null;
      let bracket = null;
      if (transition.axis === "holeRadius") {
        left = rowsByCoord.get(phaseDiagramCoordKey(transition.backlash, transition.leftHoleRadius));
        right = rowsByCoord.get(phaseDiagramCoordKey(transition.backlash, transition.rightHoleRadius));
        const lowerHoleRadius = finiteNumber(transition.leftHoleRadius, 0);
        const upperHoleRadius = finiteNumber(transition.rightHoleRadius, lowerHoleRadius);
        bracket = {
          transitionId: `TR_${String(index + 1).padStart(3, "0")}`,
          axis: "holeRadius",
          fixedBacklash: transition.backlash,
          lowerHoleRadius,
          upperHoleRadius,
          midHoleRadius: 0.5 * (lowerHoleRadius + upperHoleRadius),
          bracketWidth: Math.abs(upperHoleRadius - lowerHoleRadius),
          fromPhase: transition.fromPhase,
          toPhase: transition.toPhase,
        };
      } else {
        left = rowsByCoord.get(phaseDiagramCoordKey(transition.lowerBacklash, transition.holeRadius));
        right = rowsByCoord.get(phaseDiagramCoordKey(transition.upperBacklash, transition.holeRadius));
        const lowerBacklash = finiteNumber(transition.lowerBacklash, 0);
        const upperBacklash = finiteNumber(transition.upperBacklash, lowerBacklash);
        bracket = {
          transitionId: `TR_${String(index + 1).padStart(3, "0")}`,
          axis: "backlash",
          fixedHoleRadius: transition.holeRadius,
          lowerBacklash,
          upperBacklash,
          midBacklash: 0.5 * (lowerBacklash + upperBacklash),
          bracketWidth: Math.abs(upperBacklash - lowerBacklash),
          fromPhase: transition.fromPhase,
          toPhase: transition.toPhase,
        };
      }
      bracket.leftCaseId = left?.caseId || "";
      bracket.rightCaseId = right?.caseId || "";
      bracket.leftMatrixCaseId = left ? fidelityMatrixCaseId(left) : "";
      bracket.rightMatrixCaseId = right ? fidelityMatrixCaseId(right) : "";
      bracket.physicalAccuracyValidated = false;
      transitionBrackets.push(bracket);
      if (left) candidateRows.set(String(left.caseId), ["transition-lower-side", left]);
      if (right) candidateRows.set(String(right.caseId), ["transition-upper-side", right]);
    }
    const bySlip = [...(diagram.rows || [])].sort(
      (left, right) => Math.abs(finiteNumber(right.connectorMaxVerticalSlipMm, 0)) - Math.abs(finiteNumber(left.connectorMaxVerticalSlipMm, 0))
    );
    const byPenalty = [...(diagram.rows || [])].sort(
      (left, right) => Math.abs(finiteNumber(right.connectorTotalContactPenalty, 0)) - Math.abs(finiteNumber(left.connectorTotalContactPenalty, 0))
    );
    if (bySlip[0]) candidateRows.set(String(bySlip[0].caseId), candidateRows.get(String(bySlip[0].caseId)) || ["max-vertical-slip", bySlip[0]]);
    if (byPenalty[0]) candidateRows.set(String(byPenalty[0].caseId), candidateRows.get(String(byPenalty[0].caseId)) || ["max-contact-penalty", byPenalty[0]]);
    if (!candidateRows.size && diagram.rows?.[0]) candidateRows.set(String(diagram.rows[0].caseId), ["no-transition-baseline", diagram.rows[0]]);
    const priorityLimit = Math.max(1, Math.round(finiteNumber(options.priorityLimit, 16)));
    const measurementCases = [...candidateRows.values()]
      .slice(0, priorityLimit)
      .map(([category, row], index) => transitionCaseDigest(row, category, index + 1));
    const bracketWidths = transitionBrackets.map((item) => finiteNumber(item.bracketWidth, 0));
    const phasePairCounts = transitionBrackets.reduce((counts, item) => {
      const key = `${item.fromPhase}->${item.toPhase}`;
      counts[key] = (counts[key] || 0) + 1;
      return counts;
    }, {});
    return {
      schema: TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA,
      model: "browser-two-cell-radius-backlash-transition-bracket-calibration-report",
      diagram: {
        schema: diagram.schema,
        axes: diagram.axes,
        summary: diagram.summary,
      },
      cadReference: CAD_RAD_CELL_REFERENCE,
      transitionBrackets,
      measurementCases,
      measurementProtocol: {
        goal: "Measure the two-cell cases adjacent to predicted phase changes and fit the real pin-hole backlash boundary.",
        requiredObservables: [
          "left and right cell center xyz",
          "left and right cell alpha/theta",
          "upper, middle, and lower connector marker xyz",
          "observed contact/free-play phase",
          "whether state locks and position locks held",
        ],
        notes: [
          "Use the same pin radius, hole radius, and backlash labels recorded in each case.",
          "Do not treat the reduced-model bracket midpoint as a real threshold until measured.",
          "Repeat transition-adjacent cases after any CAD export or hardware change.",
        ],
      },
      summary: {
        status: "transition-report-ready-needs-bench-sweep-validation",
        transitionBracketCount: transitionBrackets.length,
        measurementCaseCount: measurementCases.length,
        phasePairCounts,
        minBracketWidth: bracketWidths.length ? Math.min(...bracketWidths) : 0,
        maxBracketWidth: Math.max(0, ...bracketWidths),
        physicalAccuracyValidated: false,
        remainingEvidence: [
          "measured transition-adjacent two-cell coordinates",
          "segmented one-cell and two-cell CAD contact geometry",
          "pin-hole friction and normal stiffness calibration",
          "independent holdout transition sweeps",
        ],
      },
      claimBoundary: {
        allowedClaim: "reduced-model transition brackets and prioritized bench cases",
        blockedClaim: "real RAD phase-transition law without measured transition sweeps and exact contact geometry",
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellRadiusBacklashTransitionReportJson(state, options = {}) {
    return JSON.stringify(twoCellRadiusBacklashTransitionReport(state, options), null, 2);
  }

  function exportTwoCellRadiusBacklashTransitionReportCsv(state, options = {}) {
    const report = twoCellRadiusBacklashTransitionReport(state, options);
    const columns = [
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
    ];
    return [columns.join(","), ...report.measurementCases.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function transitionMeasurementRows(input) {
    if (Array.isArray(input)) return input.filter((row) => row.caseId || row.phaseDiagramCaseId);
    if (input && Array.isArray(input.rows)) return transitionMeasurementRows(input.rows);
    if (input && Array.isArray(input.measurements)) return transitionMeasurementRows(input.measurements);
    if (input && Array.isArray(input.measurementCases)) return transitionMeasurementRows(input.measurementCases);
    if (typeof input === "string") {
      const trimmed = input.trim();
      if (!trimmed) return [];
      if (trimmed.startsWith("{") || trimmed.startsWith("[")) return transitionMeasurementRows(JSON.parse(trimmed));
      return transitionMeasurementRows(parseCsv(trimmed));
    }
    return [];
  }

  function phaseLabel(value) {
    const normalized = String(value || "")
      .trim()
      .toLowerCase()
      .replace(/\+/g, " and ")
      .replace(/[_-]+/g, " ")
      .replace(/\s+/g, " ");
    const aliases = {
      "axial vertical contact": "axial-and-vertical-contact",
      "axial and vertical contact": "axial-and-vertical-contact",
      "axial contact": "axial-contact",
      "vertical contact": "vertical-contact",
      "position locked": "position-locked",
      "state locked": "state-locked",
      "state locked vertical free": "state-locked-vertical-free",
      "free play": "free-play",
      "gravity sag": "gravity-sag",
    };
    return aliases[normalized] || normalized.replace(/\s+/g, "-");
  }

  function boolOrNull(value) {
    const normalized = String(value || "").trim().toLowerCase().replace(/[_\s]+/g, "-");
    if (!normalized) return null;
    if (["true", "t", "yes", "y", "held", "pass", "passed", "1"].includes(normalized)) return true;
    if (["false", "f", "no", "n", "fell", "fail", "failed", "released", "not-held", "unheld", "lock-fell", "0"].includes(normalized)) return false;
    return null;
  }

  function transitionShiftVote(predicted, observedPhase, brackets) {
    if (!observedPhase || observedPhase === String(predicted.phase || "")) return "within-predicted-bracket";
    for (const bracket of brackets || []) {
      if (bracket.transitionSide === "lower" && observedPhase === String(bracket.toPhase || "")) {
        return `${bracket.axis}-boundary-lower-than-predicted`;
      }
      if (bracket.transitionSide === "upper" && observedPhase === String(bracket.fromPhase || "")) {
        return `${bracket.axis}-boundary-higher-than-predicted`;
      }
    }
    return "phase-mismatch-not-adjacent";
  }

  function compareTwoCellRadiusBacklashTransitionMeasurements(state, input, options = {}) {
    const measurements = transitionMeasurementRows(input);
    const report = twoCellRadiusBacklashTransitionReport(state, options);
    const byCase = new Map((report.measurementCases || []).map((row) => [String(row.caseId || ""), row]));
    const byPhaseCase = new Map((report.measurementCases || []).map((row) => [String(row.phaseDiagramCaseId || ""), row]));
    const bracketsByPhaseCase = new Map();
    for (const bracket of report.transitionBrackets || []) {
      for (const [field, side] of [
        ["leftCaseId", "lower"],
        ["rightCaseId", "upper"],
      ]) {
        const caseId = String(bracket[field] || "");
        if (!caseId) continue;
        if (!bracketsByPhaseCase.has(caseId)) bracketsByPhaseCase.set(caseId, []);
        bracketsByPhaseCase.get(caseId).push({ ...bracket, transitionSide: side });
      }
    }
    const scalarPairs = [
      ["observedRightZ", "rightZ", "z"],
      ["observedRightAlpha", "rightAlpha", "alpha"],
      ["observedConnectorMaxVerticalSlipMm", "connectorMaxVerticalSlipMm", "verticalSlip"],
      ["observedConnectorMaxVerticalExcessMm", "connectorMaxVerticalExcessMm", "verticalExcess"],
      ["observedConnectorTotalContactPenalty", "connectorTotalContactPenalty", "penalty"],
    ];
    const comparisonRows = [];
    const shiftVotes = {};
    const zSamples = [];
    const alphaSamples = [];
    const verticalSlipSamples = [];
    const verticalExcessSamples = [];
    let matchedRowCount = 0;
    let missingPredictionCount = 0;
    let missingObservationRowCount = 0;
    let observedScalarCount = 0;
    let squaredError = 0;
    let maxAbsError = 0;
    let phaseObservationCount = 0;
    let phaseMatchCount = 0;
    let lockObservationCount = 0;
    let lockMatchCount = 0;
    for (const measurement of measurements) {
      const phaseDiagramCaseId = String(measurement.phaseDiagramCaseId || "").trim();
      const caseId = String(measurement.caseId || "").trim();
      const predicted = (phaseDiagramCaseId && byPhaseCase.get(phaseDiagramCaseId)) || byCase.get(caseId);
      if (!predicted) {
        missingPredictionCount += 1;
        comparisonRows.push({ caseId, phaseDiagramCaseId, status: "missing-prediction" });
        continue;
      }
      matchedRowCount += 1;
      const observedPhase = phaseLabel(measurement.observedPhase || measurement.observedContactPhase || "");
      const predictedPhase = String(predicted.phase || "");
      let phaseMatches = null;
      let shiftVote = "unobserved-phase";
      if (observedPhase) {
        phaseObservationCount += 1;
        phaseMatches = observedPhase === predictedPhase;
        if (phaseMatches) phaseMatchCount += 1;
        shiftVote = transitionShiftVote(predicted, observedPhase, bracketsByPhaseCase.get(String(predicted.phaseDiagramCaseId || "")) || []);
        shiftVotes[shiftVote] = (shiftVotes[shiftVote] || 0) + 1;
      }
      const residuals = {};
      let rowScalarCount = 0;
      let rowMaxAbs = 0;
      for (const [observedField, predictedField, group] of scalarPairs) {
        const observed = numericOrNull(measurement[observedField]);
        const target = numericOrNull(predicted[predictedField]);
        if (observed === null || target === null) continue;
        const residual = observed - target;
        residuals[observedField] = residual;
        rowScalarCount += 1;
        observedScalarCount += 1;
        squaredError += residual * residual;
        rowMaxAbs = Math.max(rowMaxAbs, Math.abs(residual));
        maxAbsError = Math.max(maxAbsError, Math.abs(residual));
        if (group === "z") zSamples.push([target, observed]);
        if (group === "alpha") alphaSamples.push([target - finiteNumber(state.grid.initialAlpha, 1), observed - finiteNumber(state.grid.initialAlpha, 1)]);
        if (group === "verticalSlip") verticalSlipSamples.push([target, observed]);
        if (group === "verticalExcess") verticalExcessSamples.push([target, observed]);
      }
      const observedLock = boolOrNull(measurement.observedLockHeld || measurement.lockHeldObserved || "");
      const lockExpected = ["right_state_locked", "right_position_locked"].includes(String(predicted.lockMode || ""));
      let lockHeldMatches = null;
      if (observedLock !== null && lockExpected) {
        lockObservationCount += 1;
        lockHeldMatches = observedLock === true && predicted.lockPass !== false;
        if (lockHeldMatches) lockMatchCount += 1;
      }
      if (!observedPhase && rowScalarCount === 0 && observedLock === null) missingObservationRowCount += 1;
      comparisonRows.push({
        caseId: predicted.caseId,
        phaseDiagramCaseId: predicted.phaseDiagramCaseId,
        matrixIndex: predicted.matrixIndex,
        status: "matched",
        category: predicted.category,
        predictedPhase,
        observedPhase,
        phaseMatches,
        transitionShiftVote: shiftVote,
        backlash: predicted.backlash,
        holeRadius: predicted.holeRadius,
        effectiveBacklash: predicted.effectiveBacklash ?? predicted.backlash,
        effectiveHoleRadius: predicted.effectiveHoleRadius ?? predicted.holeRadius,
        effectivePinHoleClearanceMm: predicted.effectivePinHoleClearanceMm ?? predicted.pinHoleClearanceMm,
        predictedRightZ: predicted.rightZ,
        predictedRightAlpha: predicted.rightAlpha,
        predictedConnectorMaxVerticalSlipMm: predicted.connectorMaxVerticalSlipMm,
        predictedConnectorMaxVerticalExcessMm: predicted.connectorMaxVerticalExcessMm,
        observedScalarCount: rowScalarCount,
        maxAbsResidual: rowScalarCount ? rowMaxAbs : "",
        lockExpected,
        lockHeldObserved: observedLock,
        lockHeldMatches,
        residuals,
        notes: measurement.notes || "",
      });
    }
    const holeLowerVotes = finiteNumber(shiftVotes["holeRadius-boundary-lower-than-predicted"], 0);
    const holeHigherVotes = finiteNumber(shiftVotes["holeRadius-boundary-higher-than-predicted"], 0);
    const transitionShiftDirection = holeLowerVotes && !holeHigherVotes
      ? "effective-hole-transition-lower-than-reduced-model"
      : holeHigherVotes && !holeLowerVotes
      ? "effective-hole-transition-higher-than-reduced-model"
      : holeLowerVotes || holeHigherVotes
      ? "mixed-transition-shift"
      : "not-enough-transition-phase-evidence";
    const backlashLowerVotes = finiteNumber(shiftVotes["backlash-boundary-lower-than-predicted"], 0);
    const backlashHigherVotes = finiteNumber(shiftVotes["backlash-boundary-higher-than-predicted"], 0);
    const holeWidths = (report.transitionBrackets || [])
      .filter((item) => item.axis === "holeRadius")
      .map((item) => finiteNumber(item.bracketWidth, 0));
    const backlashWidths = (report.transitionBrackets || [])
      .filter((item) => item.axis === "backlash")
      .map((item) => finiteNumber(item.bracketWidth, 0));
    const holeBias = holeLowerVotes && !holeHigherVotes
      ? 0.5 * finiteNumber(mean(holeWidths), 0)
      : holeHigherVotes && !holeLowerVotes
      ? -0.5 * finiteNumber(mean(holeWidths), 0)
      : 0;
    const backlashBias = backlashLowerVotes && !backlashHigherVotes
      ? 0.5 * finiteNumber(mean(backlashWidths), 0)
      : backlashHigherVotes && !backlashLowerVotes
      ? -0.5 * finiteNumber(mean(backlashWidths), 0)
      : 0;
    const zFit = fitScale(zSamples);
    const alphaFit = fitScale(alphaSamples);
    const verticalSlipFit = fitScale(verticalSlipSamples);
    const verticalExcessFit = fitScale(verticalExcessSamples);
    const pinRadius = Math.max(0, finiteNumber(state.grid.pinRadius, 0.18));
    const holeRadius = Math.max(pinRadius, finiteNumber(state.grid.holeRadius, 0.225));
    const currentClearance = Math.max(0, holeRadius - pinRadius);
    const proposedClearance = Math.max(0, currentClearance * bounded(verticalSlipFit.estimate, 0.25, 4, 1) + holeBias);
    const proposedReducedProxyUpdates = {
      couplingGain: bounded(finiteNumber(state.grid.couplingGain, 0.55) * bounded(alphaFit.estimate, 0.25, 4, 1), 0, 1, finiteNumber(state.grid.couplingGain, 0.55)),
      zCouplingGain: bounded(finiteNumber(state.grid.zCouplingGain, 0.32) * bounded(zFit.estimate, 0.25, 4, 1), 0, 1, finiteNumber(state.grid.zCouplingGain, 0.32)),
      backlash: Math.max(0, finiteNumber(state.grid.backlash, 0.1) + backlashBias),
      holeRadius: Math.max(pinRadius, pinRadius + proposedClearance),
      pinHoleClearance: proposedClearance,
      effectiveHoleRadiusBias: holeBias,
      effectiveBacklashBias: backlashBias,
      transitionBoundaryDirection: transitionShiftDirection,
      applyAutomatically: false,
      reason: "reduced-proxy transition estimate only; require holdout bench and segmented contact validation before physical claims",
    };
    const tolerance = Math.max(0, finiteNumber(options.tolerance, 1e-6));
    const phaseAccuracy = phaseObservationCount ? phaseMatchCount / phaseObservationCount : null;
    const rmsError = observedScalarCount ? Math.sqrt(squaredError / observedScalarCount) : null;
    const phaseMismatchCount = phaseObservationCount - phaseMatchCount;
    return {
      schema: TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA,
      model: "browser-two-cell-radius-backlash-transition-measured-vs-reduced-proxy-comparison",
      cadReference: CAD_RAD_CELL_REFERENCE,
      transitionReportSummary: report.summary,
      rows: comparisonRows,
      summary: {
        providedRowCount: measurements.length,
        matchedRowCount,
        missingPredictionCount,
        missingObservationRowCount,
        observedScalarCount,
        rmsError,
        maxAbsError: observedScalarCount ? maxAbsError : null,
        phaseObservationCount,
        phaseAccuracy,
        phaseMismatchCount,
        lockObservationCount,
        lockHeldAccuracy: lockObservationCount ? lockMatchCount / lockObservationCount : null,
        transitionShiftVoteCounts: shiftVotes,
        transitionShiftDirection,
        passesTolerance: matchedRowCount > 0 && observedScalarCount > 0 && rmsError <= tolerance && phaseMismatchCount === 0,
        physicalAccuracyValidated: false,
      },
      calibrationEstimate: {
        rightZScale: zFit,
        rightAlphaResponseScale: alphaFit,
        verticalSlipScale: verticalSlipFit,
        verticalExcessScale: verticalExcessFit,
        transitionBoundaryDirection: transitionShiftDirection,
        proposedReducedProxyUpdates,
        applyAutomatically: false,
        reason: "transition rows estimate reduced-proxy correction only; real geometry still needs segmented CAD/contact and holdout measurements",
      },
      acceptance: {
        readyForReducedProxyTransitionCalibration: observedScalarCount >= 4 && phaseObservationCount >= 2 && missingPredictionCount === 0,
        requiresMoreData: measurements.length === 0 || observedScalarCount === 0 || phaseObservationCount === 0,
        requiresSegmentedPhysics: true,
        requiresHoldoutValidation: true,
        physicalAccuracyValidated: false,
      },
      claimBoundary: {
        allowedClaim: "compares filled transition-adjacent bench rows against reduced-model boundary predictions",
        blockedClaim: "exact real RAD transition law without segmented CAD, contact parameters, and independent physical holdout sweeps",
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellRadiusBacklashTransitionComparisonJson(state, input, options = {}) {
    return JSON.stringify(compareTwoCellRadiusBacklashTransitionMeasurements(state, input, options), null, 2);
  }

  function exportTwoCellRadiusBacklashTransitionComparisonCsv(state, input, options = {}) {
    const comparison = compareTwoCellRadiusBacklashTransitionMeasurements(state, input, options);
    const columns = [
      "caseId",
      "phaseDiagramCaseId",
      "matrixIndex",
      "status",
      "category",
      "predictedPhase",
      "observedPhase",
      "phaseMatches",
      "transitionShiftVote",
      "backlash",
      "holeRadius",
      "effectiveBacklash",
      "effectiveHoleRadius",
      "effectivePinHoleClearanceMm",
      "predictedRightZ",
      "predictedRightAlpha",
      "predictedConnectorMaxVerticalSlipMm",
      "predictedConnectorMaxVerticalExcessMm",
      "observedScalarCount",
      "maxAbsResidual",
      "lockExpected",
      "lockHeldObserved",
      "lockHeldMatches",
      "notes",
    ];
    return [columns.join(","), ...comparison.rows.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function transitionMidpointValue(bracket) {
    if (String(bracket?.axis || "") === "holeRadius") return numericOrNull(bracket.midHoleRadius);
    if (String(bracket?.axis || "") === "backlash") return numericOrNull(bracket.midBacklash);
    return null;
  }

  function transitionFixedValue(bracket) {
    if (String(bracket?.axis || "") === "holeRadius") return numericOrNull(bracket.fixedBacklash);
    if (String(bracket?.axis || "") === "backlash") return numericOrNull(bracket.fixedHoleRadius);
    return null;
  }

  function transitionBoundaryMovementRows(baseline, calibrated) {
    const remaining = [...(calibrated || [])];
    const rows = [];
    for (let index = 0; index < (baseline || []).length; index += 1) {
      const before = baseline[index];
      const beforeAxis = String(before.axis || "");
      const beforePair = `${before.fromPhase || ""}->${before.toPhase || ""}`;
      const beforeMidpoint = transitionMidpointValue(before);
      const beforeFixed = transitionFixedValue(before);
      let bestIndex = -1;
      let bestScore = [Infinity, Infinity, Infinity];
      for (let candidateIndex = 0; candidateIndex < remaining.length; candidateIndex += 1) {
        const candidate = remaining[candidateIndex];
        const axisPenalty = String(candidate.axis || "") === beforeAxis ? 0 : 1;
        const pairPenalty = `${candidate.fromPhase || ""}->${candidate.toPhase || ""}` === beforePair ? 0 : 1;
        const candidateFixed = transitionFixedValue(candidate);
        const candidateMidpoint = transitionMidpointValue(candidate);
        const fixedDelta = beforeFixed !== null && candidateFixed !== null ? Math.abs(candidateFixed - beforeFixed) : 1e6;
        const midpointDelta = beforeMidpoint !== null && candidateMidpoint !== null ? Math.abs(candidateMidpoint - beforeMidpoint) : 1e6;
        const score = [axisPenalty + pairPenalty, fixedDelta, midpointDelta];
        if (
          score[0] < bestScore[0] ||
          (score[0] === bestScore[0] && score[1] < bestScore[1]) ||
          (score[0] === bestScore[0] && score[1] === bestScore[1] && score[2] < bestScore[2])
        ) {
          bestScore = score;
          bestIndex = candidateIndex;
        }
      }
      const after = bestIndex >= 0 ? remaining.splice(bestIndex, 1)[0] : null;
      const afterMidpoint = after ? transitionMidpointValue(after) : null;
      rows.push({
        index: index + 1,
        status: after ? "matched-nearest" : "missing-after-calibration",
        baselineTransitionId: before.transitionId || "",
        calibratedTransitionId: after?.transitionId || "",
        axis: beforeAxis,
        baselineFromPhase: before.fromPhase || "",
        baselineToPhase: before.toPhase || "",
        calibratedFromPhase: after?.fromPhase || "",
        calibratedToPhase: after?.toPhase || "",
        baselineMidpoint: beforeMidpoint,
        calibratedMidpoint: afterMidpoint,
        midpointDelta: beforeMidpoint !== null && afterMidpoint !== null ? afterMidpoint - beforeMidpoint : null,
        baselineFixedValue: beforeFixed,
        calibratedFixedValue: after ? transitionFixedValue(after) : null,
        physicalAccuracyValidated: false,
      });
    }
    const start = rows.length + 1;
    remaining.forEach((after, offset) => {
      rows.push({
        index: start + offset,
        status: "new-after-calibration",
        baselineTransitionId: "",
        calibratedTransitionId: after.transitionId || "",
        axis: after.axis || "",
        baselineFromPhase: "",
        baselineToPhase: "",
        calibratedFromPhase: after.fromPhase || "",
        calibratedToPhase: after.toPhase || "",
        baselineMidpoint: null,
        calibratedMidpoint: transitionMidpointValue(after),
        midpointDelta: null,
        baselineFixedValue: null,
        calibratedFixedValue: transitionFixedValue(after),
        physicalAccuracyValidated: false,
      });
    });
    return rows;
  }

  function twoCellRadiusBacklashTransitionRerun(state, input, options = {}) {
    const measurements = transitionMeasurementRows(input);
    const baselineComparison = compareTwoCellRadiusBacklashTransitionMeasurements(state, measurements, options);
    const updates = baselineComparison.calibrationEstimate?.proposedReducedProxyUpdates || {};
    const holeRadiusBias = finiteNumber(updates.effectiveHoleRadiusBias, 0);
    const backlashBias = finiteNumber(updates.effectiveBacklashBias, 0);
    const calibratedOptions = { ...options, holeRadiusBias, backlashBias };
    const baseline = twoCellRadiusBacklashTransitionReport(state, options);
    const calibrated = twoCellRadiusBacklashTransitionReport(state, calibratedOptions);
    const calibratedComparison = compareTwoCellRadiusBacklashTransitionMeasurements(state, measurements, calibratedOptions);
    const transitionMovement = transitionBoundaryMovementRows(baseline.transitionBrackets, calibrated.transitionBrackets);
    const finiteDeltas = transitionMovement
      .map((row) => numericOrNull(row.midpointDelta))
      .filter((value) => value !== null)
      .map((value) => Math.abs(value));
    return {
      schema: TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA,
      model: "browser-two-cell-radius-backlash-transition-calibrated-rerun-report",
      cadReference: CAD_RAD_CELL_REFERENCE,
      axisBiases: {
        effectiveHoleRadiusBias: holeRadiusBias,
        effectiveBacklashBias: backlashBias,
        transitionBoundaryDirection: updates.transitionBoundaryDirection || "",
      },
      proposedReducedProxyUpdates: updates,
      baseline,
      calibrated,
      baselineComparison,
      calibratedComparison,
      transitionMovement,
      summary: {
        status: "calibrated-transition-rerun-ready-needs-holdout-validation",
        providedMeasurementRowCount: measurements.length,
        baselineTransitionCount: finiteNumber(baseline.summary?.transitionBracketCount, 0),
        calibratedTransitionCount: finiteNumber(calibrated.summary?.transitionBracketCount, 0),
        transitionCountDelta: finiteNumber(calibrated.summary?.transitionBracketCount, 0) - finiteNumber(baseline.summary?.transitionBracketCount, 0),
        baselinePhaseAccuracy: baselineComparison.summary?.phaseAccuracy ?? null,
        calibratedPhaseAccuracy: calibratedComparison.summary?.phaseAccuracy ?? null,
        transitionMovementCount: transitionMovement.length,
        maxAbsMidpointDelta: Math.max(0, ...finiteDeltas),
        transitionShiftDirection: updates.transitionBoundaryDirection || "",
        physicalAccuracyValidated: false,
        remainingEvidence: [
          "independent holdout transition measurements after applying the proxy correction",
          "segmented CAD bodies with pin-hole contact surfaces",
          "external rigid-body contact results using the same nominal geometry",
          "bench repeatability across at least two physical assemblies",
        ],
      },
      claimBoundary: {
        allowedClaim: "shows how measured transition votes move the reduced-model effective boundary",
        blockedClaim: "validated physical transition law until holdout sweeps and segmented CAD/contact agree",
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellRadiusBacklashTransitionRerunJson(state, input, options = {}) {
    return JSON.stringify(twoCellRadiusBacklashTransitionRerun(state, input, options), null, 2);
  }

  function exportTwoCellRadiusBacklashTransitionRerunCsv(state, input, options = {}) {
    const report = twoCellRadiusBacklashTransitionRerun(state, input, options);
    const columns = [
      "index",
      "status",
      "baselineTransitionId",
      "calibratedTransitionId",
      "axis",
      "baselineFromPhase",
      "baselineToPhase",
      "calibratedFromPhase",
      "calibratedToPhase",
      "baselineMidpoint",
      "calibratedMidpoint",
      "midpointDelta",
      "baselineFixedValue",
      "calibratedFixedValue",
      "physicalAccuracyValidated",
    ];
    return [columns.join(","), ...report.transitionMovement.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function mean(values) {
    const finite = values.map((value) => Number(value)).filter(Number.isFinite);
    if (!finite.length) return null;
    return finite.reduce((total, value) => total + value, 0) / finite.length;
  }

  function summarizeAtlasRows(rows) {
    const lockRows = rows.filter((row) => row.lockMode === "right_state_locked" || row.lockMode === "right_position_locked");
    return {
      count: rows.length,
      lockCaseCount: lockRows.length,
      lockPassCount: lockRows.filter((row) => row.lockPass === true).length,
      stateLockVerticalMotionSamples: rows.filter(
        (row) => row.lockMode === "right_state_locked" && Math.abs(finiteNumber(row.rightZ, 0)) > 1e-9
      ).length,
      maxRightHeight: Math.max(0, ...rows.map((row) => Math.abs(finiteNumber(row.rightZ, 0)))),
      meanRightHeight: mean(rows.map((row) => row.rightZ)),
      maxVerticalSlipMm: Math.max(0, ...rows.map((row) => finiteNumber(row.connectorMaxVerticalSlipMm, 0))),
      meanVerticalSlipMm: mean(rows.map((row) => row.connectorMaxVerticalSlipMm)),
      maxContactPenalty: Math.max(0, ...rows.map((row) => finiteNumber(row.connectorTotalContactPenalty, 0))),
      meanContactPenalty: mean(rows.map((row) => row.connectorTotalContactPenalty)),
      maxTotalEnergy: Math.max(0, ...rows.map((row) => finiteNumber(row.totalEnergy, 0))),
    };
  }

  function groupAtlasRows(rows, field) {
    const groups = new Map();
    for (const row of rows) {
      const key = String(row[field]);
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(row);
    }
    return Array.from(groups.entries())
      .sort(([left], [right]) => left.localeCompare(right, undefined, { numeric: true }))
      .map(([value, groupRows]) => ({ field, value, ...summarizeAtlasRows(groupRows) }));
  }

  function atlasCaseDigest(row, category = "") {
    return {
      caseId: `${category ? `${category}:` : ""}${fidelityMatrixCaseId(row)}`,
      index: row.index,
      category,
      actuationCase: row.actuationCase,
      lockMode: row.lockMode,
      backlash: row.backlash,
      alphaCommand: row.alphaCommand,
      zCommand: row.zCommand,
      gravityForce: row.gravityForce,
      pinRadius: row.pinRadius,
      holeRadius: row.holeRadius,
      pinHoleClearanceMm: row.pinHoleClearanceMm,
      leftX: row.leftX,
      leftZ: row.leftZ,
      rightX: row.rightX,
      rightZ: row.rightZ,
      leftAlpha: row.leftAlpha,
      rightAlpha: row.rightAlpha,
      leftTheta: row.leftTheta,
      rightTheta: row.rightTheta,
      connectorMaxLateralSlipMm: row.connectorMaxLateralSlipMm,
      connectorMaxVerticalSlipMm: row.connectorMaxVerticalSlipMm,
      connectorMaxVerticalExcessMm: row.connectorMaxVerticalExcessMm,
      connectorTotalContactPenalty: row.connectorTotalContactPenalty,
      totalEnergy: row.totalEnergy,
      lockPass: row.lockPass,
      positionLocksHold: row.positionLocksHold,
      stateLocksHold: row.stateLocksHold,
      rightLockedVerticalMotionAllowed: row.rightLockedVerticalMotionAllowed,
      physicalAccuracyValidated: false,
    };
  }

  function rankedAtlasCases(rows, sortKey, category, limit, filter = () => true) {
    return rows
      .filter(filter)
      .slice()
      .sort((left, right) => finiteNumber(right[sortKey], 0) - finiteNumber(left[sortKey], 0))
      .slice(0, limit)
      .map((row) => atlasCaseDigest(row, category));
  }

  function twoCellPhysicalResponseAtlas(state, options = {}) {
    const matrix = twoCellPhysicalFidelityMatrix(state, options);
    const rows = matrix.rows || [];
    const connectorRows = matrix.connectorRows || [];
    const summary = matrix.summary || {};
    const rankLimit = Math.max(1, Math.round(finiteNumber(options.rankLimit, 12)));
    const lockRows = rows.filter((row) => row.lockMode === "right_state_locked" || row.lockMode === "right_position_locked");
    const lockFailures = lockRows.filter((row) => row.lockPass !== true);
    const stateLockRows = rows.filter((row) => row.lockMode === "right_state_locked");
    const positionLockRows = rows.filter((row) => row.lockMode === "right_position_locked");
    const contactModeCounts = connectorRows.reduce((counts, row) => {
      const mode = row.contactMode || "unknown";
      counts[mode] = (counts[mode] || 0) + 1;
      return counts;
    }, {});
    const rankedCases = {
      maxVerticalSlip: rankedAtlasCases(rows, "connectorMaxVerticalSlipMm", "maxVerticalSlip", rankLimit),
      maxContactPenalty: rankedAtlasCases(rows, "connectorTotalContactPenalty", "maxContactPenalty", rankLimit),
      maxRightHeight: rankedAtlasCases(rows, "rightZ", "maxRightHeight", rankLimit),
      lockFailures: lockFailures.slice(0, rankLimit).map((row) => atlasCaseDigest(row, "lockFailures")),
      stateLockVerticalMotion: rankedAtlasCases(
        stateLockRows,
        "rightZ",
        "stateLockVerticalMotion",
        rankLimit,
        (row) => Math.abs(finiteNumber(row.rightZ, 0)) > 1e-9
      ),
      positionLockFixtureCases: positionLockRows.slice(0, rankLimit).map((row) => atlasCaseDigest(row, "positionLockFixtureCases")),
    };
    const firstCasesToMeasure = [
      ...(rankedCases.lockFailures.length ? rankedCases.lockFailures.slice(0, 4) : rankedCases.maxContactPenalty.slice(0, 4)),
      ...rankedCases.maxVerticalSlip.slice(0, 4),
      ...rankedCases.positionLockFixtureCases.slice(0, 4),
    ];
    const uniqueFirstCases = Array.from(new Map(firstCasesToMeasure.map((item) => [item.caseId, item])).values()).slice(0, rankLimit);
    const lockPassCount = lockRows.filter((row) => row.lockPass === true).length;
    return {
      schema: TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA,
      model: "browser-two-cell-response-atlas-radius-backlash-lock-actuation-prioritizer",
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadLayout: matrix.cadLayout,
      axes: matrix.axes,
      matrixSummary: summary,
      groupSummaries: {
        byLockMode: groupAtlasRows(rows, "lockMode"),
        byHoleRadius: groupAtlasRows(rows, "holeRadius"),
        byBacklash: groupAtlasRows(rows, "backlash"),
        byActuationCase: groupAtlasRows(rows, "actuationCase"),
      },
      rankedCases,
      invariants: {
        solverAllConverged: summary.solverSuccessCount === summary.rowCount,
        lockAllPass: lockRows.length === 0 || lockPassCount === lockRows.length,
        stateLockAlphaHolds: summary.lockChecks?.stateLockAlphaHold === true,
        stateLockAllowsVerticalSamples: finiteNumber(summary.lockChecks?.stateLockVerticalMotionSamples, 0),
        positionLockFixtureHolds: summary.lockChecks?.positionLockFixtureHold === true,
        clearanceNondecreasing: summary.clearanceTrend?.clearanceNondecreasing === true,
        clearancePenaltyNonincreasing: summary.clearanceTrend?.penaltyNonincreasing === true,
        backlashAlphaResponseNonincreasing: summary.backlashTrend?.alphaResponseNonincreasing === true,
        zLiftPositive: summary.actuationPolarity?.positiveZLifts === true,
        zPushdownNegative: summary.actuationPolarity?.negativeZPushesDown === true,
        alphaContractBelowInitial: summary.actuationPolarity?.contractDecreasesAlpha === true,
        alphaExpandAboveInitial: summary.actuationPolarity?.expandIncreasesAlpha === true,
        physicalAccuracyValidated: false,
      },
      connectorSummary: {
        connectorRowCount: connectorRows.length,
        contactModeCounts,
        maxLateralSlipMm: Math.max(0, ...rows.map((row) => finiteNumber(row.connectorMaxLateralSlipMm, 0))),
        maxVerticalSlipMm: Math.max(0, ...rows.map((row) => finiteNumber(row.connectorMaxVerticalSlipMm, 0))),
        maxVerticalExcessMm: Math.max(0, ...rows.map((row) => finiteNumber(row.connectorMaxVerticalExcessMm, 0))),
        maxContactPenalty: Math.max(0, ...rows.map((row) => finiteNumber(row.connectorTotalContactPenalty, 0))),
      },
      benchPriority: {
        firstCasesToMeasure: uniqueFirstCases,
        lockCasesToMeasure: rankedCases.lockFailures.length ? rankedCases.lockFailures : rankedCases.stateLockVerticalMotion,
        clearanceCasesToMeasure: rankedCases.maxVerticalSlip,
      },
      claimBoundary: {
        physicalAccuracyValidated: false,
        remainingEvidence: summary.missingEvidence || [],
        note: "Atlas ranks reduced two-cell proxy cases for bench/external-engine measurement; it is not a fabricated-cell validation by itself.",
      },
    };
  }

  function exportTwoCellPhysicalResponseAtlasJson(state, options = {}) {
    return JSON.stringify(twoCellPhysicalResponseAtlas(state, options), null, 2);
  }

  function exportTwoCellPhysicalResponseAtlasCsv(state, options = {}) {
    const atlas = twoCellPhysicalResponseAtlas(state, options);
    const rows = [];
    for (const [category, cases] of Object.entries(atlas.rankedCases || {})) {
      for (const item of cases || []) rows.push({ category, ...item });
    }
    const columns = [
      "category",
      "caseId",
      "index",
      "actuationCase",
      "lockMode",
      "backlash",
      "alphaCommand",
      "zCommand",
      "pinRadius",
      "holeRadius",
      "pinHoleClearanceMm",
      "rightZ",
      "rightAlpha",
      "rightTheta",
      "connectorMaxVerticalSlipMm",
      "connectorMaxVerticalExcessMm",
      "connectorTotalContactPenalty",
      "totalEnergy",
      "lockPass",
      "positionLocksHold",
      "stateLocksHold",
      "physicalAccuracyValidated",
    ];
    return [columns.join(","), ...rows.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function fidelityMatrixCaseId(row) {
    return `FM_${String(Math.round(finiteNumber(row.index, 0))).padStart(3, "0")}`;
  }

  function twoCellFidelityMatrixMeasurementTemplate(state, options = {}) {
    const matrix = twoCellPhysicalFidelityMatrix(state, options);
    const byIndex = new Map(matrix.rows.map((row) => [row.index, row]));
    const rows = (matrix.connectorRows || []).map((connector) => {
      const parent = byIndex.get(connector.index) || {};
      return {
        caseId: fidelityMatrixCaseId(parent),
        matrixIndex: connector.index,
        connector: connector.connector,
        leftSite: connector.leftSite,
        rightSite: connector.rightSite,
        actuationCase: parent.actuationCase,
        lockMode: parent.lockMode,
        backlash: parent.backlash,
        alphaCommand: parent.alphaCommand,
        zCommand: parent.zCommand,
        gravityForce: parent.gravityForce,
        pinRadius: parent.pinRadius,
        holeRadius: parent.holeRadius,
        pinHoleClearanceMm: parent.pinHoleClearanceMm,
        predictedLeftCellX: parent.leftX,
        predictedLeftCellY: 0,
        predictedLeftCellZ: parent.leftZ,
        predictedRightCellX: parent.rightX,
        predictedRightCellY: 0,
        predictedRightCellZ: parent.rightZ,
        predictedLeftAlpha: parent.leftAlpha,
        predictedRightAlpha: parent.rightAlpha,
        predictedLeftTheta: parent.leftTheta,
        predictedRightTheta: parent.rightTheta,
        predictedLeftXmm: connector.leftPositionMm?.x ?? "",
        predictedLeftYmm: connector.leftPositionMm?.y ?? "",
        predictedLeftZmm: connector.leftPositionMm?.z ?? "",
        predictedRightXmm: connector.rightPositionMm?.x ?? "",
        predictedRightYmm: connector.rightPositionMm?.y ?? "",
        predictedRightZmm: connector.rightPositionMm?.z ?? "",
        predictedSlipXmm: connector.slipMm?.x ?? "",
        predictedSlipYmm: connector.slipMm?.y ?? "",
        predictedSlipZmm: connector.slipMm?.z ?? "",
        predictedLateralSlipMm: connector.lateralSlipMm,
        predictedVerticalSlipMm: connector.verticalSlipMm,
        predictedTotalSlipMm: connector.totalSlipMm,
        predictedLateralExcessMm: connector.lateralExcessMm,
        predictedVerticalExcessMm: connector.verticalExcessMm,
        predictedContactMode: connector.contactMode,
        predictedLockPass: parent.lockPass,
        observedLeftCellX: "",
        observedLeftCellY: "",
        observedLeftCellZ: "",
        observedRightCellX: "",
        observedRightCellY: "",
        observedRightCellZ: "",
        observedLeftAlpha: "",
        observedRightAlpha: "",
        observedLeftTheta: "",
        observedRightTheta: "",
        observedLeftXmm: "",
        observedLeftYmm: "",
        observedLeftZmm: "",
        observedRightXmm: "",
        observedRightYmm: "",
        observedRightZmm: "",
        observedLateralSlipMm: "",
        observedVerticalSlipMm: "",
        observedTotalSlipMm: "",
        observedContactMode: "",
        lockHeldObserved: "",
        measurementSource: "",
        notes: "",
      };
    });
    return {
      schema: TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA,
      model: "browser-two-cell-fidelity-matrix-connector-and-cell-fillable-template",
      cadReference: CAD_RAD_CELL_REFERENCE,
      matrixSummary: matrix.summary,
      rows,
      summary: {
        matrixRowCount: matrix.summary.rowCount,
        connectorRowCount: rows.length,
        connectorNames: [...new Set(rows.map((row) => row.connector))].sort(),
        measurementRequired: true,
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellFidelityMatrixMeasurementTemplateJson(state, options = {}) {
    return JSON.stringify(twoCellFidelityMatrixMeasurementTemplate(state, options), null, 2);
  }

  function exportTwoCellFidelityMatrixMeasurementTemplateCsv(state, options = {}) {
    const template = twoCellFidelityMatrixMeasurementTemplate(state, options);
    const columns = template.rows.length ? Object.keys(template.rows[0]) : [];
    return [columns.join(","), ...template.rows.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function manifestLockFlags(lockMode) {
    return {
      leftPositionLocked: lockMode !== "left_free",
      rightStateLocked: lockMode === "right_state_locked",
      rightPositionLocked: lockMode === "right_position_locked",
    };
  }

  function twoCellExternalFidelityMatrixManifestRows(state, options = {}) {
    const matrix = twoCellPhysicalFidelityMatrix(state, options);
    return (matrix.rows || []).map((row) => {
      const caseId = fidelityMatrixCaseId(row);
      return {
        caseId,
        matrixIndex: row.index,
        actuationCase: row.actuationCase,
        lockMode: row.lockMode,
        backlash: row.backlash,
        alphaCommand: row.alphaCommand,
        zCommand: row.zCommand,
        gravityForce: row.gravityForce,
        pinRadius: row.pinRadius,
        holeRadius: row.holeRadius,
        pinHoleClearanceMm: row.pinHoleClearanceMm,
        ...manifestLockFlags(row.lockMode),
        expectedLeftX: row.leftX,
        expectedLeftY: 0,
        expectedLeftZ: row.leftZ,
        expectedRightX: row.rightX,
        expectedRightY: 0,
        expectedRightZ: row.rightZ,
        expectedLeftAlpha: row.leftAlpha,
        expectedRightAlpha: row.rightAlpha,
        expectedLeftTheta: row.leftTheta,
        expectedRightTheta: row.rightTheta,
        expectedContactMode: row.contactMode,
        expectedConnectorMaxLateralSlipMm: row.connectorMaxLateralSlipMm,
        expectedConnectorMaxVerticalSlipMm: row.connectorMaxVerticalSlipMm,
        expectedConnectorMaxVerticalExcessMm: row.connectorMaxVerticalExcessMm,
        expectedConnectorTotalContactPenalty: row.connectorTotalContactPenalty,
        mjcfProxyPath: `external_fidelity_matrix_mjcf/${caseId}.xml`,
        engineResultPath: `external_fidelity_matrix_results/${caseId}.json`,
        measurementRows: 3,
        runStatus: "pending-external-engine",
      };
    });
  }

  function twoCellExternalFidelityMatrixManifest(state, options = {}) {
    const matrix = twoCellPhysicalFidelityMatrix(state, options);
    const template = twoCellFidelityMatrixMeasurementTemplate(state, options);
    const caseRows = twoCellExternalFidelityMatrixManifestRows(state, options);
    const audit = cadRadCellArchiveAudit();
    const missingEvidence = Array.from(
      new Set([
        ...(audit.summary?.missingEvidence || []),
        "segmentedCADBodies",
        "filledSegmentedCadIntake",
        "externalEngineRunResults",
        "filledFidelityMatrixMeasurements",
        "benchCoordinateHoldout",
      ])
    ).sort();
    return {
      schema: TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA,
      model: "browser-dense-two-cell-external-engine-fidelity-matrix-handoff",
      engine: {
        requested: String(options.engine || "mujoco").toLowerCase(),
        proxyEngineReady: true,
        exactSegmentedCadReady: false,
        runScope: "252 two-cell cases with three connector measurement rows per case",
      },
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadArchiveAuditSummary: audit.summary,
      axes: matrix.axes,
      caseRows,
      connectorMeasurementRows: template.rows,
      resultSchema: {
        caseResultJson: {
          caseId: "FM_000",
          sourceEngine: "mujoco",
          finalBodies: [
            { bodyId: "left", x: "float model units", y: "float model units", z: "float model units", alpha: "float", theta: "float degrees" },
            { bodyId: "right", x: "float model units", y: "float model units", z: "float model units", alpha: "float", theta: "float degrees" },
          ],
          connectors: [
            {
              connector: "upper|middle|lower",
              observedLeftXmm: "float",
              observedLeftYmm: "float",
              observedLeftZmm: "float",
              observedRightXmm: "float",
              observedRightYmm: "float",
              observedRightZmm: "float",
              observedLateralSlipMm: "float",
              observedVerticalSlipMm: "float",
              observedTotalSlipMm: "float",
              observedContactMode: "inside_clearance|clearance_edge|contacting",
              lockHeldObserved: "held|released|failed",
            },
          ],
        },
        fillableCsv: "two_cell_fidelity_matrix_measurement_template.csv",
      },
      runInstructions: [
        "Generate one external-engine case per caseId using the listed controls, locks, backlash, and hole radius.",
        "Run each case to static equilibrium under gravityForce and actuator commands.",
        "Record final body centers, alpha, theta, connector markers, slip, contact mode, and lock-held status.",
        "Fill the dense matrix template and load it back into the browser or Python comparison tool.",
      ],
      summary: {
        status: "external-fidelity-manifest-ready-needs-engine-results",
        caseCount: caseRows.length,
        connectorMeasurementRowCount: template.rows.length,
        matrixRowCount: matrix.summary?.rowCount || 0,
        proxyEngineReady: true,
        exactSegmentedCadReady: false,
        physicalAccuracyValidated: false,
        missingEvidence,
        missingEvidenceCount: missingEvidence.length,
      },
      claimLabels: {
        handoff: "external engine run manifest aligned to the dense two-cell fidelity matrix",
        geometry: "proxy-ready now; exact geometry requires segmented CAD intake",
        physicalAccuracy: "false until external solver and bench holdout measurements agree",
      },
    };
  }

  function exportTwoCellExternalFidelityMatrixManifestJson(state, options = {}) {
    return JSON.stringify(twoCellExternalFidelityMatrixManifest(state, options), null, 2);
  }

  function exportTwoCellExternalFidelityMatrixManifestCsv(state, options = {}) {
    const rows = twoCellExternalFidelityMatrixManifestRows(state, options);
    const columns = rows.length ? Object.keys(rows[0]) : [];
    return [columns.join(","), ...rows.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function handoffCaseKey(row, source, rank) {
    return String(row.caseId || row.phaseDiagramCaseId || row.matrixCaseId || row.index || `${source}-${rank}`);
  }

  function handoffLockFlags(row) {
    const lockMode = String(row.lockMode || "free");
    return {
      leftPositionLocked: row.leftPositionLocked ?? lockMode !== "left_free",
      rightStateLocked: row.rightStateLocked ?? lockMode === "right_state_locked",
      rightPositionLocked: row.rightPositionLocked ?? lockMode === "right_position_locked",
    };
  }

  function exactContactHandoffCaseRow(row, source, category, rank) {
    return {
      rank,
      caseId: handoffCaseKey(row, source, rank),
      source,
      category: String(row.category || category),
      phase: String(row.phase || ""),
      actuationCase: String(row.actuationCase || ""),
      lockMode: String(row.lockMode || "free"),
      backlash: row.backlash ?? "",
      holeRadius: row.holeRadius ?? "",
      pinHoleClearanceMm: row.pinHoleClearanceMm ?? "",
      effectiveBacklash: row.effectiveBacklash ?? row.backlash ?? "",
      effectiveHoleRadius: row.effectiveHoleRadius ?? row.holeRadius ?? "",
      effectivePinHoleClearanceMm: row.effectivePinHoleClearanceMm ?? row.pinHoleClearanceMm ?? "",
      alphaCommand: row.alphaCommand ?? "",
      zCommand: row.zCommand ?? "",
      gravityForce: row.gravityForce ?? "",
      leftZ: row.leftZ ?? "",
      rightZ: row.rightZ ?? "",
      rightAlpha: row.rightAlpha ?? "",
      connectorActiveContactCount: row.connectorActiveContactCount ?? "",
      connectorMaxVerticalSlipMm: row.connectorMaxVerticalSlipMm ?? "",
      connectorMaxVerticalExcessMm: row.connectorMaxVerticalExcessMm ?? "",
      connectorTotalContactPenalty: row.connectorTotalContactPenalty ?? "",
      ...handoffLockFlags(row),
      expectedConnectorCount: 3,
      requiredCadAssets: [...EXACT_CONTACT_HANDOFF_REQUIRED_CAD_ASSETS],
      requiredMeasuredInputs: [...EXACT_CONTACT_HANDOFF_REQUIRED_MEASURED_INPUTS],
      engineHandoffReady: false,
      canRunExactContact: false,
      physicalAccuracyValidated: false,
      notes:
        "Run with segmented rigid bodies, explicit pin-hole contacts, gravity, locks, and actuator commands; compare final cell centers and upper/middle/lower connector markers against bench data.",
    };
  }

  function shiftedCadSiteCenter(site, originMm) {
    return [
      finiteNumber(originMm[0], 0) + finiteNumber(site.xMm, 0),
      finiteNumber(originMm[1], 0) + finiteNumber(site.yMm, 0),
      finiteNumber(originMm[2], 0),
    ];
  }

  function twoCellCadContactDecompositionSpec(state, options = {}) {
    const layout = cadRadCellLayout(state, options);
    const dimensions = layout.dimensionsMm;
    const ringSegments = Math.max(4, Math.round(finiteNumber(options.ringSegments, 12)));
    const instances = {
      left: { cell: "left", originMm: [0, 0, 0], upperBody: "left_upper_free_cell_body", lowerBody: "left_lower_cell_body", pinBody: "left_screw_pin_body" },
      right: { cell: "right", originMm: [finiteNumber(dimensions.connectorPitchMm, 0), 0, 0], upperBody: "right_upper_free_cell_body", lowerBody: "right_lower_cell_body", pinBody: "right_screw_pin_body" },
    };
    const segmentStep = 360 / ringSegments;
    const holeDecomposition = [];
    for (const instance of Object.values(instances)) {
      for (const site of layout.padSites) {
        const centerMm = shiftedCadSiteCenter(site, instance.originMm);
        holeDecomposition.push({
          cell: instance.cell,
          site: site.name,
          centerMm,
          axis: [0, 0, 1],
          holeWallBody: instance.upperBody,
          pinBody: instance.pinBody,
          holeRadiusMm: site.holeRadiusMm,
          pinRadiusMm: site.pinRadiusMm,
          clearanceMm: dimensions.pinHoleClearanceMm,
          ringSegmentCount: ringSegments,
          recommendedCollisionPrimitive: "convex_ring_sector_or_capsule_wall",
          ringSegments: Array.from({ length: ringSegments }, (_, index) => ({
            segmentIndex: index,
            angleStartDeg: index * segmentStep,
            angleEndDeg: (index + 1) * segmentStep,
            centerMm,
            axis: [0, 0, 1],
            bodyA: instance.pinBody,
            bodyB: instance.upperBody,
            primitive: "convex_ring_sector_or_capsule_wall",
            source: "layout-derived-placeholder-needs-CAD-surface",
            physicalAccuracyValidated: false,
          })),
          source: "A360 layout-derived site; replace with segmented radial_pad_hole_surfaces evidence before exact claims",
          physicalAccuracyValidated: false,
        });
      }
    }
    const siteByName = Object.fromEntries(layout.padSites.map((site) => [site.name, site]));
    const markerZ = {
      lower: 0.5 * finiteNumber(dimensions.lowerBodyThicknessMm, 0),
      middle: finiteNumber(dimensions.lowerBodyThicknessMm, 0),
      upper: finiteNumber(dimensions.bodyThicknessMm, 0),
    };
    const twoCellContactPairs = CAD_RAD_X_CONNECTOR_PAIRS.map(([name, leftSite, rightSite]) => ({
      name,
      leftCell: "left",
      rightCell: "right",
      leftSite,
      rightSite,
      leftCenterMm: shiftedCadSiteCenter(siteByName[leftSite], instances.left.originMm),
      rightCenterMm: shiftedCadSiteCenter(siteByName[rightSite], instances.right.originMm),
      axis: [1, 0, 0],
      leftPinBody: `left_${leftSite}_connector_pin_or_screw`,
      rightHoleBody: instances.right.upperBody,
      holeRadiusMm: dimensions.holeRadiusMm,
      pinRadiusMm: dimensions.pinRadiusMm,
      clearanceMm: dimensions.pinHoleClearanceMm,
      requiredSurfacePair: "two_cell_connector_pairs.json plus radial_pad_hole_surfaces.json",
      observableMarkers: Object.entries(markerZ).map(([marker, zMm]) => ({ marker: `${name}_${marker}`, zMm })),
      source: "layout-derived X-neighbor connector pair; replace with measured connector origins before exact claims",
      physicalAccuracyValidated: false,
    }));
    const missingEvidence = [
      "segmentedCADBodies",
      "filledSegmentedCadIntake",
      "exactJointAxes",
      "exactPinHoleContactSurfaces",
      "convex-or-analytic-hole-wall-contact-primitives",
      "body-collision-exclusion-mask",
      "joint-limit-stop-surfaces",
      "external-engine-contact-run",
      "bench-coordinate-holdout-comparison",
    ].sort();
    return {
      schema: TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA,
      model: "browser-two-cell-cad-contact-decomposition-contract",
      cadReference: CAD_RAD_CELL_REFERENCE,
      layout,
      cellInstances: Object.values(instances),
      bodyRoles: [
        { role: "upper_free_cell_body", instances: [instances.left.upperBody, instances.right.upperBody], expectedAssets: ["upper_free_cell_body.step", "upper_free_cell_body.stl", "upper_free_cell_body.obj"] },
        { role: "lower_cell_body", instances: [instances.left.lowerBody, instances.right.lowerBody], expectedAssets: ["lower_cell_body.step", "lower_cell_body.stl", "lower_cell_body.obj"] },
        { role: "screw_pin_body", instances: [instances.left.pinBody, instances.right.pinBody], expectedAssets: ["screw_pin_body.step", "screw_pin_body.stl", "screw_pin_body.obj"] },
      ],
      jointSpec: {
        alphaRotation: { axis: [0, 0, 1], stateVariable: "alpha", thetaMapping: "theta_degrees = 70 * alpha - 60 in the reduced model; replace by measured crown/arm relation when available" },
        verticalSlide: { axis: [0, 0, 1], stateVariable: "z", contactRole: "pin-hole radius difference permits vertical residual motion under actuation and gravity" },
        lockStops: { note: "Crown angle is not equal to lattice theta because the arm has nonzero width; use measured effective stop theta." },
      },
      holeDecomposition,
      twoCellContactPairs,
      collisionPrimitiveContract: {
        holeWallPrimitive: "split each annular hole wall into convex sectors or analytic capsule/cylinder contacts",
        avoid: "single concave hole mesh as the active collision body",
        ringSegmentCount: ringSegments,
        minimumRecommendedRingSegments: 8,
        collisionPrimitiveCount: holeDecomposition.length * ringSegments + twoCellContactPairs.length,
      },
      summary: {
        status: "needs-cad-contact-decomposition-inputs",
        cellCount: 2,
        holeCount: holeDecomposition.length,
        twoCellConnectorCount: twoCellContactPairs.length,
        ringSegments,
        collisionPrimitiveCount: holeDecomposition.length * ringSegments + twoCellContactPairs.length,
        canBuildExactContactModel: false,
        physicalAccuracyValidated: false,
        missingEvidence,
        missingEvidenceCount: missingEvidence.length,
      },
      claimBoundary: {
        allowedClaim: "solver-ready decomposition contract derived from the A360 one-cell reference and current two-cell topology",
        blockedClaim: "fabrication-accurate contact mechanics until segmented CAD primitives and bench holdouts validate the model",
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellCadContactDecompositionJson(state, options = {}) {
    return JSON.stringify(twoCellCadContactDecompositionSpec(state, options), null, 2);
  }

  function exportTwoCellCadContactDecompositionCsv(state, options = {}) {
    const spec = twoCellCadContactDecompositionSpec(state, options);
    const columns = ["section", "cell", "connector", "site", "segmentIndex", "angleStartDeg", "angleEndDeg", "bodyA", "bodyB", "primitive", "centerXmm", "centerYmm", "centerZmm", "axisX", "axisY", "axisZ", "pinRadiusMm", "holeRadiusMm", "clearanceMm", "requiredAsset", "source", "physicalAccuracyValidated"];
    const rows = [];
    for (const hole of spec.holeDecomposition) {
      for (const segment of hole.ringSegments) {
        rows.push({
          section: "holeSegment",
          cell: hole.cell,
          site: hole.site,
          segmentIndex: segment.segmentIndex,
          angleStartDeg: segment.angleStartDeg,
          angleEndDeg: segment.angleEndDeg,
          bodyA: segment.bodyA,
          bodyB: segment.bodyB,
          primitive: segment.primitive,
          centerXmm: segment.centerMm[0],
          centerYmm: segment.centerMm[1],
          centerZmm: segment.centerMm[2],
          axisX: segment.axis[0],
          axisY: segment.axis[1],
          axisZ: segment.axis[2],
          pinRadiusMm: hole.pinRadiusMm,
          holeRadiusMm: hole.holeRadiusMm,
          clearanceMm: hole.clearanceMm,
          requiredAsset: "radial_pad_hole_surfaces.json",
          source: segment.source,
          physicalAccuracyValidated: false,
        });
      }
    }
    for (const pair of spec.twoCellContactPairs) {
      rows.push({
        section: "connectorPair",
        cell: `${pair.leftCell}-${pair.rightCell}`,
        connector: pair.name,
        site: `${pair.leftSite}-${pair.rightSite}`,
        bodyA: pair.leftPinBody,
        bodyB: pair.rightHoleBody,
        primitive: "explicit_pin_hole_connector_contact",
        centerXmm: pair.leftCenterMm[0],
        centerYmm: pair.leftCenterMm[1],
        centerZmm: pair.leftCenterMm[2],
        axisX: pair.axis[0],
        axisY: pair.axis[1],
        axisZ: pair.axis[2],
        pinRadiusMm: pair.pinRadiusMm,
        holeRadiusMm: pair.holeRadiusMm,
        clearanceMm: pair.clearanceMm,
        requiredAsset: pair.requiredSurfacePair,
        source: pair.source,
        physicalAccuracyValidated: false,
      });
    }
    return [columns.join(","), ...rows.map((row) => columns.map((column) => csvValue(row[column])).join(","))].join("\n");
  }

  function twoCellExactContactHandoffPlan(state, options = {}) {
    const transitionReport = twoCellRadiusBacklashTransitionReport(state, options);
    const atlas = twoCellPhysicalResponseAtlas(state, options);
    const manifest = twoCellExternalFidelityMatrixManifest(state, options);
    const audit = cadRadCellArchiveAudit();
    const contactDecomposition = twoCellCadContactDecompositionSpec(state, options);
    const missingEvidence = Array.from(
      new Set([
        ...(audit.summary?.missingEvidence || []),
        ...(manifest.summary?.missingEvidence || []),
        ...(contactDecomposition.summary?.missingEvidence || []),
        "segmentedCADBodies",
        "filledSegmentedCadIntake",
        "exactJointAxes",
        "exactPinHoleContactSurfaces",
        "externalRigidBodyContactResults",
        "benchCoordinateHoldout",
      ])
    ).sort();
    const caseRows = [];
    const byCase = new Map();
    const addCases = (source, category, rows) => {
      for (const input of rows || []) {
        const key = handoffCaseKey(input, source, caseRows.length + 1);
        const existing = byCase.get(key);
        if (existing) {
          existing.source = Array.from(new Set(String(existing.source).split(";").concat(source))).sort().join(";");
          existing.category = Array.from(new Set(String(existing.category).split(";").concat(String(input.category || category)))).sort().join(";");
          continue;
        }
        const out = exactContactHandoffCaseRow(input, source, category, caseRows.length + 1);
        byCase.set(key, out);
        caseRows.push(out);
      }
    };
    addCases("transitionReport", "transition-boundary", transitionReport.measurementCases);
    addCases("responseAtlas.first", "high-response", atlas.benchPriority?.firstCasesToMeasure);
    addCases("responseAtlas.lock", "lock-diagnostic", atlas.benchPriority?.lockCasesToMeasure);
    addCases("responseAtlas.clearance", "clearance-dieoff", atlas.benchPriority?.clearanceCasesToMeasure);
    const limit = Math.max(1, Math.round(finiteNumber(options.priorityLimit, 24)));
    const limitedRows = caseRows.slice(0, limit).map((row, index) => ({ ...row, rank: index + 1 }));
    return {
      schema: TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA,
      model: "browser-two-cell-exact-rigid-body-contact-handoff-plan",
      cadReference: CAD_RAD_CELL_REFERENCE,
      cadArchiveAuditSummary: audit.summary,
      externalFidelityManifest: manifest.summary,
      cadContactDecomposition: contactDecomposition.summary,
      engineTargets: [
        { engine: "mujoco", status: "browser-file-status-unknown", expectedAssets: ["mujoco/rad_two_cell.xml", "rad_two_cell.xml"], detectedAssets: [] },
        { engine: "gazebo", status: "browser-file-status-unknown", expectedAssets: ["gazebo/rad_two_cell.sdf", "rad_two_cell.sdf", "rad_two_cell.urdf"], detectedAssets: [] },
        { engine: "isaac-sim", status: "browser-file-status-unknown", expectedAssets: ["isaac/rad_two_cell.usd", "rad_two_cell.usd"], detectedAssets: [] },
      ],
      caseRows: limitedRows,
      protocol: {
        startFrom: "one-cell A360 CAD reference plus two connected-cell X-neighbor assembly",
        requiredRigidBodies: [
          "left upper moving cell body",
          "left lower fixture body",
          "right upper moving cell body",
          "right lower fixture body",
          "connector pins or screw bodies",
        ],
        requiredContacts: [
          "pin-hole radial clearance contacts",
          "vertical slide/free-play limits",
          "lock crown contact/stop angles",
          "ground/support contacts under gravity",
        ],
        collisionPrimitiveContract: contactDecomposition.collisionPrimitiveContract,
        observables: [
          "left/right cell center xyz",
          "left/right alpha and theta",
          "upper/middle/lower connector marker xyz",
          "contact mode",
          "lock held/released/failed",
        ],
      },
      summary: {
        status: "needs-segmented-cad-contact-inputs",
        caseCount: limitedRows.length,
        transitionCaseCount: limitedRows.filter((row) => String(row.source).includes("transitionReport")).length,
        atlasCaseCount: limitedRows.filter((row) => String(row.source).includes("responseAtlas")).length,
        engineHandoffReady: false,
        intakeValidationReady: false,
        canAttemptExactRigidBodyContact: false,
        canRunExactContact: false,
        physicalAccuracyValidated: false,
        missingEvidence,
        missingEvidenceCount: missingEvidence.length,
      },
      claimBoundary: {
        allowedClaim: "ranked exact-contact and bench handoff plan derived from current reduced two-cell diagnostics",
        blockedClaim: "exact real RAD mechanics until segmented CAD, external contact results, and bench holdouts pass",
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellExactContactHandoffPlanJson(state, options = {}) {
    return JSON.stringify(twoCellExactContactHandoffPlan(state, options), null, 2);
  }

  function exportTwoCellExactContactHandoffPlanCsv(state, options = {}) {
    const plan = twoCellExactContactHandoffPlan(state, options);
    const columns = [
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
    ];
    return [
      columns.join(","),
      ...plan.caseRows.map((row) => columns.map((column) => csvValue(row[column])).join(",")),
    ].join("\n");
  }

  function twoCellFidelityMatrixMeasurementsFromCsv(text) {
    return {
      schema: TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA,
      rows: parseCsv(text),
    };
  }

  function measurementRows(input) {
    if (Array.isArray(input)) return input;
    if (input && Array.isArray(input.rows)) return input.rows;
    if (typeof input === "string") {
      const trimmed = input.trim();
      if (!trimmed) return [];
      if (trimmed.startsWith("{") || trimmed.startsWith("[")) {
        const parsed = JSON.parse(trimmed);
        return measurementRows(parsed);
      }
      return twoCellFidelityMatrixMeasurementsFromCsv(trimmed).rows;
    }
    return [];
  }

  function compareTwoCellFidelityMatrixMeasurements(state, input, options = {}) {
    const rows = measurementRows(input);
    const template = twoCellFidelityMatrixMeasurementTemplate(state, options);
    const expected = new Map(template.rows.map((row) => [`${row.caseId}|${row.connector}`, row]));
    const scalarPairs = [
      ["observedLeftCellX", "predictedLeftCellX", "cell"],
      ["observedLeftCellY", "predictedLeftCellY", "cell"],
      ["observedLeftCellZ", "predictedLeftCellZ", "cell"],
      ["observedRightCellX", "predictedRightCellX", "cell"],
      ["observedRightCellY", "predictedRightCellY", "cell"],
      ["observedRightCellZ", "predictedRightCellZ", "cell"],
      ["observedLeftAlpha", "predictedLeftAlpha", "cell"],
      ["observedRightAlpha", "predictedRightAlpha", "cell"],
      ["observedLeftTheta", "predictedLeftTheta", "cell"],
      ["observedRightTheta", "predictedRightTheta", "cell"],
      ["observedLeftXmm", "predictedLeftXmm", "connector"],
      ["observedLeftYmm", "predictedLeftYmm", "connector"],
      ["observedLeftZmm", "predictedLeftZmm", "connector"],
      ["observedRightXmm", "predictedRightXmm", "connector"],
      ["observedRightYmm", "predictedRightYmm", "connector"],
      ["observedRightZmm", "predictedRightZmm", "connector"],
      ["observedLateralSlipMm", "predictedLateralSlipMm", "connector"],
      ["observedVerticalSlipMm", "predictedVerticalSlipMm", "connector"],
      ["observedTotalSlipMm", "predictedTotalSlipMm", "connector"],
    ];
    const tolerance = Math.max(0, finiteNumber(options.tolerance, 1e-6));
    const comparisonRows = [];
    let matchedRowCount = 0;
    let unmatchedRowCount = 0;
    let observedScalarCount = 0;
    let missingObservationCount = 0;
    let squaredError = 0;
    let maxAbsError = 0;
    let cellSquaredError = 0;
    let cellCount = 0;
    let connectorSquaredError = 0;
    let connectorCount = 0;
    let contactModeCount = 0;
    let contactModeMatchCount = 0;
    let lockHeldCount = 0;
    let lockHeldMatchCount = 0;
    for (const row of rows) {
      const key = `${row.caseId}|${row.connector}`;
      const predicted = expected.get(key);
      if (!predicted) {
        unmatchedRowCount += 1;
        comparisonRows.push({ caseId: row.caseId || "", connector: row.connector || "", matched: false, maxAbsError: "", observedScalarCount: 0 });
        continue;
      }
      matchedRowCount += 1;
      let rowMax = 0;
      let rowScalarCount = 0;
      for (const [observedField, predictedField, group] of scalarPairs) {
        const observed = numericOrNull(row[observedField]);
        const target = numericOrNull(predicted[predictedField] ?? row[predictedField]);
        if (observed === null || target === null) {
          missingObservationCount += 1;
          continue;
        }
        const delta = observed - target;
        const absDelta = Math.abs(delta);
        squaredError += delta * delta;
        maxAbsError = Math.max(maxAbsError, absDelta);
        rowMax = Math.max(rowMax, absDelta);
        observedScalarCount += 1;
        rowScalarCount += 1;
        if (group === "cell") {
          cellSquaredError += delta * delta;
          cellCount += 1;
        } else {
          connectorSquaredError += delta * delta;
          connectorCount += 1;
        }
      }
      const observedMode = String(row.observedContactMode || "").trim();
      const predictedMode = String(predicted.predictedContactMode || row.predictedContactMode || "").trim();
      if (observedMode && predictedMode) {
        contactModeCount += 1;
        if (observedMode === predictedMode) contactModeMatchCount += 1;
      }
      const observedLockHeld = String(row.lockHeldObserved || "").trim().toLowerCase();
      if (observedLockHeld) {
        lockHeldCount += 1;
        const expectedHeld = String(predicted.lockMode || row.lockMode || "").includes("locked") && predicted.predictedLockPass !== false;
        const observedHeld = ["held", "true", "yes", "1", "pass"].includes(observedLockHeld);
        if (observedHeld === expectedHeld) lockHeldMatchCount += 1;
      }
      comparisonRows.push({
        caseId: row.caseId,
        connector: row.connector,
        matched: true,
        observedScalarCount: rowScalarCount,
        maxAbsError: rowScalarCount ? rowMax : "",
        observedContactMode: observedMode,
        predictedContactMode: predictedMode,
        lockHeldObserved: row.lockHeldObserved || "",
      });
    }
    const rmsError = observedScalarCount ? Math.sqrt(squaredError / observedScalarCount) : null;
    const cellRmsError = cellCount ? Math.sqrt(cellSquaredError / cellCount) : null;
    const connectorRmsError = connectorCount ? Math.sqrt(connectorSquaredError / connectorCount) : null;
    const contactModeAccuracy = contactModeCount ? contactModeMatchCount / contactModeCount : null;
    const lockHeldAccuracy = lockHeldCount ? lockHeldMatchCount / lockHeldCount : null;
    const passesTolerance = observedScalarCount > 0 && rmsError <= tolerance && maxAbsError <= tolerance;
    const comparisonReady = observedScalarCount > 0 && unmatchedRowCount === 0;
    return {
      schema: TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA,
      model: "browser-two-cell-fidelity-matrix-measured-vs-reduced-proxy-comparison",
      cadReference: CAD_RAD_CELL_REFERENCE,
      rows: comparisonRows,
      summary: {
        comparisonReady,
        passesTolerance,
        requiresMoreData: observedScalarCount === 0,
        matchedRowCount,
        unmatchedRowCount,
        expectedRowCount: template.rows.length,
        providedRowCount: rows.length,
        observedScalarCount,
        missingObservationCount,
        rmsError,
        cellRmsError,
        connectorRmsError,
        maxAbsError: observedScalarCount ? maxAbsError : null,
        contactModeAccuracy,
        lockHeldAccuracy,
        tolerance,
        physicalAccuracyValidated: false,
        claimLabel: "validates reduced proxy against supplied measurements; exact real-life physics still requires segmented CAD/contact evidence",
      },
    };
  }

  function exportTwoCellFidelityMatrixMeasurementComparisonJson(state, input, options = {}) {
    return JSON.stringify(compareTwoCellFidelityMatrixMeasurements(state, input, options), null, 2);
  }

  function mean(values) {
    return values.length ? values.reduce((total, value) => total + value, 0) / values.length : null;
  }

  function fitScale(samples) {
    const usable = samples.filter(([predicted]) => Math.abs(finiteNumber(predicted, 0)) > 1e-12);
    if (!usable.length) {
      return {
        estimate: null,
        sampleCount: 0,
        status: "insufficient-excitation",
        meanResidual: null,
        rmsResidual: null,
      };
    }
    const numerator = usable.reduce((total, [predicted, observed]) => total + predicted * observed, 0);
    const denominator = usable.reduce((total, [predicted]) => total + predicted * predicted, 0);
    const estimate = denominator > 0 ? numerator / denominator : null;
    const residuals = usable.map(([predicted, observed]) => observed - predicted);
    const squared = residuals.reduce((total, value) => total + value * value, 0);
    return {
      estimate,
      sampleCount: usable.length,
      status: estimate === null ? "insufficient-excitation" : "estimated",
      meanResidual: mean(residuals),
      rmsResidual: Math.sqrt(squared / usable.length),
    };
  }

  function bounded(value, lower, upper, fallback) {
    const numeric = Number(value);
    return Number.isFinite(numeric) ? clamp(numeric, lower, upper) : fallback;
  }

  function axisResidualGroups(rows, expected, axis) {
    const buckets = new Map();
    for (const row of rows) {
      const predicted = expected.get(`${row.caseId}|${row.connector}`);
      if (!predicted) continue;
      const axisValue = numericOrNull(row[axis] ?? predicted[axis]);
      if (axisValue === null) continue;
      if (!buckets.has(axisValue)) {
        buckets.set(axisValue, {
          [axis]: axisValue,
          observedScalarCount: 0,
          rightZResiduals: [],
          rightAlphaResiduals: [],
          verticalSlipResiduals: [],
          lateralSlipResiduals: [],
        });
      }
      const bucket = buckets.get(axisValue);
      for (const [observedField, predictedField, residualField] of [
        ["observedRightCellZ", "predictedRightCellZ", "rightZResiduals"],
        ["observedRightAlpha", "predictedRightAlpha", "rightAlphaResiduals"],
        ["observedVerticalSlipMm", "predictedVerticalSlipMm", "verticalSlipResiduals"],
        ["observedLateralSlipMm", "predictedLateralSlipMm", "lateralSlipResiduals"],
      ]) {
        const observed = numericOrNull(row[observedField]);
        const target = numericOrNull(predicted[predictedField] ?? row[predictedField]);
        if (observed === null || target === null) continue;
        bucket[residualField].push(observed - target);
        bucket.observedScalarCount += 1;
      }
    }
    return [...buckets.values()]
      .sort((a, b) => finiteNumber(a[axis], 0) - finiteNumber(b[axis], 0))
      .map((bucket) => ({
        [axis]: bucket[axis],
        observedScalarCount: bucket.observedScalarCount,
        meanRightZResidual: mean(bucket.rightZResiduals),
        meanRightAlphaResidual: mean(bucket.rightAlphaResiduals),
        meanVerticalSlipResidualMm: mean(bucket.verticalSlipResiduals),
        meanLateralSlipResidualMm: mean(bucket.lateralSlipResiduals),
      }));
  }

  function calibrateTwoCellFidelityMatrixParameters(state, input, options = {}) {
    const rows = measurementRows(input);
    const template = twoCellFidelityMatrixMeasurementTemplate(state, options);
    const expected = new Map(template.rows.map((row) => [`${row.caseId}|${row.connector}`, row]));
    const comparison = compareTwoCellFidelityMatrixMeasurements(state, rows, options);
    const initialAlpha = finiteNumber(state.grid.initialAlpha, 1);
    const alphaSamples = [];
    const thetaSamples = [];
    const zSamples = [];
    const lateralSlipSamples = [];
    const verticalSlipSamples = [];
    const totalSlipSamples = [];
    const alphaBiases = [];
    const zBiases = [];
    const verticalSlipBiases = [];
    for (const row of rows) {
      const predicted = expected.get(`${row.caseId}|${row.connector}`);
      if (!predicted) continue;
      for (const [observedField, predictedField] of [
        ["observedLeftAlpha", "predictedLeftAlpha"],
        ["observedRightAlpha", "predictedRightAlpha"],
      ]) {
        const observed = numericOrNull(row[observedField]);
        const target = numericOrNull(predicted[predictedField] ?? row[predictedField]);
        if (observed === null || target === null) continue;
        alphaSamples.push([target - initialAlpha, observed - initialAlpha]);
        alphaBiases.push(observed - target);
      }
      for (const [observedField, predictedField] of [
        ["observedLeftTheta", "predictedLeftTheta"],
        ["observedRightTheta", "predictedRightTheta"],
      ]) {
        const observed = numericOrNull(row[observedField]);
        const target = numericOrNull(predicted[predictedField] ?? row[predictedField]);
        if (observed !== null && target !== null) thetaSamples.push([target, observed]);
      }
      for (const [observedField, predictedField] of [
        ["observedLeftCellZ", "predictedLeftCellZ"],
        ["observedRightCellZ", "predictedRightCellZ"],
      ]) {
        const observed = numericOrNull(row[observedField]);
        const target = numericOrNull(predicted[predictedField] ?? row[predictedField]);
        if (observed === null || target === null) continue;
        zSamples.push([target, observed]);
        zBiases.push(observed - target);
      }
      for (const [observedField, predictedField, samples] of [
        ["observedLateralSlipMm", "predictedLateralSlipMm", lateralSlipSamples],
        ["observedVerticalSlipMm", "predictedVerticalSlipMm", verticalSlipSamples],
        ["observedTotalSlipMm", "predictedTotalSlipMm", totalSlipSamples],
      ]) {
        const observed = numericOrNull(row[observedField]);
        const target = numericOrNull(predicted[predictedField] ?? row[predictedField]);
        if (observed === null || target === null) continue;
        samples.push([target, observed]);
        if (observedField === "observedVerticalSlipMm") verticalSlipBiases.push(observed - target);
      }
    }
    const alphaFit = fitScale(alphaSamples);
    const thetaFit = fitScale(thetaSamples);
    const zFit = fitScale(zSamples);
    const lateralSlipFit = fitScale(lateralSlipSamples);
    const verticalSlipFit = fitScale(verticalSlipSamples);
    const totalSlipFit = fitScale(totalSlipSamples);
    const alphaScale = bounded(alphaFit.estimate, 0.25, 4, 1);
    const zScale = bounded(zFit.estimate, 0.25, 4, 1);
    const verticalSlipScale = bounded(verticalSlipFit.estimate, 0.25, 4, 1);
    const pinRadius = Math.max(0, finiteNumber(state.grid.pinRadius, 0.18));
    const holeRadius = Math.max(pinRadius, finiteNumber(state.grid.holeRadius, 0.225));
    const currentClearance = Math.max(0, holeRadius - pinRadius);
    const proposedClearance = Math.max(0, currentClearance * verticalSlipScale);
    const proposed = {
      couplingGain: bounded(finiteNumber(state.grid.couplingGain, 0.55) * alphaScale, 0, 1, finiteNumber(state.grid.couplingGain, 0.55)),
      zCouplingGain: bounded(finiteNumber(state.grid.zCouplingGain, 0.32) * zScale, 0, 1, finiteNumber(state.grid.zCouplingGain, 0.32)),
      holeRadius: Math.max(pinRadius, pinRadius + proposedClearance),
      pinHoleClearance: proposedClearance,
      applyAutomatically: false,
      reason: "reduced-proxy estimate only; require holdout bench and segmented contact validation before applying to physical claims",
    };
    const ready = comparison.summary.observedScalarCount >= 24 && comparison.summary.matchedRowCount > 0;
    return {
      schema: TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA,
      model: "browser-two-cell-fidelity-matrix-reduced-proxy-calibration-estimate",
      cadReference: CAD_RAD_CELL_REFERENCE,
      comparisonSummary: comparison.summary,
      estimates: {
        alphaResponseScale: alphaFit,
        thetaScale: thetaFit,
        zResponseScale: zFit,
        lateralSlipScale: lateralSlipFit,
        verticalSlipScale: verticalSlipFit,
        totalSlipScale: totalSlipFit,
        meanAlphaBias: mean(alphaBiases),
        meanZBias: mean(zBiases),
        meanVerticalSlipBiasMm: mean(verticalSlipBiases),
      },
      axisDiagnostics: {
        byHoleRadius: axisResidualGroups(rows, expected, "holeRadius"),
        byBacklash: axisResidualGroups(rows, expected, "backlash"),
      },
      proposedReducedProxyUpdates: proposed,
      summary: {
        readyForReducedProxyCalibration: ready,
        requiresMoreData: !ready,
        requiresHoldoutValidation: true,
        requiresSegmentedPhysics: true,
        physicalAccuracyValidated: false,
        zResponseScale: zFit.estimate,
        verticalSlipScale: verticalSlipFit.estimate,
        proposedZCouplingGain: proposed.zCouplingGain,
        proposedHoleRadius: proposed.holeRadius,
      },
      claimLabels: {
        calibration: "fits scale/bias terms for the browser reduced two-cell proxy only",
        holeRadius: "vertical slip scale suggests an effective clearance update, not an exact CAD hole dimension",
        physicalAccuracy: "not validated until segmented CAD/contact and independent bench measurements agree",
      },
    };
  }

  function exportTwoCellFidelityMatrixParameterCalibrationJson(state, input, options = {}) {
    return JSON.stringify(calibrateTwoCellFidelityMatrixParameters(state, input, options), null, 2);
  }

  function twoCellConnectorMeasurementTemplate(state, options = {}) {
    const specs = twoCellPhysicalSuiteCaseSpecs(state, options);
    const rows = [];
    for (const spec of specs) {
      const caseState = stateForTwoCellCase(state, spec);
      const caseControls = controlsForLockMode(
        twoCellBenchControls(state, {
          ...options,
          alphaCommand: spec.alphaCommand,
          zCommand: spec.zCommand,
        }),
        spec.lockMode
      );
      const contact = twoCellConnectorContactReport(caseState, caseControls);
      for (const connector of contact.connectors || []) {
        rows.push({
          caseId: spec.caseId,
          connector: connector.connector,
          leftSite: connector.leftSite,
          rightSite: connector.rightSite,
          lockMode: spec.lockMode,
          alphaCommand: caseControls.alphaCommand,
          zCommand: caseControls.zCommand,
          backlash: caseState.grid.backlash,
          pinRadius: contact.dimensions.pinRadius,
          holeRadius: contact.dimensions.holeRadius,
          pinHoleClearanceMm: contact.dimensions.pinHoleClearanceMm,
          predictedLeftXmm: connector.leftPositionMm.x,
          predictedLeftYmm: connector.leftPositionMm.y,
          predictedLeftZmm: connector.leftPositionMm.z,
          predictedRightXmm: connector.rightPositionMm.x,
          predictedRightYmm: connector.rightPositionMm.y,
          predictedRightZmm: connector.rightPositionMm.z,
          predictedSlipXmm: connector.slipMm.x,
          predictedSlipYmm: connector.slipMm.y,
          predictedSlipZmm: connector.slipMm.z,
          predictedLateralSlipMm: connector.lateralSlipMm,
          predictedVerticalSlipMm: connector.verticalSlipMm,
          predictedTotalSlipMm: connector.totalSlipMm,
          predictedLateralExcessMm: connector.lateralExcessMm,
          predictedVerticalExcessMm: connector.verticalExcessMm,
          predictedContactMode: connector.contactMode,
          observedLeftXmm: "",
          observedLeftYmm: "",
          observedLeftZmm: "",
          observedRightXmm: "",
          observedRightYmm: "",
          observedRightZmm: "",
          observedLateralSlipMm: "",
          observedVerticalSlipMm: "",
          observedTotalSlipMm: "",
          observedContactMode: "",
          lockHeldObserved: "",
          measurementSource: "",
          notes: "",
        });
      }
    }
    return {
      schema: "rad-sim.two-cell-connector-measurement-template.v1",
      model: "browser-two-cell-suite-connector-slip-fillable-template",
      cadReference: CAD_RAD_CELL_REFERENCE,
      rows,
      summary: {
        caseCount: new Set(rows.map((row) => row.caseId)).size,
        connectorRowCount: rows.length,
        connectorNames: [...new Set(rows.map((row) => row.connector))].sort(),
        measurementRequired: true,
        physicalAccuracyValidated: false,
      },
    };
  }

  function exportTwoCellConnectorMeasurementTemplateJson(state, options = {}) {
    return JSON.stringify(twoCellConnectorMeasurementTemplate(state, options), null, 2);
  }

  function exportTwoCellConnectorMeasurementTemplateCsv(state, options = {}) {
    const template = twoCellConnectorMeasurementTemplate(state, options);
    const columns = [
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
    ];
    return [
      columns.join(","),
      ...template.rows.map((row) =>
        columns
          .map((column) => {
            const value = row[column];
            return typeof value === "string" ? `"${value.replace(/"/g, '""')}"` : String(value);
          })
          .join(",")
      ),
    ].join("\n");
  }

  RAD.TWO_CELL_BENCH_SCHEMA = TWO_CELL_BENCH_SCHEMA;
  RAD.TWO_CELL_SWEEP_SCHEMA = TWO_CELL_SWEEP_SCHEMA;
  RAD.CAD_RAD_CELL_LAYOUT_SCHEMA = CAD_RAD_CELL_LAYOUT_SCHEMA;
  RAD.CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA = CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA;
  RAD.TWO_CELL_CONNECTOR_CONTACT_SCHEMA = TWO_CELL_CONNECTOR_CONTACT_SCHEMA;
  RAD.TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA = TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA;
  RAD.TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA = TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA;
  RAD.TWO_CELL_CONTACT_PHASE_MAP_SCHEMA = TWO_CELL_CONTACT_PHASE_MAP_SCHEMA;
  RAD.TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA = TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA;
  RAD.TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA = TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA;
  RAD.TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA = TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA;
  RAD.TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA = TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA;
  RAD.TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA = TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA;
  RAD.TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA = TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA;
  RAD.TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA = TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA;
  RAD.TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA = TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA;
  RAD.TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA = TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA;
  RAD.TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA = TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA;
  RAD.TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA = TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA;
  RAD.CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA = CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA;
  RAD.TWO_CELL_PHYSICAL_FIDELITY_STATUS_SCHEMA = TWO_CELL_PHYSICAL_FIDELITY_STATUS_SCHEMA;
  RAD.TWO_CELL_EXTERNAL_FIDELITY_WEB_SUMMARY_SCHEMA = TWO_CELL_EXTERNAL_FIDELITY_WEB_SUMMARY_SCHEMA;
  RAD.CAD_RAD_CELL_REFERENCE = CAD_RAD_CELL_REFERENCE;
  RAD.CAD_RAD_CELL_ARCHIVE_AUDIT = CAD_RAD_CELL_ARCHIVE_AUDIT;
  RAD.DEFAULT_TWO_CELL_CONTROLS = DEFAULT_TWO_CELL_CONTROLS;
  RAD.cadRadCellArchiveAudit = cadRadCellArchiveAudit;
  RAD.twoCellExternalFidelityWebSummary = twoCellExternalFidelityWebSummary;
  RAD.twoCellExternalFidelityCaseOptions = twoCellExternalFidelityCaseOptions;
  RAD.twoCellExternalFidelityCase = twoCellExternalFidelityCase;
  RAD.twoCellExternalFidelityCaseControls = twoCellExternalFidelityCaseControls;
  RAD.exportTwoCellExternalFidelityWebSummaryJson = exportTwoCellExternalFidelityWebSummaryJson;
  RAD.cadRadCellLayout = cadRadCellLayout;
  RAD.cadRadCellReferenceProfile = cadRadCellReferenceProfile;
  RAD.exportCadRadCellReferenceProfileJson = exportCadRadCellReferenceProfileJson;
  RAD.twoCellPhysicalFidelityStatus = twoCellPhysicalFidelityStatus;
  RAD.twoCellConnectorContactReport = twoCellConnectorContactReport;
  RAD.twoCellBenchControls = twoCellBenchControls;
  RAD.twoCellPhysicalSuiteCaseSpecs = twoCellPhysicalSuiteCaseSpecs;
  RAD.twoCellPhysicalSuiteCaseOptions = twoCellPhysicalSuiteCaseOptions;
  RAD.twoCellPhysicalSimulationSuite = twoCellPhysicalSimulationSuite;
  RAD.twoCellPhysicalFidelityMatrix = twoCellPhysicalFidelityMatrix;
  RAD.twoCellContactPhaseMap = twoCellContactPhaseMap;
  RAD.twoCellRadiusBacklashPhaseDiagram = twoCellRadiusBacklashPhaseDiagram;
  RAD.twoCellRadiusBacklashTransitionReport = twoCellRadiusBacklashTransitionReport;
  RAD.compareTwoCellRadiusBacklashTransitionMeasurements = compareTwoCellRadiusBacklashTransitionMeasurements;
  RAD.twoCellRadiusBacklashTransitionRerun = twoCellRadiusBacklashTransitionRerun;
  RAD.twoCellPhysicalResponseAtlas = twoCellPhysicalResponseAtlas;
  RAD.twoCellCadContactDecompositionSpec = twoCellCadContactDecompositionSpec;
  RAD.twoCellExactContactHandoffPlan = twoCellExactContactHandoffPlan;
  RAD.twoCellFidelityMatrixMeasurementTemplate = twoCellFidelityMatrixMeasurementTemplate;
  RAD.twoCellExternalFidelityMatrixManifest = twoCellExternalFidelityMatrixManifest;
  RAD.twoCellExternalFidelityMatrixManifestRows = twoCellExternalFidelityMatrixManifestRows;
  RAD.twoCellFidelityMatrixMeasurementsFromCsv = twoCellFidelityMatrixMeasurementsFromCsv;
  RAD.compareTwoCellFidelityMatrixMeasurements = compareTwoCellFidelityMatrixMeasurements;
  RAD.calibrateTwoCellFidelityMatrixParameters = calibrateTwoCellFidelityMatrixParameters;
  RAD.twoCellConnectorMeasurementTemplate = twoCellConnectorMeasurementTemplate;
  RAD.simulateTwoCellBench = simulateTwoCellBench;
  RAD.sweepTwoCellBacklash = sweepTwoCellBacklash;
  RAD.exportTwoCellBenchJson = exportTwoCellBenchJson;
  RAD.exportCadRadCellArchiveAuditJson = exportCadRadCellArchiveAuditJson;
  RAD.exportTwoCellPhysicalFidelityStatusJson = exportTwoCellPhysicalFidelityStatusJson;
  RAD.exportTwoCellBacklashSweepCsv = exportTwoCellBacklashSweepCsv;
  RAD.exportTwoCellPhysicalSuiteJson = exportTwoCellPhysicalSuiteJson;
  RAD.exportTwoCellPhysicalSuiteCsv = exportTwoCellPhysicalSuiteCsv;
  RAD.exportTwoCellPhysicalFidelityMatrixJson = exportTwoCellPhysicalFidelityMatrixJson;
  RAD.exportTwoCellPhysicalFidelityMatrixCsv = exportTwoCellPhysicalFidelityMatrixCsv;
  RAD.exportTwoCellContactPhaseMapJson = exportTwoCellContactPhaseMapJson;
  RAD.exportTwoCellContactPhaseMapCsv = exportTwoCellContactPhaseMapCsv;
  RAD.exportTwoCellRadiusBacklashPhaseDiagramJson = exportTwoCellRadiusBacklashPhaseDiagramJson;
  RAD.exportTwoCellRadiusBacklashPhaseDiagramCsv = exportTwoCellRadiusBacklashPhaseDiagramCsv;
  RAD.exportTwoCellRadiusBacklashTransitionReportJson = exportTwoCellRadiusBacklashTransitionReportJson;
  RAD.exportTwoCellRadiusBacklashTransitionReportCsv = exportTwoCellRadiusBacklashTransitionReportCsv;
  RAD.exportTwoCellRadiusBacklashTransitionComparisonJson = exportTwoCellRadiusBacklashTransitionComparisonJson;
  RAD.exportTwoCellRadiusBacklashTransitionComparisonCsv = exportTwoCellRadiusBacklashTransitionComparisonCsv;
  RAD.exportTwoCellRadiusBacklashTransitionRerunJson = exportTwoCellRadiusBacklashTransitionRerunJson;
  RAD.exportTwoCellRadiusBacklashTransitionRerunCsv = exportTwoCellRadiusBacklashTransitionRerunCsv;
  RAD.exportTwoCellPhysicalResponseAtlasJson = exportTwoCellPhysicalResponseAtlasJson;
  RAD.exportTwoCellPhysicalResponseAtlasCsv = exportTwoCellPhysicalResponseAtlasCsv;
  RAD.exportTwoCellCadContactDecompositionJson = exportTwoCellCadContactDecompositionJson;
  RAD.exportTwoCellCadContactDecompositionCsv = exportTwoCellCadContactDecompositionCsv;
  RAD.exportTwoCellExactContactHandoffPlanJson = exportTwoCellExactContactHandoffPlanJson;
  RAD.exportTwoCellExactContactHandoffPlanCsv = exportTwoCellExactContactHandoffPlanCsv;
  RAD.exportTwoCellFidelityMatrixMeasurementTemplateJson = exportTwoCellFidelityMatrixMeasurementTemplateJson;
  RAD.exportTwoCellFidelityMatrixMeasurementTemplateCsv = exportTwoCellFidelityMatrixMeasurementTemplateCsv;
  RAD.exportTwoCellExternalFidelityMatrixManifestJson = exportTwoCellExternalFidelityMatrixManifestJson;
  RAD.exportTwoCellExternalFidelityMatrixManifestCsv = exportTwoCellExternalFidelityMatrixManifestCsv;
  RAD.exportTwoCellFidelityMatrixMeasurementComparisonJson = exportTwoCellFidelityMatrixMeasurementComparisonJson;
  RAD.exportTwoCellFidelityMatrixParameterCalibrationJson = exportTwoCellFidelityMatrixParameterCalibrationJson;
  RAD.exportTwoCellConnectorMeasurementTemplateJson = exportTwoCellConnectorMeasurementTemplateJson;
  RAD.exportTwoCellConnectorMeasurementTemplateCsv = exportTwoCellConnectorMeasurementTemplateCsv;
})();
