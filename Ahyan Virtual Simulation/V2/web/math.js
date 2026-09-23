(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function backlashActivation(x, backlash) {
    return Math.max(0, x - backlash) + Math.min(x + backlash, 0);
  }

  function alphaToTheta(alpha) {
    return 70 * alpha - 60;
  }

  function normalizedBacklashToThetaDeadZone(backlash) {
    return (Math.asin(clamp(Math.max(0, Number(backlash) || 0), 0, 1)) * 180) / Math.PI;
  }

  function normalizedBacklashToAlphaDeadZone(backlash) {
    return normalizedBacklashToThetaDeadZone(backlash) / 70;
  }

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function commandLimits(state) {
    return {
      z: Math.max(0.05, Number(state.grid.zTravelLimit || 0.8)),
      alphaContract: Math.max(0.05, Number(state.grid.alphaContractLimit || 0.55)),
      alphaExpand: Math.max(0.05, Number(state.grid.alphaExpandLimit || 0.35)),
    };
  }

  function clampCommandZ(state, value) {
    const limit = commandLimits(state).z;
    return clamp(Number(value) || 0, -limit, limit);
  }

  function clampCommandAlpha(state, value) {
    const limits = commandLimits(state);
    return clamp(Number(value) || 0, -limits.alphaContract, limits.alphaExpand);
  }

  function cellRemoved(state, r, c) {
    return Boolean(state.cells.removed?.[r]?.[c]);
  }

  function cellPresent(state, r, c) {
    return !cellRemoved(state, r, c);
  }

  function commandSaturation(state, alphaCommand, zCommand) {
    const limits = commandLimits(state);
    const alphaLimit = alphaCommand < 0 ? limits.alphaContract : limits.alphaExpand;
    const alphaSat = Math.abs(alphaCommand) / Math.max(1e-9, alphaLimit);
    const zSat = Math.abs(zCommand) / Math.max(1e-9, limits.z);
    return Math.max(alphaSat, zSat);
  }

  function pinHoleClearance(state) {
    const pinRadius = Math.max(0, Number(state.grid.pinRadius ?? 0.18));
    const holeRadius = Math.max(pinRadius, Number(state.grid.holeRadius ?? 0.225));
    return Math.max(0, holeRadius - pinRadius);
  }

  const HARDWARE_PROFILE_DIMENSIONS = Object.freeze([
    "pinRadiusMm",
    "holeRadiusMm",
    "plateThicknessMm",
    "jointStackHeightMm",
    "bossRadiusMm",
  ]);
  const HARDWARE_PROFILE_SCHEMA = "rad-sim.hardware-profile.v1";
  const CALIBRATION_VISUAL_FIELDS = Object.freeze([
    "pinRadiusMm",
    "holeRadiusMm",
    "plateThicknessMm",
    "jointStackHeightMm",
  ]);
  const CALIBRATION_MESH_FIELDS = HARDWARE_PROFILE_DIMENSIONS;
  const CALIBRATION_SOLVER_GAPS = Object.freeze([
    "measured axial and hinge stiffness",
    "actuator force/stroke calibration",
    "friction and contact characterization",
    "measured single, pair, and cluster response data",
  ]);
  const CALIBRATION_FIELD_LABELS = Object.freeze({
    pinRadiusMm: "hinge pin radius",
    holeRadiusMm: "mating hole radius",
    plateThicknessMm: "plate thickness",
    jointStackHeightMm: "joint stack height",
    bossRadiusMm: "pivot boss radius",
  });
  const CALIBRATION_SOLVER_TASKS = Object.freeze([
    {
      id: "solver_axial_hinge_stiffness",
      label: "Measure axial and hinge stiffness",
      notes: "Run isolated bar and hinge displacement tests before fitting spring-hinge constants.",
    },
    {
      id: "solver_actuator_force_stroke",
      label: "Measure actuator force and stroke",
      notes: "Record commanded stroke, realized travel, and force limits for the installed linear actuators.",
    },
    {
      id: "solver_friction_contact",
      label: "Characterize friction and contact",
      notes: "Measure pin-hole slip, contact onset, and hysteresis under repeated actuation.",
    },
    {
      id: "solver_response_data",
      label: "Measure single, pair, and cluster response",
      notes: "Capture alpha and vertical displacement fields for one cell, adjacent pairs, and small clusters.",
    },
  ]);

  function nonNegativeNumberOrNull(value) {
    if (value === null || value === undefined || value === "") return null;
    const numeric = Number(value);
    return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
  }

  function hardwareProfile(state) {
    const raw = state.grid.hardwareProfile || {};
    const sideLengthMm = Math.max(1e-9, Number(raw.sideLengthMm ?? state.grid.paperSideLengthMm ?? 35));
    const fabricationHoleToleranceMm = Math.max(0, Number(raw.fabricationHoleToleranceMm ?? state.grid.paperHoleToleranceMm ?? 0.1));
    const pinRadiusMm = nonNegativeNumberOrNull(raw.pinRadiusMm);
    let holeRadiusMm = nonNegativeNumberOrNull(raw.holeRadiusMm);
    if (pinRadiusMm !== null && holeRadiusMm !== null && holeRadiusMm < pinRadiusMm) holeRadiusMm = pinRadiusMm;
    return {
      name: String(raw.name || "paper-reference"),
      source: String(raw.source || "RAD preprint defaults"),
      sideLengthMm,
      fabricationHoleToleranceMm,
      backlashMm: nonNegativeNumberOrNull(raw.backlashMm),
      pinRadiusMm,
      holeRadiusMm,
      plateThicknessMm: nonNegativeNumberOrNull(raw.plateThicknessMm),
      jointStackHeightMm: nonNegativeNumberOrNull(raw.jointStackHeightMm),
      bossRadiusMm: nonNegativeNumberOrNull(raw.bossRadiusMm),
      notes: String(raw.notes || ""),
    };
  }

  function hardwareProfileSchemaPayload(state) {
    const profile = hardwareProfile(state);
    const measuredFields = hardwareMeasuredFields(profile);
    const missingFields = hardwareMissingFields(profile);
    const pinHoleClearanceMm =
      profile.pinRadiusMm !== null && profile.holeRadiusMm !== null ? profile.holeRadiusMm - profile.pinRadiusMm : null;
    return {
      schema: HARDWARE_PROFILE_SCHEMA,
      name: profile.name,
      source: profile.source,
      units: "mm",
      dimensionsMm: {
        sideLengthMm: profile.sideLengthMm,
        fabricationHoleToleranceMm: profile.fabricationHoleToleranceMm,
        backlashMm: profile.backlashMm,
        pinRadiusMm: profile.pinRadiusMm,
        holeRadiusMm: profile.holeRadiusMm,
        plateThicknessMm: profile.plateThicknessMm,
        jointStackHeightMm: profile.jointStackHeightMm,
        bossRadiusMm: profile.bossRadiusMm,
      },
      derived: {
        pinHoleClearanceMm,
      },
      measuredFields,
      missingFields,
      coverageRatio: measuredFields.length / HARDWARE_PROFILE_DIMENSIONS.length,
      notes: profile.notes,
    };
  }

  function profileValue(payload, dimensions, camelKey, snakeKey, fallback) {
    const value = dimensions?.[camelKey] ?? payload?.[camelKey] ?? payload?.[snakeKey] ?? fallback;
    return Number.isFinite(Number(value)) ? Number(value) : fallback;
  }

  function optionalProfileValue(payload, dimensions, camelKey, snakeKey) {
    return nonNegativeNumberOrNull(dimensions?.[camelKey] ?? payload?.[camelKey] ?? payload?.[snakeKey]);
  }

  function hardwareProfileFromObject(payload, state = null) {
    if (!payload || typeof payload !== "object" || Array.isArray(payload)) {
      throw new Error("Hardware profile JSON must contain an object.");
    }
    if (payload.schema === "rad-sim.cad-mesh-audit.v1" && payload.hardwareProfile && typeof payload.hardwareProfile === "object") {
      return hardwareProfileFromObject(payload.hardwareProfile, state);
    }
    if (
      payload.schema === "rad-sim.cad-rad-cell-reference-profile.v1" &&
      payload.hardwareProfile &&
      typeof payload.hardwareProfile === "object"
    ) {
      return hardwareProfileFromObject(payload.hardwareProfile, state);
    }
    if (payload.schema && payload.schema !== HARDWARE_PROFILE_SCHEMA) {
      throw new Error(`Unsupported hardware profile schema: ${payload.schema}`);
    }
    const base = state ? hardwareProfile(state) : RAD.defaultHardwareProfile?.() || {};
    const dimensions = payload.dimensionsMm && typeof payload.dimensionsMm === "object" ? payload.dimensionsMm : payload;
    const profile = {
      name: String(payload.name || base.name || "measured-hardware-profile"),
      source: String(payload.source || base.source || "external JSON"),
      sideLengthMm: Math.max(1e-9, profileValue(payload, dimensions, "sideLengthMm", "side_length_mm", base.sideLengthMm || 35)),
      fabricationHoleToleranceMm: Math.max(
        0,
        profileValue(
          payload,
          dimensions,
          "fabricationHoleToleranceMm",
          "fabrication_hole_tolerance_mm",
          base.fabricationHoleToleranceMm || 0.1
        )
      ),
      backlashMm: optionalProfileValue(payload, dimensions, "backlashMm", "backlash_mm"),
      pinRadiusMm: optionalProfileValue(payload, dimensions, "pinRadiusMm", "pin_radius_mm"),
      holeRadiusMm: optionalProfileValue(payload, dimensions, "holeRadiusMm", "hole_radius_mm"),
      plateThicknessMm: optionalProfileValue(payload, dimensions, "plateThicknessMm", "plate_thickness_mm"),
      jointStackHeightMm: optionalProfileValue(payload, dimensions, "jointStackHeightMm", "joint_stack_height_mm"),
      bossRadiusMm: optionalProfileValue(payload, dimensions, "bossRadiusMm", "boss_radius_mm"),
      notes: String(payload.notes || base.notes || ""),
    };
    if (profile.pinRadiusMm !== null && profile.holeRadiusMm !== null && profile.holeRadiusMm < profile.pinRadiusMm) {
      profile.holeRadiusMm = profile.pinRadiusMm;
    }
    return profile;
  }

  function importHardwareProfileJson(text, state = null) {
    const parsed = JSON.parse(String(text || "").replace(/^\uFEFF/, ""));
    return hardwareProfileFromObject(parsed, state);
  }

  function exportHardwareProfileJson(state) {
    return JSON.stringify(hardwareProfileSchemaPayload(state), null, 2);
  }

  function hardwareMeasuredFields(profile) {
    const p = profile || {};
    return HARDWARE_PROFILE_DIMENSIONS.filter((field) => p[field] !== null && p[field] !== undefined);
  }

  function hardwareMissingFields(profile) {
    const p = profile || {};
    return HARDWARE_PROFILE_DIMENSIONS.filter((field) => p[field] === null || p[field] === undefined);
  }

  function calibrationProfileSummary(state) {
    const profile = hardwareProfile(state);
    const measuredFields = hardwareMeasuredFields(profile);
    const missingFields = hardwareMissingFields(profile);
    const modelScale = Math.max(1e-9, Number(state.grid.cellSize || 1)) / profile.sideLengthMm;
    const pinHoleClearanceMm =
      profile.pinRadiusMm !== null && profile.holeRadiusMm !== null ? profile.holeRadiusMm - profile.pinRadiusMm : null;
    return {
      profile,
      measuredFields,
      missingFields,
      measuredCount: measuredFields.length,
      totalCount: HARDWARE_PROFILE_DIMENSIONS.length,
      coverageRatio: measuredFields.length / HARDWARE_PROFILE_DIMENSIONS.length,
      pinHoleClearanceMm,
      pinRadiusModel: profile.pinRadiusMm === null ? null : profile.pinRadiusMm * modelScale,
      holeRadiusModel: profile.holeRadiusMm === null ? null : profile.holeRadiusMm * modelScale,
      backlashModel: profile.backlashMm === null ? null : profile.backlashMm / profile.sideLengthMm,
    };
  }

  function calibrationReadiness(state) {
    const summary = calibrationProfileSummary(state);
    const hasValue = (field) => summary.profile[field] !== null && summary.profile[field] !== undefined;
    const visualMissingFields = CALIBRATION_VISUAL_FIELDS.filter((field) => !hasValue(field));
    const meshMissingFields = CALIBRATION_MESH_FIELDS.filter((field) => !hasValue(field));
    const visualReady = visualMissingFields.length === 0;
    const meshReady = meshMissingFields.length === 0;
    let level = "paper-scale";
    if (meshReady) level = "mesh-calibrated";
    else if (visualReady) level = "visual-calibrated";
    else if (summary.measuredCount > 0) level = "partial-measured";
    const summaryText = meshReady
      ? "mesh-calibrated geometry; solver still needs physical response calibration"
      : visualReady
        ? "visual-calibrated geometry; mesh export still has missing dimensions"
        : summary.measuredCount > 0
          ? "partial measured geometry; solver still needs geometry and response calibration"
          : "paper-scale defaults only; solver still lacks measured hardware geometry";
    return {
      profileName: summary.profile.name,
      level,
      measuredFields: summary.measuredFields,
      missingFields: summary.missingFields,
      visualMissingFields,
      meshMissingFields,
      solverGaps: [...CALIBRATION_SOLVER_GAPS],
      coverageRatio: summary.coverageRatio,
      visualReady,
      meshReady,
      solverReady: false,
      summary: summaryText,
    };
  }

  function calibrationMeasurementPlan(state) {
    const profile = hardwareProfile(state);
    const tasks = [];
    for (const [index, field] of HARDWARE_PROFILE_DIMENSIONS.entries()) {
      const measured = profile[field] !== null && profile[field] !== undefined;
      const label = CALIBRATION_FIELD_LABELS[field] || field;
      tasks.push({
        id: `geometry_${field}`,
        label: `Measure ${label}`,
        category: "geometry",
        status: measured ? "done" : "missing",
        priority: measured ? 101 + index : 1 + index,
        evidenceField: field,
        notes: measured
          ? "Available in the active hardware profile."
          : "Required for real-cell visual/export geometry before solver calibration.",
      });
    }
    for (const [index, task] of CALIBRATION_SOLVER_TASKS.entries()) {
      tasks.push({
        id: task.id,
        label: task.label,
        category: "solver",
        status: "missing",
        priority: 51 + index,
        evidenceField: null,
        notes: task.notes,
      });
    }
    return tasks.sort((a, b) => {
      if (a.status !== b.status) return a.status === "done" ? 1 : -1;
      if (a.priority !== b.priority) return a.priority - b.priority;
      return a.id.localeCompare(b.id);
    });
  }

  function exportCalibrationMeasurementPlan(state) {
    return JSON.stringify(
      {
        schema: "rad-sim.calibration-plan.v1",
        profile: hardwareProfile(state),
        readiness: calibrationReadiness(state),
        tasks: calibrationMeasurementPlan(state),
      },
      null,
      2
    );
  }

  function applyHardwareProfileToGrid(state) {
    const summary = calibrationProfileSummary(state);
    const profile = summary.profile;
    state.grid.hardwareProfile = profile;
    state.grid.paperSideLengthMm = profile.sideLengthMm;
    state.grid.paperHoleToleranceMm = profile.fabricationHoleToleranceMm;
    if (summary.backlashModel !== null) state.grid.backlash = summary.backlashModel;
    if (summary.pinRadiusModel !== null) state.grid.pinRadius = summary.pinRadiusModel;
    if (summary.holeRadiusModel !== null) state.grid.holeRadius = Math.max(summary.holeRadiusModel, state.grid.pinRadius);
    return summary;
  }

  function paperRadReference(state) {
    const profile = state.grid.hardwareProfile || {};
    const sideLengthMm = Math.max(1e-9, Number(profile.sideLengthMm ?? state.grid.paperSideLengthMm ?? 35));
    const normalizedBacklash = 0.1;
    const fabricationHoleToleranceMm = Math.max(0, Number(profile.fabricationHoleToleranceMm ?? state.grid.paperHoleToleranceMm ?? 0.1));
    return {
      sideLengthMm,
      normalizedBacklash,
      poissonRatio: -0.4,
      fabricationHoleToleranceMm,
      referenceBacklashMm: normalizedBacklash * sideLengthMm,
      concentricParts: 2,
      jointsPerPart: 4,
    };
  }

  function modelLengthToMm(state, value) {
    const reference = paperRadReference(state);
    const cellSize = Math.max(1e-9, Number(state.grid.cellSize || 1));
    return (Number(value) || 0) * reference.sideLengthMm / cellSize;
  }

  function mmToModelLength(state, valueMm) {
    const reference = paperRadReference(state);
    const cellSize = Math.max(1e-9, Number(state.grid.cellSize || 1));
    return (Number(valueMm) || 0) * cellSize / reference.sideLengthMm;
  }

  function paperRadCalibration(state) {
    const reference = paperRadReference(state);
    const backlashGap = (Number(state.grid.backlash) || 0) * (Number(state.grid.cellSize) || 1);
    return {
      ...reference,
      mmPerModelUnit: reference.sideLengthMm / Math.max(1e-9, Number(state.grid.cellSize || 1)),
      configuredBacklashMm: modelLengthToMm(state, backlashGap),
      pinRadiusMm: modelLengthToMm(state, Number(state.grid.pinRadius ?? 0.18)),
      holeRadiusMm: modelLengthToMm(state, Number(state.grid.holeRadius ?? 0.225)),
      pinHoleClearanceMm: modelLengthToMm(state, pinHoleClearance(state)),
      fabricationHoleToleranceModel: mmToModelLength(state, reference.fabricationHoleToleranceMm),
    };
  }

  function verticalDeadZone(state) {
    return Math.min(commandLimits(state).z, pinHoleClearance(state));
  }

  function clampAllCommands(state) {
    const { rows, cols } = state.grid;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) {
          state.cells.commandAlpha[r][c] = 0;
          state.cells.commandZ[r][c] = 0;
          continue;
        }
        state.cells.commandAlpha[r][c] = clampCommandAlpha(state, state.cells.commandAlpha[r][c]);
        state.cells.commandZ[r][c] = clampCommandZ(state, state.cells.commandZ[r][c]);
      }
    }
  }

  function computeLinkStrain(state, centers) {
    const { rows, cols, cellSize } = state.grid;
    const horizontal = RAD.matrix(rows, Math.max(0, cols - 1), 0);
    const vertical = RAD.matrix(Math.max(0, rows - 1), cols, 0);
    let totalAbs = 0;
    let maxAbs = 0;
    let count = 0;
    const restLength = Math.max(1e-9, cellSize);

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c + 1 < cols; c += 1) {
        if (cellRemoved(state, r, c) || cellRemoved(state, r, c + 1)) continue;
        const a = centers[r][c];
        const b = centers[r][c + 1];
        const strain = Math.hypot(b.x - a.x, b.y - a.y, b.z - a.z) / restLength - 1;
        horizontal[r][c] = strain;
        totalAbs += Math.abs(strain);
        maxAbs = Math.max(maxAbs, Math.abs(strain));
        count += 1;
      }
    }
    for (let r = 0; r + 1 < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c) || cellRemoved(state, r + 1, c)) continue;
        const a = centers[r][c];
        const b = centers[r + 1][c];
        const strain = Math.hypot(b.x - a.x, b.y - a.y, b.z - a.z) / restLength - 1;
        vertical[r][c] = strain;
        totalAbs += Math.abs(strain);
        maxAbs = Math.max(maxAbs, Math.abs(strain));
        count += 1;
      }
    }
    return { horizontal, vertical, meanAbs: totalAbs / Math.max(1, count), maxAbs };
  }

  function referenceCenter(state, r, c) {
    const { rows, cols, cellSize } = state.grid;
    return {
      x: (c - (cols - 1) / 2) * cellSize,
      y: (r - (rows - 1) / 2) * cellSize,
      z: 0,
    };
  }

  function alphaPitchScale(state, alphaValue) {
    const initialAlpha = Math.max(1e-9, Number(state.grid.initialAlpha) || 1);
    return Math.max(1e-9, Number(alphaValue) || initialAlpha) / initialAlpha;
  }

  function sharedEdgePitch(state, firstAlpha, secondAlpha) {
    const cellSize = Math.max(1e-9, Number(state.grid.cellSize) || 1);
    return 0.5 * cellSize * (alphaPitchScale(state, firstAlpha) + alphaPitchScale(state, secondAlpha));
  }

  function finiteOrFallback(value, fallback) {
    if (value === null || value === undefined || value === "") return fallback;
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : fallback;
  }

  function positionLockTarget(state, r, c) {
    const reference = referenceCenter(state, r, c);
    return {
      x: finiteOrFallback(state.cells.positionLockX?.[r]?.[c], reference.x),
      y: finiteOrFallback(state.cells.positionLockY?.[r]?.[c], reference.y),
      z: finiteOrFallback(state.cells.positionLockZ?.[r]?.[c], reference.z),
    };
  }

  function computeSurfaceSlope(state, height) {
    const { rows, cols, cellSize } = state.grid;
    const magnitude = RAD.matrix(rows, cols, 0);
    const x = RAD.matrix(rows, cols, 0);
    const y = RAD.matrix(rows, cols, 0);
    const normal = RAD.matrix(rows, cols, null);
    const tilt = RAD.matrix(rows, cols, 0);
    let mean = 0;
    let max = 0;
    let meanTilt = 0;
    let maxTilt = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const left = height[r][Math.max(0, c - 1)];
        const right = height[r][Math.min(cols - 1, c + 1)];
        const down = height[Math.max(0, r - 1)][c];
        const up = height[Math.min(rows - 1, r + 1)][c];
        const dxDenom = (c === 0 || c === cols - 1 ? 1 : 2) * cellSize;
        const dyDenom = (r === 0 || r === rows - 1 ? 1 : 2) * cellSize;
        x[r][c] = (right - left) / Math.max(1e-9, dxDenom);
        y[r][c] = (up - down) / Math.max(1e-9, dyDenom);
        magnitude[r][c] = Math.hypot(x[r][c], y[r][c]);
        const normalLength = Math.hypot(x[r][c], y[r][c], 1);
        normal[r][c] = { x: -x[r][c] / normalLength, y: -y[r][c] / normalLength, z: 1 / normalLength };
        tilt[r][c] = (Math.atan(magnitude[r][c]) * 180) / Math.PI;
        mean += magnitude[r][c];
        max = Math.max(max, magnitude[r][c]);
        meanTilt += tilt[r][c];
        maxTilt = Math.max(maxTilt, tilt[r][c]);
      }
    }
    return {
      x,
      y,
      magnitude,
      normal,
      tilt,
      mean: mean / Math.max(1, rows * cols),
      max,
      meanTilt: meanTilt / Math.max(1, rows * cols),
      maxTilt,
    };
  }

  const CUSTOM_TARGET_DEFAULT = "0.35 * exp(-3 * (nr*nr + nc*nc))";
  const CUSTOM_ALLOWED_IDENTIFIERS = new Set([
    "nr",
    "nc",
    "r",
    "c",
    "d",
    "rows",
    "cols",
    "amplitude",
    "frequency",
    "pi",
    "PI",
    "Math",
    "sin",
    "cos",
    "tan",
    "asin",
    "acos",
    "atan",
    "atan2",
    "sqrt",
    "abs",
    "min",
    "max",
    "pow",
    "exp",
    "log",
    "floor",
    "ceil",
    "round",
  ]);

  function compileCustomTargetExpression(expression) {
    const source = String(expression || CUSTOM_TARGET_DEFAULT).trim();
    if (!source || source.length > 240) throw new Error("Custom target expression must be 1-240 characters.");
    if (!/^[\w\s+\-*/%.(),<>!=?:|&]*$/.test(source)) throw new Error("Custom target expression contains unsupported characters.");
    const identifiers = source.match(/[A-Za-z_]\w*/g) || [];
    for (const name of identifiers) {
      if (!CUSTOM_ALLOWED_IDENTIFIERS.has(name)) throw new Error(`Unsupported target identifier: ${name}`);
    }
    return new Function(
      "nr",
      "nc",
      "r",
      "c",
      "d",
      "rows",
      "cols",
      "amplitude",
      "frequency",
      `"use strict"; const {sin,cos,tan,asin,acos,atan,atan2,sqrt,abs,min,max,pow,exp,log,floor,ceil,round,PI}=Math; const pi=PI; return (${source});`
    );
  }

  function neighbors(rows, cols, r, c) {
    const out = [];
    if (r > 0) out.push([r - 1, c]);
    if (r + 1 < rows) out.push([r + 1, c]);
    if (c > 0) out.push([r, c - 1]);
    if (c + 1 < cols) out.push([r, c + 1]);
    return out;
  }

  function activeNeighbors(state, r, c) {
    return neighbors(state.grid.rows, state.grid.cols, r, c).filter(([nr, nc]) => cellPresent(state, nr, nc));
  }

  function topologyDiagnostics(state) {
    const { rows, cols } = state.grid;
    const labels = RAD.matrix(rows, cols, -1);
    const componentSizes = [];
    let component = 0;
    let presentCells = 0;
    let removedCells = 0;
    for (let r0 = 0; r0 < rows; r0 += 1) {
      for (let c0 = 0; c0 < cols; c0 += 1) {
        if (cellRemoved(state, r0, c0)) {
          removedCells += 1;
          continue;
        }
        presentCells += 1;
        if (labels[r0][c0] >= 0) continue;
        const queue = [[r0, c0]];
        labels[r0][c0] = component;
        let size = 0;
        while (queue.length) {
          const [r, c] = queue.shift();
          size += 1;
          for (const [nr, nc] of activeNeighbors(state, r, c)) {
            if (labels[nr][nc] >= 0) continue;
            labels[nr][nc] = component;
            queue.push([nr, nc]);
          }
        }
        componentSizes.push(size);
        component += 1;
      }
    }

    const totalEdges = rows * Math.max(0, cols - 1) + Math.max(0, rows - 1) * cols;
    let activeEdges = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (c + 1 < cols && cellPresent(state, r, c) && cellPresent(state, r, c + 1)) activeEdges += 1;
        if (r + 1 < rows && cellPresent(state, r, c) && cellPresent(state, r + 1, c)) activeEdges += 1;
      }
    }

    return {
      componentLabels: labels,
      componentCount: component,
      componentSizes,
      largestComponentSize: componentSizes.reduce((max, value) => Math.max(max, value), 0),
      presentCells,
      removedCells,
      totalEdges,
      activeEdges,
      deletedEdges: totalEdges - activeEdges,
    };
  }

  function computeInfluence(state) {
    const { rows, cols, backlash, couplingGain } = state.grid;
    const alphaDeadZone = normalizedBacklashToAlphaDeadZone(backlash);
    const influence = RAD.matrix(rows, cols, 0);
    const dieOff = RAD.matrix(rows, cols, Infinity);
    for (let sr = 0; sr < rows; sr += 1) {
      for (let sc = 0; sc < cols; sc += 1) {
        if (cellRemoved(state, sr, sc)) continue;
        const source = state.cells.commandAlpha[sr][sc];
        if (Math.abs(source) < 1e-9) continue;
        const queue = [[sr, sc, source, 0]];
        const visited = new Map();
        while (queue.length) {
          const [r, c, signal, dist] = queue.shift();
          if (cellRemoved(state, r, c)) continue;
          if (dist > rows + cols) continue;
          const key = `${r},${c}`;
          if (Math.abs(signal) <= Math.abs(visited.get(key) || 0)) continue;
          visited.set(key, signal);
          influence[r][c] += signal;
          dieOff[r][c] = Math.min(dieOff[r][c], dist);
          const next = backlashActivation(signal, alphaDeadZone) * couplingGain;
          if (Math.abs(next) < 1e-9) continue;
          for (const [nr, nc] of activeNeighbors(state, r, c)) queue.push([nr, nc, next, dist + 1]);
        }
      }
    }
    return { influence, dieOff };
  }

  function computeVerticalResidual(state) {
    const { rows, cols } = state.grid;
    const zCouplingGain = Math.max(0, Math.min(1, Number(state.grid.zCouplingGain ?? 0.32)));
    const zResidual = RAD.matrix(rows, cols, 0);
    const zDieOff = RAD.matrix(rows, cols, Infinity);
    const deadZone = verticalDeadZone(state);
    for (let sr = 0; sr < rows; sr += 1) {
      for (let sc = 0; sc < cols; sc += 1) {
        if (cellRemoved(state, sr, sc)) continue;
        const source = state.cells.commandZ[sr][sc];
        if (Math.abs(source) < 1e-9) continue;
        const queue = [[sr, sc, source, 0]];
        const visited = new Map();
        while (queue.length) {
          const [r, c, signal, dist] = queue.shift();
          if (cellRemoved(state, r, c)) continue;
          if (dist > rows + cols) continue;
          const key = `${r},${c}`;
          if (Math.abs(signal) <= Math.abs(visited.get(key) || 0)) continue;
          visited.set(key, signal);
          zResidual[r][c] += signal;
          zDieOff[r][c] = Math.min(zDieOff[r][c], dist);
          const next = backlashActivation(signal, deadZone) * zCouplingGain;
          if (Math.abs(next) < 1e-9) continue;
          for (const [nr, nc] of activeNeighbors(state, r, c)) queue.push([nr, nc, next, dist + 1]);
        }
      }
    }
    return { zResidual, zDieOff };
  }

  function propagateSingleSource(state, source, deadZone, gain) {
    const { rows, cols } = state.grid;
    const field = RAD.matrix(rows, cols, 0);
    const dieOff = RAD.matrix(rows, cols, Infinity);
    const [sr, sc, value] = source;
    if (cellRemoved(state, sr, sc)) return { field, dieOff };
    if (Math.abs(value) < 1e-9) return { field, dieOff };
    const queue = [[sr, sc, value, 0]];
    const visited = new Map();
    while (queue.length) {
      const [r, c, signal, dist] = queue.shift();
      if (cellRemoved(state, r, c)) continue;
      if (dist > rows + cols) continue;
      const key = `${r},${c}`;
      if (Math.abs(signal) <= Math.abs(visited.get(key) || 0)) continue;
      visited.set(key, signal);
      field[r][c] += signal;
      dieOff[r][c] = Math.min(dieOff[r][c], dist);
      const next = backlashActivation(signal, deadZone) * gain;
      if (Math.abs(next) < 1e-9) continue;
      for (const [nr, nc] of activeNeighbors(state, r, c)) queue.push([nr, nc, next, dist + 1]);
    }
    return { field, dieOff };
  }

  function selectedCellFootprint(state, r, c) {
    const { backlash, couplingGain } = state.grid;
    const zCouplingGain = Math.max(0, Math.min(1, Number(state.grid.zCouplingGain ?? 0.32)));
    const zDeadZone = verticalDeadZone(state);
    const alphaDeadZone = normalizedBacklashToAlphaDeadZone(backlash);
    const alpha = propagateSingleSource(state, [r, c, state.cells.commandAlpha[r][c]], alphaDeadZone, couplingGain);
    const z = propagateSingleSource(state, [r, c, state.cells.commandZ[r][c]], zDeadZone, zCouplingGain);
    return {
      alpha: alpha.field,
      alphaDieOff: alpha.dieOff,
      zResidual: z.field,
      zDieOff: z.dieOff,
    };
  }

  function selectedCellCouplingMetrics(state, r, c) {
    const { backlash, couplingGain } = state.grid;
    const zCouplingGain = Math.max(0, Math.min(1, Number(state.grid.zCouplingGain ?? 0.32)));
    const zDeadZone = verticalDeadZone(state);
    const alphaDeadZone = normalizedBacklashToAlphaDeadZone(backlash);
    if (cellRemoved(state, r, c)) {
      return {
        alphaCommand: 0,
        zCommand: 0,
        alphaDeadZone,
        zDeadZone,
        alphaCoupled: false,
        zCoupled: false,
        alphaNeighborSignal: 0,
        zNeighborSignal: 0,
        alphaReachCells: 0,
        zReachCells: 0,
        maxAlphaNeighbor: 0,
        maxZNeighbor: 0,
      };
    }
    const alphaCommand = Number(state.cells.commandAlpha[r][c]) || 0;
    const zCommand = Number(state.cells.commandZ[r][c]) || 0;
    const alphaNeighborSignal = backlashActivation(alphaCommand, alphaDeadZone) * couplingGain;
    const zNeighborSignal = backlashActivation(zCommand, zDeadZone) * zCouplingGain;
    const footprint = selectedCellFootprint(state, r, c);
    let alphaReachCells = 0;
    let zReachCells = 0;
    let maxAlphaNeighbor = 0;
    let maxZNeighbor = 0;
    for (let rr = 0; rr < state.grid.rows; rr += 1) {
      for (let cc = 0; cc < state.grid.cols; cc += 1) {
        if (rr === r && cc === c) continue;
        const alpha = Math.abs(footprint.alpha[rr][cc] || 0);
        const z = Math.abs(footprint.zResidual[rr][cc] || 0);
        if (alpha > 1e-9) alphaReachCells += 1;
        if (z > 1e-9) zReachCells += 1;
        maxAlphaNeighbor = Math.max(maxAlphaNeighbor, alpha);
        maxZNeighbor = Math.max(maxZNeighbor, z);
      }
    }
    return {
      alphaCommand,
      zCommand,
      alphaDeadZone,
      zDeadZone,
      alphaCoupled: Math.abs(alphaNeighborSignal) > 1e-9,
      zCoupled: Math.abs(zNeighborSignal) > 1e-9,
      alphaNeighborSignal,
      zNeighborSignal,
      alphaReachCells,
      zReachCells,
      maxAlphaNeighbor,
      maxZNeighbor,
    };
  }

  function brushCells(state, r, c, radiusValue = state.view?.paintRadius || 0) {
    const radius = Math.max(0, Math.min(2, Math.round(Number(radiusValue || 0))));
    const cells = [];
    for (let rr = Math.max(0, r - radius); rr <= Math.min(state.grid.rows - 1, r + radius); rr += 1) {
      for (let cc = Math.max(0, c - radius); cc <= Math.min(state.grid.cols - 1, c + radius); cc += 1) {
        if (cellRemoved(state, rr, cc)) continue;
        if (Math.abs(rr - r) + Math.abs(cc - c) <= radius) cells.push({ r: rr, c: cc });
      }
    }
    return cells;
  }

  function computeCenters(state, alpha = null) {
    const { rows, cols, cellSize } = state.grid;
    const centers = RAD.matrix(rows, cols, null);
    const xRows = RAD.matrix(rows, cols, 0);
    const yCols = RAD.matrix(rows, cols, 0);
    const alphaGrid = alpha || RAD.matrix(rows, cols, Number(state.grid.initialAlpha) || 1);
    for (let r = 0; r < rows; r += 1) {
      for (let c = 1; c < cols; c += 1) {
        xRows[r][c] = xRows[r][c - 1] + sharedEdgePitch(state, alphaGrid[r][c - 1], alphaGrid[r][c]);
      }
    }
    for (let c = 0; c < cols; c += 1) {
      for (let r = 1; r < rows; r += 1) {
        yCols[r][c] = yCols[r - 1][c] + sharedEdgePitch(state, alphaGrid[r - 1][c], alphaGrid[r][c]);
      }
    }
    let mx = 0;
    let my = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        mx += xRows[r][c];
        my += yCols[r][c];
      }
    }
    mx /= Math.max(1, rows * cols);
    my /= Math.max(1, rows * cols);
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        centers[r][c] = { x: xRows[r][c] - mx, y: yCols[r][c] - my, z: 0 };
      }
    }
    return centers;
  }

  function pairDistanceFromTheta(state, thetaA, thetaB) {
    const cellSize = Math.max(1e-9, Number(state.grid.cellSize) || 1);
    return cellSize * (Math.cos((thetaA * Math.PI) / 180) + Math.cos((thetaB * Math.PI) / 180));
  }

  function solveConstraintKinematicApprox(state, alpha, theta, preferredCenters, baseRealization) {
    const { rows, cols } = state.grid;
    const centers = RAD.matrix(rows, cols, (r, c) => ({ ...baseRealization.centers[r][c] }));
    const compressionResidual = RAD.matrix(rows, cols, 0);
    const inducedHeight = RAD.matrix(rows, cols, 0);
    const edgeCounts = RAD.matrix(rows, cols, 0);
    let maxEdgeError = 0;
    let meanEdgeError = 0;
    let edgeCount = 0;
    let contactCount = baseRealization.metrics?.constrainedCellCount || 0;
    const boundary = state.boundary || {};
    const zGain = Math.max(0, Number(boundary.zGain) || 1.15);
    const threshold = Math.max(0, Number(boundary.zThreshold) || 0);
    const power = Math.max(0.5, Number(boundary.zPower) || 1.5);

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        for (const [nr, nc] of activeNeighbors(state, r, c)) {
          if (nr < r || (nr === r && nc < c)) continue;
          const a = centers[r][c];
          const b = centers[nr][nc];
          const desired = Math.max(0, pairDistanceFromTheta(state, theta[r][c], theta[nr][nc]));
          const xy = Math.hypot(b.x - a.x, b.y - a.y);
          const compression = Math.max(0, desired - xy);
          compressionResidual[r][c] += compression;
          compressionResidual[nr][nc] += compression;
          edgeCounts[r][c] += 1;
          edgeCounts[nr][nc] += 1;
          maxEdgeError = Math.max(maxEdgeError, compression);
          meanEdgeError += compression;
          edgeCount += 1;
        }
      }
    }

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        const residual = compressionResidual[r][c] / Math.max(1, edgeCounts[r][c]);
        compressionResidual[r][c] = Math.max(residual, baseRealization.compressionResidual?.[r]?.[c] || 0);
        const positionTarget = state.cells.positionLocked?.[r]?.[c] === true ? positionLockTarget(state, r, c) : null;
        if (positionTarget) {
          inducedHeight[r][c] = 0;
          centers[r][c].z = positionTarget.z;
        } else {
          inducedHeight[r][c] = zGain * Math.pow(Math.max(0, compressionResidual[r][c] - threshold), power);
          centers[r][c].z = inducedHeight[r][c];
        }
      }
    }

    return {
      ...baseRealization,
      schema: "rad-sim.constraint-realization.v1",
      centers,
      preferredCenters,
      compressionResidual,
      inducedHeight,
      active: true,
      metrics: {
        ...(baseRealization.metrics || {}),
        constrainedCellCount: contactCount,
        meanCompressionResidual: meanEdgeError / Math.max(1, edgeCount),
        maxCompressionResidual: maxEdgeError,
        maxInducedHeight: Math.max(...inducedHeight.flat().map((value) => Math.abs(value))),
        constraintSolved: true,
        constraintMaxEdgeError: maxEdgeError,
        constraintMeanEdgeError: meanEdgeError / Math.max(1, edgeCount),
        constraintContactCount: contactCount,
      },
    };
  }

  function targetSurface(state) {
    const { rows, cols } = state.grid;
    const { type, amplitude, frequency } = state.target;
    const midR = (rows - 1) / 2;
    const midC = (cols - 1) / 2;
    const radius = Math.max(1, Math.min(rows, cols) / 2);
    const target = RAD.matrix(rows, cols, 0);
    let customTarget = null;
    state.target.expressionError = "";
    if (type === "custom") {
      try {
        customTarget = compileCustomTargetExpression(state.target.customExpression);
      } catch (error) {
        state.target.expressionError = error.message;
      }
    }
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const nr = rows === 1 ? 0 : (r - midR) / Math.max(1, midR);
        const nc = cols === 1 ? 0 : (c - midC) / Math.max(1, midC);
        const d = Math.hypot(r - midR, c - midC);
        if (type === "dome") target[r][c] = amplitude * Math.max(0, 1 - d / radius);
        else if (type === "saddle") target[r][c] = amplitude * nr * nc;
        else if (type === "wave") target[r][c] = amplitude * Math.sin(nc * Math.PI * frequency);
        else if (type === "gaussian") target[r][c] = amplitude * Math.exp(-(nr * nr + nc * nc) * 2.4);
        else if (type === "tilt") target[r][c] = amplitude * (0.65 * nc + 0.35 * nr);
        else if (type === "custom" && customTarget) {
          try {
            const value = Number(customTarget(nr, nc, r, c, d, rows, cols, amplitude, frequency));
            target[r][c] = Number.isFinite(value) ? clamp(value, -1.4, 1.4) : 0;
          } catch (error) {
            state.target.expressionError = error.message;
            target[r][c] = 0;
          }
        }
      }
    }
    return target;
  }

  function simulate(state) {
    const { rows, cols, initialAlpha, alphaMin, alphaMax } = state.grid;
    const topology = topologyDiagnostics(state);
    const { influence, dieOff } = computeInfluence(state);
    const { zResidual, zDieOff } = computeVerticalResidual(state);
    const alpha = RAD.matrix(rows, cols, 1);
    const theta = RAD.matrix(rows, cols, 0);
    const height = RAD.matrix(rows, cols, 0);
    let meanAlpha = 0;
    let maxAbsHeight = 0;
    let activeCells = 0;
    let maxFiniteDieOff = 0;
    let rmsTargetError = 0;
    let maxTargetError = 0;
    let meanSignedTargetError = 0;
    let maxPositiveTargetError = 0;
    let maxNegativeTargetError = 0;
    let meanTravel = 0;
    let meanSaturation = 0;
    let maxSaturation = 0;
    let saturatedActuators = 0;
    const target = targetSurface(state);
    const targetError = RAD.matrix(rows, cols, 0);
    const recommended = RAD.matrix(rows, cols, false);
    const saturation = RAD.matrix(rows, cols, 0);

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) {
          alpha[r][c] = initialAlpha;
          theta[r][c] = alphaToTheta(initialAlpha);
          height[r][c] = 0;
          targetError[r][c] = 0;
          recommended[r][c] = false;
          saturation[r][c] = 0;
          meanAlpha += alpha[r][c];
          continue;
        }
        const locked = state.cells.locked[r][c];
        const positionLocked = state.cells.positionLocked?.[r]?.[c] === true;
        const rawAlpha = initialAlpha + influence[r][c];
        const lockAlpha = Number(state.cells.lockAlpha?.[r]?.[c]) || initialAlpha;
        const rawLockZ = Number(state.cells.lockZ?.[r]?.[c]);
        const lockZ = Number.isFinite(rawLockZ) ? rawLockZ : 0;
        alpha[r][c] = locked ? Math.max(alphaMin, Math.min(alphaMax, lockAlpha)) : Math.max(alphaMin, Math.min(alphaMax, rawAlpha));
        theta[r][c] = alphaToTheta(alpha[r][c]);
        const positionTarget = positionLocked ? positionLockTarget(state, r, c) : null;
        height[r][c] = positionTarget ? positionTarget.z : locked ? lockZ : zResidual[r][c];
        targetError[r][c] = height[r][c] - target[r][c];
        rmsTargetError += targetError[r][c] ** 2;
        maxTargetError = Math.max(maxTargetError, Math.abs(targetError[r][c]));
        meanSignedTargetError += targetError[r][c];
        maxPositiveTargetError = Math.max(maxPositiveTargetError, targetError[r][c]);
        maxNegativeTargetError = Math.min(maxNegativeTargetError, targetError[r][c]);
        meanAlpha += alpha[r][c];
        maxAbsHeight = Math.max(maxAbsHeight, Math.abs(height[r][c]));
        const travel = Math.abs(state.cells.commandAlpha[r][c]) + Math.abs(state.cells.commandZ[r][c]);
        meanTravel += travel;
        saturation[r][c] = commandSaturation(state, state.cells.commandAlpha[r][c], state.cells.commandZ[r][c]);
        meanSaturation += saturation[r][c];
        maxSaturation = Math.max(maxSaturation, saturation[r][c]);
        if (saturation[r][c] >= 0.98 && travel > 1e-9) saturatedActuators += 1;
        recommended[r][c] = travel > 1e-9;
        if (travel > 1e-9 || locked) {
          activeCells += 1;
        }
        if (Number.isFinite(dieOff[r][c])) maxFiniteDieOff = Math.max(maxFiniteDieOff, dieOff[r][c]);
      }
    }

    const preferredCenters = computeCenters(state, alpha);
    let constraintRealization =
      typeof RAD.solveBoundaryConstraints === "function"
        ? RAD.solveBoundaryConstraints(state, alpha, preferredCenters)
        : {
            preferredCenters,
            centers: computeCenters(state, alpha),
            compressionResidual: RAD.matrix(rows, cols, 0),
            inducedHeight: RAD.matrix(rows, cols, 0),
            constraintDisplacement: RAD.matrix(rows, cols, 0),
            wallContact: RAD.matrix(rows, cols, false),
            positionLockContact: RAD.matrix(rows, cols, false),
            active: false,
            metrics: { constrainedCellCount: 0, meanCompressionResidual: 0, maxCompressionResidual: 0, maxInducedHeight: 0 },
          };
    if (state.view?.simulationMode === "constraintSolved") {
      constraintRealization = solveConstraintKinematicApprox(state, alpha, theta, preferredCenters, constraintRealization);
    }
    const centers = constraintRealization.centers;
    rmsTargetError = 0;
    maxTargetError = 0;
    meanSignedTargetError = 0;
    maxPositiveTargetError = 0;
    maxNegativeTargetError = 0;
    maxAbsHeight = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) {
          height[r][c] = 0;
          targetError[r][c] = 0;
          centers[r][c].z = 0;
          continue;
        }
        const locked = state.cells.locked[r][c];
        const positionTarget = state.cells.positionLocked?.[r]?.[c] === true ? positionLockTarget(state, r, c) : null;
        const rawLockZ = Number(state.cells.lockZ?.[r]?.[c]);
        const lockZ = Number.isFinite(rawLockZ) ? rawLockZ : 0;
        height[r][c] = positionTarget
          ? positionTarget.z
          : locked
            ? lockZ + (constraintRealization.inducedHeight?.[r]?.[c] || 0)
            : zResidual[r][c] + (constraintRealization.inducedHeight?.[r]?.[c] || 0);
        targetError[r][c] = height[r][c] - target[r][c];
        rmsTargetError += targetError[r][c] ** 2;
        maxTargetError = Math.max(maxTargetError, Math.abs(targetError[r][c]));
        meanSignedTargetError += targetError[r][c];
        maxPositiveTargetError = Math.max(maxPositiveTargetError, targetError[r][c]);
        maxNegativeTargetError = Math.min(maxNegativeTargetError, targetError[r][c]);
        maxAbsHeight = Math.max(maxAbsHeight, Math.abs(height[r][c]));
        centers[r][c].z = height[r][c];
      }
    }
    const slope = computeSurfaceSlope(state, height);
    const linkStrain = computeLinkStrain(state, centers);
    const displacement = RAD.matrix(rows, cols, 0);
    let meanReferenceDisplacement = 0;
    let maxReferenceDisplacement = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const ref = referenceCenter(state, r, c);
        const current = centers[r][c];
        displacement[r][c] = Math.hypot(current.x - ref.x, current.y - ref.y, current.z - ref.z);
        meanReferenceDisplacement += displacement[r][c];
        maxReferenceDisplacement = Math.max(maxReferenceDisplacement, displacement[r][c]);
      }
    }

    return {
      alpha,
      theta,
      height,
      target,
      targetError,
      recommended,
      saturation,
      slope,
      linkStrain,
      displacement,
      centers,
      preferredCenters,
      constraintRealization,
      compressionResidual: constraintRealization.compressionResidual,
      inducedHeight: constraintRealization.inducedHeight,
      constraintDisplacement: constraintRealization.constraintDisplacement,
      wallContact: constraintRealization.wallContact,
      positionLockContact: constraintRealization.positionLockContact,
      influence,
      dieOff,
      zResidual,
      zDieOff,
      topology,
      metrics: {
        meanAlpha: meanAlpha / (rows * cols),
        maxAbsHeight,
        activeCells,
        maxFiniteDieOff,
        rmsTargetError: Math.sqrt(rmsTargetError / (rows * cols)),
        maxTargetError,
        meanSignedTargetError: meanSignedTargetError / (rows * cols),
        maxPositiveTargetError,
        maxNegativeTargetError,
        meanTravel: meanTravel / (rows * cols),
        meanSaturation: meanSaturation / (rows * cols),
        maxSaturation,
        saturatedActuators,
        removedCells: topology.removedCells,
        topologyComponentCount: topology.componentCount,
        deletedEdges: topology.deletedEdges,
        meanAbsLinkStrain: linkStrain.meanAbs,
        maxAbsLinkStrain: linkStrain.maxAbs,
        constrainedCellCount: constraintRealization.metrics?.constrainedCellCount || 0,
        meanCompressionResidual: constraintRealization.metrics?.meanCompressionResidual || 0,
        maxCompressionResidual: constraintRealization.metrics?.maxCompressionResidual || 0,
        maxInducedHeight: constraintRealization.metrics?.maxInducedHeight || 0,
        boundaryConstraintsActive: constraintRealization.active === true,
        constraintSolved: constraintRealization.metrics?.constraintSolved === true,
        constraintMaxEdgeError: constraintRealization.metrics?.constraintMaxEdgeError || 0,
        constraintMeanEdgeError: constraintRealization.metrics?.constraintMeanEdgeError || 0,
        constraintContactCount: constraintRealization.metrics?.constraintContactCount || 0,
        meanReferenceDisplacement: meanReferenceDisplacement / (rows * cols),
        maxReferenceDisplacement,
        meanSurfaceSlope: slope.mean,
        maxSurfaceSlope: slope.max,
        meanNormalTilt: slope.meanTilt,
        maxNormalTilt: slope.maxTilt,
        recommendedActuators: recommended.flat().filter(Boolean).length,
        minTheta: Math.min(...theta.flat()),
        maxTheta: Math.max(...theta.flat()),
        minAlpha: Math.min(...alpha.flat()),
        maxAlpha: Math.max(...alpha.flat()),
      },
    };
  }

  function objectiveForSimulation(sim, state, options = {}) {
    const actuatorPenalty = options.actuatorPenalty ?? 0.012;
    const travelPenalty = options.travelPenalty ?? 0.018;
    const { rows, cols } = state.grid;
    let active = 0;
    let travel = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        const amount = Math.abs(state.cells.commandAlpha[r][c]) + Math.abs(state.cells.commandZ[r][c]);
        if (amount > 1e-9) active += 1;
        travel += amount;
      }
    }
    const n = rows * cols;
    return sim.metrics.rmsTargetError + actuatorPenalty * (active / n) + travelPenalty * (travel / n);
  }

  function estimateCommandsForTarget(state, gain = 0.82) {
    if (!state.experiment.initialSnapshot) state.experiment.initialSnapshot = RAD.snapshotState(state);
    const target = targetSurface(state);
    const { rows, cols } = state.grid;
    RAD.clearCommands(state);
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        if (state.cells.actuatorAllowed?.[r]?.[c] === false) continue;
        state.cells.commandZ[r][c] = clampCommandZ(state, target[r][c] * gain);
        state.cells.commandAlpha[r][c] = clampCommandAlpha(state, -0.22 * Math.abs(target[r][c]));
      }
    }
    RAD.recordEvent(state, { type: "target-fit-seed", target: state.target.type });
  }

  function optimizeCommandsForTarget(state, options = {}) {
    if (!state.experiment.initialSnapshot) state.experiment.initialSnapshot = RAD.snapshotState(state);
    const iterations = options.iterations ?? state.target.optimizerIterations ?? 48;
    const actuatorPenalty = options.actuatorPenalty ?? state.target.actuatorPenalty ?? 0.012;
    const travelPenalty = options.travelPenalty ?? state.target.travelPenalty ?? 0.018;
    const target = targetSurface(state);
    const { rows, cols } = state.grid;
    RAD.clearCommands(state);
    let sim = simulate(state);
    let bestScore = objectiveForSimulation(sim, state, { actuatorPenalty, travelPenalty });
    const candidates = [];
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        if (state.cells.actuatorAllowed?.[r]?.[c] !== false) candidates.push([r, c]);
      }
    }

    for (let iter = 0; iter < iterations; iter += 1) {
      sim = simulate(state);
      candidates.sort((a, b) => {
        const ea = Math.abs(target[a[0]][a[1]] - sim.height[a[0]][a[1]]);
        const eb = Math.abs(target[b[0]][b[1]] - sim.height[b[0]][b[1]]);
        return eb - ea;
      });
      let accepted = false;
      for (const [r, c] of candidates.slice(0, Math.min(10, candidates.length))) {
        const beforeZ = state.cells.commandZ[r][c];
        const beforeAlpha = state.cells.commandAlpha[r][c];
        const residual = target[r][c] - sim.height[r][c];
        const zStep = clampCommandZ(state, beforeZ + residual * 0.75);
        const alphaStep = clampCommandAlpha(state, beforeAlpha - Math.sign(residual || target[r][c]) * Math.min(0.32, Math.abs(residual) * 0.18));
        state.cells.commandZ[r][c] = zStep;
        state.cells.commandAlpha[r][c] = alphaStep;
        const trial = simulate(state);
        const score = objectiveForSimulation(trial, state, { actuatorPenalty, travelPenalty });
        if (score + 1e-6 < bestScore) {
          bestScore = score;
          accepted = true;
          break;
        }
        state.cells.commandZ[r][c] = beforeZ;
        state.cells.commandAlpha[r][c] = beforeAlpha;
      }
      if (!accepted) break;
    }

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (Math.abs(state.cells.commandZ[r][c]) < 0.025) state.cells.commandZ[r][c] = 0;
        if (Math.abs(state.cells.commandAlpha[r][c]) < 0.025) state.cells.commandAlpha[r][c] = 0;
      }
    }
    RAD.recordEvent(state, {
      type: "target-fit-optimized",
      target: state.target.type,
      iterations,
      actuatorPenalty,
      travelPenalty,
      score: objectiveForSimulation(simulate(state), state, { actuatorPenalty, travelPenalty }),
    });
  }

  function applyPreset(state, name) {
    if (!state.experiment.initialSnapshot) state.experiment.initialSnapshot = RAD.snapshotState(state);
    RAD.clearCommands(state);
    state.target.type = "none";
    const { rows, cols } = state.grid;
    const midR = Math.floor(rows / 2);
    const midC = Math.floor(cols / 2);
    const radius = Math.max(1, Math.min(rows, cols) / 2.6);
    if (name === "center") {
      state.cells.commandAlpha[midR][midC] = -0.38;
    } else if (name === "dome") {
      state.target.type = "dome";
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          const d = Math.hypot(r - midR, c - midC);
          state.cells.commandZ[r][c] = Math.max(0, 0.55 * (1 - d / radius));
          state.cells.commandAlpha[r][c] = -0.12 * Math.max(0, 1 - d / radius);
        }
      }
    } else if (name === "saddle") {
      state.target.type = "saddle";
      state.cells.commandZ[0][0] = 0.45;
      state.cells.commandZ[rows - 1][cols - 1] = 0.45;
      state.cells.commandZ[0][cols - 1] = -0.45;
      state.cells.commandZ[rows - 1][0] = -0.45;
      state.cells.commandAlpha[0][cols - 1] = -0.25;
      state.cells.commandAlpha[rows - 1][0] = -0.25;
    } else if (name === "ridge") {
      for (let r = 0; r < rows; r += 1) {
        state.cells.commandAlpha[r][midC] = -0.22;
        state.cells.commandZ[r][midC] = 0.32;
      }
      state.cells.locked[0][0] = true;
      state.cells.locked[rows - 1][cols - 1] = true;
    } else if (name === "wave") {
      state.target.type = "wave";
      for (let c = 0; c < cols; c += 1) {
        const z = 0.36 * Math.sin((c / Math.max(1, cols - 1)) * Math.PI * 2);
        for (let r = 0; r < rows; r += 1) state.cells.commandZ[r][c] = z;
      }
    } else if (name === "corner") {
      state.cells.commandZ[rows - 1][cols - 1] = 0.65;
      state.cells.commandAlpha[rows - 1][cols - 1] = -0.3;
      state.cells.locked[0][0] = true;
    } else if (name === "twist") {
      for (let i = 0; i < Math.min(rows, cols); i += 1) {
        state.cells.commandZ[i][i] = 0.4 * (i / Math.max(1, Math.min(rows, cols) - 1) - 0.5);
        state.cells.commandAlpha[i][cols - 1 - i] = i % 2 ? 0.18 : -0.18;
      }
    } else if (name === "ring") {
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          const d = Math.hypot(r - midR, c - midC);
          if (Math.abs(d - radius * 0.65) < 0.75) {
            state.cells.commandAlpha[r][c] = -0.2;
            state.cells.commandZ[r][c] = 0.22;
          }
        }
      }
    } else if (name === "checker") {
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          state.cells.commandAlpha[r][c] = (r + c) % 2 === 0 ? -0.16 : 0.16;
          state.cells.commandZ[r][c] = (r + c) % 2 === 0 ? 0.18 : -0.18;
        }
      }
    }
    state.experiment.presetName = name;
    RAD.recordEvent(state, { type: "preset", name });
  }

  RAD.backlashActivation = backlashActivation;
  RAD.alphaToTheta = alphaToTheta;
  RAD.normalizedBacklashToThetaDeadZone = normalizedBacklashToThetaDeadZone;
  RAD.normalizedBacklashToAlphaDeadZone = normalizedBacklashToAlphaDeadZone;
  RAD.commandLimits = commandLimits;
  RAD.clampCommandZ = clampCommandZ;
  RAD.clampCommandAlpha = clampCommandAlpha;
  RAD.cellRemoved = cellRemoved;
  RAD.cellPresent = cellPresent;
  RAD.commandSaturation = commandSaturation;
  RAD.pinHoleClearance = pinHoleClearance;
  RAD.alphaPitchScale = alphaPitchScale;
  RAD.sharedEdgePitch = sharedEdgePitch;
  RAD.HARDWARE_PROFILE_DIMENSIONS = HARDWARE_PROFILE_DIMENSIONS;
  RAD.HARDWARE_PROFILE_SCHEMA = HARDWARE_PROFILE_SCHEMA;
  RAD.CALIBRATION_VISUAL_FIELDS = CALIBRATION_VISUAL_FIELDS;
  RAD.CALIBRATION_MESH_FIELDS = CALIBRATION_MESH_FIELDS;
  RAD.CALIBRATION_SOLVER_GAPS = CALIBRATION_SOLVER_GAPS;
  RAD.hardwareProfile = hardwareProfile;
  RAD.hardwareProfileSchemaPayload = hardwareProfileSchemaPayload;
  RAD.hardwareProfileFromObject = hardwareProfileFromObject;
  RAD.importHardwareProfileJson = importHardwareProfileJson;
  RAD.exportHardwareProfileJson = exportHardwareProfileJson;
  RAD.hardwareMeasuredFields = hardwareMeasuredFields;
  RAD.hardwareMissingFields = hardwareMissingFields;
  RAD.calibrationProfileSummary = calibrationProfileSummary;
  RAD.calibrationReadiness = calibrationReadiness;
  RAD.calibrationMeasurementPlan = calibrationMeasurementPlan;
  RAD.exportCalibrationMeasurementPlan = exportCalibrationMeasurementPlan;
  RAD.applyHardwareProfileToGrid = applyHardwareProfileToGrid;
  RAD.paperRadReference = paperRadReference;
  RAD.modelLengthToMm = modelLengthToMm;
  RAD.mmToModelLength = mmToModelLength;
  RAD.paperRadCalibration = paperRadCalibration;
  RAD.verticalDeadZone = verticalDeadZone;
  RAD.clampAllCommands = clampAllCommands;
  RAD.computeLinkStrain = computeLinkStrain;
  RAD.computeCenters = computeCenters;
  RAD.pairDistanceFromTheta = pairDistanceFromTheta;
  RAD.solveConstraintKinematicApprox = solveConstraintKinematicApprox;
  RAD.referenceCenter = referenceCenter;
  RAD.positionLockTarget = positionLockTarget;
  RAD.computeSurfaceSlope = computeSurfaceSlope;
  RAD.compileCustomTargetExpression = compileCustomTargetExpression;
  RAD.targetSurface = targetSurface;
  RAD.topologyDiagnostics = topologyDiagnostics;
  RAD.selectedCellFootprint = selectedCellFootprint;
  RAD.selectedCellCouplingMetrics = selectedCellCouplingMetrics;
  RAD.brushCells = brushCells;
  RAD.simulate = simulate;
  RAD.estimateCommandsForTarget = estimateCommandsForTarget;
  RAD.optimizeCommandsForTarget = optimizeCommandsForTarget;
  RAD.applyPreset = applyPreset;
})();
