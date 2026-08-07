(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function matrix(rows, cols, fill) {
    return Array.from({ length: rows }, (_, r) =>
      Array.from({ length: cols }, (_, c) => (typeof fill === "function" ? fill(r, c) : fill))
    );
  }

  function defaultThetaForAlpha(alpha) {
    return 70 * alpha - 60;
  }

  function derivedCellFields(rows, cols, alpha = 1) {
    return {
      alpha: matrix(rows, cols, alpha),
      theta: matrix(rows, cols, defaultThetaForAlpha(alpha)),
      z: matrix(rows, cols, 0),
    };
  }

  function defaultHardwareProfile() {
    return {
      name: "paper-reference",
      source: "RAD preprint defaults",
      sideLengthMm: 35,
      fabricationHoleToleranceMm: 0.1,
      backlashMm: null,
      pinRadiusMm: null,
      holeRadiusMm: null,
      plateThicknessMm: null,
      jointStackHeightMm: null,
      bossRadiusMm: null,
      notes: "",
    };
  }

  function updateDerivedCells(state, sim) {
    if (!state?.cells || !sim) return state;
    const { rows, cols } = state.grid;
    state.cells.alpha = matrix(rows, cols, (r, c) => Number(sim.alpha?.[r]?.[c]) || 0);
    state.cells.theta = matrix(rows, cols, (r, c) => Number(sim.theta?.[r]?.[c]) || 0);
    state.cells.z = matrix(rows, cols, (r, c) => Number(sim.height?.[r]?.[c]) || 0);
    return state;
  }

  function ensureCellSchema(state) {
    const { rows, cols } = state.grid;
    const initialAlpha = Number(state.grid.initialAlpha) || 1;
    const defaults = derivedCellFields(rows, cols, initialAlpha);
    if (!state.cells.alpha) state.cells.alpha = defaults.alpha;
    if (!state.cells.theta) state.cells.theta = defaults.theta;
    if (!state.cells.z) state.cells.z = defaults.z;
    if (!state.cells.commandAlpha) state.cells.commandAlpha = matrix(rows, cols, 0);
    if (!state.cells.commandZ) state.cells.commandZ = matrix(rows, cols, 0);
    if (!state.cells.locked) state.cells.locked = matrix(rows, cols, false);
    if (!state.cells.lockAlpha) state.cells.lockAlpha = matrix(rows, cols, initialAlpha);
    if (!state.cells.actuatorAllowed) state.cells.actuatorAllowed = matrix(rows, cols, true);
    return state;
  }

  function ensureGridSchema(state) {
    if (state.grid.zCouplingGain === undefined) state.grid.zCouplingGain = 0.32;
    if (state.grid.pinRadius === undefined) state.grid.pinRadius = 0.18;
    if (state.grid.holeRadius === undefined) state.grid.holeRadius = 0.225;
    if (state.grid.paperSideLengthMm === undefined) state.grid.paperSideLengthMm = 35;
    if (state.grid.paperHoleToleranceMm === undefined) state.grid.paperHoleToleranceMm = 0.1;
    const rawHardwareProfile = state.grid.hardwareProfile || {};
    state.grid.hardwareProfile = {
      ...defaultHardwareProfile(),
      ...rawHardwareProfile,
    };
    if (rawHardwareProfile.sideLengthMm === undefined || rawHardwareProfile.sideLengthMm === null) {
      state.grid.hardwareProfile.sideLengthMm = state.grid.paperSideLengthMm;
    }
    if (
      rawHardwareProfile.fabricationHoleToleranceMm === undefined ||
      rawHardwareProfile.fabricationHoleToleranceMm === null
    ) {
      state.grid.hardwareProfile.fabricationHoleToleranceMm = state.grid.paperHoleToleranceMm;
    }
    if (!Number.isFinite(Number(state.grid.hardwareProfile.sideLengthMm))) {
      state.grid.hardwareProfile.sideLengthMm = state.grid.paperSideLengthMm;
    }
    if (!Number.isFinite(Number(state.grid.hardwareProfile.fabricationHoleToleranceMm))) {
      state.grid.hardwareProfile.fabricationHoleToleranceMm = state.grid.paperHoleToleranceMm;
    }
    state.grid.hardwareProfile.sideLengthMm = Math.max(1e-9, Number(state.grid.hardwareProfile.sideLengthMm));
    state.grid.hardwareProfile.fabricationHoleToleranceMm = Math.max(0, Number(state.grid.hardwareProfile.fabricationHoleToleranceMm));
    state.grid.paperSideLengthMm = state.grid.hardwareProfile.sideLengthMm;
    state.grid.paperHoleToleranceMm = state.grid.hardwareProfile.fabricationHoleToleranceMm;
    if (state.grid.holeRadius < state.grid.pinRadius) state.grid.holeRadius = state.grid.pinRadius;
    return state;
  }

  function stateWithFreshDerivedCells(state) {
    const copy = cloneData(state);
    ensureCellSchema(copy);
    if (typeof RAD.simulate === "function") updateDerivedCells(copy, RAD.simulate(copy));
    return copy;
  }

  function createState(rows = 7, cols = 7) {
    const midR = Math.floor(rows / 2);
    const midC = Math.floor(cols / 2);
    const state = {
      grid: {
        rows,
        cols,
        cellSize: 1,
        backlash: 0.1,
        couplingGain: 0.55,
        zCouplingGain: 0.32,
        pinRadius: 0.18,
        holeRadius: 0.225,
        paperSideLengthMm: 35,
        paperHoleToleranceMm: 0.1,
        hardwareProfile: defaultHardwareProfile(),
        initialAlpha: 1,
        alphaMin: 0.25,
        alphaMax: 1.75,
        zTravelLimit: 0.8,
        alphaContractLimit: 0.55,
        alphaExpandLimit: 0.35,
      },
      cells: {
        ...derivedCellFields(rows, cols, 1),
        commandAlpha: matrix(rows, cols, 0),
        commandZ: matrix(rows, cols, 0),
        locked: matrix(rows, cols, false),
        lockAlpha: matrix(rows, cols, 1),
        actuatorAllowed: matrix(rows, cols, true),
      },
      view: {
        cellVisualMode: "abstract",
        simulationMode: "kinematic",
        isolateSelected: false,
        explodedSelected: false,
        quickDockCollapsed: false,
        paintMode: false,
        paintRadius: 0,
        overlayMode: "alpha",
        membraneVisible: true,
        referenceVisible: true,
        displacementVectorsVisible: false,
        gapsVisible: true,
        stopsVisible: false,
        pivotsVisible: false,
        linkagesVisible: false,
        fastenersVisible: false,
        actuatorsVisible: false,
        actuatorDisplayMode: "selected",
        influenceFootprintVisible: true,
        measurementsVisible: true,
        measurementLabelsVisible: false,
        measurementMode: "alpha",
        targetVisible: false,
        targetErrorVectorsVisible: true,
        surfaceNormalsVisible: false,
        surfaceInterpolation: "smooth",
        surfaceSubdivisions: 5,
        surfaceContoursVisible: false,
        contourMode: "membrane",
        animationResponse: 0.18,
        camera: {
          mode: "iso",
          projection: "perspective",
          radius: 15.6,
          theta: -Math.PI / 4,
          phi: Math.acos(1 / Math.sqrt(3)),
          target: { x: 0, y: 0, z: 0 },
        },
      },
      target: {
        type: "none",
        amplitude: 0.45,
        frequency: 1,
        customExpression: "0.35 * exp(-3 * (nr*nr + nc*nc))",
        expressionError: "",
        optimizerIterations: 48,
        actuatorPenalty: 0.012,
        travelPenalty: 0.018,
      },
      inverse: {
        maxActuators: 24,
        lastScore: 0,
        plan: { candidates: [], commands: [], history: [] },
        preview: null,
        sensitivity: { candidates: [], map: matrix(rows, cols, 0), stepZ: 0.12, stepAlpha: 0.12, controllableCells: 0, meanGain: 0, maxGain: 0 },
        jacobian: { columns: [], coverageMap: matrix(rows, cols, 0), actuatorCount: 0, columnCount: 0, meanCoverage: 0, maxCoverage: 0, stepZ: 0.12, stepAlpha: 0.12, conditionEstimate: 0 },
        linearSolution: { commands: [], history: [], steps: 0, baseError: 0, predictedError: 0, projectedError: 0, projectedActuators: 0 },
        physicalValidation: null,
      },
      timeline: {
        playing: false,
        index: 0,
        speed: 1,
        smooth: true,
        transitionMs: 900,
        keyframeName: "pose",
      },
      selection: { r: midR, c: midC },
      experiment: {
        presetName: "center",
        eventList: [],
        initialSnapshot: null,
        characterizationScope: "single",
        characterization: null,
        notes: "",
      },
    };
    state.cells.commandAlpha[midR][midC] = -0.35;
    return state;
  }

  function resizeState(state, rows, cols) {
    const old = state.cells;
    const next = createState(rows, cols);
    next.grid = { ...state.grid, rows, cols };
    next.view = { ...state.view };
    next.target = { ...state.target };
    next.inverse = { ...state.inverse, plan: { candidates: [], commands: [], history: [] }, jacobian: { columns: [], coverageMap: matrix(rows, cols, 0), actuatorCount: 0, columnCount: 0, meanCoverage: 0, maxCoverage: 0, stepZ: 0.12, stepAlpha: 0.12, conditionEstimate: 0 }, physicalValidation: null };
    next.timeline = { ...state.timeline, index: Math.min(state.timeline.index, next.experiment.eventList.length) };
    next.experiment = { ...state.experiment };
    next.selection = {
      r: Math.min(state.selection.r, rows - 1),
      c: Math.min(state.selection.c, cols - 1),
    };
    for (let r = 0; r < Math.min(rows, old.commandAlpha.length); r += 1) {
      for (let c = 0; c < Math.min(cols, old.commandAlpha[r].length); c += 1) {
        next.cells.alpha[r][c] = old.alpha?.[r]?.[c] ?? next.cells.alpha[r][c];
        next.cells.theta[r][c] = old.theta?.[r]?.[c] ?? next.cells.theta[r][c];
        next.cells.z[r][c] = old.z?.[r]?.[c] ?? next.cells.z[r][c];
        next.cells.commandAlpha[r][c] = old.commandAlpha[r][c];
        next.cells.commandZ[r][c] = old.commandZ[r][c];
        next.cells.locked[r][c] = old.locked[r][c];
        next.cells.lockAlpha[r][c] = old.lockAlpha?.[r]?.[c] ?? next.grid.initialAlpha;
        next.cells.actuatorAllowed[r][c] = old.actuatorAllowed?.[r]?.[c] ?? true;
      }
    }
    return next;
  }

  function clearCommands(state) {
    const { rows, cols } = state.grid;
    const derived = derivedCellFields(rows, cols, state.grid.initialAlpha);
    state.cells.alpha = derived.alpha;
    state.cells.theta = derived.theta;
    state.cells.z = derived.z;
    state.cells.commandAlpha = matrix(rows, cols, 0);
    state.cells.commandZ = matrix(rows, cols, 0);
    state.cells.locked = matrix(rows, cols, false);
    state.cells.lockAlpha = matrix(rows, cols, state.grid.initialAlpha);
    if (!state.cells.actuatorAllowed) state.cells.actuatorAllowed = matrix(rows, cols, true);
  }

  function cloneData(value) {
    if (value === undefined) return undefined;
    return JSON.parse(JSON.stringify(value));
  }

  function snapshotState(state) {
    const snapshot = stateWithFreshDerivedCells(state);
    return {
      grid: cloneData(snapshot.grid),
      cells: cloneData(snapshot.cells),
      view: cloneData(snapshot.view),
      target: cloneData(snapshot.target),
      inverse: cloneData(snapshot.inverse),
      selection: cloneData(snapshot.selection),
      experiment: {
        presetName: snapshot.experiment.presetName,
        notes: snapshot.experiment.notes,
      },
    };
  }

  function restoreSnapshot(state, snapshot) {
    if (!snapshot) return state;
    const snapshotGrid = cloneData(snapshot.grid) || {};
    state.grid = { ...state.grid, ...snapshotGrid };
    if (!Object.prototype.hasOwnProperty.call(snapshotGrid, "hardwareProfile")) {
      delete state.grid.hardwareProfile;
    }
    ensureGridSchema(state);
    state.cells = cloneData(snapshot.cells);
    ensureCellSchema(state);
    state.view = { ...state.view, ...cloneData(snapshot.view) };
    state.target = { ...state.target, ...cloneData(snapshot.target) };
    state.inverse = { ...state.inverse, ...(cloneData(snapshot.inverse) || {}) };
    state.selection = cloneData(snapshot.selection || state.selection);
    state.experiment.presetName = snapshot.experiment?.presetName || state.experiment.presetName;
    state.experiment.notes = snapshot.experiment?.notes || state.experiment.notes;
    return state;
  }

  function lerp(a, b, t) {
    return a + (b - a) * t;
  }

  function interpolateMatrix(a, b, t, fallback) {
    if (!Array.isArray(a) || !Array.isArray(b) || a.length !== b.length || a[0]?.length !== b[0]?.length) return cloneData(t < 1 ? a || fallback : b || fallback);
    return a.map((row, r) =>
      row.map((value, c) => {
        if (typeof value === "number" && typeof b[r][c] === "number") return lerp(value, b[r][c], t);
        return t < 0.5 ? value : b[r][c];
      })
    );
  }

  function interpolateSnapshots(from, to, t) {
    if (!from) return cloneData(to);
    if (!to) return cloneData(from);
    const amount = Math.max(0, Math.min(1, t));
    const rowsMatch = from.grid?.rows === to.grid?.rows;
    const colsMatch = from.grid?.cols === to.grid?.cols;
    if (!rowsMatch || !colsMatch) return amount < 1 ? cloneData(from) : cloneData(to);
    const out = cloneData(amount < 0.5 ? from : to);
    out.grid = { ...to.grid };
    for (const key of ["cellSize", "backlash", "couplingGain", "zCouplingGain", "pinRadius", "holeRadius", "paperSideLengthMm", "paperHoleToleranceMm", "initialAlpha", "alphaMin", "alphaMax", "zTravelLimit", "alphaContractLimit", "alphaExpandLimit"]) {
      if (typeof from.grid?.[key] === "number" && typeof to.grid?.[key] === "number") out.grid[key] = lerp(from.grid[key], to.grid[key], amount);
    }
    out.cells = {
      alpha: interpolateMatrix(from.cells?.alpha, to.cells?.alpha, amount, []),
      theta: interpolateMatrix(from.cells?.theta, to.cells?.theta, amount, []),
      z: interpolateMatrix(from.cells?.z, to.cells?.z, amount, []),
      commandAlpha: interpolateMatrix(from.cells?.commandAlpha, to.cells?.commandAlpha, amount, []),
      commandZ: interpolateMatrix(from.cells?.commandZ, to.cells?.commandZ, amount, []),
      locked: interpolateMatrix(from.cells?.locked, to.cells?.locked, amount, []),
      lockAlpha: interpolateMatrix(from.cells?.lockAlpha, to.cells?.lockAlpha, amount, []),
      actuatorAllowed: interpolateMatrix(from.cells?.actuatorAllowed, to.cells?.actuatorAllowed, amount, []),
    };
    out.target = { ...to.target };
    for (const key of ["amplitude", "frequency", "optimizerIterations", "actuatorPenalty", "travelPenalty"]) {
      if (typeof from.target?.[key] === "number" && typeof to.target?.[key] === "number") out.target[key] = lerp(from.target[key], to.target[key], amount);
    }
    out.selection = amount < 0.5 ? cloneData(from.selection || to.selection) : cloneData(to.selection || from.selection);
    out.view = { ...to.view };
    out.inverse = amount < 1 ? cloneData(from.inverse || to.inverse) : cloneData(to.inverse || from.inverse);
    return out;
  }

  function recordEvent(state, event) {
    if (!state.experiment.initialSnapshot) state.experiment.initialSnapshot = snapshotState(state);
    state.experiment.eventList.push({
      ...event,
      at: event.at || new Date().toISOString(),
      snapshot: snapshotState(state),
    });
    state.timeline.index = state.experiment.eventList.length;
  }

  function serialize(state) {
    const snapshot = stateWithFreshDerivedCells(state);
    return JSON.stringify(
      {
        schema: "rad-sim.browser.v1",
        savedAt: new Date().toISOString(),
        grid: snapshot.grid,
        cells: snapshot.cells,
        view: snapshot.view,
        target: snapshot.target,
        inverse: snapshot.inverse,
        timeline: snapshot.timeline,
        selection: snapshot.selection,
        experiment: snapshot.experiment,
      },
      null,
      2
    );
  }

  function commandSummary(snapshot) {
    const commandAlpha = snapshot?.cells?.commandAlpha || [];
    const commandZ = snapshot?.cells?.commandZ || [];
    let active = 0;
    let maxAbsAlpha = 0;
    let maxAbsZ = 0;
    for (let r = 0; r < commandAlpha.length; r += 1) {
      for (let c = 0; c < commandAlpha[r].length; c += 1) {
        const alpha = Number(commandAlpha[r][c]) || 0;
        const z = Number(commandZ?.[r]?.[c]) || 0;
        if (Math.abs(alpha) > 1e-9 || Math.abs(z) > 1e-9) active += 1;
        maxAbsAlpha = Math.max(maxAbsAlpha, Math.abs(alpha));
        maxAbsZ = Math.max(maxAbsZ, Math.abs(z));
      }
    }
    return { active, maxAbsAlpha, maxAbsZ };
  }

  function exportExperimentSequence(state) {
    const events = state.experiment.eventList || [];
    const frames = events.map((event, index) => ({
      index: index + 1,
      type: event.type,
      name: event.name || event.target || event.previewType || "",
      at: event.at,
      commandSummary: commandSummary(event.snapshot),
      snapshot: cloneData(event.snapshot),
    }));
    return JSON.stringify(
      {
        schema: "rad-sim.sequence.v1",
        savedAt: new Date().toISOString(),
        summary: {
          frameCount: frames.length,
          keyframes: events.filter((event) => event.type === "keyframe" || event.type === "inverse-preview-keyframe").length,
          eventTypes: Array.from(new Set(events.map((event) => event.type))),
          finalCommandSummary: commandSummary(frames.at(-1)?.snapshot || state.experiment.initialSnapshot || snapshotState(state)),
        },
        grid: cloneData(state.grid),
        target: cloneData(state.target),
        timeline: cloneData(state.timeline),
        initialSnapshot: cloneData(state.experiment.initialSnapshot || snapshotState(state)),
        frames,
      },
      null,
      2
    );
  }

  function importExperimentSequence(state, text) {
    const parsed = JSON.parse(text);
    if (parsed.schema !== "rad-sim.sequence.v1") {
      throw new Error("Unsupported RAD sequence JSON schema.");
    }
    const frames = parsed.frames || [];
    const initialSnapshot = cloneData(parsed.initialSnapshot);
    if (!initialSnapshot?.grid?.rows || !initialSnapshot?.grid?.cols) {
      throw new Error("RAD sequence is missing its initial snapshot.");
    }
    if (!frames.length) {
      throw new Error("RAD sequence has no frames to import.");
    }
    const importedEvents = frames
      .filter((frame) => frame.snapshot)
      .map((frame, index) => ({
        type: frame.type || "sequence-frame",
        name: frame.name || "",
        at: frame.at || parsed.savedAt || new Date().toISOString(),
        sequenceIndex: frame.index || index + 1,
        commandSummary: cloneData(frame.commandSummary || commandSummary(frame.snapshot)),
        snapshot: cloneData(frame.snapshot),
      }));
    if (!importedEvents.length) {
      throw new Error("RAD sequence frames are missing snapshots.");
    }
    const lastSnapshot = importedEvents.at(-1).snapshot;
    restoreSnapshot(state, lastSnapshot);
    const importedGrid = cloneData(parsed.grid || lastSnapshot.grid) || {};
    state.grid = { ...state.grid, ...importedGrid };
    if (!Object.prototype.hasOwnProperty.call(importedGrid, "hardwareProfile")) {
      delete state.grid.hardwareProfile;
    }
    ensureGridSchema(state);
    state.target = { ...state.target, ...cloneData(parsed.target || lastSnapshot.target) };
    state.timeline = {
      ...state.timeline,
      ...(cloneData(parsed.timeline) || {}),
      index: importedEvents.length,
    };
    state.experiment.initialSnapshot = initialSnapshot;
    state.experiment.eventList = importedEvents;
    return state;
  }

  function deserialize(text) {
    const parsed = JSON.parse(text);
    if (parsed.schema !== "rad-sim.browser.v1") {
      throw new Error("Unsupported RAD simulator JSON schema.");
    }
    const rows = parsed.grid.rows;
    const cols = parsed.grid.cols;
    const state = createState(rows, cols);
    state.grid = { ...state.grid, ...parsed.grid };
    if (!Object.prototype.hasOwnProperty.call(parsed.grid, "hardwareProfile")) {
      delete state.grid.hardwareProfile;
    }
    ensureGridSchema(state);
    if (!state.grid.zTravelLimit) state.grid.zTravelLimit = 0.8;
    if (!state.grid.alphaContractLimit) state.grid.alphaContractLimit = 0.55;
    if (!state.grid.alphaExpandLimit) state.grid.alphaExpandLimit = 0.35;
    state.cells = parsed.cells;
    ensureCellSchema(state);
    state.view = { ...state.view, ...parsed.view };
    if (!state.view.cellVisualMode) state.view.cellVisualMode = "abstract";
    if (!state.view.simulationMode) state.view.simulationMode = "kinematic";
    if (state.view.isolateSelected === undefined) state.view.isolateSelected = false;
    if (state.view.explodedSelected === undefined) state.view.explodedSelected = false;
    if (state.view.quickDockCollapsed === undefined) state.view.quickDockCollapsed = false;
    if (state.view.paintMode === undefined) state.view.paintMode = false;
    if (state.view.paintRadius === undefined) state.view.paintRadius = 0;
    if (!state.view.surfaceInterpolation) state.view.surfaceInterpolation = "smooth";
    if (!state.view.surfaceSubdivisions) state.view.surfaceSubdivisions = 5;
    if (state.view.surfaceContoursVisible === undefined) state.view.surfaceContoursVisible = false;
    if (!state.view.contourMode) state.view.contourMode = "membrane";
    if (!state.view.animationResponse) state.view.animationResponse = 0.18;
    if (state.view.referenceVisible === undefined) state.view.referenceVisible = true;
    if (state.view.displacementVectorsVisible === undefined) state.view.displacementVectorsVisible = false;
    if (state.view.targetErrorVectorsVisible === undefined) state.view.targetErrorVectorsVisible = true;
    if (state.view.surfaceNormalsVisible === undefined) state.view.surfaceNormalsVisible = false;
    if (state.view.stopsVisible === undefined) state.view.stopsVisible = false;
    if (state.view.pivotsVisible === undefined) state.view.pivotsVisible = false;
    if (state.view.linkagesVisible === undefined) state.view.linkagesVisible = false;
    if (state.view.fastenersVisible === undefined) state.view.fastenersVisible = false;
    if (state.view.actuatorsVisible === undefined) state.view.actuatorsVisible = false;
    if (!state.view.actuatorDisplayMode) state.view.actuatorDisplayMode = "selected";
    if (state.view.influenceFootprintVisible === undefined) state.view.influenceFootprintVisible = true;
    if (state.view.measurementLabelsVisible === undefined) state.view.measurementLabelsVisible = false;
    if (!state.view.measurementMode) state.view.measurementMode = "alpha";
    if (!state.view.camera) state.view.camera = createState(rows, cols).view.camera;
    if (!state.view.camera.projection) state.view.camera.projection = "perspective";
    if (!state.view.camera.target) state.view.camera.target = { x: 0, y: 0, z: 0 };
    state.target = { ...state.target, ...parsed.target };
    if (!state.target.customExpression) state.target.customExpression = "0.35 * exp(-3 * (nr*nr + nc*nc))";
    if (state.target.expressionError === undefined) state.target.expressionError = "";
    state.inverse = { ...state.inverse, ...parsed.inverse };
    if (state.inverse.preview === undefined) state.inverse.preview = null;
    if (!state.inverse.sensitivity) {
      state.inverse.sensitivity = { candidates: [], map: matrix(rows, cols, 0), stepZ: 0.12, stepAlpha: 0.12, controllableCells: 0, meanGain: 0, maxGain: 0 };
    }
    if (!state.inverse.jacobian) {
      state.inverse.jacobian = { columns: [], coverageMap: matrix(rows, cols, 0), actuatorCount: 0, columnCount: 0, meanCoverage: 0, maxCoverage: 0, stepZ: 0.12, stepAlpha: 0.12, conditionEstimate: 0 };
    }
    if (!state.inverse.linearSolution) {
      state.inverse.linearSolution = { commands: [], history: [], steps: 0, baseError: 0, predictedError: 0, projectedError: 0, projectedActuators: 0 };
    }
    if (state.inverse.physicalValidation === undefined) state.inverse.physicalValidation = null;
    state.timeline = { ...state.timeline, ...parsed.timeline };
    if (state.timeline.smooth === undefined) state.timeline.smooth = true;
    if (!state.timeline.transitionMs) state.timeline.transitionMs = 900;
    if (!state.timeline.keyframeName) state.timeline.keyframeName = "pose";
    state.selection = parsed.selection || state.selection;
    state.experiment = { ...state.experiment, ...parsed.experiment };
    return state;
  }

  RAD.matrix = matrix;
  RAD.defaultHardwareProfile = defaultHardwareProfile;
  RAD.createState = createState;
  RAD.resizeState = resizeState;
  RAD.clearCommands = clearCommands;
  RAD.updateDerivedCells = updateDerivedCells;
  RAD.snapshotState = snapshotState;
  RAD.restoreSnapshot = restoreSnapshot;
  RAD.interpolateSnapshots = interpolateSnapshots;
  RAD.recordEvent = recordEvent;
  RAD.serialize = serialize;
  RAD.exportExperimentSequence = exportExperimentSequence;
  RAD.importExperimentSequence = importExperimentSequence;
  RAD.deserialize = deserialize;
})();
