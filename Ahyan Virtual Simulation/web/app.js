(function () {
  "use strict";

  const RAD = window.RAD;
  const mount = document.getElementById("threeMount");
  const status = document.getElementById("statusLine");
  const overlayLegend = {
    title: document.getElementById("overlayLegendTitle"),
    min: document.getElementById("overlayLegendMin"),
    max: document.getElementById("overlayLegendMax"),
  };
  const selectionHud = {
    cell: document.getElementById("hudCell"),
    alpha: document.getElementById("hudAlpha"),
    theta: document.getElementById("hudTheta"),
    height: document.getElementById("hudHeight"),
    command: document.getElementById("hudCommand"),
    state: document.getElementById("hudState"),
  };
  const orientationHud = {
    viewAxis: document.getElementById("cameraViewAxis"),
    upAxis: document.getElementById("cameraUpAxis"),
    projection: document.getElementById("cameraProjection"),
  };
  const overlayLegendText = {
    alpha: ["Alpha / Dilation", "contracted", "expanded"],
    height: ["Height", "low z", "high z"],
    zresidual: ["Z Residual", "negative spillover", "positive spillover"],
    influence: ["Backlash Influence", "negative/free", "positive/coupled"],
    theta: ["Rotation Theta", "closed angle", "open angle"],
    error: ["Target Error", "below target", "above target"],
    travel: ["Command Travel", "idle", "high travel"],
    saturation: ["Actuator Saturation", "inside limits", "near limit"],
    strain: ["Linkage Strain", "compressed", "stretched"],
    modelError: ["Model Disagreement", "matches kinematic", "relaxed shift"],
    calibrationError: ["Calibration Error", "matched cell", "large measured error"],
    calibrationResidual: ["Calibration Residual", "fit absorbed", "local residual"],
    displacement: ["Reference Displacement", "small shift", "large shift"],
    slope: ["Surface Slope", "flat", "steep"],
    inverse: ["Inverse Plan", "low contribution", "high contribution"],
    sensitivity: ["Sensitivity", "low response", "high response"],
    reachability: ["Reachability", "weak coverage", "strong coverage"],
    state: ["Cell State", "free/idle", "locked/active"],
  };
  let state = RAD.createState(7, 7);
  let sim = RAD.simulateActive(state);
  let renderer = null;
  let ui = null;
  let focusMode = false;
  let allMetricsVisible = false;
  let workspaceMode = "inspect";
  let hoveredCell = null;

  function setStatus(message) {
    status.textContent = message;
  }

  function setActiveView(mode) {
    document.querySelectorAll("[data-view-mode]").forEach((button) => {
      const active = button.dataset.viewMode === mode;
      button.classList.toggle("is-active", active);
      button.setAttribute("aria-pressed", String(active));
    });
  }

  function setActiveWorkspaceMode(mode) {
    document.querySelectorAll("[data-workspace-mode]").forEach((button) => {
      const active = button.dataset.workspaceMode === mode;
      button.classList.toggle("is-active", active);
      button.setAttribute("aria-pressed", String(active));
    });
  }

  function setProjectionMode(mode) {
    const projection = renderer.setProjectionMode(mode);
    syncCameraState();
    updateProjectionButton();
    updateOrientationHud();
    updateStatusLine();
    return projection;
  }

  function updateProjectionButton() {
    const button = document.getElementById("projectionMode");
    const projection = renderer?.projectionMode || state.view.camera?.projection || "perspective";
    const orthographic = projection === "orthographic";
    button.classList.toggle("is-active", orthographic);
    button.setAttribute("aria-pressed", String(orthographic));
    button.textContent = orthographic ? "Perspective" : "Ortho";
  }

  function setView(mode) {
    renderer.setView(mode);
    syncCameraState();
    setActiveView(renderer.viewMode || mode);
    updateOrientationHud();
    updateStatusLine();
  }

  function resetCameraView() {
    renderer.resetView();
    syncCameraState();
    setActiveView(renderer.viewMode || "iso");
    updateOrientationHud();
    updateStatusLine();
  }

  function frameSelectedCell() {
    renderer.frameCell(state, sim);
    syncCameraState();
    setActiveView(renderer.viewMode || "custom");
    updateOrientationHud();
    updateStatusLine();
  }

  function setIsolateCell(enabled) {
    state.view.isolateSelected = Boolean(enabled);
    if (!state.view.isolateSelected) state.view.explodedSelected = false;
    const button = document.getElementById("isolateCell");
    button.classList.toggle("is-active", state.view.isolateSelected);
    button.setAttribute("aria-pressed", String(state.view.isolateSelected));
    button.textContent = state.view.isolateSelected ? "Show Lattice" : "Isolate Cell";
    renderAll(state);
    if (state.view.isolateSelected) frameSelectedCell();
  }

  function setExplodeCell(enabled) {
    state.view.explodedSelected = Boolean(enabled);
    if (state.view.explodedSelected) {
      ui.applyCellVisualMode("mechanism");
      state.view.isolateSelected = true;
    }
    renderAll(state);
    if (state.view.explodedSelected) frameSelectedCell();
  }

  function syncCameraState() {
    if (!renderer || !state?.view) return;
    state.view.camera = renderer.getCameraState();
  }

  function setFocusMode(enabled) {
    focusMode = enabled;
    document.body.classList.toggle("is-focus-mode", focusMode);
    const button = document.getElementById("focusMode");
    button.classList.toggle("is-active", focusMode);
    button.setAttribute("aria-pressed", String(focusMode));
    button.textContent = focusMode ? "Controls" : "Focus";
    renderer.resize();
  }

  function setMetricDensity(showAll) {
    allMetricsVisible = showAll;
    document.body.classList.toggle("show-all-metrics", allMetricsVisible);
    const button = document.getElementById("metricDensity");
    button.classList.toggle("is-active", allMetricsVisible);
    button.setAttribute("aria-pressed", String(allMetricsVisible));
    button.textContent = allMetricsVisible ? "Core Metrics" : "All Metrics";
  }

  function applyWorkspaceMode(mode) {
    workspaceMode = ["inspect", "edit", "analyze"].includes(mode) ? mode : "inspect";
    const targetAvailable = state.target.type !== "none";
    const view = state.view;
    view.explodedSelected = false;
    if (workspaceMode === "inspect") {
      Object.assign(view, {
        cellVisualMode: "abstract",
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
        measurementsVisible: true,
        measurementMode: "alpha",
        targetVisible: false,
        targetErrorVectorsVisible: false,
        surfaceNormalsVisible: false,
        surfaceContoursVisible: false,
        contourMode: "membrane",
      });
      setMetricDensity(false);
    } else if (workspaceMode === "edit") {
      Object.assign(view, {
        cellVisualMode: "abstract",
        overlayMode: "state",
        membraneVisible: true,
        referenceVisible: false,
        displacementVectorsVisible: false,
        gapsVisible: true,
        stopsVisible: true,
        pivotsVisible: true,
        linkagesVisible: true,
        fastenersVisible: false,
        actuatorsVisible: true,
        actuatorDisplayMode: "selected",
        measurementsVisible: true,
        measurementMode: "all",
        targetVisible: targetAvailable,
        targetErrorVectorsVisible: false,
        surfaceNormalsVisible: false,
        surfaceContoursVisible: false,
        contourMode: "membrane",
      });
      setMetricDensity(false);
    } else {
      Object.assign(view, {
        cellVisualMode: "abstract",
        overlayMode: targetAvailable ? "error" : "height",
        membraneVisible: true,
        referenceVisible: true,
        displacementVectorsVisible: true,
        gapsVisible: true,
        stopsVisible: true,
        pivotsVisible: false,
        linkagesVisible: true,
        fastenersVisible: false,
        actuatorsVisible: true,
        actuatorDisplayMode: "planned",
        measurementsVisible: true,
        measurementMode: targetAvailable ? "height" : "all",
        targetVisible: targetAvailable,
        targetErrorVectorsVisible: targetAvailable,
        surfaceNormalsVisible: true,
        surfaceContoursVisible: true,
        contourMode: targetAvailable ? "error" : "membrane",
      });
      setMetricDensity(true);
    }
    setActiveWorkspaceMode(workspaceMode);
    renderAll(state);
  }

  function setupControlPanelSections() {
    const initiallyCollapsed = new Set(["Target Surface", "Inverse Plan", "Timeline", "Cell Grid", "Persistence"]);
    document.querySelectorAll(".control-panel section").forEach((section) => {
      const heading = section.querySelector("h2");
      if (!heading || heading.querySelector("button")) return;
      const title = heading.textContent.trim();
      section.classList.add("control-section");
      heading.classList.add("section-heading");
      const button = document.createElement("button");
      button.type = "button";
      button.className = "section-toggle";
      button.setAttribute("aria-label", `Toggle ${title} controls`);
      heading.appendChild(button);
      const setCollapsed = (collapsed) => {
        section.classList.toggle("is-collapsed", collapsed);
        button.textContent = collapsed ? "+" : "-";
        button.setAttribute("aria-expanded", String(!collapsed));
      };
      setCollapsed(initiallyCollapsed.has(title));
      heading.addEventListener("click", (event) => {
        if (event.target.closest("button") || event.target === heading) setCollapsed(!section.classList.contains("is-collapsed"));
      });
      button.addEventListener("click", (event) => {
        event.stopPropagation();
        setCollapsed(!section.classList.contains("is-collapsed"));
      });
    });
  }

  function updateStatusLine() {
    if (!renderer || !sim) return;
    const hoverText = hoveredCell ? ` Hover: r${hoveredCell.r}, c${hoveredCell.c}.` : "";
    const isolateText = state.view.isolateSelected ? " Selected cell isolated." : "";
    const explodedText = state.view.explodedSelected ? " Exploded cell detail active." : "";
    const paintText = state.view.paintMode ? " Paint clicks apply current command." : "";
    const projectionText = renderer.projectionMode === "orthographic" ? "orthographic" : "perspective";
    const modelText = sim.metrics.physicalPreview ? " Model: spring preview." : " Model: kinematic.";
    setStatus(
      `Workspace: ${workspaceMode}. View: ${renderer.viewMode || "iso"} ${projectionText}.${hoverText}${isolateText}${explodedText}${paintText}${modelText} Cell visual: ${state.view.cellVisualMode || "abstract"}. Overlay: ${state.view.overlayMode}. Theta ${sim.metrics.minTheta.toFixed(1)} to ${sim.metrics.maxTheta.toFixed(1)} deg. Target error ${sim.metrics.rmsTargetError.toFixed(3)}. Actuators ${sim.metrics.recommendedActuators}.`
    );
  }

  function updateOverlayLegend() {
    const [title, min, max] = overlayLegendText[state.view.overlayMode] || overlayLegendText.alpha;
    overlayLegend.title.textContent = title;
    overlayLegend.min.textContent = min;
    overlayLegend.max.textContent = max;
  }

  function updateSelectionHud() {
    const { r, c } = state.selection;
    const locked = state.cells.locked[r][c];
    const allowed = state.cells.actuatorAllowed?.[r]?.[c] !== false;
    const commandAlpha = state.cells.commandAlpha[r][c];
    const commandZ = state.cells.commandZ[r][c];
    selectionHud.cell.textContent = `Cell r${r} c${c}`;
    selectionHud.alpha.textContent = sim.alpha[r][c].toFixed(3);
    selectionHud.theta.textContent = `${sim.theta[r][c].toFixed(1)} deg`;
    selectionHud.height.textContent = sim.height[r][c].toFixed(3);
    selectionHud.command.textContent = `a ${commandAlpha.toFixed(2)} / z ${commandZ.toFixed(2)}`;
    selectionHud.state.textContent = locked ? "locked" : allowed ? "free" : "masked";
  }

  function updateOrientationHud() {
    if (!renderer) return;
    const orientation = renderer.getOrientationState();
    orientationHud.viewAxis.textContent = orientation.viewAxis;
    orientationHud.upAxis.textContent = orientation.upAxis;
    orientationHud.projection.textContent = orientation.projection;
  }

  function renderAll(nextState) {
    state = nextState;
    const savedCamera = state.view.camera;
    sim = RAD.simulateActive(state);
    RAD.updateDerivedCells(state, sim);
    ui.setState(state);
    ui.updateLabels(sim);
    ui.renderGrid(sim);
    renderer.renderState(state, sim);
    if (savedCamera) renderer.applyCameraState(savedCamera);
    syncCameraState();
    updateProjectionButton();
    ui.updateRendererDiagnostics(renderer.getDiagnostics());
    updateOverlayLegend();
    updateSelectionHud();
    updateOrientationHud();
    const isolateButton = document.getElementById("isolateCell");
    isolateButton.classList.toggle("is-active", state.view.isolateSelected === true);
    isolateButton.setAttribute("aria-pressed", String(state.view.isolateSelected === true));
    isolateButton.textContent = state.view.isolateSelected ? "Show Lattice" : "Isolate Cell";
    const explodeButton = document.getElementById("explodeCell");
    explodeButton.classList.toggle("is-active", state.view.explodedSelected === true);
    explodeButton.setAttribute("aria-pressed", String(state.view.explodedSelected === true));
    explodeButton.textContent = state.view.explodedSelected ? "Assemble Cell" : "Explode Cell";
    updateStatusLine();
  }

  try {
    renderer = new RAD.RadRenderer(mount, (r, c) => {
      if (state.view.paintMode) ui.paintCell(r, c);
      else ui.select(r, c);
    }, (mode) => {
      syncCameraState();
      setActiveView(mode);
      updateOrientationHud();
      updateStatusLine();
    }, (cell) => {
      hoveredCell = cell;
      updateStatusLine();
    });
    ui = new RAD.RadUI(state, renderAll);
    setupControlPanelSections();
    document.getElementById("modeInspect").addEventListener("click", () => applyWorkspaceMode("inspect"));
    document.getElementById("modeEdit").addEventListener("click", () => applyWorkspaceMode("edit"));
    document.getElementById("modeAnalyze").addEventListener("click", () => applyWorkspaceMode("analyze"));
    document.getElementById("viewIso").addEventListener("click", () => setView("iso"));
    document.getElementById("viewTop").addEventListener("click", () => setView("top"));
    document.getElementById("viewFront").addEventListener("click", () => setView("front"));
    document.getElementById("viewSide").addEventListener("click", () => setView("side"));
    document.getElementById("projectionMode").addEventListener("click", () => setProjectionMode(renderer.projectionMode === "orthographic" ? "perspective" : "orthographic"));
    document.getElementById("frameCell").addEventListener("click", () => frameSelectedCell());
    document.getElementById("isolateCell").addEventListener("click", () => setIsolateCell(!state.view.isolateSelected));
    document.getElementById("explodeCell").addEventListener("click", () => setExplodeCell(!state.view.explodedSelected));
    document.getElementById("resetView").addEventListener("click", () => resetCameraView());
    document.getElementById("focusMode").addEventListener("click", () => setFocusMode(!focusMode));
    document.getElementById("metricDensity").addEventListener("click", () => setMetricDensity(!allMetricsVisible));
    applyWorkspaceMode(workspaceMode);
    setActiveView(renderer.viewMode || "iso");
    window.RAD_APP = {
      get state() {
        return state;
      },
      get simulation() {
        return sim;
      },
      get renderer() {
        return renderer;
      },
      update: (nextState = state) => renderAll(nextState),
      getState: () => state,
      getSimulation: () => sim,
      getRendererDiagnostics: () => renderer.getDiagnostics(),
      applyWorkspaceMode,
      setProjectionMode,
      setIsolateCell,
      setExplodeCell,
      syncCameraState,
      frameSelectedCell,
      resetCameraView,
      serialize: () => RAD.serialize(state),
    };
    window.setInterval(() => ui.updateRendererDiagnostics(renderer.getDiagnostics()), 1000);
  } catch (error) {
    console.error(error);
    setStatus(`3D renderer failed to initialize: ${error.message}`);
    mount.innerHTML = '<div class="status-line">Three.js could not load. Check your internet connection or use a local server with the CDN available.</div>';
  }
})();
