(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  class RadUI {
    constructor(state, onChange) {
      this.state = state;
      this.onChange = onChange;
      this.els = {
        rows: document.getElementById("rows"),
        cols: document.getElementById("cols"),
        backlash: document.getElementById("backlash"),
        coupling: document.getElementById("coupling"),
        zCoupling: document.getElementById("zCoupling"),
        pinRadius: document.getElementById("pinRadius"),
        holeRadius: document.getElementById("holeRadius"),
        quickDock: document.querySelector(".quick-actuation-dock"),
        toggleQuickDock: document.getElementById("toggleQuickDock"),
        quickDockMini: document.getElementById("quickDockMini"),
        paintMode: document.getElementById("paintMode"),
        paintRadius: document.getElementById("paintRadius"),
        alphaCommand: document.getElementById("alphaCommand"),
        zCommand: document.getElementById("zCommand"),
        locked: document.getElementById("locked"),
        actuatorAllowed: document.getElementById("actuatorAllowed"),
        cellVisualMode: document.getElementById("cellVisualMode"),
        simulationMode: document.getElementById("simulationMode"),
        overlayMode: document.getElementById("overlayMode"),
        showMembrane: document.getElementById("showMembrane"),
        showReference: document.getElementById("showReference"),
        showDisplacementVectors: document.getElementById("showDisplacementVectors"),
        showGaps: document.getElementById("showGaps"),
        showStops: document.getElementById("showStops"),
        showPivots: document.getElementById("showPivots"),
        showLinkages: document.getElementById("showLinkages"),
        showFasteners: document.getElementById("showFasteners"),
        showActuators: document.getElementById("showActuators"),
        actuatorDisplayMode: document.getElementById("actuatorDisplayMode"),
        showInfluenceFootprint: document.getElementById("showInfluenceFootprint"),
        showMeasurements: document.getElementById("showMeasurements"),
        showMeasurementLabels: document.getElementById("showMeasurementLabels"),
        showTargetErrorVectors: document.getElementById("showTargetErrorVectors"),
        showSurfaceNormals: document.getElementById("showSurfaceNormals"),
        measurementMode: document.getElementById("measurementMode"),
        showTarget: document.getElementById("showTarget"),
        surfaceInterpolation: document.getElementById("surfaceInterpolation"),
        surfaceSubdivisions: document.getElementById("surfaceSubdivisions"),
        showSurfaceContours: document.getElementById("showSurfaceContours"),
        contourMode: document.getElementById("contourMode"),
        animationResponse: document.getElementById("animationResponse"),
        targetType: document.getElementById("targetType"),
        targetAmp: document.getElementById("targetAmp"),
        targetFreq: document.getElementById("targetFreq"),
        customTargetExpression: document.getElementById("customTargetExpression"),
        targetExpressionStatus: document.getElementById("targetExpressionStatus"),
        optimizerIterations: document.getElementById("optimizerIterations"),
        actuatorPenalty: document.getElementById("actuatorPenalty"),
        travelPenalty: document.getElementById("travelPenalty"),
        maxActuators: document.getElementById("maxActuators"),
        zTravelLimit: document.getElementById("zTravelLimit"),
        alphaContractLimit: document.getElementById("alphaContractLimit"),
        alphaExpandLimit: document.getElementById("alphaExpandLimit"),
        timelineIndex: document.getElementById("timelineIndex"),
        timelineSpeed: document.getElementById("timelineSpeed"),
        transitionMs: document.getElementById("transitionMs"),
        smoothTimeline: document.getElementById("smoothTimeline"),
        keyframeName: document.getElementById("keyframeName"),
        eventList: document.getElementById("eventList"),
        inversePlanHistory: document.getElementById("inversePlanHistory"),
        inverseCandidateList: document.getElementById("inverseCandidateList"),
        inverseStepPreview: document.getElementById("inverseStepPreview"),
        cellGrid: document.getElementById("cellGrid"),
        sequenceChart: document.getElementById("sequenceChart"),
        fileInput: document.getElementById("fileInput"),
        sequenceFileInput: document.getElementById("sequenceFileInput"),
        saveObj: document.getElementById("saveObj"),
        operatorCommitLock: document.getElementById("operatorCommitLock"),
        operatorReleaseLock: document.getElementById("operatorReleaseLock"),
        operatorCheckOrder: document.getElementById("operatorCheckOrder"),
        operatorOrderState: document.getElementById("operatorOrderState"),
        operatorAlphaError: document.getElementById("operatorAlphaError"),
        operatorHeightError: document.getElementById("operatorHeightError"),
      };
      this.playTimer = null;
      this.transitionFrame = null;
      this.timelinePlaying = false;
      this.lastSequenceAnalysis = null;
      this.lastSequenceAnalysisSignature = "";
      this.sequenceChartLayout = null;
      this.hoveredSequenceFrame = null;
      this.lastOperatorDiagnostic = null;
      this.bind();
      this.syncControls();
    }

    bind() {
      this.els.rows.addEventListener("input", () => {
        this.state = RAD.resizeState(this.state, Number(this.els.rows.value), this.state.grid.cols);
        this.onChange(this.state);
      });
      this.els.cols.addEventListener("input", () => {
        this.state = RAD.resizeState(this.state, this.state.grid.rows, Number(this.els.cols.value));
        this.onChange(this.state);
      });
      this.els.backlash.addEventListener("input", () => {
        this.state.grid.backlash = Number(this.els.backlash.value);
        this.onChange(this.state);
      });
      this.els.coupling.addEventListener("input", () => {
        this.state.grid.couplingGain = Number(this.els.coupling.value);
        this.onChange(this.state);
      });
      this.els.zCoupling.addEventListener("input", () => {
        this.state.grid.zCouplingGain = Number(this.els.zCoupling.value);
        this.onChange(this.state);
      });
      this.els.pinRadius.addEventListener("input", () => {
        this.state.grid.pinRadius = Number(this.els.pinRadius.value);
        if (this.state.grid.holeRadius < this.state.grid.pinRadius) this.state.grid.holeRadius = this.state.grid.pinRadius;
        this.onChange(this.state);
      });
      this.els.holeRadius.addEventListener("input", () => {
        this.state.grid.holeRadius = Math.max(Number(this.els.holeRadius.value), Number(this.state.grid.pinRadius || 0));
        this.onChange(this.state);
      });
      this.els.toggleQuickDock.addEventListener("click", () => {
        this.state.view.quickDockCollapsed = !this.state.view.quickDockCollapsed;
        this.syncQuickDock();
        this.onChange(this.state);
      });
      this.els.paintMode.addEventListener("click", () => {
        this.state.view.paintMode = !this.state.view.paintMode;
        this.syncPaintMode();
        this.onChange(this.state);
      });
      this.els.paintRadius.addEventListener("input", () => {
        this.state.view.paintRadius = Number(this.els.paintRadius.value);
        this.updateLabels(null);
        this.onChange(this.state);
      });
      this.els.alphaCommand.addEventListener("input", () => this.previewSelected());
      this.els.zCommand.addEventListener("input", () => this.previewSelected());
      this.els.locked.addEventListener("change", () => this.previewSelected());
      this.els.actuatorAllowed.addEventListener("change", () => this.updateSelectedActuatorMask());
      this.els.cellVisualMode.addEventListener("change", () => {
        this.applyCellVisualMode(this.els.cellVisualMode.value);
        this.onChange(this.state);
      });
      this.els.simulationMode.addEventListener("change", () => {
        this.state.view.simulationMode = this.els.simulationMode.value === "springPreview" ? "springPreview" : "kinematic";
        this.onChange(this.state);
      });
      this.els.overlayMode.addEventListener("change", () => {
        this.state.view.overlayMode = this.els.overlayMode.value;
        this.onChange(this.state);
      });
      this.els.showMembrane.addEventListener("change", () => {
        this.state.view.membraneVisible = this.els.showMembrane.checked;
        this.onChange(this.state);
      });
      this.els.showReference.addEventListener("change", () => {
        this.state.view.referenceVisible = this.els.showReference.checked;
        this.onChange(this.state);
      });
      this.els.showDisplacementVectors.addEventListener("change", () => {
        this.state.view.displacementVectorsVisible = this.els.showDisplacementVectors.checked;
        this.onChange(this.state);
      });
      this.els.showGaps.addEventListener("change", () => {
        this.state.view.gapsVisible = this.els.showGaps.checked;
        this.onChange(this.state);
      });
      this.els.showStops.addEventListener("change", () => {
        this.state.view.stopsVisible = this.els.showStops.checked;
        this.onChange(this.state);
      });
      this.els.showPivots.addEventListener("change", () => {
        this.state.view.pivotsVisible = this.els.showPivots.checked;
        this.onChange(this.state);
      });
      this.els.showLinkages.addEventListener("change", () => {
        this.state.view.linkagesVisible = this.els.showLinkages.checked;
        this.onChange(this.state);
      });
      this.els.showFasteners.addEventListener("change", () => {
        this.state.view.fastenersVisible = this.els.showFasteners.checked;
        this.onChange(this.state);
      });
      this.els.showActuators.addEventListener("change", () => {
        this.state.view.actuatorsVisible = this.els.showActuators.checked;
        this.onChange(this.state);
      });
      this.els.actuatorDisplayMode.addEventListener("change", () => {
        this.state.view.actuatorDisplayMode = this.els.actuatorDisplayMode.value;
        this.onChange(this.state);
      });
      this.els.showInfluenceFootprint.addEventListener("change", () => {
        this.state.view.influenceFootprintVisible = this.els.showInfluenceFootprint.checked;
        this.onChange(this.state);
      });
      this.els.showMeasurements.addEventListener("change", () => {
        this.state.view.measurementsVisible = this.els.showMeasurements.checked;
        this.onChange(this.state);
      });
      this.els.showMeasurementLabels.addEventListener("change", () => {
        this.state.view.measurementLabelsVisible = this.els.showMeasurementLabels.checked;
        this.onChange(this.state);
      });
      this.els.showTargetErrorVectors.addEventListener("change", () => {
        this.state.view.targetErrorVectorsVisible = this.els.showTargetErrorVectors.checked;
        this.onChange(this.state);
      });
      this.els.showSurfaceNormals.addEventListener("change", () => {
        this.state.view.surfaceNormalsVisible = this.els.showSurfaceNormals.checked;
        this.onChange(this.state);
      });
      this.els.measurementMode.addEventListener("change", () => {
        this.state.view.measurementMode = this.els.measurementMode.value;
        this.onChange(this.state);
      });
      this.els.showTarget.addEventListener("change", () => {
        this.state.view.targetVisible = this.els.showTarget.checked;
        this.onChange(this.state);
      });
      this.els.surfaceInterpolation.addEventListener("change", () => {
        this.state.view.surfaceInterpolation = this.els.surfaceInterpolation.value;
        this.onChange(this.state);
      });
      this.els.surfaceSubdivisions.addEventListener("input", () => {
        this.state.view.surfaceSubdivisions = Number(this.els.surfaceSubdivisions.value);
        this.onChange(this.state);
      });
      this.els.showSurfaceContours.addEventListener("change", () => {
        this.state.view.surfaceContoursVisible = this.els.showSurfaceContours.checked;
        this.onChange(this.state);
      });
      this.els.contourMode.addEventListener("change", () => {
        this.state.view.contourMode = this.els.contourMode.value;
        if (this.state.view.contourMode !== "membrane") this.state.view.targetVisible = this.state.target.type !== "none";
        this.syncControls();
        this.onChange(this.state);
      });
      this.els.animationResponse.addEventListener("input", () => {
        this.state.view.animationResponse = Number(this.els.animationResponse.value);
        this.updateLabels(null);
      });
      this.els.targetType.addEventListener("change", () => {
        this.state.target.type = this.els.targetType.value;
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.syncControls();
        this.onChange(this.state);
      });
      this.els.targetAmp.addEventListener("input", () => {
        this.state.target.amplitude = Number(this.els.targetAmp.value);
        this.onChange(this.state);
      });
      this.els.targetFreq.addEventListener("input", () => {
        this.state.target.frequency = Number(this.els.targetFreq.value);
        this.onChange(this.state);
      });
      this.els.customTargetExpression.addEventListener("input", () => {
        this.state.target.customExpression = this.els.customTargetExpression.value;
        this.state.target.type = "custom";
        this.state.view.targetVisible = true;
        this.syncControls();
        this.onChange(this.state);
      });
      this.els.optimizerIterations.addEventListener("input", () => {
        this.state.target.optimizerIterations = Number(this.els.optimizerIterations.value);
        this.updateLabels(null);
      });
      this.els.actuatorPenalty.addEventListener("input", () => {
        this.state.target.actuatorPenalty = Number(this.els.actuatorPenalty.value);
        this.updateLabels(null);
      });
      this.els.travelPenalty.addEventListener("input", () => {
        this.state.target.travelPenalty = Number(this.els.travelPenalty.value);
        this.updateLabels(null);
      });
      this.els.maxActuators.addEventListener("input", () => {
        this.state.inverse.maxActuators = Number(this.els.maxActuators.value);
        this.updateLabels(null);
      });
      this.els.zTravelLimit.addEventListener("input", () => this.updateTravelLimits());
      this.els.alphaContractLimit.addEventListener("input", () => this.updateTravelLimits());
      this.els.alphaExpandLimit.addEventListener("input", () => this.updateTravelLimits());
      this.els.timelineIndex.addEventListener("input", () => {
        this.stopTimelinePlayback();
        const duration = this.state.timeline.smooth === false ? 0 : Math.min(600, this.state.timeline.transitionMs || 900);
        this.animateTimelineIndex(Number(this.els.timelineIndex.value), duration);
      });
      this.els.sequenceChart.addEventListener("click", (event) => this.jumpSequenceChart(event));
      this.els.sequenceChart.addEventListener("pointermove", (event) => this.hoverSequenceChart(event));
      this.els.sequenceChart.addEventListener("pointerleave", () => this.clearSequenceChartHover());
      this.els.timelineSpeed.addEventListener("input", () => {
        this.state.timeline.speed = Number(this.els.timelineSpeed.value);
        this.updateLabels(null);
      });
      this.els.transitionMs.addEventListener("input", () => {
        this.state.timeline.transitionMs = Number(this.els.transitionMs.value);
        this.updateLabels(null);
      });
      this.els.smoothTimeline.addEventListener("change", () => {
        this.state.timeline.smooth = this.els.smoothTimeline.checked;
        this.updateLabels(null);
      });
      document.getElementById("applyCommand").addEventListener("click", () => this.applySelected());
      document.getElementById("clearCell").addEventListener("click", () => this.clearSelected());
      this.els.operatorCommitLock.addEventListener("click", () => this.commitSelectedLockEvent());
      this.els.operatorReleaseLock.addEventListener("click", () => this.releaseSelectedLockEvent());
      this.els.operatorCheckOrder.addEventListener("click", () => this.checkSelectedEventOrder());
      document.getElementById("nudgeAlphaContract").addEventListener("click", () => this.nudgeSelectedCommand(-0.08, 0));
      document.getElementById("nudgeAlphaExpand").addEventListener("click", () => this.nudgeSelectedCommand(0.08, 0));
      document.getElementById("nudgeZDown").addEventListener("click", () => this.nudgeSelectedCommand(0, -0.08));
      document.getElementById("nudgeZUp").addEventListener("click", () => this.nudgeSelectedCommand(0, 0.08));
      document.getElementById("zeroSelectedCommand").addEventListener("click", () => this.setSelectedCommand(0, 0));
      document.getElementById("allowAllActuators").addEventListener("click", () => this.applyActuatorMaskPreset("all"));
      document.getElementById("blockRimActuators").addEventListener("click", () => this.applyActuatorMaskPreset("block-rim"));
      document.getElementById("allowCenterActuators").addEventListener("click", () => this.applyActuatorMaskPreset("center"));
      document.getElementById("invertActuatorMask").addEventListener("click", () => this.applyActuatorMaskPreset("invert"));
      document.getElementById("analyzeInversePlan").addEventListener("click", () => {
        const plan = RAD.buildInverseDesignPlan(this.state);
        RAD.recordEvent(this.state, {
          type: "inverse-plan-analyzed",
          target: this.state.target.type,
          actuators: plan.commands.length,
          score: plan.baseScore,
        });
        this.state.view.overlayMode = "inverse";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("analyzeSensitivity").addEventListener("click", () => {
        const sensitivity = RAD.analyzeActuatorSensitivity(this.state);
        RAD.recordEvent(this.state, {
          type: "sensitivity-analyzed",
          cells: sensitivity.controllableCells,
          meanGain: sensitivity.meanGain,
          maxGain: sensitivity.maxGain,
        });
        this.state.view.overlayMode = "sensitivity";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("buildJacobian").addEventListener("click", () => {
        const jacobian = RAD.buildResponseJacobian(this.state);
        RAD.recordEvent(this.state, {
          type: "jacobian-built",
          columns: jacobian.columnCount,
          actuators: jacobian.actuatorCount,
          meanCoverage: jacobian.meanCoverage,
          conditionEstimate: jacobian.conditionEstimate,
        });
        this.state.view.overlayMode = "reachability";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("solveLinearFit").addEventListener("click", () => {
        const solution = RAD.solveLinearizedTargetFit(this.state);
        RAD.recordEvent(this.state, {
          type: "linear-fit-solved",
          actuators: solution.commands.length,
          steps: solution.steps,
          baseError: solution.baseError,
          projectedError: solution.projectedError,
        });
        this.state.view.overlayMode = "inverse";
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("applyInversePlan").addEventListener("click", () => {
        RAD.applyInverseDesignPlan(this.state);
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.state.view.overlayMode = "inverse";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("applyLinearFit").addEventListener("click", () => {
        RAD.applyLinearizedTargetFit(this.state);
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.state.view.overlayMode = "inverse";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("seedTargetFit").addEventListener("click", () => {
        RAD.estimateCommandsForTarget(this.state);
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("optimizeTargetFit").addEventListener("click", () => {
        RAD.optimizeCommandsForTarget(this.state);
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("clearTarget").addEventListener("click", () => {
        this.state.target.type = "none";
        this.state.view.targetVisible = false;
        this.state.inverse.plan = { candidates: [], commands: [], history: [] };
        this.state.inverse.preview = null;
        this.state.inverse.lastScore = 0;
        this.state.inverse.sensitivity = { candidates: [], map: RAD.matrix(this.state.grid.rows, this.state.grid.cols, 0), stepZ: 0.12, stepAlpha: 0.12, controllableCells: 0, meanGain: 0, maxGain: 0 };
        this.state.inverse.jacobian = { columns: [], coverageMap: RAD.matrix(this.state.grid.rows, this.state.grid.cols, 0), actuatorCount: 0, columnCount: 0, meanCoverage: 0, maxCoverage: 0, stepZ: 0.12, stepAlpha: 0.12, conditionEstimate: 0 };
        this.state.inverse.linearSolution = { commands: [], history: [], steps: 0, baseError: 0, predictedError: 0, projectedError: 0, projectedActuators: 0 };
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("clearInversePreview").addEventListener("click", () => {
        RAD.setInversePreview(this.state, null);
        this.els.inverseStepPreview.value = 0;
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("captureInversePreview").addEventListener("click", () => this.captureInversePreviewKeyframe());
      this.els.inverseStepPreview.addEventListener("input", () => {
        RAD.setInversePlanStepPreview(this.state, Number(this.els.inverseStepPreview.value));
        this.state.view.overlayMode = "inverse";
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("playTimeline").addEventListener("click", () => this.toggleTimeline());
      document.getElementById("stepTimeline").addEventListener("click", () => this.stepTimeline());
      document.getElementById("captureKeyframe").addEventListener("click", () => this.captureKeyframe());
      document.getElementById("firstKeyframe").addEventListener("click", () => this.jumpKeyframe("first"));
      document.getElementById("prevKeyframe").addEventListener("click", () => this.jumpKeyframe("prev"));
      document.getElementById("nextKeyframe").addEventListener("click", () => this.jumpKeyframe("next"));
      document.getElementById("lastKeyframe").addEventListener("click", () => this.jumpKeyframe("last"));
      document.getElementById("clearTimeline").addEventListener("click", () => {
        this.stopTimelinePlayback();
        this.state.experiment.eventList = [];
        this.state.experiment.initialSnapshot = null;
        this.state.timeline.index = 0;
        this.lastSequenceAnalysis = null;
        this.lastSequenceAnalysisSignature = "";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("saveSequence").addEventListener("click", () => this.saveSequence());
      document.getElementById("loadSequence").addEventListener("click", () => this.els.sequenceFileInput.click());
      this.els.sequenceFileInput.addEventListener("change", () => this.loadSequence());
      document.getElementById("analyzeSequence").addEventListener("click", () => this.analyzeSequence());
      document.getElementById("jumpBestFrame").addEventListener("click", () => this.jumpAnalysisFrame("best-error"));
      document.getElementById("jumpMaxSaturationFrame").addEventListener("click", () => this.jumpAnalysisFrame("max-saturation"));
      document.getElementById("selectWorstResidualCell").addEventListener("click", () => this.selectWorstResidualCell());
      document.getElementById("saveSequenceCsv").addEventListener("click", () => this.saveSequenceCsv());
      document.querySelectorAll("[data-preset]").forEach((button) => {
        button.addEventListener("click", () => {
          RAD.applyPreset(this.state, button.dataset.preset);
          this.syncControls();
          this.onChange(this.state);
        });
      });
      document.getElementById("saveJson").addEventListener("click", () => this.saveJson());
      this.els.saveObj.addEventListener("click", () => this.saveObj());
      document.getElementById("loadJson").addEventListener("click", () => this.els.fileInput.click());
      this.els.fileInput.addEventListener("change", () => this.loadJson());
    }

    setState(state) {
      this.state = state;
      this.syncControls();
    }

    select(r, c) {
      this.state.selection = { r, c };
      this.syncControls();
      this.onChange(this.state);
    }

    paintCell(r, c) {
      this.state.selection = { r, c };
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      const alpha = RAD.clampCommandAlpha(this.state, this.els.alphaCommand.value);
      const z = RAD.clampCommandZ(this.state, this.els.zCommand.value);
      const cells = this.paintBrushCells(r, c);
      for (const cell of cells) {
        this.state.cells.commandAlpha[cell.r][cell.c] = alpha;
        this.state.cells.commandZ[cell.r][cell.c] = z;
        this.state.cells.locked[cell.r][cell.c] = this.els.locked.checked;
        this.state.cells.actuatorAllowed[cell.r][cell.c] = this.els.actuatorAllowed.checked;
      }
      RAD.recordEvent(this.state, {
        type: "cell-paint-command",
        r,
        c,
        alpha,
        z,
        radius: Number(this.state.view.paintRadius || 0),
        affected: cells.length,
        locked: this.els.locked.checked,
        actuatorAllowed: this.els.actuatorAllowed.checked,
      });
      this.syncControls();
      this.onChange(this.state);
    }

    paintBrushCells(r, c) {
      return RAD.brushCells(this.state, r, c, this.state.view.paintRadius);
    }

    applySelected() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      this.state.cells.commandAlpha[r][c] = RAD.clampCommandAlpha(this.state, this.els.alphaCommand.value);
      this.state.cells.commandZ[r][c] = RAD.clampCommandZ(this.state, this.els.zCommand.value);
      this.state.cells.locked[r][c] = this.els.locked.checked;
      this.state.cells.actuatorAllowed[r][c] = this.els.actuatorAllowed.checked;
      RAD.recordEvent(this.state, {
        type: "cell-command",
        r,
        c,
        alpha: this.state.cells.commandAlpha[r][c],
        z: this.state.cells.commandZ[r][c],
        locked: this.state.cells.locked[r][c],
        actuatorAllowed: this.state.cells.actuatorAllowed[r][c],
      });
      this.onChange(this.state);
    }

    previewSelected() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      this.state.cells.commandAlpha[r][c] = RAD.clampCommandAlpha(this.state, this.els.alphaCommand.value);
      this.state.cells.commandZ[r][c] = RAD.clampCommandZ(this.state, this.els.zCommand.value);
      this.state.cells.locked[r][c] = this.els.locked.checked;
      this.onChange(this.state);
    }

    nudgeSelectedCommand(deltaAlpha, deltaZ) {
      const nextAlpha = Number(this.els.alphaCommand.value || 0) + deltaAlpha;
      const nextZ = Number(this.els.zCommand.value || 0) + deltaZ;
      this.setSelectedCommand(nextAlpha, nextZ);
    }

    setSelectedCommand(alpha, z) {
      this.els.alphaCommand.value = RAD.clampCommandAlpha(this.state, alpha);
      this.els.zCommand.value = RAD.clampCommandZ(this.state, z);
      this.previewSelected();
      this.updateLabels(null);
    }

    updateSelectedActuatorMask() {
      const { r, c } = this.state.selection;
      this.state.cells.actuatorAllowed[r][c] = this.els.actuatorAllowed.checked;
      this.state.inverse.plan = { candidates: [], commands: [], history: [] };
      this.state.inverse.preview = null;
      this.state.inverse.sensitivity = { candidates: [], map: RAD.matrix(this.state.grid.rows, this.state.grid.cols, 0), stepZ: 0.12, stepAlpha: 0.12, controllableCells: 0, meanGain: 0, maxGain: 0 };
      this.onChange(this.state);
    }

    updateTravelLimits() {
      this.state.grid.zTravelLimit = Number(this.els.zTravelLimit.value);
      this.state.grid.alphaContractLimit = Number(this.els.alphaContractLimit.value);
      this.state.grid.alphaExpandLimit = Number(this.els.alphaExpandLimit.value);
      RAD.clampAllCommands(this.state);
      this.state.inverse.plan = { candidates: [], commands: [], history: [] };
      this.state.inverse.preview = null;
      this.state.inverse.sensitivity = { candidates: [], map: RAD.matrix(this.state.grid.rows, this.state.grid.cols, 0), stepZ: 0.12, stepAlpha: 0.12, controllableCells: 0, meanGain: 0, maxGain: 0 };
      this.syncControls();
      this.onChange(this.state);
    }

    applyActuatorMaskPreset(name) {
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      const { rows, cols } = this.state.grid;
      const midR = (rows - 1) / 2;
      const midC = (cols - 1) / 2;
      const radius = Math.max(1, Math.min(rows, cols) * 0.32);
      let allowed = 0;
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          let value = true;
          if (name === "block-rim") value = !(r === 0 || c === 0 || r === rows - 1 || c === cols - 1);
          else if (name === "center") value = Math.hypot(r - midR, c - midC) <= radius;
          else if (name === "invert") value = this.state.cells.actuatorAllowed?.[r]?.[c] === false;
          this.state.cells.actuatorAllowed[r][c] = value;
          if (value) allowed += 1;
        }
      }
      this.state.inverse.plan = { candidates: [], commands: [], history: [] };
      this.state.inverse.preview = null;
      this.state.inverse.sensitivity = { candidates: [], map: RAD.matrix(rows, cols, 0), stepZ: 0.12, stepAlpha: 0.12, controllableCells: 0, meanGain: 0, maxGain: 0 };
      RAD.recordEvent(this.state, { type: "actuator-mask", name, allowed });
      this.syncControls();
      this.onChange(this.state);
    }

    clearSelected() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      this.state.cells.commandAlpha[r][c] = 0;
      this.state.cells.commandZ[r][c] = 0;
      this.state.cells.locked[r][c] = false;
      if (this.state.cells.lockAlpha) this.state.cells.lockAlpha[r][c] = this.state.grid.initialAlpha;
      RAD.recordEvent(this.state, { type: "cell-clear", r, c });
      this.syncControls();
      this.onChange(this.state);
    }

    commitSelectedLockEvent() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      this.state.cells.commandAlpha[r][c] = RAD.clampCommandAlpha(this.state, this.els.alphaCommand.value);
      this.state.cells.commandZ[r][c] = RAD.clampCommandZ(this.state, this.els.zCommand.value);
      this.state = RAD.applyProgrammableEvent(this.state, RAD.lockEvent({ r, c }));
      RAD.recordEvent(this.state, {
        type: "operator-lock",
        r,
        c,
        lockAlpha: this.state.cells.lockAlpha?.[r]?.[c] ?? this.state.grid.initialAlpha,
      });
      this.syncControls();
      this.onChange(this.state);
    }

    releaseSelectedLockEvent() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      this.state = RAD.applyProgrammableEvent(this.state, RAD.releaseEvent({ r, c }));
      RAD.recordEvent(this.state, { type: "operator-release", r, c });
      this.syncControls();
      this.onChange(this.state);
    }

    selectedEventBaseline(r, c) {
      const base = JSON.parse(JSON.stringify(this.state));
      if (!base.cells.lockAlpha) {
        base.cells.lockAlpha = RAD.matrix(base.grid.rows, base.grid.cols, base.grid.initialAlpha);
      }
      base.cells.commandAlpha[r][c] = 0;
      base.cells.commandZ[r][c] = 0;
      base.cells.locked[r][c] = false;
      base.cells.lockAlpha[r][c] = base.grid.initialAlpha;
      return base;
    }

    checkSelectedEventOrder() {
      const { r, c } = this.state.selection;
      const alpha = RAD.clampCommandAlpha(this.state, this.els.alphaCommand.value);
      const z = RAD.clampCommandZ(this.state, this.els.zCommand.value);
      const diagnostic = RAD.compareEventOrder(
        this.selectedEventBaseline(r, c),
        RAD.localActuationEvent({ r, c }, alpha, z),
        RAD.lockEvent({ r, c })
      );
      this.lastOperatorDiagnostic = diagnostic;
      RAD.recordEvent(this.state, {
        type: "operator-order-check",
        r,
        c,
        alpha,
        z,
        finalAlphaError: diagnostic.finalAlphaError,
        finalHeightError: diagnostic.finalHeightError,
      });
      this.updateOperatorInspector(diagnostic);
      this.onChange(this.state);
    }

    updateOperatorInspector(diagnostic) {
      if (!diagnostic) {
        this.els.operatorOrderState.textContent = "not checked";
        this.els.operatorAlphaError.textContent = "0.000";
        this.els.operatorHeightError.textContent = "0.000";
        return;
      }
      const commutes = diagnostic.modeCommutes && diagnostic.commandCommutes && diagnostic.lockAlphaCommutes && diagnostic.finalAlphaError < 1e-9 && diagnostic.finalHeightError < 1e-9;
      this.els.operatorOrderState.textContent = commutes ? "commutes" : "path dependent";
      this.els.operatorAlphaError.textContent = diagnostic.finalAlphaError.toFixed(3);
      this.els.operatorHeightError.textContent = diagnostic.finalHeightError.toFixed(3);
      this.els.operatorOrderState.title = `mode ${diagnostic.modeCommutes ? "same" : "diff"}, commands ${diagnostic.commandCommutes ? "same" : "diff"}, lock alpha ${diagnostic.lockAlphaCommutes ? "same" : "diff"}`;
    }

    applyCellVisualMode(mode) {
      const visualMode = ["abstract", "paperRad", "mechanism"].includes(mode) ? mode : "abstract";
      this.state.view.cellVisualMode = visualMode;
      if (visualMode === "mechanism") {
        Object.assign(this.state.view, {
          gapsVisible: true,
          stopsVisible: true,
          pivotsVisible: true,
          linkagesVisible: true,
          fastenersVisible: true,
          actuatorsVisible: true,
          measurementsVisible: true,
          measurementMode: "all",
        });
      } else if (visualMode === "paperRad") {
        this.state.view.explodedSelected = false;
        Object.assign(this.state.view, {
          gapsVisible: true,
          stopsVisible: true,
          pivotsVisible: true,
          linkagesVisible: true,
          fastenersVisible: false,
          actuatorsVisible: true,
          measurementsVisible: true,
          measurementMode: "backlash",
        });
      } else {
        this.state.view.explodedSelected = false;
        Object.assign(this.state.view, {
          gapsVisible: true,
          stopsVisible: false,
          pivotsVisible: false,
          linkagesVisible: false,
          fastenersVisible: false,
          actuatorsVisible: false,
          measurementsVisible: true,
          measurementMode: "alpha",
        });
      }
      this.syncControls();
    }

    syncControls() {
      const s = this.state;
      this.els.rows.value = s.grid.rows;
      this.els.cols.value = s.grid.cols;
      this.els.backlash.value = s.grid.backlash;
      this.els.coupling.value = s.grid.couplingGain;
      this.els.zCoupling.value = s.grid.zCouplingGain ?? 0.32;
      this.els.pinRadius.value = s.grid.pinRadius ?? 0.18;
      this.els.holeRadius.value = s.grid.holeRadius ?? 0.225;
      this.els.paintRadius.value = Math.max(0, Math.min(2, Number(s.view.paintRadius || 0)));
      this.els.cellVisualMode.value = s.view.cellVisualMode || "abstract";
      this.els.simulationMode.value = s.view.simulationMode || "kinematic";
      this.els.overlayMode.value = s.view.overlayMode;
      this.els.showMembrane.checked = s.view.membraneVisible;
      this.els.showReference.checked = s.view.referenceVisible !== false;
      this.els.showDisplacementVectors.checked = s.view.displacementVectorsVisible !== false;
      this.els.showGaps.checked = s.view.gapsVisible;
      this.els.showStops.checked = s.view.stopsVisible !== false;
      this.els.showPivots.checked = s.view.pivotsVisible !== false;
      this.els.showLinkages.checked = s.view.linkagesVisible !== false;
      this.els.showFasteners.checked = s.view.fastenersVisible !== false;
      this.els.showActuators.checked = s.view.actuatorsVisible;
      this.els.actuatorDisplayMode.value = s.view.actuatorDisplayMode || "selected";
      this.els.showInfluenceFootprint.checked = s.view.influenceFootprintVisible !== false;
      this.els.showMeasurements.checked = s.view.measurementsVisible;
      this.els.showMeasurementLabels.checked = s.view.measurementLabelsVisible === true;
      this.els.showTargetErrorVectors.checked = s.view.targetErrorVectorsVisible !== false;
      this.els.showSurfaceNormals.checked = s.view.surfaceNormalsVisible === true;
      this.els.measurementMode.value = s.view.measurementMode || "all";
      this.els.showTarget.checked = s.view.targetVisible;
      this.els.surfaceInterpolation.value = s.view.surfaceInterpolation || "smooth";
      this.els.surfaceSubdivisions.value = s.view.surfaceSubdivisions || 5;
      this.els.showSurfaceContours.checked = s.view.surfaceContoursVisible === true;
      this.els.contourMode.value = s.view.contourMode || "membrane";
      this.els.animationResponse.value = s.view.animationResponse || 0.18;
      this.els.targetType.value = s.target.type;
      this.els.targetAmp.value = s.target.amplitude;
      this.els.targetFreq.value = s.target.frequency;
      this.els.customTargetExpression.value = s.target.customExpression || "";
      this.els.optimizerIterations.value = s.target.optimizerIterations;
      this.els.actuatorPenalty.value = s.target.actuatorPenalty;
      this.els.travelPenalty.value = s.target.travelPenalty;
      this.els.maxActuators.value = s.inverse?.maxActuators || 24;
      this.els.zTravelLimit.value = s.grid.zTravelLimit || 0.8;
      this.els.alphaContractLimit.value = s.grid.alphaContractLimit || 0.55;
      this.els.alphaExpandLimit.value = s.grid.alphaExpandLimit || 0.35;
      this.els.timelineIndex.max = Math.max(0, s.experiment.eventList.length);
      this.els.timelineIndex.value = Math.min(s.timeline.index, s.experiment.eventList.length);
      this.els.timelineSpeed.value = s.timeline.speed;
      this.els.transitionMs.value = s.timeline.transitionMs || 900;
      this.els.smoothTimeline.checked = s.timeline.smooth !== false;
      this.els.keyframeName.value = s.timeline.keyframeName || "pose";
      const { r, c } = s.selection;
      const limits = RAD.commandLimits(s);
      this.els.alphaCommand.min = -limits.alphaContract;
      this.els.alphaCommand.max = limits.alphaExpand;
      this.els.zCommand.min = -limits.z;
      this.els.zCommand.max = limits.z;
      this.els.alphaCommand.value = s.cells.commandAlpha[r][c] ?? 0;
      this.els.zCommand.value = s.cells.commandZ[r][c] ?? 0;
      this.els.locked.checked = s.cells.locked[r][c];
      this.els.actuatorAllowed.checked = s.cells.actuatorAllowed?.[r]?.[c] !== false;
      this.updateLabels(null);
      this.syncQuickDock();
      this.syncPaintMode();
    }

    syncQuickDock() {
      const collapsed = this.state.view.quickDockCollapsed === true;
      this.els.quickDock.classList.toggle("is-collapsed", collapsed);
      this.els.toggleQuickDock.textContent = collapsed ? "Show" : "Hide";
      this.els.toggleQuickDock.setAttribute("aria-pressed", String(collapsed));
    }

    syncPaintMode() {
      const enabled = this.state.view.paintMode === true;
      this.els.paintMode.classList.toggle("is-active", enabled);
      this.els.paintMode.setAttribute("aria-pressed", String(enabled));
      this.els.paintMode.textContent = enabled ? "Paint Clicks On" : "Paint Clicks Off";
    }

    updateLabels(sim) {
      const s = this.state;
      const { r, c } = s.selection;
      document.getElementById("rowsOut").textContent = s.grid.rows;
      document.getElementById("colsOut").textContent = s.grid.cols;
      document.getElementById("backlashOut").textContent = Number(s.grid.backlash).toFixed(2);
      document.getElementById("couplingOut").textContent = Number(s.grid.couplingGain).toFixed(2);
      document.getElementById("zCouplingOut").textContent = Number(s.grid.zCouplingGain ?? 0.32).toFixed(2);
      document.getElementById("pinRadiusOut").textContent = Number(s.grid.pinRadius ?? 0.18).toFixed(3);
      document.getElementById("holeRadiusOut").textContent = Number(s.grid.holeRadius ?? 0.225).toFixed(3);
      document.getElementById("pinHoleClearanceOut").textContent = RAD.pinHoleClearance(s).toFixed(3);
      document.getElementById("alphaCommandOut").textContent = Number(this.els.alphaCommand.value).toFixed(2);
      document.getElementById("zCommandOut").textContent = Number(this.els.zCommand.value).toFixed(2);
      document.getElementById("targetAmpOut").textContent = Number(s.target.amplitude).toFixed(2);
      document.getElementById("targetFreqOut").textContent = Number(s.target.frequency).toFixed(2);
      this.els.targetExpressionStatus.textContent = s.target.expressionError || "Variables: nr, nc, r, c, d, amplitude, frequency";
      this.els.targetExpressionStatus.classList.toggle("is-error", Boolean(s.target.expressionError));
      document.getElementById("surfaceSubdivisionsOut").textContent = String(s.view.surfaceSubdivisions || 5);
      document.getElementById("animationResponseOut").textContent = Number(s.view.animationResponse || 0.18).toFixed(2);
      document.getElementById("paintRadiusOut").textContent = String(Math.max(0, Math.min(2, Number(s.view.paintRadius || 0))));
      document.getElementById("optimizerIterationsOut").textContent = String(s.target.optimizerIterations);
      document.getElementById("actuatorPenaltyOut").textContent = Number(s.target.actuatorPenalty).toFixed(3);
      document.getElementById("travelPenaltyOut").textContent = Number(s.target.travelPenalty).toFixed(3);
      document.getElementById("maxActuatorsOut").textContent = String(s.inverse?.maxActuators || 24);
      document.getElementById("zTravelLimitOut").textContent = Number(s.grid.zTravelLimit || 0.8).toFixed(2);
      document.getElementById("alphaContractLimitOut").textContent = Number(s.grid.alphaContractLimit || 0.55).toFixed(2);
      document.getElementById("alphaExpandLimitOut").textContent = Number(s.grid.alphaExpandLimit || 0.35).toFixed(2);
      document.getElementById("timelineOut").textContent = `${Math.min(s.timeline.index, s.experiment.eventList.length)} / ${s.experiment.eventList.length}`;
      document.getElementById("timelineSpeedOut").textContent = `${Number(s.timeline.speed).toFixed(1)}x`;
      document.getElementById("transitionMsOut").textContent = `${Number(s.timeline.transitionMs || 900).toFixed(0)} ms`;
      document.getElementById("sequenceSummary").textContent = this.sequenceSummaryText();
      document.getElementById("sequenceAnalysisSummary").textContent = this.sequenceAnalysisText();
      document.getElementById("sequenceFrameDetail").textContent = this.sequenceFrameDetailText();
      this.renderSequenceChart();
      document.getElementById("selectedCell").textContent = `row ${r}, col ${c}`;
      this.els.quickDockMini.textContent = `r${r} c${c} / a ${Number(this.els.alphaCommand.value).toFixed(2)} / z ${Number(this.els.zCommand.value).toFixed(2)}`;
      this.updateCouplingInspector(r, c);
      if (sim) {
        const residual = sim.target[r][c] - sim.height[r][c];
        const localStrain = this.localLinkStrain(sim, r, c);
        const zResidual = sim.zResidual?.[r]?.[c] || 0;
        document.getElementById("selectedValues").textContent = `alpha ${sim.alpha[r][c].toFixed(3)}, theta ${sim.theta[r][c].toFixed(1)}, z ${sim.height[r][c].toFixed(3)}, z residual ${zResidual.toFixed(3)}, slope ${sim.slope.magnitude[r][c].toFixed(3)}, normal tilt ${sim.slope.tilt[r][c].toFixed(1)}, residual ${residual.toFixed(3)}, link strain ${localStrain.toFixed(3)}, ref disp ${sim.displacement[r][c].toFixed(3)}`;
        document.getElementById("meanAlpha").textContent = sim.metrics.meanAlpha.toFixed(3);
        document.getElementById("maxHeight").textContent = sim.metrics.maxAbsHeight.toFixed(3);
        document.getElementById("dieOff").textContent = String(sim.metrics.maxFiniteDieOff);
        document.getElementById("activeCells").textContent = String(sim.metrics.activeCells);
        document.getElementById("targetError").textContent = sim.metrics.rmsTargetError.toFixed(3);
        document.getElementById("signedTargetError").textContent = sim.metrics.meanSignedTargetError.toFixed(3);
        document.getElementById("targetErrorRange").textContent = `${sim.metrics.maxNegativeTargetError.toFixed(2)}:${sim.metrics.maxPositiveTargetError.toFixed(2)}`;
        document.getElementById("thetaRange").textContent = `${sim.metrics.minTheta.toFixed(1)}-${sim.metrics.maxTheta.toFixed(1)}`;
        document.getElementById("recommendedActuators").textContent = String(sim.metrics.recommendedActuators);
        document.getElementById("candidateActuators").textContent = String(s.inverse?.plan?.commands?.length || 0);
        document.getElementById("designScore").textContent = Number(s.inverse?.plan?.baseScore || 0).toFixed(3);
        document.getElementById("projectedScore").textContent = Number(s.inverse?.plan?.projectedScore || 0).toFixed(3);
        document.getElementById("plannerSteps").textContent = String(s.inverse?.plan?.steps || 0);
        document.getElementById("meanTravel").textContent = sim.metrics.meanTravel.toFixed(3);
        document.getElementById("maxSaturation").textContent = sim.metrics.maxSaturation.toFixed(3);
        document.getElementById("saturatedActuators").textContent = String(sim.metrics.saturatedActuators);
        document.getElementById("maxLinkStrain").textContent = sim.metrics.maxAbsLinkStrain.toFixed(3);
        document.getElementById("meanLinkStrain").textContent = sim.metrics.meanAbsLinkStrain.toFixed(3);
        document.getElementById("maxReferenceDisplacement").textContent = sim.metrics.maxReferenceDisplacement.toFixed(3);
        document.getElementById("maxSurfaceSlope").textContent = sim.metrics.maxSurfaceSlope.toFixed(3);
        document.getElementById("meanSurfaceSlope").textContent = sim.metrics.meanSurfaceSlope.toFixed(3);
        document.getElementById("maxNormalTilt").textContent = sim.metrics.maxNormalTilt.toFixed(1);
        document.getElementById("meanNormalTilt").textContent = sim.metrics.meanNormalTilt.toFixed(1);
        document.getElementById("sensitivityCells").textContent = String(s.inverse?.sensitivity?.controllableCells || 0);
        document.getElementById("meanSensitivity").textContent = Number(s.inverse?.sensitivity?.meanGain || 0).toFixed(3);
        document.getElementById("jacobianColumns").textContent = String(s.inverse?.jacobian?.columnCount || 0);
        document.getElementById("meanReachability").textContent = Number(s.inverse?.jacobian?.meanCoverage || 0).toFixed(3);
        document.getElementById("linearFitSteps").textContent = String(s.inverse?.linearSolution?.steps || 0);
        document.getElementById("linearFitError").textContent = Number(s.inverse?.linearSolution?.projectedError || 0).toFixed(3);
      }
      this.renderInversePlan();
      this.updateOperatorInspector(this.lastOperatorDiagnostic);
      this.renderEvents();
    }

    updateCouplingInspector(r, c) {
      if (typeof RAD.selectedCellCouplingMetrics !== "function") return;
      const metrics = RAD.selectedCellCouplingMetrics(this.state, r, c);
      document.getElementById("alphaNeighborSignal").textContent = metrics.alphaNeighborSignal.toFixed(3);
      document.getElementById("zNeighborSignal").textContent = metrics.zNeighborSignal.toFixed(3);
      document.getElementById("alphaReachCells").textContent = String(metrics.alphaReachCells);
      document.getElementById("zReachCells").textContent = String(metrics.zReachCells);
      document.getElementById("alphaNeighborSignal").title = `dead-zone +/-${metrics.alphaDeadZone.toFixed(3)}, max neighbor ${metrics.maxAlphaNeighbor.toFixed(3)}`;
      document.getElementById("zNeighborSignal").title = `dead-zone +/-${metrics.zDeadZone.toFixed(3)}, max neighbor ${metrics.maxZNeighbor.toFixed(3)}`;
      this.updateDeadZoneBadge(document.getElementById("alphaDeadZoneState"), metrics.alphaCoupled, "alpha");
      this.updateDeadZoneBadge(document.getElementById("zDeadZoneState"), metrics.zCoupled, "z");
    }

    updateDeadZoneBadge(element, coupled, label) {
      element.textContent = coupled ? `${label} coupled` : `${label} free gap`;
      element.classList.toggle("is-coupled", coupled);
      element.classList.toggle("is-gap", !coupled);
    }

    localLinkStrain(sim, r, c) {
      const values = [];
      if (sim.linkStrain?.horizontal?.[r]?.[c] !== undefined) values.push(Math.abs(sim.linkStrain.horizontal[r][c]));
      if (sim.linkStrain?.horizontal?.[r]?.[c - 1] !== undefined) values.push(Math.abs(sim.linkStrain.horizontal[r][c - 1]));
      if (sim.linkStrain?.vertical?.[r]?.[c] !== undefined) values.push(Math.abs(sim.linkStrain.vertical[r][c]));
      if (sim.linkStrain?.vertical?.[r - 1]?.[c] !== undefined) values.push(Math.abs(sim.linkStrain.vertical[r - 1][c]));
      return values.length ? Math.max(...values) : 0;
    }

    updateRendererDiagnostics(diagnostics) {
      if (!diagnostics) return;
      document.getElementById("fps").textContent = String(diagnostics.fps);
      document.getElementById("drawCalls").textContent = String(diagnostics.drawCalls);
      document.getElementById("triangles").textContent = String(diagnostics.triangles);
    }

    renderGrid(sim) {
      const s = this.state;
      this.els.cellGrid.style.gridTemplateColumns = `repeat(${s.grid.cols}, minmax(18px, 1fr))`;
      this.els.cellGrid.innerHTML = "";
      for (let r = s.grid.rows - 1; r >= 0; r -= 1) {
        for (let c = 0; c < s.grid.cols; c += 1) {
          const button = document.createElement("button");
          button.type = "button";
          button.className = "cell";
          if (s.selection.r === r && s.selection.c === c) button.classList.add("is-selected");
          if (Math.abs(s.cells.commandAlpha[r][c]) > 1e-9 || Math.abs(s.cells.commandZ[r][c]) > 1e-9) button.classList.add("is-actuated");
          if (s.cells.locked[r][c]) button.classList.add("is-locked");
          if (s.cells.actuatorAllowed?.[r]?.[c] === false) button.classList.add("is-disallowed");
          button.title = `row ${r}, col ${c}, alpha ${sim.alpha[r][c].toFixed(3)}, z ${sim.height[r][c].toFixed(3)}`;
          button.addEventListener("click", () => this.select(r, c));
          this.els.cellGrid.appendChild(button);
        }
      }
    }

    renderInversePlan() {
      const plan = this.state.inverse?.plan;
      const commands = plan?.commands || [];
      const candidates = plan?.candidates || [];
      const linearCommands = this.state.inverse?.linearSolution?.commands || [];
      const sensitivityCandidates = this.state.inverse?.sensitivity?.candidates || [];
      const preview = this.state.inverse?.preview;
      this.els.inverseStepPreview.max = commands.length;
      this.els.inverseStepPreview.value = preview?.type === "plan-step" ? preview.step : 0;
      document.getElementById("inverseStepPreviewOut").textContent = `${this.els.inverseStepPreview.value} / ${commands.length}`;
      document.getElementById("inversePlanSummary").textContent = commands.length
        ? `${commands.length} actuators, error ${Number(plan.projectedError || 0).toFixed(3)}`
        : "No plan analyzed";
      document.getElementById("inversePlanDelta").textContent = `delta ${Number(plan?.totalImprovement || 0).toFixed(3)}`;
      document.getElementById("inversePreviewSummary").textContent = preview
        ? this.formatInversePreview(preview)
        : linearCommands.length
          ? `linear fit: ${linearCommands.length} actuators, err ${Number(this.state.inverse.linearSolution.projectedError || 0).toFixed(3)}`
          : "No actuator preview";
      this.els.inversePlanHistory.innerHTML = "";
      for (const step of (plan?.history || []).slice(-6)) {
        const item = document.createElement("div");
        item.className = "history-row";
        item.textContent = `${step.step}: score ${step.score.toFixed(3)}, err ${step.rmsError.toFixed(3)}, gain ${step.improvement.toFixed(3)}`;
        this.els.inversePlanHistory.appendChild(item);
      }
      this.els.inverseCandidateList.innerHTML = "";
      const rows = commands.length ? commands.slice(0, 10) : linearCommands.length ? linearCommands.slice(0, 10) : candidates.length ? candidates.slice(0, 10) : sensitivityCandidates.slice(0, 10);
      for (const candidate of rows) {
        const button = document.createElement("button");
        button.type = "button";
        button.className = "candidate-row";
        if (candidate.step) button.classList.add("is-accepted");
        if (this.state.selection.r === candidate.r && this.state.selection.c === candidate.c) button.classList.add("is-selected");
        if (this.isPreviewedCandidate(preview, candidate)) button.classList.add("is-previewed");
        const label = candidate.step ? `${candidate.step}. r${candidate.r}, c${candidate.c}` : `r${candidate.r}, c${candidate.c}`;
        const gain = candidate.combinedGain !== undefined ? candidate.combinedGain : candidate.improvement || 0;
        const reach = candidate.reach !== undefined ? `reach ${candidate.reach}` : candidate.predictedGain !== undefined ? `gain ${candidate.predictedGain.toFixed(3)}` : `gain ${gain.toFixed(3)}`;
        button.innerHTML = `<span>${label}</span><span>dz ${candidate.commandZ.toFixed(2)}</span><span>da ${candidate.commandAlpha.toFixed(2)}</span><span>${reach}</span>`;
        button.addEventListener("click", () => this.previewInverseCandidate(candidate));
        this.els.inverseCandidateList.appendChild(button);
      }
    }

    isPreviewedCandidate(preview, candidate) {
      if (!preview) return false;
      if (preview.type === "plan-step") {
        return (preview.commands || []).some((command) => command.r === candidate.r && command.c === candidate.c);
      }
      return preview.r === candidate.r && preview.c === candidate.c;
    }

    formatInversePreview(preview) {
      if (preview.type === "plan-step") {
        return `preview first ${preview.step}: err ${preview.previewError.toFixed(3)}, d ${preview.errorDelta.toFixed(3)}, max dz ${preview.maxContribution.toFixed(3)}`;
      }
      return `preview r${preview.r}, c${preview.c}: err ${preview.previewError.toFixed(3)}, d ${preview.errorDelta.toFixed(3)}, max dz ${preview.maxContribution.toFixed(3)}`;
    }

    previewInverseCandidate(candidate) {
      this.state.selection = { r: candidate.r, c: candidate.c };
      RAD.setInversePreview(this.state, candidate);
      this.state.view.overlayMode = "inverse";
      this.state.view.targetVisible = this.state.target.type !== "none";
      this.syncControls();
      this.onChange(this.state);
    }

    saveJson() {
      const blob = new Blob([RAD.serialize(this.state)], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-state.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveObj() {
      const mesh = RAD.buildPaperRadMesh(this.state, { sim: RAD.simulateActive(this.state) });
      const blob = new Blob([RAD.exportPaperRadMeshObj(mesh)], { type: "text/plain" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-paper-rad.obj";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveSequence() {
      const blob = new Blob([RAD.exportExperimentSequence(this.state)], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-sequence.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    analyzeSequence() {
      this.lastSequenceAnalysis = RAD.analyzeExperimentSequence(this.state);
      this.lastSequenceAnalysisSignature = this.sequenceAnalysisSignature();
      this.updateLabels(null);
    }

    saveSequenceCsv() {
      const blob = new Blob([RAD.exportSequenceMetricsCsv(this.state)], { type: "text/csv" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-sequence-metrics.csv";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
      this.lastSequenceAnalysis = RAD.analyzeExperimentSequence(this.state);
      this.lastSequenceAnalysisSignature = this.sequenceAnalysisSignature();
      this.updateLabels(null);
    }

    jumpAnalysisFrame(kind) {
      const analysis = this.currentSequenceAnalysis();
      if (!analysis.frames.length) return;
      const frame =
        kind === "max-saturation"
          ? analysis.frames.reduce((best, candidate) => (candidate.maxSaturation > best.maxSaturation ? candidate : best), analysis.frames[0])
          : analysis.bestFrame;
      this.stopTimelinePlayback();
      this.hoveredSequenceFrame = null;
      const duration = this.state.timeline.smooth === false ? 0 : Math.min(600, this.state.timeline.transitionMs || 900);
      this.animateTimelineIndex(frame.index, duration);
    }

    analysisFrameForCurrentView() {
      const analysis = this.currentSequenceAnalysis();
      const targetIndex = this.hoveredSequenceFrame === null ? this.state.timeline.index || 0 : this.hoveredSequenceFrame;
      const index = Math.max(0, Math.min(analysis.frames.length - 1, targetIndex));
      return analysis.frames[index] || null;
    }

    selectWorstResidualCell() {
      const frame = this.analysisFrameForCurrentView();
      const hotspot = frame?.worstResidual;
      if (!hotspot) return;
      this.stopTimelinePlayback();
      this.state.selection = { r: hotspot.row, c: hotspot.col };
      this.state.view.overlayMode = "error";
      this.state.view.targetVisible = this.state.target.type !== "none";
      this.state.view.targetErrorVectorsVisible = true;
      this.syncControls();
      this.onChange(this.state);
    }

    loadSequence() {
      const file = this.els.sequenceFileInput.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        try {
          this.state = RAD.importExperimentSequence(this.state, String(reader.result));
          this.lastSequenceAnalysis = RAD.analyzeExperimentSequence(this.state);
          this.lastSequenceAnalysisSignature = this.sequenceAnalysisSignature();
          this.syncControls();
          this.onChange(this.state);
        } catch (error) {
          window.alert(error.message);
        }
      };
      reader.readAsText(file);
      this.els.sequenceFileInput.value = "";
    }

    loadJson() {
      const file = this.els.fileInput.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        try {
          this.state = RAD.deserialize(String(reader.result));
          this.syncControls();
          this.onChange(this.state);
        } catch (error) {
          window.alert(error.message);
        }
      };
      reader.readAsText(file);
      this.els.fileInput.value = "";
    }

    renderEvents() {
      const events = this.state.experiment.eventList.slice(-8).reverse();
      this.els.eventList.innerHTML = "";
      for (const event of events) {
        const item = document.createElement("li");
        item.textContent = this.formatEvent(event);
        this.els.eventList.appendChild(item);
      }
    }

    sequenceSummaryText() {
      const events = this.state.experiment.eventList || [];
      if (!events.length) return "No sequence frames";
      const keyframes = events.filter((event) => event.type === "keyframe" || event.type === "inverse-preview-keyframe").length;
      const last = events[events.length - 1];
      const active = this.countSnapshotCommands(last.snapshot);
      return `${events.length} frames, ${keyframes} keyframes, final ${active} actuators`;
    }

    sequenceAnalysisSignature() {
      const events = this.state.experiment.eventList || [];
      const last = events.at(-1);
      return [
        this.state.grid.rows,
        this.state.grid.cols,
        this.state.grid.backlash,
        this.state.grid.couplingGain,
        this.state.timeline.index,
        events.length,
        last?.at || "",
        last?.type || "",
        last?.name || last?.target || "",
      ].join("|");
    }

    sequenceAnalysisText() {
      const signature = this.sequenceAnalysisSignature();
      if (!this.lastSequenceAnalysis || this.lastSequenceAnalysisSignature !== signature) {
        this.lastSequenceAnalysis = RAD.analyzeExperimentSequence(this.state);
        this.lastSequenceAnalysisSignature = signature;
      }
      const analysis = this.lastSequenceAnalysis;
      if (!analysis.frameCount) return "Initial frame only";
      const final = analysis.finalFrame;
      const best = analysis.bestFrame;
      return `RMS ${final.rmsTargetError.toFixed(3)}, best frame ${best.index} (${best.rmsTargetError.toFixed(3)}), max sat ${final.maxSaturation.toFixed(2)}, delta ${final.commandDelta.toFixed(3)}`;
    }

    currentSequenceAnalysis() {
      const signature = this.sequenceAnalysisSignature();
      if (!this.lastSequenceAnalysis || this.lastSequenceAnalysisSignature !== signature) {
        this.lastSequenceAnalysis = RAD.analyzeExperimentSequence(this.state);
        this.lastSequenceAnalysisSignature = signature;
      }
      return this.lastSequenceAnalysis;
    }

    sequenceFrameDetailText() {
      const frame = this.analysisFrameForCurrentView();
      if (!frame) return "Frame 0: initial";
      const prefix = this.hoveredSequenceFrame === null ? "Frame" : "Hover";
      const label = frame.name && frame.name !== frame.type ? `${frame.type} ${frame.name}` : frame.type;
      const hotspot = frame.worstResidual;
      const worst = hotspot ? `, worst r${hotspot.row} c${hotspot.col} ${hotspot.value.toFixed(3)}` : "";
      return `${prefix} ${frame.index}: ${label}, RMS ${frame.rmsTargetError.toFixed(3)}, sat ${frame.maxSaturation.toFixed(2)}, delta ${frame.commandDelta.toFixed(3)}, active ${frame.active}${worst}`;
    }

    renderSequenceChart() {
      const canvas = this.els.sequenceChart;
      if (!canvas) return;
      const rect = canvas.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      const width = Math.max(280, Math.round((rect.width || canvas.clientWidth || 320) * dpr));
      const height = Math.max(130, Math.round((rect.height || canvas.clientHeight || 150) * dpr));
      if (canvas.width !== width || canvas.height !== height) {
        canvas.width = width;
        canvas.height = height;
      }
      const ctx = canvas.getContext("2d");
      const analysis = this.currentSequenceAnalysis();
      ctx.clearRect(0, 0, width, height);
      ctx.fillStyle = "#f7fafc";
      ctx.fillRect(0, 0, width, height);
      const pad = { left: 34 * dpr, right: 12 * dpr, top: 14 * dpr, bottom: 24 * dpr };
      const plotW = Math.max(1, width - pad.left - pad.right);
      const plotH = Math.max(1, height - pad.top - pad.bottom);
      this.sequenceChartLayout = { pad, plotW, plotH, width, height, dpr, frameCount: analysis.frames.length };
      ctx.strokeStyle = "#dbe3ec";
      ctx.lineWidth = 1 * dpr;
      ctx.beginPath();
      for (let i = 0; i <= 3; i += 1) {
        const y = pad.top + (plotH * i) / 3;
        ctx.moveTo(pad.left, y);
        ctx.lineTo(width - pad.right, y);
      }
      ctx.stroke();
      ctx.fillStyle = "#6b7785";
      ctx.font = `${10 * dpr}px system-ui, sans-serif`;
      ctx.fillText("0", 8 * dpr, height - pad.bottom + 3 * dpr);
      if (!analysis.frames.length) return;
      const fields = [
        { key: "rmsTargetError", color: "#c94848" },
        { key: "maxSaturation", color: "#2f6fcb" },
        { key: "commandDelta", color: "#2c8c5b" },
      ];
      const maxValue = Math.max(
        1e-6,
        ...analysis.frames.flatMap((frame) => fields.map((field) => Number(frame[field.key]) || 0))
      );
      ctx.fillText(maxValue.toFixed(2), 8 * dpr, pad.top + 3 * dpr);
      const point = (frameIndex, value) => {
        const x = pad.left + (analysis.frames.length === 1 ? 0 : (plotW * frameIndex) / (analysis.frames.length - 1));
        const y = pad.top + plotH - (plotH * Math.max(0, Number(value) || 0)) / maxValue;
        return { x, y };
      };
      for (const field of fields) {
        ctx.strokeStyle = field.color;
        ctx.lineWidth = 2 * dpr;
        ctx.beginPath();
        analysis.frames.forEach((frame, i) => {
          const p = point(i, frame[field.key]);
          if (i === 0) ctx.moveTo(p.x, p.y);
          else ctx.lineTo(p.x, p.y);
        });
        ctx.stroke();
      }
      const current = Math.max(0, Math.min(analysis.frames.length - 1, this.state.timeline.index || 0));
      const currentX = point(current, 0).x;
      ctx.strokeStyle = "#252c33";
      ctx.lineWidth = 1.5 * dpr;
      ctx.beginPath();
      ctx.moveTo(currentX, pad.top);
      ctx.lineTo(currentX, pad.top + plotH);
      ctx.stroke();
      if (this.hoveredSequenceFrame !== null) {
        const hoverIndex = Math.max(0, Math.min(analysis.frames.length - 1, this.hoveredSequenceFrame));
        const hoverX = point(hoverIndex, 0).x;
        ctx.strokeStyle = "#f0a027";
        ctx.lineWidth = 1.5 * dpr;
        ctx.setLineDash([4 * dpr, 4 * dpr]);
        ctx.beginPath();
        ctx.moveTo(hoverX, pad.top);
        ctx.lineTo(hoverX, pad.top + plotH);
        ctx.stroke();
        ctx.setLineDash([]);
      }
      const bestIndex = Math.max(0, analysis.frames.findIndex((frame) => frame.index === analysis.bestFrame.index));
      const best = point(bestIndex, analysis.bestFrame.rmsTargetError);
      ctx.fillStyle = "#c94848";
      ctx.beginPath();
      ctx.arc(best.x, best.y, 4 * dpr, 0, Math.PI * 2);
      ctx.fill();
      const maxSaturationFrame = analysis.frames.reduce((bestFrame, frame) => (frame.maxSaturation > bestFrame.maxSaturation ? frame : bestFrame), analysis.frames[0]);
      const maxSatIndex = Math.max(0, analysis.frames.findIndex((frame) => frame.index === maxSaturationFrame.index));
      const maxSat = point(maxSatIndex, maxSaturationFrame.maxSaturation);
      ctx.strokeStyle = "#2f6fcb";
      ctx.lineWidth = 2 * dpr;
      ctx.beginPath();
      ctx.arc(maxSat.x, maxSat.y, 5 * dpr, 0, Math.PI * 2);
      ctx.stroke();
      ctx.fillStyle = "#4f5c69";
      ctx.fillText(`frames ${analysis.frameCount}`, pad.left, height - 7 * dpr);
    }

    sequenceChartFrameFromEvent(event) {
      const canvas = this.els.sequenceChart;
      const analysis = this.currentSequenceAnalysis();
      const layout = this.sequenceChartLayout;
      if (!canvas || !layout || analysis.frames.length <= 1) return null;
      const rect = canvas.getBoundingClientRect();
      const x = (event.clientX - rect.left) * layout.dpr;
      const t = Math.max(0, Math.min(1, (x - layout.pad.left) / layout.plotW));
      return Math.round(t * (analysis.frames.length - 1));
    }

    hoverSequenceChart(event) {
      const frameIndex = this.sequenceChartFrameFromEvent(event);
      if (frameIndex === null || frameIndex === this.hoveredSequenceFrame) return;
      this.hoveredSequenceFrame = frameIndex;
      document.getElementById("sequenceFrameDetail").textContent = this.sequenceFrameDetailText();
      this.renderSequenceChart();
    }

    clearSequenceChartHover() {
      if (this.hoveredSequenceFrame === null) return;
      this.hoveredSequenceFrame = null;
      document.getElementById("sequenceFrameDetail").textContent = this.sequenceFrameDetailText();
      this.renderSequenceChart();
    }

    jumpSequenceChart(event) {
      const analysis = this.currentSequenceAnalysis();
      const frameIndex = this.sequenceChartFrameFromEvent(event);
      if (frameIndex === null) return;
      this.stopTimelinePlayback();
      const target = analysis.frames[frameIndex]?.index ?? frameIndex;
      const duration = this.state.timeline.smooth === false ? 0 : Math.min(600, this.state.timeline.transitionMs || 900);
      this.animateTimelineIndex(target, duration);
    }

    countSnapshotCommands(snapshot) {
      const commandAlpha = snapshot?.cells?.commandAlpha || [];
      const commandZ = snapshot?.cells?.commandZ || [];
      let active = 0;
      for (let r = 0; r < commandAlpha.length; r += 1) {
        for (let c = 0; c < commandAlpha[r].length; c += 1) {
          if (Math.abs(commandAlpha[r][c] || 0) > 1e-9 || Math.abs(commandZ?.[r]?.[c] || 0) > 1e-9) active += 1;
        }
      }
      return active;
    }

    formatEvent(event) {
      if (event.type === "preset") return `preset: ${event.name}`;
      if (event.type === "cell-command") {
        const alpha = typeof event.alpha === "number" ? ` a ${event.alpha.toFixed(2)}` : "";
        const z = typeof event.z === "number" ? ` z ${event.z.toFixed(2)}` : "";
        const locked = event.locked ? " locked" : "";
        return `cell command: r${event.r}, c${event.c}${alpha}${z}${locked}`;
      }
      if (event.type === "actuator-mask") return `actuator mask: ${event.name}, ${event.allowed} allowed`;
      if (event.type === "cell-clear") return `cell clear: r${event.r}, c${event.c}`;
      if (event.type === "target-fit-seed") return `seed fit: ${event.target}`;
      if (event.type === "target-fit-optimized") {
        const score = typeof event.score === "number" ? ` score ${event.score.toFixed(3)}` : "";
        return `optimize: ${event.target}, ${event.iterations} iter${score}`;
      }
      if (event.type === "inverse-plan-analyzed") return `analyze: ${event.target}, ${event.actuators} cells`;
      if (event.type === "inverse-plan-applied") return `apply plan: ${event.target}, ${event.actuators} cells`;
      if (event.type === "sensitivity-analyzed") return `sensitivity: ${event.cells || 0} cells, mean ${Number(event.meanGain || 0).toFixed(3)}`;
      if (event.type === "jacobian-built") return `jacobian: ${event.columns || 0} columns, reach ${Number(event.meanCoverage || 0).toFixed(3)}`;
      if (event.type === "linear-fit-solved") return `linear solve: ${event.actuators || 0} actuators, err ${Number(event.projectedError || 0).toFixed(3)}`;
      if (event.type === "linear-fit-applied") return `linear apply: ${event.actuators || 0} actuators, err ${Number(event.projectedError || 0).toFixed(3)}`;
      if (event.type === "operator-lock") return `event lock: r${event.r}, c${event.c} a ${Number(event.lockAlpha || 0).toFixed(3)}`;
      if (event.type === "operator-release") return `event release: r${event.r}, c${event.c}`;
      if (event.type === "operator-order-check") return `order check: r${event.r}, c${event.c} da ${Number(event.finalAlphaError || 0).toFixed(3)} dz ${Number(event.finalHeightError || 0).toFixed(3)}`;
      if (event.type === "inverse-preview-keyframe") return `preview keyframe: ${event.name || "plan preview"}, ${event.actuators || 0} cells`;
      if (event.type === "keyframe") return `keyframe: ${event.name || "pose"} (${event.index || 0})`;
      return event.name ? `${event.type}: ${event.name}` : event.type;
    }

    captureInversePreviewKeyframe() {
      const preview = this.state.inverse?.preview;
      const commands = preview?.commands || [];
      if (!commands.length) return;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      const locked = JSON.parse(JSON.stringify(this.state.cells.locked));
      RAD.clearCommands(this.state);
      this.state.cells.locked = locked;
      for (const command of commands) {
        this.state.cells.commandZ[command.r][command.c] = command.commandZ;
        this.state.cells.commandAlpha[command.r][command.c] = command.commandAlpha;
      }
      const name = preview.type === "plan-step" ? `plan first ${preview.step}` : `candidate r${preview.r}, c${preview.c}`;
      this.state.timeline.keyframeName = name;
      RAD.recordEvent(this.state, {
        type: "inverse-preview-keyframe",
        name,
        previewType: preview.type || "candidate",
        step: preview.step || 1,
        actuators: commands.length,
        previewError: preview.previewError,
        errorDelta: preview.errorDelta,
      });
      this.syncControls();
      this.onChange(this.state);
    }

    captureKeyframe() {
      const name = (this.els.keyframeName.value || "pose").trim() || "pose";
      this.state.timeline.keyframeName = name;
      RAD.recordEvent(this.state, {
        type: "keyframe",
        name,
        index: this.keyframeEvents().length + 1,
      });
      this.syncControls();
      this.onChange(this.state);
    }

    keyframeEvents() {
      return this.state.experiment.eventList
        .map((event, index) => ({ event, index: index + 1 }))
        .filter((item) => item.event.type === "keyframe");
    }

    jumpKeyframe(direction) {
      const frames = this.keyframeEvents();
      if (!frames.length) return;
      const current = this.state.timeline.index || 0;
      let target = frames[0].index;
      if (direction === "last") target = frames[frames.length - 1].index;
      else if (direction === "prev") {
        const before = frames.filter((frame) => frame.index < current);
        target = (before.at(-1) || frames[frames.length - 1]).index;
      } else if (direction === "next") {
        const after = frames.find((frame) => frame.index > current);
        target = (after || frames[0]).index;
      }
      const duration = this.state.timeline.smooth === false ? 0 : this.state.timeline.transitionMs || 900;
      this.animateTimelineIndex(target, duration);
    }

    toggleTimeline() {
      if (this.timelinePlaying) {
        this.stopTimelinePlayback();
        return;
      }
      this.timelinePlaying = true;
      document.getElementById("playTimeline").textContent = "Pause";
      this.playTimelineLoop();
    }

    stopTimelinePlayback() {
      this.timelinePlaying = false;
      if (this.playTimer) clearTimeout(this.playTimer);
      if (this.transitionFrame) cancelAnimationFrame(this.transitionFrame);
      this.playTimer = null;
      this.transitionFrame = null;
      document.getElementById("playTimeline").textContent = "Play";
    }

    playTimelineLoop() {
      if (!this.timelinePlaying) return;
      const max = this.state.experiment.eventList.length;
      const next = max === 0 ? 0 : (this.state.timeline.index + 1) % (max + 1);
      const duration = this.state.timeline.smooth === false ? 0 : (this.state.timeline.transitionMs || 900) / Math.max(0.25, this.state.timeline.speed || 1);
      this.animateTimelineIndex(next, duration, () => {
        if (!this.timelinePlaying) return;
        this.playTimer = setTimeout(() => this.playTimelineLoop(), Math.max(120, 220 / Math.max(0.25, this.state.timeline.speed || 1)));
      });
    }

    stepTimeline() {
      this.stopTimelinePlayback();
      const max = this.state.experiment.eventList.length;
      const next = max === 0 ? 0 : (this.state.timeline.index + 1) % (max + 1);
      const duration = this.state.timeline.smooth === false ? 0 : this.state.timeline.transitionMs || 900;
      this.animateTimelineIndex(next, duration);
    }

    applyTimelineIndex(index) {
      const max = this.state.experiment.eventList.length;
      this.state.timeline.index = Math.max(0, Math.min(max, index));
      if (this.state.timeline.index > 0) {
        const event = this.state.experiment.eventList[this.state.timeline.index - 1];
        RAD.restoreSnapshot(this.state, event.snapshot);
      } else if (this.state.experiment.initialSnapshot) {
        RAD.restoreSnapshot(this.state, this.state.experiment.initialSnapshot);
      }
      this.state.timeline.index = Math.max(0, Math.min(max, index));
    }

    snapshotForTimelineIndex(index) {
      const max = this.state.experiment.eventList.length;
      const clamped = Math.max(0, Math.min(max, index));
      if (clamped > 0) return this.state.experiment.eventList[clamped - 1]?.snapshot;
      return this.state.experiment.initialSnapshot || RAD.snapshotState(this.state);
    }

    animateTimelineIndex(index, duration, done = null) {
      if (this.transitionFrame) cancelAnimationFrame(this.transitionFrame);
      const max = this.state.experiment.eventList.length;
      const clamped = Math.max(0, Math.min(max, index));
      const start = RAD.snapshotState(this.state);
      const end = this.snapshotForTimelineIndex(clamped);
      const startTimeline = { ...this.state.timeline };
      const startTime = performance.now();
      const run = (now) => {
        const t = duration <= 0 ? 1 : Math.min(1, (now - startTime) / duration);
        const eased = t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
        RAD.restoreSnapshot(this.state, RAD.interpolateSnapshots(start, end, eased));
        this.state.timeline = { ...startTimeline, index: clamped };
        this.syncControls();
        this.onChange(this.state);
        if (t < 1) {
          this.transitionFrame = requestAnimationFrame(run);
        } else {
          this.transitionFrame = null;
          done?.();
        }
      };
      run(startTime);
    }
  }

  RAD.RadUI = RadUI;
})();
