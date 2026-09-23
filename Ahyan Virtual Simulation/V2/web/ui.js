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
        paperSideLengthMm: document.getElementById("paperSideLengthMm"),
        paperHoleToleranceMm: document.getElementById("paperHoleToleranceMm"),
        hardwareProfileName: document.getElementById("hardwareProfileName"),
        hardwareProfileSource: document.getElementById("hardwareProfileSource"),
        hardwareBacklashMm: document.getElementById("hardwareBacklashMm"),
        hardwarePinRadiusMm: document.getElementById("hardwarePinRadiusMm"),
        hardwareHoleRadiusMm: document.getElementById("hardwareHoleRadiusMm"),
        hardwarePlateThicknessMm: document.getElementById("hardwarePlateThicknessMm"),
        hardwareJointStackHeightMm: document.getElementById("hardwareJointStackHeightMm"),
        hardwareBossRadiusMm: document.getElementById("hardwareBossRadiusMm"),
        applyHardwareProfile: document.getElementById("applyHardwareProfile"),
        applyA360CadProfile: document.getElementById("applyA360CadProfile"),
        clearHardwareProfile: document.getElementById("clearHardwareProfile"),
        saveHardwareProfile: document.getElementById("saveHardwareProfile"),
        loadHardwareProfile: document.getElementById("loadHardwareProfile"),
        hardwareProfileFileInput: document.getElementById("hardwareProfileFileInput"),
        saveCalibrationPlan: document.getElementById("saveCalibrationPlan"),
        boundaryMode: document.getElementById("boundaryMode"),
        boundaryWallType: document.getElementById("boundaryWallType"),
        boundaryXMin: document.getElementById("boundaryXMin"),
        boundaryXMax: document.getElementById("boundaryXMax"),
        boundaryYMin: document.getElementById("boundaryYMin"),
        boundaryYMax: document.getElementById("boundaryYMax"),
        boundaryChannelAxis: document.getElementById("boundaryChannelAxis"),
        boundaryChannelWidth: document.getElementById("boundaryChannelWidth"),
        boundaryZGain: document.getElementById("boundaryZGain"),
        boundaryZThreshold: document.getElementById("boundaryZThreshold"),
        boundaryZPower: document.getElementById("boundaryZPower"),
        fitBoundaryWalls: document.getElementById("fitBoundaryWalls"),
        clearBoundaryWalls: document.getElementById("clearBoundaryWalls"),
        twoCellAlpha: document.getElementById("twoCellAlpha"),
        twoCellZ: document.getElementById("twoCellZ"),
        twoCellHoleMax: document.getElementById("twoCellHoleMax"),
        twoCellSweepSteps: document.getElementById("twoCellSweepSteps"),
        twoCellLeftPositionLocked: document.getElementById("twoCellLeftPositionLocked"),
        twoCellRightLocked: document.getElementById("twoCellRightLocked"),
        twoCellRightPositionLocked: document.getElementById("twoCellRightPositionLocked"),
        twoCellSuiteCase: document.getElementById("twoCellSuiteCase"),
        twoCellExternalCase: document.getElementById("twoCellExternalCase"),
        twoCellAtlasCase: document.getElementById("twoCellAtlasCase"),
        runTwoCellBench: document.getElementById("runTwoCellBench"),
        runTwoCellSuite: document.getElementById("runTwoCellSuite"),
        runTwoCellFidelityMatrix: document.getElementById("runTwoCellFidelityMatrix"),
        runTwoCellPhaseMap: document.getElementById("runTwoCellPhaseMap"),
        runTwoCellPhaseDiagram: document.getElementById("runTwoCellPhaseDiagram"),
        runTwoCellTransitionReport: document.getElementById("runTwoCellTransitionReport"),
        runTwoCellResponseAtlas: document.getElementById("runTwoCellResponseAtlas"),
        applyTwoCellBench: document.getElementById("applyTwoCellBench"),
        applyTwoCellSuiteCase: document.getElementById("applyTwoCellSuiteCase"),
        applyTwoCellExternalCase: document.getElementById("applyTwoCellExternalCase"),
        applyTwoCellAtlasCase: document.getElementById("applyTwoCellAtlasCase"),
        saveTwoCellBench: document.getElementById("saveTwoCellBench"),
        saveTwoCellSweep: document.getElementById("saveTwoCellSweep"),
        saveTwoCellSuite: document.getElementById("saveTwoCellSuite"),
        saveTwoCellSuiteCsv: document.getElementById("saveTwoCellSuiteCsv"),
        saveTwoCellFidelityMatrix: document.getElementById("saveTwoCellFidelityMatrix"),
        saveTwoCellFidelityMatrixCsv: document.getElementById("saveTwoCellFidelityMatrixCsv"),
        saveTwoCellPhaseMap: document.getElementById("saveTwoCellPhaseMap"),
        saveTwoCellPhaseMapCsv: document.getElementById("saveTwoCellPhaseMapCsv"),
        saveTwoCellPhaseDiagram: document.getElementById("saveTwoCellPhaseDiagram"),
        saveTwoCellPhaseDiagramCsv: document.getElementById("saveTwoCellPhaseDiagramCsv"),
        saveTwoCellTransitionReport: document.getElementById("saveTwoCellTransitionReport"),
        saveTwoCellTransitionReportCsv: document.getElementById("saveTwoCellTransitionReportCsv"),
        loadTwoCellTransitionMeasurements: document.getElementById("loadTwoCellTransitionMeasurements"),
        saveTwoCellTransitionComparison: document.getElementById("saveTwoCellTransitionComparison"),
        saveTwoCellTransitionRerun: document.getElementById("saveTwoCellTransitionRerun"),
        applyTwoCellTransitionCalibration: document.getElementById("applyTwoCellTransitionCalibration"),
        saveTwoCellResponseAtlas: document.getElementById("saveTwoCellResponseAtlas"),
        saveTwoCellResponseAtlasCsv: document.getElementById("saveTwoCellResponseAtlasCsv"),
        saveTwoCellFidelityMatrixTemplateCsv: document.getElementById("saveTwoCellFidelityMatrixTemplateCsv"),
        saveTwoCellExternalFidelityManifest: document.getElementById("saveTwoCellExternalFidelityManifest"),
        saveTwoCellCadContactDecomposition: document.getElementById("saveTwoCellCadContactDecomposition"),
        saveTwoCellCadContactDecompositionCsv: document.getElementById("saveTwoCellCadContactDecompositionCsv"),
        saveTwoCellExactContactHandoffPlan: document.getElementById("saveTwoCellExactContactHandoffPlan"),
        saveTwoCellExactContactHandoffPlanCsv: document.getElementById("saveTwoCellExactContactHandoffPlanCsv"),
        loadTwoCellFidelityMatrixMeasurements: document.getElementById("loadTwoCellFidelityMatrixMeasurements"),
        saveTwoCellFidelityMatrixComparison: document.getElementById("saveTwoCellFidelityMatrixComparison"),
        saveTwoCellFidelityMatrixCalibration: document.getElementById("saveTwoCellFidelityMatrixCalibration"),
        applyTwoCellFidelityMatrixCalibration: document.getElementById("applyTwoCellFidelityMatrixCalibration"),
        saveTwoCellExternalFidelitySummary: document.getElementById("saveTwoCellExternalFidelitySummary"),
        twoCellFidelityMatrixMeasurementFileInput: document.getElementById("twoCellFidelityMatrixMeasurementFileInput"),
        twoCellTransitionMeasurementFileInput: document.getElementById("twoCellTransitionMeasurementFileInput"),
        saveTwoCellConnectorTemplateCsv: document.getElementById("saveTwoCellConnectorTemplateCsv"),
        saveTwoCellPhysicalStatus: document.getElementById("saveTwoCellPhysicalStatus"),
        twoCellBenchSummary: document.getElementById("twoCellBenchSummary"),
        twoCellContactSummary: document.getElementById("twoCellContactSummary"),
        twoCellEnergySummary: document.getElementById("twoCellEnergySummary"),
        twoCellSweepSummary: document.getElementById("twoCellSweepSummary"),
        twoCellSuiteSummary: document.getElementById("twoCellSuiteSummary"),
        twoCellSuiteCaseSummary: document.getElementById("twoCellSuiteCaseSummary"),
        twoCellFidelityMatrixSummary: document.getElementById("twoCellFidelityMatrixSummary"),
        twoCellPhaseMapSummary: document.getElementById("twoCellPhaseMapSummary"),
        twoCellPhaseDiagramSummary: document.getElementById("twoCellPhaseDiagramSummary"),
        twoCellTransitionReportSummary: document.getElementById("twoCellTransitionReportSummary"),
        twoCellTransitionComparisonSummary: document.getElementById("twoCellTransitionComparisonSummary"),
        twoCellTransitionRerunSummary: document.getElementById("twoCellTransitionRerunSummary"),
        twoCellPhaseDiagramGrid: document.getElementById("twoCellPhaseDiagramGrid"),
        twoCellResponseAtlasSummary: document.getElementById("twoCellResponseAtlasSummary"),
        twoCellResponseAtlasCaseSummary: document.getElementById("twoCellResponseAtlasCaseSummary"),
        twoCellMatrixComparisonSummary: document.getElementById("twoCellMatrixComparisonSummary"),
        twoCellMatrixCalibrationSummary: document.getElementById("twoCellMatrixCalibrationSummary"),
        twoCellExternalFidelitySummary: document.getElementById("twoCellExternalFidelitySummary"),
        twoCellCorrectionSummary: document.getElementById("twoCellCorrectionSummary"),
        twoCellExternalCaseSummary: document.getElementById("twoCellExternalCaseSummary"),
        twoCellExternalCaseDetail: document.getElementById("twoCellExternalCaseDetail"),
        twoCellPhysicalStatusSummary: document.getElementById("twoCellPhysicalStatusSummary"),
        twoCellArchiveSummary: document.getElementById("twoCellArchiveSummary"),
        twoCellCadContactDecompositionSummary: document.getElementById("twoCellCadContactDecompositionSummary"),
        twoCellExactContactPlanSummary: document.getElementById("twoCellExactContactPlanSummary"),
        quickDock: document.querySelector(".quick-actuation-dock"),
        viewportTitle: document.getElementById("viewportTitle"),
        toggleQuickDock: document.getElementById("toggleQuickDock"),
        quickDockMini: document.getElementById("quickDockMini"),
        paintMode: document.getElementById("paintMode"),
        paintRadius: document.getElementById("paintRadius"),
        alphaCommand: document.getElementById("alphaCommand"),
        zCommand: document.getElementById("zCommand"),
        locked: document.getElementById("locked"),
        positionLocked: document.getElementById("positionLocked"),
        actuatorAllowed: document.getElementById("actuatorAllowed"),
        removed: document.getElementById("removed"),
        cellVisualMode: document.getElementById("cellVisualMode"),
        simulationMode: document.getElementById("simulationMode"),
        empiricalLockModel: document.getElementById("empiricalLockModel"),
        lockStateIndex: document.getElementById("lockStateIndex"),
        overlayMode: document.getElementById("overlayMode"),
        showMembrane: document.getElementById("showMembrane"),
        showReference: document.getElementById("showReference"),
        showDisplacementVectors: document.getElementById("showDisplacementVectors"),
        showExternalCase: document.getElementById("showExternalCase"),
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
        characterizationScope: document.getElementById("characterizationScope"),
        cellGrid: document.getElementById("cellGrid"),
        sequenceChart: document.getElementById("sequenceChart"),
        fileInput: document.getElementById("fileInput"),
        sequenceFileInput: document.getElementById("sequenceFileInput"),
        calibrationResultsFileInput: document.getElementById("calibrationResultsFileInput"),
        calibrationHoldoutFileInput: document.getElementById("calibrationHoldoutFileInput"),
        saveObj: document.getElementById("saveObj"),
        provenancePaper: document.getElementById("provenancePaper"),
        provenanceAssumptions: document.getElementById("provenanceAssumptions"),
        provenanceDiagnostics: document.getElementById("provenanceDiagnostics"),
        provenanceGaps: document.getElementById("provenanceGaps"),
        modelProvenanceList: document.getElementById("modelProvenanceList"),
        operatorCommitLock: document.getElementById("operatorCommitLock"),
        operatorReleaseLock: document.getElementById("operatorReleaseLock"),
        operatorRemoveCell: document.getElementById("operatorRemoveCell"),
        operatorRestoreCell: document.getElementById("operatorRestoreCell"),
        operatorCheckOrder: document.getElementById("operatorCheckOrder"),
        positionLockRowEnds: document.getElementById("positionLockRowEnds"),
        operatorOrderState: document.getElementById("operatorOrderState"),
        operatorAlphaError: document.getElementById("operatorAlphaError"),
        operatorHeightError: document.getElementById("operatorHeightError"),
        operatorSequenceState: document.getElementById("operatorSequenceState"),
        operatorMaxOrderError: document.getElementById("operatorMaxOrderError"),
        compareModelProfile: document.getElementById("compareModelProfile"),
        saveModelProfileComparison: document.getElementById("saveModelProfileComparison"),
        checkModelProfileHoldout: document.getElementById("checkModelProfileHoldout"),
        saveModelProfileHoldout: document.getElementById("saveModelProfileHoldout"),
        saveModelProfileHoldoutCsv: document.getElementById("saveModelProfileHoldoutCsv"),
        saveBenchExecutionValidation: document.getElementById("saveBenchExecutionValidation"),
        saveBenchExecutionValidationCsv: document.getElementById("saveBenchExecutionValidationCsv"),
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
      this.renderModelProvenance();
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
      this.els.paperSideLengthMm.addEventListener("input", () => {
        this.state.grid.paperSideLengthMm = Math.max(1e-9, Number(this.els.paperSideLengthMm.value));
        this.state.grid.hardwareProfile = {
          ...(this.state.grid.hardwareProfile || {}),
          sideLengthMm: this.state.grid.paperSideLengthMm,
        };
        this.onChange(this.state);
      });
      this.els.paperHoleToleranceMm.addEventListener("input", () => {
        this.state.grid.paperHoleToleranceMm = Math.max(0, Number(this.els.paperHoleToleranceMm.value));
        this.state.grid.hardwareProfile = {
          ...(this.state.grid.hardwareProfile || {}),
          fabricationHoleToleranceMm: this.state.grid.paperHoleToleranceMm,
        };
        this.onChange(this.state);
      });
      this.els.hardwareProfileName.addEventListener("input", () => this.updateHardwareProfileText("name", this.els.hardwareProfileName.value));
      this.els.hardwareProfileSource.addEventListener("input", () => this.updateHardwareProfileText("source", this.els.hardwareProfileSource.value));
      for (const [field, element] of this.hardwareProfileNumberInputs()) {
        element.addEventListener("input", () => this.updateHardwareProfileNumber(field, element.value));
      }
      this.els.applyHardwareProfile.addEventListener("click", () => this.applyHardwareProfile());
      this.els.applyA360CadProfile.addEventListener("click", () => this.applyA360CadProfile());
      this.els.clearHardwareProfile.addEventListener("click", () => this.clearHardwareProfileMeasurements());
      this.els.saveHardwareProfile.addEventListener("click", () => this.saveHardwareProfile());
      this.els.loadHardwareProfile.addEventListener("click", () => this.els.hardwareProfileFileInput.click());
      this.els.hardwareProfileFileInput.addEventListener("change", () => this.loadHardwareProfile());
      this.els.saveCalibrationPlan.addEventListener("click", () => this.saveCalibrationPlan());
      for (const element of [
        this.els.boundaryMode,
        this.els.boundaryWallType,
        this.els.boundaryXMin,
        this.els.boundaryXMax,
        this.els.boundaryYMin,
        this.els.boundaryYMax,
        this.els.boundaryChannelAxis,
        this.els.boundaryChannelWidth,
        this.els.boundaryZGain,
        this.els.boundaryZThreshold,
        this.els.boundaryZPower,
      ]) {
        element.addEventListener("input", () => this.updateBoundaryControls());
        element.addEventListener("change", () => this.updateBoundaryControls());
      }
      this.els.fitBoundaryWalls.addEventListener("click", () => {
        if (typeof RAD.fitBoundaryToReference === "function") RAD.fitBoundaryToReference(this.state);
        this.syncControls();
        this.onChange(this.state);
      });
      this.els.clearBoundaryWalls.addEventListener("click", () => {
        if (typeof RAD.clearBoundaryConstraints === "function") RAD.clearBoundaryConstraints(this.state);
        this.syncControls();
        this.onChange(this.state);
      });
      for (const element of [
        this.els.twoCellAlpha,
        this.els.twoCellZ,
        this.els.twoCellHoleMax,
        this.els.twoCellSweepSteps,
        this.els.twoCellLeftPositionLocked,
        this.els.twoCellRightLocked,
        this.els.twoCellRightPositionLocked,
      ]) {
        element.addEventListener("input", () => this.updateTwoCellBenchControls());
        element.addEventListener("change", () => this.updateTwoCellBenchControls());
      }
      this.els.runTwoCellBench.addEventListener("click", () => this.runTwoCellBench());
      this.els.runTwoCellSuite.addEventListener("click", () => this.runTwoCellSuite());
      this.els.runTwoCellFidelityMatrix.addEventListener("click", () => this.runTwoCellFidelityMatrix());
      this.els.runTwoCellPhaseMap.addEventListener("click", () => this.runTwoCellPhaseMap());
      this.els.runTwoCellPhaseDiagram.addEventListener("click", () => this.runTwoCellPhaseDiagram());
      this.els.runTwoCellTransitionReport.addEventListener("click", () => this.runTwoCellTransitionReport());
      this.els.runTwoCellResponseAtlas.addEventListener("click", () => this.runTwoCellResponseAtlas());
      this.els.applyTwoCellBench.addEventListener("click", () => this.applyTwoCellBench());
      this.els.applyTwoCellSuiteCase.addEventListener("click", () => this.applyTwoCellSuiteCase());
      this.els.applyTwoCellExternalCase.addEventListener("click", () => this.applyTwoCellExternalCase());
      this.els.applyTwoCellAtlasCase.addEventListener("click", () => this.applyTwoCellAtlasCase());
      this.els.saveTwoCellBench.addEventListener("click", () => this.saveTwoCellBench());
      this.els.saveTwoCellSweep.addEventListener("click", () => this.saveTwoCellSweep());
      this.els.saveTwoCellSuite.addEventListener("click", () => this.saveTwoCellSuite());
      this.els.saveTwoCellSuiteCsv.addEventListener("click", () => this.saveTwoCellSuiteCsv());
      this.els.saveTwoCellFidelityMatrix.addEventListener("click", () => this.saveTwoCellFidelityMatrix());
      this.els.saveTwoCellFidelityMatrixCsv.addEventListener("click", () => this.saveTwoCellFidelityMatrixCsv());
      this.els.saveTwoCellPhaseMap.addEventListener("click", () => this.saveTwoCellPhaseMap());
      this.els.saveTwoCellPhaseMapCsv.addEventListener("click", () => this.saveTwoCellPhaseMapCsv());
      this.els.saveTwoCellPhaseDiagram.addEventListener("click", () => this.saveTwoCellPhaseDiagram());
      this.els.saveTwoCellPhaseDiagramCsv.addEventListener("click", () => this.saveTwoCellPhaseDiagramCsv());
      this.els.saveTwoCellTransitionReport.addEventListener("click", () => this.saveTwoCellTransitionReport());
      this.els.saveTwoCellTransitionReportCsv.addEventListener("click", () => this.saveTwoCellTransitionReportCsv());
      this.els.loadTwoCellTransitionMeasurements.addEventListener("click", () => this.els.twoCellTransitionMeasurementFileInput.click());
      this.els.twoCellTransitionMeasurementFileInput.addEventListener("change", () => this.loadTwoCellTransitionMeasurements());
      this.els.saveTwoCellTransitionComparison.addEventListener("click", () => this.saveTwoCellTransitionComparison());
      this.els.saveTwoCellTransitionRerun.addEventListener("click", () => this.saveTwoCellTransitionRerun());
      this.els.applyTwoCellTransitionCalibration.addEventListener("click", () => this.applyTwoCellTransitionCalibration());
      this.els.saveTwoCellResponseAtlas.addEventListener("click", () => this.saveTwoCellResponseAtlas());
      this.els.saveTwoCellResponseAtlasCsv.addEventListener("click", () => this.saveTwoCellResponseAtlasCsv());
      this.els.saveTwoCellFidelityMatrixTemplateCsv.addEventListener("click", () => this.saveTwoCellFidelityMatrixTemplateCsv());
      this.els.saveTwoCellExternalFidelityManifest.addEventListener("click", () => this.saveTwoCellExternalFidelityManifest());
      this.els.saveTwoCellCadContactDecomposition.addEventListener("click", () => this.saveTwoCellCadContactDecomposition());
      this.els.saveTwoCellCadContactDecompositionCsv.addEventListener("click", () => this.saveTwoCellCadContactDecompositionCsv());
      this.els.saveTwoCellExactContactHandoffPlan.addEventListener("click", () => this.saveTwoCellExactContactHandoffPlan());
      this.els.saveTwoCellExactContactHandoffPlanCsv.addEventListener("click", () => this.saveTwoCellExactContactHandoffPlanCsv());
      this.els.loadTwoCellFidelityMatrixMeasurements.addEventListener("click", () => this.els.twoCellFidelityMatrixMeasurementFileInput.click());
      this.els.twoCellFidelityMatrixMeasurementFileInput.addEventListener("change", () => this.loadTwoCellFidelityMatrixMeasurements());
      this.els.saveTwoCellFidelityMatrixComparison.addEventListener("click", () => this.saveTwoCellFidelityMatrixComparison());
      this.els.saveTwoCellFidelityMatrixCalibration.addEventListener("click", () => this.saveTwoCellFidelityMatrixCalibration());
      this.els.applyTwoCellFidelityMatrixCalibration.addEventListener("click", () => this.applyTwoCellFidelityMatrixCalibration());
      this.els.saveTwoCellExternalFidelitySummary.addEventListener("click", () => this.saveTwoCellExternalFidelitySummary());
      this.els.saveTwoCellConnectorTemplateCsv.addEventListener("click", () => this.saveTwoCellConnectorTemplateCsv());
      this.els.saveTwoCellPhysicalStatus.addEventListener("click", () => this.saveTwoCellPhysicalStatus());
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
      this.els.positionLocked.addEventListener("change", () => this.previewSelected());
      this.els.removed.addEventListener("change", () => this.previewSelected());
      this.els.actuatorAllowed.addEventListener("change", () => this.updateSelectedActuatorMask());
      this.els.cellVisualMode.addEventListener("change", () => {
        this.applyCellVisualMode(this.els.cellVisualMode.value);
        this.onChange(this.state);
      });
      this.els.simulationMode.addEventListener("change", () => {
        this.state.view.simulationMode = ["springPreview", "constraintSolved"].includes(this.els.simulationMode.value) ? this.els.simulationMode.value : "kinematic";
        this.onChange(this.state);
      });
      this.els.empiricalLockModel.addEventListener("change", () => {
        this.state.grid.empiricalLockModel = this.els.empiricalLockModel.checked;
        this.onChange(this.state);
      });
      this.els.lockStateIndex.addEventListener("input", () => {
        this.state.grid.lockStateIndex = Math.max(1, Math.min(3, Math.round(Number(this.els.lockStateIndex.value) || 1)));
        this.onChange(this.state);
      });
      this.els.overlayMode.addEventListener("change", () => {
        this.state.view.overlayMode = this.els.overlayMode.value;
        if (this.state.view.overlayMode === "modelError") {
          this.state.view.simulationMode = "springPreview";
          this.els.simulationMode.value = "springPreview";
        }
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
      this.els.showExternalCase.addEventListener("change", () => {
        this.state.view.externalCaseVisible = this.els.showExternalCase.checked;
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
      this.els.operatorRemoveCell.addEventListener("click", () => this.removeSelectedCellEvent());
      this.els.operatorRestoreCell.addEventListener("click", () => this.restoreSelectedCellEvent());
      this.els.operatorCheckOrder.addEventListener("click", () => this.checkSelectedEventOrder());
      this.els.positionLockRowEnds.addEventListener("click", () => this.positionLockSelectedRowEnds());
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
          underactuatedTargets: jacobian.targetReachability?.underactuatedHeightCells || 0,
        });
        this.state.view.overlayMode = jacobian.targetReachability?.underactuatedHeightCells ? "underactuated" : "reachability";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("selectUnderactuatedTarget").addEventListener("click", () => this.selectUnderactuatedTarget());
      document.getElementById("solveLinearFit").addEventListener("click", () => {
        const solution = RAD.solveLinearizedTargetFit(this.state);
        RAD.recordEvent(this.state, {
          type: "linear-fit-solved",
          actuators: solution.commands.length,
          steps: solution.steps,
          baseError: solution.baseError,
          projectedError: solution.projectedError,
          underactuatedTargets: solution.underactuatedHeightCells || 0,
        });
        this.state.view.overlayMode = "inverse";
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("validatePhysicalFit").addEventListener("click", () => {
        const validation = RAD.validateInversePlanPhysical(this.state);
        RAD.recordEvent(this.state, {
          type: "inverse-physical-validated",
          source: validation.source,
          actuators: validation.commandCount,
          physicalError: validation.physicalProjectedError,
          modelError: validation.centerModelRms,
        });
        this.state.view.simulationMode = "springPreview";
        this.els.simulationMode.value = "springPreview";
        this.state.view.overlayMode = "modelError";
        this.state.view.targetVisible = this.state.target.type !== "none";
        this.syncControls();
        this.onChange(this.state);
      });
      document.getElementById("saveInverseReport").addEventListener("click", () => this.saveInverseReport());
      document.getElementById("buildPreviewPacket").addEventListener("click", () => this.buildProfileInversePreviewPacket());
      document.getElementById("runPreviewReplay").addEventListener("click", () => this.runProfileInversePreviewReplay());
      document.getElementById("runPreviewPhysical").addEventListener("click", () => this.runProfileInversePreviewPhysical());
      document.getElementById("savePreviewPacket").addEventListener("click", () => this.saveProfileInversePreviewPacket());
      document.getElementById("savePreviewReplay").addEventListener("click", () => this.saveProfileInversePreviewReplay());
      document.getElementById("savePreviewPhysical").addEventListener("click", () => this.saveProfileInversePreviewPhysical());
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
        this.state.inverse.jacobian = { columns: [], coverageMap: RAD.matrix(this.state.grid.rows, this.state.grid.cols, 0), actuatorCount: 0, columnCount: 0, meanCoverage: 0, maxCoverage: 0, stepZ: 0.12, stepAlpha: 0.12, conditionEstimate: 0, targetReachability: null };
        this.state.inverse.linearSolution = { commands: [], history: [], steps: 0, baseError: 0, predictedError: 0, projectedError: 0, projectedActuators: 0, targetReachability: null };
        this.state.inverse.physicalValidation = null;
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
      this.els.characterizationScope.addEventListener("change", () => {
        this.state.experiment.characterizationScope = this.els.characterizationScope.value;
        this.updateLabels(null);
      });
      document.getElementById("runCharacterization").addEventListener("click", () => this.runCharacterization());
      document.getElementById("saveResponseMatrix").addEventListener("click", () => this.saveResponseMatrix());
      document.getElementById("saveTopologyReport").addEventListener("click", () => this.saveTopologyReport());
      document.getElementById("saveProgrammableReport").addEventListener("click", () => this.saveProgrammableReport());
      document.getElementById("saveFormalizationTargets").addEventListener("click", () => this.saveFormalizationTargets());
      document.getElementById("selectInteractionHotspot").addEventListener("click", () => this.selectInteractionHotspot());
      document.getElementById("selectCalibrationHotspot").addEventListener("click", () => this.selectCalibrationHotspot());
      document.getElementById("nextCalibrationHotspot").addEventListener("click", () => this.selectCalibrationHotspot(1));
      document.getElementById("saveExperimentProtocol").addEventListener("click", () => this.saveExperimentProtocol());
      document.getElementById("saveResponseAtlas").addEventListener("click", () => this.saveResponseAtlas());
      document.getElementById("runResponseAtlasSweep").addEventListener("click", () => this.runResponseAtlasSweep());
      document.getElementById("saveResponseAtlasSweep").addEventListener("click", () => this.saveResponseAtlasSweep());
      document.getElementById("saveResultsTemplate").addEventListener("click", () => this.saveResultsTemplate());
      document.getElementById("saveBenchPacket").addEventListener("click", () => this.saveBenchPacket());
      document.getElementById("saveBenchNotebook").addEventListener("click", () => this.saveBenchNotebook());
      document.getElementById("saveBenchNotebookCsv").addEventListener("click", () => this.saveBenchNotebookCsv());
      document.getElementById("loadResultsJson").addEventListener("click", () => this.els.calibrationResultsFileInput.click());
      this.els.calibrationResultsFileInput.addEventListener("change", () => this.loadCalibrationResults());
      document.getElementById("loadHoldoutJson").addEventListener("click", () => this.els.calibrationHoldoutFileInput.click());
      this.els.calibrationHoldoutFileInput.addEventListener("change", () => this.loadCalibrationHoldoutResults());
      document.getElementById("saveComparisonReport").addEventListener("click", () => this.saveComparisonReport());
      document.getElementById("saveModelProfile").addEventListener("click", () => this.saveModelProfile());
      this.els.compareModelProfile.addEventListener("click", () => this.compareModelProfile());
      this.els.saveModelProfileComparison.addEventListener("click", () => this.saveModelProfileComparison());
      this.els.checkModelProfileHoldout.addEventListener("click", () => this.checkModelProfileHoldout());
      this.els.saveModelProfileHoldout.addEventListener("click", () => this.saveModelProfileHoldout());
      this.els.saveModelProfileHoldoutCsv.addEventListener("click", () => this.saveModelProfileHoldoutCsv());
      this.els.saveBenchExecutionValidation.addEventListener("click", () => this.saveBenchExecutionValidation());
      this.els.saveBenchExecutionValidationCsv.addEventListener("click", () => this.saveBenchExecutionValidationCsv());
      document.getElementById("applyModelProfile").addEventListener("click", () => this.applyModelProfile());
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
      const removed = this.els.removed.checked;
      for (const cell of cells) {
        this.state.cells.removed[cell.r][cell.c] = removed;
        if (removed) {
          this.state.cells.commandAlpha[cell.r][cell.c] = 0;
          this.state.cells.commandZ[cell.r][cell.c] = 0;
          this.state.cells.locked[cell.r][cell.c] = false;
          this.state.cells.positionLocked[cell.r][cell.c] = false;
          if (this.state.cells.lockZ) this.state.cells.lockZ[cell.r][cell.c] = 0;
          this.clearPositionLockTarget(cell.r, cell.c);
          this.state.cells.actuatorAllowed[cell.r][cell.c] = false;
        } else {
          this.state.cells.commandAlpha[cell.r][cell.c] = alpha;
          this.state.cells.commandZ[cell.r][cell.c] = z;
          if (this.els.locked.checked) this.commitDirectLockState(cell.r, cell.c);
          if (this.els.positionLocked.checked && !this.state.cells.positionLocked?.[cell.r]?.[cell.c]) {
            this.commitPositionLockState(cell.r, cell.c);
          }
          this.state.cells.locked[cell.r][cell.c] = this.els.locked.checked;
          this.state.cells.positionLocked[cell.r][cell.c] = this.els.positionLocked.checked;
          this.state.cells.actuatorAllowed[cell.r][cell.c] = this.els.actuatorAllowed.checked;
        }
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
        positionLocked: this.els.positionLocked.checked,
        actuatorAllowed: this.els.actuatorAllowed.checked,
        removed,
        operator: removed ? null : RAD.groupActuationEvent(cells, alpha, z),
      });
      this.syncControls();
      this.onChange(this.state);
    }

    paintBrushCells(r, c) {
      return RAD.brushCells(this.state, r, c, this.state.view.paintRadius);
    }

    commitDirectLockState(r, c) {
      if (!this.state.cells.lockAlpha) {
        this.state.cells.lockAlpha = RAD.matrix(this.state.grid.rows, this.state.grid.cols, this.state.grid.initialAlpha);
      }
      if (!this.state.cells.lockZ) {
        this.state.cells.lockZ = RAD.matrix(this.state.grid.rows, this.state.grid.cols, 0);
      }
      const preview = JSON.parse(JSON.stringify(this.state));
      preview.cells.locked[r][c] = false;
      const sim = RAD.simulate(preview);
      this.state.cells.lockAlpha[r][c] = sim.alpha[r][c];
      this.state.cells.lockZ[r][c] = sim.height[r][c];
    }

    ensurePositionLockTargetMatrices() {
      const { rows, cols } = this.state.grid;
      if (!this.state.cells.positionLockX) this.state.cells.positionLockX = RAD.matrix(rows, cols, null);
      if (!this.state.cells.positionLockY) this.state.cells.positionLockY = RAD.matrix(rows, cols, null);
      if (!this.state.cells.positionLockZ) this.state.cells.positionLockZ = RAD.matrix(rows, cols, null);
    }

    clearPositionLockTarget(r, c) {
      this.ensurePositionLockTargetMatrices();
      this.state.cells.positionLockX[r][c] = null;
      this.state.cells.positionLockY[r][c] = null;
      this.state.cells.positionLockZ[r][c] = null;
    }

    commitPositionLockState(r, c) {
      this.ensurePositionLockTargetMatrices();
      const preview = JSON.parse(JSON.stringify(this.state));
      preview.cells.positionLocked[r][c] = false;
      const sim = typeof RAD.simulateActive === "function" ? RAD.simulateActive(preview) : RAD.simulate(preview);
      const center = sim.centers?.[r]?.[c] || (typeof RAD.referenceCenter === "function" ? RAD.referenceCenter(this.state, r, c) : { x: 0, y: 0, z: 0 });
      this.state.cells.positionLockX[r][c] = Number(center.x) || 0;
      this.state.cells.positionLockY[r][c] = Number(center.y) || 0;
      this.state.cells.positionLockZ[r][c] = Number(center.z) || 0;
    }

    currentTwoCellBenchControls(overrides = {}) {
      const uiValues = {
        alphaCommand: Number(this.els.twoCellAlpha.value),
        zCommand: Number(this.els.twoCellZ.value),
        holeSweepMax: Number(this.els.twoCellHoleMax.value),
        holeSweepSteps: Number(this.els.twoCellSweepSteps.value),
        leftPositionLocked: this.els.twoCellLeftPositionLocked.checked,
        rightLocked: this.els.twoCellRightLocked.checked,
        rightPositionLocked: this.els.twoCellRightPositionLocked.checked,
        ...overrides,
      };
      return typeof RAD.twoCellBenchControls === "function"
        ? RAD.twoCellBenchControls(this.state, uiValues)
        : { ...(this.state.experiment?.twoCellBenchControls || {}), ...uiValues };
    }

    updateTwoCellBenchControls() {
      const controls = this.currentTwoCellBenchControls();
      this.state.experiment.twoCellBenchControls = controls;
      this.renderTwoCellBench();
      return controls;
    }

    runTwoCellBench() {
      if (typeof RAD.simulateTwoCellBench !== "function" || typeof RAD.sweepTwoCellBacklash !== "function") return null;
      const controls = this.updateTwoCellBenchControls();
      const bench = RAD.simulateTwoCellBench(this.state, controls);
      const sweep = RAD.sweepTwoCellBacklash(this.state, controls);
      const contact = typeof RAD.twoCellConnectorContactReport === "function" ? RAD.twoCellConnectorContactReport(this.state, controls) : null;
      const suite = typeof RAD.twoCellPhysicalSimulationSuite === "function" ? RAD.twoCellPhysicalSimulationSuite(this.state, controls) : null;
      this.state.experiment.twoCellBench = { bench, sweep, contact, suite, ranAt: new Date().toISOString() };
      if (suite) this.state.experiment.twoCellPhysicalSuite = suite;
      RAD.recordEvent(this.state, {
        type: "two-cell-bench-run",
        alphaCommand: controls.alphaCommand,
        zCommand: controls.zCommand,
        contactMode: bench.connector.contactMode,
        totalEnergy: bench.energy.totalEnergy,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return this.state.experiment.twoCellBench;
    }

    runTwoCellSuite() {
      if (typeof RAD.twoCellPhysicalSimulationSuite !== "function") return null;
      const controls = this.updateTwoCellBenchControls();
      const suite = RAD.twoCellPhysicalSimulationSuite(this.state, controls);
      this.state.experiment.twoCellPhysicalSuite = suite;
      this.state.experiment.twoCellBench = {
        ...(this.state.experiment.twoCellBench || {}),
        suite,
        ranAt: new Date().toISOString(),
      };
      RAD.recordEvent(this.state, {
        type: "two-cell-physical-suite-run",
        cases: suite.summary?.caseCount || 0,
        ready: suite.summary?.internalSuiteReady === true,
        missing: suite.summary?.missingEvidence?.length || 0,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return suite;
    }

    runTwoCellFidelityMatrix() {
      if (typeof RAD.twoCellPhysicalFidelityMatrix !== "function") return null;
      const controls = this.updateTwoCellBenchControls();
      const matrix = RAD.twoCellPhysicalFidelityMatrix(this.state, controls);
      this.state.experiment.twoCellPhysicalFidelityMatrix = matrix;
      this.state.experiment.twoCellBench = {
        ...(this.state.experiment.twoCellBench || {}),
        fidelityMatrix: matrix,
        ranAt: new Date().toISOString(),
      };
      RAD.recordEvent(this.state, {
        type: "two-cell-fidelity-matrix-run",
        rows: matrix.summary?.rowCount || 0,
        connectorRows: matrix.summary?.connectorRowCount || 0,
        ready: matrix.summary?.internalMatrixReady === true,
        missing: matrix.summary?.missingEvidence?.length || 0,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return matrix;
    }

    runTwoCellPhaseMap() {
      if (typeof RAD.twoCellContactPhaseMap !== "function") return null;
      const controls = this.updateTwoCellBenchControls();
      const phaseMap = RAD.twoCellContactPhaseMap(this.state, controls);
      this.state.experiment.twoCellContactPhaseMap = phaseMap;
      this.state.experiment.twoCellBench = {
        ...(this.state.experiment.twoCellBench || {}),
        contactPhaseMap: phaseMap,
        fidelityMatrix: this.state.experiment.twoCellPhysicalFidelityMatrix || this.state.experiment.twoCellBench?.fidelityMatrix,
        ranAt: new Date().toISOString(),
      };
      RAD.recordEvent(this.state, {
        type: "two-cell-contact-phase-map-run",
        rows: phaseMap.summary?.rowCount || 0,
        dominantPhase: phaseMap.summary?.dominantPhase || "none",
        activeContactCases: phaseMap.summary?.activeContactCaseCount || 0,
        priorityCases: phaseMap.summary?.measurementPriority?.length || 0,
        physicalAccuracyValidated: phaseMap.summary?.physicalAccuracyValidated === true,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return phaseMap;
    }

    runTwoCellPhaseDiagram() {
      if (typeof RAD.twoCellRadiusBacklashPhaseDiagram !== "function") return null;
      const controls = this.updateTwoCellBenchControls();
      const diagram = RAD.twoCellRadiusBacklashPhaseDiagram(this.state, controls);
      this.state.experiment.twoCellRadiusBacklashPhaseDiagram = diagram;
      this.state.experiment.twoCellBench = {
        ...(this.state.experiment.twoCellBench || {}),
        radiusBacklashPhaseDiagram: diagram,
        ranAt: new Date().toISOString(),
      };
      RAD.recordEvent(this.state, {
        type: "two-cell-radius-backlash-phase-diagram-run",
        rows: diagram.summary?.rowCount || 0,
        holeSteps: diagram.summary?.holeStepCount || 0,
        backlashSteps: diagram.summary?.backlashStepCount || 0,
        dominantPhase: diagram.summary?.dominantPhase || "none",
        transitions: diagram.summary?.transitionCount || 0,
        physicalAccuracyValidated: diagram.summary?.physicalAccuracyValidated === true,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return diagram;
    }

    runTwoCellTransitionReport() {
      if (typeof RAD.twoCellRadiusBacklashTransitionReport !== "function") return null;
      const controls = this.updateTwoCellBenchControls();
      const report = RAD.twoCellRadiusBacklashTransitionReport(this.state, controls);
      this.state.experiment.twoCellRadiusBacklashTransitionReport = report;
      this.state.experiment.twoCellBench = {
        ...(this.state.experiment.twoCellBench || {}),
        radiusBacklashTransitionReport: report,
        ranAt: new Date().toISOString(),
      };
      RAD.recordEvent(this.state, {
        type: "two-cell-radius-backlash-transition-report-run",
        brackets: report.summary?.transitionBracketCount || 0,
        measurementCases: report.summary?.measurementCaseCount || 0,
        physicalAccuracyValidated: report.summary?.physicalAccuracyValidated === true,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return report;
    }

    runTwoCellResponseAtlas() {
      if (typeof RAD.twoCellPhysicalResponseAtlas !== "function") return null;
      const controls = this.updateTwoCellBenchControls();
      const atlas = RAD.twoCellPhysicalResponseAtlas(this.state, controls);
      this.state.experiment.twoCellPhysicalResponseAtlas = atlas;
      this.state.experiment.twoCellBench = {
        ...(this.state.experiment.twoCellBench || {}),
        responseAtlas: atlas,
        fidelityMatrix: this.state.experiment.twoCellPhysicalFidelityMatrix || this.state.experiment.twoCellBench?.fidelityMatrix,
        ranAt: new Date().toISOString(),
      };
      RAD.recordEvent(this.state, {
        type: "two-cell-response-atlas-run",
        rows: atlas.matrixSummary?.rowCount || 0,
        connectorRows: atlas.connectorSummary?.connectorRowCount || 0,
        physicalAccuracyValidated: atlas.claimBoundary?.physicalAccuracyValidated === true,
        priorityCases: atlas.benchPriority?.firstCasesToMeasure?.length || 0,
      });
      this.syncTwoCellAtlasCases(atlas);
      this.renderTwoCellBench();
      this.onChange(this.state);
      return atlas;
    }

    applyTwoCellBench() {
      const controls = this.updateTwoCellBenchControls();
      this.state = RAD.resizeState(this.state, 1, 2);
      RAD.clearCommands(this.state);
      this.state.selection = { r: 0, c: 1 };
      this.state.experiment.twoCellBenchControls = controls;
      this.state.cells.positionLocked = RAD.matrix(1, 2, false);
      this.state.cells.positionLockX = RAD.matrix(1, 2, null);
      this.state.cells.positionLockY = RAD.matrix(1, 2, null);
      this.state.cells.positionLockZ = RAD.matrix(1, 2, null);
      this.state.cells.locked = RAD.matrix(1, 2, false);
      this.state.cells.lockAlpha = RAD.matrix(1, 2, this.state.grid.initialAlpha);
      this.state.cells.lockZ = RAD.matrix(1, 2, 0);
      this.state.cells.removed = RAD.matrix(1, 2, false);
      this.state.cells.actuatorAllowed = RAD.matrix(1, 2, true);
      this.state.cells.commandAlpha[0][1] = controls.alphaCommand;
      this.state.cells.commandZ[0][1] = controls.zCommand;
      this.state.cells.positionLocked[0][0] = controls.leftPositionLocked;
      this.state.cells.locked[0][1] = controls.rightLocked;
      this.state.cells.positionLocked[0][1] = controls.rightPositionLocked;
      if (controls.leftPositionLocked) this.commitPositionLockState(0, 0);
      if (controls.rightPositionLocked) this.commitPositionLockState(0, 1);
      this.applyCellVisualMode("cadRad");
      Object.assign(this.state.view, {
        simulationMode: "springPreview",
        overlayMode: "state",
        membraneVisible: false,
        referenceVisible: true,
        displacementVectorsVisible: true,
        gapsVisible: true,
        stopsVisible: true,
        pivotsVisible: true,
        linkagesVisible: true,
        fastenersVisible: true,
        actuatorsVisible: true,
        measurementMode: "hardware",
      });
      if (typeof RAD.simulateTwoCellBench === "function" && typeof RAD.sweepTwoCellBacklash === "function") {
        const bench = RAD.simulateTwoCellBench(this.state, controls);
        const sweep = RAD.sweepTwoCellBacklash(this.state, controls);
        const contact = typeof RAD.twoCellConnectorContactReport === "function" ? RAD.twoCellConnectorContactReport(this.state, controls) : null;
        this.state.experiment.twoCellBench = { bench, sweep, contact, ranAt: new Date().toISOString(), appliedToLattice: true };
      }
      RAD.recordEvent(this.state, {
        type: "two-cell-bench-applied",
        alphaCommand: controls.alphaCommand,
        zCommand: controls.zCommand,
        leftPositionLocked: controls.leftPositionLocked,
        rightLocked: controls.rightLocked,
        rightPositionLocked: controls.rightPositionLocked,
      });
      this.syncControls();
      this.onChange(this.state);
      return this.state;
    }

    applyTwoCellSuiteCase() {
      if (typeof RAD.twoCellPhysicalSuiteCaseOptions !== "function") return this.applyTwoCellBench();
      const caseId = this.els.twoCellSuiteCase?.value || "free_contract_lift_nominal";
      const controls = this.updateTwoCellBenchControls();
      const resolved = RAD.twoCellPhysicalSuiteCaseOptions(this.state, caseId, controls);
      this.state.grid.backlash = Number(resolved.backlash);
      this.state.grid.holeRadius = Math.max(Number(this.state.grid.pinRadius || 0), Number(resolved.holeRadius));
      this.state.experiment.twoCellSuiteCaseId = resolved.caseId;
      this.state.experiment.twoCellBenchControls = resolved.controls;
      this.els.twoCellAlpha.value = resolved.controls.alphaCommand;
      this.els.twoCellZ.value = resolved.controls.zCommand;
      this.els.twoCellLeftPositionLocked.checked = resolved.controls.leftPositionLocked !== false;
      this.els.twoCellRightLocked.checked = resolved.controls.rightLocked === true;
      this.els.twoCellRightPositionLocked.checked = resolved.controls.rightPositionLocked === true;
      RAD.recordEvent(this.state, {
        type: "two-cell-suite-case-applied",
        caseId: resolved.caseId,
        backlash: resolved.backlash,
        holeRadius: resolved.holeRadius,
      });
      return this.applyTwoCellBench();
    }

    applyTwoCellExternalCase() {
      if (typeof RAD.twoCellExternalFidelityWebSummary !== "function" || typeof RAD.twoCellExternalFidelityCaseControls !== "function") {
        window.alert("External fidelity summary is not loaded.");
        return null;
      }
      const externalSummary = RAD.twoCellExternalFidelityWebSummary();
      const caseId = this.els.twoCellExternalCase?.value || this.state.experiment?.twoCellExternalCaseId;
      const preview = RAD.twoCellExternalFidelityCaseControls(caseId, externalSummary, this.state.experiment?.twoCellBenchControls || {});
      if (!preview) {
        window.alert("Select an external MuJoCo proxy case first.");
        return null;
      }
      this.state.grid.backlash = Number(preview.grid.backlash);
      this.state.grid.pinRadius = Math.max(0, Number(preview.grid.pinRadius));
      this.state.grid.holeRadius = Math.max(this.state.grid.pinRadius, Number(preview.grid.holeRadius));
      this.els.twoCellAlpha.value = preview.controls.alphaCommand;
      this.els.twoCellZ.value = preview.controls.zCommand;
      this.els.twoCellHoleMax.value = Math.max(preview.grid.holeRadius, Number(this.els.twoCellHoleMax.value || 0.5));
      this.els.twoCellLeftPositionLocked.checked = preview.controls.leftPositionLocked !== false;
      this.els.twoCellRightLocked.checked = preview.controls.rightLocked === true;
      this.els.twoCellRightPositionLocked.checked = preview.controls.rightPositionLocked === true;
      this.state.experiment.twoCellExternalCaseId = preview.caseId;
      this.state.experiment.twoCellExternalCasePreview = preview;
      this.applyTwoCellBench();
      this.state.experiment.twoCellExternalCaseId = preview.caseId;
      this.state.experiment.twoCellExternalCasePreview = preview;
      this.state.view.externalCaseVisible = true;
      RAD.recordEvent(this.state, {
        type: "two-cell-external-case-applied",
        caseId: preview.caseId,
        backlash: preview.grid.backlash,
        pinRadius: preview.grid.pinRadius,
        holeRadius: preview.grid.holeRadius,
        alphaCommand: preview.controls.alphaCommand,
        zCommand: preview.controls.zCommand,
        physicalAccuracyValidated: false,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return preview;
    }

    twoCellAtlasCaseOptions(atlas) {
      const source = atlas || this.state.experiment?.twoCellPhysicalResponseAtlas || this.state.experiment?.twoCellBench?.responseAtlas;
      const entries = [];
      const seen = new Set();
      const add = (bucket, label, cases) => {
        for (const item of cases || []) {
          const caseId = item.caseId || `${bucket}:${item.index}`;
          if (seen.has(caseId)) continue;
          seen.add(caseId);
          entries.push({ ...item, caseId, bucket, label: `${label} ${String(item.index).padStart(3, "0")} ${item.actuationCase} ${item.lockMode}` });
        }
      };
      add("firstCasesToMeasure", "Measure", source?.benchPriority?.firstCasesToMeasure);
      add("lockCasesToMeasure", "Lock", source?.benchPriority?.lockCasesToMeasure);
      add("clearanceCasesToMeasure", "Slip", source?.benchPriority?.clearanceCasesToMeasure);
      add("maxContactPenalty", "Contact", source?.rankedCases?.maxContactPenalty);
      add("maxRightHeight", "Height", source?.rankedCases?.maxRightHeight);
      return entries;
    }

    syncTwoCellAtlasCases(atlas) {
      if (!this.els.twoCellAtlasCase) return null;
      const options = this.twoCellAtlasCaseOptions(atlas);
      const previous = this.state.experiment?.twoCellResponseAtlasCaseId || this.els.twoCellAtlasCase.value || options[0]?.caseId || "";
      if (this.els.twoCellAtlasCase.options.length !== options.length || options.some((item, index) => this.els.twoCellAtlasCase.options[index]?.value !== item.caseId)) {
        this.els.twoCellAtlasCase.replaceChildren(
          ...options.map((item) => {
            const option = document.createElement("option");
            option.value = item.caseId;
            option.textContent = item.label;
            return option;
          })
        );
      }
      const selected = options.some((item) => item.caseId === previous) ? previous : options[0]?.caseId || "";
      this.els.twoCellAtlasCase.value = selected;
      return options.find((item) => item.caseId === selected) || null;
    }

    applyTwoCellAtlasCase() {
      let atlas = this.state.experiment?.twoCellPhysicalResponseAtlas || this.state.experiment?.twoCellBench?.responseAtlas;
      if (!atlas && typeof RAD.twoCellPhysicalResponseAtlas === "function") {
        atlas = this.runTwoCellResponseAtlas();
      }
      const options = this.twoCellAtlasCaseOptions(atlas);
      const caseId = this.els.twoCellAtlasCase?.value || options[0]?.caseId;
      const selected = options.find((item) => item.caseId === caseId);
      if (!selected) {
        window.alert("Run the response atlas first.");
        return null;
      }
      this.state.grid.backlash = Number(selected.backlash);
      this.state.grid.pinRadius = Math.max(0, Number(selected.pinRadius || this.state.grid.pinRadius || 0));
      this.state.grid.holeRadius = Math.max(this.state.grid.pinRadius, Number(selected.holeRadius || this.state.grid.holeRadius || 0));
      this.els.twoCellAlpha.value = selected.alphaCommand;
      this.els.twoCellZ.value = selected.zCommand;
      this.els.twoCellHoleMax.value = Math.max(Number(this.els.twoCellHoleMax.value || 0.5), Number(selected.holeRadius || 0));
      this.els.twoCellLeftPositionLocked.checked = selected.lockMode !== "left_free";
      this.els.twoCellRightLocked.checked = selected.lockMode === "right_state_locked";
      this.els.twoCellRightPositionLocked.checked = selected.lockMode === "right_position_locked";
      this.state.experiment.twoCellResponseAtlasCaseId = selected.caseId;
      this.state.experiment.twoCellResponseAtlasCasePreview = selected;
      this.applyTwoCellBench();
      this.state.experiment.twoCellResponseAtlasCaseId = selected.caseId;
      this.state.experiment.twoCellResponseAtlasCasePreview = selected;
      RAD.recordEvent(this.state, {
        type: "two-cell-response-atlas-case-applied",
        caseId: selected.caseId,
        lockMode: selected.lockMode,
        backlash: selected.backlash,
        holeRadius: selected.holeRadius,
        alphaCommand: selected.alphaCommand,
        zCommand: selected.zCommand,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return selected;
    }

    syncTwoCellSuiteCases(controls) {
      if (!this.els.twoCellSuiteCase || typeof RAD.twoCellPhysicalSuiteCaseSpecs !== "function") return;
      const previous = this.state.experiment?.twoCellSuiteCaseId || this.els.twoCellSuiteCase.value;
      const specs = RAD.twoCellPhysicalSuiteCaseSpecs(this.state, controls);
      if (this.els.twoCellSuiteCase.options.length !== specs.length) {
        this.els.twoCellSuiteCase.replaceChildren(
          ...specs.map((spec) => {
            const option = document.createElement("option");
            option.value = spec.caseId;
            option.textContent = spec.caseId.replace(/_/g, " ");
            return option;
          })
        );
      }
      const selected = specs.some((spec) => spec.caseId === previous) ? previous : specs[0]?.caseId;
      if (selected) this.els.twoCellSuiteCase.value = selected;
      const spec = specs.find((item) => item.caseId === this.els.twoCellSuiteCase.value);
      if (this.els.twoCellSuiteCaseSummary) {
        this.els.twoCellSuiteCaseSummary.textContent = spec
          ? `case ${spec.lockMode}, hole ${Number(spec.holeRadius).toFixed(3)}, b ${Number(spec.backlash).toFixed(2)}`
          : "case --";
      }
    }

    syncTwoCellExternalCases(externalSummary) {
      if (!this.els.twoCellExternalCase || typeof RAD.twoCellExternalFidelityCaseOptions !== "function") return null;
      const options = RAD.twoCellExternalFidelityCaseOptions(externalSummary);
      const previous = this.state.experiment?.twoCellExternalCaseId || this.els.twoCellExternalCase.value || options[0]?.caseId || "";
      if (this.els.twoCellExternalCase.options.length !== options.length) {
        this.els.twoCellExternalCase.replaceChildren(
          ...options.map((item) => {
            const option = document.createElement("option");
            option.value = item.caseId;
            option.textContent = item.label;
            return option;
          })
        );
      }
      const selected = options.some((item) => item.caseId === previous) ? previous : options[0]?.caseId || "";
      this.els.twoCellExternalCase.value = selected;
      return options.find((item) => item.caseId === selected) || null;
    }

    phaseDiagramClassName(phase) {
      return String(phase || "unknown").replace(/[^a-z0-9]+/gi, "-").replace(/^-|-$/g, "").toLowerCase() || "unknown";
    }

    phaseDiagramLabel(phase) {
      const labels = {
        "axial-and-vertical-contact": "AV",
        "axial-contact": "A",
        "vertical-contact": "V",
        "position-locked": "P",
        "state-locked-vertical-free": "SZ",
        "state-locked": "S",
        "gravity-sag": "G",
        "free-play": "F",
      };
      return labels[String(phase)] || "?";
    }

    renderTwoCellPhaseDiagramGrid(diagram) {
      const container = this.els.twoCellPhaseDiagramGrid;
      if (!container) return;
      container.replaceChildren();
      const holes = diagram?.axes?.holeRadii || [];
      const backlashes = diagram?.axes?.backlashValues || [];
      const rows = diagram?.rows || [];
      const grid = diagram?.phaseGrid || [];
      if (!holes.length || !backlashes.length || !rows.length || !grid.length) {
        container.hidden = true;
        return;
      }
      container.hidden = false;
      container.style.gridTemplateColumns = `repeat(${holes.length}, minmax(18px, 1fr))`;
      const selected = this.state.experiment?.twoCellPhaseDiagramCasePreview;
      for (let r = 0; r < backlashes.length; r += 1) {
        for (let c = 0; c < holes.length; c += 1) {
          const phase = grid[r]?.[c] || "unknown";
          const row = rows.find(
            (item) =>
              Math.abs(Number(item.backlash || 0) - Number(backlashes[r] || 0)) <= 1e-9 &&
              Math.abs(Number(item.holeRadius || 0) - Number(holes[c] || 0)) <= 1e-9
          );
          const button = document.createElement("button");
          button.type = "button";
          button.className = `phase-diagram-cell phase-${this.phaseDiagramClassName(phase)}`;
          if (selected && row && String(selected.caseId) === String(row.caseId)) button.classList.add("is-active");
          button.textContent = this.phaseDiagramLabel(phase);
          button.title = row
            ? `b ${Number(row.backlash || 0).toFixed(3)}, hole ${Number(row.holeRadius || 0).toFixed(3)}, ${phase}, right z ${Number(row.rightZ || 0).toFixed(3)}`
            : `${phase}`;
          button.setAttribute("aria-label", button.title);
          button.addEventListener("click", () => this.applyTwoCellPhaseDiagramGridCell(r, c));
          container.append(button);
        }
      }
    }

    applyTwoCellPhaseDiagramGridCell(rowIndex, columnIndex) {
      const diagram = this.state.experiment?.twoCellRadiusBacklashPhaseDiagram || this.state.experiment?.twoCellBench?.radiusBacklashPhaseDiagram;
      if (!diagram) return null;
      const backlash = Number(diagram.axes?.backlashValues?.[rowIndex]);
      const holeRadius = Number(diagram.axes?.holeRadii?.[columnIndex]);
      const selected = (diagram.rows || []).find(
        (row) =>
          Math.abs(Number(row.backlash || 0) - backlash) <= 1e-9 &&
          Math.abs(Number(row.holeRadius || 0) - holeRadius) <= 1e-9
      );
      if (!selected) return null;
      this.state.grid.backlash = Math.max(0, backlash);
      this.state.grid.pinRadius = Math.max(0, Number(selected.pinRadius || this.state.grid.pinRadius || 0));
      this.state.grid.holeRadius = Math.max(this.state.grid.pinRadius, holeRadius);
      this.els.twoCellAlpha.value = Number(selected.alphaCommand || 0);
      this.els.twoCellZ.value = Number(selected.zCommand || 0);
      this.els.twoCellHoleMax.value = Math.max(Number(this.els.twoCellHoleMax.value || 0), holeRadius);
      this.els.twoCellLeftPositionLocked.checked = selected.lockMode !== "left_free";
      this.els.twoCellRightLocked.checked = selected.lockMode === "right_state_locked";
      this.els.twoCellRightPositionLocked.checked = selected.lockMode === "right_position_locked";
      this.state.experiment.twoCellPhaseDiagramCasePreview = selected;
      this.applyTwoCellBench();
      this.state.experiment.twoCellRadiusBacklashPhaseDiagram = diagram;
      this.state.experiment.twoCellPhaseDiagramCasePreview = selected;
      this.state.experiment.twoCellBench = {
        ...(this.state.experiment.twoCellBench || {}),
        radiusBacklashPhaseDiagram: diagram,
      };
      RAD.recordEvent(this.state, {
        type: "two-cell-radius-backlash-phase-diagram-case-applied",
        caseId: selected.caseId,
        phase: selected.phase,
        backlash,
        holeRadius,
        alphaCommand: selected.alphaCommand,
        zCommand: selected.zCommand,
      });
      this.syncControls();
      this.renderTwoCellBench();
      this.onChange(this.state);
      return selected;
    }

    renderTwoCellBench() {
      if (!this.els.twoCellAlpha) return;
      const controls =
        typeof RAD.twoCellBenchControls === "function"
          ? RAD.twoCellBenchControls(this.state, this.state.experiment?.twoCellBenchControls || {})
          : this.state.experiment?.twoCellBenchControls || {};
      this.els.twoCellAlpha.value = controls.alphaCommand ?? -0.35;
      this.els.twoCellZ.value = controls.zCommand ?? 0.35;
      this.els.twoCellHoleMax.value = controls.holeSweepMax ?? 0.5;
      this.els.twoCellSweepSteps.value = controls.holeSweepSteps ?? 9;
      this.els.twoCellLeftPositionLocked.checked = controls.leftPositionLocked !== false;
      this.els.twoCellRightLocked.checked = controls.rightLocked === true;
      this.els.twoCellRightPositionLocked.checked = controls.rightPositionLocked === true;
      this.syncTwoCellSuiteCases(controls);
      document.getElementById("twoCellAlphaOut").textContent = Number(controls.alphaCommand ?? -0.35).toFixed(2);
      document.getElementById("twoCellZOut").textContent = Number(controls.zCommand ?? 0.35).toFixed(2);
      document.getElementById("twoCellHoleMaxOut").textContent = Number(controls.holeSweepMax ?? 0.5).toFixed(3);
      document.getElementById("twoCellSweepStepsOut").textContent = String(Math.round(Number(controls.holeSweepSteps ?? 9)));
      const bench = this.state.experiment?.twoCellBench?.bench;
      const sweep = this.state.experiment?.twoCellBench?.sweep;
      const contact = this.state.experiment?.twoCellBench?.contact;
      const suite = this.state.experiment?.twoCellPhysicalSuite || this.state.experiment?.twoCellBench?.suite;
      const matrix = this.state.experiment?.twoCellPhysicalFidelityMatrix || this.state.experiment?.twoCellBench?.fidelityMatrix;
      const phaseMap = this.state.experiment?.twoCellContactPhaseMap || this.state.experiment?.twoCellBench?.contactPhaseMap;
      const phaseDiagram = this.state.experiment?.twoCellRadiusBacklashPhaseDiagram || this.state.experiment?.twoCellBench?.radiusBacklashPhaseDiagram;
      const transitionReport = this.state.experiment?.twoCellRadiusBacklashTransitionReport || this.state.experiment?.twoCellBench?.radiusBacklashTransitionReport;
      const transitionComparison = this.state.experiment?.twoCellRadiusBacklashTransitionComparison;
      const transitionRerun = this.state.experiment?.twoCellRadiusBacklashTransitionRerun;
      const cadContactDecomposition = this.state.experiment?.twoCellCadContactDecomposition;
      const exactContactPlan = this.state.experiment?.twoCellExactContactHandoffPlan;
      this.renderTwoCellPhaseDiagramGrid(phaseDiagram);
      const atlas = this.state.experiment?.twoCellPhysicalResponseAtlas || this.state.experiment?.twoCellBench?.responseAtlas;
      const comparison = this.state.experiment?.twoCellFidelityMatrixMeasurementComparison;
      const calibration = this.state.experiment?.twoCellFidelityMatrixParameterCalibration;
      const externalSummary =
        typeof RAD.twoCellExternalFidelityWebSummary === "function" ? RAD.twoCellExternalFidelityWebSummary() : null;
      const externalCase = this.syncTwoCellExternalCases(externalSummary);
      const atlasCase = this.syncTwoCellAtlasCases(atlas);
      const status =
        typeof RAD.twoCellPhysicalFidelityStatus === "function"
          ? RAD.twoCellPhysicalFidelityStatus(this.state, controls, { suite, matrix, comparison, calibration })
          : null;
      if (this.els.twoCellExternalFidelitySummary) {
        const ext = externalSummary?.summary;
        this.els.twoCellExternalFidelitySummary.textContent =
          ext?.available
            ? `MuJoCo ${ext.caseCount || 0} cases/${ext.connectorRowCount || 0} connectors, contact ${Number(100 * (ext.contactModeAccuracy || 0)).toFixed(1)}%, lock ${Number(100 * (ext.lockHeldAccuracy || 0)).toFixed(1)}%`
            : "external MuJoCo summary not exported";
      }
      if (this.els.twoCellCorrectionSummary) {
        const ext = externalSummary?.summary;
        this.els.twoCellCorrectionSummary.textContent =
          ext?.available
            ? `corrected RMS ${Number(ext.baselineRms || 0).toFixed(3)} -> ${Number(ext.correctedRms || 0).toFixed(3)}, improvement ${Number(100 * (ext.improvement || 0)).toFixed(1)}%, physical exact ${ext.physicalAccuracyValidated ? "yes" : "no"}`
            : "corrected preview --";
      }
      if (this.els.twoCellExternalCaseSummary) {
        const active = this.state.experiment?.twoCellExternalCasePreview || externalCase;
        this.els.twoCellExternalCaseSummary.textContent = active
          ? `case ${active.caseId}, ${active.actuationCase || ""} ${active.lockMode || ""}, b ${Number(active.grid?.backlash ?? active.backlash ?? 0).toFixed(2)}, hole ${Number(active.grid?.holeRadius ?? active.holeRadius ?? 0).toFixed(3)}`
          : "external case --";
      }
      if (this.els.twoCellExternalCaseDetail) {
        const active = this.state.experiment?.twoCellExternalCasePreview;
        const selected = active || externalCase;
        const correctedZ = active?.correctedPreview?.rightCellZ ?? selected?.correctedRightCellZ;
        const observedZ = active?.observedProxy?.rightCellZ ?? selected?.observedRightCellZ;
        const correctedSlip = active?.correctedPreview?.verticalSlipMm ?? selected?.correctedVerticalSlipMm;
        const observedSlip = active?.observedProxy?.verticalSlipMm ?? selected?.observedVerticalSlipMm;
        this.els.twoCellExternalCaseDetail.textContent = selected
          ? `right z corr ${Number(correctedZ || 0).toFixed(3)} / obs ${Number(observedZ || 0).toFixed(3)}, v-slip corr ${Number(correctedSlip || 0).toFixed(2)} / obs ${Number(observedSlip || 0).toFixed(2)} mm`
          : "external response --";
      }
      if (this.els.twoCellPhysicalStatusSummary) {
        this.els.twoCellPhysicalStatusSummary.textContent = status
          ? `exact ${status.summary.exactGeometryReady ? "ready" : "not ready"}, proxy ${status.summary.reducedProxyReady ? "ready" : "pending"}, missing ${status.summary.missingEvidenceCount}`
          : "exact CAD status --";
      }
      if (this.els.twoCellArchiveSummary) {
        this.els.twoCellArchiveSummary.textContent = status
          ? `archive BREP ${status.cadArchiveAudit.summary.brepEntryCount}, preview ${status.cadArchiveAudit.summary.previewEntryCount}, rows ${status.summary.matrixMeasurementRowsExpected}`
          : "archive --";
      }
      if (this.els.twoCellCadContactDecompositionSummary) {
        const summary = cadContactDecomposition?.summary;
        this.els.twoCellCadContactDecompositionSummary.textContent = summary
          ? `CAD contact ${summary.holeCount || 0} holes/${summary.twoCellConnectorCount || 0} connectors, ${summary.collisionPrimitiveCount || 0} primitives`
          : "CAD contact not built";
      }
      if (this.els.twoCellExactContactPlanSummary) {
        this.els.twoCellExactContactPlanSummary.textContent = exactContactPlan
          ? `contact plan ${exactContactPlan.summary?.caseCount || 0} cases, exact ${exactContactPlan.summary?.canRunExactContact ? "ready" : "blocked"}, missing ${exactContactPlan.summary?.missingEvidenceCount || 0}`
          : "contact plan not built";
      }
      if (this.els.twoCellMatrixComparisonSummary) {
        this.els.twoCellMatrixComparisonSummary.textContent = comparison
          ? `measured ${comparison.summary?.observedScalarCount || 0} values, rms ${comparison.summary?.rmsError === null ? "--" : Number(comparison.summary?.rmsError || 0).toFixed(4)}, ${comparison.summary?.passesTolerance ? "pass" : "review"}`
          : "measurements not loaded";
      }
      if (this.els.twoCellMatrixCalibrationSummary) {
        const zScale = calibration?.estimates?.zResponseScale?.estimate;
        const slipScale = calibration?.estimates?.verticalSlipScale?.estimate;
        const updates = calibration?.proposedReducedProxyUpdates;
        this.els.twoCellMatrixCalibrationSummary.textContent = calibration
          ? `cal z-scale ${zScale === null || zScale === undefined ? "--" : Number(zScale).toFixed(2)}, slip-scale ${slipScale === null || slipScale === undefined ? "--" : Number(slipScale).toFixed(2)}, proposed zGain ${Number(updates?.zCouplingGain || 0).toFixed(2)}, hole ${Number(updates?.holeRadius || 0).toFixed(3)}`
          : "calibration not estimated";
      }
      if (this.els.twoCellResponseAtlasSummary) {
        this.els.twoCellResponseAtlasSummary.textContent = atlas
          ? `atlas ${atlas.matrixSummary?.rowCount || 0} cases, priority ${atlas.benchPriority?.firstCasesToMeasure?.length || 0}, max slip ${Number(atlas.connectorSummary?.maxVerticalSlipMm || 0).toFixed(2)} mm, exact ${atlas.claimBoundary?.physicalAccuracyValidated ? "yes" : "no"}`
          : "atlas not run";
      }
      if (this.els.twoCellResponseAtlasCaseSummary) {
        const active = this.state.experiment?.twoCellResponseAtlasCasePreview || atlasCase;
        this.els.twoCellResponseAtlasCaseSummary.textContent = active
          ? `atlas ${active.caseId}, ${active.actuationCase} ${active.lockMode}, b ${Number(active.backlash || 0).toFixed(2)}, hole ${Number(active.holeRadius || 0).toFixed(3)}, z ${Number(active.rightZ || 0).toFixed(3)}`
          : "atlas case --";
      }
      if (!bench) {
        this.els.twoCellBenchSummary.textContent = "two-cell bench not run";
        this.els.twoCellContactSummary.textContent = "contact --";
        this.els.twoCellEnergySummary.textContent = "energy --";
        this.els.twoCellSweepSummary.textContent = "sweep --";
        if (this.els.twoCellSuiteSummary) {
          this.els.twoCellSuiteSummary.textContent = suite
            ? `suite ${suite.summary?.caseCount || 0} cases, ${suite.summary?.internalSuiteReady ? "internal ready" : "review"}`
            : "suite not run";
        }
        if (this.els.twoCellFidelityMatrixSummary) {
          this.els.twoCellFidelityMatrixSummary.textContent = matrix
            ? `matrix ${matrix.summary?.rowCount || 0} rows, ${matrix.summary?.connectorRowCount || 0} connector rows, missing ${matrix.summary?.missingEvidence?.length || 0}`
            : "matrix not run";
        }
        if (this.els.twoCellPhaseMapSummary) {
          this.els.twoCellPhaseMapSummary.textContent = phaseMap
            ? `phase ${phaseMap.summary?.dominantPhase || "none"}, contacts ${phaseMap.summary?.activeContactCaseCount || 0}/${phaseMap.summary?.rowCount || 0}, free ${phaseMap.summary?.freePlayCaseCount || 0}, exact ${phaseMap.summary?.physicalAccuracyValidated ? "yes" : "no"}`
            : "phase map not run";
        }
        if (this.els.twoCellPhaseDiagramSummary) {
          this.els.twoCellPhaseDiagramSummary.textContent = phaseDiagram
            ? `diagram ${phaseDiagram.summary?.holeStepCount || 0}x${phaseDiagram.summary?.backlashStepCount || 0}, phase ${phaseDiagram.summary?.dominantPhase || "none"}, transitions ${phaseDiagram.summary?.transitionCount || 0}, exact ${phaseDiagram.summary?.physicalAccuracyValidated ? "yes" : "no"}`
            : "phase diagram not run";
        }
        if (this.els.twoCellTransitionReportSummary) {
          this.els.twoCellTransitionReportSummary.textContent = transitionReport
            ? `transition brackets ${transitionReport.summary?.transitionBracketCount || 0}, measure ${transitionReport.summary?.measurementCaseCount || 0}, exact ${transitionReport.summary?.physicalAccuracyValidated ? "yes" : "no"}`
            : "transition report not run";
        }
        if (this.els.twoCellTransitionComparisonSummary) {
          this.els.twoCellTransitionComparisonSummary.textContent = transitionComparison
            ? `transition measured ${transitionComparison.summary?.matchedRowCount || 0}/${transitionComparison.summary?.providedRowCount || 0}, phase ${transitionComparison.summary?.phaseAccuracy === null ? "--" : Number(100 * (transitionComparison.summary?.phaseAccuracy || 0)).toFixed(1)}%, shift ${transitionComparison.summary?.transitionShiftDirection || "--"}`
            : "transition measurements not loaded";
        }
        if (this.els.twoCellTransitionRerunSummary) {
          this.els.twoCellTransitionRerunSummary.textContent = transitionRerun
            ? `rerun bias h ${Number(transitionRerun.axisBiases?.effectiveHoleRadiusBias || 0).toFixed(4)}, b ${Number(transitionRerun.axisBiases?.effectiveBacklashBias || 0).toFixed(4)}, moved ${Number(transitionRerun.summary?.maxAbsMidpointDelta || 0).toFixed(4)}`
            : "calibrated rerun not built";
        }
        return;
      }
      const right = bench.cells?.[1] || {};
      const left = bench.cells?.[0] || {};
      this.els.twoCellBenchSummary.textContent = `left z ${Number(left.residualZ || 0).toFixed(3)}, right z ${Number(right.center?.z || 0).toFixed(3)}`;
      this.els.twoCellContactSummary.textContent = contact
        ? `connectors ${contact.summary?.activeContactCount || 0}/${contact.summary?.connectorCount || 0}, max slip ${Number(contact.summary?.maxVerticalSlipMm || 0).toFixed(2)} mm`
        : `${bench.connector.contactMode}, clearance ${Number(bench.dimensions.pinHoleClearance || 0).toFixed(3)}`;
      this.els.twoCellEnergySummary.textContent = `energy ${Number(bench.energy.totalEnergy || 0).toFixed(4)}, strain ${Number(bench.connector.axialStrain || 0).toFixed(4)}`;
      this.els.twoCellSweepSummary.textContent = sweep
        ? `sweep ${sweep.rows.length} holes, neighbor z ${Number(sweep.trend.neighborZStart || 0).toFixed(3)} -> ${Number(sweep.trend.neighborZEnd || 0).toFixed(3)}`
        : "sweep --";
      if (this.els.twoCellSuiteSummary) {
        this.els.twoCellSuiteSummary.textContent = suite
          ? `suite ${suite.summary?.caseCount || 0} cases, locks ${suite.summary?.lockPassCount || 0}/${suite.summary?.lockCaseCount || 0}, missing ${suite.summary?.missingEvidence?.length || 0}`
          : "suite not run";
      }
      if (this.els.twoCellFidelityMatrixSummary) {
        this.els.twoCellFidelityMatrixSummary.textContent = matrix
          ? `matrix ${matrix.summary?.rowCount || 0} rows, locks ${matrix.summary?.lockChecks?.stateLockVerticalMotionSamples || 0} z-free, missing ${matrix.summary?.missingEvidence?.length || 0}`
          : "matrix not run";
      }
      if (this.els.twoCellPhaseMapSummary) {
        this.els.twoCellPhaseMapSummary.textContent = phaseMap
          ? `phase ${phaseMap.summary?.dominantPhase || "none"}, contacts ${phaseMap.summary?.activeContactCaseCount || 0}/${phaseMap.summary?.rowCount || 0}, free ${phaseMap.summary?.freePlayCaseCount || 0}, exact ${phaseMap.summary?.physicalAccuracyValidated ? "yes" : "no"}`
          : "phase map not run";
      }
      if (this.els.twoCellPhaseDiagramSummary) {
        this.els.twoCellPhaseDiagramSummary.textContent = phaseDiagram
          ? `diagram ${phaseDiagram.summary?.holeStepCount || 0}x${phaseDiagram.summary?.backlashStepCount || 0}, phase ${phaseDiagram.summary?.dominantPhase || "none"}, transitions ${phaseDiagram.summary?.transitionCount || 0}, exact ${phaseDiagram.summary?.physicalAccuracyValidated ? "yes" : "no"}`
          : "phase diagram not run";
      }
      if (this.els.twoCellTransitionReportSummary) {
        this.els.twoCellTransitionReportSummary.textContent = transitionReport
          ? `transition brackets ${transitionReport.summary?.transitionBracketCount || 0}, measure ${transitionReport.summary?.measurementCaseCount || 0}, exact ${transitionReport.summary?.physicalAccuracyValidated ? "yes" : "no"}`
          : "transition report not run";
      }
      if (this.els.twoCellTransitionComparisonSummary) {
        this.els.twoCellTransitionComparisonSummary.textContent = transitionComparison
          ? `transition measured ${transitionComparison.summary?.matchedRowCount || 0}/${transitionComparison.summary?.providedRowCount || 0}, phase ${transitionComparison.summary?.phaseAccuracy === null ? "--" : Number(100 * (transitionComparison.summary?.phaseAccuracy || 0)).toFixed(1)}%, shift ${transitionComparison.summary?.transitionShiftDirection || "--"}`
          : "transition measurements not loaded";
      }
      if (this.els.twoCellTransitionRerunSummary) {
        this.els.twoCellTransitionRerunSummary.textContent = transitionRerun
          ? `rerun bias h ${Number(transitionRerun.axisBiases?.effectiveHoleRadiusBias || 0).toFixed(4)}, b ${Number(transitionRerun.axisBiases?.effectiveBacklashBias || 0).toFixed(4)}, moved ${Number(transitionRerun.summary?.maxAbsMidpointDelta || 0).toFixed(4)}`
          : "calibrated rerun not built";
      }
      if (this.els.twoCellMatrixComparisonSummary) {
        this.els.twoCellMatrixComparisonSummary.textContent = comparison
          ? `measured ${comparison.summary?.observedScalarCount || 0} values, rms ${comparison.summary?.rmsError === null ? "--" : Number(comparison.summary?.rmsError || 0).toFixed(4)}, ${comparison.summary?.passesTolerance ? "pass" : "review"}`
          : "measurements not loaded";
      }
    }

    saveTwoCellBench() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellBenchJson === "function"
          ? RAD.exportTwoCellBenchJson(this.state, controls)
          : JSON.stringify({ schema: "rad-sim.two-cell-physical-bench.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-bench.json");
    }

    saveTwoCellSweep() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellBacklashSweepCsv === "function"
          ? RAD.exportTwoCellBacklashSweepCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-backlash-sweep.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-backlash-sweep.csv", "text/csv");
    }

    saveTwoCellSuite() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellPhysicalSuiteJson === "function"
          ? RAD.exportTwoCellPhysicalSuiteJson(this.state, controls)
          : JSON.stringify({ schema: "rad-sim.two-cell-physical-simulation-suite.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-physical-suite.json");
    }

    saveTwoCellSuiteCsv() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellPhysicalSuiteCsv === "function"
          ? RAD.exportTwoCellPhysicalSuiteCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-physical-simulation-suite.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-physical-suite.csv", "text/csv");
    }

    saveTwoCellFidelityMatrix() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellPhysicalFidelityMatrixJson === "function"
          ? RAD.exportTwoCellPhysicalFidelityMatrixJson(this.state, controls)
          : JSON.stringify({ schema: "rad-sim.two-cell-physical-fidelity-matrix.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-physical-fidelity-matrix.json");
    }

    saveTwoCellFidelityMatrixCsv() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellPhysicalFidelityMatrixCsv === "function"
          ? RAD.exportTwoCellPhysicalFidelityMatrixCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-physical-fidelity-matrix.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-physical-fidelity-matrix.csv", "text/csv");
    }

    saveTwoCellPhaseMap() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellContactPhaseMapJson === "function"
          ? RAD.exportTwoCellContactPhaseMapJson(this.state, controls)
          : JSON.stringify({ schema: "rad-sim.two-cell-contact-phase-map.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-contact-phase-map.json");
    }

    saveTwoCellPhaseMapCsv() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellContactPhaseMapCsv === "function"
          ? RAD.exportTwoCellContactPhaseMapCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-contact-phase-map.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-contact-phase-map.csv", "text/csv");
    }

    saveTwoCellPhaseDiagram() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellRadiusBacklashPhaseDiagramJson === "function"
          ? RAD.exportTwoCellRadiusBacklashPhaseDiagramJson(this.state, controls)
          : JSON.stringify({ schema: "rad-sim.two-cell-radius-backlash-phase-diagram.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-radius-backlash-phase-diagram.json");
    }

    saveTwoCellPhaseDiagramCsv() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellRadiusBacklashPhaseDiagramCsv === "function"
          ? RAD.exportTwoCellRadiusBacklashPhaseDiagramCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-radius-backlash-phase-diagram.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-radius-backlash-phase-diagram.csv", "text/csv");
    }

    saveTwoCellTransitionReport() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellRadiusBacklashTransitionReportJson === "function"
          ? RAD.exportTwoCellRadiusBacklashTransitionReportJson(this.state, controls)
          : JSON.stringify({ schema: "rad-sim.two-cell-radius-backlash-transition-report.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-radius-backlash-transition-report.json");
    }

    saveTwoCellTransitionReportCsv() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellRadiusBacklashTransitionReportCsv === "function"
          ? RAD.exportTwoCellRadiusBacklashTransitionReportCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-radius-backlash-transition-report.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-radius-backlash-transition-report.csv", "text/csv");
    }

    loadTwoCellTransitionMeasurements() {
      const file = this.els.twoCellTransitionMeasurementFileInput.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        try {
          const controls = this.updateTwoCellBenchControls();
          const comparison =
            typeof RAD.compareTwoCellRadiusBacklashTransitionMeasurements === "function"
              ? RAD.compareTwoCellRadiusBacklashTransitionMeasurements(this.state, String(reader.result), controls)
              : null;
          const rerun =
            typeof RAD.twoCellRadiusBacklashTransitionRerun === "function"
              ? RAD.twoCellRadiusBacklashTransitionRerun(this.state, String(reader.result), controls)
              : null;
          this.state.experiment.twoCellRadiusBacklashTransitionComparison = comparison;
          this.state.experiment.twoCellRadiusBacklashTransitionRerun = rerun;
          RAD.recordEvent(this.state, {
            type: "two-cell-radius-backlash-transition-measurements-load",
            rows: comparison?.summary?.providedRowCount || 0,
            matched: comparison?.summary?.matchedRowCount || 0,
            phaseAccuracy: comparison?.summary?.phaseAccuracy ?? null,
            shift: comparison?.summary?.transitionShiftDirection || "unknown",
            maxAbsRerunMidpointDelta: rerun?.summary?.maxAbsMidpointDelta ?? null,
          });
          this.renderTwoCellBench();
          this.onChange(this.state);
        } catch (error) {
          window.alert(error.message);
        }
      };
      reader.readAsText(file);
      this.els.twoCellTransitionMeasurementFileInput.value = "";
    }

    saveTwoCellTransitionComparison() {
      const comparison = this.state.experiment?.twoCellRadiusBacklashTransitionComparison;
      if (!comparison) {
        window.alert("Load filled transition measurements before saving a comparison.");
        return;
      }
      this.downloadText(JSON.stringify(comparison, null, 2), "rad-sim-two-cell-radius-backlash-transition-comparison.json");
    }

    saveTwoCellTransitionRerun() {
      const rerun = this.state.experiment?.twoCellRadiusBacklashTransitionRerun;
      if (!rerun) {
        window.alert("Load filled transition measurements before saving a calibrated rerun.");
        return;
      }
      this.downloadText(JSON.stringify(rerun, null, 2), "rad-sim-two-cell-radius-backlash-transition-rerun.json");
    }

    applyTwoCellTransitionCalibration() {
      const comparison = this.state.experiment?.twoCellRadiusBacklashTransitionComparison;
      const updates = comparison?.calibrationEstimate?.proposedReducedProxyUpdates;
      if (!updates) {
        window.alert("Load filled transition measurements before applying a transition calibration estimate.");
        return;
      }
      this.state.grid.couplingGain = Math.max(0, Math.min(1, Number(updates.couplingGain ?? this.state.grid.couplingGain)));
      this.state.grid.zCouplingGain = Math.max(0, Math.min(1, Number(updates.zCouplingGain ?? this.state.grid.zCouplingGain)));
      this.state.grid.backlash = Math.max(0, Number(updates.backlash ?? this.state.grid.backlash));
      this.state.grid.holeRadius = Math.max(Number(this.state.grid.pinRadius || 0), Number(updates.holeRadius ?? this.state.grid.holeRadius));
      RAD.recordEvent(this.state, {
        type: "two-cell-radius-backlash-transition-calibration-applied",
        couplingGain: this.state.grid.couplingGain,
        zCouplingGain: this.state.grid.zCouplingGain,
        backlash: this.state.grid.backlash,
        holeRadius: this.state.grid.holeRadius,
        transitionBoundaryDirection: updates.transitionBoundaryDirection || comparison.summary?.transitionShiftDirection || "unknown",
        physicalAccuracyValidated: false,
      });
      this.renderTwoCellBench();
      this.syncControls();
      this.onChange(this.state);
    }

    saveTwoCellResponseAtlas() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellPhysicalResponseAtlasJson === "function"
          ? RAD.exportTwoCellPhysicalResponseAtlasJson(this.state, controls)
          : JSON.stringify({ schema: "rad-sim.two-cell-physical-response-atlas.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-physical-response-atlas.json");
    }

    saveTwoCellResponseAtlasCsv() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellPhysicalResponseAtlasCsv === "function"
          ? RAD.exportTwoCellPhysicalResponseAtlasCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-physical-response-atlas.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-physical-response-atlas.csv", "text/csv");
    }

    saveTwoCellFidelityMatrixTemplateCsv() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellFidelityMatrixMeasurementTemplateCsv === "function"
          ? RAD.exportTwoCellFidelityMatrixMeasurementTemplateCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-fidelity-matrix-measurement-template.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-fidelity-matrix-measurement-template.csv", "text/csv");
    }

    saveTwoCellExternalFidelityManifest() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellExternalFidelityMatrixManifestJson === "function"
          ? RAD.exportTwoCellExternalFidelityMatrixManifestJson(this.state, controls)
          : JSON.stringify({ schema: "rad-sim.browser-two-cell-external-fidelity-matrix-manifest.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-external-fidelity-matrix-manifest.json");
    }

    buildTwoCellCadContactDecomposition() {
      const controls = this.updateTwoCellBenchControls();
      const spec =
        typeof RAD.twoCellCadContactDecompositionSpec === "function"
          ? RAD.twoCellCadContactDecompositionSpec(this.state, controls)
          : { schema: "rad-sim.two-cell-cad-contact-decomposition.missing", summary: { holeCount: 0, twoCellConnectorCount: 0, collisionPrimitiveCount: 0 } };
      this.state.experiment.twoCellCadContactDecomposition = spec;
      RAD.recordEvent(this.state, {
        type: "two-cell-cad-contact-decomposition-built",
        holes: spec.summary?.holeCount || 0,
        connectors: spec.summary?.twoCellConnectorCount || 0,
        primitives: spec.summary?.collisionPrimitiveCount || 0,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return { controls, spec };
    }

    saveTwoCellCadContactDecomposition() {
      const { spec } = this.buildTwoCellCadContactDecomposition();
      this.downloadText(JSON.stringify(spec, null, 2), "rad-sim-two-cell-cad-contact-decomposition.json");
    }

    saveTwoCellCadContactDecompositionCsv() {
      const { controls } = this.buildTwoCellCadContactDecomposition();
      const payload =
        typeof RAD.exportTwoCellCadContactDecompositionCsv === "function"
          ? RAD.exportTwoCellCadContactDecompositionCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-cad-contact-decomposition.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-cad-contact-decomposition.csv", "text/csv");
    }

    buildTwoCellExactContactHandoffPlan() {
      const controls = this.updateTwoCellBenchControls();
      const plan =
        typeof RAD.twoCellExactContactHandoffPlan === "function"
          ? RAD.twoCellExactContactHandoffPlan(this.state, controls)
          : { schema: "rad-sim.two-cell-exact-contact-handoff-plan.missing", summary: { caseCount: 0, missingEvidenceCount: 1 } };
      this.state.experiment.twoCellExactContactHandoffPlan = plan;
      RAD.recordEvent(this.state, {
        type: "two-cell-exact-contact-handoff-plan-built",
        cases: plan.summary?.caseCount || 0,
        canRunExactContact: plan.summary?.canRunExactContact || false,
        missingEvidenceCount: plan.summary?.missingEvidenceCount || 0,
      });
      this.renderTwoCellBench();
      this.onChange(this.state);
      return { controls, plan };
    }

    saveTwoCellExactContactHandoffPlan() {
      const { plan } = this.buildTwoCellExactContactHandoffPlan();
      this.downloadText(JSON.stringify(plan, null, 2), "rad-sim-two-cell-exact-contact-handoff-plan.json");
    }

    saveTwoCellExactContactHandoffPlanCsv() {
      const { controls } = this.buildTwoCellExactContactHandoffPlan();
      const payload =
        typeof RAD.exportTwoCellExactContactHandoffPlanCsv === "function"
          ? RAD.exportTwoCellExactContactHandoffPlanCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-exact-contact-handoff-plan.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-exact-contact-handoff-plan.csv", "text/csv");
    }

    saveTwoCellExternalFidelitySummary() {
      const payload =
        typeof RAD.exportTwoCellExternalFidelityWebSummaryJson === "function"
          ? RAD.exportTwoCellExternalFidelityWebSummaryJson()
          : JSON.stringify({ schema: "rad-sim.web-two-cell-external-fidelity-summary.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-external-fidelity-summary.json");
    }

    saveTwoCellConnectorTemplateCsv() {
      const controls = this.updateTwoCellBenchControls();
      const payload =
        typeof RAD.exportTwoCellConnectorMeasurementTemplateCsv === "function"
          ? RAD.exportTwoCellConnectorMeasurementTemplateCsv(this.state, controls)
          : "schema\nrad-sim.two-cell-connector-measurement-template.missing\n";
      this.downloadText(payload, "rad-sim-two-cell-connector-measurement-template.csv", "text/csv");
    }

    loadTwoCellFidelityMatrixMeasurements() {
      const file = this.els.twoCellFidelityMatrixMeasurementFileInput.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        try {
          const controls = this.updateTwoCellBenchControls();
          const comparison =
            typeof RAD.compareTwoCellFidelityMatrixMeasurements === "function"
              ? RAD.compareTwoCellFidelityMatrixMeasurements(this.state, String(reader.result), controls)
              : null;
          const calibration =
            typeof RAD.calibrateTwoCellFidelityMatrixParameters === "function"
              ? RAD.calibrateTwoCellFidelityMatrixParameters(this.state, String(reader.result), controls)
              : null;
          this.state.experiment.twoCellFidelityMatrixMeasurementComparison = comparison;
          this.state.experiment.twoCellFidelityMatrixParameterCalibration = calibration;
          RAD.recordEvent(this.state, {
            type: "two-cell-fidelity-matrix-measurements-load",
            rows: comparison?.summary?.providedRowCount || 0,
            observed: comparison?.summary?.observedScalarCount || 0,
            pass: comparison?.summary?.passesTolerance || false,
            calibrationReady: calibration?.summary?.readyForReducedProxyCalibration || false,
          });
          this.renderTwoCellBench();
          this.onChange(this.state);
        } catch (error) {
          window.alert(error.message);
        }
      };
      reader.readAsText(file);
      this.els.twoCellFidelityMatrixMeasurementFileInput.value = "";
    }

    saveTwoCellFidelityMatrixComparison() {
      const comparison = this.state.experiment?.twoCellFidelityMatrixMeasurementComparison;
      if (!comparison) {
        window.alert("Load filled matrix measurements before saving a comparison.");
        return;
      }
      this.downloadText(JSON.stringify(comparison, null, 2), "rad-sim-two-cell-fidelity-matrix-measurement-comparison.json");
    }

    saveTwoCellFidelityMatrixCalibration() {
      const calibration = this.state.experiment?.twoCellFidelityMatrixParameterCalibration;
      if (!calibration) {
        window.alert("Load filled matrix measurements before saving a calibration estimate.");
        return;
      }
      this.downloadText(JSON.stringify(calibration, null, 2), "rad-sim-two-cell-fidelity-matrix-parameter-calibration.json");
    }

    applyTwoCellFidelityMatrixCalibration() {
      const calibration = this.state.experiment?.twoCellFidelityMatrixParameterCalibration;
      const updates = calibration?.proposedReducedProxyUpdates;
      if (!updates) {
        window.alert("Load filled matrix measurements before applying a calibration estimate.");
        return;
      }
      this.state.grid.couplingGain = Math.max(0, Math.min(1, Number(updates.couplingGain ?? this.state.grid.couplingGain)));
      this.state.grid.zCouplingGain = Math.max(0, Math.min(1, Number(updates.zCouplingGain ?? this.state.grid.zCouplingGain)));
      this.state.grid.holeRadius = Math.max(Number(this.state.grid.pinRadius || 0), Number(updates.holeRadius ?? this.state.grid.holeRadius));
      this.state.grid.backlash = Math.max(0, Number(this.state.grid.backlash || 0));
      RAD.recordEvent(this.state, {
        type: "two-cell-fidelity-matrix-calibration-applied",
        couplingGain: this.state.grid.couplingGain,
        zCouplingGain: this.state.grid.zCouplingGain,
        holeRadius: this.state.grid.holeRadius,
        physicalAccuracyValidated: false,
      });
      this.renderTwoCellBench();
      this.syncControls();
      this.onChange(this.state);
    }

    saveTwoCellPhysicalStatus() {
      const controls = this.updateTwoCellBenchControls();
      const suite = this.state.experiment?.twoCellPhysicalSuite || this.state.experiment?.twoCellBench?.suite;
      const matrix = this.state.experiment?.twoCellPhysicalFidelityMatrix || this.state.experiment?.twoCellBench?.fidelityMatrix;
      const comparison = this.state.experiment?.twoCellFidelityMatrixMeasurementComparison;
      const calibration = this.state.experiment?.twoCellFidelityMatrixParameterCalibration;
      const payload =
        typeof RAD.exportTwoCellPhysicalFidelityStatusJson === "function"
          ? RAD.exportTwoCellPhysicalFidelityStatusJson(this.state, controls, { suite, matrix, comparison, calibration })
          : JSON.stringify({ schema: "rad-sim.browser-two-cell-physical-fidelity-status.missing" }, null, 2);
      this.downloadText(payload, "rad-sim-two-cell-physical-fidelity-status.json");
    }

    hardwareProfileNumberInputs() {
      return [
        ["backlashMm", this.els.hardwareBacklashMm],
        ["pinRadiusMm", this.els.hardwarePinRadiusMm],
        ["holeRadiusMm", this.els.hardwareHoleRadiusMm],
        ["plateThicknessMm", this.els.hardwarePlateThicknessMm],
        ["jointStackHeightMm", this.els.hardwareJointStackHeightMm],
        ["bossRadiusMm", this.els.hardwareBossRadiusMm],
      ];
    }

    ensureHardwareProfile() {
      const defaults = typeof RAD.defaultHardwareProfile === "function" ? RAD.defaultHardwareProfile() : {};
      this.state.grid.hardwareProfile = {
        ...defaults,
        ...(this.state.grid.hardwareProfile || {}),
        sideLengthMm: Math.max(1e-9, Number(this.state.grid.paperSideLengthMm ?? 35)),
        fabricationHoleToleranceMm: Math.max(0, Number(this.state.grid.paperHoleToleranceMm ?? 0.1)),
      };
      return this.state.grid.hardwareProfile;
    }

    optionalProfileNumber(value) {
      if (value === "" || value === null || value === undefined) return null;
      const numeric = Number(value);
      return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
    }

    profileNumberInputValue(value) {
      return value === null || value === undefined ? "" : String(value);
    }

    enforceHardwareProfileRadii(profile) {
      if (
        profile.pinRadiusMm !== null &&
        profile.pinRadiusMm !== undefined &&
        profile.holeRadiusMm !== null &&
        profile.holeRadiusMm !== undefined &&
        profile.holeRadiusMm < profile.pinRadiusMm
      ) {
        profile.holeRadiusMm = profile.pinRadiusMm;
        this.els.hardwareHoleRadiusMm.value = this.profileNumberInputValue(profile.holeRadiusMm);
      }
    }

    updateHardwareProfileText(field, value) {
      const profile = this.ensureHardwareProfile();
      profile[field] = String(value || "").trim() || (field === "name" ? "paper-reference" : "unspecified");
      this.updateLabels(null);
      this.onChange(this.state);
    }

    updateHardwareProfileNumber(field, value) {
      const profile = this.ensureHardwareProfile();
      profile[field] = this.optionalProfileNumber(value);
      this.enforceHardwareProfileRadii(profile);
      this.updateLabels(null);
      this.onChange(this.state);
    }

    boundaryNumberInputValue(value) {
      return value === null || value === undefined || !Number.isFinite(Number(value)) ? "" : String(Number(value));
    }

    readBoundaryNumber(element) {
      const value = String(element.value || "").trim();
      if (!value) return null;
      const numeric = Number(value);
      return Number.isFinite(numeric) ? numeric : null;
    }

    updateBoundaryControls() {
      if (typeof RAD.ensureBoundarySchema === "function") RAD.ensureBoundarySchema(this.state);
      const boundary = this.state.boundary;
      boundary.mode = this.els.boundaryMode.value;
      boundary.wallType = this.els.boundaryWallType.value;
      boundary.xMin = this.readBoundaryNumber(this.els.boundaryXMin);
      boundary.xMax = this.readBoundaryNumber(this.els.boundaryXMax);
      boundary.yMin = this.readBoundaryNumber(this.els.boundaryYMin);
      boundary.yMax = this.readBoundaryNumber(this.els.boundaryYMax);
      boundary.channelAxis = this.els.boundaryChannelAxis.value === "y" ? "y" : "x";
      boundary.channelWidth = Math.max(0, Number(this.els.boundaryChannelWidth.value) || 0);
      boundary.zGain = Math.max(0, Number(this.els.boundaryZGain.value) || 0);
      boundary.zThreshold = Math.max(0, Number(this.els.boundaryZThreshold.value) || 0);
      boundary.zPower = Math.max(0.5, Number(this.els.boundaryZPower.value) || 1.5);
      this.onChange(this.state);
    }

    applyHardwareProfile() {
      this.ensureHardwareProfile();
      if (typeof RAD.applyHardwareProfileToGrid === "function") RAD.applyHardwareProfileToGrid(this.state);
      if (this.state.grid.holeRadius < this.state.grid.pinRadius) this.state.grid.holeRadius = this.state.grid.pinRadius;
      this.syncControls();
      this.onChange(this.state);
    }

    applyA360CadProfile() {
      if (typeof RAD.cadRadCellReferenceProfile !== "function" || typeof RAD.hardwareProfileFromObject !== "function") return;
      const referenceProfile = RAD.cadRadCellReferenceProfile(this.state);
      this.state.grid.hardwareProfile = RAD.hardwareProfileFromObject(referenceProfile.hardwareProfile, this.state);
      if (typeof RAD.applyHardwareProfileToGrid === "function") RAD.applyHardwareProfileToGrid(this.state);
      if (this.state.grid.holeRadius < this.state.grid.pinRadius) this.state.grid.holeRadius = this.state.grid.pinRadius;
      this.state.view.cellVisualMode = "cadRad";
      this.state.view.gapsVisible = true;
      this.state.view.fastenersVisible = true;
      RAD.recordEvent(this.state, {
        type: "a360-cad-profile-applied",
        schema: referenceProfile.schema,
        source: referenceProfile.cadReference?.source,
        sideLengthMm: referenceProfile.hardwareProfile?.dimensionsMm?.sideLengthMm,
        holeRadiusMm: referenceProfile.hardwareProfile?.dimensionsMm?.holeRadiusMm,
      });
      this.syncControls();
      this.onChange(this.state);
    }

    clearHardwareProfileMeasurements() {
      const profile = this.ensureHardwareProfile();
      for (const [field, element] of this.hardwareProfileNumberInputs()) {
        profile[field] = null;
        element.value = "";
      }
      this.updateLabels(null);
      this.onChange(this.state);
    }

    saveHardwareProfile() {
      this.ensureHardwareProfile();
      const payload =
        typeof RAD.exportHardwareProfileJson === "function"
          ? RAD.exportHardwareProfileJson(this.state)
          : JSON.stringify(this.state.grid.hardwareProfile || {}, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "hardware_profile.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    loadHardwareProfile() {
      const file = this.els.hardwareProfileFileInput.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        try {
          const profile =
            typeof RAD.importHardwareProfileJson === "function"
              ? RAD.importHardwareProfileJson(String(reader.result), this.state)
              : JSON.parse(String(reader.result));
          this.state.grid.hardwareProfile = profile;
          this.state.grid.paperSideLengthMm = profile.sideLengthMm;
          this.state.grid.paperHoleToleranceMm = profile.fabricationHoleToleranceMm;
          if (typeof RAD.applyHardwareProfileToGrid === "function") RAD.applyHardwareProfileToGrid(this.state);
          if (this.state.grid.holeRadius < this.state.grid.pinRadius) this.state.grid.holeRadius = this.state.grid.pinRadius;
          RAD.recordEvent(this.state, {
            type: "hardware-profile-import",
            name: profile.name,
            measuredFields:
              typeof RAD.hardwareMeasuredFields === "function"
                ? RAD.hardwareMeasuredFields(profile)
                : [],
          });
          this.syncControls();
          this.onChange(this.state);
        } catch (error) {
          window.alert(error.message);
        }
      };
      reader.readAsText(file);
      this.els.hardwareProfileFileInput.value = "";
    }

    applySelected() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      const removed = this.els.removed.checked;
      this.state.cells.removed[r][c] = removed;
      if (removed) {
        this.state.cells.commandAlpha[r][c] = 0;
        this.state.cells.commandZ[r][c] = 0;
        this.state.cells.locked[r][c] = false;
        this.state.cells.positionLocked[r][c] = false;
        if (this.state.cells.lockZ) this.state.cells.lockZ[r][c] = 0;
        this.clearPositionLockTarget(r, c);
        this.state.cells.actuatorAllowed[r][c] = false;
      } else {
        this.state.cells.commandAlpha[r][c] = RAD.clampCommandAlpha(this.state, this.els.alphaCommand.value);
        this.state.cells.commandZ[r][c] = RAD.clampCommandZ(this.state, this.els.zCommand.value);
        if (this.els.locked.checked) this.commitDirectLockState(r, c);
        if (this.els.positionLocked.checked && !this.state.cells.positionLocked?.[r]?.[c]) this.commitPositionLockState(r, c);
        if (!this.els.positionLocked.checked) this.clearPositionLockTarget(r, c);
        this.state.cells.locked[r][c] = this.els.locked.checked;
        this.state.cells.positionLocked[r][c] = this.els.positionLocked.checked;
        this.state.cells.actuatorAllowed[r][c] = this.els.actuatorAllowed.checked;
      }
      RAD.recordEvent(this.state, {
        type: "cell-command",
        r,
        c,
        alpha: this.state.cells.commandAlpha[r][c],
        z: this.state.cells.commandZ[r][c],
        locked: this.state.cells.locked[r][c],
        positionLocked: this.state.cells.positionLocked[r][c],
        actuatorAllowed: this.state.cells.actuatorAllowed[r][c],
        removed: this.state.cells.removed[r][c],
      });
      this.onChange(this.state);
    }

    previewSelected() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      const removed = this.els.removed.checked;
      this.state.cells.removed[r][c] = removed;
      if (removed) {
        this.state.cells.commandAlpha[r][c] = 0;
        this.state.cells.commandZ[r][c] = 0;
        this.state.cells.locked[r][c] = false;
        this.state.cells.positionLocked[r][c] = false;
        if (this.state.cells.lockZ) this.state.cells.lockZ[r][c] = 0;
        this.clearPositionLockTarget(r, c);
        this.state.cells.actuatorAllowed[r][c] = false;
      } else {
        this.state.cells.commandAlpha[r][c] = RAD.clampCommandAlpha(this.state, this.els.alphaCommand.value);
        this.state.cells.commandZ[r][c] = RAD.clampCommandZ(this.state, this.els.zCommand.value);
        if (this.els.locked.checked) this.commitDirectLockState(r, c);
        if (this.els.positionLocked.checked && !this.state.cells.positionLocked?.[r]?.[c]) this.commitPositionLockState(r, c);
        if (!this.els.positionLocked.checked) this.clearPositionLockTarget(r, c);
        this.state.cells.locked[r][c] = this.els.locked.checked;
        this.state.cells.positionLocked[r][c] = this.els.positionLocked.checked;
      }
      this.state.inverse.physicalValidation = null;
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
      if (this.state.cells.removed?.[r]?.[c]) {
        this.state.cells.actuatorAllowed[r][c] = false;
        this.els.actuatorAllowed.checked = false;
      } else {
        this.state.cells.actuatorAllowed[r][c] = this.els.actuatorAllowed.checked;
      }
      this.state.inverse.plan = { candidates: [], commands: [], history: [] };
      this.state.inverse.preview = null;
      this.state.inverse.sensitivity = { candidates: [], map: RAD.matrix(this.state.grid.rows, this.state.grid.cols, 0), stepZ: 0.12, stepAlpha: 0.12, controllableCells: 0, meanGain: 0, maxGain: 0 };
      this.state.inverse.physicalValidation = null;
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
      this.state.inverse.physicalValidation = null;
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
          if (this.state.cells.removed?.[r]?.[c]) {
            this.state.cells.actuatorAllowed[r][c] = false;
            continue;
          }
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
      this.state.inverse.physicalValidation = null;
      RAD.recordEvent(this.state, { type: "actuator-mask", name, allowed });
      this.syncControls();
      this.onChange(this.state);
    }

    positionLockSelectedRowEnds() {
      const row = this.state.selection.r;
      const lastCol = this.state.grid.cols - 1;
      if (lastCol < 0) return;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      for (const c of [0, lastCol]) {
        if (this.state.cells.removed?.[row]?.[c] === true) {
          this.state.cells.positionLocked[row][c] = false;
          this.clearPositionLockTarget(row, c);
          continue;
        }
        if (!this.state.cells.positionLocked[row][c]) this.commitPositionLockState(row, c);
        this.state.cells.positionLocked[row][c] = true;
      }
      this.state.inverse.physicalValidation = null;
      RAD.recordEvent(this.state, { type: "position-lock-row-ends", row, first: 0, last: lastCol });
      this.syncControls();
      this.onChange(this.state);
    }

    clearSelected() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      this.state.cells.commandAlpha[r][c] = 0;
      this.state.cells.commandZ[r][c] = 0;
      this.state.cells.locked[r][c] = false;
      this.state.cells.positionLocked[r][c] = false;
      this.state.cells.removed[r][c] = false;
      this.state.cells.actuatorAllowed[r][c] = true;
      if (this.state.cells.lockAlpha) this.state.cells.lockAlpha[r][c] = this.state.grid.initialAlpha;
      if (this.state.cells.lockZ) this.state.cells.lockZ[r][c] = 0;
      this.clearPositionLockTarget(r, c);
      this.state.inverse.physicalValidation = null;
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
        lockZ: this.state.cells.lockZ?.[r]?.[c] ?? 0,
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

    removeSelectedCellEvent() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      this.state = RAD.applyProgrammableEvent(this.state, RAD.removeCellEvent({ r, c }));
      RAD.recordEvent(this.state, { type: "operator-remove-cell", r, c });
      this.syncControls();
      this.onChange(this.state);
    }

    restoreSelectedCellEvent() {
      const { r, c } = this.state.selection;
      if (!this.state.experiment.initialSnapshot) this.state.experiment.initialSnapshot = RAD.snapshotState(this.state);
      this.state = RAD.applyProgrammableEvent(this.state, RAD.restoreCellEvent({ r, c }));
      RAD.recordEvent(this.state, { type: "operator-restore-cell", r, c });
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
      if (base.cells.removed) base.cells.removed[r][c] = false;
      base.cells.lockAlpha[r][c] = base.grid.initialAlpha;
      if (!base.cells.lockZ) base.cells.lockZ = RAD.matrix(base.grid.rows, base.grid.cols, 0);
      base.cells.lockZ[r][c] = 0;
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
      const sequence = RAD.compareSequenceOrder(this.selectedEventBaseline(r, c), [
        RAD.localActuationEvent({ r, c }, alpha, z),
        RAD.lockEvent({ r, c }),
        RAD.clearActuationEvent({ r, c }),
      ]);
      diagnostic.sequence = sequence;
      this.lastOperatorDiagnostic = diagnostic;
      RAD.recordEvent(this.state, {
        type: "operator-order-check",
        r,
        c,
        alpha,
        z,
        finalAlphaError: diagnostic.finalAlphaError,
        finalHeightError: diagnostic.finalHeightError,
        maxOrderError: sequence.maxOrderError,
        noncommutingAdjacentPairs: sequence.noncommutingAdjacentPairs,
      });
      this.updateOperatorInspector(diagnostic);
      this.onChange(this.state);
    }

    updateOperatorInspector(diagnostic) {
      if (!diagnostic) {
        this.els.operatorOrderState.textContent = "not checked";
        this.els.operatorAlphaError.textContent = "0.000";
        this.els.operatorHeightError.textContent = "0.000";
        this.els.operatorSequenceState.textContent = "0/0";
        this.els.operatorMaxOrderError.textContent = "0.000";
        return;
      }
      const commutes = diagnostic.modeCommutes && diagnostic.commandCommutes && diagnostic.lockAlphaCommutes && diagnostic.lockZCommutes && diagnostic.finalAlphaError < 1e-9 && diagnostic.finalHeightError < 1e-9;
      const sequence = diagnostic.sequence;
      this.els.operatorOrderState.textContent = commutes ? "commutes" : "path dependent";
      this.els.operatorAlphaError.textContent = diagnostic.finalAlphaError.toFixed(3);
      this.els.operatorHeightError.textContent = diagnostic.finalHeightError.toFixed(3);
      this.els.operatorOrderState.title = `mode ${diagnostic.modeCommutes ? "same" : "diff"}, commands ${diagnostic.commandCommutes ? "same" : "diff"}, lock alpha ${diagnostic.lockAlphaCommutes ? "same" : "diff"}, lock z ${diagnostic.lockZCommutes ? "same" : "diff"}`;
      this.els.operatorSequenceState.textContent = sequence
        ? `${sequence.noncommutingAdjacentPairs}/${sequence.adjacentPairCount}`
        : "0/0";
      this.els.operatorMaxOrderError.textContent = Number(sequence?.maxOrderError || 0).toFixed(3);
      this.els.operatorSequenceState.title = sequence
        ? `reverse da ${sequence.reverseAlphaError.toFixed(3)}, dz ${sequence.reverseHeightError.toFixed(3)}, mode changes ${sequence.reverseModeChanges}`
        : "";
    }

    applyCellVisualMode(mode) {
      const visualMode = ["abstract", "graph", "sheetOnly", "paperRad", "calibratedRad", "cadRad", "mechanism"].includes(mode) ? mode : "abstract";
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
      } else if (visualMode === "calibratedRad") {
        this.state.view.explodedSelected = false;
        Object.assign(this.state.view, {
          gapsVisible: true,
          stopsVisible: true,
          pivotsVisible: true,
          linkagesVisible: true,
          fastenersVisible: false,
          actuatorsVisible: true,
          measurementsVisible: true,
          measurementMode: "hardware",
        });
      } else if (visualMode === "cadRad") {
        this.state.view.explodedSelected = false;
        Object.assign(this.state.view, {
          gapsVisible: true,
          stopsVisible: true,
          pivotsVisible: true,
          linkagesVisible: true,
          fastenersVisible: true,
          actuatorsVisible: true,
          measurementsVisible: true,
          measurementMode: "hardware",
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
      } else if (visualMode === "graph") {
        this.state.view.explodedSelected = false;
        Object.assign(this.state.view, {
          membraneVisible: false,
          referenceVisible: false,
          displacementVectorsVisible: false,
          gapsVisible: false,
          stopsVisible: false,
          pivotsVisible: false,
          linkagesVisible: false,
          fastenersVisible: false,
          actuatorsVisible: false,
          measurementsVisible: false,
          measurementLabelsVisible: false,
          targetVisible: false,
          surfaceContoursVisible: false,
          overlayMode: "state",
        });
      } else if (visualMode === "sheetOnly") {
        this.state.view.explodedSelected = false;
        Object.assign(this.state.view, {
          membraneVisible: true,
          referenceVisible: false,
          displacementVectorsVisible: false,
          gapsVisible: false,
          stopsVisible: false,
          pivotsVisible: false,
          linkagesVisible: false,
          fastenersVisible: false,
          actuatorsVisible: false,
          measurementsVisible: false,
          measurementLabelsVisible: false,
          targetVisible: this.state.target.type !== "none",
          surfaceContoursVisible: true,
          contourMode: this.state.target.type === "none" ? "membrane" : "error",
          overlayMode: this.state.inverse?.jacobian?.targetReachability ? "topologyBlocked" : "error",
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
      this.els.paperSideLengthMm.value = s.grid.paperSideLengthMm ?? 35;
      this.els.paperHoleToleranceMm.value = s.grid.paperHoleToleranceMm ?? 0.1;
      const profile = typeof RAD.hardwareProfile === "function" ? RAD.hardwareProfile(s) : this.ensureHardwareProfile();
      this.els.hardwareProfileName.value = profile.name || "paper-reference";
      this.els.hardwareProfileSource.value = profile.source || "unspecified";
      for (const [field, element] of this.hardwareProfileNumberInputs()) {
        element.value = this.profileNumberInputValue(profile[field]);
      }
      if (typeof RAD.ensureBoundarySchema === "function") RAD.ensureBoundarySchema(s);
      const boundary = s.boundary || {};
      this.els.boundaryMode.value = boundary.mode || "free";
      this.els.boundaryWallType.value = boundary.wallType || "rigid";
      this.els.boundaryXMin.value = this.boundaryNumberInputValue(boundary.xMin);
      this.els.boundaryXMax.value = this.boundaryNumberInputValue(boundary.xMax);
      this.els.boundaryYMin.value = this.boundaryNumberInputValue(boundary.yMin);
      this.els.boundaryYMax.value = this.boundaryNumberInputValue(boundary.yMax);
      this.els.boundaryChannelAxis.value = boundary.channelAxis || "x";
      this.els.boundaryChannelWidth.value = Number(boundary.channelWidth || 0);
      this.els.boundaryZGain.value = Number(boundary.zGain ?? 1.15);
      this.els.boundaryZThreshold.value = Number(boundary.zThreshold ?? 0.02);
      this.els.boundaryZPower.value = Number(boundary.zPower ?? 1.5);
      this.els.paintRadius.value = Math.max(0, Math.min(2, Number(s.view.paintRadius || 0)));
      this.els.cellVisualMode.value = s.view.cellVisualMode || "abstract";
      this.syncViewportTitle(s.view.cellVisualMode || "abstract");
      this.els.simulationMode.value = s.view.simulationMode || "kinematic";
      this.els.empiricalLockModel.checked = s.grid.empiricalLockModel !== false;
      this.els.lockStateIndex.value = Math.max(1, Math.min(3, Number(s.grid.lockStateIndex || 1)));
      this.renderTwoCellBench();
      this.els.overlayMode.value = s.view.overlayMode;
      this.els.showMembrane.checked = s.view.membraneVisible;
      this.els.showReference.checked = s.view.referenceVisible !== false;
      this.els.showDisplacementVectors.checked = s.view.displacementVectorsVisible !== false;
      this.els.showExternalCase.checked = s.view.externalCaseVisible !== false;
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
      this.els.characterizationScope.value = s.experiment.characterizationScope || "single";
      const { r, c } = s.selection;
      const limits = RAD.commandLimits(s);
      this.els.alphaCommand.min = -limits.alphaContract;
      this.els.alphaCommand.max = limits.alphaExpand;
      this.els.zCommand.min = -limits.z;
      this.els.zCommand.max = limits.z;
      this.els.alphaCommand.value = s.cells.commandAlpha[r][c] ?? 0;
      this.els.zCommand.value = s.cells.commandZ[r][c] ?? 0;
      this.els.locked.checked = s.cells.locked[r][c];
      this.els.positionLocked.checked = s.cells.positionLocked?.[r]?.[c] === true;
      this.els.removed.checked = s.cells.removed?.[r]?.[c] === true;
      this.els.actuatorAllowed.checked = s.cells.removed?.[r]?.[c] === true ? false : s.cells.actuatorAllowed?.[r]?.[c] !== false;
      this.updateLabels(null);
      this.syncQuickDock();
      this.syncPaintMode();
    }

    syncViewportTitle(mode) {
      if (!this.els.viewportTitle) return;
      const titles = {
        abstract: "3D Cell Abstraction",
        graph: "Topology Graph",
        sheetOnly: "Sheet Contour",
        paperRad: "Paper RAD Cell",
        calibratedRad: "Calibrated RAD Cell",
        cadRad: "CAD RAD Cell",
        mechanism: "Mechanism Detail",
      };
      this.els.viewportTitle.textContent = titles[mode] || titles.abstract;
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

    renderModelProvenance() {
      if (typeof RAD.modelProvenance !== "function" || !this.els.modelProvenanceList) return;
      const summary = typeof RAD.provenanceSummary === "function" ? RAD.provenanceSummary() : {};
      this.els.provenancePaper.textContent = `paper ${summary["paper-supported"] || 0}`;
      this.els.provenanceAssumptions.textContent = `assumptions ${summary["implementation-assumption"] || 0}`;
      this.els.provenanceDiagnostics.textContent = `diagnostics ${summary["simulator-diagnostic"] || 0}`;
      this.els.provenanceGaps.textContent = `gaps ${summary["calibration-gap"] || 0}`;
      this.els.modelProvenanceList.replaceChildren();
      for (const item of RAD.modelProvenance()) {
        const row = document.createElement("div");
        row.className = `provenance-item provenance-${item.status}`;
        row.title = `${item.source}: ${item.evidence} Limitation: ${item.limitation}`;
        const label = document.createElement("span");
        label.textContent = item.label;
        const status = document.createElement("small");
        status.textContent = item.status;
        row.append(label, status);
        this.els.modelProvenanceList.append(row);
      }
    }

    updateLabels(sim) {
      const s = this.state;
      const { r, c } = s.selection;
      document.getElementById("rowsOut").textContent = s.grid.rows;
      document.getElementById("colsOut").textContent = s.grid.cols;
      document.getElementById("backlashOut").textContent = Number(s.grid.backlash).toFixed(2);
      document.getElementById("couplingOut").textContent = Number(s.grid.couplingGain).toFixed(2);
      document.getElementById("zCouplingOut").textContent = Number(s.grid.zCouplingGain ?? 0.32).toFixed(2);
      document.getElementById("lockStateIndexOut").textContent = String(Math.max(1, Math.min(3, Number(s.grid.lockStateIndex || 1))));
      if (typeof RAD.lockDatasetSummary === "function") {
        const lockSummary = RAD.lockDatasetSummary();
        const applied = sim?.metrics?.empiricalLockDatasetApplied ? `, applied ${sim.metrics.empiricalLockNearestFile}` : "";
        document.getElementById("lockDatasetOut").textContent = lockSummary.available
          ? `${lockSummary.configurationCount} configs, ${lockSummary.outputDimension} xyz outputs${applied}`
          : "not loaded";
      }
      document.getElementById("pinRadiusOut").textContent = Number(s.grid.pinRadius ?? 0.18).toFixed(3);
      document.getElementById("holeRadiusOut").textContent = Number(s.grid.holeRadius ?? 0.225).toFixed(3);
      const calibration = RAD.paperRadCalibration(s);
      document.getElementById("pinHoleClearanceOut").textContent = RAD.pinHoleClearance(s).toFixed(3);
      document.getElementById("paperSideLengthMmOut").textContent = calibration.sideLengthMm.toFixed(1);
      document.getElementById("paperHoleToleranceMmOut").textContent = calibration.fabricationHoleToleranceMm.toFixed(3);
      document.getElementById("backlashMmOut").textContent = calibration.configuredBacklashMm.toFixed(3);
      document.getElementById("pinHoleClearanceMmOut").textContent = calibration.pinHoleClearanceMm.toFixed(3);
      document.getElementById("holeToleranceModelOut").textContent = calibration.fabricationHoleToleranceModel.toFixed(4);
      if (typeof RAD.calibrationProfileSummary === "function") {
        const profileSummary = RAD.calibrationProfileSummary(s);
        document.getElementById("hardwareProfileOut").textContent = profileSummary.profile.name;
        document.getElementById("hardwareCoverageOut").textContent = `${profileSummary.measuredCount}/${profileSummary.totalCount}`;
        const labels = {
          pinRadiusMm: "pin",
          holeRadiusMm: "hole",
          plateThicknessMm: "plate",
          jointStackHeightMm: "stack",
          bossRadiusMm: "boss",
        };
        document.getElementById("hardwareMissingOut").textContent =
          profileSummary.missingFields.map((field) => labels[field] || field).join(", ") || "none";
        if (typeof RAD.calibrationReadiness === "function") {
          const readiness = RAD.calibrationReadiness(s);
          document.getElementById("hardwareReadinessOut").textContent = readiness.level;
          document.getElementById("hardwareSolverGapOut").textContent =
            readiness.solverGaps.map((gap) => gap.split(" ")[0]).join(", ");
        }
        if (typeof RAD.calibrationMeasurementPlan === "function") {
          const missingTasks = RAD.calibrationMeasurementPlan(s).filter((task) => task.status === "missing");
          document.getElementById("hardwareMeasurementPlanOut").textContent =
            missingTasks.slice(0, 3).map((task) => task.label.replace(/^Measure /, "")).join(", ") || "none";
        }
      }
      document.getElementById("alphaCommandOut").textContent = Number(this.els.alphaCommand.value).toFixed(2);
      document.getElementById("zCommandOut").textContent = Number(this.els.zCommand.value).toFixed(2);
      document.getElementById("boundaryChannelWidthOut").textContent = Number(s.boundary?.channelWidth || 0).toFixed(2);
      document.getElementById("boundaryZGainOut").textContent = Number(s.boundary?.zGain ?? 1.15).toFixed(2);
      document.getElementById("boundaryZThresholdOut").textContent = Number(s.boundary?.zThreshold ?? 0.02).toFixed(3);
      document.getElementById("boundaryZPowerOut").textContent = Number(s.boundary?.zPower ?? 1.5).toFixed(2);
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
      this.renderCharacterization();
      this.renderSequenceChart();
      document.getElementById("selectedCell").textContent = `row ${r}, col ${c}`;
      this.els.quickDockMini.textContent = `r${r} c${c} / a ${Number(this.els.alphaCommand.value).toFixed(2)} / z ${Number(this.els.zCommand.value).toFixed(2)}`;
      this.updateCouplingInspector(r, c);
      if (sim) {
        const residual = sim.target[r][c] - sim.height[r][c];
        const localStrain = this.localLinkStrain(sim, r, c);
        const zResidual = sim.zResidual?.[r]?.[c] || 0;
        const compressionResidual = sim.compressionResidual?.[r]?.[c] || 0;
        const inducedHeight = sim.inducedHeight?.[r]?.[c] || 0;
        const topology = s.cells.removed?.[r]?.[c] ? "removed, " : "";
        document.getElementById("selectedValues").textContent = `${topology}alpha ${sim.alpha[r][c].toFixed(3)}, theta ${sim.theta[r][c].toFixed(1)}, z ${sim.height[r][c].toFixed(3)}, z residual ${zResidual.toFixed(3)}, compression ${compressionResidual.toFixed(3)}, induced z ${inducedHeight.toFixed(3)}, slope ${sim.slope.magnitude[r][c].toFixed(3)}, normal tilt ${sim.slope.tilt[r][c].toFixed(1)}, residual ${residual.toFixed(3)}, link strain ${localStrain.toFixed(3)}, ref disp ${sim.displacement[r][c].toFixed(3)}`;
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
        document.getElementById("maxCompressionResidualOut").textContent = Number(sim.metrics.maxCompressionResidual || 0).toFixed(3);
        document.getElementById("maxInducedHeightOut").textContent = Number(sim.metrics.maxInducedHeight || 0).toFixed(3);
        document.getElementById("constrainedCellsOut").textContent = String(sim.metrics.constrainedCellCount || 0);
        document.getElementById("constraintSolverStatusOut").textContent = sim.metrics.constraintSolved ? "constraint solved" : sim.metrics.physicalPreview ? "spring preview" : "fast preview";
        document.getElementById("constraintMaxEdgeErrorOut").textContent = Number(sim.metrics.constraintMaxEdgeError || 0).toFixed(3);
        document.getElementById("constraintMeanEdgeErrorOut").textContent = Number(sim.metrics.constraintMeanEdgeError || 0).toFixed(3);
        document.getElementById("constraintContactCountOut").textContent = String(sim.metrics.constraintContactCount || 0);
        document.getElementById("sensitivityCells").textContent = String(s.inverse?.sensitivity?.controllableCells || 0);
        document.getElementById("meanSensitivity").textContent = Number(s.inverse?.sensitivity?.meanGain || 0).toFixed(3);
        document.getElementById("jacobianColumns").textContent = String(s.inverse?.jacobian?.columnCount || 0);
        document.getElementById("meanReachability").textContent = Number(s.inverse?.jacobian?.meanCoverage || 0).toFixed(3);
        document.getElementById("underTargetCells").textContent = String(s.inverse?.jacobian?.targetReachability?.underactuatedHeightCells || 0);
        document.getElementById("linearFitSteps").textContent = String(s.inverse?.linearSolution?.steps || 0);
        document.getElementById("linearFitError").textContent = Number(s.inverse?.linearSolution?.projectedError || 0).toFixed(3);
      }
      this.renderInversePlan();
      this.updateOperatorInspector(this.lastOperatorDiagnostic);
      this.renderEvents();
    }

    runCharacterization() {
      const scope = this.els.characterizationScope.value || "single";
      this.state.experiment.characterizationScope = scope;
      const result = RAD.characterizeLocalResponse(this.state, { scope, ...this.state.selection });
      this.state.experiment.characterization = result;
      const framework = this.currentFrameworkReportMetadata(scope);
      this.state.experiment.frameworkLawCandidates = framework.operatorLawCandidates;
      this.state.experiment.frameworkFormalizationTargets = framework.formalizationTargets;
      this.state.view.overlayMode = "operatorInteraction";
      RAD.recordEvent(this.state, {
        type: "characterization",
        scope: result.scope,
        cells: result.regionCellCount,
        activeSources: result.activeSources,
        responseCells: result.responseCells,
        superpositionError: result.superpositionError,
        pairwisePairs: result.pairwiseEvaluatedPairs,
        pairwiseNonadditive: result.pairwiseNonadditivePairs,
        pairwiseMaxError: result.pairwiseMaxInteractionError,
        pairwiseMaxDegree: result.pairwiseInteractionDegreeMax,
        pairwiseDensity: result.pairwiseInteractionDensity,
      });
      this.syncControls();
      this.onChange(this.state);
    }

    renderCharacterization() {
      const result = this.state.experiment.characterization;
      if (!result) {
        document.getElementById("characterizationSummary").textContent = "No response run";
        document.getElementById("characterizationDetail").textContent = "reach a0 z0";
        document.getElementById("characterizationZSign").textContent = "z sign +0/-0";
        document.getElementById("characterizationZExtrema").textContent = "z max +0.000/-0.000";
        document.getElementById("characterizationSuperposition").textContent = "superposition 0.000";
        document.getElementById("characterizationScale").textContent = "clearance 0.000 mm";
        document.getElementById("characterizationPairwise").textContent = "pairs 0/0 nonadd 0";
        document.getElementById("characterizationPairwiseMax").textContent = "pair max 0.000";
        document.getElementById("characterizationHotspot").textContent = "hotspot none";
        document.getElementById("characterizationRank").textContent = "rank a0 z0";
        document.getElementById("characterizationUnderactuated").textContent = "under a0 z0";
        document.getElementById("characterizationPhysical").textContent = "phys rms 0.000";
        document.getElementById("characterizationPhysicalMax").textContent = "phys max 0.000";
        document.getElementById("characterizationDecay").textContent = "decay a0.00 z0.00";
        document.getElementById("characterizationDecayLength").textContent = "len a0.0 z0.0";
        this.renderFrameworkLawCandidates();
        this.renderFormalizationTargets();
        this.renderCalibrationResults();
        this.renderResponseAtlasSweep();
        return;
      }
      const superposition = result.superpositionSkipped
        ? `superposition skipped (${result.superpositionSources} sources)`
        : `superposition ${Number(result.superpositionError || 0).toFixed(3)}`;
      document.getElementById("characterizationSummary").textContent = `${result.scope}: ${result.activeSources}/${result.regionCellCount} active, ${result.responseCells} cells`;
      document.getElementById("characterizationDetail").textContent = `reach a${result.alphaReachCells} z${result.zReachCells}, die ${result.alphaDieOff}/${result.zDieOff}`;
      document.getElementById("characterizationZSign").textContent = `z sign +${result.positiveZReachCells || 0}/-${result.negativeZReachCells || 0}`;
      document.getElementById("characterizationZExtrema").textContent = `z max +${Number(result.maxPositiveHeightDelta || 0).toFixed(3)}/-${Math.abs(Number(result.maxNegativeHeightDelta || 0)).toFixed(3)}`;
      document.getElementById("characterizationSuperposition").textContent = superposition;
      document.getElementById("characterizationScale").textContent = `clearance ${Number(result.pinHoleClearanceMm || 0).toFixed(3)} mm`;
      document.getElementById("characterizationPairwise").textContent = `pairs ${result.pairwiseEvaluatedPairs || 0}/${result.pairwiseTotalPairs || 0} nonadd ${result.pairwiseNonadditivePairs || 0}${result.pairwiseTruncated ? " trunc" : ""}`;
      document.getElementById("characterizationPairwiseMax").textContent = `pair max ${Number(result.pairwiseMaxInteractionError || 0).toFixed(3)} deg ${result.pairwiseInteractionDegreeMax || 0}`;
      const hotspot = this.strongestInteractionHotspot(result);
      document.getElementById("characterizationHotspot").textContent = hotspot ? `hotspot r${hotspot.r} c${hotspot.c} ${hotspot.value.toFixed(3)} d${hotspot.degree}` : "hotspot none";
      document.getElementById("characterizationRank").textContent = `rank a${result.responseRankAlpha || 0} z${result.responseRankHeight || 0}`;
      document.getElementById("characterizationUnderactuated").textContent = `under a${result.alphaUnderactuatedCells || 0} z${result.heightUnderactuatedCells || 0}`;
      const physicalLabel = result.physicalPreviewAvailable
        ? `phys rms ${Number(result.physicalHeightRmsError || 0).toFixed(3)}`
        : "phys unavailable";
      const physicalMaxLabel = result.physicalPreviewAvailable
        ? `phys max ${Number(result.physicalHeightMaxError || 0).toFixed(3)}`
        : "phys max --";
      document.getElementById("characterizationPhysical").textContent = physicalLabel;
      document.getElementById("characterizationPhysicalMax").textContent = physicalMaxLabel;
      document.getElementById("characterizationDecay").textContent = `decay a${Number(result.alphaDecayRatio || 0).toFixed(2)} z${Number(result.zDecayRatio || 0).toFixed(2)}`;
      document.getElementById("characterizationDecayLength").textContent = `len a${Number(result.alphaDecayLength || 0).toFixed(1)} z${Number(result.zDecayLength || 0).toFixed(1)}`;
      this.renderFrameworkLawCandidates();
      this.renderFormalizationTargets();
      this.renderCalibrationResults();
      this.renderResponseAtlasSweep();
    }

    currentFrameworkReportMetadata(scope) {
      const empty = { operatorLawCandidates: null, formalizationTargets: null };
      if (typeof RAD.programmableDiscontinuityReport !== "function") return empty;
      try {
        const report = RAD.programmableDiscontinuityReport(this.state, {
          scope,
          ...this.state.selection,
          includeResponseMatrix: false,
          includeFields: false,
        });
        return {
          operatorLawCandidates: report.operatorLawCandidates || null,
          formalizationTargets: report.formalizationTargets || null,
        };
      } catch {
        return empty;
      }
    }

    formatFrameworkLawEvidence(law) {
      if (!law) return "";
      const evidence = law.evidence || {};
      if (law.id === "composition_nonadditivity") {
        return ` err ${Number(evidence.maxSuperpositionError || 0).toFixed(3)}`;
      }
      if (law.id === "event_order_noncommutativity") {
        return ` err ${Number(evidence.maxOrderError || 0).toFixed(3)}`;
      }
      if (law.id === "rank_limited_reachability") {
        return ` rank a${evidence.alphaRank || 0} z${evidence.heightRank || 0}`;
      }
      if (law.id === "bounded_locality") {
        return ` die a${evidence.alphaLocalityRadius || 0} z${evidence.zLocalityRadius || 0}`;
      }
      return "";
    }

    renderFrameworkLawCandidates() {
      const candidates = this.state.experiment.frameworkLawCandidates;
      if (!candidates || !Array.isArray(candidates.laws)) {
        document.getElementById("frameworkLawSummary").textContent = "laws not run";
        document.getElementById("frameworkLawDetail").textContent = "candidate --";
        return;
      }
      const laws = candidates.laws;
      const supported = laws.filter((law) => law.supportedByDiagnostic);
      const preferred = [
        "composition_nonadditivity",
        "event_order_noncommutativity",
        "rank_limited_reachability",
        "bounded_locality",
      ];
      const primary =
        preferred.map((id) => supported.find((law) => law.id === id)).find(Boolean) ||
        supported[0] ||
        laws[0];
      const status = primary?.supportedByDiagnostic ? "supported" : "not supported";
      const label = primary ? primary.id.replace(/_/g, " ") : "none";
      document.getElementById("frameworkLawSummary").textContent = `laws ${supported.length}/${laws.length} supported`;
      document.getElementById("frameworkLawDetail").textContent = `candidate ${label} ${status}${this.formatFrameworkLawEvidence(primary)}`;
    }

    renderFormalizationTargets() {
      const manifest = this.state.experiment.frameworkFormalizationTargets;
      if (!manifest || !Array.isArray(manifest.targets)) {
        document.getElementById("formalizationSummary").textContent = "formal targets not run";
        document.getElementById("formalizationDetail").textContent = "proof --";
        return;
      }
      const targets = manifest.targets;
      const ready = targets.filter((target) => target.readyForLean);
      const tooling = manifest.tooling || {};
      const status = tooling.status || "unknown";
      const preferred = [
        "dead_zone_zero_inside_backlash",
        "dead_zone_piecewise_linear_outside_gap",
        "lock_projection_idempotent",
        "noncommutativity_witness_from_order_error",
        "bounded_locality_witness",
      ];
      const primary =
        preferred.map((id) => targets.find((target) => target.id === id)).find(Boolean) ||
        targets[0];
      const label = primary ? primary.id.replace(/_/g, " ") : "none";
      const targetStatus = primary?.status || "unknown";
      document.getElementById("formalizationSummary").textContent = `formal ${ready.length}/${targets.length} ready, Lean ${status}`;
      document.getElementById("formalizationDetail").textContent = `proof ${label} ${targetStatus}`;
    }

    renderCalibrationResults() {
      const summary = this.state.experiment.calibrationComparisonSummary;
      if (!summary) {
        document.getElementById("calibrationResultsSummary").textContent = "calibration results none";
        document.getElementById("calibrationResultsError").textContent = "height rmse --";
        document.getElementById("modelProfileSummary").textContent = "model profile unchecked";
        document.getElementById("modelProfileAudit").textContent = "profile audit --";
        document.getElementById("modelProfileHoldoutSummary").textContent = "holdout not loaded";
        document.getElementById("modelProfileHoldoutDetail").textContent = "holdout validation --";
        return;
      }
      const measured = Number(summary.measuredCellCount || 0);
      const missing = Number(summary.missingObservationCount || 0);
      const stepCount = Number(summary.stepCount || 0);
      const height = summary.heightRmseMean === null || summary.heightRmseMean === undefined ? "--" : Number(summary.heightRmseMean).toFixed(4);
      const alpha = summary.alphaRmseMean === null || summary.alphaRmseMean === undefined ? "--" : Number(summary.alphaRmseMean).toFixed(4);
      const heightBias = summary.meanSignedHeightError === null || summary.meanSignedHeightError === undefined ? "--" : Number(summary.meanSignedHeightError).toFixed(4);
      const alphaBias = summary.meanSignedAlphaError === null || summary.meanSignedAlphaError === undefined ? "--" : Number(summary.meanSignedAlphaError).toFixed(4);
      const max = summary.maxCombinedError === null || summary.maxCombinedError === undefined ? "--" : Number(summary.maxCombinedError).toFixed(4);
      const heightFit = summary.fit?.height;
      const alphaFit = summary.fit?.alpha;
      const heightGain = heightFit?.suggestedGain === null || heightFit?.suggestedGain === undefined ? "--" : Number(heightFit.suggestedGain).toFixed(3);
      const alphaGain = alphaFit?.suggestedGain === null || alphaFit?.suggestedGain === undefined ? "--" : Number(alphaFit.suggestedGain).toFixed(3);
      const parameterEstimates =
        this.state.experiment.calibrationComparison && typeof RAD.calibrationParameterEstimates === "function"
          ? RAD.calibrationParameterEstimates(
              this.state,
              this.state.experiment.calibrationResults || {},
              this.state.experiment.calibrationComparison,
              summary
            )
          : null;
      const zEstimate = parameterEstimates?.estimates?.zCouplingGain?.estimate;
      const zCoupling = zEstimate === null || zEstimate === undefined ? "--" : Number(zEstimate).toFixed(3);
      const slipEstimate = parameterEstimates?.estimates?.verticalFreePlay?.meanMeasuredSlipMm;
      const slip = slipEstimate === null || slipEstimate === undefined ? "--" : `${Number(slipEstimate).toFixed(3)}mm`;
      const contactStatus = parameterEstimates?.estimates?.contactProxy?.status || "contact --";
      const residualMax = summary.fitResidualMaxCombinedError === null || summary.fitResidualMaxCombinedError === undefined ? "--" : Number(summary.fitResidualMaxCombinedError).toFixed(4);
      const cell = summary.worstCell ? ` cell r${summary.worstCell.row} c${summary.worstCell.col}` : "";
      const worst = summary.worstStepId ? ` worst ${summary.worstStepId}` : "";
      document.getElementById("calibrationResultsSummary").textContent = `cal results ${stepCount} steps, ${measured} cells, missing ${missing}, fit h ${heightGain} a ${alphaGain}, z ${zCoupling}`;
      document.getElementById("calibrationResultsError").textContent = `h ${height}/${heightBias}, a ${alpha}/${alphaBias}, raw ${max}, resid ${residualMax}, slip ${slip}, ${contactStatus}${cell}${worst}`;
      const profile = this.state.experiment.calibrationModelProfile;
      const comparison = this.state.experiment.calibrationModelProfileComparison;
      const selection = this.state.experiment.calibrationModelProfileSelection;
      const application = this.state.experiment.calibrationModelProfileApplication;
      const historyCount = Array.isArray(this.state.experiment.calibrationModelProfileHistory)
        ? this.state.experiment.calibrationModelProfileHistory.length
        : 0;
      if (comparison) {
        const beforeScore = comparison.before?.metrics?.residualScore;
        const afterScore = comparison.after?.metrics?.residualScore;
        const deltaScore = comparison.delta?.residualScore;
        const beforeLabel = beforeScore === null || beforeScore === undefined ? "--" : Number(beforeScore).toFixed(4);
        const afterLabel = afterScore === null || afterScore === undefined ? "--" : Number(afterScore).toFixed(4);
        const deltaLabel = deltaScore === null || deltaScore === undefined ? "--" : Number(deltaScore).toFixed(4);
        const selectedLabel =
          selection?.selectedIndex === null || selection?.selectedIndex === undefined
            ? "no selected profile"
            : `selected ${selection.selectedIndex}/${selection.candidateCount}`;
        document.getElementById("modelProfileSummary").textContent = `profile ${comparison.improved ? "improved" : "not improved"} ${beforeLabel} -> ${afterLabel}, d ${deltaLabel}, ${selectedLabel}`;
      } else {
        const safeCount = Number(profile?.safeUpdateCount || 0);
        const diagnosticCount = Number(profile?.diagnosticOnlyCount || 0);
        document.getElementById("modelProfileSummary").textContent = `profile unchecked, safe ${safeCount}, diagnostic ${diagnosticCount}`;
      }
      if (application) {
        const applied = Array.isArray(application.appliedUpdates) ? application.appliedUpdates.length : 0;
        const skipped = Array.isArray(application.skippedUpdates) ? application.skippedUpdates.length : 0;
        const z = application.configAfter?.zCouplingGain;
        const zLabel = z === null || z === undefined ? "--" : Number(z).toFixed(3);
        document.getElementById("modelProfileAudit").textContent = `audit applied ${applied}, skipped ${skipped}, z ${zLabel}, history ${historyCount}`;
      } else {
        document.getElementById("modelProfileAudit").textContent = `profile audit none, history ${historyCount}`;
      }
      const holdoutValidation = this.state.experiment.calibrationModelProfileHoldoutValidation;
      const executionValidation = this.state.experiment.calibrationBenchExecutionValidation;
      const holdoutSummary = this.state.experiment.calibrationHoldoutComparisonSummary;
      if (holdoutValidation) {
        const beforeScore = holdoutValidation.holdout?.before?.metrics?.residualScore;
        const afterScore = holdoutValidation.holdout?.after?.metrics?.residualScore;
        const deltaScore = holdoutValidation.holdout?.delta?.residualScore;
        const beforeLabel = beforeScore === null || beforeScore === undefined ? "--" : Number(beforeScore).toFixed(4);
        const afterLabel = afterScore === null || afterScore === undefined ? "--" : Number(afterScore).toFixed(4);
        const deltaLabel = deltaScore === null || deltaScore === undefined ? "--" : Number(deltaScore).toFixed(4);
        document.getElementById("modelProfileHoldoutSummary").textContent =
          `holdout ${holdoutValidation.holdoutPass ? "passed" : "not passed"} ${beforeLabel} -> ${afterLabel}, d ${deltaLabel}`;
        const missingEvidence = holdoutValidation.splitMetadata?.missingEvidence || [];
        const missingLabel = Array.isArray(missingEvidence) && missingEvidence.length
          ? `missing ${missingEvidence.slice(0, 2).join("/")}`
          : holdoutValidation.splitMetadata?.independenceStatus || "split ok";
        const executionLabel = executionValidation?.summary
          ? `, exec ${executionValidation.summary.executionValidationPass ? "ready" : executionValidation.summary.status || "pending"}`
          : "";
        document.getElementById("modelProfileHoldoutDetail").textContent =
          `fit ${holdoutValidation.fitPass ? "pass" : "fail"}, residual ${holdoutValidation.residualValidationPass ? "pass" : "fail"}, independent ${holdoutValidation.independentValidationPass ? "pass" : "pending"}, ${missingLabel}${executionLabel}`;
      } else if (holdoutSummary) {
        const holdoutSteps = Number(holdoutSummary.stepCount || 0);
        const holdoutMissing = Number(holdoutSummary.missingObservationCount || 0);
        const holdoutResidual =
          holdoutSummary.fitResidualMaxCombinedError === null || holdoutSummary.fitResidualMaxCombinedError === undefined
            ? "--"
            : Number(holdoutSummary.fitResidualMaxCombinedError).toFixed(4);
        document.getElementById("modelProfileHoldoutSummary").textContent =
          `holdout loaded ${holdoutSteps} steps, missing ${holdoutMissing}, resid ${holdoutResidual}`;
        document.getElementById("modelProfileHoldoutDetail").textContent = "run Check Holdout to replay profile";
      } else {
        document.getElementById("modelProfileHoldoutSummary").textContent = "holdout not loaded";
        document.getElementById("modelProfileHoldoutDetail").textContent = "holdout validation --";
      }
    }

    runResponseAtlasSweep() {
      if (typeof RAD.responseAtlasSweep !== "function") return;
      const sweep = RAD.responseAtlasSweep(this.state, { ...this.state.selection });
      this.state.experiment.responseAtlasSweep = sweep;
      RAD.recordEvent(this.state, {
        type: "response-atlas-sweep",
        samples: sweep.summary.sampleCount,
        maxNeighborResidual: sweep.summary.maxObservedNeighborZResidual,
        maxSuperpositionError: sweep.summary.maxSuperpositionError,
      });
      this.renderResponseAtlasSweep();
      this.syncControls();
      this.onChange(this.state);
    }

    renderResponseAtlasSweep() {
      const sweep = this.state.experiment.responseAtlasSweep;
      if (!sweep) {
        document.getElementById("sweepSummary").textContent = "sweep not run";
        document.getElementById("sweepTrend").textContent = "neighbor z --";
        return;
      }
      const summary = sweep.summary || {};
      const trend = sweep.trends?.byPinHoleClearance || [];
      const low = trend[0];
      const high = trend[trend.length - 1];
      const lowResidual = low ? Number(low.maxObservedNeighborZResidual || 0) : 0;
      const highResidual = high ? Number(high.maxObservedNeighborZResidual || 0) : 0;
      const drop = lowResidual > 1e-12 ? 100 * (1 - highResidual / lowResidual) : 0;
      const dominant = sweep.sensitivity?.dominant;
      const dominantLabel = dominant
        ? `, sens ${dominant.parameter}/${dominant.metric} ${Number(dominant.slope || 0).toFixed(2)}`
        : "";
      const residualLaw = (sweep.operatorLawCandidates?.laws || []).find(
        (law) => law.parameter === "pinHoleClearance" && law.metric === "maxObservedNeighborZResidual"
      );
      const lawLabel = residualLaw ? `, law ${residualLaw.parameter}/${residualLaw.monotonicity}` : "";
      document.getElementById("sweepSummary").textContent = `sweep ${summary.sampleCount || 0} samples, reach a${summary.maxAlphaReach || 0} z${summary.maxZReach || 0}`;
      document.getElementById("sweepTrend").textContent = `neighbor z ${Number(summary.maxObservedNeighborZResidual || 0).toFixed(3)}, clear drop ${drop.toFixed(0)}%, super ${Number(summary.maxSuperpositionError || 0).toFixed(3)}${dominantLabel}${lawLabel}`;
    }

    strongestInteractionHotspot(result) {
      const map = result?.pairwiseInteractionMap;
      if (!Array.isArray(map)) return null;
      let best = null;
      for (let r = 0; r < map.length; r += 1) {
        const row = map[r] || [];
        for (let c = 0; c < row.length; c += 1) {
          const value = Math.abs(Number(row[c]) || 0);
          const degree = Number(result?.pairwiseInteractionDegreeMap?.[r]?.[c]) || 0;
          if (!best || value > best.value) best = { r, c, value, degree };
        }
      }
      return best && best.value > 1e-12 ? best : null;
    }

    selectInteractionHotspot() {
      const hotspot = this.strongestInteractionHotspot(this.state.experiment.characterization);
      if (!hotspot) return;
      this.state.selection = { r: hotspot.r, c: hotspot.c };
      this.state.view.overlayMode = "operatorInteraction";
      this.syncControls();
      this.onChange(this.state);
    }

    calibrationHotspotRanking(residualMode) {
      const summary = this.state.experiment.calibrationComparisonSummary;
      const comparison = this.state.experiment.calibrationComparison;
      const ranked = residualMode
        ? summary?.fitResidualTopCells || comparison?.fitResidualField?.topCells
        : summary?.topCells || comparison?.field?.topCells;
      if (Array.isArray(ranked) && ranked.length) return ranked;
      const fallback =
        (residualMode
          ? summary?.fitResidualWorstCell || comparison?.fitResidualField?.worstCell
          : summary?.worstCell || comparison?.field?.worstCell) ||
        summary?.worstCell ||
        comparison?.field?.worstCell;
      return fallback ? [fallback] : [];
    }

    selectCalibrationHotspot(step = 0) {
      const residualMode = this.state.view.overlayMode === "calibrationResidual";
      const ranking = this.calibrationHotspotRanking(residualMode);
      if (!ranking.length) return;
      let index = 0;
      if (step !== 0 && this.state.selection) {
        const currentIndex = ranking.findIndex((cell) => cell.row === this.state.selection.r && cell.col === this.state.selection.c);
        index = currentIndex >= 0 ? (currentIndex + step + ranking.length) % ranking.length : 0;
      }
      const cell = ranking[index] || ranking[0];
      if (!cell) return;
      this.state.selection = { r: cell.row, c: cell.col };
      this.state.view.overlayMode = residualMode ? "calibrationResidual" : "calibrationError";
      this.syncControls();
      this.onChange(this.state);
    }

    selectUnderactuatedTarget() {
      if (!this.state.inverse?.jacobian?.targetReachability && typeof RAD.buildResponseJacobian === "function") {
        RAD.buildResponseJacobian(this.state);
      }
      const report = this.state.inverse?.jacobian?.targetReachability || this.state.inverse?.linearSolution?.targetReachability;
      const cell = report?.worstUnderactuatedCell;
      if (!cell) return;
      this.state.selection = { r: cell.row, c: cell.col };
      this.state.view.overlayMode = "underactuated";
      this.state.view.targetVisible = this.state.target.type !== "none";
      this.syncControls();
      this.onChange(this.state);
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
          if (s.cells.positionLocked?.[r]?.[c]) button.classList.add("is-position-locked");
          if (s.cells.actuatorAllowed?.[r]?.[c] === false) button.classList.add("is-disallowed");
          if (s.cells.removed?.[r]?.[c]) button.classList.add("is-removed");
          const topology = s.cells.removed?.[r]?.[c] ? "removed, " : "";
          const fixture = s.cells.positionLocked?.[r]?.[c] ? "position locked, " : "";
          button.title = `${topology}${fixture}row ${r}, col ${c}, alpha ${sim.alpha[r][c].toFixed(3)}, z ${sim.height[r][c].toFixed(3)}`;
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
      const validation = this.state.inverse?.physicalValidation;
      this.els.inverseStepPreview.max = commands.length;
      this.els.inverseStepPreview.value = preview?.type === "plan-step" ? preview.step : 0;
      document.getElementById("inverseStepPreviewOut").textContent = `${this.els.inverseStepPreview.value} / ${commands.length}`;
      document.getElementById("inversePlanSummary").textContent = commands.length
        ? `${commands.length} actuators, error ${Number(plan.projectedError || 0).toFixed(3)}`
        : "No plan analyzed";
      document.getElementById("inversePlanDelta").textContent = `delta ${Number(plan?.totalImprovement || 0).toFixed(3)}`;
      const reachability = this.state.inverse?.jacobian?.targetReachability || this.state.inverse?.linearSolution?.targetReachability;
      document.getElementById("inverseReachabilitySummary").textContent = reachability
        ? `under targets ${reachability.underactuatedHeightCells || 0} up ${reachability.positiveUnderactuatedHeightCells || 0} down ${reachability.negativeUnderactuatedHeightCells || 0} topo ${reachability.topologyBlockedHeightCells || 0} comps ${reachability.topologyActuatedComponentCount || 0}/${reachability.topologyComponentCount || 0}, rms ${Number(reachability.unreachableHeightRms || 0).toFixed(3)}`
        : "under targets 0";
      document.getElementById("inversePreviewSummary").textContent = preview
        ? this.formatInversePreview(preview)
        : linearCommands.length
          ? `linear fit: ${linearCommands.length} actuators, err ${Number(this.state.inverse.linearSolution.projectedError || 0).toFixed(3)}`
          : "No actuator preview";
      const packet = this.state.experiment?.reachableEquilibriumProfileInversePreviewPacket;
      const replay = this.state.experiment?.reachableEquilibriumProfileInversePreviewReplay;
      const packetPhysical = this.state.experiment?.reachableEquilibriumProfileInversePreviewPhysical;
      document.getElementById("inversePacketSummary").textContent = packet
        ? `packet ${packet.summary?.profileInversePreviewPacketReady ? "ready" : "review"}: ${packet.summary?.commandCount || 0} cmds, missing ${packet.summary?.missingEvidenceCount || 0}`
        : "packet not built";
      document.getElementById("inversePhysicalSummary").textContent = validation
        ? `${validation.source} physical err ${Number(validation.physicalProjectedError || 0).toFixed(3)}, d ${Number(validation.physicalErrorDelta || 0).toFixed(3)}`
        : "physical validation not run";
      document.getElementById("inversePhysicalModel").textContent = validation
        ? `model rms h${Number(validation.heightModelRms || 0).toFixed(3)} c${Number(validation.centerModelRms || 0).toFixed(3)}`
        : "model delta 0.000";
      document.getElementById("inversePacketReplaySummary").textContent = replay
        ? `replay ${replay.summary?.profileInversePreviewReplayReady ? "ready" : "review"}: ${replay.commands?.replayedCommandCount || 0} cmds, miss ${replay.summary?.missingEvidenceCount || 0}`
        : "packet replay not run";
      document.getElementById("inversePacketPhysicalSummary").textContent = packetPhysical
        ? `packet phys ${packetPhysical.summary?.profileInversePreviewPhysicalReady ? "ready" : "review"}: h${Number(packetPhysical.comparison?.heightRmsModelError || 0).toFixed(3)} c${Number(packetPhysical.comparison?.centerRmsModelError || 0).toFixed(3)}`
        : "packet physical not run";
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

    saveCalibrationPlan() {
      const payload =
        typeof RAD.exportCalibrationMeasurementPlan === "function"
          ? RAD.exportCalibrationMeasurementPlan(this.state)
          : JSON.stringify({ schema: "rad-sim.calibration-plan.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-plan.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveExperimentProtocol() {
      const payload =
        typeof RAD.exportCalibrationExperimentProtocol === "function"
          ? RAD.exportCalibrationExperimentProtocol(this.state)
          : JSON.stringify({ schema: "rad-sim.calibration-experiment-protocol.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-experiment-protocol.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveResponseMatrix() {
      const scope = this.els.characterizationScope.value || "single";
      const payload =
        typeof RAD.exportResponseMatrix === "function"
          ? RAD.exportResponseMatrix(this.state, { scope, ...this.state.selection })
          : JSON.stringify({ schema: "rad-sim.response-matrix.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = `rad-sim-response-matrix-${scope}.json`;
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveTopologyReport() {
      const payload =
        typeof RAD.exportTopologyExperimentReport === "function"
          ? RAD.exportTopologyExperimentReport(this.state, null, { center: this.state.selection })
          : JSON.stringify({ schema: "rad-sim.browser-topology-experiment-report.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-topology-experiment-report.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveResponseAtlas() {
      const payload =
        typeof RAD.exportResponseAtlas === "function"
          ? RAD.exportResponseAtlas(this.state, { ...this.state.selection })
          : JSON.stringify({ schema: "rad-sim.response-atlas.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-response-atlas.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveResponseAtlasSweep() {
      const payload =
        typeof RAD.exportResponseAtlasSweep === "function"
          ? RAD.exportResponseAtlasSweep(this.state, { ...this.state.selection })
          : JSON.stringify({ schema: "rad-sim.response-atlas-sweep.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-response-atlas-sweep.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveProgrammableReport() {
      const scope = this.els.characterizationScope.value || "single";
      const payload =
        typeof RAD.exportProgrammableDiscontinuityReport === "function"
          ? RAD.exportProgrammableDiscontinuityReport(this.state, { scope, ...this.state.selection })
          : JSON.stringify({ schema: "rad-sim.programmable-discontinuity-report.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = `rad-sim-programmable-discontinuity-report-${scope}.json`;
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveFormalizationTargets() {
      const scope = this.els.characterizationScope.value || "single";
      const payload =
        typeof RAD.exportFormalizationTargetManifest === "function"
          ? RAD.exportFormalizationTargetManifest(this.state, { scope, ...this.state.selection })
          : JSON.stringify({ schema: "rad-sim.formalization-targets.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = `rad-sim-formalization-targets-${scope}.json`;
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveInverseReport() {
      const payload =
        typeof RAD.exportInverseDesignReport === "function"
          ? RAD.exportInverseDesignReport(this.state)
          : JSON.stringify({ schema: "rad-sim.inverse-design-report.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-inverse-design-report.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    downloadText(payload, filename, type = "application/json") {
      const blob = new Blob([payload], { type });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    currentProfileInverseReport() {
      if (typeof RAD.reachableEquilibriumProfileInverseReport !== "function") return null;
      return RAD.reachableEquilibriumProfileInverseReport(
        this.state,
        this.state.experiment?.reachableEquilibriumEmpiricalProfile || null,
        { store: true }
      );
    }

    currentProfileInverseAcceptance(report = null) {
      if (typeof RAD.reachableEquilibriumProfileInverseAcceptanceReport !== "function") return null;
      const profileInverse = report || this.state.experiment?.reachableEquilibriumProfileInverse || this.currentProfileInverseReport();
      const acceptance = RAD.reachableEquilibriumProfileInverseAcceptanceReport(profileInverse, {
        maxActiveActuators: this.state.inverse?.maxActuators || 24,
      });
      this.state.experiment.reachableEquilibriumProfileInverseAcceptance = JSON.parse(JSON.stringify(acceptance));
      return acceptance;
    }

    buildProfileInversePreviewPacket() {
      if (typeof RAD.reachableEquilibriumProfileInversePreviewPacket !== "function") return null;
      const profileInverse = this.currentProfileInverseReport();
      const acceptance = this.currentProfileInverseAcceptance(profileInverse);
      const packet = RAD.reachableEquilibriumProfileInversePreviewPacket(this.state, profileInverse, acceptance, {
        store: true,
        packetId: `browser-preview-${Date.now()}`,
      });
      RAD.recordEvent(this.state, {
        type: "profile-inverse-preview-packet-built",
        ready: packet.summary?.profileInversePreviewPacketReady,
        commands: packet.summary?.commandCount || 0,
        missing: packet.summary?.missingEvidenceCount || 0,
      });
      this.syncControls();
      this.onChange(this.state);
      return packet;
    }

    runProfileInversePreviewReplay() {
      if (typeof RAD.reachableEquilibriumProfileInversePreviewReplayReport !== "function") return null;
      const packet = this.state.experiment?.reachableEquilibriumProfileInversePreviewPacket || this.buildProfileInversePreviewPacket();
      const replay = RAD.reachableEquilibriumProfileInversePreviewReplayReport(this.state, packet, {
        store: true,
        residualAgreementTolerance: 1e-8,
      });
      RAD.recordEvent(this.state, {
        type: "profile-inverse-preview-replayed",
        ready: replay.summary?.profileInversePreviewReplayReady,
        commands: replay.commands?.replayedCommandCount || 0,
        missing: replay.summary?.missingEvidenceCount || 0,
      });
      this.state.view.overlayMode = "targetError";
      this.state.view.targetVisible = this.state.target.type !== "none";
      this.syncControls();
      this.onChange(this.state);
      return replay;
    }

    runProfileInversePreviewPhysical() {
      if (typeof RAD.reachableEquilibriumProfileInversePreviewPhysicalReport !== "function") return null;
      const packet = this.state.experiment?.reachableEquilibriumProfileInversePreviewPacket || this.buildProfileInversePreviewPacket();
      const physical = RAD.reachableEquilibriumProfileInversePreviewPhysicalReport(this.state, packet, {
        store: true,
        residualAgreementTolerance: 1e-8,
      });
      RAD.recordEvent(this.state, {
        type: "profile-inverse-preview-physical-checked",
        ready: physical.summary?.profileInversePreviewPhysicalReady,
        commands: physical.commands?.replayedCommandCount || 0,
        modelError: physical.comparison?.centerRmsModelError || 0,
        missing: physical.summary?.missingEvidenceCount || 0,
      });
      this.state.view.simulationMode = "springPreview";
      this.els.simulationMode.value = "springPreview";
      this.state.view.overlayMode = "modelError";
      this.state.view.targetVisible = this.state.target.type !== "none";
      this.syncControls();
      this.onChange(this.state);
      return physical;
    }

    saveProfileInversePreviewPacket() {
      const packet = this.state.experiment?.reachableEquilibriumProfileInversePreviewPacket || this.buildProfileInversePreviewPacket();
      const payload =
        typeof RAD.exportReachableEquilibriumProfileInversePreviewPacket === "function"
          ? RAD.exportReachableEquilibriumProfileInversePreviewPacket(
              this.state,
              this.state.experiment?.reachableEquilibriumProfileInverse || null,
              this.state.experiment?.reachableEquilibriumProfileInverseAcceptance || null,
              { packetId: packet?.packetId || "profile-inverse-preview" }
            )
          : JSON.stringify(packet || { schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1" }, null, 2);
      this.downloadText(payload, "rad-sim-profile-inverse-preview-packet.json");
    }

    saveProfileInversePreviewReplay() {
      const replay = this.state.experiment?.reachableEquilibriumProfileInversePreviewReplay || this.runProfileInversePreviewReplay();
      const payload =
        typeof RAD.exportReachableEquilibriumProfileInversePreviewReplay === "function"
          ? RAD.exportReachableEquilibriumProfileInversePreviewReplay(
              this.state,
              this.state.experiment?.reachableEquilibriumProfileInversePreviewPacket || null,
              { residualAgreementTolerance: 1e-8 }
            )
          : JSON.stringify(replay || { schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1" }, null, 2);
      this.downloadText(payload, "rad-sim-profile-inverse-preview-replay.json");
    }

    saveProfileInversePreviewPhysical() {
      const physical = this.state.experiment?.reachableEquilibriumProfileInversePreviewPhysical || this.runProfileInversePreviewPhysical();
      const payload =
        typeof RAD.exportReachableEquilibriumProfileInversePreviewPhysical === "function"
          ? RAD.exportReachableEquilibriumProfileInversePreviewPhysical(
              this.state,
              this.state.experiment?.reachableEquilibriumProfileInversePreviewPacket || null,
              { residualAgreementTolerance: 1e-8 }
            )
          : JSON.stringify(physical || { schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1" }, null, 2);
      this.downloadText(payload, "rad-sim-profile-inverse-preview-physical.json");
    }

    saveResultsTemplate() {
      const payload =
        typeof RAD.exportCalibrationExperimentResultsTemplate === "function"
          ? RAD.exportCalibrationExperimentResultsTemplate(this.state)
          : JSON.stringify({ schema: "rad-sim.calibration-experiment-results.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-results-template.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveBenchNotebook() {
      const payload =
        typeof RAD.exportCalibrationBenchNotebook === "function"
          ? RAD.exportCalibrationBenchNotebook(this.state)
          : JSON.stringify({ schema: "rad-sim.calibration-bench-notebook.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-bench-notebook.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveBenchPacket() {
      const payload =
        typeof RAD.exportCalibrationBenchPacket === "function"
          ? RAD.exportCalibrationBenchPacket(this.state)
          : JSON.stringify({ schema: "rad-sim.calibration-bench-packet.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-bench-packet.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveBenchNotebookCsv() {
      const payload =
        typeof RAD.exportCalibrationBenchNotebookCsv === "function"
          ? RAD.exportCalibrationBenchNotebookCsv(this.state)
          : "schema,step_id\nrad-sim.calibration-bench-notebook.v1,";
      const blob = new Blob([payload], { type: "text/csv" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-bench-notebook.csv";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveComparisonReport() {
      if (!this.state.experiment.calibrationComparison) {
        window.alert("Load calibration results before saving a comparison report.");
        return;
      }
      const payload =
        typeof RAD.exportCalibrationComparisonReport === "function"
          ? RAD.exportCalibrationComparisonReport(this.state)
          : JSON.stringify(this.state.experiment.calibrationComparison, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-comparison-report.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveModelProfile() {
      if (!this.state.experiment.calibrationComparison) {
        window.alert("Load calibration results before saving a model profile.");
        return;
      }
      const payload =
        typeof RAD.exportCalibrationModelProfile === "function"
          ? RAD.exportCalibrationModelProfile(this.state)
          : JSON.stringify({ schema: "rad-sim.calibration-model-profile.v1" }, null, 2);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-model-profile.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    compareModelProfile() {
      if (!this.state.experiment.calibrationComparison || typeof RAD.calibrationModelProfileResidualComparison !== "function") {
        window.alert("Load calibration results before checking a model profile.");
        return;
      }
      const profile =
        typeof RAD.calibrationModelProfile === "function"
          ? RAD.calibrationModelProfile(this.state)
          : this.state.experiment.calibrationModelProfile;
      const comparison = RAD.calibrationModelProfileResidualComparison(this.state, profile);
      this.state.experiment.calibrationModelProfile = profile;
      this.state.experiment.calibrationModelProfileComparison = comparison;
      const priorComparisons = (this.state.experiment.calibrationModelProfileHistory || [])
        .map((entry) => entry.residualComparison)
        .filter(Boolean);
      this.state.experiment.calibrationModelProfileSelection =
        typeof RAD.selectCalibrationModelProfile === "function"
          ? RAD.selectCalibrationModelProfile([...priorComparisons, comparison])
          : null;
      if (!Array.isArray(this.state.experiment.calibrationModelProfileHistory)) this.state.experiment.calibrationModelProfileHistory = [];
      this.state.experiment.calibrationModelProfileHistory.push({
        at: new Date().toISOString(),
        kind: "residual-comparison",
        residualComparison: comparison,
      });
      if (this.state.experiment.calibrationModelProfileHistory.length > 50) {
        this.state.experiment.calibrationModelProfileHistory.splice(0, this.state.experiment.calibrationModelProfileHistory.length - 50);
      }
      RAD.recordEvent(this.state, {
        type: "calibration-model-profile-check",
        improved: comparison.improved,
        residualScoreDelta: comparison.delta?.residualScore ?? null,
        selectedIndex: this.state.experiment.calibrationModelProfileSelection?.selectedIndex ?? null,
      });
      this.renderCalibrationResults();
      this.onChange(this.state);
    }

    saveModelProfileComparison() {
      if (!this.state.experiment.calibrationComparison || typeof RAD.exportCalibrationModelProfileResidualComparison !== "function") {
        window.alert("Load calibration results before saving a profile check.");
        return;
      }
      const payload = this.state.experiment.calibrationModelProfileComparison
        ? JSON.stringify(this.state.experiment.calibrationModelProfileComparison, null, 2)
        : RAD.exportCalibrationModelProfileResidualComparison(this.state, this.state.experiment.calibrationModelProfile);
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-model-profile-residual-comparison.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    checkModelProfileHoldout() {
      if (!this.state.experiment.calibrationComparison || !this.state.experiment.calibrationHoldoutResults) {
        window.alert("Load both fit calibration results and holdout results before checking holdout validation.");
        return;
      }
      if (typeof RAD.calibrationModelProfileHoldoutValidation !== "function") return;
      const profile =
        typeof RAD.calibrationModelProfile === "function"
          ? RAD.calibrationModelProfile(this.state)
          : this.state.experiment.calibrationModelProfile;
      const validation = RAD.calibrationModelProfileHoldoutValidation(this.state, this.state.experiment.calibrationHoldoutResults, profile);
      this.state.experiment.calibrationModelProfile = profile;
      this.state.experiment.calibrationModelProfileHoldoutValidation = validation;
      this.state.experiment.calibrationBenchExecutionValidation =
        typeof RAD.calibrationBenchExecutionValidation === "function"
          ? RAD.calibrationBenchExecutionValidation(this.state, this.state.experiment.calibrationHoldoutResults, profile)
          : null;
      this.state.experiment.calibrationModelProfileComparison = validation.fit;
      if (!Array.isArray(this.state.experiment.calibrationModelProfileHistory)) this.state.experiment.calibrationModelProfileHistory = [];
      this.state.experiment.calibrationModelProfileHistory.push({
        at: new Date().toISOString(),
        kind: "holdout-validation",
        holdoutValidation: validation,
        executionValidation: this.state.experiment.calibrationBenchExecutionValidation,
        residualComparison: validation.holdout,
      });
      if (this.state.experiment.calibrationModelProfileHistory.length > 50) {
        this.state.experiment.calibrationModelProfileHistory.splice(0, this.state.experiment.calibrationModelProfileHistory.length - 50);
      }
      RAD.recordEvent(this.state, {
        type: "calibration-model-profile-holdout",
        fitPass: validation.fitPass,
        holdoutPass: validation.holdoutPass,
        independentValidationPass: validation.independentValidationPass,
        holdoutResidualScoreDelta: validation.holdout?.delta?.residualScore ?? null,
      });
      this.renderCalibrationResults();
      this.onChange(this.state);
    }

    saveModelProfileHoldout() {
      if (!this.state.experiment.calibrationComparison || !this.state.experiment.calibrationHoldoutResults) {
        window.alert("Load both fit calibration results and holdout results before saving holdout validation.");
        return;
      }
      if (typeof RAD.exportCalibrationModelProfileHoldoutValidation !== "function") return;
      const payload = this.state.experiment.calibrationModelProfileHoldoutValidation
        ? JSON.stringify(this.state.experiment.calibrationModelProfileHoldoutValidation, null, 2)
        : RAD.exportCalibrationModelProfileHoldoutValidation(
            this.state,
            this.state.experiment.calibrationHoldoutResults,
            this.state.experiment.calibrationModelProfile
          );
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-model-profile-holdout-validation.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveModelProfileHoldoutCsv() {
      if (!this.state.experiment.calibrationModelProfileHoldoutValidation) {
        window.alert("Run Check Holdout before saving a holdout CSV table.");
        return;
      }
      if (typeof RAD.exportCalibrationModelProfileHoldoutValidationCsv !== "function") return;
      const payload = RAD.exportCalibrationModelProfileHoldoutValidationCsv(
        this.state.experiment.calibrationModelProfileHoldoutValidation
      );
      const blob = new Blob([payload], { type: "text/csv" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-model-profile-holdout-validation.csv";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveBenchExecutionValidation() {
      if (!this.state.experiment.calibrationComparison || !this.state.experiment.calibrationHoldoutResults) {
        window.alert("Load both fit calibration results and holdout results before saving an execution gate.");
        return;
      }
      if (typeof RAD.exportCalibrationBenchExecutionValidation !== "function") return;
      const payload = this.state.experiment.calibrationBenchExecutionValidation
        ? JSON.stringify(this.state.experiment.calibrationBenchExecutionValidation, null, 2)
        : RAD.exportCalibrationBenchExecutionValidation(
            this.state,
            this.state.experiment.calibrationHoldoutResults,
            this.state.experiment.calibrationModelProfile
          );
      const blob = new Blob([payload], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-bench-execution-validation.json";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    saveBenchExecutionValidationCsv() {
      if (!this.state.experiment.calibrationBenchExecutionValidation) {
        window.alert("Run Check Holdout before saving an execution-gate CSV table.");
        return;
      }
      if (typeof RAD.exportCalibrationBenchExecutionValidationCsv !== "function") return;
      const payload = RAD.exportCalibrationBenchExecutionValidationCsv(
        this.state.experiment.calibrationBenchExecutionValidation
      );
      const blob = new Blob([payload], { type: "text/csv" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "rad-sim-calibration-bench-execution-validation.csv";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    applyModelProfile() {
      if (!this.state.experiment.calibrationComparison || typeof RAD.applyCalibrationModelProfile !== "function") {
        window.alert("Load calibration results before applying a model profile.");
        return;
      }
      const profile =
        typeof RAD.calibrationModelProfile === "function"
          ? RAD.calibrationModelProfile(this.state)
          : null;
      if (typeof RAD.calibrationModelProfileResidualComparison === "function") {
        const comparison = RAD.calibrationModelProfileResidualComparison(this.state, profile);
        this.state.experiment.calibrationModelProfileComparison = comparison;
        const priorComparisons = (this.state.experiment.calibrationModelProfileHistory || [])
          .map((entry) => entry.residualComparison)
          .filter(Boolean);
        if (typeof RAD.selectCalibrationModelProfile === "function") {
          this.state.experiment.calibrationModelProfileSelection = RAD.selectCalibrationModelProfile([...priorComparisons, comparison]);
        }
      }
      const audit = RAD.applyCalibrationModelProfile(this.state, profile);
      RAD.recordEvent(this.state, {
        type: "calibration-model-profile-apply",
        applied: audit.appliedUpdates.length,
        skipped: audit.skippedUpdates.length,
        improved: this.state.experiment.calibrationModelProfileComparison?.improved ?? null,
        zCouplingGain: this.state.grid.zCouplingGain,
      });
      this.syncControls();
      this.renderCalibrationResults();
      this.onChange(this.state);
    }

    loadCalibrationResults() {
      const file = this.els.calibrationResultsFileInput.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        try {
          const results = JSON.parse(String(reader.result));
          results.provenance = {
            ...(results.provenance || {}),
            schema: "rad-sim.calibration-dataset-provenance.v1",
            sourceFileId: results.provenance?.sourceFileId || file.name,
            datasetRole: results.provenance?.datasetRole || "fit",
          };
          const comparison = RAD.compareCalibrationExperimentResults(this.state, results);
          const summary =
            typeof RAD.summarizeCalibrationComparison === "function"
              ? RAD.summarizeCalibrationComparison(comparison)
              : { stepCount: comparison.comparisons?.length || 0 };
          this.state.experiment.calibrationResults = results;
          this.state.experiment.calibrationComparison = comparison;
          this.state.experiment.calibrationComparisonSummary = summary;
          if (typeof RAD.calibrationModelProfile === "function") {
            this.state.experiment.calibrationModelProfile = RAD.calibrationModelProfile(this.state);
          }
          this.state.experiment.calibrationModelProfileComparison = null;
          this.state.experiment.calibrationModelProfileSelection = null;
          this.state.experiment.calibrationModelProfileHoldoutValidation = null;
          this.state.experiment.calibrationBenchExecutionValidation = null;
          this.state.view.overlayMode = "calibrationError";
          this.syncControls();
          this.renderCalibrationResults();
          this.onChange(this.state);
        } catch (error) {
          window.alert(error.message);
        }
      };
      reader.readAsText(file);
      this.els.calibrationResultsFileInput.value = "";
    }

    loadCalibrationHoldoutResults() {
      const file = this.els.calibrationHoldoutFileInput.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        try {
          const results = JSON.parse(String(reader.result));
          results.provenance = {
            ...(results.provenance || {}),
            schema: "rad-sim.calibration-dataset-provenance.v1",
            sourceFileId: results.provenance?.sourceFileId || file.name,
            datasetRole: results.provenance?.datasetRole || "holdout",
          };
          const comparison = RAD.compareCalibrationExperimentResults(this.state, results);
          const summary =
            typeof RAD.summarizeCalibrationComparison === "function"
              ? RAD.summarizeCalibrationComparison(comparison)
              : { stepCount: comparison.comparisons?.length || 0 };
          this.state.experiment.calibrationHoldoutResults = results;
          this.state.experiment.calibrationHoldoutComparison = comparison;
          this.state.experiment.calibrationHoldoutComparisonSummary = summary;
          this.state.experiment.calibrationModelProfileHoldoutValidation = null;
          this.state.experiment.calibrationBenchExecutionValidation = null;
          RAD.recordEvent(this.state, {
            type: "calibration-holdout-results-load",
            steps: summary.stepCount || 0,
            missing: summary.missingObservationCount || 0,
          });
          this.syncControls();
          this.renderCalibrationResults();
          this.onChange(this.state);
        } catch (error) {
          window.alert(error.message);
        }
      };
      reader.readAsText(file);
      this.els.calibrationHoldoutFileInput.value = "";
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
        const positionLocked = event.positionLocked ? " fixed" : "";
        return `cell command: r${event.r}, c${event.c}${alpha}${z}${locked}${positionLocked}`;
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
      if (event.type === "jacobian-built") return `jacobian: ${event.columns || 0} columns, reach ${Number(event.meanCoverage || 0).toFixed(3)}, under ${event.underactuatedTargets || 0}`;
      if (event.type === "linear-fit-solved") return `linear solve: ${event.actuators || 0} actuators, err ${Number(event.projectedError || 0).toFixed(3)}, under ${event.underactuatedTargets || 0}`;
      if (event.type === "linear-fit-applied") return `linear apply: ${event.actuators || 0} actuators, err ${Number(event.projectedError || 0).toFixed(3)}`;
      if (event.type === "inverse-physical-validated") return `physical validate: ${event.source || "plan"}, err ${Number(event.physicalError || 0).toFixed(3)}, model ${Number(event.modelError || 0).toFixed(3)}`;
      if (event.type === "profile-inverse-preview-packet-built") return `packet gate: ${event.ready ? "ready" : "review"}, ${event.commands || 0} cmds, missing ${event.missing || 0}`;
      if (event.type === "profile-inverse-preview-replayed") return `packet replay: ${event.ready ? "ready" : "review"}, ${event.commands || 0} cmds, missing ${event.missing || 0}`;
      if (event.type === "profile-inverse-preview-physical-checked") return `packet physical: ${event.ready ? "ready" : "review"}, model ${Number(event.modelError || 0).toFixed(3)}, missing ${event.missing || 0}`;
      if (event.type === "two-cell-radius-backlash-transition-report-run") return `transition report: ${event.brackets || 0} brackets, ${event.measurementCases || 0} cases`;
      if (event.type === "two-cell-radius-backlash-transition-measurements-load") return `transition data: ${event.matched || 0}/${event.rows || 0} rows, ${event.shift || "unknown"}`;
      if (event.type === "two-cell-radius-backlash-transition-calibration-applied") return `transition apply: hole ${Number(event.holeRadius || 0).toFixed(3)}, b ${Number(event.backlash || 0).toFixed(3)}`;
      if (event.type === "characterization") return `characterize: ${event.scope}, ${event.responseCells || 0} cells, sup ${Number(event.superpositionError || 0).toFixed(3)}, pairs ${event.pairwiseNonadditive || 0}/${event.pairwisePairs || 0}, deg ${event.pairwiseMaxDegree || 0}`;
      if (event.type === "operator-lock") return `event lock: r${event.r}, c${event.c} a ${Number(event.lockAlpha || 0).toFixed(3)} z ${Number(event.lockZ || 0).toFixed(3)}`;
      if (event.type === "operator-release") return `event release: r${event.r}, c${event.c}`;
      if (event.type === "operator-remove-cell") return `event remove: r${event.r}, c${event.c}`;
      if (event.type === "operator-restore-cell") return `event restore: r${event.r}, c${event.c}`;
      if (event.type === "operator-order-check") return `order check: r${event.r}, c${event.c} da ${Number(event.finalAlphaError || 0).toFixed(3)} dz ${Number(event.finalHeightError || 0).toFixed(3)} max ${Number(event.maxOrderError || 0).toFixed(3)}`;
      if (event.type === "position-lock-row-ends") return `position lock row ends: r${event.row}`;
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
