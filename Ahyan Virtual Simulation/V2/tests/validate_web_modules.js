const assert = require("assert");
const fs = require("fs");
const http = require("http");
const path = require("path");
const vm = require("vm");

const root = path.resolve(__dirname, "..");
const web = path.join(root, "web");
const context = {
  console,
  window: {},
  performance,
};
context.window.window = context.window;
context.window.console = console;
context.window.performance = performance;
context.document = {
  body: { classList: { toggle: () => {} } },
  getElementById: (id) => ({
    id,
    classList: { toggle: () => {} },
    setAttribute: () => {},
    addEventListener: () => {},
    replaceChildren: () => {},
    append: () => {},
    textContent: "",
    title: "",
    hidden: false,
  }),
  querySelectorAll: () => [],
  createElement: (tag) => ({
    tagName: String(tag).toUpperCase(),
    classList: { toggle: () => {} },
    setAttribute: () => {},
    addEventListener: () => {},
    append: () => {},
    replaceChildren: () => {},
    textContent: "",
  }),
};
context.window.document = context.document;
vm.createContext(context);

for (const filename of ["state.js", "provenance.js", "math.js", "constraints.js", "operators.js", "inverse.js", "analysis.js", "data/lock_dataset.js", "data/two_cell_external_fidelity_summary.js", "lock_surrogate.js", "physics.js", "two_cell_bench.js", "mesh_export.js", "primitives.js"]) {
  const source = fs.readFileSync(path.join(web, filename), "utf8");
  vm.runInContext(source, context, { filename });
}

const RAD = context.window.RAD;
assert.ok(RAD, "RAD namespace should load");
assert.strictEqual(typeof RAD.modelProvenance, "function", "provenance module should expose model provenance");
assert.strictEqual(typeof RAD.provenanceSummary, "function", "provenance module should expose provenance counts");
assert.strictEqual(typeof RAD.calibrationProfileSummary, "function", "math module should expose calibration profile summary");
assert.strictEqual(typeof RAD.calibrationReadiness, "function", "math module should expose calibration readiness");
assert.strictEqual(typeof RAD.calibrationMeasurementPlan, "function", "math module should expose calibration measurement plan");
assert.strictEqual(typeof RAD.exportCalibrationMeasurementPlan, "function", "math module should export calibration measurement plan");
assert.strictEqual(typeof RAD.applyHardwareProfileToGrid, "function", "math module should apply hardware profiles");
assert.strictEqual(typeof RAD.exportHardwareProfileJson, "function", "math module should export hardware profiles");
assert.strictEqual(typeof RAD.importHardwareProfileJson, "function", "math module should import hardware profiles");
assert.strictEqual(typeof RAD.solveBoundaryConstraints, "function", "constraints module should expose boundary solver");
assert.strictEqual(typeof RAD.fitBoundaryToReference, "function", "constraints module should fit wall bounds");
const freeConstraintState = RAD.createState(5, 5);
RAD.clearCommands(freeConstraintState);
freeConstraintState.cells.commandAlpha[2][2] = 0.35;
let constrainedRule = RAD.simulate(freeConstraintState);
assert.ok(constrainedRule.metrics.maxInducedHeight < 1e-9, "free in-plane alpha should not create induced height");
assert.ok(constrainedRule.metrics.maxCompressionResidual < 1e-9, "free in-plane alpha should not create compression residual");
const wallConstraintState = RAD.createState(5, 5);
RAD.clearCommands(wallConstraintState);
wallConstraintState.cells.commandAlpha[2][2] = 0.35;
RAD.fitBoundaryToReference(wallConstraintState, 0);
constrainedRule = RAD.simulate(wallConstraintState);
assert.ok(constrainedRule.metrics.maxCompressionResidual > 0, "fitted rigid walls should create blocked XY residual");
assert.ok(constrainedRule.metrics.maxInducedHeight > 0, "blocked XY residual should create induced height");
assert.strictEqual(typeof RAD.calibratedMeshDimensions, "function", "mesh exporter should expose calibrated mesh dimensions");
assert.strictEqual(typeof RAD.simulateTwoCellBench, "function", "two-cell bench module should expose bench simulator");
assert.strictEqual(typeof RAD.sweepTwoCellBacklash, "function", "two-cell bench module should expose backlash sweep");
assert.strictEqual(typeof RAD.exportTwoCellBenchJson, "function", "two-cell bench module should export JSON bench packets");
assert.strictEqual(typeof RAD.exportTwoCellBacklashSweepCsv, "function", "two-cell bench module should export CSV sweeps");
assert.strictEqual(typeof RAD.cadRadCellLayout, "function", "two-cell bench module should expose CAD-derived one-cell layout");
assert.strictEqual(typeof RAD.twoCellConnectorContactReport, "function", "two-cell bench module should expose connector contact report");
assert.strictEqual(typeof RAD.twoCellPhysicalSimulationSuite, "function", "two-cell bench module should expose the named physical suite");
assert.strictEqual(typeof RAD.twoCellPhysicalSuiteCaseOptions, "function", "two-cell bench module should expose suite case controls");
assert.strictEqual(typeof RAD.exportTwoCellPhysicalSuiteCsv, "function", "two-cell bench module should export suite CSV");
assert.strictEqual(typeof RAD.twoCellPhysicalFidelityMatrix, "function", "two-cell bench module should expose the physical fidelity matrix");
assert.strictEqual(typeof RAD.exportTwoCellPhysicalFidelityMatrixJson, "function", "two-cell bench module should export matrix JSON");
assert.strictEqual(typeof RAD.exportTwoCellPhysicalFidelityMatrixCsv, "function", "two-cell bench module should export matrix CSV");
assert.strictEqual(typeof RAD.twoCellContactPhaseMap, "function", "two-cell bench module should expose the contact phase map");
assert.strictEqual(typeof RAD.exportTwoCellContactPhaseMapJson, "function", "two-cell bench module should export phase-map JSON");
assert.strictEqual(typeof RAD.exportTwoCellContactPhaseMapCsv, "function", "two-cell bench module should export phase-map CSV");
assert.strictEqual(typeof RAD.twoCellRadiusBacklashPhaseDiagram, "function", "two-cell bench module should expose the radius/backlash phase diagram");
assert.strictEqual(typeof RAD.exportTwoCellRadiusBacklashPhaseDiagramJson, "function", "two-cell bench module should export phase-diagram JSON");
assert.strictEqual(typeof RAD.exportTwoCellRadiusBacklashPhaseDiagramCsv, "function", "two-cell bench module should export phase-diagram CSV");
assert.strictEqual(typeof RAD.twoCellRadiusBacklashTransitionReport, "function", "two-cell bench module should expose the transition calibration report");
assert.strictEqual(typeof RAD.exportTwoCellRadiusBacklashTransitionReportJson, "function", "two-cell bench module should export transition-report JSON");
assert.strictEqual(typeof RAD.exportTwoCellRadiusBacklashTransitionReportCsv, "function", "two-cell bench module should export transition-report CSV");
assert.strictEqual(typeof RAD.compareTwoCellRadiusBacklashTransitionMeasurements, "function", "two-cell bench module should compare transition measurements");
assert.strictEqual(typeof RAD.exportTwoCellRadiusBacklashTransitionComparisonJson, "function", "two-cell bench module should export transition comparison JSON");
assert.strictEqual(typeof RAD.exportTwoCellRadiusBacklashTransitionComparisonCsv, "function", "two-cell bench module should export transition comparison CSV");
assert.strictEqual(typeof RAD.twoCellRadiusBacklashTransitionRerun, "function", "two-cell bench module should rerun calibrated transition diagrams");
assert.strictEqual(typeof RAD.exportTwoCellRadiusBacklashTransitionRerunJson, "function", "two-cell bench module should export transition rerun JSON");
assert.strictEqual(typeof RAD.exportTwoCellRadiusBacklashTransitionRerunCsv, "function", "two-cell bench module should export transition rerun CSV");
assert.strictEqual(typeof RAD.twoCellPhysicalResponseAtlas, "function", "two-cell bench module should expose the physical response atlas");
assert.strictEqual(typeof RAD.exportTwoCellPhysicalResponseAtlasJson, "function", "two-cell bench module should export atlas JSON");
assert.strictEqual(typeof RAD.exportTwoCellPhysicalResponseAtlasCsv, "function", "two-cell bench module should export atlas CSV");
assert.strictEqual(typeof RAD.twoCellCadContactDecompositionSpec, "function", "two-cell bench module should expose CAD contact decomposition specs");
assert.strictEqual(typeof RAD.exportTwoCellCadContactDecompositionJson, "function", "two-cell bench module should export CAD contact JSON");
assert.strictEqual(typeof RAD.exportTwoCellCadContactDecompositionCsv, "function", "two-cell bench module should export CAD contact CSV");
assert.strictEqual(typeof RAD.twoCellExactContactHandoffPlan, "function", "two-cell bench module should expose exact-contact handoff plans");
assert.strictEqual(typeof RAD.exportTwoCellExactContactHandoffPlanJson, "function", "two-cell bench module should export exact-contact handoff JSON");
assert.strictEqual(typeof RAD.exportTwoCellExactContactHandoffPlanCsv, "function", "two-cell bench module should export exact-contact handoff CSV");
assert.strictEqual(typeof RAD.twoCellFidelityMatrixMeasurementTemplate, "function", "two-cell bench module should expose matrix measurement templates");
assert.strictEqual(typeof RAD.exportTwoCellFidelityMatrixMeasurementTemplateCsv, "function", "two-cell bench module should export matrix measurement CSV");
assert.strictEqual(typeof RAD.twoCellExternalFidelityMatrixManifest, "function", "two-cell bench module should expose external fidelity manifests");
assert.strictEqual(typeof RAD.exportTwoCellExternalFidelityMatrixManifestJson, "function", "two-cell bench module should export external fidelity manifest JSON");
assert.strictEqual(typeof RAD.exportTwoCellExternalFidelityMatrixManifestCsv, "function", "two-cell bench module should export external fidelity manifest CSV");
assert.strictEqual(typeof RAD.twoCellFidelityMatrixMeasurementsFromCsv, "function", "two-cell bench module should parse matrix measurements");
assert.strictEqual(typeof RAD.compareTwoCellFidelityMatrixMeasurements, "function", "two-cell bench module should compare matrix measurements");
assert.strictEqual(typeof RAD.exportTwoCellFidelityMatrixMeasurementComparisonJson, "function", "two-cell bench module should export matrix measurement comparisons");
assert.strictEqual(typeof RAD.calibrateTwoCellFidelityMatrixParameters, "function", "two-cell bench module should calibrate matrix parameters");
assert.strictEqual(typeof RAD.exportTwoCellFidelityMatrixParameterCalibrationJson, "function", "two-cell bench module should export matrix parameter calibration");
assert.strictEqual(typeof RAD.cadRadCellArchiveAudit, "function", "two-cell bench module should expose the CAD archive audit");
assert.strictEqual(typeof RAD.twoCellExternalFidelityWebSummary, "function", "two-cell bench module should expose compact external-fidelity browser evidence");
assert.strictEqual(typeof RAD.twoCellExternalFidelityCaseOptions, "function", "two-cell bench module should expose external-fidelity representative case options");
assert.strictEqual(typeof RAD.twoCellExternalFidelityCase, "function", "two-cell bench module should expose external-fidelity representative case lookup");
assert.strictEqual(typeof RAD.twoCellExternalFidelityCaseControls, "function", "two-cell bench module should expose external-fidelity case controls");
assert.strictEqual(typeof RAD.exportTwoCellExternalFidelityWebSummaryJson, "function", "two-cell bench module should export compact external-fidelity browser evidence");
assert.strictEqual(typeof RAD.twoCellPhysicalFidelityStatus, "function", "two-cell bench module should expose physical fidelity status");
assert.strictEqual(typeof RAD.exportTwoCellPhysicalFidelityStatusJson, "function", "two-cell bench module should export physical fidelity status");
assert.strictEqual(typeof RAD.twoCellConnectorMeasurementTemplate, "function", "two-cell bench module should expose connector measurement templates");
assert.strictEqual(typeof RAD.exportTwoCellConnectorMeasurementTemplateCsv, "function", "two-cell bench module should export connector measurement CSV");
const twoCellBenchState = RAD.createState(1, 2);
const cadLayout = RAD.cadRadCellLayout(twoCellBenchState);
assert.strictEqual(cadLayout.schema, "rad-sim.cad-rad-cell-layout.v1", "CAD layout should expose stable schema");
assert.strictEqual(cadLayout.padSites.length, 8, "CAD layout should contain eight radial pad sites");
assert.strictEqual(cadLayout.envelopeCheck.matchesBoundingBox, true, "CAD layout should match the A360 bounding-box envelope");
assert.strictEqual(cadLayout.twoCellConnectorPairs.length, 3, "CAD layout should expose the three two-cell connector pairs");
assert.strictEqual(typeof RAD.cadRadCellReferenceProfile, "function", "two-cell bench module should expose the A360 CAD reference profile");
assert.strictEqual(typeof RAD.exportCadRadCellReferenceProfileJson, "function", "two-cell bench module should export the A360 CAD reference profile");
const cadReferenceProfile = RAD.cadRadCellReferenceProfile(twoCellBenchState);
assert.strictEqual(cadReferenceProfile.schema, "rad-sim.cad-rad-cell-reference-profile.v1", "CAD reference profile should expose stable schema");
assert.strictEqual(cadReferenceProfile.hardwareProfile.schema, "rad-sim.hardware-profile.v1", "CAD reference profile should wrap a hardware profile");
assert.strictEqual(Number(cadReferenceProfile.hardwareProfile.dimensionsMm.holeRadiusMm.toFixed(6)), 1.7, "CAD reference profile should preserve the 3.4 mm nominal hole");
assert.ok(cadReferenceProfile.claimBoundary.blockedClaim.includes("exact real-cell dynamics"), "CAD reference profile should preserve physical claim boundary");
const cadAudit = RAD.cadRadCellArchiveAudit();
assert.strictEqual(cadAudit.schema, "rad-sim.cad-rad-cell-archive-audit.v1", "CAD archive audit should expose stable schema");
assert.strictEqual(cadAudit.summary.brepEntryCount, 2, "CAD archive audit should preserve inspected BREP count");
assert.strictEqual(cadAudit.summary.physicalAccuracyValidated, false, "CAD archive audit should not claim physical validation");
assert.ok(cadAudit.summary.missingEvidence.includes("exactPinHoleContactSurfaces"), "CAD audit should require exact contact surfaces");
const externalFidelitySummary = RAD.twoCellExternalFidelityWebSummary();
assert.strictEqual(externalFidelitySummary.schema, "rad-sim.web-two-cell-external-fidelity-summary.v1", "external-fidelity web summary should expose stable schema");
assert.strictEqual(externalFidelitySummary.summary.externalRunComplete, true, "external-fidelity web summary should route completed MuJoCo matrix status");
assert.strictEqual(externalFidelitySummary.summary.physicalAccuracyValidated, false, "external-fidelity web summary should not claim hardware validation");
assert.ok(externalFidelitySummary.summary.correctedRms < externalFidelitySummary.summary.baselineRms, "external correction should improve compact preview RMS");
assert.ok(RAD.exportTwoCellExternalFidelityWebSummaryJson().includes("segmentedCadContactModel"), "external-fidelity web export should preserve remaining physical blockers");
const externalCaseOptions = RAD.twoCellExternalFidelityCaseOptions(externalFidelitySummary);
assert.ok(externalCaseOptions.length >= 3, "external-fidelity web summary should expose representative cases");
assert.strictEqual(externalCaseOptions[0].connectorCount, 3, "external representative cases should group upper/middle/lower connectors");
const externalCasePreview = RAD.twoCellExternalFidelityCaseControls(externalCaseOptions[0].caseId, externalFidelitySummary);
assert.strictEqual(externalCasePreview.schema, "rad-sim.browser-two-cell-external-fidelity-case-preview.v1", "external case preview should expose stable schema");
assert.strictEqual(externalCasePreview.grid.holeRadius, externalCaseOptions[0].holeRadius, "external case preview should route hole radius into the 1x2 bench");
assert.strictEqual(externalCasePreview.physicalAccuracyValidated, false, "external case preview should not claim hardware validation");
assert.ok(Number.isFinite(externalCasePreview.correctedPreview.rightCellZ), "external case preview should expose corrected right-cell height");
assert.ok(Array.isArray(externalCasePreview.connectorRows), "external case preview should expose connector rows for 3D overlay");
assert.strictEqual(externalCasePreview.connectorRows.length, externalCasePreview.connectorCount, "external connector rows should match connector count");
const twoCellBench = RAD.simulateTwoCellBench(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(twoCellBench.schema, "rad-sim.two-cell-physical-bench.v1", "two-cell bench should expose stable schema");
assert.strictEqual(twoCellBench.cadLayout.schema, cadLayout.schema, "two-cell bench should include CAD layout metadata");
assert.strictEqual(twoCellBench.cells.length, 2, "two-cell bench should report both connected cells");
const connectorContact = RAD.twoCellConnectorContactReport(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(connectorContact.schema, "rad-sim.two-cell-connector-contact-report.v1", "connector contact should expose stable schema");
assert.strictEqual(connectorContact.connectors.length, 3, "connector contact should report all three CAD connector pairs");
assert.ok(connectorContact.summary.maxVerticalSlipMm > 0, "vertical actuation should create connector-level slip");
const twoCellSweep = RAD.sweepTwoCellBacklash(twoCellBenchState, { holeSweepSteps: 5, holeSweepMax: 0.5 });
assert.strictEqual(twoCellSweep.schema, "rad-sim.two-cell-backlash-sweep.v1", "two-cell sweep should expose stable schema");
assert.strictEqual(twoCellSweep.rows.length, 5, "two-cell sweep should honor requested sweep count");
assert.strictEqual(twoCellSweep.trend.clearanceIncreasesBacklash, true, "larger holes should increase effective backlash");
const stateLockedBench = RAD.simulateTwoCellBench(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4, rightLocked: true });
assert.strictEqual(stateLockedBench.cells[1].alpha, twoCellBenchState.grid.initialAlpha, "state lock should freeze alpha");
assert.ok(stateLockedBench.cells[1].center.z > 0, "state lock should not suppress vertical motion");
const browserSuite = RAD.twoCellPhysicalSimulationSuite(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(browserSuite.schema, "rad-sim.two-cell-physical-simulation-suite.v1", "two-cell suite should expose stable schema");
assert.strictEqual(browserSuite.summary.caseCount, 12, "two-cell suite should include the named physical cases");
assert.strictEqual(browserSuite.summary.internalSuiteReady, true, "browser two-cell suite should be internally ready");
assert.strictEqual(browserSuite.summary.physicalAccuracyValidated, false, "browser suite should not claim physical validation");
assert.strictEqual(browserSuite.summary.backlashTrend.highBacklashReducesAlphaResponse, true, "high backlash should reduce alpha response");
assert.ok(RAD.exportTwoCellPhysicalSuiteCsv(twoCellBenchState).includes("leftAlphaDeltaFromInitial"), "suite CSV should expose alpha die-off fields");
const browserMatrix = RAD.twoCellPhysicalFidelityMatrix(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(browserMatrix.schema, "rad-sim.two-cell-physical-fidelity-matrix.v1", "fidelity matrix should expose stable schema");
assert.strictEqual(browserMatrix.summary.rowCount, 252, "fidelity matrix should cover 3 backlash x 3 holes x 4 locks x 7 actuations");
assert.strictEqual(browserMatrix.summary.connectorRowCount, 756, "fidelity matrix should include three connector rows per case");
assert.strictEqual(browserMatrix.summary.internalMatrixReady, true, "fidelity matrix internal checks should pass");
assert.strictEqual(browserMatrix.summary.physicalAccuracyValidated, false, "fidelity matrix should not claim physical validation");
assert.strictEqual(browserMatrix.summary.lockChecks.stateLockAlphaHold, true, "state locks should hold alpha in the matrix");
assert.strictEqual(browserMatrix.summary.lockChecks.positionLockFixtureHold, true, "position locks should hold fixture coordinates in the matrix");
assert.ok(browserMatrix.summary.missingEvidence.includes("externalRigidBodyContactRun"), "matrix should require external contact evidence");
assert.ok(RAD.exportTwoCellPhysicalFidelityMatrixCsv(twoCellBenchState).includes("connectorMaxVerticalSlipMm"), "matrix CSV should expose connector slip fields");
const browserPhaseMap = RAD.twoCellContactPhaseMap(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(browserPhaseMap.schema, "rad-sim.two-cell-contact-phase-map.v1", "phase map should expose stable schema");
assert.strictEqual(browserPhaseMap.summary.rowCount, browserMatrix.summary.rowCount, "phase map should classify every matrix row");
assert.strictEqual(browserPhaseMap.summary.physicalAccuracyValidated, false, "phase map should not claim hardware validation");
assert.ok(browserPhaseMap.summary.activeContactCaseCount > 0, "phase map should report active contact cases");
assert.ok(browserPhaseMap.summary.phaseCounts["position-locked"] > 0, "phase map should expose position-lock cases");
assert.ok(browserPhaseMap.summary.phaseCounts["state-locked-vertical-free"] > 0, "phase map should expose state-lock vertical-free cases");
assert.ok(browserPhaseMap.summary.measurementPriority.length > 0, "phase map should prioritize cases for measurement");
assert.ok(RAD.exportTwoCellContactPhaseMapCsv(twoCellBenchState).includes("phase"), "phase-map CSV should expose phase labels");
const browserPhaseDiagram = RAD.twoCellRadiusBacklashPhaseDiagram(twoCellBenchState, { holeSteps: 5, backlashSteps: 4 });
assert.strictEqual(browserPhaseDiagram.schema, "rad-sim.two-cell-radius-backlash-phase-diagram.v1", "phase diagram should expose stable schema");
assert.strictEqual(browserPhaseDiagram.summary.rowCount, 20, "phase diagram should sweep hole radius against backlash");
assert.strictEqual(browserPhaseDiagram.phaseGrid.length, 4, "phase diagram should expose backlash rows");
assert.strictEqual(browserPhaseDiagram.phaseGrid[0].length, 5, "phase diagram should expose hole-radius columns");
assert.strictEqual(browserPhaseDiagram.summary.physicalAccuracyValidated, false, "phase diagram should not claim hardware validation");
assert.ok(Object.keys(browserPhaseDiagram.summary.phaseCounts).length > 0, "phase diagram should classify phases");
assert.strictEqual(browserPhaseDiagram.axes.effectiveHoleRadii.length, 5, "phase diagram should expose effective hole radii");
assert.ok(RAD.exportTwoCellRadiusBacklashPhaseDiagramCsv(twoCellBenchState, { holeSteps: 3, backlashSteps: 3 }).includes("effectiveHoleRadius"), "phase-diagram CSV should expose effective radius columns");
const browserTransitionReport = RAD.twoCellRadiusBacklashTransitionReport(twoCellBenchState, { holeSteps: 5, backlashSteps: 4 });
assert.strictEqual(browserTransitionReport.schema, "rad-sim.two-cell-radius-backlash-transition-report.v1", "transition report should expose stable schema");
assert.strictEqual(browserTransitionReport.summary.physicalAccuracyValidated, false, "transition report should not claim hardware validation");
assert.ok(browserTransitionReport.summary.measurementCaseCount > 0, "transition report should prioritize measurement cases");
assert.ok(browserTransitionReport.measurementCases[0].measurementTarget.includes("connector marker"), "transition report should identify bench observables");
assert.ok(RAD.exportTwoCellRadiusBacklashTransitionReportCsv(twoCellBenchState, { holeSteps: 3, backlashSteps: 3 }).includes("observedPhase"), "transition CSV should be lab-fillable");
const transitionMeasurement = {
  ...browserTransitionReport.measurementCases[0],
  observedPhase: browserTransitionReport.measurementCases[0].phase,
  observedRightZ: browserTransitionReport.measurementCases[0].rightZ,
  observedRightAlpha: browserTransitionReport.measurementCases[0].rightAlpha,
  observedConnectorMaxVerticalSlipMm: browserTransitionReport.measurementCases[0].connectorMaxVerticalSlipMm,
};
const browserTransitionComparison = RAD.compareTwoCellRadiusBacklashTransitionMeasurements(
  twoCellBenchState,
  [transitionMeasurement],
  { holeSteps: 5, backlashSteps: 4 }
);
assert.strictEqual(browserTransitionComparison.schema, "rad-sim.two-cell-radius-backlash-transition-comparison.v1", "transition comparison should expose stable schema");
assert.strictEqual(browserTransitionComparison.summary.matchedRowCount, 1, "transition comparison should match filled rows");
assert.strictEqual(browserTransitionComparison.summary.phaseAccuracy, 1, "transition comparison should score matching phases");
assert.ok(RAD.exportTwoCellRadiusBacklashTransitionComparisonCsv(twoCellBenchState, [transitionMeasurement], { holeSteps: 5, backlashSteps: 4 }).includes("transitionShiftVote"), "transition comparison CSV should expose shift diagnostics");
const shiftedTransitionMeasurement = {
  ...browserTransitionReport.measurementCases[0],
  observedPhase: "axial-contact",
};
const shiftedTransitionComparison = RAD.compareTwoCellRadiusBacklashTransitionMeasurements(
  twoCellBenchState,
  [transitionMeasurement, shiftedTransitionMeasurement],
  { holeSteps: 5, backlashSteps: 4 }
);
assert.strictEqual(
  shiftedTransitionComparison.summary.transitionShiftDirection,
  "effective-hole-transition-lower-than-reduced-model",
  "transition comparison should diagnose lower-than-predicted hole boundary"
);
assert.ok(
  shiftedTransitionComparison.calibrationEstimate.proposedReducedProxyUpdates.holeRadius > twoCellBenchState.grid.holeRadius,
  "lower transition boundary should propose larger effective hole radius"
);
const shiftedTransitionRerun = RAD.twoCellRadiusBacklashTransitionRerun(
  twoCellBenchState,
  [transitionMeasurement, shiftedTransitionMeasurement],
  { holeSteps: 5, backlashSteps: 4 }
);
assert.strictEqual(shiftedTransitionRerun.schema, "rad-sim.two-cell-radius-backlash-transition-rerun.v1", "transition rerun should expose stable schema");
assert.ok(shiftedTransitionRerun.axisBiases.effectiveHoleRadiusBias > 0, "transition rerun should apply the effective hole-radius bias");
assert.strictEqual(
  shiftedTransitionRerun.calibrated.diagram.axes.holeRadiusBias,
  shiftedTransitionRerun.axisBiases.effectiveHoleRadiusBias,
  "calibrated rerun should feed the estimated bias into the phase diagram"
);
assert.ok(
  RAD.exportTwoCellRadiusBacklashTransitionRerunCsv(twoCellBenchState, [transitionMeasurement, shiftedTransitionMeasurement], { holeSteps: 5, backlashSteps: 4 }).includes("midpointDelta"),
  "transition rerun CSV should expose boundary movement columns"
);
const browserTwoCellAtlas = RAD.twoCellPhysicalResponseAtlas(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(browserTwoCellAtlas.schema, "rad-sim.two-cell-physical-response-atlas.v1", "response atlas should expose stable schema");
assert.strictEqual(browserTwoCellAtlas.matrixSummary.rowCount, browserMatrix.summary.rowCount, "response atlas should summarize the full matrix");
assert.ok(browserTwoCellAtlas.benchPriority.firstCasesToMeasure.length > 0, "response atlas should rank first bench cases");
assert.ok(browserTwoCellAtlas.rankedCases.maxVerticalSlip.length > 0, "response atlas should rank high-slip cases");
assert.strictEqual(browserTwoCellAtlas.invariants.physicalAccuracyValidated, false, "response atlas should not claim hardware validation");
assert.ok(RAD.exportTwoCellPhysicalResponseAtlasCsv(twoCellBenchState).includes("maxVerticalSlip"), "atlas CSV should expose ranked case categories");
const browserCadContactSpec = RAD.twoCellCadContactDecompositionSpec(twoCellBenchState, { ringSegments: 8 });
assert.strictEqual(browserCadContactSpec.schema, "rad-sim.two-cell-cad-contact-decomposition.v1", "CAD contact decomposition should expose stable schema");
assert.strictEqual(browserCadContactSpec.summary.holeCount, 16, "CAD contact decomposition should cover two cells with eight pad holes each");
assert.strictEqual(browserCadContactSpec.summary.twoCellConnectorCount, 3, "CAD contact decomposition should cover the three X-neighbor connector pairs");
assert.strictEqual(browserCadContactSpec.holeDecomposition[0].ringSegments.length, 8, "CAD contact decomposition should segment hole-wall contacts");
assert.ok(RAD.exportTwoCellCadContactDecompositionCsv(twoCellBenchState, { ringSegments: 8 }).includes("convex_ring_sector_or_capsule_wall"), "CAD contact CSV should expose collision primitives");
const browserExactContactPlan = RAD.twoCellExactContactHandoffPlan(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4, priorityLimit: 8 });
assert.strictEqual(browserExactContactPlan.schema, "rad-sim.two-cell-exact-contact-handoff-plan.v1", "exact-contact plan should expose stable schema");
assert.strictEqual(browserExactContactPlan.summary.canRunExactContact, false, "exact-contact plan should not claim readiness from one-cell CAD alone");
assert.ok(browserExactContactPlan.summary.caseCount > 0, "exact-contact plan should prioritize cases");
assert.ok(browserExactContactPlan.summary.missingEvidence.includes("exactPinHoleContactSurfaces"), "exact-contact plan should require pin-hole contact surfaces");
assert.ok(browserExactContactPlan.caseRows[0].requiredCadAssets.includes("upper_free_cell_body"), "exact-contact plan should list segmented body assets");
assert.ok(RAD.exportTwoCellExactContactHandoffPlanCsv(twoCellBenchState, { priorityLimit: 4 }).includes("requiredCadAssets"), "exact-contact CSV should expose CAD requirements");
const browserFidelityStatus = RAD.twoCellPhysicalFidelityStatus(twoCellBenchState, {}, { suite: browserSuite, matrix: browserMatrix });
assert.strictEqual(browserFidelityStatus.schema, "rad-sim.browser-two-cell-physical-fidelity-status.v1", "physical fidelity status should expose stable schema");
assert.strictEqual(browserFidelityStatus.summary.exactGeometryReady, false, "status should not claim exact CAD readiness");
assert.strictEqual(browserFidelityStatus.summary.matrixRunReady, true, "status should see the current matrix as internally ready");
assert.strictEqual(browserFidelityStatus.summary.matrixMeasurementRowsExpected, 756, "status should route the dense measurement row count");
assert.ok(RAD.exportTwoCellPhysicalFidelityStatusJson(twoCellBenchState, {}, { suite: browserSuite, matrix: browserMatrix }).includes("segmentedBodyMeshExport"), "status export should include exact-geometry blockers");
const browserMatrixTemplate = RAD.twoCellFidelityMatrixMeasurementTemplate(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(browserMatrixTemplate.schema, "rad-sim.two-cell-fidelity-matrix-measurement-template.v1", "matrix template should expose stable schema");
assert.strictEqual(browserMatrixTemplate.summary.connectorRowCount, 756, "matrix template should cover all matrix connectors");
assert.strictEqual(browserMatrixTemplate.rows[0].caseId, "FM_000", "matrix template should use stable matrix case ids");
assert.strictEqual(browserMatrixTemplate.rows[0].observedRightCellZ, "", "matrix template observed fields should be blank");
assert.ok(RAD.exportTwoCellFidelityMatrixMeasurementTemplateCsv(twoCellBenchState).includes("observedVerticalSlipMm"), "matrix template CSV should expose observed slip fields");
const browserExternalManifest = RAD.twoCellExternalFidelityMatrixManifest(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(browserExternalManifest.schema, "rad-sim.browser-two-cell-external-fidelity-matrix-manifest.v1", "external fidelity manifest should expose stable schema");
assert.strictEqual(browserExternalManifest.summary.caseCount, 252, "external manifest should cover every matrix case");
assert.strictEqual(browserExternalManifest.summary.connectorMeasurementRowCount, 756, "external manifest should include fillable connector rows");
assert.strictEqual(browserExternalManifest.summary.physicalAccuracyValidated, false, "external manifest should not claim physical validation");
assert.ok(browserExternalManifest.summary.missingEvidence.includes("externalEngineRunResults"), "external manifest should require engine results");
assert.ok(RAD.exportTwoCellExternalFidelityMatrixManifestCsv(twoCellBenchState).includes("FM_000"), "external manifest CSV should include stable case IDs");
assert.ok(RAD.exportTwoCellExternalFidelityMatrixManifestJson(twoCellBenchState).includes("external_fidelity_matrix_results"), "external manifest JSON should include result paths");
const parsedBlankMatrixMeasurements = RAD.twoCellFidelityMatrixMeasurementsFromCsv(
  RAD.exportTwoCellFidelityMatrixMeasurementTemplateCsv(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 })
);
assert.strictEqual(parsedBlankMatrixMeasurements.rows.length, 756, "browser should parse the full blank matrix template CSV");
const blankMatrixComparison = RAD.compareTwoCellFidelityMatrixMeasurements(twoCellBenchState, parsedBlankMatrixMeasurements, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(blankMatrixComparison.schema, "rad-sim.browser-two-cell-fidelity-matrix-measurement-comparison.v1", "matrix comparison should expose stable schema");
assert.strictEqual(blankMatrixComparison.summary.requiresMoreData, true, "blank matrix measurements should require more data");
assert.strictEqual(blankMatrixComparison.summary.observedScalarCount, 0, "blank matrix measurements should not fake observed values");
const filledMatrixRows = browserMatrixTemplate.rows.map((row) => {
  const filled = { ...row };
  for (const [observed, predicted] of Object.entries({
    observedLeftCellX: "predictedLeftCellX",
    observedLeftCellY: "predictedLeftCellY",
    observedLeftCellZ: "predictedLeftCellZ",
    observedRightCellX: "predictedRightCellX",
    observedRightCellY: "predictedRightCellY",
    observedRightCellZ: "predictedRightCellZ",
    observedLeftAlpha: "predictedLeftAlpha",
    observedRightAlpha: "predictedRightAlpha",
    observedLeftTheta: "predictedLeftTheta",
    observedRightTheta: "predictedRightTheta",
    observedLeftXmm: "predictedLeftXmm",
    observedLeftYmm: "predictedLeftYmm",
    observedLeftZmm: "predictedLeftZmm",
    observedRightXmm: "predictedRightXmm",
    observedRightYmm: "predictedRightYmm",
    observedRightZmm: "predictedRightZmm",
    observedLateralSlipMm: "predictedLateralSlipMm",
    observedVerticalSlipMm: "predictedVerticalSlipMm",
    observedTotalSlipMm: "predictedTotalSlipMm",
  })) {
    filled[observed] = filled[predicted];
  }
  filled.observedContactMode = filled.predictedContactMode;
  filled.lockHeldObserved = String(filled.lockMode || "").includes("locked") ? "held" : "released";
  return filled;
});
const filledMatrixComparison = RAD.compareTwoCellFidelityMatrixMeasurements(twoCellBenchState, filledMatrixRows, { alphaCommand: -0.42, zCommand: 0.4, tolerance: 1e-9 });
assert.strictEqual(filledMatrixComparison.summary.passesTolerance, true, "synthetic filled matrix measurements should roundtrip at zero residual");
assert.ok(filledMatrixComparison.summary.observedScalarCount > 10000, "filled matrix comparison should count observed scalar fields");
assert.strictEqual(filledMatrixComparison.summary.maxAbsError, 0, "filled matrix comparison should report zero max residual");
assert.ok(RAD.exportTwoCellFidelityMatrixMeasurementComparisonJson(twoCellBenchState, filledMatrixRows).includes("observedScalarCount"), "matrix comparison export should include residual summary");
const filledMatrixCalibration = RAD.calibrateTwoCellFidelityMatrixParameters(twoCellBenchState, filledMatrixRows, { alphaCommand: -0.42, zCommand: 0.4, tolerance: 1e-9 });
assert.strictEqual(filledMatrixCalibration.schema, "rad-sim.browser-two-cell-fidelity-matrix-parameter-calibration.v1", "matrix parameter calibration should expose stable schema");
assert.strictEqual(filledMatrixCalibration.summary.readyForReducedProxyCalibration, true, "filled matrix should be ready for reduced-proxy calibration");
assert.strictEqual(filledMatrixCalibration.summary.physicalAccuracyValidated, false, "calibration should not claim exact physical validation");
assert.strictEqual(Number(filledMatrixCalibration.estimates.zResponseScale.estimate.toFixed(8)), 1, "synthetic filled matrix should estimate z scale one");
assert.ok(RAD.exportTwoCellFidelityMatrixParameterCalibrationJson(twoCellBenchState, filledMatrixRows).includes("proposedReducedProxyUpdates"), "calibration export should include proposed proxy updates");
const verticallyScaledRows = filledMatrixRows.map((row) => {
  const scaled = { ...row };
  for (const field of ["observedLeftCellZ", "observedRightCellZ", "observedLeftZmm", "observedRightZmm", "observedVerticalSlipMm", "observedTotalSlipMm"]) {
    scaled[field] = Number(scaled[field]) * 1.5;
  }
  return scaled;
});
const scaledMatrixCalibration = RAD.calibrateTwoCellFidelityMatrixParameters(twoCellBenchState, verticallyScaledRows, { alphaCommand: -0.42, zCommand: 0.4 });
assert.ok(scaledMatrixCalibration.estimates.zResponseScale.estimate > 1.2, "larger measured z should increase z response scale");
assert.ok(scaledMatrixCalibration.estimates.verticalSlipScale.estimate > 1.2, "larger measured connector slip should increase vertical slip scale");
assert.ok(scaledMatrixCalibration.proposedReducedProxyUpdates.zCouplingGain > twoCellBenchState.grid.zCouplingGain, "larger measured z should propose stronger z coupling");
assert.ok(scaledMatrixCalibration.proposedReducedProxyUpdates.holeRadius > twoCellBenchState.grid.holeRadius, "larger measured slip should propose larger effective hole radius");
const measuredFidelityStatus = RAD.twoCellPhysicalFidelityStatus(twoCellBenchState, {}, { suite: browserSuite, matrix: browserMatrix, comparison: filledMatrixComparison, calibration: filledMatrixCalibration });
assert.strictEqual(measuredFidelityStatus.summary.measurementComparisonReady, true, "status should detect loaded matrix measurements");
assert.strictEqual(measuredFidelityStatus.summary.parameterCalibrationReady, true, "status should detect calibrated matrix measurements");
assert.strictEqual(measuredFidelityStatus.summary.reducedProxyMeasurementValidated, true, "status should mark the reduced proxy validated by matching measurements");
assert.strictEqual(measuredFidelityStatus.summary.physicalAccuracyValidated, false, "status should still not claim exact physical validation");
const connectorTemplate = RAD.twoCellConnectorMeasurementTemplate(twoCellBenchState, { alphaCommand: -0.42, zCommand: 0.4 });
assert.strictEqual(connectorTemplate.schema, "rad-sim.two-cell-connector-measurement-template.v1", "connector template should expose stable schema");
assert.strictEqual(connectorTemplate.summary.connectorRowCount, 36, "connector template should cover 12 cases x 3 connectors");
assert.deepStrictEqual(new Set(connectorTemplate.rows.map((row) => row.connector)), new Set(["upper", "middle", "lower"]), "connector template should cover all connector names");
assert.strictEqual(connectorTemplate.rows[0].observedVerticalSlipMm, "", "connector template observed fields should be fillable blanks");
assert.ok(RAD.exportTwoCellConnectorMeasurementTemplateCsv(twoCellBenchState).includes("observedRightZmm"), "connector template CSV should expose observed marker fields");
assert.strictEqual(typeof RAD.lockDatasetSummary, "function", "lock surrogate should expose dataset summary");
assert.strictEqual(typeof RAD.predictEmpiricalLockCoordinates, "function", "lock surrogate should expose coordinate prediction");
assert.strictEqual(typeof RAD.applyEmpiricalLockSurrogate, "function", "lock surrogate should expose simulation postprocessor");
assert.strictEqual(typeof RAD.Primitives.runSelfCheck, "function", "primitive library should expose browser self-check");
const primitiveSelfCheck = RAD.Primitives.runSelfCheck({ renderer: { getDiagnostics: () => ({ drawCalls: 1 }) } });
assert.strictEqual(primitiveSelfCheck.schema, "rad-sim.browser-self-check.v1");
assert.strictEqual(primitiveSelfCheck.ok, true, primitiveSelfCheck.checks.filter((check) => !check.ok).map((check) => `${check.name}: ${check.detail}`).join("; "));
assert.strictEqual(context.window.RAD_LOCK_DATASET.schema, "rad-sim.browser-lock-coordinate-dataset.v1", "RAD_LOCK_DATASET should load processed measured coordinates");
assert.strictEqual(RAD.lockDatasetSummary().configurationCount, 37, "lock dataset should expose measured configurations");
assert.strictEqual(typeof RAD.physicalPreviewComparison, "function", "analysis module should expose physical preview comparison");
assert.strictEqual(typeof RAD.responseDecayProfile, "function", "analysis module should expose response decay profile");
assert.strictEqual(typeof RAD.calibrationExperimentProtocol, "function", "analysis module should expose calibration experiment protocol");
assert.strictEqual(typeof RAD.exportCalibrationExperimentProtocol, "function", "analysis module should export calibration experiment protocol");
assert.strictEqual(typeof RAD.responseAtlas, "function", "analysis module should expose response atlas");
assert.strictEqual(typeof RAD.exportResponseAtlas, "function", "analysis module should export response atlas");
assert.strictEqual(typeof RAD.responseAtlasSweep, "function", "analysis module should expose response atlas sweeps");
assert.strictEqual(typeof RAD.exportResponseAtlasSweep, "function", "analysis module should export response atlas sweeps");
assert.strictEqual(typeof RAD.springHingePhysicalPreviewReportDescriptor, "function", "analysis module should describe Python spring-hinge physical preview reports");
assert.strictEqual(typeof RAD.buildResponseMatrix, "function", "analysis module should expose response matrix export");
assert.strictEqual(typeof RAD.exportResponseMatrix, "function", "analysis module should serialize response matrices");
assert.strictEqual(typeof RAD.topologyExperimentReport, "function", "analysis module should expose topology experiment reports");
assert.strictEqual(typeof RAD.exportTopologyExperimentReport, "function", "analysis module should serialize topology experiment reports");
assert.strictEqual(typeof RAD.cellRemoved, "function", "math module should expose removed-cell topology checks");
assert.strictEqual(typeof RAD.topologyDiagnostics, "function", "math module should expose topology diagnostics");
assert.strictEqual(typeof RAD.groupActuationEvent, "function", "operator module should expose group actuation events");
assert.strictEqual(typeof RAD.compareGroupActuationDecomposition, "function", "operator module should expose group decomposition diagnostics");
assert.strictEqual(typeof RAD.compareGroupActuationUnderRemoval, "function", "operator module should expose group removal diagnostics");
assert.strictEqual(typeof RAD.compareVerticalResidualUnderRemoval, "function", "operator module should expose vertical removal diagnostics");
assert.strictEqual(typeof RAD.removeCellEvent, "function", "operator module should expose remove-cell events");
assert.strictEqual(typeof RAD.restoreCellEvent, "function", "operator module should expose restore-cell events");
assert.strictEqual(typeof RAD.programmableDiscontinuityReport, "function", "analysis module should expose programmable-discontinuity reports");
assert.strictEqual(typeof RAD.exportProgrammableDiscontinuityReport, "function", "analysis module should export programmable-discontinuity reports");
assert.strictEqual(typeof RAD.formalizationTargetManifest, "function", "analysis module should expose formalization target manifests");
assert.strictEqual(typeof RAD.exportFormalizationTargetManifest, "function", "analysis module should export formalization target manifests");
assert.strictEqual(typeof RAD.calibrationExperimentResultsTemplate, "function", "analysis module should expose calibration results template");
assert.strictEqual(typeof RAD.exportCalibrationExperimentResultsTemplate, "function", "analysis module should export calibration results template");
assert.strictEqual(typeof RAD.calibrationBenchNotebook, "function", "analysis module should expose calibration bench notebooks");
assert.strictEqual(typeof RAD.exportCalibrationBenchNotebook, "function", "analysis module should export calibration bench notebooks");
assert.strictEqual(typeof RAD.exportCalibrationBenchNotebookCsv, "function", "analysis module should export calibration bench notebook CSV tables");
assert.strictEqual(typeof RAD.calibrationBenchPacket, "function", "analysis module should expose calibration bench packets");
assert.strictEqual(typeof RAD.exportCalibrationBenchPacket, "function", "analysis module should export calibration bench packets");
assert.strictEqual(typeof RAD.calibrationBenchExecutionValidation, "function", "analysis module should expose executed calibration bench validations");
assert.strictEqual(typeof RAD.exportCalibrationBenchExecutionValidation, "function", "analysis module should export executed calibration bench validations");
assert.strictEqual(typeof RAD.exportCalibrationBenchExecutionValidationCsv, "function", "analysis module should export executed calibration bench validation CSV tables");
assert.strictEqual(typeof RAD.physicalValidationReadinessReport, "function", "analysis module should expose physical validation readiness reports");
assert.strictEqual(typeof RAD.exportPhysicalValidationReadiness, "function", "analysis module should export physical validation readiness reports");
assert.strictEqual(typeof RAD.exportPhysicalValidationReadinessCsv, "function", "analysis module should export physical validation readiness CSV tables");
assert.strictEqual(typeof RAD.contactStateAbstractionReport, "function", "analysis module should expose contact-state abstraction reports");
assert.strictEqual(typeof RAD.exportContactStateAbstraction, "function", "analysis module should export contact-state abstraction reports");
assert.strictEqual(typeof RAD.exportContactStateAbstractionCsv, "function", "analysis module should export contact-state abstraction CSV tables");
assert.strictEqual(typeof RAD.contactGraphConsistencyReport, "function", "analysis module should expose contact graph consistency reports");
assert.strictEqual(typeof RAD.exportContactGraphConsistency, "function", "analysis module should export contact graph consistency reports");
assert.strictEqual(typeof RAD.exportContactGraphConsistencyCsv, "function", "analysis module should export contact graph consistency CSV tables");
assert.strictEqual(typeof RAD.physicalRealizationMapReport, "function", "analysis module should expose physical realization map reports");
assert.strictEqual(typeof RAD.exportPhysicalRealizationMap, "function", "analysis module should export physical realization map reports");
assert.strictEqual(typeof RAD.exportPhysicalRealizationMapCsv, "function", "analysis module should export physical realization map CSV tables");
assert.strictEqual(typeof RAD.externalPhysicsEngineAuditReport, "function", "analysis module should expose external physics engine audits");
assert.strictEqual(typeof RAD.exportExternalPhysicsEngineAudit, "function", "analysis module should export external physics engine audits");
assert.strictEqual(typeof RAD.exportExternalPhysicsEngineAuditCsv, "function", "analysis module should export external physics engine audit CSV tables");
assert.strictEqual(typeof RAD.mujocoModelExportReport, "function", "analysis module should expose MuJoCo model exports");
assert.strictEqual(typeof RAD.exportMujocoModelXml, "function", "analysis module should export MuJoCo XML");
assert.strictEqual(typeof RAD.exportMujocoModelReport, "function", "analysis module should export MuJoCo model reports");
assert.strictEqual(typeof RAD.mujocoPinHoleContactGeometryReport, "function", "analysis module should expose MuJoCo pin-hole contact geometry");
assert.strictEqual(typeof RAD.exportMujocoPinHoleContactGeometry, "function", "analysis module should export MuJoCo pin-hole contact geometry");
assert.strictEqual(typeof RAD.exportMujocoPinHoleContactGeometryCsv, "function", "analysis module should export MuJoCo pin-hole contact geometry CSV tables");
assert.strictEqual(typeof RAD.mujocoContactParameterReport, "function", "analysis module should expose MuJoCo contact parameter profiles");
assert.strictEqual(typeof RAD.exportMujocoContactParameter, "function", "analysis module should export MuJoCo contact parameter profiles");
assert.strictEqual(typeof RAD.exportMujocoContactParameterCsv, "function", "analysis module should export MuJoCo contact parameter CSV tables");
assert.strictEqual(typeof RAD.contactParameterCalibrationPacket, "function", "analysis module should expose contact parameter calibration packets");
assert.strictEqual(typeof RAD.exportContactParameterCalibrationPacket, "function", "analysis module should export contact parameter calibration packets");
assert.strictEqual(typeof RAD.exportContactParameterCalibrationPacketCsv, "function", "analysis module should export contact parameter calibration CSV tables");
assert.strictEqual(typeof RAD.contactParameterCalibrationResultsTemplate, "function", "analysis module should expose contact parameter calibration result templates");
assert.strictEqual(typeof RAD.exportContactParameterCalibrationResultsTemplate, "function", "analysis module should export contact parameter calibration result templates");
assert.strictEqual(typeof RAD.contactParameterCalibrationResultsFromJson, "function", "analysis module should parse contact parameter calibration result artifacts");
assert.strictEqual(typeof RAD.compareContactParameterCalibrationResults, "function", "analysis module should compare contact parameter calibration results");
assert.strictEqual(typeof RAD.exportContactParameterBenchValidation, "function", "analysis module should export contact parameter bench validations");
assert.strictEqual(typeof RAD.exportContactParameterBenchValidationCsv, "function", "analysis module should export contact parameter bench validation CSV tables");
assert.strictEqual(typeof RAD.contactParameterIntervalCalibrationReport, "function", "analysis module should expose contact parameter interval calibrations");
assert.strictEqual(typeof RAD.exportContactParameterIntervalCalibration, "function", "analysis module should export contact parameter interval calibrations");
assert.strictEqual(typeof RAD.exportContactParameterIntervalCalibrationCsv, "function", "analysis module should export contact parameter interval calibration CSV tables");
assert.strictEqual(typeof RAD.mujocoExternalRunReport, "function", "analysis module should expose MuJoCo run reports");
assert.strictEqual(typeof RAD.exportMujocoExternalRun, "function", "analysis module should export MuJoCo run reports");
assert.strictEqual(typeof RAD.mujocoExternalComparisonReport, "function", "analysis module should expose MuJoCo comparison reports");
assert.strictEqual(typeof RAD.exportMujocoExternalComparison, "function", "analysis module should export MuJoCo comparison reports");
assert.strictEqual(typeof RAD.exportMujocoExternalComparisonCsv, "function", "analysis module should export MuJoCo comparison CSV tables");
assert.strictEqual(typeof RAD.equilibriumRelationReport, "function", "analysis module should expose equilibrium relation reports");
assert.strictEqual(typeof RAD.exportEquilibriumRelation, "function", "analysis module should export equilibrium relation reports");
assert.strictEqual(typeof RAD.exportEquilibriumRelationCsv, "function", "analysis module should export equilibrium relation CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumControllabilityReport, "function", "analysis module should expose reachable equilibrium controllability reports");
assert.strictEqual(typeof RAD.exportReachableEquilibriumControllability, "function", "analysis module should export reachable equilibrium controllability reports");
assert.strictEqual(typeof RAD.exportReachableEquilibriumControllabilityCsv, "function", "analysis module should export reachable equilibrium controllability CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumBenchProtocol, "function", "analysis module should expose reachable equilibrium bench protocols");
assert.strictEqual(typeof RAD.exportReachableEquilibriumBenchProtocol, "function", "analysis module should export reachable equilibrium bench protocols");
assert.strictEqual(typeof RAD.exportReachableEquilibriumBenchProtocolCsv, "function", "analysis module should export reachable equilibrium bench protocol CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumBenchResultsTemplate, "function", "analysis module should expose reachable equilibrium bench result templates");
assert.strictEqual(typeof RAD.exportReachableEquilibriumBenchResultsTemplate, "function", "analysis module should export reachable equilibrium bench result templates");
assert.strictEqual(typeof RAD.reachableEquilibriumBenchResultsFromJson, "function", "analysis module should parse reachable equilibrium bench results");
assert.strictEqual(typeof RAD.compareReachableEquilibriumBenchResults, "function", "analysis module should compare reachable equilibrium bench results");
assert.strictEqual(typeof RAD.exportReachableEquilibriumBenchComparison, "function", "analysis module should export reachable equilibrium bench comparisons");
assert.strictEqual(typeof RAD.exportReachableEquilibriumBenchComparisonCsv, "function", "analysis module should export reachable equilibrium bench comparison CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumAmplitudeCalibrationReport, "function", "analysis module should expose reachable equilibrium amplitude calibration reports");
assert.strictEqual(typeof RAD.exportReachableEquilibriumAmplitudeCalibration, "function", "analysis module should export reachable equilibrium amplitude calibration reports");
assert.strictEqual(typeof RAD.exportReachableEquilibriumAmplitudeCalibrationCsv, "function", "analysis module should export reachable equilibrium amplitude calibration CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumEmpiricalProfileFromAmplitude, "function", "analysis module should expose reachable equilibrium empirical profiles");
assert.strictEqual(typeof RAD.exportReachableEquilibriumEmpiricalProfile, "function", "analysis module should export reachable equilibrium empirical profiles");
assert.strictEqual(typeof RAD.exportReachableEquilibriumEmpiricalProfileCsv, "function", "analysis module should export reachable equilibrium empirical profile CSV tables");
assert.strictEqual(typeof RAD.compareCalibrationExperimentResults, "function", "analysis module should compare calibration results");
assert.strictEqual(typeof RAD.summarizeCalibrationComparison, "function", "analysis module should summarize calibration result comparisons");
assert.strictEqual(typeof RAD.calibrationParameterEstimates, "function", "analysis module should expose calibration parameter estimates");
assert.strictEqual(typeof RAD.calibrationModelProfile, "function", "analysis module should expose calibration model profiles");
assert.strictEqual(typeof RAD.applyCalibrationModelProfile, "function", "analysis module should apply safe calibration model profile updates");
assert.strictEqual(typeof RAD.exportCalibrationModelProfile, "function", "analysis module should export calibration model profiles");
assert.strictEqual(typeof RAD.calibrationModelProfileResidualComparison, "function", "analysis module should compare calibration model-profile residuals");
assert.strictEqual(typeof RAD.calibrationModelProfileHoldoutValidation, "function", "analysis module should validate model profiles on holdout calibration results");
assert.strictEqual(typeof RAD.selectCalibrationModelProfile, "function", "analysis module should select calibration model-profile candidates");
assert.strictEqual(typeof RAD.exportCalibrationModelProfileResidualComparison, "function", "analysis module should export calibration model-profile residual comparisons");
assert.strictEqual(typeof RAD.exportCalibrationModelProfileHoldoutValidation, "function", "analysis module should export calibration model-profile holdout validations");
assert.strictEqual(typeof RAD.exportCalibrationModelProfileHoldoutValidationCsv, "function", "analysis module should export calibration model-profile holdout validation CSV tables");
assert.strictEqual(typeof RAD.exportCalibrationModelProfileSelection, "function", "analysis module should export calibration model-profile selections");
assert.strictEqual(typeof RAD.calibrationComparisonReport, "function", "analysis module should expose calibration comparison reports");
assert.strictEqual(typeof RAD.exportCalibrationComparisonReport, "function", "analysis module should export calibration comparison reports");
assert.strictEqual(typeof RAD.validateInversePlanPhysical, "function", "inverse module should expose physical inverse validation");
assert.strictEqual(typeof RAD.inverseDesignReport, "function", "inverse module should expose inverse design reports");
assert.strictEqual(typeof RAD.exportInverseDesignReport, "function", "inverse module should export inverse design reports");
assert.strictEqual(typeof RAD.reachableEquilibriumProfileInverseReport, "function", "inverse module should expose profile-aware inverse diagnostics");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInverse, "function", "inverse module should export profile-aware inverse diagnostics");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInverseCsv, "function", "inverse module should export profile-aware inverse CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumProfileInverseAcceptanceReport, "function", "inverse module should expose profile-aware inverse acceptance gates");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInverseAcceptance, "function", "inverse module should export profile-aware inverse acceptance gates");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInverseAcceptanceCsv, "function", "inverse module should export profile-aware inverse acceptance CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumProfileInversePreviewPacket, "function", "inverse module should expose profile-aware inverse preview packets");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInversePreviewPacket, "function", "inverse module should export profile-aware inverse preview packets");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInversePreviewPacketCsv, "function", "inverse module should export profile-aware inverse preview-packet CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumProfileInversePreviewReplayReport, "function", "inverse module should expose profile-aware inverse preview replay reports");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInversePreviewReplay, "function", "inverse module should export profile-aware inverse preview replay reports");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInversePreviewReplayCsv, "function", "inverse module should export profile-aware inverse preview replay CSV tables");
assert.strictEqual(typeof RAD.reachableEquilibriumProfileInversePreviewPhysicalReport, "function", "inverse module should expose profile-aware inverse preview physical reports");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInversePreviewPhysical, "function", "inverse module should export profile-aware inverse preview physical reports");
assert.strictEqual(typeof RAD.exportReachableEquilibriumProfileInversePreviewPhysicalCsv, "function", "inverse module should export profile-aware inverse preview physical CSV tables");
const provenance = RAD.modelProvenance();
assert.ok(
  provenance.some((item) => item.id === "backlash_dead_zone" && item.status === "paper-supported"),
  "backlash dead-zone should be marked paper-supported"
);
assert.ok(provenance.some((item) => item.status === "implementation-assumption"), "provenance should expose implementation assumptions");
assert.ok(provenance.some((item) => item.status === "calibration-gap"), "provenance should expose calibration gaps");
const provenanceSummary = RAD.provenanceSummary();
assert.ok(provenanceSummary["paper-supported"] >= 1, "provenance summary should count paper-supported items");

const html = fs.readFileSync(path.join(web, "index.html"), "utf8");
const localRefs = Array.from(html.matchAll(/(?:src|href)="(\.\/[^"]+)"/g)).map((match) => match[1]);
assert.ok(localRefs.includes("./vendor/three.min.js"), "local Three.js asset should be referenced");
assert.ok(localRefs.includes("./primitives.js"), "local primitive library should be referenced");
assert.ok(localRefs.includes("./data/lock_dataset.js"), "processed measured lock dataset should be referenced");
assert.ok(localRefs.includes("./lock_surrogate.js"), "browser lock surrogate should be referenced");
assert.ok(localRefs.includes("./two_cell_bench.js"), "browser two-cell bench module should be referenced");
for (const ref of localRefs) {
  assert.ok(fs.existsSync(path.join(web, ref.replace("./", ""))), `missing local browser asset: ${ref}`);
}
assert.ok(!/https?:\/\//.test(html), "browser entry should not require external scripts or styles");
assert.ok(html.includes('value="operatorInteraction"'), "operator interaction overlay should be available in the browser UI");
assert.ok(html.includes('value="calibrationError"'), "calibration error overlay should be available in the browser UI");
assert.ok(html.includes('value="calibrationResidual"'), "calibration residual overlay should be available in the browser UI");
assert.ok(html.includes('value="underactuated"'), "underactuated target overlay should be available in the browser UI");
assert.ok(html.includes('value="topologyBlocked"'), "topology-blocked target overlay should be available in the browser UI");
assert.ok(html.includes('value="graph"'), "topology graph display mode should be available in the browser UI");
assert.ok(html.includes('value="sheetOnly"'), "sheet-only contour display mode should be available in the browser UI");
assert.ok(html.includes('id="modeSheet"'), "sheet workspace mode should be available from the viewport toolbar");
assert.ok(html.includes('data-workspace-mode="sheet"'), "sheet workspace mode should be first-class, not only a dropdown value");
assert.ok(html.includes('id="selectInteractionHotspot"'), "response panel should expose a hotspot selection button");
assert.ok(html.includes('id="saveResponseMatrix"'), "response panel should expose a response-matrix export button");
assert.ok(html.includes('id="saveTopologyReport"'), "response panel should expose a topology report export button");
assert.ok(html.includes('id="saveProgrammableReport"'), "response panel should expose a programmable-discontinuity report export button");
assert.ok(html.includes('id="saveFormalizationTargets"'), "response panel should expose a formalization-target export button");
assert.ok(html.includes('id="selectCalibrationHotspot"'), "response panel should expose a calibration-error selection button");
assert.ok(html.includes('id="nextCalibrationHotspot"'), "response panel should expose ranked calibration-error navigation");
assert.ok(html.includes('id="selectUnderactuatedTarget"'), "inverse panel should expose underactuated target selection");
assert.ok(html.includes('id="saveInverseReport"'), "inverse panel should expose inverse report export");
assert.ok(html.includes('id="buildPreviewPacket"'), "inverse panel should expose profile-aware preview-packet build control");
assert.ok(html.includes('id="runPreviewReplay"'), "inverse panel should expose preview-packet replay control");
assert.ok(html.includes('id="runPreviewPhysical"'), "inverse panel should expose preview-packet physical preview control");
assert.ok(html.includes('id="savePreviewPacket"'), "inverse panel should expose preview-packet export");
assert.ok(html.includes('id="savePreviewReplay"'), "inverse panel should expose preview-replay export");
assert.ok(html.includes('id="savePreviewPhysical"'), "inverse panel should expose preview-physical export");
assert.ok(html.includes('id="underTargetCells"'), "metric strip should count underactuated targets");
assert.ok(html.includes('id="inverseReachabilitySummary"'), "inverse plan readout should summarize underactuated targets");
assert.ok(html.includes('id="inversePacketSummary"'), "inverse plan readout should summarize preview-packet readiness");
assert.ok(html.includes('id="inversePacketReplaySummary"'), "inverse plan readout should summarize preview replay readiness");
assert.ok(html.includes('id="inversePacketPhysicalSummary"'), "inverse plan readout should summarize preview physical readiness");
assert.ok(html.includes('id="saveExperimentProtocol"'), "response panel should expose a protocol export button");
assert.ok(html.includes('id="saveResponseAtlas"'), "response panel should expose a response-atlas export button");
assert.ok(html.includes('id="runResponseAtlasSweep"'), "response panel should expose a sweep run button");
assert.ok(html.includes('id="saveResponseAtlasSweep"'), "response panel should expose a response-atlas sweep export button");
assert.ok(html.includes('id="saveResultsTemplate"'), "response panel should expose a results-template export button");
assert.ok(html.includes('id="saveBenchPacket"'), "response panel should expose a bench-packet export button");
assert.ok(html.includes('id="saveBenchNotebook"'), "response panel should expose a bench-notebook export button");
assert.ok(html.includes('id="saveBenchNotebookCsv"'), "response panel should expose a bench-notebook CSV export button");
assert.ok(html.includes('id="loadResultsJson"'), "response panel should expose a results import button");
assert.ok(html.includes('id="loadHoldoutJson"'), "response panel should expose a holdout results import button");
assert.ok(html.includes('id="saveComparisonReport"'), "response panel should expose a comparison-report export button");
assert.ok(html.includes('id="saveModelProfile"'), "response panel should expose a calibration model-profile export button");
assert.ok(html.includes('id="compareModelProfile"'), "response panel should expose a calibration model-profile residual-check button");
assert.ok(html.includes('id="saveModelProfileComparison"'), "response panel should expose a calibration model-profile residual-check export button");
assert.ok(html.includes('id="checkModelProfileHoldout"'), "response panel should expose a holdout validation button");
assert.ok(html.includes('id="saveModelProfileHoldout"'), "response panel should expose a holdout validation export button");
assert.ok(html.includes('id="saveModelProfileHoldoutCsv"'), "response panel should expose a holdout validation CSV export button");
assert.ok(html.includes('id="applyModelProfile"'), "response panel should expose a safe model-profile apply button");
assert.ok(html.includes('id="calibrationResultsFileInput"'), "response panel should include hidden calibration results file input");
assert.ok(html.includes('id="calibrationHoldoutFileInput"'), "response panel should include hidden holdout results file input");
assert.ok(html.includes('id="calibrationResultsSummary"'), "response panel should expose calibration result summary");
assert.ok(html.includes('id="calibrationResultsError"'), "response panel should expose calibration result error readout");
assert.ok(html.includes('id="modelProfileSummary"'), "response panel should expose calibration model-profile residual summary");
assert.ok(html.includes('id="modelProfileAudit"'), "response panel should expose calibration model-profile audit summary");
assert.ok(html.includes('id="modelProfileHoldoutSummary"'), "response panel should expose holdout validation summary");
assert.ok(html.includes('id="modelProfileHoldoutDetail"'), "response panel should expose holdout validation detail");
assert.ok(html.includes('id="frameworkLawSummary"'), "response panel should expose framework law summary");
assert.ok(html.includes('id="frameworkLawDetail"'), "response panel should expose framework law detail");
assert.ok(html.includes('id="formalizationSummary"'), "response panel should expose formalization summary");
assert.ok(html.includes('id="formalizationDetail"'), "response panel should expose formalization detail");
assert.ok(html.includes('id="sweepSummary"'), "response panel should expose sweep summary readout");
assert.ok(html.includes('id="sweepTrend"'), "response panel should expose sweep trend readout");
assert.ok(html.includes('id="characterizationHotspot"'), "response panel should expose a hotspot readout");
assert.ok(html.includes("./provenance.js"), "provenance module should be loaded by the browser entry");
assert.ok(html.includes('id="modelProvenanceList"'), "browser UI should expose model provenance list");
assert.ok(html.includes('id="provenancePaper"'), "browser UI should expose paper-supported provenance count");
assert.ok(html.includes('id="hardwareCoverageOut"'), "browser UI should expose hardware profile coverage");
assert.ok(html.includes('id="hardwareMissingOut"'), "browser UI should expose missing calibration fields");
assert.ok(html.includes('id="hardwareReadinessOut"'), "browser UI should expose calibration readiness");
assert.ok(html.includes('id="hardwareSolverGapOut"'), "browser UI should expose solver calibration gaps");
assert.ok(html.includes('id="hardwareMeasurementPlanOut"'), "browser UI should expose next calibration measurements");
assert.ok(html.includes('id="hardwareProfileName"'), "browser UI should expose measured profile name input");
assert.ok(html.includes('id="hardwarePinRadiusMm"'), "browser UI should expose measured pin-radius input");
assert.ok(html.includes('id="hardwareHoleRadiusMm"'), "browser UI should expose measured hole-radius input");
assert.ok(html.includes('id="applyHardwareProfile"'), "browser UI should expose hardware profile apply action");
assert.ok(html.includes('id="saveHardwareProfile"'), "browser UI should expose hardware profile export action");
assert.ok(html.includes('id="loadHardwareProfile"'), "browser UI should expose hardware profile import action");
assert.ok(html.includes('id="hardwareProfileFileInput"'), "browser UI should include hidden hardware profile file input");
assert.ok(html.includes('id="saveCalibrationPlan"'), "browser UI should expose calibration plan export action");
assert.ok(html.includes('value="calibratedRad"'), "browser UI should expose calibrated RAD visual mode");
assert.ok(html.includes('value="hardware"'), "browser UI should expose hardware profile measurement mode");
assert.ok(html.includes('id="empiricalLockModel"'), "browser UI should expose measured lock-data toggle");
assert.ok(html.includes('id="lockStateIndex"'), "browser UI should expose stable-state selection");
assert.ok(html.includes('id="runSystemCheck"'), "browser UI should expose primitive self-check control near the viewport");
assert.ok(html.includes('id="systemCheckPanel"'), "browser UI should expose primitive self-check result panel");
assert.ok(html.includes('value="cadRad"'), "browser UI should expose CAD RAD cell visual mode");
assert.ok(html.includes('id="runTwoCellBench"'), "browser UI should expose two-cell bench run action");
assert.ok(html.includes('id="applyTwoCellBench"'), "browser UI should expose 1x2 two-cell bench action");
assert.ok(html.includes('id="twoCellExternalCase"'), "browser UI should expose external two-cell proxy case selection");
assert.ok(html.includes('id="applyTwoCellExternalCase"'), "browser UI should expose external proxy case application");
assert.ok(html.includes('id="twoCellPhaseDiagramGrid"'), "browser UI should expose the clickable radius/backlash phase diagram grid");
assert.ok(html.includes('id="runTwoCellTransitionReport"'), "browser UI should expose transition-report run action");
assert.ok(html.includes('id="saveTwoCellTransitionReportCsv"'), "browser UI should expose transition-report CSV export");
assert.ok(html.includes('id="loadTwoCellTransitionMeasurements"'), "browser UI should expose transition measurement import");
assert.ok(html.includes('id="twoCellTransitionComparisonSummary"'), "browser UI should summarize transition measurement comparison");
assert.ok(html.includes('id="saveTwoCellTransitionRerun"'), "browser UI should expose calibrated transition rerun export");
assert.ok(html.includes('id="twoCellTransitionRerunSummary"'), "browser UI should summarize calibrated transition reruns");
assert.ok(html.includes('id="applyTwoCellTransitionCalibration"'), "browser UI should expose explicit transition-calibration apply");
assert.ok(html.includes('id="saveTwoCellExactContactHandoffPlan"'), "browser UI should expose exact-contact handoff export");
assert.ok(html.includes('id="saveTwoCellExactContactHandoffPlanCsv"'), "browser UI should expose exact-contact handoff CSV export");
assert.ok(html.includes('id="twoCellExactContactPlanSummary"'), "browser UI should summarize exact-contact handoff readiness");
const scriptOrder = [
  "./vendor/three.min.js",
  "./state.js",
  "./provenance.js",
  "./math.js",
  "./operators.js",
  "./inverse.js",
  "./analysis.js",
  "./data/lock_dataset.js",
  "./lock_surrogate.js",
  "./physics.js",
  "./two_cell_bench.js",
  "./mesh_export.js",
  "./renderer.js",
  "./ui.js",
  "./primitives.js",
  "./app.js",
];
for (let i = 1; i < scriptOrder.length; i += 1) {
  assert.ok(html.indexOf(scriptOrder[i - 1]) < html.indexOf(scriptOrder[i]), `${scriptOrder[i - 1]} should load before ${scriptOrder[i]}`);
}
for (const filename of ["renderer.js", "ui.js", "primitives.js", "app.js"]) {
  const source = fs.readFileSync(path.join(web, filename), "utf8");
  assert.doesNotThrow(() => new Function(source), `${filename} should compile`);
  if (filename === "ui.js") {
    assert.ok(source.includes("selectCalibrationHotspot(1)"), "UI should cycle through ranked calibration hotspots");
    assert.ok(source.includes("fitResidualTopCells"), "UI should support ranked residual calibration hotspots");
    assert.ok(source.includes("RAD.calibrationParameterEstimates"), "UI should surface calibration parameter estimates");
    assert.ok(source.includes("contactStatus"), "UI calibration readout should expose contact proxy status");
    assert.ok(source.includes("saveModelProfile"), "UI should save calibration model profiles");
    assert.ok(source.includes("compareModelProfile"), "UI should compare model-profile residuals before/after");
    assert.ok(source.includes("saveModelProfileComparison"), "UI should export model-profile residual comparisons");
    assert.ok(source.includes("checkModelProfileHoldout"), "UI should run holdout model-profile validation");
    assert.ok(source.includes("saveModelProfileHoldout"), "UI should export holdout model-profile validation");
    assert.ok(source.includes("saveModelProfileHoldoutCsv"), "UI should export holdout model-profile validation CSV");
    assert.ok(source.includes("applyModelProfile"), "UI should apply safe calibration model-profile updates");
  }
}

const threeContext = { console: { ...console, warn: () => {} }, window: {} };
threeContext.self = threeContext.window;
threeContext.globalThis = threeContext;
threeContext.window.window = threeContext.window;
vm.createContext(threeContext);
vm.runInContext(fs.readFileSync(path.join(web, "vendor", "three.min.js"), "utf8"), threeContext, { filename: "three.min.js" });
const THREE = threeContext.THREE;
assert.ok(THREE, "vendored Three.js should expose THREE");
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(45, 1, 0.1, 200);
const orthoCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.1, 200);
const group = new THREE.Group();
const mesh = new THREE.Mesh(new THREE.BoxGeometry(1, 1, 0.1), new THREE.MeshStandardMaterial({ color: 0x14799d }));
group.add(mesh);
scene.add(group);
assert.strictEqual(scene.children.length, 1);
assert.strictEqual(group.children.length, 1);
assert.strictEqual(camera.isPerspectiveCamera, true);
assert.strictEqual(orthoCamera.isOrthographicCamera, true);
assert.strictEqual(typeof THREE.Raycaster, "function");
assert.strictEqual(typeof THREE.Vector2, "function");

assert.strictEqual(RAD.backlashActivation(-0.05, 0.1), 0);
assert.strictEqual(RAD.backlashActivation(0.05, 0.1), 0);
assert.strictEqual(Number(RAD.backlashActivation(0.2, 0.1).toFixed(6)), 0.1);
assert.strictEqual(Number(RAD.backlashActivation(-0.2, 0.1).toFixed(6)), -0.1);
assert.strictEqual(Number(RAD.normalizedBacklashToThetaDeadZone(0.1).toFixed(6)), Number(((Math.asin(0.1) * 180) / Math.PI).toFixed(6)));
assert.strictEqual(Number(RAD.normalizedBacklashToAlphaDeadZone(0.1).toFixed(6)), Number((((Math.asin(0.1) * 180) / Math.PI) / 70).toFixed(6)));

const state = RAD.createState(5, 6);
assert.strictEqual(state.cells.alpha.length, 5);
assert.strictEqual(state.cells.theta[0].length, 6);
assert.strictEqual(state.cells.z[2][3], 0);
assert.strictEqual(state.grid.empiricalLockModel, true);
assert.strictEqual(state.grid.lockStateIndex, 1);
assert.strictEqual(state.view.externalCaseVisible, true);
state.grid.cellSize = 1.25;
state.grid.backlash = 0.18;
state.grid.couplingGain = 0.42;
state.grid.zCouplingGain = 0.37;
state.grid.pinRadius = 0.16;
state.grid.holeRadius = 0.22;
state.grid.paperSideLengthMm = 42;
state.grid.paperHoleToleranceMm = 0.12;
state.grid.hardwareProfile = {
  name: "bench-v1",
  source: "caliper",
  sideLengthMm: 42,
  fabricationHoleToleranceMm: 0.12,
  backlashMm: 4.2,
  pinRadiusMm: 7.56,
  holeRadiusMm: 9.45,
  plateThicknessMm: 2.4,
  jointStackHeightMm: null,
  bossRadiusMm: null,
  notes: "validator profile",
};
state.view.camera = {
  mode: "custom",
  projection: "orthographic",
  radius: 18.25,
  theta: -0.7,
  phi: 0.9,
  target: { x: 0.2, y: -0.1, z: 0.05 },
};
state.view.overlayMode = "height";
state.view.simulationMode = "springPreview";
state.view.membraneVisible = false;
state.view.externalCaseVisible = false;
state.view.paintRadius = 2;
state.cells.commandAlpha[2][3] = -0.24;
state.cells.commandZ[2][3] = 0.31;
state.cells.locked[0][0] = true;
state.cells.lockAlpha[0][0] = 0.82;
state.cells.lockZ[0][0] = 0.14;
state.cells.actuatorAllowed[4][5] = false;
state.cells.removed[1][2] = true;
state.selection = { r: 2, c: 3 };
state.experiment.presetName = "custom-roundtrip";
state.experiment.notes = "roundtrip validation";

const serialized = RAD.serialize(state);
const parsedSerialized = JSON.parse(serialized);
assert.ok(Array.isArray(parsedSerialized.cells.alpha), "serialized cells should include derived alpha");
assert.ok(Array.isArray(parsedSerialized.cells.theta), "serialized cells should include derived theta");
assert.ok(Array.isArray(parsedSerialized.cells.z), "serialized cells should include derived z");
assert.ok(Array.isArray(parsedSerialized.cells.lockAlpha), "serialized cells should include committed lock alpha");
assert.ok(Array.isArray(parsedSerialized.cells.lockZ), "serialized cells should include committed lock height");
assert.ok(Array.isArray(parsedSerialized.cells.removed), "serialized cells should include removed-cell topology");
assert.strictEqual(parsedSerialized.experiment.calibrationHoldoutResults, null);
assert.strictEqual(parsedSerialized.experiment.calibrationModelProfileHoldoutValidation, null);
const restored = RAD.deserialize(serialized);
assert.deepStrictEqual(restored.grid.rows, 5);
assert.deepStrictEqual(restored.grid.cols, 6);
assert.strictEqual(restored.grid.cellSize, 1.25);
assert.strictEqual(restored.grid.backlash, 0.18);
assert.strictEqual(restored.grid.couplingGain, 0.42);
assert.strictEqual(restored.grid.zCouplingGain, 0.37);
assert.strictEqual(restored.grid.pinRadius, 0.16);
assert.strictEqual(restored.grid.holeRadius, 0.22);
assert.strictEqual(restored.grid.paperSideLengthMm, 42);
assert.strictEqual(restored.grid.paperHoleToleranceMm, 0.12);
assert.strictEqual(Number(RAD.pinHoleClearance(restored).toFixed(6)), 0.06);
const calibration = RAD.paperRadCalibration(restored);
assert.strictEqual(Number(calibration.mmPerModelUnit.toFixed(6)), Number((42 / 1.25).toFixed(6)));
assert.strictEqual(Number(calibration.configuredBacklashMm.toFixed(6)), Number((0.18 * 42).toFixed(6)));
assert.strictEqual(Number(calibration.referenceBacklashMm.toFixed(6)), 4.2);
assert.strictEqual(Number(calibration.pinHoleClearanceMm.toFixed(6)), Number(((0.22 - 0.16) * 42 / 1.25).toFixed(6)));
assert.strictEqual(Number(calibration.fabricationHoleToleranceModel.toFixed(6)), Number((0.12 * 1.25 / 42).toFixed(6)));
assert.strictEqual(Number(RAD.modelLengthToMm(restored, 1.25).toFixed(6)), 42);
assert.strictEqual(Number(RAD.mmToModelLength(restored, 42).toFixed(6)), 1.25);
assert.strictEqual(restored.grid.hardwareProfile.name, "bench-v1");
assert.strictEqual(restored.grid.hardwareProfile.source, "caliper");
const profileSummary = RAD.calibrationProfileSummary(restored);
assert.strictEqual(profileSummary.measuredCount, 3);
assert.strictEqual(profileSummary.totalCount, 5);
assert.strictEqual(Number(profileSummary.pinHoleClearanceMm.toFixed(6)), 1.89);
assert.strictEqual(Number(profileSummary.pinRadiusModel.toFixed(6)), Number((7.56 * 1.25 / 42).toFixed(6)));
assert.ok(profileSummary.missingFields.includes("jointStackHeightMm"), "profile summary should list missing stack measurement");
const hardwareProfileExport = JSON.parse(RAD.exportHardwareProfileJson(restored));
assert.strictEqual(hardwareProfileExport.schema, "rad-sim.hardware-profile.v1");
assert.strictEqual(hardwareProfileExport.name, "bench-v1");
assert.strictEqual(hardwareProfileExport.dimensionsMm.pinRadiusMm, 7.56);
assert.strictEqual(Number(hardwareProfileExport.derived.pinHoleClearanceMm.toFixed(6)), 1.89);
assert.deepStrictEqual(Array.from(hardwareProfileExport.measuredFields), [
  "pinRadiusMm",
  "holeRadiusMm",
  "plateThicknessMm",
]);
const importedHardwareProfile = RAD.importHardwareProfileJson("\ufeff" + RAD.exportHardwareProfileJson(restored), restored);
assert.strictEqual(importedHardwareProfile.name, "bench-v1");
assert.strictEqual(importedHardwareProfile.source, "caliper");
assert.strictEqual(importedHardwareProfile.holeRadiusMm, 9.45);
const importedCadReferenceHardwareProfile = RAD.importHardwareProfileJson(
  RAD.exportCadRadCellReferenceProfileJson(twoCellBenchState),
  restored
);
assert.strictEqual(importedCadReferenceHardwareProfile.name, "a360-rads-unit-cell-visual-profile");
assert.strictEqual(Number(importedCadReferenceHardwareProfile.holeRadiusMm.toFixed(6)), 1.7);
const importedCadMeshAuditProfile = RAD.importHardwareProfileJson(
  JSON.stringify({
    schema: "rad-sim.cad-mesh-audit.v1",
    hardwareProfile: hardwareProfileExport,
  }),
  restored
);
assert.strictEqual(importedCadMeshAuditProfile.name, "bench-v1");
assert.strictEqual(importedCadMeshAuditProfile.plateThicknessMm, 2.4);
const importedFlatHardwareProfile = RAD.importHardwareProfileJson(
  JSON.stringify({
    name: "flat-profile",
    source: "browser-flat",
    sideLengthMm: 40,
    fabricationHoleToleranceMm: 0.08,
    backlashMm: 4,
    pinRadiusMm: 5,
    holeRadiusMm: 5.8,
  }),
  restored
);
assert.strictEqual(importedFlatHardwareProfile.name, "flat-profile");
assert.strictEqual(importedFlatHardwareProfile.sideLengthMm, 40);
assert.strictEqual(importedFlatHardwareProfile.pinRadiusMm, 5);
const partialReadiness = RAD.calibrationReadiness(restored);
assert.strictEqual(partialReadiness.level, "partial-measured");
assert.strictEqual(partialReadiness.visualReady, false);
assert.strictEqual(partialReadiness.meshReady, false);
assert.strictEqual(partialReadiness.solverReady, false);
assert.ok(partialReadiness.solverGaps.includes("friction and contact characterization"), "readiness should preserve solver calibration gaps");
const partialMeasurementPlan = RAD.calibrationMeasurementPlan(restored);
const missingPlanTasks = partialMeasurementPlan.filter((task) => task.status === "missing");
const donePlanTasks = partialMeasurementPlan.filter((task) => task.status === "done");
assert.strictEqual(missingPlanTasks[0].id, "geometry_jointStackHeightMm");
assert.strictEqual(missingPlanTasks[1].id, "geometry_bossRadiusMm");
assert.strictEqual(donePlanTasks[0].id, "geometry_pinRadiusMm");
assert.ok(
  missingPlanTasks.some((task) => task.id === "solver_response_data" && task.category === "solver"),
  "measurement plan should include response-data solver task"
);
const exportedPlan = JSON.parse(RAD.exportCalibrationMeasurementPlan(restored));
assert.strictEqual(exportedPlan.schema, "rad-sim.calibration-plan.v1");
assert.strictEqual(exportedPlan.profile.name, "bench-v1");
assert.strictEqual(exportedPlan.readiness.level, "partial-measured");
assert.strictEqual(exportedPlan.tasks.length, partialMeasurementPlan.length);
const readyState = RAD.createState(2, 2);
readyState.grid.hardwareProfile = {
  ...state.grid.hardwareProfile,
  jointStackHeightMm: 3.8,
  bossRadiusMm: 1.6,
};
const readyReadiness = RAD.calibrationReadiness(readyState);
assert.strictEqual(readyReadiness.level, "mesh-calibrated");
assert.strictEqual(readyReadiness.visualReady, true);
assert.strictEqual(readyReadiness.meshReady, true);
assert.strictEqual(readyReadiness.solverReady, false);
const readyMeasurementPlan = RAD.calibrationMeasurementPlan(readyState);
assert.ok(readyMeasurementPlan.filter((task) => task.category === "geometry").every((task) => task.status === "done"));
assert.ok(readyMeasurementPlan.filter((task) => task.category === "solver").every((task) => task.status === "missing"));
readyState.selection = { r: 1, c: 1 };
const experimentProtocol = RAD.calibrationExperimentProtocol(readyState);
assert.strictEqual(experimentProtocol.schema, "rad-sim.calibration-experiment-protocol.v1");
assert.strictEqual(experimentProtocol.hardwareProfile, "bench-v1");
assert.strictEqual(experimentProtocol.steps.length, 7);
assert.ok(experimentProtocol.steps.some((step) => step.id === "pair_z_residual" && step.scope === "pair"));
assert.ok(experimentProtocol.steps.some((step) => step.id === "locked_cell_control" && step.lockedCells.length === 1));
assert.ok(experimentProtocol.measurementFields.includes("pin_hole_slip_mm"));
const exportedExperimentProtocol = JSON.parse(RAD.exportCalibrationExperimentProtocol(readyState));
assert.strictEqual(exportedExperimentProtocol.schema, experimentProtocol.schema);
assert.strictEqual(exportedExperimentProtocol.steps[0].id, "single_alpha_contract");
const browserAtlas = RAD.responseAtlas(readyState);
assert.strictEqual(browserAtlas.schema, "rad-sim.response-atlas.v1");
assert.strictEqual(browserAtlas.summary.entryCount, experimentProtocol.steps.length);
assert.strictEqual(browserAtlas.summary.scopeCounts.single, 3);
assert.strictEqual(browserAtlas.summary.scopeCounts.pair, 2);
assert.strictEqual(browserAtlas.summary.scopeCounts.cluster, 1);
assert.strictEqual(browserAtlas.summary.scopeCounts.lock, 1);
assert.ok(browserAtlas.entries.some((entry) => entry.stepId === "single_z_lift" && entry.observationCells.length >= 2), "response atlas should include residual z observation cells");
assert.strictEqual(browserAtlas.entries[0].physicalValidation.schema, "rad-sim.browser-physical-preview.v1");
assert.strictEqual(JSON.parse(RAD.exportResponseAtlas(readyState)).schema, browserAtlas.schema);
const browserSweep = RAD.responseAtlasSweep(readyState, {
  backlashValues: [0.02, 0.14],
  clearanceValues: [0.02, 0.12],
});
assert.strictEqual(browserSweep.schema, "rad-sim.response-atlas-sweep.v1");
assert.strictEqual(browserSweep.summary.sampleCount, 4);
assert.deepStrictEqual(Array.from(browserSweep.parameters.backlashValues), [0.02, 0.14]);
assert.deepStrictEqual(Array.from(browserSweep.parameters.pinHoleClearanceValues), [0.02, 0.12]);
assert.strictEqual(browserSweep.samples[0].atlas.schema, "rad-sim.response-atlas.v1");
assert.ok(browserSweep.summary.maxObservedNeighborZResidual > 0, "response atlas sweep should track residual neighbor z motion");
assert.ok(
  browserSweep.trends.byPinHoleClearance[0].maxObservedNeighborZResidual >
    browserSweep.trends.byPinHoleClearance[1].maxObservedNeighborZResidual,
  "response atlas sweep should show clearance-gated neighbor residual trend"
);
assert.strictEqual(browserSweep.sensitivity.method, "endpoint finite difference over each parameter trend");
assert.ok(
  browserSweep.sensitivity.metrics.pinHoleClearance.maxObservedNeighborZResidual < 0,
  "response atlas sweep should report negative clearance sensitivity for neighbor residual"
);
assert.ok(browserSweep.sensitivity.dominant, "response atlas sweep should identify a dominant sensitivity");
assert.strictEqual(
  browserSweep.operatorLawCandidates.method,
  "adjacent monotonicity over sampled parameter trend plus endpoint sensitivity"
);
const clearanceResidualLaw = browserSweep.operatorLawCandidates.laws.find(
  (law) => law.parameter === "pinHoleClearance" && law.metric === "maxObservedNeighborZResidual"
);
assert.ok(clearanceResidualLaw, "response atlas sweep should propose a clearance residual law");
assert.strictEqual(clearanceResidualLaw.monotonicity, "decreasing");
assert.strictEqual(clearanceResidualLaw.supportedBySweep, true);
assert.strictEqual(clearanceResidualLaw.status, "simulator-diagnostic");
assert.strictEqual(JSON.parse(RAD.exportResponseAtlasSweep(readyState, { backlashValues: [0.02], clearanceValues: [0.02] })).schema, browserSweep.schema);
const resultsTemplate = RAD.calibrationExperimentResultsTemplate(readyState, { protocol: experimentProtocol });
assert.strictEqual(resultsTemplate.schema, "rad-sim.calibration-experiment-results.v1");
assert.strictEqual(resultsTemplate.provenance.schema, "rad-sim.calibration-dataset-provenance.v1");
assert.strictEqual(resultsTemplate.provenance.datasetRole, "unassigned");
assert.ok(resultsTemplate.provenance.notes.includes("datasetRole"), "results template should prompt fit/holdout provenance");
assert.strictEqual(resultsTemplate.steps.length, experimentProtocol.steps.reduce((sum, step) => sum + step.repeatCount, 0));
const exportedResultsTemplate = JSON.parse(RAD.exportCalibrationExperimentResultsTemplate(readyState, { protocol: experimentProtocol }));
assert.strictEqual(exportedResultsTemplate.schema, resultsTemplate.schema);
const benchNotebook = RAD.calibrationBenchNotebook(readyState, { protocol: experimentProtocol });
assert.strictEqual(benchNotebook.schema, "rad-sim.calibration-bench-notebook.v1");
assert.strictEqual(benchNotebook.datasetPlan.fit.role, "fit");
assert.strictEqual(benchNotebook.datasetPlan.holdout.role, "holdout");
assert.strictEqual(benchNotebook.protocol.stepCount, experimentProtocol.steps.length);
assert.strictEqual(benchNotebook.scenarios.length, experimentProtocol.steps.length);
assert.ok(benchNotebook.instruments.some((instrument) => instrument.id === "motion-tracking-camera"), "bench notebook should list response measurement instruments");
assert.ok(benchNotebook.outputs.some((output) => output.id === "holdout-validation-csv"), "bench notebook should point to holdout CSV output");
assert.ok(benchNotebook.scenarios.some((scenario) => scenario.stepId === "pair_z_residual"), "bench notebook should include residual z scenario");
assert.strictEqual(JSON.parse(RAD.exportCalibrationBenchNotebook(readyState, { protocol: experimentProtocol })).schema, benchNotebook.schema);
const benchCsv = RAD.exportCalibrationBenchNotebookCsv(benchNotebook);
assert.ok(benchCsv.split("\n")[0].includes("step_id"), "bench notebook CSV should expose scenario IDs");
assert.ok(benchCsv.includes("pair_z_residual"), "bench notebook CSV should include the residual z scenario");
assert.ok(benchCsv.includes("pin_hole_slip_mm"), "bench notebook CSV should include measurement fields");
const benchPacket = RAD.calibrationBenchPacket(readyState, {
  protocol: experimentProtocol,
  fitDatasetId: "browser-fit",
  holdoutDatasetId: "browser-holdout",
  profileId: "browser-profile-v1",
  profileFrozenAt: "2026-08-09T12:30:00Z",
});
assert.strictEqual(benchPacket.schema, "rad-sim.calibration-bench-packet.v1");
assert.strictEqual(benchPacket.benchNotebook.schema, benchNotebook.schema);
assert.strictEqual(benchPacket.fitResultsTemplate.provenance.datasetRole, "fit");
assert.strictEqual(benchPacket.holdoutResultsTemplate.provenance.datasetRole, "holdout");
assert.strictEqual(benchPacket.holdoutResultsTemplate.provenance.profileId, "browser-profile-v1");
assert.ok(benchPacket.benchNotebookCsv.includes("pair_z_residual"), "bench packet should embed the notebook CSV table");
assert.ok(benchPacket.artifactManifest.some((item) => item.id === "holdout-template"), "bench packet should include a holdout template manifest entry");
assert.strictEqual(
  JSON.parse(RAD.exportCalibrationBenchPacket(readyState, { protocol: experimentProtocol })).schema,
  benchPacket.schema
);
const liftStep = experimentProtocol.steps.find((step) => step.id === "single_z_lift");
const liftResult = resultsTemplate.steps.find((step) => step.stepId === "single_z_lift" && step.repeatIndex === 1);
function stateForProtocolStep(baseState, step, includeCommands) {
  const temp = RAD.deserialize(RAD.serialize(baseState));
  for (let r = 0; r < temp.grid.rows; r += 1) {
    for (let c = 0; c < temp.grid.cols; c += 1) {
      temp.cells.commandAlpha[r][c] = 0;
      temp.cells.commandZ[r][c] = 0;
      temp.cells.locked[r][c] = false;
    }
  }
  for (const locked of step.lockedCells || []) temp.cells.locked[locked.row][locked.col] = true;
  if (includeCommands) {
    for (const command of step.commands || []) {
      temp.cells.commandAlpha[command.row][command.col] += Number(command.alpha) || 0;
      temp.cells.commandZ[command.row][command.col] += Number(command.z) || 0;
    }
  }
  return temp;
}
const liftCommandState = stateForProtocolStep(readyState, liftStep, true);
const liftBaselineState = stateForProtocolStep(readyState, liftStep, false);
const liftSim = RAD.simulate(liftCommandState);
const liftBaseline = RAD.simulate(liftBaselineState);
for (const cell of liftResult.cells) {
  cell.alphaDelta = (liftSim.alpha[cell.row][cell.col] || 0) - (liftBaseline.alpha[cell.row][cell.col] || 0);
  cell.heightDelta = (liftSim.height[cell.row][cell.col] || 0) - (liftBaseline.height[cell.row][cell.col] || 0);
  cell.pinHoleSlipMm = 0.11;
  cell.actuatorForceN = 2.5;
}
const comparison = RAD.compareCalibrationExperimentResults(readyState, resultsTemplate, { protocol: experimentProtocol });
const liftComparison = comparison.comparisons.find((item) => item.stepId === "single_z_lift" && item.repeatIndex === 1);
assert.strictEqual(comparison.schema, "rad-sim.calibration-experiment-comparison.v1");
assert.strictEqual(liftComparison.missingObservationCount, 0);
assert.strictEqual(Number(liftComparison.heightRmse.toFixed(8)), 0);
assert.strictEqual(Number(liftComparison.meanPinHoleSlipMm.toFixed(8)), 0.11);
assert.strictEqual(Number(liftComparison.meanSignedHeightError.toFixed(8)), 0);
assert.strictEqual(Number(liftComparison.meanAbsHeightError.toFixed(8)), 0);
assert.ok(Array.isArray(comparison.field.combinedError), "calibration comparison should include per-cell error field");
assert.ok(comparison.field.sampleCount.some((row) => row.some((value) => value > 0)), "calibration field should count measured cells");
assert.strictEqual(Number(comparison.field.maxCombinedError.toFixed(8)), 0);
assert.ok(comparison.field.worstCell, "calibration comparison should track the strongest measured cell");
assert.ok(comparison.field.topCells.length > 0, "calibration comparison should rank measured cells");
assert.deepStrictEqual(comparison.field.topCells[0], comparison.field.worstCell);
assert.ok(Array.isArray(comparison.fitResidualField.combinedError), "calibration comparison should include fitted residual field");
assert.strictEqual(Number(comparison.fitResidualField.maxCombinedError.toFixed(8)), 0);
assert.ok(comparison.fitResidualField.topCells.length > 0, "calibration residual field should rank measured cells");
assert.strictEqual(comparison.fit.height.sampleCount, liftResult.cells.length);
assert.ok(comparison.fit.height.gainIdentifiable, "calibration fit should identify height gain when predicted samples vary");
assert.strictEqual(Number(comparison.fit.height.suggestedGain.toFixed(8)), 1);
assert.strictEqual(Number(comparison.fit.height.suggestedBias.toFixed(8)), 0);
const comparisonSummary = RAD.summarizeCalibrationComparison(comparison);
assert.strictEqual(comparisonSummary.schema, "rad-sim.calibration-experiment-comparison-summary.v1");
assert.strictEqual(comparisonSummary.stepCount, comparison.comparisons.length);
assert.ok(comparisonSummary.measuredCellCount >= liftResult.cells.length);
assert.strictEqual(Number(comparisonSummary.heightRmseMean.toFixed(8)), 0);
assert.strictEqual(Number(comparisonSummary.maxCombinedError.toFixed(8)), 0);
assert.strictEqual(Number(comparisonSummary.meanSignedHeightError.toFixed(8)), 0);
assert.strictEqual(Number(comparisonSummary.fit.height.rmsResidual.toFixed(8)), 0);
assert.strictEqual(Number(comparisonSummary.fitResidualMaxCombinedError.toFixed(8)), 0);
const perturbedResults = JSON.parse(JSON.stringify(resultsTemplate));
const perturbedLift = perturbedResults.steps.find((step) => step.stepId === "single_z_lift" && step.repeatIndex === 1);
perturbedLift.cells[0].heightDelta += 0.05;
const perturbedComparison = RAD.compareCalibrationExperimentResults(readyState, perturbedResults, { protocol: experimentProtocol });
assert.ok(perturbedComparison.field.maxAbsHeightError > 0, "perturbed height measurement should produce a height-error field");
assert.ok(perturbedComparison.field.maxCombinedError > 0, "perturbed height measurement should produce a visible calibration-error overlay scale");
assert.ok(perturbedComparison.fit.height.rmsRawError > 0, "perturbed height measurement should produce raw fit error");
assert.ok(perturbedComparison.fit.height.sampleCount > 0, "perturbed height fit should keep sample count");
assert.ok(
  perturbedComparison.fitResidualField.sampleCount.some((row) => row.some((value) => value > 0)),
  "perturbed height measurement should populate a fitted residual field"
);
assert.strictEqual(perturbedComparison.field.worstCell.row, perturbedLift.cells[0].row);
assert.strictEqual(perturbedComparison.field.worstCell.col, perturbedLift.cells[0].col);
const perturbedSummary = RAD.summarizeCalibrationComparison(perturbedComparison);
assert.deepStrictEqual(perturbedSummary.worstCell, perturbedComparison.field.worstCell);
assert.deepStrictEqual(perturbedSummary.fitResidualWorstCell, perturbedComparison.fitResidualField.worstCell);
assert.deepStrictEqual(perturbedSummary.topCells[0], perturbedComparison.field.topCells[0]);
assert.deepStrictEqual(perturbedSummary.fitResidualTopCells[0], perturbedComparison.fitResidualField.topCells[0]);
readyState.experiment.calibrationResults = perturbedResults;
readyState.experiment.calibrationComparison = perturbedComparison;
readyState.experiment.calibrationComparisonSummary = perturbedSummary;
const report = RAD.calibrationComparisonReport(readyState);
assert.strictEqual(report.schema, "rad-sim.calibration-comparison-report.v1");
assert.strictEqual(report.parameterEstimates.schema, "rad-sim.calibration-parameter-estimates.v1");
assert.strictEqual(report.parameterEstimates.estimates.alphaResponse.status, "estimated-from-results");
assert.ok(
  report.parameterEstimates.estimates.zCouplingGain.sampleCount > 0,
  "parameter report should estimate z-coupling from neighbor/direct height ratios"
);
assert.strictEqual(report.parameterEstimates.estimates.backlash.status, "not-identifiable-from-current-results");
assert.ok(
  report.parameterEstimates.estimates.contactProxy.forceHeightSlopeNPerModelUnit > 0,
  "parameter report should expose a force-height proxy when force data exists"
);
assert.strictEqual(report.modelProfile.schema, "rad-sim.calibration-model-profile.v1");
const zProfileUpdate = report.modelProfile.recommendedUpdates.find((update) => update.name === "zCouplingGain");
assert.ok(zProfileUpdate.safeToApply, "bounded measured z coupling should be marked safe to apply");
assert.strictEqual(zProfileUpdate.configField, "zCouplingGain");
const alphaProfileUpdate = report.modelProfile.recommendedUpdates.find((update) => update.name === "alphaResponseGain");
assert.strictEqual(alphaProfileUpdate.safeToApply, false, "alpha response gain should remain diagnostic in v1");
const exportedModelProfile = JSON.parse(RAD.exportCalibrationModelProfile(readyState));
assert.strictEqual(exportedModelProfile.schema, "rad-sim.calibration-model-profile.v1");
const lowZProfileState = RAD.createState(3, 3);
lowZProfileState.grid.zCouplingGain = 0.05;
lowZProfileState.selection = { r: 1, c: 1 };
const highZProfileState = RAD.createState(3, 3);
highZProfileState.grid.zCouplingGain = 0.5;
highZProfileState.selection = { r: 1, c: 1 };
const profileProtocol = RAD.calibrationExperimentProtocol(lowZProfileState);
const profileResults = RAD.calibrationExperimentResultsTemplate(lowZProfileState, { protocol: profileProtocol });
const profileLiftStep = profileProtocol.steps.find((step) => step.id === "single_z_lift");
const profileLiftResult = profileResults.steps.find((step) => step.stepId === "single_z_lift" && step.repeatIndex === 1);
const profileHighCommand = stateForProtocolStep(highZProfileState, profileLiftStep, true);
const profileHighBaseline = stateForProtocolStep(highZProfileState, profileLiftStep, false);
const profileHighSim = RAD.simulate(profileHighCommand);
const profileHighBaseSim = RAD.simulate(profileHighBaseline);
for (const cell of profileLiftResult.cells) {
  cell.alphaDelta = (profileHighSim.alpha[cell.row][cell.col] || 0) - (profileHighBaseSim.alpha[cell.row][cell.col] || 0);
  cell.heightDelta = (profileHighSim.height[cell.row][cell.col] || 0) - (profileHighBaseSim.height[cell.row][cell.col] || 0);
}
lowZProfileState.experiment.calibrationResults = profileResults;
lowZProfileState.experiment.calibrationComparison = RAD.compareCalibrationExperimentResults(lowZProfileState, profileResults, { protocol: profileProtocol });
lowZProfileState.experiment.calibrationComparisonSummary = RAD.summarizeCalibrationComparison(lowZProfileState.experiment.calibrationComparison);
const profileResidualComparison = RAD.calibrationModelProfileResidualComparison(lowZProfileState, null, { protocol: profileProtocol });
assert.strictEqual(profileResidualComparison.schema, "rad-sim.calibration-model-profile-residual-comparison.v1");
assert.ok(profileResidualComparison.application.appliedUpdates.length > 0, "profile residual comparison should apply bounded safe updates to the trial state");
assert.ok(
  profileResidualComparison.after.metrics.residualScore < profileResidualComparison.before.metrics.residualScore,
  "profile residual comparison should detect lower residual after applying the fitted z coupling"
);
assert.strictEqual(
  JSON.parse(RAD.exportCalibrationModelProfileResidualComparison(lowZProfileState, null, { protocol: profileProtocol })).delta.residualScore,
  profileResidualComparison.delta.residualScore
);
const rejectedProfileComparison = JSON.parse(JSON.stringify(profileResidualComparison));
rejectedProfileComparison.application.appliedUpdates = [];
rejectedProfileComparison.delta.residualScore = 0;
const profileSelection = RAD.selectCalibrationModelProfile([rejectedProfileComparison, profileResidualComparison]);
assert.strictEqual(profileSelection.schema, "rad-sim.calibration-model-profile-selection.v1");
assert.strictEqual(profileSelection.candidateCount, 2);
assert.strictEqual(profileSelection.eligibleCount, 1);
assert.strictEqual(profileSelection.selectedIndex, 1);
assert.strictEqual(JSON.parse(RAD.exportCalibrationModelProfileSelection([rejectedProfileComparison, profileResidualComparison])).selectedIndex, 1);
const holdoutProfileState = RAD.createState(3, 3);
holdoutProfileState.grid.zCouplingGain = 0.45;
holdoutProfileState.selection = { r: 1, c: 1 };
const holdoutResults = RAD.calibrationExperimentResultsTemplate(lowZProfileState, { protocol: profileProtocol });
for (const measuredStep of holdoutResults.steps) {
  const protocolStep = profileProtocol.steps.find((step) => step.id === measuredStep.stepId);
  const holdoutCommand = stateForProtocolStep(holdoutProfileState, protocolStep, true);
  const holdoutBaseline = stateForProtocolStep(holdoutProfileState, protocolStep, false);
  const holdoutSim = RAD.simulate(holdoutCommand);
  const holdoutBaseSim = RAD.simulate(holdoutBaseline);
  for (const cell of measuredStep.cells) {
    cell.alphaDelta = (holdoutSim.alpha[cell.row][cell.col] || 0) - (holdoutBaseSim.alpha[cell.row][cell.col] || 0);
    cell.heightDelta = (holdoutSim.height[cell.row][cell.col] || 0) - (holdoutBaseSim.height[cell.row][cell.col] || 0);
  }
}
lowZProfileState.experiment.calibrationHoldoutResults = holdoutResults;
lowZProfileState.experiment.calibrationHoldoutComparison = RAD.compareCalibrationExperimentResults(lowZProfileState, holdoutResults, { protocol: profileProtocol });
lowZProfileState.experiment.calibrationHoldoutComparisonSummary = RAD.summarizeCalibrationComparison(lowZProfileState.experiment.calibrationHoldoutComparison);
const profileHoldoutValidation = RAD.calibrationModelProfileHoldoutValidation(lowZProfileState, holdoutResults, null, { protocol: profileProtocol });
assert.strictEqual(profileHoldoutValidation.schema, "rad-sim.calibration-model-profile-holdout-validation.v1");
assert.strictEqual(profileHoldoutValidation.fitPass, true);
assert.strictEqual(profileHoldoutValidation.holdoutPass, true);
assert.strictEqual(profileHoldoutValidation.residualValidationPass, true);
assert.strictEqual(profileHoldoutValidation.independentValidationPass, false);
assert.strictEqual(profileHoldoutValidation.splitMetadata.schema, "rad-sim.calibration-train-holdout-split.v1");
assert.strictEqual(profileHoldoutValidation.splitMetadata.independenceStatus, "requires-external-bench-protocol");
assert.ok(profileHoldoutValidation.splitMetadata.missingEvidence.includes("fitDatasetId"), "holdout validation should warn when fit dataset ID is missing");
assert.ok(
  profileHoldoutValidation.holdout.after.metrics.residualScore < profileHoldoutValidation.holdout.before.metrics.residualScore,
  "holdout validation should detect residual reduction on separately supplied holdout data"
);
const holdoutCsv = RAD.exportCalibrationModelProfileHoldoutValidationCsv(profileHoldoutValidation);
assert.ok(holdoutCsv.split("\n")[0].includes("independence_status"), "holdout CSV should include split independence status");
assert.ok(holdoutCsv.includes("holdout"), "holdout CSV should include a holdout row");
assert.ok(holdoutCsv.includes("requires-external-bench-protocol"), "holdout CSV should preserve split caveat");
assert.ok(holdoutCsv.includes("independent_validation_pass"), "holdout CSV should report independent validation separately from residual validation");
const documentedProfile = JSON.parse(JSON.stringify(profileHoldoutValidation.modelProfile));
documentedProfile.profileId = "browser-synthetic-profile-v1";
const documentedFitResults = JSON.parse(JSON.stringify(profileResults));
documentedFitResults.provenance = {
  schema: "rad-sim.calibration-dataset-provenance.v1",
  datasetId: "browser-fit-synthetic-z-050",
  datasetRole: "fit",
  sourceFileId: "browser-fit-synthetic-z-050.json",
  collectedAt: "2026-08-09T10:00:00Z",
  operator: "validator",
};
const documentedHoldoutResults = JSON.parse(JSON.stringify(holdoutResults));
documentedHoldoutResults.provenance = {
  schema: "rad-sim.calibration-dataset-provenance.v1",
  datasetId: "browser-holdout-synthetic-z-045",
  datasetRole: "holdout",
  sourceFileId: "browser-holdout-synthetic-z-045.json",
  collectedAt: "2026-08-09T11:00:00Z",
  operator: "validator",
  profileId: "browser-synthetic-profile-v1",
  profileFrozenAt: "2026-08-09T10:30:00Z",
};
lowZProfileState.experiment.calibrationResults = documentedFitResults;
lowZProfileState.experiment.calibrationHoldoutResults = documentedHoldoutResults;
const documentedHoldoutValidation = RAD.calibrationModelProfileHoldoutValidation(
  lowZProfileState,
  documentedHoldoutResults,
  documentedProfile,
  { protocol: profileProtocol }
);
assert.strictEqual(documentedHoldoutValidation.residualValidationPass, true);
assert.strictEqual(documentedHoldoutValidation.independentValidationPass, true);
assert.strictEqual(documentedHoldoutValidation.provenanceWarnings.length, 0);
assert.strictEqual(documentedHoldoutValidation.splitMetadata.independenceStatus, "documented-independent-holdout");
assert.strictEqual(documentedHoldoutValidation.splitMetadata.fitDataset.datasetId, "browser-fit-synthetic-z-050");
const documentedHoldoutCsv = RAD.exportCalibrationModelProfileHoldoutValidationCsv(documentedHoldoutValidation);
assert.ok(documentedHoldoutCsv.includes("documented-independent-holdout"), "documented holdout CSV should preserve independence status");
const benchExecutionValidation = RAD.calibrationBenchExecutionValidation(
  lowZProfileState,
  documentedHoldoutResults,
  documentedProfile,
  { protocol: profileProtocol }
);
assert.strictEqual(benchExecutionValidation.schema, "rad-sim.calibration-bench-execution-validation.v1");
assert.strictEqual(benchExecutionValidation.summary.status, "documented-independent-validation");
assert.strictEqual(benchExecutionValidation.summary.residualValidationPass, true);
assert.strictEqual(benchExecutionValidation.summary.independentValidationPass, true);
assert.strictEqual(benchExecutionValidation.summary.executionValidationPass, true);
assert.strictEqual(benchExecutionValidation.summary.missingEvidenceCount, 0);
assert.strictEqual(benchExecutionValidation.fitDataset.datasetId, "browser-fit-synthetic-z-050");
assert.strictEqual(benchExecutionValidation.holdoutDataset.datasetId, "browser-holdout-synthetic-z-045");
assert.strictEqual(benchExecutionValidation.formalization.targetId, "calibration_bench_executed_validation_gate");
const benchExecutionCsv = RAD.exportCalibrationBenchExecutionValidationCsv(benchExecutionValidation);
assert.ok(benchExecutionCsv.includes("execution_validation_pass"), "execution CSV should expose the gate status");
assert.ok(benchExecutionCsv.includes("browser-fit-synthetic-z-050"), "execution CSV should include fit dataset provenance");
assert.ok(benchExecutionCsv.includes("browser-holdout-synthetic-z-045"), "execution CSV should include holdout dataset provenance");
assert.strictEqual(
  JSON.parse(RAD.exportCalibrationBenchExecutionValidation(lowZProfileState, documentedHoldoutResults, documentedProfile, { protocol: profileProtocol })).summary.executionValidationPass,
  true
);
const syntheticVerticalLoadComparison = {
  schema: "rad-sim.vertical-load-energy-comparison-report.v1",
  summary: {
    scenarioCount: 1,
    missingMeasurementCount: 0,
    allScenariosPassTolerance: true,
    maxAbsWorkOrContactError: 0,
  },
  scenarios: [
    {
      name: "browser-readiness-smoke",
      validation: {
        intact: {
          errors: { signedLoadWork: 0, loadWorkMagnitude: 0, heightContactPenalty: 0 },
          measured: { signedLoadWork: 0, loadWorkMagnitude: 0, heightContactPenalty: 0 },
          simulatedObserved: { signedLoadWork: 0, loadWorkMagnitude: 0, heightContactPenalty: 0 },
        },
        removed: {
          errors: { signedLoadWork: 0, loadWorkMagnitude: 0, heightContactPenalty: 0 },
          measured: { signedLoadWork: 0, loadWorkMagnitude: 0, heightContactPenalty: 0 },
          simulatedObserved: { signedLoadWork: 0, loadWorkMagnitude: 0, heightContactPenalty: 0 },
        },
      },
    },
  ],
};
const physicalReadiness = RAD.physicalValidationReadinessReport(lowZProfileState, {
  calibrationExecutionValidation: benchExecutionValidation,
  verticalLoadComparisonReport: syntheticVerticalLoadComparison,
});
assert.strictEqual(physicalReadiness.schema, "rad-sim.physical-validation-readiness.v1");
assert.strictEqual(physicalReadiness.summary.physicalValidationReady, true);
assert.strictEqual(physicalReadiness.summary.missingEvidenceCount, 0);
assert.ok(physicalReadiness.evidence.loadProxyTermCount > 0, "physical readiness should count load-work proxy terms");
assert.ok(physicalReadiness.evidence.contactProxyTermCount > 0, "physical readiness should count contact proxy terms");
assert.strictEqual(physicalReadiness.formalization.targetId, "physical_validation_readiness_gate");
const physicalReadinessCsv = RAD.exportPhysicalValidationReadinessCsv(physicalReadiness);
assert.ok(physicalReadinessCsv.includes("physical_validation_ready"), "physical readiness CSV should expose the gate status");
assert.ok(physicalReadinessCsv.includes("ready-for-physical-claim-review"), "physical readiness CSV should preserve readiness status");
assert.strictEqual(
  JSON.parse(RAD.exportPhysicalValidationReadiness(lowZProfileState, {
    calibrationExecutionValidation: benchExecutionValidation,
    verticalLoadComparisonReport: syntheticVerticalLoadComparison,
  })).summary.physicalValidationReady,
  true
);
const missingPhysicalReadiness = RAD.physicalValidationReadinessReport(lowZProfileState);
assert.strictEqual(missingPhysicalReadiness.summary.physicalValidationReady, false);
assert.ok(missingPhysicalReadiness.summary.missingEvidence.includes("verticalLoadComparisonReport"), "physical readiness should fail without vertical-load evidence");
const contactState = RAD.createState(1, 3);
contactState.grid.pinRadius = 0.10;
contactState.grid.holeRadius = 0.12;
contactState.cells.commandZ[0][0] = 0.01;
contactState.cells.commandZ[0][1] = 0.08;
contactState.cells.removed[0][2] = true;
const contactReport = RAD.contactStateAbstractionReport(contactState, { contactStiffness: 2 });
assert.strictEqual(contactReport.schema, "rad-sim.contact-state-abstraction.v1");
assert.strictEqual(contactReport.summary.contactStateAbstractionReady, true);
assert.strictEqual(contactReport.topology.activeBodyCount, 2);
assert.strictEqual(contactReport.topology.removedBodyCount, 1);
assert.strictEqual(contactReport.topology.penaltyTermCount, 2);
const contactModes = new Map(contactReport.cells.map((cell) => [`${cell.row},${cell.col}`, cell.mode]));
assert.strictEqual(contactModes.get("0,0"), "free-clearance");
assert.strictEqual(contactModes.get("0,1"), "engaged");
assert.strictEqual(contactModes.get("0,2"), "removed");
assert.strictEqual(contactReport.formalization.targetId, "contact_state_abstraction_gate");
const contactCsv = RAD.exportContactStateAbstractionCsv(contactReport);
assert.ok(contactCsv.includes("contact_penalty"), "contact-state CSV should expose penalty terms");
assert.ok(contactCsv.includes("engaged"), "contact-state CSV should include contact modes");
assert.strictEqual(JSON.parse(RAD.exportContactStateAbstraction(contactState, { contactStiffness: 2 })).summary.engagedContactCount, 1);
const graphContactState = RAD.createState(2, 2);
graphContactState.grid.pinRadius = 0.10;
graphContactState.grid.holeRadius = 0.12;
graphContactState.cells.removed[0][1] = true;
graphContactState.cells.commandZ[1][1] = 0.08;
graphContactState.cells.commandAlpha[0][0] = 0.2;
const graphContactReport = RAD.contactStateAbstractionReport(graphContactState);
const graphConsistency = RAD.contactGraphConsistencyReport(graphContactState, {
  contactReport: graphContactReport,
  groupSupport: [{ row: 0, col: 1 }, { row: 1, col: 1 }],
});
assert.strictEqual(graphConsistency.schema, "rad-sim.contact-graph-consistency.v1");
assert.strictEqual(graphConsistency.summary.contactGraphConsistent, true);
assert.strictEqual(graphConsistency.formalization.targetId, "contact_graph_consistency_gate");
assert.strictEqual(graphConsistency.graph.activeBodyCount, 3);
assert.strictEqual(graphConsistency.graph.removedBodyCount, 1);
assert.strictEqual(graphConsistency.graph.totalEdgeCount, 4);
assert.strictEqual(graphConsistency.graph.activeEdgeCount, 2);
assert.strictEqual(graphConsistency.graph.deletedEdgeCount, 2);
assert.strictEqual(graphConsistency.graph.removedIncidentActiveEdgeCount, 0);
assert.strictEqual(graphConsistency.contact.removedActiveContactCount, 0);
assert.strictEqual(graphConsistency.support.supportCellCount, 2);
assert.strictEqual(graphConsistency.support.removedSupportCellCount, 1);
assert.strictEqual(graphConsistency.support.supportContactRecordCount, 2);
const graphConsistencyCsv = RAD.exportContactGraphConsistencyCsv(graphConsistency);
assert.ok(graphConsistencyCsv.includes("contact_graph_consistent"), "contact graph CSV should expose the gate status");
assert.strictEqual(JSON.parse(RAD.exportContactGraphConsistency(graphContactState, {
  contactReport: graphContactReport,
  groupSupport: [{ row: 0, col: 1 }, { row: 1, col: 1 }],
})).summary.contactGraphConsistent, true);
const inconsistentGraphContactReport = JSON.parse(JSON.stringify(graphContactReport));
const removedGraphContactCell = inconsistentGraphContactReport.cells.find((cell) => cell.row === 0 && cell.col === 1);
removedGraphContactCell.mode = "engaged";
removedGraphContactCell.contactState.unilateralContactActive = true;
const inconsistentGraph = RAD.contactGraphConsistencyReport(graphContactState, {
  contactReport: inconsistentGraphContactReport,
  groupSupport: [{ row: 0, col: 1 }, { row: 1, col: 1 }],
});
assert.strictEqual(inconsistentGraph.summary.contactGraphConsistent, false);
assert.ok(
  inconsistentGraph.summary.missingEvidence.includes("removedActiveContacts"),
  "removed active contacts should fail the contact graph consistency gate"
);
const realizationState = RAD.createState(2, 2);
realizationState.grid.pinRadius = 0.10;
realizationState.grid.holeRadius = 0.12;
const realizationEvents = [
  RAD.groupActuationEvent([{ r: 0, c: 0 }, { r: 0, c: 1 }], 0.2, 0.05),
  RAD.lockEvent({ r: 1, c: 1 }),
  RAD.removeCellEvent({ r: 0, c: 1 }),
];
const realizationReport = RAD.physicalRealizationMapReport(realizationState, { eventSequence: realizationEvents });
assert.strictEqual(realizationReport.schema, "rad-sim.physical-realization-map.v1");
assert.strictEqual(realizationReport.summary.physicalRealizationMapReady, true);
assert.strictEqual(realizationReport.summary.abstractOperatorCount, 3);
assert.strictEqual(realizationReport.summary.realizedOperatorCount, 3);
assert.strictEqual(realizationReport.summary.supportRecordCount, 3);
assert.strictEqual(realizationReport.support.supportCellCount, 3);
assert.strictEqual(realizationReport.formalization.targetId, "physical_realization_map_gate");
assert.strictEqual(realizationReport.operators[0].operatorClass, "group actuation operator");
assert.strictEqual(realizationReport.operators[2].operatorClass, "graph deletion operator");
assert.strictEqual(realizationReport.operators[2].stateEffect.removedDeltaCount, 1);
assert.strictEqual(realizationReport.contactGraph.contactGraphConsistent, true);
const realizationCsv = RAD.exportPhysicalRealizationMapCsv(realizationReport);
assert.ok(realizationCsv.includes("hardware_channel"), "realization map CSV should expose mechanism channels");
assert.ok(realizationCsv.includes("graph deletion operator"), "realization map CSV should include graph deletion rows");
assert.strictEqual(
  JSON.parse(RAD.exportPhysicalRealizationMap(realizationState, { eventSequence: realizationEvents })).summary.physicalRealizationMapReady,
  true
);
const noEffectRealization = RAD.physicalRealizationMapReport(realizationState, {
  eventSequence: [RAD.releaseEvent({ r: 0, c: 0 })],
});
assert.strictEqual(noEffectRealization.summary.physicalRealizationMapReady, false);
assert.ok(
  noEffectRealization.summary.missingEvidence.includes("stateEffectRecords"),
  "operators with no state effect should fail the realization map gate"
);
const externalAudit = RAD.externalPhysicsEngineAuditReport(realizationState, {
  eventSequence: realizationEvents,
  engineAvailability: { mujoco: true },
  store: true,
});
assert.strictEqual(externalAudit.schema, "rad-sim.external-physics-engine-audit.v1");
assert.strictEqual(externalAudit.summary.externalPhysicsEngineAuditReady, true);
assert.strictEqual(externalAudit.summary.availableEngineCount, 1);
assert.ok(externalAudit.summary.feasibleEngineCount >= 1, "audit should count an available independent engine");
assert.ok(externalAudit.summary.contactModelRecordCount > 0, "audit should count contact-model records");
assert.strictEqual(externalAudit.formalization.targetId, "external_physics_engine_audit_gate");
assert.strictEqual(realizationState.experiment.externalPhysicsEngineAudit.schema, externalAudit.schema);
const externalAuditCsv = RAD.exportExternalPhysicsEngineAuditCsv(externalAudit);
assert.ok(externalAuditCsv.includes("MuJoCo"), "external audit CSV should include MuJoCo candidate");
assert.strictEqual(
  JSON.parse(RAD.exportExternalPhysicsEngineAudit(realizationState, {
    eventSequence: realizationEvents,
    engineAvailability: { mujoco: true },
  })).summary.externalPhysicsEngineAuditReady,
  true
);
const missingExternalAudit = RAD.externalPhysicsEngineAuditReport(realizationState, {
  eventSequence: realizationEvents,
  engineAvailability: { mujoco: false, pybullet: false, pychrono: false },
});
assert.strictEqual(missingExternalAudit.summary.externalPhysicsEngineAuditReady, false);
assert.ok(
  missingExternalAudit.summary.missingEvidence.includes("availableExternalEngine"),
  "external audit should fail without an available independent engine"
);
assert.ok(
  missingExternalAudit.summary.missingEvidence.includes("independentToolRecords"),
  "external audit should fail without feasible independent tool records"
);
const mujocoExport = RAD.mujocoModelExportReport(realizationState, {
  eventSequence: realizationEvents,
  fixedCells: [{ row: 0, col: 0 }],
  externalForces: [{ cell: { row: 1, col: 1 }, force: [0, 0, -0.25] }],
  store: true,
});
assert.strictEqual(mujocoExport.schema, "rad-sim.mujoco-model-export.v1");
assert.strictEqual(mujocoExport.summary.mujocoModelExportReady, true);
assert.ok(mujocoExport.summary.bodyRecordCount > 0, "MuJoCo export should include active body records");
assert.strictEqual(mujocoExport.summary.fixedBodyCount, 1);
assert.strictEqual(mujocoExport.summary.loadRecordCount, 1);
assert.strictEqual(mujocoExport.summary.pinRecordCount, 12);
assert.strictEqual(mujocoExport.summary.holeRecordCount, 12);
assert.strictEqual(mujocoExport.summary.clearanceRecordCount, 12);
assert.strictEqual(mujocoExport.summary.contactPairRecordCount, 12);
assert.strictEqual(mujocoExport.summary.contactParameterRecordCount, 12);
assert.strictEqual(mujocoExport.summary.frictionRecordCount, 12);
assert.strictEqual(mujocoExport.summary.solverParameterRecordCount, 12);
assert.strictEqual(mujocoExport.formalization.targetId, "mujoco_model_export_gate");
assert.ok(mujocoExport.xml.includes("<mujoco"), "MuJoCo export should include MJCF XML");
assert.ok(mujocoExport.xml.includes("rad_cell_0_0_nw_pin"), "MuJoCo export should include pin proxy geometry");
assert.ok(mujocoExport.xml.includes("rad_cell_0_0_nw_hole_clearance"), "MuJoCo export should include clearance proxy geometry");
assert.ok(mujocoExport.xml.includes('friction="0.4 0.02 0.001"'), "MuJoCo export should include contact friction attributes");
assert.ok(RAD.exportMujocoModelXml(realizationState, { eventSequence: realizationEvents }).includes("rad_cell_0_0"), "MuJoCo XML export should include body names");
assert.strictEqual(JSON.parse(RAD.exportMujocoModelReport(realizationState, { eventSequence: realizationEvents })).schema, mujocoExport.schema);
const mujocoContactGeometry = RAD.mujocoPinHoleContactGeometryReport(realizationState, { exportReport: mujocoExport });
assert.strictEqual(mujocoContactGeometry.schema, "rad-sim.mujoco-pin-hole-contact-geometry.v1");
assert.strictEqual(mujocoContactGeometry.summary.mujocoPinHoleContactGeometryReady, true);
assert.strictEqual(mujocoContactGeometry.summary.pinRecordCount, 12);
assert.strictEqual(mujocoContactGeometry.summary.holeRecordCount, 12);
assert.strictEqual(mujocoContactGeometry.summary.activeContactPairCount, 12);
assert.strictEqual(mujocoContactGeometry.summary.contactParameterRecordCount, 12);
assert.ok(mujocoContactGeometry.summary.minClearance > 0, "pin-hole contact geometry should preserve positive clearance");
assert.strictEqual(mujocoContactGeometry.formalization.targetId, "mujoco_pin_hole_contact_geometry_gate");
assert.ok(RAD.exportMujocoPinHoleContactGeometryCsv(mujocoContactGeometry).includes("clearance"), "contact geometry CSV should include clearance fields");
assert.strictEqual(
  JSON.parse(RAD.exportMujocoPinHoleContactGeometry(realizationState, { exportReport: mujocoExport })).summary.pinRecordCount,
  12
);
const mujocoParameters = RAD.mujocoContactParameterReport(realizationState, {
  exportReport: mujocoExport,
  contactStiffness: 1500,
  contactDamping: 3,
  friction: [0.6, 0.03, 0.002],
  calibratedContact: true,
});
assert.strictEqual(mujocoParameters.schema, "rad-sim.mujoco-contact-parameter-profile.v1");
assert.strictEqual(mujocoParameters.summary.mujocoContactParameterProfileReady, true);
assert.strictEqual(mujocoParameters.summary.parameterRecordCount, 12);
assert.strictEqual(mujocoParameters.summary.calibratedRecordCount, 12);
assert.strictEqual(mujocoParameters.formalization.targetId, "mujoco_contact_parameter_profile_gate");
assert.strictEqual(mujocoParameters.parameters[0].contactStiffness, 1500);
assert.strictEqual(JSON.stringify(mujocoParameters.parameters[0].friction), JSON.stringify([0.6, 0.03, 0.002]));
assert.ok(RAD.exportMujocoContactParameterCsv(mujocoParameters).includes("contact_stiffness"), "contact parameter CSV should include stiffness fields");
assert.strictEqual(
  JSON.parse(RAD.exportMujocoContactParameter(realizationState, {
    exportReport: mujocoExport,
    contactStiffness: 1500,
    contactDamping: 3,
    friction: [0.6, 0.03, 0.002],
    calibratedContact: true,
  })).summary.parameterRecordCount,
  12
);
const contactCalibrationPacket = RAD.contactParameterCalibrationPacket(realizationState, {
  exportReport: mujocoExport,
  contactParameterReport: mujocoParameters,
  repeatCount: 2,
  fitDatasetId: "fit-contact",
  holdoutDatasetId: "holdout-contact",
});
assert.strictEqual(contactCalibrationPacket.schema, "rad-sim.contact-parameter-calibration-packet.v1");
assert.strictEqual(contactCalibrationPacket.summary.contactParameterCalibrationPacketReady, true);
assert.strictEqual(contactCalibrationPacket.summary.contactPairRecordCount, 12);
assert.strictEqual(contactCalibrationPacket.summary.fitTemplateRowCount, 24);
assert.strictEqual(contactCalibrationPacket.summary.holdoutTemplateRowCount, 24);
assert.ok(
  contactCalibrationPacket.measurementColumns.includes("measured_pin_hole_slip_mm"),
  "contact calibration packet should request slip measurements"
);
assert.strictEqual(contactCalibrationPacket.formalization.targetId, "contact_parameter_calibration_packet_completeness");
assert.ok(
  RAD.exportContactParameterCalibrationPacketCsv(contactCalibrationPacket).includes("fitted_contact_stiffness"),
  "contact calibration CSV should include fitted stiffness fields"
);
assert.strictEqual(
  JSON.parse(RAD.exportContactParameterCalibrationPacket(realizationState, {
    exportReport: mujocoExport,
    contactParameterReport: mujocoParameters,
    repeatCount: 2,
  })).summary.fitTemplateRowCount,
  24
);
const blankContactResults = RAD.contactParameterCalibrationResultsTemplate(contactCalibrationPacket);
assert.strictEqual(blankContactResults.schema, "rad-sim.contact-parameter-calibration-results.v1");
assert.strictEqual(blankContactResults.measurements.length, 48);
assert.strictEqual(
  RAD.contactParameterCalibrationResultsFromJson(JSON.stringify(blankContactResults)).schema,
  "rad-sim.contact-parameter-calibration-results.v1"
);
assert.strictEqual(
  JSON.parse(RAD.exportContactParameterCalibrationResultsTemplate(contactCalibrationPacket)).sourcePacketSchema,
  "rad-sim.contact-parameter-calibration-packet.v1"
);
const blankContactValidation = RAD.compareContactParameterCalibrationResults(contactCalibrationPacket, blankContactResults);
assert.strictEqual(blankContactValidation.summary.contactParameterBenchValidationPass, false);
assert.ok(blankContactValidation.metrics.missingMeasurementCount > 0, "blank contact validation should require filled measurements");
const blankIntervalCalibration = RAD.contactParameterIntervalCalibrationReport(contactCalibrationPacket, blankContactResults, {
  benchValidation: blankContactValidation,
});
assert.strictEqual(blankIntervalCalibration.summary.contactParameterIntervalCalibrationReady, false);
assert.ok(
  blankIntervalCalibration.summary.missingEvidence.includes("contactParameterBenchValidationPass"),
  "blank interval calibration should depend on a passing bench validation"
);
const filledContactResults = JSON.parse(JSON.stringify(blankContactResults));
for (const resultRow of filledContactResults.measurements) {
  const measured = resultRow.requiredMeasurements;
  measured.pinRadiusMm = contactCalibrationPacket.grid.pinRadius;
  measured.holeRadiusMm = contactCalibrationPacket.grid.holeRadius;
  measured.clearanceMm = contactCalibrationPacket.grid.pinHoleClearance;
  measured.normalLoadN = 1.0;
  measured.tangentialLoadN = 0.2;
  measured.imposedZMm = 0.1;
  measured.measuredPinHoleSlipMm = 0.02;
  measured.measuredNormalForceN = 1.0;
  measured.measuredTangentForceN = 0.6;
  measured.measuredReboundRatio = 0.5;
  measured.measuredContactDurationS = 0.1;
  measured.measuredStaticFrictionCoeff = resultRow.assumedFriction[0];
  measured.measuredDynamicFrictionCoeff = resultRow.assumedFriction[1];
  measured.fittedContactStiffness = resultRow.assumedContactStiffness;
  measured.fittedContactDamping = resultRow.assumedContactDamping;
  measured.fittedSolrefTimeconst = resultRow.assumedSolref[0];
  measured.fittedSolrefDampingRatio = resultRow.assumedSolref[1];
  measured.fittedSolimpWidth = resultRow.assumedSolimp[0];
  measured.fixtureNotes = "synthetic validation row";
}
const contactBenchValidation = RAD.compareContactParameterCalibrationResults(contactCalibrationPacket, filledContactResults, {
  tolerance: 1e-9,
  holdoutTolerance: 1e-9,
});
assert.strictEqual(contactBenchValidation.schema, "rad-sim.contact-parameter-bench-validation.v1");
assert.strictEqual(contactBenchValidation.summary.contactParameterBenchValidationPass, true);
assert.strictEqual(contactBenchValidation.summary.fitParameterPass, true);
assert.strictEqual(contactBenchValidation.summary.holdoutParameterPass, true);
assert.strictEqual(contactBenchValidation.summary.independentHoldoutPass, true);
assert.strictEqual(contactBenchValidation.metrics.missingMeasurementCount, 0);
assert.strictEqual(contactBenchValidation.formalization.targetId, "contact_parameter_bench_validation_gate");
assert.ok(
  RAD.exportContactParameterBenchValidationCsv(contactBenchValidation).includes("max_parameter_residual"),
  "contact bench validation CSV should expose residual fields"
);
assert.strictEqual(
  JSON.parse(RAD.exportContactParameterBenchValidation(contactCalibrationPacket, filledContactResults, {
    tolerance: 1e-9,
    holdoutTolerance: 1e-9,
  })).summary.contactParameterBenchValidationPass,
  true
);
const contactIntervalCalibration = RAD.contactParameterIntervalCalibrationReport(contactCalibrationPacket, filledContactResults, {
  benchValidation: contactBenchValidation,
});
assert.strictEqual(contactIntervalCalibration.schema, "rad-sim.contact-parameter-interval-calibration.v1");
assert.strictEqual(contactIntervalCalibration.summary.contactParameterIntervalCalibrationReady, true);
assert.strictEqual(contactIntervalCalibration.summary.parameterIntervalCount, 10);
assert.strictEqual(contactIntervalCalibration.summary.acceptedParameterIntervalCount, 10);
assert.strictEqual(contactIntervalCalibration.summary.simulatorParametersInsideBounds, true);
assert.strictEqual(contactIntervalCalibration.formalization.targetId, "contact_parameter_interval_calibration_gate");
assert.ok(
  RAD.exportContactParameterIntervalCalibrationCsv(contactIntervalCalibration).includes("lower_bound"),
  "contact interval calibration CSV should expose interval bounds"
);
assert.strictEqual(
  JSON.parse(RAD.exportContactParameterIntervalCalibration(contactCalibrationPacket, filledContactResults, {
    benchValidation: contactBenchValidation,
  })).summary.contactParameterIntervalCalibrationReady,
  true
);
const browserMujocoRun = RAD.mujocoExternalRunReport(realizationState, { exportReport: mujocoExport });
assert.strictEqual(browserMujocoRun.schema, "rad-sim.mujoco-external-run.v1");
assert.strictEqual(browserMujocoRun.summary.mujocoExternalRunComplete, false);
assert.ok(browserMujocoRun.summary.missingEvidence.includes("browserCannotExecuteMujoco"), "browser run report should not claim MuJoCo execution");
assert.strictEqual(JSON.parse(RAD.exportMujocoExternalRun(realizationState, { exportReport: mujocoExport })).summary.mujocoExternalRunComplete, false);
const realizationFinalState = RAD.applyEventSequence(realizationState, realizationEvents);
const realizationFinalSim = RAD.simulate(realizationFinalState);
const syntheticMujocoBodies = mujocoExport.bodies.map((body) => {
  const center = realizationFinalSim.centers[body.row][body.col];
  return {
    row: body.row,
    col: body.col,
    name: body.name,
    initialPosition: body.position,
    finalPosition: [center.x, center.y, center.z],
    displacement: [0, 0, 0],
    fixed: body.fixed,
  };
});
const syntheticMujocoRun = {
  schema: "rad-sim.mujoco-external-run.v1",
  summary: {
    mujocoExternalRunComplete: true,
    bodyResultCount: syntheticMujocoBodies.length,
    expectedBodyResultCount: syntheticMujocoBodies.length,
    missingEvidence: [],
  },
  results: { bodies: syntheticMujocoBodies },
  modelExport: mujocoExport,
};
const mujocoComparison = RAD.mujocoExternalComparisonReport(realizationState, syntheticMujocoRun, {
  eventSequence: realizationEvents,
  tolerance: 1e-9,
});
assert.strictEqual(mujocoComparison.schema, "rad-sim.mujoco-external-comparison.v1");
assert.strictEqual(mujocoComparison.summary.mujocoExternalComparisonReady, true);
assert.strictEqual(mujocoComparison.summary.maxPositionError, 0);
assert.strictEqual(mujocoComparison.formalization.targetId, "mujoco_external_comparison_gate");
assert.ok(RAD.exportMujocoExternalComparisonCsv(mujocoComparison).includes("position_error"), "MuJoCo comparison CSV should expose position errors");
assert.strictEqual(
  JSON.parse(RAD.exportMujocoExternalComparison(realizationState, syntheticMujocoRun, { eventSequence: realizationEvents })).summary.mujocoExternalComparisonReady,
  true
);
const missingMujocoComparison = RAD.mujocoExternalComparisonReport(realizationState, browserMujocoRun, {
  eventSequence: realizationEvents,
});
assert.strictEqual(missingMujocoComparison.summary.mujocoExternalComparisonReady, false);
assert.ok(missingMujocoComparison.summary.missingEvidence.includes("mujocoExternalRun"), "comparison should fail without a complete external run");
const equilibriumReport = RAD.equilibriumRelationReport(realizationState, {
  eventSequence: realizationEvents,
  tolerance: 1e-8,
  residualTolerance: 1.0,
});
assert.strictEqual(equilibriumReport.schema, "rad-sim.equilibrium-relation.v1");
assert.strictEqual(equilibriumReport.summary.equilibriumRelationReady, true);
assert.strictEqual(equilibriumReport.solver.success, true);
assert.strictEqual(equilibriumReport.residual.passesTolerance, true);
assert.strictEqual(
  equilibriumReport.energy.nonnegativeTermCount,
  equilibriumReport.energy.requiredNonnegativeTermCount
);
assert.strictEqual(equilibriumReport.formalization.targetId, "equilibrium_relation_gate");
const equilibriumCsv = RAD.exportEquilibriumRelationCsv(equilibriumReport);
assert.ok(equilibriumCsv.includes("equilibrium_relation_ready"), "equilibrium CSV should expose readiness");
assert.ok(equilibriumCsv.includes("stored_energy"), "equilibrium CSV should expose energy terms");
assert.strictEqual(
  JSON.parse(RAD.exportEquilibriumRelation(realizationState, {
    eventSequence: realizationEvents,
    tolerance: 1e-8,
    residualTolerance: 1.0,
  })).summary.equilibriumRelationReady,
  true
);
const missingRealizationEquilibrium = RAD.equilibriumRelationReport(realizationState, {
  eventSequence: realizationEvents,
  realizationReport: {
    schema: "rad-sim.physical-realization-map.v1",
    summary: { physicalRealizationMapReady: false },
  },
  residualTolerance: 1.0,
});
assert.strictEqual(missingRealizationEquilibrium.summary.equilibriumRelationReady, false);
assert.ok(
  missingRealizationEquilibrium.summary.missingEvidence.includes("physicalRealizationMap"),
  "equilibrium relation should fail without ready realization evidence"
);
const reachableEquilibriumState = RAD.createState(1, 3);
reachableEquilibriumState.grid.backlash = 0;
reachableEquilibriumState.grid.couplingGain = 1;
reachableEquilibriumState.grid.zCouplingGain = 1;
reachableEquilibriumState.grid.pinRadius = 0.10;
reachableEquilibriumState.grid.holeRadius = 0.12;
RAD.clearCommands(reachableEquilibriumState);
const reachableEquilibriumEvents = [
  RAD.groupActuationEvent([{ r: 0, c: 0 }], 0.2, 0.05),
  RAD.removeCellEvent({ r: 0, c: 1 }),
];
const reachableEquilibrium = RAD.reachableEquilibriumControllabilityReport(reachableEquilibriumState, {
  eventSequence: reachableEquilibriumEvents,
  targetCells: [{ row: 0, col: 2 }],
  residualTolerance: 1.0,
});
assert.strictEqual(reachableEquilibrium.schema, "rad-sim.reachable-equilibrium-controllability.v1");
assert.strictEqual(reachableEquilibrium.summary.reachableEquilibriumControllabilityReady, true);
assert.strictEqual(reachableEquilibrium.summary.targetFullyReachable, false);
assert.strictEqual(reachableEquilibrium.actuatorBasis.cellCount, 1);
assert.strictEqual(reachableEquilibrium.targets.targetCellCount, 1);
assert.strictEqual(reachableEquilibrium.topology.topologyBlockedTargetCells, 1);
assert.strictEqual(reachableEquilibrium.formalization.targetId, "reachable_equilibrium_controllability_gate");
const reachableEquilibriumCsv = RAD.exportReachableEquilibriumControllabilityCsv(reachableEquilibrium);
assert.ok(
  reachableEquilibriumCsv.includes("reachable_equilibrium_controllability_ready"),
  "reachable equilibrium CSV should expose readiness"
);
assert.ok(
  reachableEquilibriumCsv.includes("topology_blocked_target_cells"),
  "reachable equilibrium CSV should expose topology-blocked targets"
);
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumControllability(reachableEquilibriumState, {
    eventSequence: reachableEquilibriumEvents,
    targetCells: [{ row: 0, col: 2 }],
    residualTolerance: 1.0,
  })).summary.reachableEquilibriumControllabilityReady,
  true
);
const strictReachableEquilibrium = RAD.reachableEquilibriumControllabilityReport(reachableEquilibriumState, {
  eventSequence: reachableEquilibriumEvents,
  targetCells: [{ row: 0, col: 2 }],
  requireFullTargetReachability: true,
  residualTolerance: 1.0,
});
assert.strictEqual(strictReachableEquilibrium.summary.reachableEquilibriumControllabilityReady, false);
assert.ok(
  strictReachableEquilibrium.summary.missingEvidence.includes("targetReachability"),
  "strict reachable equilibrium should fail when the target is topology-blocked"
);
assert.strictEqual(reachableEquilibrium.response.alphaReachableMap[0][0], true);
assert.strictEqual(reachableEquilibrium.response.heightReachableMap[0][0], true);
const reachableBenchProtocol = RAD.reachableEquilibriumBenchProtocol(reachableEquilibriumState, {
  controllabilityReport: reachableEquilibrium,
  repeatCount: 2,
});
assert.strictEqual(reachableBenchProtocol.schema, "rad-sim.reachable-equilibrium-bench-protocol.v1");
assert.strictEqual(reachableBenchProtocol.summary.benchProtocolReady, true);
assert.strictEqual(reachableBenchProtocol.summary.topologyBlockedControlCount, 1);
assert.strictEqual(reachableBenchProtocol.topology.topologyPolicySatisfied, true);
assert.strictEqual(reachableBenchProtocol.targetSummaries[0].blockedByTopology, true);
assert.ok(
  reachableBenchProtocol.steps.some((step) => step.id === "topology_blocked_target_control"),
  "reachable bench protocol should include a topology-blocked target control"
);
assert.strictEqual(reachableBenchProtocol.formalization.targetId, "reachable_equilibrium_bench_protocol_gate");
const reachableBenchCsv = RAD.exportReachableEquilibriumBenchProtocolCsv(reachableBenchProtocol);
assert.ok(reachableBenchCsv.includes("protocol_ready"), "reachable bench CSV should expose readiness");
assert.ok(reachableBenchCsv.includes("topology_blocked_target_control"), "reachable bench CSV should expose blocked controls");
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumBenchProtocol(reachableEquilibriumState, {
    controllabilityReport: reachableEquilibrium,
    repeatCount: 2,
  })).summary.benchProtocolReady,
  true
);
const reachableBenchResults = RAD.reachableEquilibriumBenchResultsTemplate(reachableBenchProtocol, {
  datasetId: "blocked-target-run",
});
assert.strictEqual(reachableBenchResults.schema, "rad-sim.reachable-equilibrium-bench-results.v1");
assert.ok(reachableBenchResults.measurements.length > 0, "reachable bench template should contain fillable rows");
for (const row of reachableBenchResults.measurements) {
  const modes = new Set(row.measurementMode);
  for (const measurement of row.targetMeasurements) {
    measurement.measuredTopologyComponentLabel = measurement.expectedTopologyComponentLabel;
    measurement.measuredAlphaDelta = modes.has("alpha") && measurement.predictedAlphaReachable ? 0.02 : 0;
    measurement.measuredHeightDelta = modes.has("height") && measurement.predictedHeightReachable ? 0.02 : 0;
  }
}
const parsedReachableBenchResults = RAD.reachableEquilibriumBenchResultsFromJson(JSON.stringify(reachableBenchResults));
const reachableBenchComparison = RAD.compareReachableEquilibriumBenchResults(reachableBenchProtocol, parsedReachableBenchResults, {
  responseTolerance: 1e-6,
  blockedTolerance: 1e-6,
  groupSequenceTolerance: 1e-6,
});
assert.strictEqual(reachableBenchComparison.schema, "rad-sim.reachable-equilibrium-bench-comparison.v1");
assert.strictEqual(reachableBenchComparison.summary.benchComparisonPass, true);
assert.strictEqual(reachableBenchComparison.summary.topologyBlockedLeakagePass, true);
assert.strictEqual(reachableBenchComparison.summary.groupSequencePass, true);
assert.strictEqual(reachableBenchComparison.metrics.topologyLeakageCount, 0);
assert.strictEqual(reachableBenchComparison.formalization.targetId, "reachable_equilibrium_bench_validation_gate");
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumBenchResultsTemplate(reachableBenchProtocol, {
    datasetId: "blocked-target-run",
  })).schema,
  "rad-sim.reachable-equilibrium-bench-results.v1"
);
const reachableBenchComparisonCsv = RAD.exportReachableEquilibriumBenchComparisonCsv(reachableBenchComparison);
assert.ok(reachableBenchComparisonCsv.includes("comparison_pass"), "reachable bench comparison CSV should expose pass status");
assert.ok(reachableBenchComparisonCsv.includes("group_target_reachability_probe"), "reachable bench comparison CSV should expose group rows");
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumBenchComparison(reachableBenchProtocol, parsedReachableBenchResults, {
    responseTolerance: 1e-6,
    blockedTolerance: 1e-6,
    groupSequenceTolerance: 1e-6,
  })).summary.benchComparisonPass,
  true
);
const reachableAmplitude = RAD.reachableEquilibriumAmplitudeCalibrationReport(reachableBenchProtocol, parsedReachableBenchResults, {
  comparisonReport: reachableBenchComparison,
  responseTolerance: 1e-6,
  blockedTolerance: 1e-6,
  groupSequenceTolerance: 1e-6,
  confidenceSigma: 2.0,
});
assert.strictEqual(reachableAmplitude.schema, "rad-sim.reachable-equilibrium-amplitude-calibration.v1");
assert.strictEqual(reachableAmplitude.summary.amplitudeCalibrationReady, true);
assert.strictEqual(reachableAmplitude.summary.topologyLeakageBandsPass, true);
assert.strictEqual(reachableAmplitude.summary.groupSequenceResidualPass, true);
assert.ok(reachableAmplitude.metrics.amplitudeEstimateCount > 0, "amplitude calibration should include estimates");
assert.ok(reachableAmplitude.metrics.repeatedTrialGroupCount > 0, "amplitude calibration should include repeated trials");
assert.strictEqual(reachableAmplitude.metrics.failedTopologyBandCount, 0);
assert.strictEqual(reachableAmplitude.formalization.targetId, "reachable_equilibrium_amplitude_calibration_gate");
const reachableAmplitudeCsv = RAD.exportReachableEquilibriumAmplitudeCalibrationCsv(reachableAmplitude);
assert.ok(reachableAmplitudeCsv.includes("calibration_ready"), "amplitude calibration CSV should expose readiness");
assert.ok(reachableAmplitudeCsv.includes("alpha_uncertainty"), "amplitude calibration CSV should expose uncertainty");
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumAmplitudeCalibration(reachableBenchProtocol, parsedReachableBenchResults, {
    comparisonReport: reachableBenchComparison,
    responseTolerance: 1e-6,
    blockedTolerance: 1e-6,
    groupSequenceTolerance: 1e-6,
    confidenceSigma: 2.0,
  })).summary.amplitudeCalibrationReady,
  true
);
const profileReachableEquilibrium = RAD.reachableEquilibriumControllabilityReport(reachableEquilibriumState, {
  eventSequence: [RAD.groupActuationEvent([{ r: 0, c: 0 }], 0.2, 0.05)],
  targetCells: [{ row: 0, col: 0 }],
  residualTolerance: 1.0,
});
const profileBenchProtocol = RAD.reachableEquilibriumBenchProtocol(reachableEquilibriumState, {
  controllabilityReport: profileReachableEquilibrium,
  repeatCount: 2,
});
const profileBenchResults = RAD.reachableEquilibriumBenchResultsTemplate(profileBenchProtocol, {
  datasetId: "reachable-target-run",
});
for (const row of profileBenchResults.measurements) {
  const modes = new Set(row.measurementMode);
  for (const measurement of row.targetMeasurements) {
    measurement.measuredTopologyComponentLabel = measurement.expectedTopologyComponentLabel;
    measurement.measuredAlphaDelta = modes.has("alpha") && measurement.predictedAlphaReachable ? 0.02 : 0;
    measurement.measuredHeightDelta = modes.has("height") && measurement.predictedHeightReachable ? 0.02 : 0;
  }
}
const profileBenchComparison = RAD.compareReachableEquilibriumBenchResults(profileBenchProtocol, profileBenchResults, {
  responseTolerance: 1e-6,
  blockedTolerance: 1e-6,
  groupSequenceTolerance: 1e-6,
});
const profileAmplitude = RAD.reachableEquilibriumAmplitudeCalibrationReport(profileBenchProtocol, profileBenchResults, {
  comparisonReport: profileBenchComparison,
  responseTolerance: 1e-6,
  blockedTolerance: 1e-6,
  groupSequenceTolerance: 1e-6,
  confidenceSigma: 2.0,
});
const reachableEmpiricalProfile = RAD.reachableEquilibriumEmpiricalProfileFromAmplitude(profileAmplitude, {
  safetyFactor: 1.5,
  minSamples: 1,
});
assert.strictEqual(reachableEmpiricalProfile.schema, "rad-sim.reachable-equilibrium-empirical-profile.v1");
assert.strictEqual(reachableEmpiricalProfile.summary.empiricalProfileReady, true);
assert.ok(reachableEmpiricalProfile.metrics.safeProposalCount > 0, "empirical profile should include safe proposals");
assert.strictEqual(reachableEmpiricalProfile.holdoutValidation.holdoutPass, true);
assert.ok(
  reachableEmpiricalProfile.recommendedUpdates.some((update) => update.name === "alphaResponseScale"),
  "empirical profile should include alpha response scale"
);
assert.strictEqual(reachableEmpiricalProfile.formalization.targetId, "reachable_equilibrium_empirical_profile_gate");
const reachableEmpiricalProfileCsv = RAD.exportReachableEquilibriumEmpiricalProfileCsv(reachableEmpiricalProfile);
assert.ok(reachableEmpiricalProfileCsv.includes("profile_ready"), "empirical profile CSV should expose readiness");
assert.ok(reachableEmpiricalProfileCsv.includes("alphaResponseScale"), "empirical profile CSV should expose alpha profile rows");
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumEmpiricalProfile(profileAmplitude, {
    safetyFactor: 1.5,
    minSamples: 1,
  })).summary.empiricalProfileReady,
  true
);
const profileInverseState = RAD.createState(1, 3);
RAD.clearCommands(profileInverseState);
profileInverseState.target.type = "custom";
profileInverseState.target.amplitude = 0.03;
profileInverseState.target.customExpression = "r==0&&c==0?amplitude:0";
profileInverseState.cells.commandZ[0][0] = 0.03;
const profileInverse = RAD.reachableEquilibriumProfileInverseReport(profileInverseState, reachableEmpiricalProfile, {
  store: true,
});
assert.strictEqual(profileInverse.schema, "rad-sim.reachable-equilibrium-profile-inverse.v1");
assert.strictEqual(profileInverse.summary.profileInverseReady, true);
assert.strictEqual(profileInverse.target.targetCellCount, 1);
assert.strictEqual(profileInverse.formalization.targetId, "reachable_equilibrium_profile_inverse_gate");
assert.strictEqual(profileInverseState.experiment.reachableEquilibriumProfileInverse.schema, profileInverse.schema);
const profileInverseCsv = RAD.exportReachableEquilibriumProfileInverseCsv(profileInverse);
assert.ok(profileInverseCsv.includes("profile_inverse_ready"), "profile-aware inverse CSV should expose readiness");
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumProfileInverse(profileInverseState, reachableEmpiricalProfile)).summary.profileInverseReady,
  true
);
const profileInverseAcceptance = RAD.reachableEquilibriumProfileInverseAcceptanceReport(profileInverse, {
  maxWeightedResidualScore: 1.0,
  maxBandFailures: 1,
  maxActiveActuators: 1,
});
assert.strictEqual(profileInverseAcceptance.schema, "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1");
assert.strictEqual(profileInverseAcceptance.summary.profileInverseAcceptanceReady, true);
assert.strictEqual(profileInverseAcceptance.decision.decision, "accept-for-preview");
assert.strictEqual(profileInverseAcceptance.formalization.targetId, "reachable_equilibrium_profile_inverse_acceptance_gate");
assert.ok(
  RAD.exportReachableEquilibriumProfileInverseAcceptanceCsv(profileInverseAcceptance).includes("acceptance_ready"),
  "profile-aware inverse acceptance CSV should expose readiness"
);
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumProfileInverseAcceptance(profileInverse, {
    maxWeightedResidualScore: 1.0,
    maxBandFailures: 1,
    maxActiveActuators: 1,
  })).summary.profileInverseAcceptanceReady,
  true
);
const reviewProfileInverse = JSON.parse(JSON.stringify(profileInverse));
reviewProfileInverse.summary.profileWeightedResidualScore = 5.0;
const reviewProfileInverseAcceptance = RAD.reachableEquilibriumProfileInverseAcceptanceReport(reviewProfileInverse, {
  maxWeightedResidualScore: 1.0,
});
assert.strictEqual(reviewProfileInverseAcceptance.summary.profileInverseAcceptanceReady, false);
assert.strictEqual(reviewProfileInverseAcceptance.decision.decision, "review-required");
assert.ok(reviewProfileInverseAcceptance.decision.failedCriteria.includes("profileWeightedResidualScore"));
const profileInversePreviewPacket = RAD.reachableEquilibriumProfileInversePreviewPacket(
  profileInverseState,
  profileInverse,
  profileInverseAcceptance,
  { store: true, packetId: "browser-profile-inverse-preview" }
);
assert.strictEqual(profileInversePreviewPacket.schema, "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1");
assert.strictEqual(profileInversePreviewPacket.summary.profileInversePreviewPacketReady, true);
assert.strictEqual(profileInversePreviewPacket.summary.commandCount, 1);
assert.strictEqual(profileInversePreviewPacket.formalization.targetId, "reachable_equilibrium_profile_inverse_preview_packet_gate");
assert.strictEqual(profileInverseState.experiment.reachableEquilibriumProfileInversePreviewPacket.schema, profileInversePreviewPacket.schema);
assert.ok(
  RAD.exportReachableEquilibriumProfileInversePreviewPacketCsv(profileInversePreviewPacket).includes("packet_ready"),
  "profile-aware inverse preview packet CSV should expose readiness"
);
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumProfileInversePreviewPacket(profileInverseState, profileInverse, profileInverseAcceptance)).summary.profileInversePreviewPacketReady,
  true
);
const rejectedProfileInversePreviewPacket = RAD.reachableEquilibriumProfileInversePreviewPacket(
  profileInverseState,
  profileInverse,
  reviewProfileInverseAcceptance
);
assert.strictEqual(rejectedProfileInversePreviewPacket.summary.profileInversePreviewPacketReady, false);
assert.ok(rejectedProfileInversePreviewPacket.summary.missingEvidence.includes("profileInverseAcceptanceReady"));
const profileInversePreviewReplay = RAD.reachableEquilibriumProfileInversePreviewReplayReport(
  profileInverseState,
  profileInversePreviewPacket,
  { store: true, residualAgreementTolerance: 1e-8 }
);
assert.strictEqual(profileInversePreviewReplay.schema, "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1");
assert.strictEqual(profileInversePreviewReplay.summary.profileInversePreviewReplayReady, true);
assert.strictEqual(profileInversePreviewReplay.commands.replayedCommandCount, 1);
assert.strictEqual(profileInversePreviewReplay.commands.invalidCommandCount, 0);
assert.strictEqual(profileInversePreviewReplay.residuals.residualAgreementPass, true);
assert.strictEqual(profileInversePreviewReplay.formalization.targetId, "reachable_equilibrium_profile_inverse_preview_replay_gate");
assert.strictEqual(profileInverseState.experiment.reachableEquilibriumProfileInversePreviewReplay.schema, profileInversePreviewReplay.schema);
assert.ok(
  RAD.exportReachableEquilibriumProfileInversePreviewReplayCsv(profileInversePreviewReplay).includes("replay_ready"),
  "profile-aware inverse preview replay CSV should expose readiness"
);
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumProfileInversePreviewReplay(profileInverseState, profileInversePreviewPacket, {
    residualAgreementTolerance: 1e-8,
  })).summary.profileInversePreviewReplayReady,
  true
);
const profileInversePreviewPhysical = RAD.reachableEquilibriumProfileInversePreviewPhysicalReport(
  profileInverseState,
  profileInversePreviewPacket,
  { store: true, residualAgreementTolerance: 1e-8 }
);
assert.strictEqual(profileInversePreviewPhysical.schema, "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1");
assert.strictEqual(profileInversePreviewPhysical.summary.profileInversePreviewPhysicalReady, true);
assert.strictEqual(profileInversePreviewPhysical.physical.physicalSuccess, true);
assert.strictEqual(profileInversePreviewPhysical.commands.replayedCommandCount, 1);
assert.ok(Number.isFinite(profileInversePreviewPhysical.comparison.heightRmsModelError));
assert.ok(Number.isFinite(profileInversePreviewPhysical.comparison.centerRmsModelError));
assert.strictEqual(profileInversePreviewPhysical.formalization.targetId, "reachable_equilibrium_profile_inverse_preview_physical_gate");
assert.strictEqual(profileInverseState.experiment.reachableEquilibriumProfileInversePreviewPhysical.schema, profileInversePreviewPhysical.schema);
assert.ok(
  RAD.exportReachableEquilibriumProfileInversePreviewPhysicalCsv(profileInversePreviewPhysical).includes("physical_ready"),
  "profile-aware inverse preview physical CSV should expose readiness"
);
assert.strictEqual(
  JSON.parse(RAD.exportReachableEquilibriumProfileInversePreviewPhysical(profileInverseState, profileInversePreviewPacket, {
    residualAgreementTolerance: 1e-8,
  })).summary.profileInversePreviewPhysicalReady,
  true
);
const invalidReplayPacket = JSON.parse(JSON.stringify(profileInversePreviewPacket));
invalidReplayPacket.commands[0].row = profileInverseState.grid.rows + 10;
const invalidProfileInversePreviewReplay = RAD.reachableEquilibriumProfileInversePreviewReplayReport(
  profileInverseState,
  invalidReplayPacket
);
assert.strictEqual(invalidProfileInversePreviewReplay.summary.profileInversePreviewReplayReady, false);
assert.ok(invalidProfileInversePreviewReplay.summary.missingEvidence.includes("validCommandRecords"));
const needsHoldoutProfile = RAD.reachableEquilibriumEmpiricalProfileFromAmplitude(profileAmplitude, {
  requireHoldout: true,
});
assert.strictEqual(needsHoldoutProfile.summary.empiricalProfileReady, false);
assert.ok(needsHoldoutProfile.summary.missingEvidence.includes("holdoutValidation"));
const leakedReachableBenchResults = JSON.parse(JSON.stringify(reachableBenchResults));
for (const row of leakedReachableBenchResults.measurements) {
  if (row.stepId === "topology_blocked_target_control") {
    row.targetMeasurements[0].measuredHeightDelta = 0.01;
    break;
  }
}
const leakedReachableBenchComparison = RAD.compareReachableEquilibriumBenchResults(reachableBenchProtocol, leakedReachableBenchResults, {
  blockedTolerance: 1e-6,
});
assert.strictEqual(leakedReachableBenchComparison.summary.benchComparisonPass, false);
assert.strictEqual(leakedReachableBenchComparison.summary.topologyBlockedLeakagePass, false);
assert.ok(leakedReachableBenchComparison.metrics.topologyLeakageCount > 0, "leaked blocked target should fail topology leakage");
const leakedReachableAmplitude = RAD.reachableEquilibriumAmplitudeCalibrationReport(reachableBenchProtocol, leakedReachableBenchResults, {
  comparisonReport: leakedReachableBenchComparison,
  blockedTolerance: 1e-6,
});
assert.strictEqual(leakedReachableAmplitude.summary.amplitudeCalibrationReady, false);
assert.ok(
  leakedReachableAmplitude.summary.missingEvidence.includes("topologyLeakageBand"),
  "amplitude calibration should fail when blocked leakage bands exceed tolerance"
);
assert.strictEqual(
  JSON.parse(RAD.exportCalibrationModelProfileHoldoutValidation(lowZProfileState, holdoutResults, null, { protocol: profileProtocol })).holdoutPass,
  true
);
const applicationAudit = RAD.applyCalibrationModelProfile(readyState, exportedModelProfile);
assert.strictEqual(applicationAudit.schema, "rad-sim.calibration-model-profile-application.v1");
assert.strictEqual(applicationAudit.appliedUpdates[0].configField, "zCouplingGain");
assert.strictEqual(Number(readyState.grid.zCouplingGain.toFixed(6)), Number(zProfileUpdate.proposed.toFixed(6)));
assert.strictEqual(readyState.experiment.calibrationModelProfile.schema, "rad-sim.calibration-model-profile.v1");
assert.strictEqual(readyState.experiment.calibrationModelProfileApplication.schema, "rad-sim.calibration-model-profile-application.v1");
assert.ok(Array.isArray(readyState.experiment.calibrationModelProfileHistory), "profile applications should be audited in state history");
assert.ok(readyState.experiment.calibrationModelProfileHistory.length >= 1, "profile application audit history should not be empty");
assert.strictEqual(report.summary.worstCell.row, perturbedLift.cells[0].row);
assert.strictEqual(report.comparison.fit.height.sampleCount, perturbedComparison.fit.height.sampleCount);
assert.ok(Array.isArray(report.comparison.fitResidualField.combinedError), "comparison report should include fit residual field");
assert.ok(report.summary.topCells.length > 0, "comparison report should include ranked raw cells");
assert.ok(report.summary.fitResidualTopCells.length > 0, "comparison report should include ranked residual cells");
assert.strictEqual(JSON.parse(RAD.exportCalibrationComparisonReport(readyState)).schema, report.schema);
const appliedProfileState = RAD.createState(2, 2);
appliedProfileState.grid.cellSize = 1.25;
appliedProfileState.grid.hardwareProfile = JSON.parse(JSON.stringify(state.grid.hardwareProfile));
const appliedProfile = RAD.applyHardwareProfileToGrid(appliedProfileState);
assert.strictEqual(Number(appliedProfileState.grid.backlash.toFixed(6)), 0.1);
assert.strictEqual(Number(appliedProfileState.grid.pinRadius.toFixed(6)), Number((7.56 * 1.25 / 42).toFixed(6)));
assert.strictEqual(Number(appliedProfileState.grid.holeRadius.toFixed(6)), Number((9.45 * 1.25 / 42).toFixed(6)));
assert.strictEqual(appliedProfile.measuredCount, 3);
const legacyProfileJson = JSON.parse(RAD.serialize(RAD.createState(2, 2)));
legacyProfileJson.grid.paperSideLengthMm = 44;
legacyProfileJson.grid.paperHoleToleranceMm = 0.18;
delete legacyProfileJson.grid.hardwareProfile;
const migratedLegacy = RAD.deserialize(JSON.stringify(legacyProfileJson));
assert.strictEqual(migratedLegacy.grid.hardwareProfile.sideLengthMm, 44);
assert.strictEqual(migratedLegacy.grid.hardwareProfile.fabricationHoleToleranceMm, 0.18);
assert.deepStrictEqual(JSON.parse(JSON.stringify(restored.view.camera)), state.view.camera);
assert.strictEqual(restored.view.overlayMode, "height");
assert.strictEqual(restored.view.simulationMode, "springPreview");
assert.strictEqual(restored.view.membraneVisible, false);
assert.strictEqual(restored.view.externalCaseVisible, false);
assert.deepStrictEqual(JSON.parse(JSON.stringify(restored.selection)), state.selection);
assert.strictEqual(restored.view.paintRadius, 2);
assert.strictEqual(restored.cells.commandAlpha[2][3], -0.24);
assert.strictEqual(restored.cells.commandZ[2][3], 0.31);
assert.strictEqual(restored.cells.locked[0][0], true);
assert.strictEqual(restored.cells.lockAlpha[0][0], 0.82);
assert.strictEqual(restored.cells.lockZ[0][0], 0.14);
assert.strictEqual(restored.cells.actuatorAllowed[4][5], false);
assert.strictEqual(restored.cells.removed[1][2], true);
assert.strictEqual(restored.experiment.presetName, "custom-roundtrip");
assert.strictEqual(restored.experiment.notes, "roundtrip validation");

const sim = RAD.simulate(restored);
RAD.updateDerivedCells(restored, sim);
assert.strictEqual(sim.alpha.length, 5);
assert.strictEqual(sim.alpha[0].length, 6);
assert.strictEqual(Number(restored.cells.alpha[2][3].toFixed(6)), Number(sim.alpha[2][3].toFixed(6)));
assert.strictEqual(Number(restored.cells.theta[2][3].toFixed(6)), Number(sim.theta[2][3].toFixed(6)));
assert.strictEqual(Number(restored.cells.z[2][3].toFixed(6)), Number(sim.height[2][3].toFixed(6)));
assert.ok(Number.isFinite(sim.metrics.meanAlpha));
assert.ok(Number.isFinite(sim.metrics.rmsTargetError));
assert.strictEqual(Number(sim.alpha[0][0].toFixed(6)), 0.82);
assert.strictEqual(Number(sim.height[0][0].toFixed(6)), 0.14);
assert.strictEqual(typeof RAD.pairDistanceFromTheta, "function");
assert.strictEqual(typeof RAD.solveConstraintKinematicApprox, "function");
assert.ok(RAD.pairDistanceFromTheta(restored, 10, 10) > 1);
const constraintState = RAD.deserialize(RAD.serialize(restored));
constraintState.view.simulationMode = "constraintSolved";
const constraintSim = RAD.simulate(constraintState);
assert.strictEqual(constraintSim.metrics.constraintSolved, true);
assert.ok(Number.isFinite(constraintSim.metrics.constraintMaxEdgeError));
assert.ok(Number.isFinite(constraintSim.metrics.constraintMeanEdgeError));
const activeSim = RAD.simulateActive(restored);
assert.strictEqual(activeSim.metrics.model, "spring-preview");
assert.strictEqual(activeSim.metrics.physicalPreview, true);
assert.ok(Number.isFinite(activeSim.metrics.physicalRmsHeightDelta));
assert.ok(Number.isFinite(activeSim.metrics.physicalMaxHeightDelta));
assert.ok(Number.isFinite(activeSim.metrics.physicalRmsCenterDelta));
assert.ok(Number.isFinite(activeSim.metrics.physicalMaxCenterDelta));
assert.ok(activeSim.metrics.physicalIterations > 0);
assert.strictEqual(activeSim.metrics.physicalActiveCells, restored.grid.rows * restored.grid.cols - 1);
assert.strictEqual(activeSim.metrics.physicalSkippedSpringEdges, 4);
assert.strictEqual(activeSim.metrics.physicalActiveSpringEdges, 45);
assert.strictEqual(activeSim.linkStrain.horizontal[1][1], null);
assert.strictEqual(activeSim.linkStrain.horizontal[1][2], null);
assert.strictEqual(activeSim.linkStrain.vertical[0][2], null);
assert.strictEqual(activeSim.linkStrain.vertical[1][2], null);
assert.strictEqual(activeSim.slope.magnitude[1][2], null);
assert.strictEqual(activeSim.slope.tilt[1][2], null);
assert.strictEqual(activeSim.modelErrorHeight.length, restored.grid.rows);
assert.strictEqual(activeSim.modelErrorCenter[0].length, restored.grid.cols);
assert.strictEqual(
  Number(activeSim.modelErrorHeight[2][3].toFixed(6)),
  Number((activeSim.height[2][3] - sim.height[2][3]).toFixed(6)),
  "spring preview should expose per-cell height disagreement from the kinematic state"
);
const positionFixtureState = RAD.createState(3, 3);
positionFixtureState.view.simulationMode = "springPreview";
positionFixtureState.cells.positionLocked[1][1] = true;
positionFixtureState.cells.commandAlpha[1][1] = -0.4;
positionFixtureState.cells.commandZ[1][1] = 0.5;
const positionFixtureRef = RAD.referenceCenter(positionFixtureState, 1, 1);
const positionFixtureSim = RAD.simulateActive(positionFixtureState);
assert.strictEqual(
  positionFixtureSim.metrics.physicalPositionLockedCells,
  1,
  "spring preview should count hard position fixtures"
);
assert.ok(Math.abs(positionFixtureSim.centers[1][1].x - positionFixtureRef.x) < 1e-9, "position lock should preserve x");
assert.ok(Math.abs(positionFixtureSim.centers[1][1].y - positionFixtureRef.y) < 1e-9, "position lock should preserve y");
assert.ok(Math.abs(positionFixtureSim.centers[1][1].z - positionFixtureRef.z) < 1e-9, "position lock should preserve z");
const browserMesh = RAD.buildPaperRadMesh(restored, { sim: activeSim, includePins: false });
assert.ok(browserMesh.vertexCount > 0, "browser OBJ mesh should include vertices");
assert.ok(browserMesh.faceCount > 0, "browser OBJ mesh should include faces");
assert.ok(browserMesh.components.some((part) => part.kind === "connector"), "browser OBJ mesh should include connector components");
assert.strictEqual(browserMesh.source.calibrationProfile, "bench-v1");
assert.strictEqual(browserMesh.source.measuredCalibrationFields, 3);
assert.strictEqual(Number(browserMesh.source.plateThickness.toFixed(6)), Number((2.4 * 1.25 / 42).toFixed(6)));
assert.strictEqual(Number(browserMesh.source.pinRadius.toFixed(6)), Number((7.56 * 1.25 / 42).toFixed(6)));
const browserObj = RAD.exportPaperRadMeshObj(browserMesh);
assert.ok(browserObj.includes("o cell_2_3_outer_plate"), "browser OBJ should include cell object names");
assert.ok(browserObj.includes("# kind connector"), "browser OBJ should include connector metadata");
assert.ok(browserObj.includes("# calibrationProfile bench-v1"), "browser OBJ should include calibration profile metadata");
assert.ok(browserObj.includes("# measuredCalibrationFields 3"), "browser OBJ should include measured calibration count");
assert.strictEqual((browserObj.match(/^v /gm) || []).length, browserMesh.vertexCount);
assert.strictEqual((browserObj.match(/^f /gm) || []).length, browserMesh.faceCount);

const inverseValidationState = RAD.createState(3, 3);
RAD.clearCommands(inverseValidationState);
inverseValidationState.target.type = "gaussian";
inverseValidationState.target.amplitude = 0.28;
inverseValidationState.grid.zCouplingGain = 0.0;
const linearFit = RAD.solveLinearizedTargetFit(inverseValidationState, { maxActuators: 3, maxColumns: 4 });
const physicalValidation = RAD.validateInversePlanPhysical(inverseValidationState);
assert.strictEqual(physicalValidation.strategy, "spring-preview-inverse-validation");
assert.strictEqual(physicalValidation.source, linearFit.commands.length ? "linear" : "plan");
assert.strictEqual(physicalValidation.physicalAvailable, true);
assert.ok(Number.isFinite(physicalValidation.physicalProjectedError), "physical inverse validation should report finite target error");
assert.ok(Number.isFinite(physicalValidation.centerModelRms), "physical inverse validation should report finite model disagreement");
assert.ok(physicalValidation.modelAgreementScore > 0 && physicalValidation.modelAgreementScore <= 1, "model agreement score should be normalized");
assert.deepStrictEqual(
  JSON.parse(JSON.stringify(inverseValidationState.inverse.physicalValidation)),
  JSON.parse(JSON.stringify(physicalValidation)),
  "physical inverse validation should be stored on state.inverse"
);
const inverseReport = RAD.inverseDesignReport(inverseValidationState);
assert.strictEqual(inverseReport.schema, "rad-sim.inverse-design-report.v1");
assert.strictEqual(inverseReport.inverse.linearSolution.strategy, "linearized-jacobian-greedy-fit");
assert.strictEqual(inverseReport.inverse.physicalValidation.strategy, "spring-preview-inverse-validation");
assert.ok(inverseReport.inverse.targetReachability.model, "inverse report should include target reachability diagnostics");
assert.strictEqual(JSON.parse(RAD.exportInverseDesignReport(inverseValidationState)).schema, inverseReport.schema);
const underactuatedTargetState = RAD.createState(3, 3);
RAD.clearCommands(underactuatedTargetState);
underactuatedTargetState.grid.zCouplingGain = 0;
underactuatedTargetState.target.type = "custom";
underactuatedTargetState.target.amplitude = 0.3;
underactuatedTargetState.target.customExpression = "r==2&&c==2?amplitude:0";
for (let r = 0; r < underactuatedTargetState.grid.rows; r += 1) {
  for (let c = 0; c < underactuatedTargetState.grid.cols; c += 1) {
    underactuatedTargetState.cells.actuatorAllowed[r][c] = r === 0 && c === 0;
  }
}
const underactuatedJacobian = RAD.buildResponseJacobian(underactuatedTargetState, { responseThreshold: 0.01 });
assert.strictEqual(underactuatedJacobian.targetReachability.model, "finite-response-height-reachability");
assert.strictEqual(underactuatedJacobian.targetReachability.directionalModel, "sign-compatible finite-response-height-reachability");
assert.strictEqual(underactuatedJacobian.targetReachability.targetHeightCells, 1);
assert.strictEqual(underactuatedJacobian.targetReachability.upwardTargetHeightCells, 1);
assert.strictEqual(underactuatedJacobian.targetReachability.downwardTargetHeightCells, 0);
assert.strictEqual(underactuatedJacobian.targetReachability.underactuatedHeightCells, 1);
assert.strictEqual(underactuatedJacobian.targetReachability.positiveUnderactuatedHeightCells, 1);
assert.strictEqual(underactuatedJacobian.targetReachability.negativeUnderactuatedHeightCells, 0);
assert.strictEqual(underactuatedJacobian.targetReachability.worstUnderactuatedCell.row, 2);
assert.strictEqual(underactuatedJacobian.targetReachability.worstUnderactuatedCell.col, 2);
assert.ok(underactuatedJacobian.targetReachability.maxUnreachableHeightResidual > 0, "underactuated target report should expose residual magnitude");
const underactuatedFit = RAD.solveLinearizedTargetFit(underactuatedTargetState, { maxActuators: 1, maxColumns: 2 });
assert.strictEqual(underactuatedFit.underactuatedHeightCells, 1, "linear fit should carry underactuated target count");
assert.strictEqual(underactuatedFit.targetReachability.worstUnderactuatedCell.col, 2);
const downwardUnderactuatedTargetState = RAD.createState(3, 3);
RAD.clearCommands(downwardUnderactuatedTargetState);
downwardUnderactuatedTargetState.grid.zCouplingGain = 0;
downwardUnderactuatedTargetState.target.type = "custom";
downwardUnderactuatedTargetState.target.amplitude = 0.3;
downwardUnderactuatedTargetState.target.customExpression = "r==2&&c==2?-amplitude:0";
for (let r = 0; r < downwardUnderactuatedTargetState.grid.rows; r += 1) {
  for (let c = 0; c < downwardUnderactuatedTargetState.grid.cols; c += 1) {
    downwardUnderactuatedTargetState.cells.actuatorAllowed[r][c] = r === 0 && c === 0;
  }
}
const downwardUnderactuatedJacobian = RAD.buildResponseJacobian(downwardUnderactuatedTargetState, { responseThreshold: 0.01 });
assert.strictEqual(downwardUnderactuatedJacobian.targetReachability.upwardTargetHeightCells, 0);
assert.strictEqual(downwardUnderactuatedJacobian.targetReachability.downwardTargetHeightCells, 1);
assert.strictEqual(downwardUnderactuatedJacobian.targetReachability.positiveUnderactuatedHeightCells, 0);
assert.strictEqual(downwardUnderactuatedJacobian.targetReachability.negativeUnderactuatedHeightCells, 1);

const topologyBlockedTargetState = RAD.createState(1, 3);
RAD.clearCommands(topologyBlockedTargetState);
topologyBlockedTargetState.grid.backlash = 0;
topologyBlockedTargetState.grid.couplingGain = 1;
topologyBlockedTargetState.target.type = "custom";
topologyBlockedTargetState.target.amplitude = 0.3;
topologyBlockedTargetState.target.customExpression = "r==0&&c==2?amplitude:0";
topologyBlockedTargetState.cells.removed[0][1] = true;
for (let r = 0; r < topologyBlockedTargetState.grid.rows; r += 1) {
  for (let c = 0; c < topologyBlockedTargetState.grid.cols; c += 1) {
    topologyBlockedTargetState.cells.actuatorAllowed[r][c] = r === 0 && c === 0;
  }
}
const topologyBlockedJacobian = RAD.buildResponseJacobian(topologyBlockedTargetState, { responseThreshold: 0.01 });
assert.strictEqual(topologyBlockedJacobian.targetReachability.topologyComponentCount, 2);
assert.strictEqual(topologyBlockedJacobian.targetReachability.topologyActuatedComponentCount, 1);
assert.strictEqual(topologyBlockedJacobian.targetReachability.topologyReachableMap[0][0], true);
assert.strictEqual(topologyBlockedJacobian.targetReachability.topologyReachableMap[0][1], false);
assert.strictEqual(topologyBlockedJacobian.targetReachability.topologyReachableMap[0][2], false);
assert.strictEqual(topologyBlockedJacobian.targetReachability.topologyBlockedHeightCells, 1);
assert.strictEqual(topologyBlockedJacobian.targetReachability.worstTopologyBlockedCell.row, 0);
assert.strictEqual(topologyBlockedJacobian.targetReachability.worstTopologyBlockedCell.col, 2);
assert.ok(topologyBlockedJacobian.targetReachability.maxTopologyBlockedHeightResidual > 0, "removed-cell topology should expose blocked target residual magnitude");
const topologyReport = RAD.topologyExperimentReport(topologyBlockedTargetState, [
  {
    name: "intact-left-source",
    commands: [{ r: 0, c: 0, alpha: -0.2, z: 0.2 }],
    actuatorCells: [{ r: 0, c: 0 }],
  },
  {
    name: "removed-middle-left-source",
    commands: [{ r: 0, c: 0, alpha: -0.2, z: 0.2 }],
    removedCells: [{ r: 0, c: 1 }],
    actuatorCells: [{ r: 0, c: 0 }],
  },
]);
assert.strictEqual(topologyReport.schema, "rad-sim.browser-topology-experiment-report.v1");
assert.strictEqual(topologyReport.summary.scenarioCount, 2);
assert.strictEqual(topologyReport.scenarios[1].topology.componentCount, 2);
assert.strictEqual(topologyReport.scenarios[1].metrics.componentReachableCells, 1);
assert.strictEqual(topologyReport.scenarios[1].metrics.componentBlockedCells, 1);
assert.strictEqual(topologyReport.scenarios[1].metrics.componentResponseRank.length, 2);
assert.strictEqual(topologyReport.scenarios[1].metrics.componentResponseRank[1].blockedHeightCells, 1);
assert.strictEqual(topologyReport.scenarioComparisons[0].componentCountDelta, 1);
assert.ok(topologyReport.scenarioComparisons[0].reachableHeightCellDelta < 0, "removed topology should reduce reachable height cells");
assert.ok(topologyReport.scenarios[1].metrics.maxAbsTargetHeightResidual > 0, "browser topology report should expose target residual");
assert.strictEqual(JSON.parse(RAD.exportTopologyExperimentReport(topologyBlockedTargetState)).schema, topologyReport.schema);

const operatorState = RAD.createState(3, 3);
RAD.clearCommands(operatorState);
const preLockState = RAD.applyEventSequence(operatorState, [
  RAD.localActuationEvent({ r: 1, c: 1 }, -0.3, 0.2),
]);
const expectedLockHeight = RAD.simulate(preLockState).height[1][1];
const lockedEventState = RAD.applyEventSequence(operatorState, [
  RAD.localActuationEvent({ r: 1, c: 1 }, -0.3, 0.2),
  RAD.lockEvent({ r: 1, c: 1 }),
  RAD.clearActuationEvent(),
]);
assert.strictEqual(lockedEventState.cells.locked[1][1], true);
assert.strictEqual(Number(lockedEventState.cells.lockAlpha[1][1].toFixed(6)), 0.7);
assert.strictEqual(Number(lockedEventState.cells.lockZ[1][1].toFixed(6)), Number(expectedLockHeight.toFixed(6)));
assert.strictEqual(Number(lockedEventState.cells.commandAlpha[1][1].toFixed(6)), 0);
assert.strictEqual(Number(lockedEventState.cells.commandZ[1][1].toFixed(6)), 0);
assert.strictEqual(Number(RAD.simulate(lockedEventState).alpha[1][1].toFixed(6)), 0.7);
assert.strictEqual(Number(RAD.simulate(lockedEventState).height[1][1].toFixed(6)), Number(expectedLockHeight.toFixed(6)));
assert.notStrictEqual(Number(RAD.simulate(lockedEventState).height[1][1].toFixed(6)), 0);
assert.strictEqual(RAD.finiteDieOffRadius(RAD.simulate(lockedEventState).dieOff), 0);
const removedPathState = RAD.createState(1, 3);
RAD.clearCommands(removedPathState);
removedPathState.grid.backlash = 0;
removedPathState.grid.couplingGain = 1;
const removedPathFinal = RAD.applyEventSequence(removedPathState, [
  RAD.removeCellEvent({ r: 0, c: 1 }),
  RAD.localActuationEvent({ r: 0, c: 0 }, 0.3, 0.2),
  RAD.localActuationEvent({ r: 0, c: 1 }, 0.3, 0.2),
]);
const removedPathSim = RAD.simulate(removedPathFinal);
assert.strictEqual(removedPathFinal.cells.removed[0][1], true);
assert.strictEqual(Number(removedPathFinal.cells.commandAlpha[0][1].toFixed(6)), 0);
assert.strictEqual(Number(removedPathSim.influence[0][2].toFixed(6)), 0);
assert.strictEqual(removedPathSim.dieOff[0][2], Infinity);
assert.strictEqual(removedPathSim.topology.componentCount, 2);
assert.strictEqual(removedPathSim.topology.deletedEdges, 2);
assert.strictEqual(removedPathSim.metrics.topologyComponentCount, 2);
const restoredPathFinal = RAD.applyProgrammableEvent(removedPathFinal, RAD.restoreCellEvent({ r: 0, c: 1 }));
assert.strictEqual(restoredPathFinal.cells.removed[0][1], false);
const groupState = RAD.applyEventSequence(RAD.createState(2, 2), [
  RAD.clearActuationEvent(),
  RAD.removeCellEvent({ r: 1, c: 1 }),
  RAD.groupActuationEvent([{ r: 0, c: 0 }, { r: 1, c: 0 }, { r: 1, c: 1 }], -0.2, 0.15),
]);
assert.strictEqual(Number(groupState.cells.commandAlpha[0][0].toFixed(6)), -0.2);
assert.strictEqual(Number(groupState.cells.commandZ[1][0].toFixed(6)), 0.15);
assert.strictEqual(Number(groupState.cells.commandAlpha[1][1].toFixed(6)), 0);
const groupDecomposition = RAD.compareGroupActuationDecomposition(
  groupState,
  [{ r: 0, c: 0 }, { r: 1, c: 0 }, { r: 1, c: 1 }],
  -0.2,
  0.15
);
assert.strictEqual(groupDecomposition.decomposes, true);
assert.strictEqual(groupDecomposition.cells.length, 3);
assert.strictEqual(groupDecomposition.appliedCells.length, 2);
assert.strictEqual(groupDecomposition.skippedRemovedCells.length, 1);
assert.strictEqual(groupDecomposition.commandError, 0);
assert.strictEqual(groupDecomposition.finalError, 0);
const groupRemovalState = RAD.createState(1, 3);
RAD.clearCommands(groupRemovalState);
groupRemovalState.grid.backlash = 0;
groupRemovalState.grid.couplingGain = 0.8;
const groupRemoval = RAD.compareGroupActuationUnderRemoval(
  groupRemovalState,
  [{ r: 0, c: 0 }, { r: 0, c: 1 }, { r: 0, c: 2 }],
  [{ r: 0, c: 1 }],
  -0.2,
  0.15
);
assert.strictEqual(groupRemoval.intact.decomposes, true);
assert.strictEqual(groupRemoval.removed.decomposes, true);
assert.strictEqual(groupRemoval.lostAppliedCells.map((cell) => `${cell.r},${cell.c}`).join("|"), "0,1");
assert.strictEqual(groupRemoval.newlySkippedCells.map((cell) => `${cell.r},${cell.c}`).join("|"), "0,1");
assert.strictEqual(groupRemoval.componentCountDelta, 1);
assert.strictEqual(groupRemoval.deletedEdgeDelta, 2);
assert.ok(groupRemoval.commandError > 0, "removed group cell should create command support loss");
const verticalRemovalState = RAD.createState(1, 5);
RAD.clearCommands(verticalRemovalState);
verticalRemovalState.grid.pinRadius = 0.10;
verticalRemovalState.grid.holeRadius = 0.12;
verticalRemovalState.grid.zCouplingGain = 0.5;
const verticalRemoval = RAD.compareVerticalResidualUnderRemoval(
  verticalRemovalState,
  [{ r: 0, c: 0 }, { r: 0, c: 2 }],
  [{ r: 0, c: 2 }],
  0.4,
  0,
  1e-9,
  1,
  { fixedCells: [{ r: 0, c: 0 }], externalZLoad: -0.2 }
);
assert.strictEqual(verticalRemoval.fixedCells.map((cell) => `${cell.r},${cell.c}`).join("|"), "0,0");
assert.strictEqual(Number(verticalRemoval.externalZLoad.toFixed(6)), -0.2);
assert.strictEqual(verticalRemoval.removedSourceCells.map((cell) => `${cell.r},${cell.c}`).join("|"), "0,2");
assert.strictEqual(verticalRemoval.affectedNeighborCells.map((cell) => `${cell.r},${cell.c}`).join("|"), "0,1");
assert.ok(
  verticalRemoval.topologyBlockedCells.some((cell) => cell.r === 0 && cell.c === 3),
  "removed cell should block right-side z residual propagation"
);
assert.ok(
  verticalRemoval.topologyBlockedCells.some((cell) => cell.r === 0 && cell.c === 4),
  "removed cell should block farther right-side z residual propagation"
);
assert.strictEqual(verticalRemoval.componentCountDelta, 1);
assert.strictEqual(verticalRemoval.deletedEdgeDelta, 2);
assert.ok(
  verticalRemoval.intactDieOffRadius > verticalRemoval.removedDieOffRadius,
  "removed topology should reduce vertical residual die-off radius"
);
assert.ok(verticalRemoval.maxAbsResidualDelta > 0, "removed topology should change the residual field");
assert.strictEqual(Number(verticalRemoval.clearance.toFixed(6)), 0.02);
assert.strictEqual(verticalRemoval.intactContactEngagedCells.length, 2);
assert.strictEqual(verticalRemoval.removedContactEngagedCells.length, 1);
assert.ok(
  verticalRemoval.contactPenaltyDelta < 0,
  "removing a z-commanded source should reduce the uncalibrated contact penalty"
);
assert.ok(
  !verticalRemoval.intactLoadActiveCells.some((cell) => cell.r === 0 && cell.c === 0),
  "fixed cells should be excluded from load-active diagnostic cells"
);
assert.ok(
  verticalRemoval.intactLoadActiveCells.some((cell) => cell.r === 0 && cell.c === 1),
  "neighbor residual cell should remain load-active"
);
assert.ok(
  !verticalRemoval.removedLoadActiveCells.some((cell) => cell.r === 0 && cell.c === 2),
  "removed cells should be excluded from load-active diagnostic cells"
);
assert.ok(verticalRemoval.loadWorkMagnitudeDelta < 0, "removing topology should reduce load-work magnitude in this fixture");
assert.ok(
  !verticalRemoval.intactHeightContactEngagedCells.some((cell) => cell.r === 0 && cell.c === 0),
  "fixed cells should be excluded from height-contact engagement"
);
assert.ok(
  verticalRemoval.heightContactPenaltyDelta < 0,
  "removing topology should reduce height-contact penalty in this fixture"
);
const orderState = RAD.createState(3, 3);
RAD.clearCommands(orderState);
const orderDiagnostic = RAD.compareEventOrder(
  orderState,
  RAD.localActuationEvent({ r: 1, c: 1 }, -0.3, 0.2),
  RAD.lockEvent({ r: 1, c: 1 })
);
assert.strictEqual(orderDiagnostic.modeCommutes, true);
assert.strictEqual(orderDiagnostic.commandCommutes, false);
assert.strictEqual(orderDiagnostic.lockAlphaCommutes, false);
assert.strictEqual(orderDiagnostic.lockZCommutes, false);
assert.ok(orderDiagnostic.finalAlphaError > 0.1, "actuation and lock events should not commute in alpha");
assert.ok(orderDiagnostic.finalHeightError > 0.1, "actuation and lock events should not commute in height");
const sequenceDiagnostic = RAD.compareSequenceOrder(orderState, [
  RAD.localActuationEvent({ r: 1, c: 1 }, -0.3, 0.2),
  RAD.lockEvent({ r: 1, c: 1 }),
  RAD.clearActuationEvent({ r: 1, c: 1 }),
]);
assert.strictEqual(sequenceDiagnostic.eventCount, 3);
assert.strictEqual(sequenceDiagnostic.adjacentPairCount, 2);
assert.ok(sequenceDiagnostic.noncommutingAdjacentPairs > 0, "local operator sequence should identify noncommuting adjacent swaps");
assert.strictEqual(sequenceDiagnostic.orderSensitive, true, "local operator sequence should be order sensitive");
assert.ok(sequenceDiagnostic.maxOrderError > 0.1, "sequence order sensitivity should report a nonzero max error");
assert.ok(Number.isFinite(sequenceDiagnostic.reverseAlphaError), "sequence order sensitivity should report finite reverse alpha error");
assert.ok(Number.isFinite(sequenceDiagnostic.reverseHeightError), "sequence order sensitivity should report finite reverse height error");
assert.strictEqual(sequenceDiagnostic.adjacent.length, 2);

const zResidualState = RAD.createState(5, 5);
zResidualState.grid.backlash = 0.02;
zResidualState.grid.zCouplingGain = 0.35;
zResidualState.grid.pinRadius = 0.18;
zResidualState.grid.holeRadius = 0.225;
zResidualState.cells.commandAlpha = RAD.matrix(5, 5, 0);
zResidualState.cells.commandZ[2][2] = 0.4;
const zResidualSim = RAD.simulate(zResidualState);
assert.strictEqual(Number(zResidualSim.height[2][2].toFixed(6)), 0.4);
assert.ok(zResidualSim.height[2][3] > 0, "vertical actuation should leave residual motion in neighboring cells");
assert.ok(Math.abs(zResidualSim.height[2][3]) < Math.abs(zResidualSim.height[2][2]), "neighbor residual should decay");
const zFootprint = RAD.selectedCellFootprint(zResidualState, 2, 2);
assert.strictEqual(Number(zFootprint.zResidual[2][2].toFixed(6)), 0.4);
assert.ok(zFootprint.zResidual[2][3] > 0, "selected footprint should expose vertical spillover neighbors");
assert.ok(Number.isFinite(zFootprint.zDieOff[2][3]), "selected footprint should expose z die-off distance");
const zMetrics = RAD.selectedCellCouplingMetrics(zResidualState, 2, 2);
assert.strictEqual(Number(zMetrics.zNeighborSignal.toFixed(6)), Number(zFootprint.zResidual[2][3].toFixed(6)));
assert.ok(zMetrics.zReachCells > 0, "selected coupling metrics should count z spillover reach");
assert.strictEqual(zMetrics.zCoupled, true, "selected coupling metrics should report z coupled outside dead-zone");
const upZCharacterization = RAD.characterizeLocalResponse(zResidualState, { r: 2, c: 2, scope: "single" });
assert.ok(upZCharacterization.positiveZReachCells > 0, "characterization should count upward z response cells");
assert.strictEqual(upZCharacterization.negativeZReachCells, 0);
assert.ok(upZCharacterization.maxPositiveHeightDelta > 0);
const downZState = RAD.createState(5, 5);
downZState.grid.backlash = 0.02;
downZState.grid.zCouplingGain = 0.35;
downZState.grid.pinRadius = 0.18;
downZState.grid.holeRadius = 0.225;
downZState.cells.commandAlpha = RAD.matrix(5, 5, 0);
downZState.cells.commandZ[2][2] = -0.4;
const downZCharacterization = RAD.characterizeLocalResponse(downZState, { r: 2, c: 2, scope: "single" });
assert.ok(downZCharacterization.negativeZReachCells > 0, "characterization should count downward z response cells");
assert.strictEqual(downZCharacterization.positiveZReachCells, 0);
assert.ok(downZCharacterization.maxNegativeHeightDelta < 0);
const characterizationState = RAD.createState(5, 5);
characterizationState.grid.backlash = 0.02;
characterizationState.grid.zCouplingGain = 0.35;
characterizationState.cells.commandAlpha[2][2] = -0.35;
characterizationState.cells.commandZ[2][2] = 0.4;
characterizationState.cells.commandAlpha[2][3] = -0.18;
characterizationState.cells.commandZ[2][3] = 0.2;
const singleCharacterization = RAD.characterizeLocalResponse(characterizationState, { r: 2, c: 2, scope: "single" });
const pairCharacterization = RAD.characterizeLocalResponse(characterizationState, { r: 2, c: 2, scope: "pair" });
assert.strictEqual(singleCharacterization.scope, "single");
assert.strictEqual(singleCharacterization.regionCellCount, 1);
assert.strictEqual(pairCharacterization.scope, "pair");
assert.strictEqual(pairCharacterization.regionCellCount, 2);
assert.strictEqual(pairCharacterization.activeSources, 2);
assert.ok(pairCharacterization.responseCells >= singleCharacterization.responseCells, "pair characterization should include at least the single-cell response footprint");
assert.ok(Number.isFinite(pairCharacterization.superpositionError), "pair characterization should report finite superposition residual");
assert.strictEqual(pairCharacterization.pairwiseInteractionModel, "pairwise superposition residual", "pair characterization should report pairwise interaction model");
assert.strictEqual(pairCharacterization.pairwiseTotalPairs, 1, "selected pair should have one source pair");
assert.strictEqual(pairCharacterization.pairwiseEvaluatedPairs, 1, "selected pair should evaluate one source pair");
assert.ok(Number.isFinite(pairCharacterization.pairwiseMaxInteractionError), "pair characterization should report finite pairwise max interaction");
assert.ok(pairCharacterization.pairwiseInteractions.length <= 1, "pair characterization should store bounded pairwise rows");
assert.ok(Array.isArray(pairCharacterization.pairwiseInteractionMap), "pair characterization should include a pairwise hotspot map");
assert.strictEqual(pairCharacterization.pairwiseInteractionMap.length, characterizationState.grid.rows, "pairwise hotspot map should match lattice rows");
assert.strictEqual(pairCharacterization.pairwiseInteractionMap[0].length, characterizationState.grid.cols, "pairwise hotspot map should match lattice columns");
assert.ok(Number.isFinite(pairCharacterization.pairwiseInteractionMapMax), "pair characterization should report finite pairwise hotspot scale");
assert.ok(pairCharacterization.pairwiseInteractionMapMax >= 0, "pairwise hotspot scale should be nonnegative");
assert.ok(Array.isArray(pairCharacterization.pairwiseInteractionDegreeMap), "pair characterization should include a nonadditive degree map");
assert.strictEqual(pairCharacterization.pairwiseInteractionDegreeMap.length, characterizationState.grid.rows, "pairwise degree map should match lattice rows");
assert.strictEqual(pairCharacterization.pairwiseInteractionDegreeMap[0].length, characterizationState.grid.cols, "pairwise degree map should match lattice columns");
assert.ok(Number.isFinite(pairCharacterization.pairwiseInteractionDegreeMax), "pair characterization should report finite max degree");
assert.ok(Number.isFinite(pairCharacterization.pairwiseInteractionDensity), "pair characterization should report finite interaction density");
assert.ok(pairCharacterization.pairwiseInteractionDensity >= 0 && pairCharacterization.pairwiseInteractionDensity <= 1, "interaction density should be normalized");
assert.deepStrictEqual(
  JSON.parse(JSON.stringify(pairCharacterization.pairwiseAlphaErrorMatrix)),
  JSON.parse(JSON.stringify(pairCharacterization.pairwiseAlphaErrorMatrix.map((row, r) => row.map((_, c) => pairCharacterization.pairwiseAlphaErrorMatrix[c][r])))),
  "pairwise alpha interaction matrix should be symmetric"
);
assert.ok(pairCharacterization.pinHoleClearanceMm > 0, "characterization should carry paper-scale clearance");
assert.ok(pairCharacterization.diagnosticColumnCount >= 4, "pair characterization should build alpha and z response columns");
assert.ok(pairCharacterization.responseRankAlpha > 0, "pair characterization should report alpha response rank");
assert.ok(pairCharacterization.responseRankHeight > 0, "pair characterization should report height response rank");
assert.ok(pairCharacterization.reachableAlphaCells > 0, "pair characterization should report reachable alpha cells");
assert.ok(pairCharacterization.reachableHeightCells > 0, "pair characterization should report reachable height cells");
assert.ok(pairCharacterization.alphaUnderactuatedCells >= 0, "pair characterization should report alpha underactuation");
assert.ok(pairCharacterization.heightUnderactuatedCells >= 0, "pair characterization should report height underactuation");
const responseMatrix = RAD.buildResponseMatrix(characterizationState, { r: 2, c: 2, scope: "pair" });
assert.strictEqual(responseMatrix.schema, "rad-sim.response-matrix.v1");
assert.strictEqual(responseMatrix.commands.length, responseMatrix.diagnostics.columnCount);
assert.ok(responseMatrix.commands.length >= 4, "pair response matrix should include alpha and z columns");
assert.strictEqual(responseMatrix.alpha.length, characterizationState.grid.rows * characterizationState.grid.cols);
assert.strictEqual(responseMatrix.height.length, characterizationState.grid.rows * characterizationState.grid.cols);
assert.strictEqual(responseMatrix.alpha[0].length, responseMatrix.commands.length);
assert.ok(responseMatrix.diagnostics.alphaRank > 0, "response matrix should report alpha rank");
assert.ok(responseMatrix.diagnostics.reachableHeightCells > 0, "response matrix should report height reachability");
assert.strictEqual(JSON.parse(RAD.exportResponseMatrix(characterizationState, { r: 2, c: 2, scope: "pair" })).schema, responseMatrix.schema);
const frameworkReport = RAD.programmableDiscontinuityReport(characterizationState, { r: 2, c: 2, scope: "pair", includeResponseMatrix: false, includeFields: false });
assert.strictEqual(frameworkReport.schema, "rad-sim.programmable-discontinuity-report.v1");
assert.strictEqual(frameworkReport.operators.activeOperatorCount, 2, "framework report should count active command operators");
assert.strictEqual(frameworkReport.paperSupportedAssumptions[0].formula, "f(x)=max(0,x-b)+min(x+b,0)");
assert.strictEqual(frameworkReport.simulatorDiagnostics[1].name, "superposition residual");
assert.strictEqual(frameworkReport.operatorLawCandidates.schema, "rad-sim.framework-law-candidates.v1");
assert.strictEqual(
  frameworkReport.operatorLawCandidates.method,
  "thresholded diagnostic predicates over locality, reachability, composition, and event-order metrics"
);
const frameworkLaws = new Map(frameworkReport.operatorLawCandidates.laws.map((law) => [law.id, law]));
assert.ok(frameworkLaws.has("composition_nonadditivity"), "framework report should include a composition law candidate");
assert.ok(frameworkLaws.has("event_order_noncommutativity"), "framework report should include an order law candidate");
assert.strictEqual(frameworkLaws.get("event_order_noncommutativity").status, "simulator-diagnostic");
assert.strictEqual(frameworkReport.formalizationTargets.schema, "rad-sim.formalization-targets.v1");
assert.strictEqual(frameworkReport.formalizationTargets.tooling.engine, "Lean");
assert.strictEqual(frameworkReport.formalizationTargets.tooling.status, "browser-cannot-inspect-path");
assert.deepStrictEqual(RAD.formalizationTargetManifest(characterizationState, { r: 2, c: 2, scope: "pair" }), frameworkReport.formalizationTargets);
assert.strictEqual(JSON.parse(RAD.exportFormalizationTargetManifest(characterizationState, { r: 2, c: 2, scope: "pair" })).schema, "rad-sim.formalization-targets.v1");
const formalTargets = new Map(frameworkReport.formalizationTargets.targets.map((target) => [target.id, target]));
assert.strictEqual(formalTargets.get("dead_zone_zero_inside_backlash").source, "paper-supported");
assert.strictEqual(formalTargets.get("group_operator_support_decomposition").status, "lean-proved-discrete");
assert.ok(
  formalTargets.get("group_operator_support_decomposition").evidence.leanTheorems.includes("LocalModeOperator.groupOperators_with_disjoint_lists_commute"),
  "framework report should expose the Lean-proven group support theorem"
);
assert.strictEqual(formalTargets.get("removed_cell_clears_group_supported_constraints").status, "lean-proved-discrete");
assert.ok(
  formalTargets.get("removed_cell_clears_group_supported_constraints").evidence.leanTheorems.includes("CellGraph.removed_cell_constraints_clear_group_operator"),
  "framework report should expose the Lean-proven removal/group theorem"
);
assert.strictEqual(formalTargets.get("vertical_clearance_gates_residual_contact").status, "lean-proved-discrete");
assert.ok(
  formalTargets.get("vertical_clearance_gates_residual_contact").evidence.leanTheorems.includes("verticalResidualStepNat_zero_inside"),
  "framework report should expose the Lean-proven clearance residual theorem"
);
assert.ok(
  formalTargets.get("vertical_clearance_gates_residual_contact").evidence.leanTheorems.includes("contactPenaltyFromClearanceNat_nonnegative"),
  "framework report should expose the Lean-proven clearance contact theorem"
);
assert.strictEqual(formalTargets.get("fixed_cell_load_work_proxy_zero").status, "lean-proved-discrete");
assert.ok(
  formalTargets.get("fixed_cell_load_work_proxy_zero").evidence.leanTheorems.includes("loadWorkMagnitudeNat_zero_fixed"),
  "framework report should expose the Lean-proven fixed-cell load-work theorem"
);
assert.strictEqual(formalTargets.get("signed_vertical_load_work_residual_int").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("signed_vertical_load_work_residual_int").evidence.cliModule,
  "rad_sim.compare_vertical_load_measurements"
);
assert.ok(
  formalTargets.get("signed_vertical_load_work_residual_int").evidence.leanTheorems.includes("signedLoadWorkInt_zero_displacement"),
  "framework report should expose the Lean-proven signed load-work zero-displacement theorem"
);
assert.ok(
  formalTargets.get("signed_vertical_load_work_residual_int").evidence.leanTheorems.includes("signedEnergyResidualTripleInt_zero_when_equal"),
  "framework report should expose the Lean-proven signed residual triple theorem"
);
assert.strictEqual(formalTargets.get("integer_scaled_mechanics_scaffold").status, "lean-proved-discrete");
assert.ok(
  formalTargets.get("integer_scaled_mechanics_scaffold").evidence.leanTheorems.includes("scaledSpringEnergyNat_preserves_denominator"),
  "framework report should expose the Lean-proven scaled spring denominator theorem"
);
assert.ok(
  formalTargets.get("integer_scaled_mechanics_scaffold").evidence.leanTheorems.includes("scaledSignedEnergyResidualTripleInt_zero_when_equal"),
  "framework report should expose the Lean-proven scaled signed residual theorem"
);
assert.ok(
  formalTargets.get("integer_scaled_mechanics_scaffold").evidence.nextFormalStep.includes("Rat/Real"),
  "framework report should identify the rational/real next step"
);
assert.strictEqual(formalTargets.get("measurement_unit_scale_invariants").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("measurement_unit_scale_invariants").evidence.pythonFunction,
  "calibrate_paper_rad_config"
);
assert.strictEqual(
  formalTargets.get("measurement_unit_scale_invariants").evidence.unitScaleSchema,
  "rad-sim.physical-unit-scale-metadata.v1"
);
assert.strictEqual(
  formalTargets.get("measurement_unit_scale_invariants").evidence.unitScaleFunction,
  "physical_unit_scale_metadata"
);
assert.strictEqual(
  formalTargets.get("measurement_unit_scale_invariants").evidence.hardwareProfileSchema,
  "rad-sim.hardware-profile.v1"
);
assert.ok(
  formalTargets
    .get("measurement_unit_scale_invariants")
    .evidence.hardwareProfileFunctions.includes("hardware_profile_from_json"),
  "framework report should expose measured hardware profile JSON loading"
);
assert.ok(
  formalTargets.get("measurement_unit_scale_invariants").evidence.leanTheorems.includes("measurementUnitScaleNat_preserves_denominator"),
  "framework report should expose the Lean-proven unit-scale denominator theorem"
);
assert.ok(
  formalTargets.get("measurement_unit_scale_invariants").evidence.leanTheorems.includes("measurementUnitScaleResidualInt_zero_when_equal"),
  "framework report should expose the Lean-proven signed unit-scale residual theorem"
);
assert.ok(
  formalTargets.get("measurement_unit_scale_invariants").evidence.leanTheorems.includes("hardwareProfileCoverageMissing_zero_when_complete"),
  "framework report should expose the Lean-proven hardware-profile coverage theorem"
);
assert.strictEqual(formalTargets.get("calibration_parameter_estimate_residual_bookkeeping").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("calibration_parameter_estimate_residual_bookkeeping").evidence.parameterEstimateSchema,
  "rad-sim.calibration-parameter-estimates.v1"
);
assert.ok(
  formalTargets
    .get("calibration_parameter_estimate_residual_bookkeeping")
    .evidence.leanTheorems.includes("calibrationFitResidualPassNat_zero"),
  "framework report should expose calibration fit residual theorem"
);
assert.strictEqual(formalTargets.get("calibration_model_profile_safe_update_bounds").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("calibration_model_profile_safe_update_bounds").evidence.profileSchema,
  "rad-sim.calibration-model-profile.v1"
);
assert.strictEqual(
  formalTargets.get("calibration_model_profile_safe_update_bounds").evidence.applicationSchema,
  "rad-sim.calibration-model-profile-application.v1"
);
assert.ok(
  formalTargets
    .get("calibration_model_profile_safe_update_bounds")
    .evidence.leanTheorems.includes("calibrationModelProfileUpdateSafeNat_intro"),
  "framework report should expose calibration model-profile safety theorem"
);
assert.strictEqual(formalTargets.get("calibration_bench_protocol_coverage").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("calibration_bench_protocol_coverage").evidence.schema,
  "rad-sim.calibration-bench-notebook.v1"
);
assert.strictEqual(
  formalTargets.get("calibration_bench_protocol_coverage").evidence.csvExportFunction,
  "export_calibration_bench_notebook_csv"
);
assert.ok(
  formalTargets
    .get("calibration_bench_protocol_coverage")
    .evidence.leanTheorems.includes("calibrationBenchProtocolCoverageReadyNat_has_two_dataset_roles"),
  "framework report should expose calibration bench protocol coverage theorem"
);
assert.strictEqual(formalTargets.get("calibration_bench_packet_completeness").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("calibration_bench_packet_completeness").evidence.schema,
  "rad-sim.calibration-bench-packet.v1"
);
assert.strictEqual(
  formalTargets.get("calibration_bench_packet_completeness").evidence.writerFunction,
  "write_calibration_bench_packet_artifacts"
);
assert.ok(
  formalTargets
    .get("calibration_bench_packet_completeness")
    .evidence.leanTheorems.includes("calibrationBenchPacketCompleteNat_has_manifest"),
  "framework report should expose calibration bench packet completeness theorem"
);
assert.strictEqual(formalTargets.get("calibration_bench_executed_validation_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("calibration_bench_executed_validation_gate").evidence.schema,
  "rad-sim.calibration-bench-execution-validation.v1"
);
assert.strictEqual(
  formalTargets.get("calibration_bench_executed_validation_gate").evidence.browserFunction,
  "calibrationBenchExecutionValidation"
);
assert.ok(
  formalTargets
    .get("calibration_bench_executed_validation_gate")
    .evidence.leanTheorems.includes("calibrationBenchExecutedValidationReadyNat_zero_missing_evidence"),
  "framework report should expose calibration bench execution validation theorem"
);
assert.strictEqual(formalTargets.get("mechanics_energy_certificate_nonnegative_proxy").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("mechanics_energy_certificate_nonnegative_proxy").evidence.schema,
  "rad-sim.mechanics-energy-certificate.v1"
);
assert.ok(
  formalTargets.get("mechanics_energy_certificate_nonnegative_proxy").evidence.leanTheorems.includes("mechanicalStoredEnergyNat_nonnegative"),
  "framework report should expose the Lean-proven stored-energy certificate theorem"
);
assert.strictEqual(formalTargets.get("measured_vertical_load_work_validation_zero_residual").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("measured_vertical_load_work_validation_zero_residual").evidence.schema,
  "rad-sim.vertical-load-energy-validation.v1"
);
assert.strictEqual(
  formalTargets.get("measured_vertical_load_work_validation_zero_residual").evidence.resultsSchema,
  "rad-sim.vertical-load-energy-measurement-results.v1"
);
assert.strictEqual(
  formalTargets.get("measured_vertical_load_work_validation_zero_residual").evidence.protocolSchema,
  "rad-sim.vertical-load-energy-experiment-protocol.v1"
);
assert.strictEqual(
  formalTargets.get("measured_vertical_load_work_validation_zero_residual").evidence.packetSchema,
  "rad-sim.vertical-load-bench-packet.v1"
);
assert.strictEqual(
  formalTargets.get("measured_vertical_load_work_validation_zero_residual").evidence.reportSchema,
  "rad-sim.vertical-load-energy-comparison-report.v1"
);
assert.ok(
  formalTargets.get("measured_vertical_load_work_validation_zero_residual").evidence.leanTheorems.includes("measuredWorkResidualNat_zero_when_equal"),
  "framework report should expose the Lean-proven measured work residual theorem"
);
assert.strictEqual(formalTargets.get("vertical_load_comparison_pass_predicate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("vertical_load_comparison_pass_predicate").evidence.cliModule,
  "rad_sim.compare_vertical_load_measurements"
);
assert.strictEqual(
  formalTargets.get("vertical_load_comparison_pass_predicate").evidence.reportSchema,
  "rad-sim.vertical-load-energy-comparison-report.v1"
);
assert.ok(
  formalTargets.get("vertical_load_comparison_pass_predicate").evidence.leanTheorems.includes("verticalLoadScenarioPassNat_zero_errors"),
  "framework report should expose the Lean-proven zero-error comparison predicate theorem"
);
assert.ok(
  formalTargets.get("vertical_load_comparison_pass_predicate").evidence.leanTheorems.includes("verticalLoadScenarioPassNat_false_when_missing"),
  "framework report should expose the Lean-proven missing-measurement comparison theorem"
);
assert.strictEqual(formalTargets.get("physical_validation_readiness_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("physical_validation_readiness_gate").evidence.schema,
  "rad-sim.physical-validation-readiness.v1"
);
assert.strictEqual(
  formalTargets.get("physical_validation_readiness_gate").evidence.browserFunction,
  "physicalValidationReadinessReport"
);
assert.ok(
  formalTargets.get("physical_validation_readiness_gate").evidence.leanTheorems.includes("physicalValidationReadyNat_zero_missing_evidence"),
  "framework report should expose physical validation readiness theorem"
);
assert.strictEqual(formalTargets.get("contact_state_abstraction_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("contact_state_abstraction_gate").evidence.schema,
  "rad-sim.contact-state-abstraction.v1"
);
assert.strictEqual(
  formalTargets.get("contact_state_abstraction_gate").evidence.browserFunction,
  "contactStateAbstractionReport"
);
assert.ok(
  formalTargets.get("contact_state_abstraction_gate").evidence.leanTheorems.includes("contactStateAbstractionReadyNat_zero_missing_evidence"),
  "framework report should expose contact-state abstraction theorem"
);
assert.strictEqual(formalTargets.get("contact_graph_consistency_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("contact_graph_consistency_gate").evidence.schema,
  "rad-sim.contact-graph-consistency.v1"
);
assert.strictEqual(
  formalTargets.get("contact_graph_consistency_gate").evidence.browserFunction,
  "contactGraphConsistencyReport"
);
assert.strictEqual(
  formalTargets.get("contact_graph_consistency_gate").evidence.pythonFunction,
  "contact_graph_consistency_report"
);
assert.strictEqual(
  formalTargets.get("contact_graph_consistency_gate").evidence.csvExportFunction,
  "export_contact_graph_consistency_csv"
);
assert.ok(
  formalTargets.get("contact_graph_consistency_gate").evidence.leanTheorems.includes("contactGraphConsistentNat_zero_missing_evidence"),
  "framework report should expose contact graph consistency theorem"
);
assert.strictEqual(formalTargets.get("physical_realization_map_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("physical_realization_map_gate").evidence.schema,
  "rad-sim.physical-realization-map.v1"
);
assert.strictEqual(
  formalTargets.get("physical_realization_map_gate").evidence.browserFunction,
  "physicalRealizationMapReport"
);
assert.strictEqual(
  formalTargets.get("physical_realization_map_gate").evidence.pythonFunction,
  "physical_realization_map_report"
);
assert.strictEqual(
  formalTargets.get("physical_realization_map_gate").evidence.csvExportFunction,
  "export_physical_realization_map_csv"
);
assert.ok(
  formalTargets.get("physical_realization_map_gate").evidence.leanTheorems.includes("physicalRealizationMapReadyNat_zero_missing_evidence"),
  "framework report should expose physical realization map theorem"
);
assert.strictEqual(formalTargets.get("external_physics_engine_audit_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("external_physics_engine_audit_gate").evidence.schema,
  "rad-sim.external-physics-engine-audit.v1"
);
assert.strictEqual(
  formalTargets.get("external_physics_engine_audit_gate").evidence.browserFunction,
  "externalPhysicsEngineAuditReport"
);
assert.strictEqual(
  formalTargets.get("external_physics_engine_audit_gate").evidence.pythonFunction,
  "external_physics_engine_audit_report"
);
assert.strictEqual(
  formalTargets.get("external_physics_engine_audit_gate").evidence.csvExportFunction,
  "export_external_physics_engine_audit_csv"
);
assert.ok(
  formalTargets.get("external_physics_engine_audit_gate").evidence.leanTheorems.includes("externalPhysicsEngineAuditReadyNat_zero_missing_evidence"),
  "framework report should expose external physics engine audit theorem"
);
assert.strictEqual(formalTargets.get("mujoco_model_export_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("mujoco_model_export_gate").evidence.schema,
  "rad-sim.mujoco-model-export.v1"
);
assert.strictEqual(
  formalTargets.get("mujoco_model_export_gate").evidence.browserFunction,
  "mujocoModelExportReport"
);
assert.strictEqual(
  formalTargets.get("mujoco_model_export_gate").evidence.pythonFunction,
  "mujoco_model_export_report"
);
assert.strictEqual(
  formalTargets.get("mujoco_model_export_gate").evidence.xmlExportFunction,
  "export_mujoco_model_xml"
);
assert.ok(
  formalTargets.get("mujoco_model_export_gate").evidence.leanTheorems.includes("externalPhysicsModelExportReadyNat_zero_missing_evidence"),
  "framework report should expose MuJoCo model export theorem"
);
assert.strictEqual(formalTargets.get("mujoco_pin_hole_contact_geometry_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("mujoco_pin_hole_contact_geometry_gate").evidence.schema,
  "rad-sim.mujoco-pin-hole-contact-geometry.v1"
);
assert.strictEqual(
  formalTargets.get("mujoco_pin_hole_contact_geometry_gate").evidence.browserFunction,
  "mujocoPinHoleContactGeometryReport"
);
assert.strictEqual(
  formalTargets.get("mujoco_pin_hole_contact_geometry_gate").evidence.pythonFunction,
  "mujoco_pin_hole_contact_geometry_report"
);
assert.strictEqual(
  formalTargets.get("mujoco_pin_hole_contact_geometry_gate").evidence.csvExportFunction,
  "export_mujoco_pin_hole_contact_geometry_csv"
);
assert.ok(
  formalTargets.get("mujoco_pin_hole_contact_geometry_gate").evidence.leanTheorems.includes("externalContactGeometryReadyNat_zero_missing_evidence"),
  "framework report should expose MuJoCo pin-hole contact geometry theorem"
);
assert.strictEqual(formalTargets.get("mujoco_contact_parameter_profile_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("mujoco_contact_parameter_profile_gate").evidence.schema,
  "rad-sim.mujoco-contact-parameter-profile.v1"
);
assert.strictEqual(
  formalTargets.get("mujoco_contact_parameter_profile_gate").evidence.browserFunction,
  "mujocoContactParameterReport"
);
assert.strictEqual(
  formalTargets.get("mujoco_contact_parameter_profile_gate").evidence.pythonFunction,
  "mujoco_contact_parameter_report"
);
assert.strictEqual(
  formalTargets.get("mujoco_contact_parameter_profile_gate").evidence.csvExportFunction,
  "export_mujoco_contact_parameter_csv"
);
assert.ok(
  formalTargets.get("mujoco_contact_parameter_profile_gate").evidence.leanTheorems.includes("externalContactParameterProfileReadyNat_zero_missing_evidence"),
  "framework report should expose MuJoCo contact parameter theorem"
);
assert.strictEqual(formalTargets.get("contact_parameter_calibration_packet_completeness").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("contact_parameter_calibration_packet_completeness").evidence.schema,
  "rad-sim.contact-parameter-calibration-packet.v1"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_calibration_packet_completeness").evidence.browserFunction,
  "contactParameterCalibrationPacket"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_calibration_packet_completeness").evidence.pythonFunction,
  "contact_parameter_calibration_packet"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_calibration_packet_completeness").evidence.csvExportFunction,
  "export_contact_parameter_calibration_packet_csv"
);
assert.ok(
  formalTargets
    .get("contact_parameter_calibration_packet_completeness")
    .evidence.leanTheorems.includes("contactParameterCalibrationPacketCompleteNat_zero_missing_evidence"),
  "framework report should expose contact parameter calibration packet theorem"
);
assert.strictEqual(formalTargets.get("contact_parameter_bench_validation_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("contact_parameter_bench_validation_gate").evidence.schema,
  "rad-sim.contact-parameter-bench-validation.v1"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_bench_validation_gate").evidence.browserFunction,
  "compareContactParameterCalibrationResults"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_bench_validation_gate").evidence.pythonFunction,
  "compare_contact_parameter_calibration_results"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_bench_validation_gate").evidence.csvExportFunction,
  "export_contact_parameter_bench_validation_csv"
);
assert.ok(
  formalTargets
    .get("contact_parameter_bench_validation_gate")
    .evidence.leanTheorems.includes("contactParameterBenchValidationReadyNat_zero_missing_evidence"),
  "framework report should expose contact parameter bench validation theorem"
);
assert.strictEqual(formalTargets.get("contact_parameter_interval_calibration_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("contact_parameter_interval_calibration_gate").evidence.schema,
  "rad-sim.contact-parameter-interval-calibration.v1"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_interval_calibration_gate").evidence.browserFunction,
  "contactParameterIntervalCalibrationReport"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_interval_calibration_gate").evidence.pythonFunction,
  "contact_parameter_interval_calibration_report"
);
assert.strictEqual(
  formalTargets.get("contact_parameter_interval_calibration_gate").evidence.csvExportFunction,
  "export_contact_parameter_interval_calibration_csv"
);
assert.ok(
  formalTargets
    .get("contact_parameter_interval_calibration_gate")
    .evidence.leanTheorems.includes("contactParameterIntervalCalibrationReadyNat_zero_missing_evidence"),
  "framework report should expose contact parameter interval calibration theorem"
);
assert.strictEqual(formalTargets.get("mujoco_external_run_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("mujoco_external_run_gate").evidence.schema,
  "rad-sim.mujoco-external-run.v1"
);
assert.strictEqual(
  formalTargets.get("mujoco_external_run_gate").evidence.browserFunction,
  "mujocoExternalRunReport"
);
assert.strictEqual(
  formalTargets.get("mujoco_external_run_gate").evidence.pythonFunction,
  "mujoco_external_run_report"
);
assert.ok(
  formalTargets.get("mujoco_external_run_gate").evidence.leanTheorems.includes("externalPhysicsRunReadyNat_zero_missing_evidence"),
  "framework report should expose MuJoCo external run theorem"
);
assert.strictEqual(formalTargets.get("mujoco_external_comparison_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("mujoco_external_comparison_gate").evidence.schema,
  "rad-sim.mujoco-external-comparison.v1"
);
assert.strictEqual(
  formalTargets.get("mujoco_external_comparison_gate").evidence.browserFunction,
  "mujocoExternalComparisonReport"
);
assert.strictEqual(
  formalTargets.get("mujoco_external_comparison_gate").evidence.pythonFunction,
  "mujoco_external_comparison_report"
);
assert.strictEqual(
  formalTargets.get("mujoco_external_comparison_gate").evidence.csvExportFunction,
  "export_mujoco_external_comparison_csv"
);
assert.ok(
  formalTargets.get("mujoco_external_comparison_gate").evidence.leanTheorems.includes("externalPhysicsComparisonReadyNat_zero_missing_evidence"),
  "framework report should expose MuJoCo external comparison theorem"
);
assert.strictEqual(formalTargets.get("equilibrium_relation_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("equilibrium_relation_gate").evidence.schema,
  "rad-sim.equilibrium-relation.v1"
);
assert.strictEqual(
  formalTargets.get("equilibrium_relation_gate").evidence.browserFunction,
  "equilibriumRelationReport"
);
assert.strictEqual(
  formalTargets.get("equilibrium_relation_gate").evidence.pythonFunction,
  "equilibrium_relation_report"
);
assert.strictEqual(
  formalTargets.get("equilibrium_relation_gate").evidence.csvExportFunction,
  "export_equilibrium_relation_csv"
);
assert.ok(
  formalTargets.get("equilibrium_relation_gate").evidence.leanTheorems.includes("equilibriumRelationReadyNat_zero_missing_evidence"),
  "framework report should expose equilibrium relation theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_controllability_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_controllability_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-controllability.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_controllability_gate").evidence.browserFunction,
  "reachableEquilibriumControllabilityReport"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_controllability_gate").evidence.pythonFunction,
  "reachable_equilibrium_controllability_report"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_controllability_gate").evidence.csvExportFunction,
  "export_reachable_equilibrium_controllability_csv"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_controllability_gate").evidence.leanTheorems.includes("reachableEquilibriumControllabilityReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_bench_protocol_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_protocol_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-bench-protocol.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_protocol_gate").evidence.browserFunction,
  "reachableEquilibriumBenchProtocol"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_protocol_gate").evidence.pythonFunction,
  "reachable_equilibrium_bench_protocol"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_protocol_gate").evidence.csvExportFunction,
  "export_reachable_equilibrium_bench_protocol_csv"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_bench_protocol_gate").evidence.leanTheorems.includes("reachableEquilibriumBenchProtocolReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium bench theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_bench_validation_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_validation_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-bench-comparison.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_validation_gate").evidence.browserFunction,
  "compareReachableEquilibriumBenchResults"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_validation_gate").evidence.pythonFunction,
  "compare_reachable_equilibrium_bench_results"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_validation_gate").evidence.csvExportFunction,
  "export_reachable_equilibrium_bench_comparison_csv"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_bench_validation_gate").evidence.resultsSchema,
  "rad-sim.reachable-equilibrium-bench-results.v1"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_bench_validation_gate").evidence.leanTheorems.includes("reachableEquilibriumBenchValidationReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium bench validation theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_amplitude_calibration_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_amplitude_calibration_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-amplitude-calibration.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_amplitude_calibration_gate").evidence.browserFunction,
  "reachableEquilibriumAmplitudeCalibrationReport"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_amplitude_calibration_gate").evidence.pythonFunction,
  "reachable_equilibrium_amplitude_calibration_report"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_amplitude_calibration_gate").evidence.csvExportFunction,
  "export_reachable_equilibrium_amplitude_calibration_csv"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_amplitude_calibration_gate").evidence.leanTheorems.includes("reachableEquilibriumAmplitudeCalibrationReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium amplitude theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_empirical_profile_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_empirical_profile_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-empirical-profile.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_empirical_profile_gate").evidence.browserFunction,
  "reachableEquilibriumEmpiricalProfileFromAmplitude"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_empirical_profile_gate").evidence.pythonFunction,
  "reachable_equilibrium_empirical_profile_from_amplitude"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_empirical_profile_gate").evidence.csvExportFunction,
  "export_reachable_equilibrium_empirical_profile_csv"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_empirical_profile_gate").evidence.leanTheorems.includes("reachableEquilibriumEmpiricalProfileReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium empirical profile theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_profile_inverse_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-profile-inverse.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_gate").evidence.browserFunction,
  "reachableEquilibriumProfileInverseReport"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_profile_inverse_gate").evidence.leanTheorems.includes("reachableEquilibriumProfileInverseReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium profile inverse theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_profile_inverse_acceptance_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_acceptance_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_acceptance_gate").evidence.browserFunction,
  "reachableEquilibriumProfileInverseAcceptanceReport"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_profile_inverse_acceptance_gate").evidence.leanTheorems.includes("reachableEquilibriumProfileInverseAcceptanceReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium profile inverse acceptance theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_profile_inverse_preview_packet_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_packet_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_packet_gate").evidence.browserFunction,
  "reachableEquilibriumProfileInversePreviewPacket"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_packet_gate").evidence.leanTheorems.includes("reachableEquilibriumProfileInversePreviewPacketReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium profile inverse preview-packet theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_profile_inverse_preview_replay_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_replay_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_replay_gate").evidence.browserFunction,
  "reachableEquilibriumProfileInversePreviewReplayReport"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_replay_gate").evidence.leanTheorems.includes("reachableEquilibriumProfileInversePreviewReplayReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium profile inverse preview-replay theorem"
);
assert.strictEqual(formalTargets.get("reachable_equilibrium_profile_inverse_preview_physical_gate").status, "lean-proved-discrete");
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_physical_gate").evidence.schema,
  "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1"
);
assert.strictEqual(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_physical_gate").evidence.browserFunction,
  "reachableEquilibriumProfileInversePreviewPhysicalReport"
);
assert.ok(
  formalTargets.get("reachable_equilibrium_profile_inverse_preview_physical_gate").evidence.leanTheorems.includes("reachableEquilibriumProfileInversePreviewPhysicalReadyNat_zero_missing_evidence"),
  "framework report should expose reachable equilibrium profile inverse preview-physical theorem"
);
assert.strictEqual(formalTargets.get("spring_hinge_removed_topology_load_comparison").status, "simulator-diagnostic");
assert.strictEqual(
  formalTargets.get("spring_hinge_removed_topology_load_comparison").evidence.pythonFunction,
  "compare_vertical_residual_spring_hinge_3d"
);
assert.strictEqual(
  formalTargets.get("spring_hinge_removed_topology_load_comparison").evidence.reportSchema,
  "rad-sim.vertical-load-physical-preview-report.v1"
);
assert.strictEqual(
  formalTargets.get("spring_hinge_removed_topology_load_comparison").evidence.csvExportFunction,
  "export_vertical_load_physical_preview_report_csv"
);
assert.strictEqual(formalTargets.get("browser_spring_preview_removed_edge_deletion").status, "simulator-diagnostic");
assert.strictEqual(
  formalTargets.get("browser_spring_preview_removed_edge_deletion").evidence.browserFunction,
  "RAD.simulatePhysicalRelaxation"
);
assert.ok(
  formalTargets.get("browser_spring_preview_removed_edge_deletion").evidence.browserMetrics.includes("physicalSkippedSpringEdges"),
  "framework report should expose removed-edge browser preview metrics"
);
assert.ok(
  formalTargets.get("browser_spring_preview_removed_edge_deletion").evidence.leanTheorems.includes("CellGraph.removal_deletes_one_step_to_removed"),
  "browser preview diagnostic should reference graph-deletion Lean evidence"
);
assert.ok(
  formalTargets.get("noncommutativity_witness_from_order_error").evidence.orderSensitive,
  "framework report should expose a noncommutativity formalization witness"
);
assert.ok(Number.isFinite(frameworkReport.locality.alphaDecayRatio), "framework report should include locality decay");
assert.ok(frameworkReport.reachability.reachableHeightCells > 0, "framework report should include reachable height cells");
assert.ok(Number.isFinite(frameworkReport.composition.maxPairwiseInteractionError), "framework report should include pairwise composition metrics");
assert.strictEqual(frameworkReport.physicalValidation.schema, "rad-sim.browser-physical-preview.v1", "framework report should include browser physical preview validation");
assert.strictEqual(frameworkReport.physicalValidation.physicalSuccess, true, "framework report physical preview should complete");
assert.ok(Number.isFinite(frameworkReport.physicalValidation.heightRmsModelError), "framework report should include physical height model error");
assert.ok("positiveZReachCells" in frameworkReport.combinedResponse, "framework report should include signed z reach");
assert.ok("maxNegativeHeightDelta" in frameworkReport.combinedResponse, "framework report should include signed z extrema");
assert.ok(!("fields" in frameworkReport.combinedResponse), "framework report should support compact exports without response fields");
assert.ok(!("alpha" in frameworkReport.responseMatrix), "framework report should support compact response-matrix diagnostics");
assert.strictEqual(JSON.parse(RAD.exportProgrammableDiscontinuityReport(characterizationState, { r: 2, c: 2, scope: "pair", includeResponseMatrix: false, includeFields: false })).schema, frameworkReport.schema);
const springHingeDescriptor = RAD.springHingePhysicalPreviewReportDescriptor();
assert.strictEqual(springHingeDescriptor.schema, "rad-sim.vertical-load-physical-preview-report.v1");
assert.strictEqual(springHingeDescriptor.comparisonSchema, "rad-sim.vertical-removal-physical-comparison.v1");
assert.strictEqual(springHingeDescriptor.mechanicsCertificateSchema, "rad-sim.mechanics-energy-certificate.v1");
assert.strictEqual(springHingeDescriptor.energyValidationSchema, "rad-sim.vertical-load-energy-validation.v1");
assert.strictEqual(springHingeDescriptor.energyMeasurementSchema, "rad-sim.vertical-load-energy-measurement-results.v1");
assert.strictEqual(springHingeDescriptor.energyProtocolSchema, "rad-sim.vertical-load-energy-experiment-protocol.v1");
assert.strictEqual(springHingeDescriptor.benchPacketSchema, "rad-sim.vertical-load-bench-packet.v1");
assert.strictEqual(springHingeDescriptor.energyComparisonSchema, "rad-sim.vertical-load-energy-comparison-report.v1");
assert.strictEqual(springHingeDescriptor.browserRunsSolver, false);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("build_vertical_load_physical_preview_report"),
  "browser descriptor should point to the Python spring-hinge report builder"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("export_vertical_load_physical_preview_report_csv"),
  "browser descriptor should point to the Python spring-hinge CSV exporter"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("mechanics_energy_certificate_to_dict"),
  "browser descriptor should point to the Python mechanics certificate exporter"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("validate_vertical_load_energy_measurements"),
  "browser descriptor should point to the Python measured energy validation helper"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("vertical_load_energy_measurement_template"),
  "browser descriptor should point to the Python measured energy template helper"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("vertical_load_energy_experiment_protocol"),
  "browser descriptor should point to the Python measured energy protocol helper"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("export_vertical_load_energy_experiment_protocol_json"),
  "browser descriptor should point to the Python measured energy protocol exporter"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("vertical_load_bench_packet"),
  "browser descriptor should point to the Python bench packet helper"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("export_vertical_load_bench_packet_json"),
  "browser descriptor should point to the Python bench packet exporter"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("vertical_load_energy_measurement_results_from_json"),
  "browser descriptor should point to the Python measured energy result parser"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("compare_vertical_load_energy_measurement_results"),
  "browser descriptor should point to the Python measured energy comparison helper"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("export_vertical_load_energy_measurement_template_json"),
  "browser descriptor should point to the Python measured energy template exporter"
);
assert.ok(
  springHingeDescriptor.pythonFunctions.includes("export_vertical_load_energy_validation_json"),
  "browser descriptor should point to the Python measured energy validation exporter"
);
assert.strictEqual(pairCharacterization.decayModel, "log-linear shell max", "characterization should report the decay fitting model");
assert.ok(Number.isFinite(pairCharacterization.alphaDecayRatio), "pair characterization should report finite alpha decay ratio");
assert.ok(Number.isFinite(pairCharacterization.zDecayRatio), "pair characterization should report finite z decay ratio");
assert.ok(Number.isFinite(pairCharacterization.alphaDecayLength), "pair characterization should report finite alpha decay length");
assert.ok(Number.isFinite(pairCharacterization.zDecayLength), "pair characterization should report finite z decay length");
assert.ok(pairCharacterization.alphaDecayShells > 0, "pair characterization should count alpha decay shells");
assert.ok(pairCharacterization.zDecayShells > 0, "pair characterization should count z decay shells");
assert.ok(pairCharacterization.alphaDecayReach >= 0, "pair characterization should report alpha decay reach");
assert.ok(pairCharacterization.zDecayReach >= 0, "pair characterization should report z decay reach");
assert.strictEqual(pairCharacterization.physicalPreviewAvailable, true, "characterization should compare against the browser physical preview");
assert.strictEqual(pairCharacterization.physicalPreviewSuccess, true, "physical preview comparison should complete");
assert.ok(Number.isFinite(pairCharacterization.physicalHeightRmsError), "physical preview comparison should report finite height RMS error");
assert.ok(Number.isFinite(pairCharacterization.physicalHeightMaxError), "physical preview comparison should report finite height max error");
assert.ok(Number.isFinite(pairCharacterization.physicalCenterRmsError), "physical preview comparison should report finite center RMS error");
assert.ok(Number.isFinite(pairCharacterization.physicalCenterMaxError), "physical preview comparison should report finite center max error");
assert.ok(pairCharacterization.physicalPreviewIterations > 0, "physical preview comparison should report relaxation iterations");
characterizationState.experiment.characterizationScope = "pair";
characterizationState.experiment.characterization = pairCharacterization;
const characterizationRoundtrip = RAD.deserialize(RAD.serialize(characterizationState));
assert.strictEqual(characterizationRoundtrip.experiment.characterizationScope, "pair");
assert.strictEqual(characterizationRoundtrip.experiment.characterization.scope, "pair");
assert.strictEqual(characterizationRoundtrip.experiment.characterization.activeSources, 2);
assert.strictEqual(characterizationRoundtrip.experiment.characterization.responseRankHeight, pairCharacterization.responseRankHeight);
const tightClearanceState = RAD.createState(5, 5);
tightClearanceState.grid.zCouplingGain = 0.35;
tightClearanceState.grid.pinRadius = 0.18;
tightClearanceState.grid.holeRadius = 0.20;
tightClearanceState.cells.commandZ[2][2] = 0.4;
const looseClearanceState = RAD.createState(5, 5);
looseClearanceState.grid.zCouplingGain = 0.35;
looseClearanceState.grid.pinRadius = 0.18;
looseClearanceState.grid.holeRadius = 0.32;
looseClearanceState.cells.commandZ[2][2] = 0.4;
assert.ok(
  RAD.simulate(tightClearanceState).zResidual[2][3] > RAD.simulate(looseClearanceState).zResidual[2][3],
  "larger pin-hole clearance should increase z dead-zone and reduce neighbor residual"
);
zResidualState.grid.zCouplingGain = 0;
const localZSim = RAD.simulate(zResidualState);
assert.strictEqual(Number(localZSim.height[2][3].toFixed(6)), 0);
const localFootprint = RAD.selectedCellFootprint(zResidualState, 2, 2);
assert.strictEqual(Number(localFootprint.zResidual[2][3].toFixed(6)), 0);
const localMetrics = RAD.selectedCellCouplingMetrics(zResidualState, 2, 2);
assert.strictEqual(localMetrics.zReachCells, 0);
assert.strictEqual(localMetrics.zCoupled, false, "selected coupling metrics should report z free when coupling is zero");

const alphaFootprintState = RAD.createState(5, 5);
alphaFootprintState.grid.backlash = 0.05;
alphaFootprintState.grid.couplingGain = 0.5;
alphaFootprintState.cells.commandAlpha[2][2] = -0.35;
const alphaFootprint = RAD.selectedCellFootprint(alphaFootprintState, 2, 2);
assert.ok(alphaFootprint.alpha[2][3] < 0, "selected footprint should expose alpha coupling neighbors");
assert.ok(Math.abs(alphaFootprint.alpha[2][3]) < Math.abs(alphaFootprint.alpha[2][2]), "alpha footprint should decay from selected cell");
const alphaMetrics = RAD.selectedCellCouplingMetrics(alphaFootprintState, 2, 2);
assert.strictEqual(Number(alphaMetrics.alphaNeighborSignal.toFixed(6)), Number(alphaFootprint.alpha[2][3].toFixed(6)));
assert.ok(alphaMetrics.alphaReachCells > 0, "selected coupling metrics should count alpha footprint reach");
assert.strictEqual(alphaMetrics.alphaCoupled, true, "selected coupling metrics should report alpha coupled outside dead-zone");
alphaFootprintState.cells.commandAlpha[2][2] = -0.02;
const alphaGapMetrics = RAD.selectedCellCouplingMetrics(alphaFootprintState, 2, 2);
assert.strictEqual(alphaGapMetrics.alphaCoupled, false, "selected coupling metrics should report alpha free inside backlash gap");

const lockDataState = RAD.createState(2, 12);
assert.strictEqual(lockDataState.cells.positionLocked[0][0], false, "state should initialize positional locks");
for (let r = 0; r < lockDataState.grid.rows; r += 1) {
  lockDataState.cells.locked[r][0] = true;
  lockDataState.cells.locked[r][1] = true;
  lockDataState.cells.locked[r][2] = true;
  lockDataState.cells.positionLocked[r][0] = true;
  lockDataState.cells.positionLocked[r][11] = true;
}
lockDataState.grid.lockStateIndex = 2;
const empiricalPrediction = RAD.predictEmpiricalLockCoordinates(lockDataState);
assert.ok(empiricalPrediction, "measured lock-data surrogate should predict repeated 12-cell rows");
assert.strictEqual(JSON.stringify(empiricalPrediction.failedLockCells), JSON.stringify([2]));
const lockBaseSim = RAD.simulate(lockDataState);
const lockDataSim = RAD.simulateActive(lockDataState);
assert.strictEqual(lockDataSim.metrics.empiricalLockDatasetApplied, true);
assert.ok(lockDataSim.metrics.empiricalLockNearestFile.includes("AOM_L3_1_2_3_303030_2"));
assert.ok(Math.abs(lockDataSim.centers[0][0].z) > 0 || Math.abs(lockDataSim.centers[0][1].z) > 0);
assert.strictEqual(lockDataSim.metrics.empiricalLockEndpointAnchoredRows, 2, "fixed row endpoints should anchor empirical rows");
for (let r = 0; r < lockDataState.grid.rows; r += 1) {
  for (const c of [0, 11]) {
    const anchored = lockDataSim.centers[r][c];
    const fixed = lockBaseSim.centers[r][c];
    assert.ok(Math.abs(anchored.x - fixed.x) < 1e-9, "position lock should preserve endpoint x");
    assert.ok(Math.abs(anchored.y - fixed.y) < 1e-9, "position lock should preserve endpoint y");
    assert.ok(Math.abs(anchored.z - fixed.z) < 1e-9, "position lock should preserve endpoint z");
  }
}

const brushState = RAD.createState(5, 5);
const brushPairs = (radius) => JSON.stringify(RAD.brushCells(brushState, 2, 2, radius).map((cell) => [cell.r, cell.c]));
assert.strictEqual(brushPairs(0), JSON.stringify([[2, 2]]));
assert.strictEqual(brushPairs(1), JSON.stringify([[1, 2], [2, 1], [2, 2], [2, 3], [3, 2]]));
assert.strictEqual(RAD.brushCells(brushState, 2, 2, 2).length, 13, "radius 2 brush should include a Manhattan diamond");
brushState.cells.removed[2][3] = true;
assert.ok(!RAD.brushCells(brushState, 2, 2, 1).some((cell) => cell.r === 2 && cell.c === 3), "brush should skip removed cells");

const presetSignatures = ["center", "dome", "saddle", "ridge", "wave", "corner", "twist", "ring", "checker"].map((preset) => {
  const presetState = RAD.createState(7, 7);
  RAD.applyPreset(presetState, preset);
  const presetSim = RAD.simulate(presetState);
  return JSON.stringify({
    preset,
    target: presetState.target.type,
    commandAlpha: presetState.cells.commandAlpha,
    commandZ: presetState.cells.commandZ,
    locked: presetState.cells.locked,
    meanAlpha: Number(presetSim.metrics.meanAlpha.toFixed(4)),
    maxHeight: Number(presetSim.metrics.maxAbsHeight.toFixed(4)),
  });
});
assert.ok(new Set(presetSignatures).size > 6, "presets should produce varied simulator states");

const sequenceState = RAD.createState(4, 4);
RAD.applyPreset(sequenceState, "center");
RAD.recordEvent(sequenceState, { type: "keyframe", name: "center pose" });
const sequenceJson = RAD.exportExperimentSequence(sequenceState);
const imported = RAD.createState(2, 2);
RAD.importExperimentSequence(imported, sequenceJson);
assert.strictEqual(imported.grid.rows, 4);
assert.strictEqual(imported.grid.cols, 4);
assert.strictEqual(imported.experiment.eventList.length, 2);
assert.deepStrictEqual(JSON.parse(JSON.stringify(imported.experiment.eventList.map((event) => event.type))), ["preset", "keyframe"]);

async function validateLocalHttpEntry() {
  const contentTypes = {
    ".css": "text/css",
    ".html": "text/html",
    ".js": "application/javascript",
  };
  const server = http.createServer((request, response) => {
    const urlPath = request.url === "/" ? "/index.html" : request.url.split("?")[0];
    const resolved = path.normalize(path.join(web, urlPath));
    if (!resolved.startsWith(web)) {
      response.writeHead(403);
      response.end("Forbidden");
      return;
    }
    if (!fs.existsSync(resolved) || fs.statSync(resolved).isDirectory()) {
      response.writeHead(404);
      response.end("Not found");
      return;
    }
    response.writeHead(200, { "Content-Type": contentTypes[path.extname(resolved)] || "application/octet-stream" });
    response.end(fs.readFileSync(resolved));
  });
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  const { port } = server.address();
  try {
    const base = `http://127.0.0.1:${port}`;
    const page = await fetch(`${base}/index.html`);
    assert.strictEqual(page.status, 200);
    assert.ok((await page.text()).includes("RAD CAD-Like Lattice Simulator"));
    for (const ref of localRefs) {
      const asset = await fetch(`${base}/${ref.replace("./", "")}`);
      assert.strictEqual(asset.status, 200, `asset should load over local HTTP: ${ref}`);
    }
  } finally {
    await new Promise((resolve) => server.close(resolve));
  }
}

validateLocalHttpEntry()
  .then(() => console.log("web module validation passed"))
  .catch((error) => {
    console.error(error);
    process.exitCode = 1;
  });
