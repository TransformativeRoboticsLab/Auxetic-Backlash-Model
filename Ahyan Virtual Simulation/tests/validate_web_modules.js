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
};
context.window.window = context.window;
context.window.console = console;
vm.createContext(context);

for (const filename of ["state.js", "provenance.js", "math.js", "operators.js", "inverse.js", "analysis.js", "physics.js", "mesh_export.js"]) {
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
assert.strictEqual(typeof RAD.calibratedMeshDimensions, "function", "mesh exporter should expose calibrated mesh dimensions");
assert.strictEqual(typeof RAD.physicalPreviewComparison, "function", "analysis module should expose physical preview comparison");
assert.strictEqual(typeof RAD.responseDecayProfile, "function", "analysis module should expose response decay profile");
assert.strictEqual(typeof RAD.calibrationExperimentProtocol, "function", "analysis module should expose calibration experiment protocol");
assert.strictEqual(typeof RAD.exportCalibrationExperimentProtocol, "function", "analysis module should export calibration experiment protocol");
assert.strictEqual(typeof RAD.responseAtlas, "function", "analysis module should expose response atlas");
assert.strictEqual(typeof RAD.exportResponseAtlas, "function", "analysis module should export response atlas");
assert.strictEqual(typeof RAD.responseAtlasSweep, "function", "analysis module should expose response atlas sweeps");
assert.strictEqual(typeof RAD.exportResponseAtlasSweep, "function", "analysis module should export response atlas sweeps");
assert.strictEqual(typeof RAD.buildResponseMatrix, "function", "analysis module should expose response matrix export");
assert.strictEqual(typeof RAD.exportResponseMatrix, "function", "analysis module should serialize response matrices");
assert.strictEqual(typeof RAD.programmableDiscontinuityReport, "function", "analysis module should expose programmable-discontinuity reports");
assert.strictEqual(typeof RAD.exportProgrammableDiscontinuityReport, "function", "analysis module should export programmable-discontinuity reports");
assert.strictEqual(typeof RAD.calibrationExperimentResultsTemplate, "function", "analysis module should expose calibration results template");
assert.strictEqual(typeof RAD.exportCalibrationExperimentResultsTemplate, "function", "analysis module should export calibration results template");
assert.strictEqual(typeof RAD.compareCalibrationExperimentResults, "function", "analysis module should compare calibration results");
assert.strictEqual(typeof RAD.summarizeCalibrationComparison, "function", "analysis module should summarize calibration result comparisons");
assert.strictEqual(typeof RAD.calibrationComparisonReport, "function", "analysis module should expose calibration comparison reports");
assert.strictEqual(typeof RAD.exportCalibrationComparisonReport, "function", "analysis module should export calibration comparison reports");
assert.strictEqual(typeof RAD.validateInversePlanPhysical, "function", "inverse module should expose physical inverse validation");
assert.strictEqual(typeof RAD.inverseDesignReport, "function", "inverse module should expose inverse design reports");
assert.strictEqual(typeof RAD.exportInverseDesignReport, "function", "inverse module should export inverse design reports");
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
for (const ref of localRefs) {
  assert.ok(fs.existsSync(path.join(web, ref.replace("./", ""))), `missing local browser asset: ${ref}`);
}
assert.ok(!/https?:\/\//.test(html), "browser entry should not require external scripts or styles");
assert.ok(html.includes('value="operatorInteraction"'), "operator interaction overlay should be available in the browser UI");
assert.ok(html.includes('value="calibrationError"'), "calibration error overlay should be available in the browser UI");
assert.ok(html.includes('value="calibrationResidual"'), "calibration residual overlay should be available in the browser UI");
assert.ok(html.includes('value="underactuated"'), "underactuated target overlay should be available in the browser UI");
assert.ok(html.includes('id="selectInteractionHotspot"'), "response panel should expose a hotspot selection button");
assert.ok(html.includes('id="saveResponseMatrix"'), "response panel should expose a response-matrix export button");
assert.ok(html.includes('id="saveProgrammableReport"'), "response panel should expose a programmable-discontinuity report export button");
assert.ok(html.includes('id="selectCalibrationHotspot"'), "response panel should expose a calibration-error selection button");
assert.ok(html.includes('id="nextCalibrationHotspot"'), "response panel should expose ranked calibration-error navigation");
assert.ok(html.includes('id="selectUnderactuatedTarget"'), "inverse panel should expose underactuated target selection");
assert.ok(html.includes('id="saveInverseReport"'), "inverse panel should expose inverse report export");
assert.ok(html.includes('id="underTargetCells"'), "metric strip should count underactuated targets");
assert.ok(html.includes('id="inverseReachabilitySummary"'), "inverse plan readout should summarize underactuated targets");
assert.ok(html.includes('id="saveExperimentProtocol"'), "response panel should expose a protocol export button");
assert.ok(html.includes('id="saveResponseAtlas"'), "response panel should expose a response-atlas export button");
assert.ok(html.includes('id="runResponseAtlasSweep"'), "response panel should expose a sweep run button");
assert.ok(html.includes('id="saveResponseAtlasSweep"'), "response panel should expose a response-atlas sweep export button");
assert.ok(html.includes('id="saveResultsTemplate"'), "response panel should expose a results-template export button");
assert.ok(html.includes('id="loadResultsJson"'), "response panel should expose a results import button");
assert.ok(html.includes('id="saveComparisonReport"'), "response panel should expose a comparison-report export button");
assert.ok(html.includes('id="calibrationResultsFileInput"'), "response panel should include hidden calibration results file input");
assert.ok(html.includes('id="calibrationResultsSummary"'), "response panel should expose calibration result summary");
assert.ok(html.includes('id="calibrationResultsError"'), "response panel should expose calibration result error readout");
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
assert.ok(html.includes('id="saveCalibrationPlan"'), "browser UI should expose calibration plan export action");
assert.ok(html.includes('value="calibratedRad"'), "browser UI should expose calibrated RAD visual mode");
const scriptOrder = [
  "./vendor/three.min.js",
  "./state.js",
  "./provenance.js",
  "./math.js",
  "./operators.js",
  "./inverse.js",
  "./analysis.js",
  "./physics.js",
  "./mesh_export.js",
  "./renderer.js",
  "./ui.js",
  "./app.js",
];
for (let i = 1; i < scriptOrder.length; i += 1) {
  assert.ok(html.indexOf(scriptOrder[i - 1]) < html.indexOf(scriptOrder[i]), `${scriptOrder[i - 1]} should load before ${scriptOrder[i]}`);
}
for (const filename of ["renderer.js", "ui.js", "app.js"]) {
  const source = fs.readFileSync(path.join(web, filename), "utf8");
  assert.doesNotThrow(() => new Function(source), `${filename} should compile`);
  if (filename === "ui.js") {
    assert.ok(source.includes("selectCalibrationHotspot(1)"), "UI should cycle through ranked calibration hotspots");
    assert.ok(source.includes("fitResidualTopCells"), "UI should support ranked residual calibration hotspots");
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

const state = RAD.createState(5, 6);
assert.strictEqual(state.cells.alpha.length, 5);
assert.strictEqual(state.cells.theta[0].length, 6);
assert.strictEqual(state.cells.z[2][3], 0);
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
state.view.paintRadius = 2;
state.cells.commandAlpha[2][3] = -0.24;
state.cells.commandZ[2][3] = 0.31;
state.cells.locked[0][0] = true;
state.cells.lockAlpha[0][0] = 0.82;
state.cells.actuatorAllowed[4][5] = false;
state.selection = { r: 2, c: 3 };
state.experiment.presetName = "custom-roundtrip";
state.experiment.notes = "roundtrip validation";

const serialized = RAD.serialize(state);
const parsedSerialized = JSON.parse(serialized);
assert.ok(Array.isArray(parsedSerialized.cells.alpha), "serialized cells should include derived alpha");
assert.ok(Array.isArray(parsedSerialized.cells.theta), "serialized cells should include derived theta");
assert.ok(Array.isArray(parsedSerialized.cells.z), "serialized cells should include derived z");
assert.ok(Array.isArray(parsedSerialized.cells.lockAlpha), "serialized cells should include committed lock alpha");
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
assert.strictEqual(JSON.parse(RAD.exportResponseAtlasSweep(readyState, { backlashValues: [0.02], clearanceValues: [0.02] })).schema, browserSweep.schema);
const resultsTemplate = RAD.calibrationExperimentResultsTemplate(readyState, { protocol: experimentProtocol });
assert.strictEqual(resultsTemplate.schema, "rad-sim.calibration-experiment-results.v1");
assert.strictEqual(resultsTemplate.steps.length, experimentProtocol.steps.reduce((sum, step) => sum + step.repeatCount, 0));
const exportedResultsTemplate = JSON.parse(RAD.exportCalibrationExperimentResultsTemplate(readyState, { protocol: experimentProtocol }));
assert.strictEqual(exportedResultsTemplate.schema, resultsTemplate.schema);
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
assert.deepStrictEqual(JSON.parse(JSON.stringify(restored.selection)), state.selection);
assert.strictEqual(restored.view.paintRadius, 2);
assert.strictEqual(restored.cells.commandAlpha[2][3], -0.24);
assert.strictEqual(restored.cells.commandZ[2][3], 0.31);
assert.strictEqual(restored.cells.locked[0][0], true);
assert.strictEqual(restored.cells.lockAlpha[0][0], 0.82);
assert.strictEqual(restored.cells.actuatorAllowed[4][5], false);
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
const activeSim = RAD.simulateActive(restored);
assert.strictEqual(activeSim.metrics.model, "spring-preview");
assert.strictEqual(activeSim.metrics.physicalPreview, true);
assert.ok(Number.isFinite(activeSim.metrics.physicalRmsHeightDelta));
assert.ok(Number.isFinite(activeSim.metrics.physicalMaxHeightDelta));
assert.ok(Number.isFinite(activeSim.metrics.physicalRmsCenterDelta));
assert.ok(Number.isFinite(activeSim.metrics.physicalMaxCenterDelta));
assert.ok(activeSim.metrics.physicalIterations > 0);
assert.strictEqual(activeSim.modelErrorHeight.length, restored.grid.rows);
assert.strictEqual(activeSim.modelErrorCenter[0].length, restored.grid.cols);
assert.strictEqual(
  Number(activeSim.modelErrorHeight[2][3].toFixed(6)),
  Number((activeSim.height[2][3] - sim.height[2][3]).toFixed(6)),
  "spring preview should expose per-cell height disagreement from the kinematic state"
);
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
assert.strictEqual(underactuatedJacobian.targetReachability.targetHeightCells, 1);
assert.strictEqual(underactuatedJacobian.targetReachability.underactuatedHeightCells, 1);
assert.strictEqual(underactuatedJacobian.targetReachability.worstUnderactuatedCell.row, 2);
assert.strictEqual(underactuatedJacobian.targetReachability.worstUnderactuatedCell.col, 2);
assert.ok(underactuatedJacobian.targetReachability.maxUnreachableHeightResidual > 0, "underactuated target report should expose residual magnitude");
const underactuatedFit = RAD.solveLinearizedTargetFit(underactuatedTargetState, { maxActuators: 1, maxColumns: 2 });
assert.strictEqual(underactuatedFit.underactuatedHeightCells, 1, "linear fit should carry underactuated target count");
assert.strictEqual(underactuatedFit.targetReachability.worstUnderactuatedCell.col, 2);

const operatorState = RAD.createState(3, 3);
RAD.clearCommands(operatorState);
const lockedEventState = RAD.applyEventSequence(operatorState, [
  RAD.localActuationEvent({ r: 1, c: 1 }, -0.3, 0.2),
  RAD.lockEvent({ r: 1, c: 1 }),
  RAD.clearActuationEvent(),
]);
assert.strictEqual(lockedEventState.cells.locked[1][1], true);
assert.strictEqual(Number(lockedEventState.cells.lockAlpha[1][1].toFixed(6)), 0.7);
assert.strictEqual(Number(lockedEventState.cells.commandAlpha[1][1].toFixed(6)), 0);
assert.strictEqual(Number(RAD.simulate(lockedEventState).alpha[1][1].toFixed(6)), 0.7);
assert.strictEqual(RAD.finiteDieOffRadius(RAD.simulate(lockedEventState).dieOff), 0);
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
assert.ok(Number.isFinite(frameworkReport.locality.alphaDecayRatio), "framework report should include locality decay");
assert.ok(frameworkReport.reachability.reachableHeightCells > 0, "framework report should include reachable height cells");
assert.ok(Number.isFinite(frameworkReport.composition.maxPairwiseInteractionError), "framework report should include pairwise composition metrics");
assert.strictEqual(frameworkReport.physicalValidation.schema, "rad-sim.browser-physical-preview.v1", "framework report should include browser physical preview validation");
assert.strictEqual(frameworkReport.physicalValidation.physicalSuccess, true, "framework report physical preview should complete");
assert.ok(Number.isFinite(frameworkReport.physicalValidation.heightRmsModelError), "framework report should include physical height model error");
assert.ok(!("fields" in frameworkReport.combinedResponse), "framework report should support compact exports without response fields");
assert.ok(!("alpha" in frameworkReport.responseMatrix), "framework report should support compact response-matrix diagnostics");
assert.strictEqual(JSON.parse(RAD.exportProgrammableDiscontinuityReport(characterizationState, { r: 2, c: 2, scope: "pair", includeResponseMatrix: false, includeFields: false })).schema, frameworkReport.schema);
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

const brushState = RAD.createState(5, 5);
const brushPairs = (radius) => JSON.stringify(RAD.brushCells(brushState, 2, 2, radius).map((cell) => [cell.r, cell.c]));
assert.strictEqual(brushPairs(0), JSON.stringify([[2, 2]]));
assert.strictEqual(brushPairs(1), JSON.stringify([[1, 2], [2, 1], [2, 2], [2, 3], [3, 2]]));
assert.strictEqual(RAD.brushCells(brushState, 2, 2, 2).length, 13, "radius 2 brush should include a Manhattan diamond");

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
