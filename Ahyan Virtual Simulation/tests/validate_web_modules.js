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

for (const filename of ["state.js", "math.js", "operators.js", "inverse.js", "analysis.js", "physics.js", "mesh_export.js"]) {
  const source = fs.readFileSync(path.join(web, filename), "utf8");
  vm.runInContext(source, context, { filename });
}

const RAD = context.window.RAD;
assert.ok(RAD, "RAD namespace should load");
assert.strictEqual(typeof RAD.physicalPreviewComparison, "function", "analysis module should expose physical preview comparison");
assert.strictEqual(typeof RAD.responseDecayProfile, "function", "analysis module should expose response decay profile");
assert.strictEqual(typeof RAD.validateInversePlanPhysical, "function", "inverse module should expose physical inverse validation");

const html = fs.readFileSync(path.join(web, "index.html"), "utf8");
const localRefs = Array.from(html.matchAll(/(?:src|href)="(\.\/[^"]+)"/g)).map((match) => match[1]);
assert.ok(localRefs.includes("./vendor/three.min.js"), "local Three.js asset should be referenced");
for (const ref of localRefs) {
  assert.ok(fs.existsSync(path.join(web, ref.replace("./", ""))), `missing local browser asset: ${ref}`);
}
assert.ok(!/https?:\/\//.test(html), "browser entry should not require external scripts or styles");
assert.ok(html.includes('value="operatorInteraction"'), "operator interaction overlay should be available in the browser UI");
assert.ok(html.includes('id="selectInteractionHotspot"'), "response panel should expose a hotspot selection button");
assert.ok(html.includes('id="characterizationHotspot"'), "response panel should expose a hotspot readout");
const scriptOrder = [
  "./vendor/three.min.js",
  "./state.js",
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
const browserObj = RAD.exportPaperRadMeshObj(browserMesh);
assert.ok(browserObj.includes("o cell_2_3_outer_plate"), "browser OBJ should include cell object names");
assert.ok(browserObj.includes("# kind connector"), "browser OBJ should include connector metadata");
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
