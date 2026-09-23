import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / "web"


class WebStaticTests(unittest.TestCase):
    def test_index_references_existing_local_assets(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        refs = re.findall(r'(?:src|href)="(\./[^"]+)"', html)
        self.assertIn("./styles.css", refs)
        for ref in refs:
            self.assertTrue((WEB / ref.removeprefix("./")).exists(), ref)

    def test_one_cell_rotation_page_references_existing_local_assets(self):
        page = WEB / "one-cell-rotation" / "index.html"
        html = page.read_text(encoding="utf-8")
        refs = re.findall(r'(?:src|href)="([^"]+)"', html)
        self.assertIn("../vendor/three.min.js", refs)
        self.assertIn("./styles.css", refs)
        self.assertIn("./one_cell_rotation.js", refs)
        self.assertNotIn("cdn.jsdelivr.net", html)
        for ref in refs:
            if ref.startswith("data:"):
                continue
            self.assertTrue((page.parent / ref).resolve().exists(), ref)

    def test_two_cell_attachment_page_references_existing_local_assets(self):
        page = WEB / "two-cell-attachment" / "index.html"
        html = page.read_text(encoding="utf-8")
        refs = re.findall(r'(?:src|href)="([^"]+)"', html)
        self.assertIn("../vendor/three.min.js", refs)
        self.assertIn("./styles.css", refs)
        self.assertIn("./two_cell_attachment.js", refs)
        self.assertNotIn("cdn.jsdelivr.net", html)
        for ref in refs:
            if ref.startswith("data:"):
                continue
            self.assertTrue((page.parent / ref).resolve().exists(), ref)

    def test_two_cell_attachment_page_has_configurable_lattice_controls(self):
        page = WEB / "two-cell-attachment" / "index.html"
        html = page.read_text(encoding="utf-8")
        script = (page.parent / "two_cell_attachment.js").read_text(encoding="utf-8")
        for control_id in [
            "rowCount",
            "colCount",
            "applyGrid",
            "showPins",
            "showCenterLines",
            "fitView",
            "actuationMode",
            "selectedCellMetric",
            "dieOffRadiusMetric",
            "maxPropagatedMetric",
            "latticeSizeMetric",
            "cellCountMetric",
            "pinRadiusRatio",
            "contactTestCommand",
            "pinRadiusMetric",
            "radialClearanceMetric",
            "relativeBacklashMetric",
            "contactStateMetric",
            "transmittedCommandMetric",
        ]:
            self.assertIn(f'id="{control_id}"', html)
        for symbol in [
            "function rebuildLattice",
            "rowCountInput",
            "colCountInput",
            "connections.push",
            "updateLatticeFromPose",
            "fitCurrentView",
            "function annularDisk",
            "new THREE.ExtrudeGeometry",
            "shape.holes.push(hole)",
            "function addCircularHole",
            "function extrudedShapeMesh",
            "addCircularHole(shape, CAD.siteRadiusMm, 0, holeRadius)",
            "function pinRadiusMm",
            "function pinHoleClearanceMm",
            "function normalizedBacklash",
            "function researchAngularBacklashRad",
            "function researchAngularBacklashDeg",
            "function updatePinGeometry",
            "function updateContactReadout",
            "function contactTransmit",
            "function contactDriveField",
            "function explicitConnectionPinPairs",
            "function contactPinPointForPair",
            "function contactPinPoints",
            "function buildLatticeSnapshot",
            "function snapshotHoleRefPoint",
            "function selectCellFromPointer",
            "raycaster.intersectObjects(selectableMeshes",
            "actuationModeInput.value = \"contact\"",
            "const contactMode = actuationModeInput.value === \"contact\" && latestContactField",
            "latestRestSnapshot = contactField ? buildLatticeSnapshot(uniformDriveField(contactField.restDrive), requestedBLower, requestedBTop) : null",
            "const restPinPoint = restFirstPoint.add(restSecondPoint).multiplyScalar(0.5)",
            "const contactMode = actuationModeInput.value === \"contact\"",
            "const range = widestRange(driveEnvelope)",
            "contactMode ? contactDriveField(aTopDeg, driveEnvelope) : null",
            "pinHoleClearanceMm()",
            "Math.asin(normalizedBacklash())",
            "return restPinPoint.add(displacement.multiplyScalar((separation - clearance) / separation))",
            "Math.sign(command) * Math.max(0, Math.abs(command) - backlash)",
            "function collisionReportForBodies",
            "function latticeCollisionReportForDrive",
            "function latticeCollisionReportFromRenderedCells",
            "function mergeCollisionReports",
            "return envelope.ranges.some",
            "LONG_SINGLE_STRAND_MIN_CELLS = 15",
            "LONG_SINGLE_STRAND_MIN_DRIVE_DEG = 33",
            "function latticeMinimumDriveDeg",
            "MIN_CONTINUOUS_CENTER_PITCH_RATIO = 0.75",
            "function isContinuousAssemblyPose",
            "candidateFeasible = continuous",
            "best.continuous && clear && attached",
            "function scanStringContinuity",
            "window.RADAttachmentDiagnostics",
            "maxSnapshotJump",
            "function continuousDriveSamples",
            "MAX_DRIVE_STEP_CENTER_JUMP_MM = 4",
            "MAX_DRIVE_STEP_ROTATION_JUMP_DEG = 5",
            "const feasible = continuousDriveSamples",
            "deltaField",
            "contactRestDriveDeg",
        ]:
            self.assertIn(symbol, script)
        self.assertIn("const floorGrid = new THREE.GridHelper", script)
        self.assertIn("function resizeFloorGridForDimensions", script)
        self.assertLess(script.index("function resizeFloorGridForDimensions"), script.index("function rebuildLattice"))
        self.assertNotIn("floorGrid.scale.setScalar", script[script.index("function updateLatticeFromPose") : script.index("function updateMechanism")])
        for removed_control_id in ["aSite", "bSite", "bTopSite", "aLowerSite"]:
            self.assertNotIn(f'id="{removed_control_id}"', html)
            self.assertNotIn(f'getElementById("{removed_control_id}")', script)
        for removed_control_id in ["enableBacklash", "driveSourceMode", "backlashDeg", "couplingGain", "boundaryMode", "boundaryResponse", "fitWalls"]:
            self.assertNotIn(f'id="{removed_control_id}"', html)
            self.assertNotIn(f'getElementById("{removed_control_id}")', script)
        for removed_symbol in ["function backlashActivation", "computeEffectiveDriveGrid", "applyBoundaryConstraints"]:
            self.assertNotIn(removed_symbol, script)
        self.assertNotIn("Primary A-B Attachment", html)
        self.assertNotIn("A-B Loop Closure", html)
        self.assertNotIn("Backlash Transmission", html)
        self.assertNotIn("Boundary Conditions", html)
        for retained_metric_id in ["residualMetric", "pitchMetric", "aCenterMetric", "bCenterMetric", "cCenterMetric", "systemCenterMetric", "secondResidualMetric", "secondTargetMetric", "secondBMetric"]:
            self.assertIn(f'id="{retained_metric_id}"', html)
        self.assertIn("FIXED_AB_SITES", script)
        self.assertIn("const displayRotation = -safeSolvedAngle", script)
        self.assertIn("centerlineAngle: 0", script)
        self.assertIn("centerPin", script)
        self.assertIn("CAD.nominalHoleDiameterMm * 0.5", script)
        self.assertIn('id="animationHalfPeriodSeconds"', html)
        self.assertIn("animationHalfPeriodSecondsInput", script)
        self.assertIn("ANIMATION_EASE_BLEND = 0.18", script)
        self.assertIn("function mildlyEasedPingPong", script)
        self.assertIn("sineEaseT", script)
        self.assertNotIn("amplitude * Math.sin(timeMs", script)
        for camera_symbol in ["function minCameraRadius", "function maxCameraRadius", "function updateCameraClipPlanes", "camera.far = far"]:
            self.assertIn(camera_symbol, script)
        self.assertIn("clamp(cameraState.radius *", script)

    def test_executable_web_module_validation_script_is_present(self):
        script = (ROOT / "tests" / "validate_web_modules.js").read_text(encoding="utf-8")
        for symbol in [
            "vm.runInContext",
            "RAD.backlashActivation",
            "RAD.serialize",
            "RAD.deserialize",
            "RAD.applyPreset",
            "RAD.exportExperimentSequence",
            "RAD.importExperimentSequence",
            "RAD.exportPaperRadMeshObj",
            "RAD.simulateActive",
            "RAD.lockDatasetSummary",
            "RAD.predictEmpiricalLockCoordinates",
            "RAD.applyEmpiricalLockSurrogate",
            "RAD.Primitives.runSelfCheck",
            "rad-sim.browser-self-check.v1",
            "RAD_LOCK_DATASET",
            "rad-sim.browser-lock-coordinate-dataset.v1",
            "RAD.compareEventOrder",
            "RAD.compareSequenceOrder",
            "RAD.paperRadCalibration",
            "RAD.exportHardwareProfileJson",
            "RAD.importHardwareProfileJson",
            "RAD.calibrationModelProfileResidualComparison",
            "RAD.calibrationModelProfileHoldoutValidation",
            "RAD.calibrationBenchExecutionValidation",
            "RAD.exportCalibrationBenchExecutionValidationCsv",
            "RAD.physicalValidationReadinessReport",
            "RAD.exportPhysicalValidationReadinessCsv",
            "RAD.contactStateAbstractionReport",
            "RAD.exportContactStateAbstractionCsv",
            "RAD.contactGraphConsistencyReport",
            "RAD.exportContactGraphConsistencyCsv",
            "RAD.physicalRealizationMapReport",
            "RAD.exportPhysicalRealizationMapCsv",
            "RAD.externalPhysicsEngineAuditReport",
            "RAD.exportExternalPhysicsEngineAuditCsv",
            "RAD.mujocoModelExportReport",
            "RAD.exportMujocoModelXml",
            "RAD.mujocoPinHoleContactGeometryReport",
            "RAD.exportMujocoPinHoleContactGeometryCsv",
            "RAD.mujocoContactParameterReport",
            "RAD.exportMujocoContactParameterCsv",
            "RAD.contactParameterCalibrationPacket",
            "RAD.exportContactParameterCalibrationPacketCsv",
            "RAD.contactParameterCalibrationResultsTemplate",
            "RAD.compareContactParameterCalibrationResults",
            "RAD.exportContactParameterBenchValidationCsv",
            "RAD.contactParameterIntervalCalibrationReport",
            "RAD.exportContactParameterIntervalCalibrationCsv",
            "RAD.mujocoExternalRunReport",
            "RAD.mujocoExternalComparisonReport",
            "RAD.exportMujocoExternalComparisonCsv",
            "RAD.equilibriumRelationReport",
            "RAD.exportEquilibriumRelationCsv",
            "RAD.reachableEquilibriumControllabilityReport",
            "RAD.exportReachableEquilibriumControllabilityCsv",
            "RAD.reachableEquilibriumBenchProtocol",
            "RAD.exportReachableEquilibriumBenchProtocolCsv",
            "RAD.reachableEquilibriumBenchResultsTemplate",
            "RAD.compareReachableEquilibriumBenchResults",
            "RAD.exportReachableEquilibriumBenchComparisonCsv",
            "RAD.reachableEquilibriumAmplitudeCalibrationReport",
            "RAD.exportReachableEquilibriumAmplitudeCalibrationCsv",
            "RAD.reachableEquilibriumEmpiricalProfileFromAmplitude",
            "RAD.exportReachableEquilibriumEmpiricalProfileCsv",
            "RAD.exportCalibrationModelProfileHoldoutValidationCsv",
            "RAD.selectCalibrationModelProfile",
            "rad-sim.calibration-model-profile-residual-comparison.v1",
            "rad-sim.calibration-model-profile-holdout-validation.v1",
            "rad-sim.calibration-train-holdout-split.v1",
            "rad-sim.calibration-model-profile-selection.v1",
            "rad-sim.hardware-profile.v1",
            "rad-sim.contact-graph-consistency.v1",
            "contact_graph_consistency_gate",
            "contactGraphConsistentNat_zero_missing_evidence",
            "contact_graph_consistency_report",
            "export_contact_graph_consistency_csv",
            "rad-sim.physical-realization-map.v1",
            "physical_realization_map_gate",
            "physicalRealizationMapReadyNat_zero_missing_evidence",
            "physical_realization_map_report",
            "export_physical_realization_map_csv",
            "rad-sim.external-physics-engine-audit.v1",
            "external_physics_engine_audit_gate",
            "externalPhysicsEngineAuditReadyNat_zero_missing_evidence",
            "external_physics_engine_audit_report",
            "export_external_physics_engine_audit_csv",
            "rad-sim.mujoco-model-export.v1",
            "mujoco_model_export_gate",
            "externalPhysicsModelExportReadyNat_zero_missing_evidence",
            "mujoco_model_export_report",
            "export_mujoco_model_xml",
            "rad-sim.mujoco-pin-hole-contact-geometry.v1",
            "mujoco_pin_hole_contact_geometry_gate",
            "externalContactGeometryReadyNat_zero_missing_evidence",
            "mujoco_pin_hole_contact_geometry_report",
            "export_mujoco_pin_hole_contact_geometry_csv",
            "rad-sim.mujoco-contact-parameter-profile.v1",
            "mujoco_contact_parameter_profile_gate",
            "externalContactParameterProfileReadyNat_zero_missing_evidence",
            "mujoco_contact_parameter_report",
            "export_mujoco_contact_parameter_csv",
            "rad-sim.contact-parameter-calibration-packet.v1",
            "contact_parameter_calibration_packet_completeness",
            "contactParameterCalibrationPacketCompleteNat_zero_missing_evidence",
            "contact_parameter_calibration_packet",
            "export_contact_parameter_calibration_packet_csv",
            "rad-sim.contact-parameter-bench-validation.v1",
            "contact_parameter_bench_validation_gate",
            "contactParameterBenchValidationReadyNat_zero_missing_evidence",
            "compare_contact_parameter_calibration_results",
            "export_contact_parameter_bench_validation_csv",
            "rad-sim.contact-parameter-interval-calibration.v1",
            "contact_parameter_interval_calibration_gate",
            "contactParameterIntervalCalibrationReadyNat_zero_missing_evidence",
            "contact_parameter_interval_calibration_report",
            "export_contact_parameter_interval_calibration_csv",
            "rad-sim.mujoco-external-run.v1",
            "mujoco_external_run_gate",
            "externalPhysicsRunReadyNat_zero_missing_evidence",
            "mujoco_external_run_report",
            "rad-sim.mujoco-external-comparison.v1",
            "mujoco_external_comparison_gate",
            "externalPhysicsComparisonReadyNat_zero_missing_evidence",
            "mujoco_external_comparison_report",
            "export_mujoco_external_comparison_csv",
            "rad-sim.equilibrium-relation.v1",
            "equilibrium_relation_gate",
            "equilibriumRelationReadyNat_zero_missing_evidence",
            "equilibrium_relation_report",
            "export_equilibrium_relation_csv",
            "rad-sim.reachable-equilibrium-controllability.v1",
            "reachable_equilibrium_controllability_gate",
            "reachableEquilibriumControllabilityReadyNat_zero_missing_evidence",
            "reachable_equilibrium_controllability_report",
            "export_reachable_equilibrium_controllability_csv",
            "rad-sim.reachable-equilibrium-bench-protocol.v1",
            "reachable_equilibrium_bench_protocol_gate",
            "reachableEquilibriumBenchProtocolReadyNat_zero_missing_evidence",
            "reachable_equilibrium_bench_protocol",
            "export_reachable_equilibrium_bench_protocol_csv",
            "rad-sim.reachable-equilibrium-bench-results.v1",
            "rad-sim.reachable-equilibrium-bench-comparison.v1",
            "reachable_equilibrium_bench_validation_gate",
            "reachableEquilibriumBenchValidationReadyNat_zero_missing_evidence",
            "compare_reachable_equilibrium_bench_results",
            "export_reachable_equilibrium_bench_comparison_csv",
            "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
            "reachable_equilibrium_amplitude_calibration_gate",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_zero_missing_evidence",
            "reachable_equilibrium_amplitude_calibration_report",
            "export_reachable_equilibrium_amplitude_calibration_csv",
            "rad-sim.reachable-equilibrium-empirical-profile.v1",
            "reachable_equilibrium_empirical_profile_gate",
            "reachableEquilibriumEmpiricalProfileReadyNat_zero_missing_evidence",
            "reachable_equilibrium_empirical_profile_from_amplitude",
            "export_reachable_equilibrium_empirical_profile_csv",
            "RAD.characterizeLocalResponse",
            "RAD.physicalPreviewComparison",
            "RAD.responseDecayProfile",
            "alphaDecayRatio",
            "zDecayRatio",
            "alphaDecayLength",
            "zDecayLength",
            "physicalHeightRmsError",
            "physicalHeightMaxError",
            "physicalPreviewAvailable",
            "modelErrorHeight",
            "physicalMaxHeightDelta",
            "physicalRmsCenterDelta",
            "localRefs",
            "validateLocalHttpEntry",
            "http.createServer",
            "fetch(`${base}/index.html`)",
            "asset should load over local HTTP",
            "vendored Three.js should expose THREE",
            "new THREE.Scene()",
            "new THREE.PerspectiveCamera",
            "new THREE.OrthographicCamera",
            "new THREE.Mesh",
            "THREE.Raycaster",
            "./vendor/three.min.js",
            "./primitives.js",
            "./data/lock_dataset.js",
            "./lock_surrogate.js",
            "browser entry should not require external scripts or styles",
            "scriptOrder",
            "state.grid.backlash",
            "state.grid.couplingGain",
            "state.grid.zCouplingGain",
            "state.grid.empiricalLockModel",
            "state.grid.lockStateIndex",
            "state.grid.paperSideLengthMm",
            "state.grid.paperHoleToleranceMm",
            "state.view.overlayMode",
            "state.view.simulationMode",
            "state.view.membraneVisible",
            "state.experiment.presetName",
            "web module validation passed",
        ]:
            self.assertIn(symbol, script)

    def test_local_threejs_is_declared_before_renderer(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn("./vendor/three.min.js", html)
        self.assertTrue((WEB / "vendor" / "three.min.js").exists())
        self.assertNotIn("cdn.jsdelivr.net", html)
        self.assertLess(html.index("three.min.js"), html.index("./renderer.js"))
        self.assertLess(html.index("./math.js"), html.index("./operators.js"))
        self.assertLess(html.index("./math.js"), html.index("./constraints.js"))
        self.assertLess(html.index("./constraints.js"), html.index("./operators.js"))
        self.assertLess(html.index("./operators.js"), html.index("./inverse.js"))
        self.assertLess(html.index("./math.js"), html.index("./analysis.js"))
        self.assertLess(html.index("./analysis.js"), html.index("./data/lock_dataset.js"))
        self.assertLess(html.index("./data/lock_dataset.js"), html.index("./data/two_cell_external_fidelity_summary.js"))
        self.assertLess(html.index("./data/two_cell_external_fidelity_summary.js"), html.index("./lock_surrogate.js"))
        self.assertLess(html.index("./data/lock_dataset.js"), html.index("./lock_surrogate.js"))
        self.assertLess(html.index("./lock_surrogate.js"), html.index("./physics.js"))
        self.assertLess(html.index("./analysis.js"), html.index("./physics.js"))
        self.assertLess(html.index("./physics.js"), html.index("./two_cell_bench.js"))
        self.assertLess(html.index("./two_cell_bench.js"), html.index("./mesh_export.js"))
        self.assertLess(html.index("./mesh_export.js"), html.index("./renderer.js"))
        self.assertLess(html.index("./ui.js"), html.index("./primitives.js"))
        self.assertLess(html.index("./primitives.js"), html.index("./app.js"))

    def test_browser_modules_expose_expected_api(self):
        expected = {
            "state.js": ["RAD.createState", "RAD.serialize", "RAD.updateDerivedCells", "RAD.exportExperimentSequence", "RAD.importExperimentSequence", "RAD.deserialize"],
            "operators.js": ["RAD.localActuationEvent", "RAD.lockEvent", "RAD.applyEventSequence", "RAD.compareEventOrder", "RAD.compareSequenceOrder", "RAD.finiteDieOffRadius"],
            "analysis.js": ["RAD.analyzeExperimentSequence", "RAD.exportSequenceMetricsCsv", "RAD.sequenceFrames", "RAD.characterizeLocalResponse", "RAD.physicalPreviewComparison", "RAD.responseDecayProfile", "RAD.calibrationParameterEstimates", "RAD.calibrationModelProfile", "RAD.applyCalibrationModelProfile", "RAD.exportCalibrationModelProfile", "RAD.calibrationModelProfileResidualComparison", "RAD.calibrationModelProfileHoldoutValidation", "RAD.calibrationBenchExecutionValidation", "RAD.exportCalibrationBenchExecutionValidation", "RAD.exportCalibrationBenchExecutionValidationCsv", "RAD.physicalValidationReadinessReport", "RAD.exportPhysicalValidationReadiness", "RAD.exportPhysicalValidationReadinessCsv", "RAD.contactStateAbstractionReport", "RAD.exportContactStateAbstraction", "RAD.exportContactStateAbstractionCsv", "RAD.contactGraphConsistencyReport", "RAD.exportContactGraphConsistency", "RAD.exportContactGraphConsistencyCsv", "RAD.physicalRealizationMapReport", "RAD.exportPhysicalRealizationMap", "RAD.exportPhysicalRealizationMapCsv", "RAD.externalPhysicsEngineAuditReport", "RAD.exportExternalPhysicsEngineAudit", "RAD.exportExternalPhysicsEngineAuditCsv", "RAD.mujocoModelExportReport", "RAD.exportMujocoModelXml", "RAD.exportMujocoModelReport", "RAD.mujocoPinHoleContactGeometryReport", "RAD.exportMujocoPinHoleContactGeometry", "RAD.exportMujocoPinHoleContactGeometryCsv", "RAD.mujocoContactParameterReport", "RAD.exportMujocoContactParameter", "RAD.exportMujocoContactParameterCsv", "RAD.contactParameterCalibrationPacket", "RAD.exportContactParameterCalibrationPacket", "RAD.exportContactParameterCalibrationPacketCsv", "RAD.contactParameterCalibrationResultsTemplate", "RAD.exportContactParameterCalibrationResultsTemplate", "RAD.contactParameterCalibrationResultsFromJson", "RAD.compareContactParameterCalibrationResults", "RAD.exportContactParameterBenchValidation", "RAD.exportContactParameterBenchValidationCsv", "RAD.contactParameterIntervalCalibrationReport", "RAD.exportContactParameterIntervalCalibration", "RAD.exportContactParameterIntervalCalibrationCsv", "RAD.mujocoExternalRunReport", "RAD.exportMujocoExternalRun", "RAD.mujocoExternalComparisonReport", "RAD.exportMujocoExternalComparison", "RAD.exportMujocoExternalComparisonCsv", "RAD.equilibriumRelationReport", "RAD.exportEquilibriumRelation", "RAD.exportEquilibriumRelationCsv", "RAD.reachableEquilibriumControllabilityReport", "RAD.exportReachableEquilibriumControllability", "RAD.exportReachableEquilibriumControllabilityCsv", "RAD.reachableEquilibriumBenchProtocol", "RAD.exportReachableEquilibriumBenchProtocol", "RAD.exportReachableEquilibriumBenchProtocolCsv", "RAD.reachableEquilibriumBenchResultsTemplate", "RAD.exportReachableEquilibriumBenchResultsTemplate", "RAD.reachableEquilibriumBenchResultsFromJson", "RAD.compareReachableEquilibriumBenchResults", "RAD.exportReachableEquilibriumBenchComparison", "RAD.exportReachableEquilibriumBenchComparisonCsv", "RAD.reachableEquilibriumAmplitudeCalibrationReport", "RAD.exportReachableEquilibriumAmplitudeCalibration", "RAD.exportReachableEquilibriumAmplitudeCalibrationCsv", "RAD.reachableEquilibriumEmpiricalProfileFromAmplitude", "RAD.exportReachableEquilibriumEmpiricalProfile", "RAD.exportReachableEquilibriumEmpiricalProfileCsv", "RAD.selectCalibrationModelProfile", "RAD.exportCalibrationModelProfileResidualComparison", "RAD.exportCalibrationModelProfileHoldoutValidation", "RAD.exportCalibrationModelProfileHoldoutValidationCsv", "RAD.exportCalibrationModelProfileSelection"],
            "lock_surrogate.js": ["RAD.lockDatasetSummary", "RAD.predictEmpiricalLockCoordinates", "RAD.applyEmpiricalLockSurrogate"],
            "physics.js": ["RAD.simulatePhysicalRelaxation", "RAD.simulateActive"],
            "constraints.js": ["RAD.solveBoundaryConstraints", "RAD.fitBoundaryToReference", "RAD.clearBoundaryConstraints", "rad-sim.constraint-realization.v1"],
            "two_cell_bench.js": ["RAD.simulateTwoCellBench", "RAD.sweepTwoCellBacklash", "RAD.exportTwoCellBenchJson", "RAD.exportTwoCellBacklashSweepCsv", "RAD.twoCellPhysicalSimulationSuite", "RAD.exportTwoCellPhysicalSuiteCsv", "RAD.twoCellPhysicalFidelityMatrix", "RAD.exportTwoCellPhysicalFidelityMatrixCsv", "RAD.twoCellContactPhaseMap", "RAD.exportTwoCellContactPhaseMapCsv", "RAD.twoCellRadiusBacklashPhaseDiagram", "RAD.exportTwoCellRadiusBacklashPhaseDiagramCsv", "RAD.twoCellRadiusBacklashTransitionReport", "RAD.exportTwoCellRadiusBacklashTransitionReportCsv", "RAD.compareTwoCellRadiusBacklashTransitionMeasurements", "RAD.exportTwoCellRadiusBacklashTransitionComparisonCsv", "RAD.twoCellRadiusBacklashTransitionRerun", "RAD.exportTwoCellRadiusBacklashTransitionRerunCsv", "RAD.twoCellPhysicalResponseAtlas", "RAD.exportTwoCellPhysicalResponseAtlasCsv", "RAD.twoCellCadContactDecompositionSpec", "RAD.exportTwoCellCadContactDecompositionJson", "RAD.exportTwoCellCadContactDecompositionCsv", "RAD.twoCellExactContactHandoffPlan", "RAD.exportTwoCellExactContactHandoffPlanJson", "RAD.exportTwoCellExactContactHandoffPlanCsv", "RAD.twoCellFidelityMatrixMeasurementTemplate", "RAD.exportTwoCellFidelityMatrixMeasurementTemplateCsv", "RAD.twoCellExternalFidelityMatrixManifest", "RAD.exportTwoCellExternalFidelityMatrixManifestJson", "RAD.exportTwoCellExternalFidelityMatrixManifestCsv", "RAD.twoCellFidelityMatrixMeasurementsFromCsv", "RAD.compareTwoCellFidelityMatrixMeasurements", "RAD.exportTwoCellFidelityMatrixMeasurementComparisonJson", "RAD.calibrateTwoCellFidelityMatrixParameters", "RAD.exportTwoCellFidelityMatrixParameterCalibrationJson", "RAD.cadRadCellArchiveAudit", "RAD.cadRadCellReferenceProfile", "RAD.exportCadRadCellReferenceProfileJson", "RAD.twoCellExternalFidelityWebSummary", "RAD.twoCellExternalFidelityCaseOptions", "RAD.twoCellExternalFidelityCaseControls", "RAD.exportTwoCellExternalFidelityWebSummaryJson", "RAD.twoCellPhysicalFidelityStatus", "RAD.exportTwoCellPhysicalFidelityStatusJson", "RAD.twoCellConnectorMeasurementTemplate", "RAD.exportTwoCellConnectorMeasurementTemplateCsv", "rad-sim.two-cell-physical-bench.v1", "rad-sim.two-cell-contact-phase-map.v1", "rad-sim.two-cell-radius-backlash-phase-diagram.v1", "rad-sim.two-cell-radius-backlash-transition-report.v1", "rad-sim.two-cell-radius-backlash-transition-comparison.v1", "rad-sim.two-cell-radius-backlash-transition-rerun.v1", "rad-sim.two-cell-physical-response-atlas.v1", "rad-sim.two-cell-cad-contact-decomposition.v1", "rad-sim.two-cell-exact-contact-handoff-plan.v1", "rad-sim.two-cell-fidelity-matrix-measurement-template.v1", "rad-sim.browser-two-cell-external-fidelity-matrix-manifest.v1", "rad-sim.browser-two-cell-fidelity-matrix-measurement-comparison.v1", "rad-sim.browser-two-cell-fidelity-matrix-parameter-calibration.v1", "rad-sim.cad-rad-cell-archive-audit.v1", "rad-sim.cad-rad-cell-reference-profile.v1", "rad-sim.browser-two-cell-physical-fidelity-status.v1", "rad-sim.web-two-cell-external-fidelity-summary.v1", "rad-sim.browser-two-cell-external-fidelity-case-preview.v1"],
            "mesh_export.js": ["RAD.buildPaperRadMesh", "RAD.exportPaperRadMeshObj"],
            "math.js": [
                "RAD.backlashActivation",
                "RAD.alphaToTheta",
                "RAD.normalizedBacklashToThetaDeadZone",
                "RAD.normalizedBacklashToAlphaDeadZone",
                "RAD.compileCustomTargetExpression",
                "RAD.simulate",
                "RAD.applyPreset",
                "RAD.optimizeCommandsForTarget",
                "RAD.pairDistanceFromTheta",
                "RAD.solveConstraintKinematicApprox",
                "RAD.paperRadCalibration",
                "RAD.exportHardwareProfileJson",
                "RAD.importHardwareProfileJson",
                "RAD.modelLengthToMm",
                "RAD.mmToModelLength",
            ],
            "inverse.js": ["RAD.buildInverseDesignPlan", "RAD.applyInverseDesignPlan"],
            "renderer.js": ["RAD.RadRenderer"],
            "ui.js": ["RAD.RadUI"],
            "primitives.js": ["RAD.Primitives", "runSelfCheck", "installSelfCheckPanel", "installCollapsibleSections", "rad-sim.browser-self-check.v1"],
            "app.js": ["RAD_APP", "runSelfCheck"],
        }
        for filename, names in expected.items():
            text = (WEB / filename).read_text(encoding="utf-8")
            for name in names:
                self.assertIn(name, text, f"{filename} should expose {name}")

    def test_app_exposes_browser_validation_handles(self):
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in ["get state()", "get simulation()", "get renderer()", "update: (nextState = state) => renderAll(nextState)", "RAD.Primitives.installSelfCheckPanel", "runSelfCheck"]:
            self.assertIn(symbol, app_js)

    def test_browser_primitive_library_and_self_check_ui_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        css = (WEB / "styles.css").read_text(encoding="utf-8")
        primitives_js = (WEB / "primitives.js").read_text(encoding="utf-8")
        for control_id in ["runSystemCheck", "systemCheckPanel", "systemCheckSummary", "systemCheckList"]:
            self.assertIn(f'id="{control_id}"', html)
        self.assertIn("./primitives.js", html)
        for symbol in ["CORE_DOM_IDS", "CORE_APIS", "installSelfCheckPanel", "runSelfCheck", "Event lock preserves height", "Position lock fixture", "Preset diversity"]:
            self.assertIn(symbol, primitives_js)
        for symbol in [".system-check-panel", ".system-check-panel.is-ok", ".system-check-panel.is-failing", "body.system-check-ok #runSystemCheck", "body.system-check-failing #runSystemCheck"]:
            self.assertIn(symbol, css)

    def test_selected_cell_actuation_controls_are_viewport_docked(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        css = (WEB / "styles.css").read_text(encoding="utf-8")
        self.assertIn('class="quick-actuation-dock"', html)
        self.assertIn('id="toggleQuickDock"', html)
        self.assertIn('id="quickDockMini"', html)
        self.assertIn('class="quick-dock-body"', html)
        self.assertIn('title="Checkerboard pattern"', html)
        self.assertIn(">Check</button>", html)
        self.assertLess(html.index('id="threeMount"'), html.index('class="quick-actuation-dock"'))
        self.assertLess(html.index('class="quick-actuation-dock"'), html.index('id="overlayLegend"'))
        self.assertIn('id="alphaCommand"', html)
        self.assertIn('id="zCommand"', html)
        self.assertIn('id="zCoupling"', html)
        self.assertIn('id="applyCommand"', html)
        self.assertIn('id="clearCell"', html)
        for control_id in ["nudgeAlphaContract", "nudgeAlphaExpand", "nudgeZDown", "nudgeZUp", "zeroSelectedCommand"]:
            self.assertIn(f'id="{control_id}"', html)
        self.assertIn('class="quick-preset-strip"', html)
        self.assertIn(".quick-actuation-dock", css)
        self.assertIn(".quick-actuation-dock.is-collapsed", css)
        self.assertIn(".quick-actuation-dock.is-collapsed .quick-dock-body", css)
        self.assertIn(".quick-dock-mini", css)
        self.assertIn("white-space: nowrap", css)
        self.assertIn("position: absolute", css)
        self.assertIn("top: calc(var(--viewport-toolbar-height) + 12px)", css)
        self.assertIn(".quick-toggle-row", css)
        self.assertIn(".quick-nudge-grid", css)
        self.assertIn(".quick-preset-strip", css)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn("quickDockCollapsed: false", state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["toggleQuickDock", "quickDockMini", "syncQuickDock", "quickDockCollapsed", "is-collapsed"]:
            self.assertIn(symbol, ui_js)

    def test_response_experiments_show_physical_preview_comparison(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        analysis_js = (WEB / "analysis.js").read_text(encoding="utf-8")
        for control_id in ["characterizationPhysical", "characterizationPhysicalMax", "characterizationDecay", "characterizationDecayLength"]:
            self.assertIn(f'id="{control_id}"', html)
        for symbol in [
            "physicalPreviewAvailable",
            "physicalHeightRmsError",
            "physicalHeightMaxError",
            "alphaDecayRatio",
            "zDecayRatio",
            "alphaDecayLength",
            "zDecayLength",
        ]:
            self.assertIn(symbol, ui_js)
            self.assertIn(symbol, analysis_js)
        for label in ["phys rms", "phys max", "decay a", "len a"]:
            self.assertIn(label, ui_js)
        for symbol in ["responseDecayProfile", "log-linear shell max", "fitShellDecay", "nearestSourceDistance"]:
            self.assertIn(symbol, analysis_js)

    def test_desktop_layout_keeps_lattice_visible_while_controls_scroll(self):
        css = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [
            "body {\n  margin: 0;\n  min-height: 100vh;\n  overflow: hidden;",
            ".app-shell {\n  display: grid;\n  grid-template-columns: minmax(0, 1fr) clamp(360px, 27vw, 430px);\n  height: 100vh;",
            ".workspace {\n  min-width: 0;\n  height: 100vh;",
            "display: flex;\n  flex-direction: column;\n  overflow: hidden;",
            ".viewport {\n  position: relative;\n  --viewport-toolbar-height: 82px;\n  flex: 1 1 auto;",
            ".three-mount {\n  width: 100%;\n  height: 100%;",
            ".control-panel {\n  height: 100vh;",
            "overflow-y: auto;\n  overscroll-behavior: contain;",
            "border-radius: 8px;",
            "box-shadow: var(--shadow);",
            "@media (max-width: 1100px)",
            "body {\n    overflow: auto;",
            ".control-panel {\n    height: auto;",
        ]:
            self.assertIn(symbol, css)

    def test_paint_clicks_apply_current_dock_command_from_3d_selection(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="paintMode"', html)
        self.assertIn("Paint Clicks Off", html)
        self.assertIn('id="paintRadius"', html)
        self.assertIn('id="paintRadiusOut"', html)
        css = (WEB / "styles.css").read_text(encoding="utf-8")
        self.assertIn(".paint-mode-button", css)
        self.assertIn(".paint-mode-button.is-active", css)
        self.assertIn(".paint-radius-control", css)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn("paintMode: false", state_js)
        self.assertIn("paintRadius: 0", state_js)
        self.assertIn("state.view.paintMode === undefined", state_js)
        self.assertIn("state.view.paintRadius === undefined", state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in [
            "paintMode",
            "paintRadius",
            "paintCell(r, c)",
            "paintBrushCells",
            "RAD.brushCells",
            "cell-paint-command",
            "affected: cells.length",
            "radius: Number(this.state.view.paintRadius || 0)",
            "syncPaintMode",
            "Paint Clicks On",
            "Paint Clicks Off",
            "RAD.clampCommandAlpha",
            "RAD.clampCommandZ",
        ]:
            self.assertIn(symbol, ui_js)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        self.assertIn("state.view.paintMode", app_js)
        self.assertIn("ui.paintCell(r, c)", app_js)
        self.assertIn("ui.select(r, c)", app_js)
        self.assertIn("Paint clicks apply current command", app_js)
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        for symbol in ["function brushCells", "RAD.brushCells", "Math.abs(rr - r) + Math.abs(cc - c) <= radius"]:
            self.assertIn(symbol, math_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in [
            "paintBrushPreviewRoot",
            "paintBrushPreview",
            "renderPaintBrushPreview",
            "addBrushPreviewRing",
            "state.view.paintMode !== true",
            "RAD.brushCells(state, hover.r, hover.c, state.view.paintRadius)",
        ]:
            self.assertIn(symbol, renderer_js)

    def test_camera_view_presets_are_axis_explicit(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for mode in ["iso", "top", "front", "side"]:
            self.assertIn(f'data-view-mode="{mode}"', html)
        self.assertIn('id="isolateCell"', html)
        self.assertIn('id="explodeCell"', html)
        self.assertIn('id="projectionMode"', html)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in ["setActiveView", "aria-pressed", "renderer.viewMode", "resetCameraView", "setFocusMode", "is-focus-mode", "syncCameraState", "frameSelectedCell", "setProjectionMode", "updateProjectionButton", "projectionMode", "orthographic", "Perspective", "Ortho", "setIsolateCell", "isolateCell", "Show Lattice", "Selected cell isolated", "setExplodeCell", "explodeCell", "Exploded cell detail active", "Assemble Cell", "hoveredCell", "Hover: r"]:
            self.assertIn(symbol, app_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in [
            "this.perspectiveCamera",
            "this.orthographicCamera",
            "new this.THREE.OrthographicCamera",
            "this.projectionMode = \"perspective\"",
            "setProjectionMode(mode)",
            "updateProjectionMatrices",
            "this.projectionMode === \"orthographic\"",
            "this.viewMode = \"iso\"",
            "setCustomView",
            "notifyViewChange",
            "resetView",
            "this.cameraTarget.set(0, 0, 0)",
            "getCameraState",
            "applyCameraState",
            "frameCell",
            "sim?.centers",
            "isSelectedCell",
            "isCellVisible",
            "explodedCellAmount",
            "explodedPlatePosition",
            "partCalloutRoot",
            "partCallout",
            "renderExplodedPartCallouts",
            "addPartCallout",
            "makePartCalloutLabel",
            "rotating plate",
            "hinge pin",
            "backlash gap",
            "z actuator",
            "alpha actuator",
            "state.view.isolateSelected !== true",
            "state.view.isolateSelected === true",
            "state.view.explodedSelected === true",
            "this.clearGroup(this.surfaceContourRoot)",
            "panCamera",
            "hoverPick",
            "setHoveredCell",
            "this.hoveredCell",
            "this.materials.hovered",
            "onHover",
            'mode === "pan"',
            "event.shiftKey",
            "event.button === 1",
            "event.button === 2",
            "contextmenu",
            "setCustomView(true)",
            "this.perspectiveCamera.up.set(0, 0, 1)",
            "this.orthographicCamera.up.set(0, 0, 1)",
            "applyCameraUpForCurrentView",
            "this.viewMode === \"top\"",
            "this.perspectiveCamera.up.set(0, 1, 0)",
            "this.orthographicCamera.up.set(0, 1, 0)",
            "setCameraPosition(9, -9, 9)",
            "setCameraPosition(0, -0.001, 13)",
            "setCameraPosition(0, -13, 0.001)",
            "setCameraPosition(13, 0, 0.001)",
            "addAxisLabel(\"X\"",
            "addAxisLabel(\"Y\"",
            "addAxisLabel(\"Z\"",
            "addScaleMarkers",
            "normalized cell units",
            "scaleMarker",
            "makeSceneLabel",
        ]:
            self.assertIn(symbol, renderer_js)
        fit_root_body = renderer_js.split("fitRoot(rows, cols) {", 1)[1].split("\n    pick(event)", 1)[0]
        self.assertNotIn("cameraTarget.set(0, 0, 0)", fit_root_body)
        self.assertIn('id="frameCell"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn("isolateSelected: false", state_js)
        self.assertIn("explodedSelected: false", state_js)
        self.assertIn("state.view.isolateSelected === undefined", state_js)
        self.assertIn("state.view.explodedSelected === undefined", state_js)

    def test_viewport_orientation_hud_reports_camera_basis(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for symbol in [
            'id="orientationHud"',
            'id="cameraViewAxis"',
            'id="cameraUpAxis"',
            'id="cameraProjection"',
            "axis-chip axis-x",
            "axis-chip axis-y",
            "axis-chip axis-z",
        ]:
            self.assertIn(symbol, html)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [
            ".orientation-hud",
            "pointer-events: none",
            ".orientation-axis-row",
            ".axis-x",
            ".axis-y",
            ".axis-z",
        ]:
            self.assertIn(symbol, styles)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in [
            "getOrientationState",
            "formatDominantAxis",
            "viewAxis",
            "upAxis",
            "this.camera.position",
            "this.cameraTarget",
            "this.projectionMode === \"orthographic\" ? \"ortho\" : \"persp\"",
        ]:
            self.assertIn(symbol, renderer_js)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in [
            "orientationHud",
            "cameraViewAxis",
            "cameraUpAxis",
            "cameraProjection",
            "updateOrientationHud",
            "renderer.getOrientationState()",
        ]:
            self.assertIn(symbol, app_js)

    def test_browser_state_serializes_full_camera_pose(self):
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in [
            "camera: {",
            'projection: "perspective"',
            "radius: 15.6",
            "theta: -Math.PI / 4",
            "Math.acos(1 / Math.sqrt(3))",
            "target: { x: 0, y: 0, z: 0 }",
            "state.view.camera",
            "state.view.camera.projection",
        ]:
            self.assertIn(symbol, state_js)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in ["state.view.camera = renderer.getCameraState()", "renderer.applyCameraState(savedCamera)", "syncCameraState", "updateProjectionButton"]:
            self.assertIn(symbol, app_js)

    def test_browser_state_schema_includes_derived_cell_fields(self):
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in [
            "derivedCellFields",
            "updateDerivedCells",
            "state.cells.alpha",
            "state.cells.theta",
            "state.cells.z",
            "sim.alpha",
            "sim.theta",
            "sim.height",
            "alpha: interpolateMatrix",
            "theta: interpolateMatrix",
            "z: interpolateMatrix",
            "stateWithFreshDerivedCells",
        ]:
            self.assertIn(symbol, state_js)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        self.assertIn("RAD.updateDerivedCells(state, sim)", app_js)
        script = (ROOT / "tests" / "validate_web_modules.js").read_text(encoding="utf-8")
        for symbol in [
            "serialized cells should include derived alpha",
            "serialized cells should include derived theta",
            "serialized cells should include derived z",
        ]:
            self.assertIn(symbol, script)

    def test_focus_mode_workspace_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="focusMode"', html)
        self.assertIn('aria-pressed="false"', html)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [
            "body.is-focus-mode .app-shell",
            "body.is-focus-mode .control-panel",
            "body.is-focus-mode .workspace",
            "body.is-focus-mode .topbar",
            "body.is-focus-mode .viewport",
        ]:
            self.assertIn(symbol, styles)

    def test_metric_density_toggle_is_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="metricDensity"', html)
        self.assertIn("advanced-metric", html)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [".advanced-metric", "body.show-all-metrics .advanced-metric", ".metric-density"]:
            self.assertIn(symbol, styles)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in ["setMetricDensity", "show-all-metrics", "allMetricsVisible", "Core Metrics", "All Metrics"]:
            self.assertIn(symbol, app_js)

    def test_workspace_mode_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for mode in ["inspect", "edit", "sheet", "analyze"]:
            self.assertIn(f'data-workspace-mode="{mode}"', html)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in [
            "workspaceMode",
            "setActiveWorkspaceMode",
            "syncViewChrome",
            "is-sheet-view",
            "applyWorkspaceMode",
            "modeInspect",
            "modeEdit",
            "modeSheet",
            "modeAnalyze",
            "Workspace:",
            "targetAvailable",
            'workspaceMode === "sheet"',
            'cellVisualMode: "sheetOnly"',
            "actuatorDisplayMode: \"planned\"",
            "overlayMode: targetAvailable ? \"error\" : \"height\"",
        ]:
            self.assertIn(symbol, app_js)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        self.assertIn("[data-workspace-mode]", styles)
        self.assertIn('body.is-sheet-view .quick-actuation-dock', styles)

    def test_control_panel_sections_are_collapsible(self):
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in [
            "setupControlPanelSections",
            "initiallyCollapsed",
            "section-toggle",
            "is-collapsed",
            "aria-expanded",
            "Target Surface",
            "Inverse Plan",
            "Timeline",
        ]:
            self.assertIn(symbol, app_js)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [".section-heading", ".section-toggle", ".control-section.is-collapsed"]:
            self.assertIn(symbol, styles)

    def test_overlay_legend_is_present_and_dynamic(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["overlayLegend", "overlayLegendTitle", "overlayLegendMin", "overlayLegendMax"]:
            self.assertIn(f'id="{control_id}"', html)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in ["overlayLegendText", "updateOverlayLegend", "Target Error", "Reachability", "Cell State", "Z Residual", "Model Disagreement"]:
            self.assertIn(symbol, app_js)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [".overlay-legend", ".legend-scale", ".legend-labels"]:
            self.assertIn(symbol, styles)

    def test_selected_cell_viewport_hud_is_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["selectionHud", "hudCell", "hudAlpha", "hudTheta", "hudHeight", "hudCommand", "hudState"]:
            self.assertIn(f'id="{control_id}"', html)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in ["selectionHud", "updateSelectionHud", "Cell r", "commandAlpha", "commandZ", "masked"]:
            self.assertIn(symbol, app_js)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [".selection-hud", ".selection-hud dl", ".selection-hud dd"]:
            self.assertIn(symbol, styles)
        self.assertIn("pointer-events: none", styles)
        self.assertIn("z-index: 6", styles)
        self.assertIn("--viewport-toolbar-height", styles)
        self.assertIn("bottom: 16px", styles)
        self.assertIn("width: min(520px, calc(100% - 390px))", styles)
        self.assertIn(".viewport-header .view-buttons", styles)
        self.assertIn("overflow-x: auto", styles)

    def test_deadzone_formula_matches_reference_shape(self):
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        self.assertIn("Math.max(0, x - backlash) + Math.min(x + backlash, 0)", math_js)

    def test_presets_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for preset in ["center", "dome", "saddle", "ridge", "wave", "corner", "twist", "ring", "checker", "clear"]:
            self.assertIn(f'data-preset="{preset}"', html)

    def test_inverse_and_timeline_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in [
            "analyzeInversePlan",
            "analyzeSensitivity",
            "buildJacobian",
            "solveLinearFit",
            "applyInversePlan",
            "applyLinearFit",
            "seedTargetFit",
            "optimizeTargetFit",
            "maxActuators",
            "optimizerIterations",
            "actuatorPenalty",
            "travelPenalty",
            "zTravelLimit",
            "alphaContractLimit",
            "alphaExpandLimit",
            "actuatorDisplayMode",
            "actuatorAllowed",
            "allowAllActuators",
            "blockRimActuators",
            "allowCenterActuators",
            "invertActuatorMask",
            "customTargetExpression",
            "targetExpressionStatus",
            "timelineIndex",
            "playTimeline",
            "stepTimeline",
            "transitionMs",
            "smoothTimeline",
            "keyframeName",
            "captureKeyframe",
            "firstKeyframe",
            "prevKeyframe",
            "nextKeyframe",
            "lastKeyframe",
            "saveSequence",
            "loadSequence",
            "analyzeSequence",
            "jumpBestFrame",
            "jumpMaxSaturationFrame",
            "selectWorstResidualCell",
            "saveSequenceCsv",
        ]:
            self.assertIn(f'id="{control_id}"', html)

    def test_custom_target_surface_controls_and_math_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('value="custom"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["customExpression", "expressionError"]:
            self.assertIn(symbol, state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("this.els.customTargetExpression.addEventListener", ui_js)
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        for symbol in ["compileCustomTargetExpression", "CUSTOM_ALLOWED_IDENTIFIERS", "Unsupported target identifier", 'type === "custom"', "clamp(value, -1.4, 1.4)"]:
            self.assertIn(symbol, math_js)

    def test_error_and_travel_overlays_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="cellVisualMode"', html)
        self.assertIn('id="simulationMode"', html)
        self.assertIn('value="abstract"', html)
        self.assertIn('value="graph"', html)
        self.assertIn('value="sheetOnly"', html)
        self.assertIn('value="paperRad"', html)
        self.assertIn('value="calibratedRad"', html)
        self.assertIn('value="mechanism"', html)
        self.assertIn('value="springPreview"', html)
        self.assertIn('value="constraintSolved"', html)
        for constraint_id in ["constraintSolverStatusOut", "constraintMaxEdgeErrorOut", "constraintMeanEdgeErrorOut", "constraintContactCountOut"]:
            self.assertIn(f'id="{constraint_id}"', html)
        for option in ['value="zresidual"', 'value="error"', 'value="travel"', 'value="saturation"', 'value="strain"', 'value="modelError"', 'value="displacement"', 'value="slope"', 'value="inverse"', 'value="sensitivity"', 'value="reachability"', 'value="topologyBlocked"']:
            self.assertIn(option, html)
        for metric_id in ["recommendedActuators", "candidateActuators", "designScore", "projectedScore", "plannerSteps", "meanTravel", "maxSaturation", "saturatedActuators", "maxLinkStrain", "meanLinkStrain", "maxReferenceDisplacement", "maxSurfaceSlope", "meanSurfaceSlope", "signedTargetError", "targetErrorRange", "sensitivityCells", "meanSensitivity", "jacobianColumns", "meanReachability", "linearFitSteps", "linearFitError", "fps", "drawCalls", "triangles"]:
            self.assertIn(f'id="{metric_id}"', html)
        for plan_id in ["inversePlanSummary", "inversePlanDelta", "inversePlanHistory", "inverseCandidateList"]:
            self.assertIn(f'id="{plan_id}"', html)
        for preview_id in ["inversePreviewSummary", "clearInversePreview", "captureInversePreview", "inverseStepPreview", "inverseStepPreviewOut"]:
            self.assertIn(f'id="{preview_id}"', html)

    def test_signed_target_error_diagnostics_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["signedTargetError", "targetErrorRange", "showTargetErrorVectors"]:
            self.assertIn(f'id="{control_id}"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn('cellVisualMode: "abstract"', state_js)
        self.assertIn('simulationMode: "kinematic"', state_js)
        self.assertIn("lockAlpha", state_js)
        self.assertIn("lockZ", state_js)
        for symbol in ["targetErrorVectorsVisible: true", "state.view.targetErrorVectorsVisible"]:
            self.assertIn(symbol, state_js)
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        for symbol in ["meanSignedTargetError", "maxPositiveTargetError", "maxNegativeTargetError"]:
            self.assertIn(symbol, math_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("this.els.cellVisualMode.addEventListener", ui_js)
        self.assertIn("this.els.simulationMode.addEventListener", ui_js)
        for symbol in [
            "applyCellVisualMode(mode)",
            '["abstract", "graph", "sheetOnly", "paperRad", "calibratedRad", "cadRad", "mechanism"].includes(mode)',
            'visualMode === "graph"',
            'visualMode === "sheetOnly"',
            "this.state.view.explodedSelected = false",
            "stopsVisible: true",
            "pivotsVisible: true",
            "linkagesVisible: true",
            "fastenersVisible: true",
            "actuatorsVisible: true",
            "measurementMode: \"all\"",
            "measurementMode: \"hardware\"",
            "measurementMode: \"backlash\"",
            "stopsVisible: false",
            "pivotsVisible: false",
            "linkagesVisible: false",
            "fastenersVisible: false",
            "actuatorsVisible: false",
            "measurementMode: \"alpha\"",
        ]:
            self.assertIn(symbol, ui_js)
        for symbol in ["showTargetErrorVectors", "signedTargetError", "targetErrorRange", "meanSignedTargetError"]:
            self.assertIn(symbol, ui_js)
        for symbol in ['this.state.view.overlayMode === "modelError"', 'this.state.view.simulationMode = "springPreview"']:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["targetErrorVectorsVisible", "renderErrorRods", "mean signed", 'mode === "modelError"', "modelErrorHeight", "modelErrorCenter", "physicalMaxHeightDelta"]:
            self.assertIn(symbol, renderer_js)

    def test_browser_programmable_discontinuity_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertLess(html.index("./math.js"), html.index("./operators.js"))
        self.assertLess(html.index("./operators.js"), html.index("./inverse.js"))
        for control_id in [
            "operatorCommitLock",
            "operatorReleaseLock",
            "operatorCheckOrder",
            "operatorOrderState",
            "operatorAlphaError",
            "operatorHeightError",
            "operatorSequenceState",
            "operatorMaxOrderError",
        ]:
            self.assertIn(f'id="{control_id}"', html)

        operators_js = (WEB / "operators.js").read_text(encoding="utf-8")
        for symbol in [
            "localActuationEvent",
            "lockEvent",
            "releaseEvent",
            "clearActuationEvent",
            "applyProgrammableEvent",
            "applyEventSequence",
            "compareEventOrder",
            "compareSequenceOrder",
            "stateDistance",
            "noncommutingAdjacentPairs",
            "maxOrderError",
            "finiteDieOffRadius",
            "lockAlphaCommutes",
            "lockZCommutes",
            "finalHeightError",
        ]:
            self.assertIn(symbol, operators_js)

        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in [
            "commitSelectedLockEvent",
            "commitDirectLockState",
            "releaseSelectedLockEvent",
            "checkSelectedEventOrder",
            "selectedEventBaseline",
            "updateOperatorInspector",
            "RAD.applyProgrammableEvent",
            "RAD.compareEventOrder",
            "RAD.compareSequenceOrder",
            "operatorSequenceState",
            "operatorMaxOrderError",
            "noncommutingAdjacentPairs",
            "operator-order-check",
        ]:
            self.assertIn(symbol, ui_js)

    def test_surface_interpolation_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["surfaceInterpolation", "surfaceSubdivisions", "showSurfaceContours", "contourMode", "animationResponse", "animationResponseOut"]:
            self.assertIn(f'id="{control_id}"', html)
        for option in ['value="membrane"', 'value="target"', 'value="error"']:
            self.assertIn(option, html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn('surfaceInterpolation: "smooth"', state_js)
        self.assertIn("surfaceSubdivisions: 5", state_js)
        self.assertIn("surfaceContoursVisible: false", state_js)
        self.assertIn('contourMode: "membrane"', state_js)
        self.assertIn("animationResponse: 0.18", state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("this.els.animationResponse.addEventListener", ui_js)

    def test_timeline_interpolation_controls_and_state_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["transitionMs", "smoothTimeline", "transitionMsOut", "sequenceSummary", "sequenceAnalysisSummary", "sequenceFrameDetail", "sequenceChart", "sequenceFileInput"]:
            self.assertIn(f'id="{control_id}"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["interpolateSnapshots", "interpolateMatrix", "exportExperimentSequence", "importExperimentSequence", "rad-sim.sequence.v1", "commandSummary", "transitionMs: 900", "smooth: true", 'keyframeName: "pose"', "Unsupported RAD sequence JSON schema"]:
            self.assertIn(symbol, state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["animateTimelineIndex", "snapshotForTimelineIndex", "playTimelineLoop", "captureKeyframe", "jumpKeyframe", "keyframeEvents", "saveSequence", "loadSequence", "analyzeSequence", "jumpAnalysisFrame", "analysisFrameForCurrentView", "selectWorstResidualCell", "max-saturation", "saveSequenceCsv", "sequenceSummaryText", "sequenceAnalysisText", "sequenceFrameDetailText", "renderSequenceChart", "sequenceChartFrameFromEvent", "hoverSequenceChart", "clearSequenceChartHover", "jumpSequenceChart", "countSnapshotCommands", "requestAnimationFrame"]:
            self.assertIn(symbol, ui_js)
        self.assertIn("this.els.timelineIndex.addEventListener", ui_js)
        self.assertIn("this.els.sequenceChart.addEventListener", ui_js)
        self.assertIn('this.els.sequenceChart.addEventListener("pointermove"', ui_js)
        self.assertIn('this.els.sequenceChart.addEventListener("pointerleave"', ui_js)
        self.assertIn("this.animateTimelineIndex(Number(this.els.timelineIndex.value)", ui_js)

    def test_selected_cell_controls_preview_live_before_event_commit(self):
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in [
            'this.els.alphaCommand.addEventListener("input", () => this.previewSelected())',
            'this.els.zCommand.addEventListener("input", () => this.previewSelected())',
            'this.els.locked.addEventListener("change", () => this.previewSelected())',
            'this.els.positionLocked.addEventListener("change", () => this.previewSelected())',
            "previewSelected()",
            "nudgeSelectedCommand(deltaAlpha, deltaZ)",
            "setSelectedCommand(alpha, z)",
            'document.getElementById("nudgeAlphaContract").addEventListener',
            'document.getElementById("nudgeZUp").addEventListener',
            "RAD.snapshotState(this.state)",
            "alpha: this.state.cells.commandAlpha[r][c]",
            "z: this.state.cells.commandZ[r][c]",
        ]:
            self.assertIn(symbol, ui_js)
        self.assertIn("s.cells.commandAlpha[r][c] ?? 0", ui_js)
        self.assertIn("s.cells.commandZ[r][c] ?? 0", ui_js)

    def test_position_lock_controls_schema_and_empirical_anchor_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["positionLocked", "positionLockRowEnds"]:
            self.assertIn(f'id="{control_id}"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in [
            "positionLocked: matrix(rows, cols, false)",
            "old.positionLocked",
            "positionLocked: interpolateMatrix",
            "state.cells.positionLocked",
            "lockZ: matrix(rows, cols, 0)",
            "state.cells.lockZ",
        ]:
            self.assertIn(symbol, state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in [
            "positionLockSelectedRowEnds",
            "this.els.positionLocked.checked",
            "is-position-locked",
            "position-lock-row-ends",
            "lockZ",
        ]:
            self.assertIn(symbol, ui_js)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        self.assertIn(".cell.is-position-locked", styles)
        surrogate_js = (WEB / "lock_surrogate.js").read_text(encoding="utf-8")
        for symbol in ["positionAnchorIndices", "empiricalLockEndpointAnchoredRows", "empiricalLockPositionAnchorCells"]:
            self.assertIn(symbol, surrogate_js)
        physics_js = (WEB / "physics.js").read_text(encoding="utf-8")
        for symbol in [
            "function positionLocked(state, r, c)",
            "function positionLockTarget(state, r, c)",
            "fixedTargets[r][c]",
            "physicalPositionLockedCells",
        ]:
            self.assertIn(symbol, physics_js)

    def test_vertical_actuation_residual_couples_to_neighbors(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn('id="zCoupling"', html)
        self.assertIn('id="zCouplingOut"', html)
        for symbol in [
            'id="pinRadius"',
            'id="holeRadius"',
            'id="pinHoleClearanceOut"',
            'id="paperSideLengthMm"',
            'id="paperHoleToleranceMm"',
            'id="backlashMmOut"',
            'id="pinHoleClearanceMmOut"',
            'id="holeToleranceModelOut"',
        ]:
            self.assertIn(symbol, html)
        for symbol in ["zCouplingGain: 0.32", "pinRadius: 0.18", "holeRadius: 0.225", "paperSideLengthMm: 35", "paperHoleToleranceMm: 0.1", "ensureGridSchema", "state.grid.zCouplingGain === undefined", "state.grid.paperSideLengthMm === undefined", "state.grid.paperHoleToleranceMm === undefined", "zCouplingGain"]:
            self.assertIn(symbol, state_js)
        for symbol in [
            "computeVerticalResidual",
            "pinHoleClearance",
            "paperRadReference",
            "paperRadCalibration",
            "modelLengthToMm",
            "mmToModelLength",
            "verticalDeadZone",
            "zResidual",
            "zDieOff",
            "function positionLockTarget(state, r, c)",
            "positionLockX",
            "positionLockZ",
            "const deadZone",
            "backlashActivation(signal, deadZone) * zCouplingGain",
            "RAD.solveBoundaryConstraints",
            "inducedHeight",
            "compressionResidual",
            "zResidual[r][c] + (constraintRealization.inducedHeight?.[r]?.[c] || 0)",
        ]:
            self.assertIn(symbol, math_js)
        self.assertIn('id="applyA360CadProfile"', html)
        for symbol in ["zCoupling: document.getElementById(\"zCoupling\")", "pinRadius: document.getElementById(\"pinRadius\")", "holeRadius: document.getElementById(\"holeRadius\")", "paperSideLengthMm: document.getElementById(\"paperSideLengthMm\")", "paperHoleToleranceMm: document.getElementById(\"paperHoleToleranceMm\")", "this.state.grid.zCouplingGain", "this.state.grid.pinRadius", "this.state.grid.holeRadius", "this.state.grid.paperSideLengthMm", "this.state.grid.paperHoleToleranceMm", "zCouplingOut", "pinHoleClearanceOut", "paperSideLengthMmOut", "paperHoleToleranceMmOut", "backlashMmOut", "pinHoleClearanceMmOut", "holeToleranceModelOut", "RAD.paperRadCalibration", "saveHardwareProfile", "loadHardwareProfile", "hardwareProfileFileInput", "applyA360CadProfile", "RAD.cadRadCellReferenceProfile", "a360-cad-profile-applied", "RAD.exportHardwareProfileJson", "RAD.importHardwareProfileJson"]:
            self.assertIn(symbol, ui_js)

    def test_selected_influence_footprint_is_visible_and_toggleable(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="showInfluenceFootprint"', html)
        self.assertIn("Show influence footprint", html)
        for symbol in ['id="couplingInspector"', 'id="alphaNeighborSignal"', 'id="zNeighborSignal"', 'id="alphaReachCells"', 'id="zReachCells"', 'id="alphaDeadZoneState"', 'id="zDeadZoneState"']:
            self.assertIn(symbol, html)
        css = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [".coupling-inspector", ".coupling-inspector dt", ".coupling-inspector dd", ".coupling-state-row", ".coupling-state.is-coupled", ".coupling-state.is-gap"]:
            self.assertIn(symbol, css)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["influenceFootprintVisible: true", "state.view.influenceFootprintVisible === undefined"]:
            self.assertIn(symbol, state_js)
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        for symbol in ["selectedCellFootprint", "selectedCellCouplingMetrics", "propagateSingleSource", "alphaDieOff", "zDieOff", "alphaNeighborSignal", "zNeighborSignal", "alphaCoupled", "zCoupled", "RAD.selectedCellFootprint", "RAD.selectedCellCouplingMetrics"]:
            self.assertIn(symbol, math_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["showInfluenceFootprint", "this.state.view.influenceFootprintVisible", "this.els.showInfluenceFootprint.checked", "updateCouplingInspector", "updateDeadZoneBadge", "RAD.selectedCellCouplingMetrics", "alphaReachCells", "zReachCells", "`${label} free gap`", "`${label} coupled`"]:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in [
            "influenceFootprintRoot",
            "footprintAlpha",
            "footprintZ",
            "renderSelectedInfluenceFootprint",
            "RAD.selectedCellFootprint",
            "addFootprintRing",
            "addFootprintLine",
            "state.view.influenceFootprintVisible === false",
        ]:
            self.assertIn(symbol, renderer_js)

    def test_actuator_allowed_mask_is_present_and_used_by_inverse_tools(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="actuatorAllowed"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["actuatorAllowed: matrix(rows, cols, true)", "old.actuatorAllowed", "state.cells.actuatorAllowed = matrix(rows, cols, true)"]:
            self.assertIn(symbol, state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["updateSelectedActuatorMask", "applyActuatorMaskPreset", "actuator-mask", "is-disallowed", "this.els.actuatorAllowed.addEventListener"]:
            self.assertIn(symbol, ui_js)
        inverse_js = (WEB / "inverse.js").read_text(encoding="utf-8")
        self.assertIn("actuatorAllowed?.[r]?.[c] === false", inverse_js)
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        self.assertIn("state.cells.actuatorAllowed?.[r]?.[c] === false", math_js)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        self.assertIn(".cell.is-disallowed", styles)

    def test_backlash_limit_stop_controls_and_geometry_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="showStops"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn("stopsVisible: false", state_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["backlashStop", "updateBacklashStops", "state.view.stopsVisible", "this.materials.stop", 'state.view.cellVisualMode || "abstract"']:
            self.assertIn(symbol, renderer_js)

    def test_pivot_and_brace_controls_and_geometry_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="showPivots"', html)
        self.assertIn('id="showLinkages"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn("pivotsVisible: false", state_js)
        self.assertIn("linkagesVisible: false", state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("this.els.showLinkages.addEventListener", ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["pivotBoss", "diagonalBrace", "pivotBosses", "updateDiagonalBraces", "state.view.pivotsVisible", "this.materials.brace"]:
            self.assertIn(symbol, renderer_js)
        for symbol in ["linkageRoot", "latticeConnectors", "createLatticeConnectors", "addLatticeConnector", "updateLatticeConnectors", "connectorEndpoint", "connectorStrain", "interCellLink", "state.view.linkagesVisible"]:
            self.assertIn(symbol, renderer_js)

    def test_linkage_strain_diagnostics_are_present(self):
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        for symbol in ["computeLinkStrain", "linkStrain", "meanAbsLinkStrain", "maxAbsLinkStrain", "horizontal", "vertical"]:
            self.assertIn(symbol, math_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["localLinkStrain", "maxLinkStrain", "meanLinkStrain", "link strain"]:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ['mode === "strain"', "localLinkStrainStrength", "connectorStrain", "overlayMaterial(\"strain\"", "link strain"]:
            self.assertIn(symbol, renderer_js)

    def test_surface_slope_diagnostics_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["maxSurfaceSlope", "meanSurfaceSlope", "maxNormalTilt", "meanNormalTilt", "showSurfaceNormals"]:
            self.assertIn(f'id="{control_id}"', html)
        self.assertIn('value="slope"', html)
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        for symbol in ["computeSurfaceSlope", "meanSurfaceSlope", "maxSurfaceSlope", "meanNormalTilt", "maxNormalTilt", "slope", "magnitude", "normal", "tilt"]:
            self.assertIn(symbol, math_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["maxSurfaceSlope", "meanSurfaceSlope", "maxNormalTilt", "meanNormalTilt", "showSurfaceNormals", "surfaceNormalsVisible", "sim.slope.magnitude", "normal tilt"]:
            self.assertIn(symbol, ui_js)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["surfaceNormalsVisible: false", "state.view.surfaceNormalsVisible"]:
            self.assertIn(symbol, state_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ['mode === "slope"', "maxSurfaceSlope", "surfaceNormalRoot", "renderSurfaceNormalVectors", "addSurfaceNormalVector", "surfaceNormalsVisible", "surfaceNormalHead", "normal tilt"]:
            self.assertIn(symbol, renderer_js)

    def test_reference_lattice_overlay_is_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="showReference"', html)
        self.assertIn('id="showDisplacementVectors"', html)
        self.assertIn('value="displacement"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["referenceVisible: true", "displacementVectorsVisible: false", "state.view.referenceVisible", "state.view.displacementVectorsVisible"]:
            self.assertIn(symbol, state_js)
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        for symbol in ["referenceCenter", "displacement", "meanReferenceDisplacement", "maxReferenceDisplacement"]:
            self.assertIn(symbol, math_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["showReference", "showDisplacementVectors", "referenceVisible", "displacementVectorsVisible", "maxReferenceDisplacement", "ref disp"]:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["referenceRoot", "displacementVectorRoot", "syncReferenceLattice", "renderDisplacementVectors", "addDisplacementVector", "this.materials.reference", "displacementHead", 'mode === "displacement"', "ref disp"]:
            self.assertIn(symbol, renderer_js)

    def test_plate_fastener_controls_and_geometry_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="showFasteners"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn("fastenersVisible: false", state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("this.els.showFasteners.addEventListener", ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["plateEdge", "screwHead", "screwSlot", "fastener", "createPlateHardware", "plateHardware", "state.view.fastenersVisible"]:
            self.assertIn(symbol, renderer_js)

    def test_abstract_cad_cell_visual_mode_is_present(self):
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        two_cell_js = (WEB / "two_cell_bench.js").read_text(encoding="utf-8")
        for symbol in [
            "abstractBody",
            "abstractInner",
            "abstractEdge",
            "abstractCorner",
            "abstractGuide",
            "abstractDatum",
            "topologyGraphRoot",
            "renderTopologyGraph",
            "topologyNode",
            "topologyDeletedEdge",
            "Topology graph",
            "../assets/cad/RADs_unit_cell_preview.png",
            "cad-reference-card",
            "cellVisualMode",
            "abstractMode",
            "paperRadMode",
            "record.abstract",
            "record.paperRad",
            "Cell abstraction",
            "Paper RAD cell",
            "paperOuter",
            "paperInner",
            "siteRadiusToPitch",
            "padSites",
            "paperClearance",
            "cadHub: this.makeCadAnnularDiskGeometry",
            "cadPad: this.makeCadAnnularDiskGeometry",
            "makeCadAnnularDiskGeometry",
            "shape.holes.push(hole)",
            "pinVisualRadius",
            "holeVisualRadius",
            "state.view.cellVisualMode || \"abstract\"",
            "state.view.cellVisualMode || \"abstract\") === \"abstract\"",
        ]:
            if symbol in {"Cell abstraction", "Paper RAD cell", "Topology graph"}:
                self.assertIn(symbol, (WEB / "index.html").read_text(encoding="utf-8"))
            elif symbol in {"../assets/cad/RADs_unit_cell_preview.png", "cad-reference-card"}:
                self.assertIn(symbol, (WEB / "index.html").read_text(encoding="utf-8") + (WEB / "styles.css").read_text(encoding="utf-8"))
            else:
                self.assertIn(symbol, renderer_js)
        for symbol in [
            "RAD.cadRadCellLayout",
            "RAD.cadRadCellReferenceProfile",
            "RAD.exportCadRadCellReferenceProfileJson",
            "twoCellConnectorContactReport",
            "connectorPitchMm",
            "siteRadiusToWidth",
            "holeToPadRadiusRatio",
            "pinToHoleRadiusRatio",
            "viewerAccess",
            "visibleFeatures",
            "a360-rads-unit-cell-visual-profile",
            "rad-sim.cad-rad-cell-reference-profile.v1",
        ]:
            self.assertIn(symbol, two_cell_js)
        for symbol in [
            "twoCellConnectorRoot",
            "renderTwoCellConnectorContacts",
            "twoCellConnectorPoint",
            "twoCellConnectorMaterial",
            "createTwoCellContactConnector",
            "twoCellConnectorVertical",
            "twoCellConnectorContact",
        ]:
            self.assertIn(symbol, renderer_js)

    def test_two_cell_physical_suite_browser_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        two_cell_js = (WEB / "two_cell_bench.js").read_text(encoding="utf-8")
        for symbol in [
            "twoCellSuiteCase",
            "twoCellExternalCase",
            "twoCellAtlasCase",
            "applyTwoCellSuiteCase",
            "applyTwoCellExternalCase",
            "applyTwoCellAtlasCase",
            "runTwoCellSuite",
            "runTwoCellFidelityMatrix",
            "runTwoCellPhaseMap",
            "runTwoCellPhaseDiagram",
            "runTwoCellTransitionReport",
            "runTwoCellResponseAtlas",
            "saveTwoCellSuiteCsv",
            "saveTwoCellFidelityMatrixCsv",
            "saveTwoCellPhaseMap",
            "saveTwoCellPhaseMapCsv",
            "saveTwoCellPhaseDiagram",
            "saveTwoCellPhaseDiagramCsv",
            "saveTwoCellTransitionReport",
            "saveTwoCellTransitionReportCsv",
            "loadTwoCellTransitionMeasurements",
            "saveTwoCellTransitionComparison",
            "saveTwoCellTransitionRerun",
            "applyTwoCellTransitionCalibration",
            "saveTwoCellResponseAtlas",
            "saveTwoCellResponseAtlasCsv",
            "saveTwoCellFidelityMatrixTemplateCsv",
            "saveTwoCellExternalFidelityManifest",
            "saveTwoCellCadContactDecomposition",
            "saveTwoCellCadContactDecompositionCsv",
            "saveTwoCellExactContactHandoffPlan",
            "saveTwoCellExactContactHandoffPlanCsv",
            "loadTwoCellFidelityMatrixMeasurements",
            "saveTwoCellFidelityMatrixComparison",
            "saveTwoCellFidelityMatrixCalibration",
            "applyTwoCellFidelityMatrixCalibration",
            "twoCellFidelityMatrixMeasurementFileInput",
            "twoCellTransitionMeasurementFileInput",
            "saveTwoCellConnectorTemplateCsv",
            "saveTwoCellPhysicalStatus",
            "twoCellSuiteSummary",
            "twoCellFidelityMatrixSummary",
            "twoCellPhaseMapSummary",
            "twoCellPhaseDiagramSummary",
            "twoCellPhaseDiagramGrid",
            "twoCellTransitionReportSummary",
            "twoCellTransitionComparisonSummary",
            "twoCellTransitionRerunSummary",
            "twoCellResponseAtlasSummary",
            "twoCellResponseAtlasCaseSummary",
            "twoCellMatrixComparisonSummary",
            "twoCellMatrixCalibrationSummary",
            "twoCellExternalCaseSummary",
            "twoCellExternalCaseDetail",
            "showExternalCase",
            "twoCellPhysicalStatusSummary",
            "twoCellArchiveSummary",
            "twoCellCadContactDecompositionSummary",
            "twoCellExactContactPlanSummary",
        ]:
            self.assertIn(symbol, html + ui_js)
        for symbol in [
            "rad-sim.two-cell-physical-simulation-suite.v1",
            "rad-sim.two-cell-physical-fidelity-matrix.v1",
            "rad-sim.two-cell-contact-phase-map.v1",
            "rad-sim.two-cell-physical-response-atlas.v1",
            "twoCellPhysicalSimulationSuite",
            "twoCellPhysicalSuiteCaseOptions",
            "twoCellPhysicalFidelityMatrix",
            "exportTwoCellPhysicalFidelityMatrixCsv",
            "twoCellContactPhaseMap",
            "exportTwoCellContactPhaseMapCsv",
            "two-cell-contact-phase-map-run",
            "twoCellRadiusBacklashPhaseDiagram",
            "exportTwoCellRadiusBacklashPhaseDiagramCsv",
            "two-cell-radius-backlash-phase-diagram-run",
            "twoCellPhaseDiagramCasePreview",
            "renderTwoCellPhaseDiagramGrid",
            "applyTwoCellPhaseDiagramGridCell",
            "phaseDiagramClassName",
            "phaseDiagramLabel",
            "two-cell-radius-backlash-phase-diagram-case-applied",
            "twoCellRadiusBacklashTransitionReport",
            "exportTwoCellRadiusBacklashTransitionReportCsv",
            "rad-sim.two-cell-radius-backlash-transition-report.v1",
            "two-cell-radius-backlash-transition-report-run",
            "compareTwoCellRadiusBacklashTransitionMeasurements",
            "exportTwoCellRadiusBacklashTransitionComparisonCsv",
            "rad-sim.two-cell-radius-backlash-transition-comparison.v1",
            "twoCellRadiusBacklashTransitionRerun",
            "exportTwoCellRadiusBacklashTransitionRerunCsv",
            "rad-sim.two-cell-radius-backlash-transition-rerun.v1",
            "two-cell-radius-backlash-transition-measurements-load",
            "proposedReducedProxyUpdates",
            "two-cell-radius-backlash-transition-calibration-applied",
            "twoCellPhysicalResponseAtlas",
            "exportTwoCellPhysicalResponseAtlasCsv",
            "two-cell-response-atlas-run",
            "two-cell-response-atlas-case-applied",
            "twoCellFidelityMatrixMeasurementTemplate",
            "exportTwoCellFidelityMatrixMeasurementTemplateCsv",
            "rad-sim.two-cell-fidelity-matrix-measurement-template.v1",
            "twoCellExternalFidelityMatrixManifest",
            "exportTwoCellExternalFidelityMatrixManifestJson",
            "exportTwoCellExternalFidelityMatrixManifestCsv",
            "rad-sim.browser-two-cell-external-fidelity-matrix-manifest.v1",
            "twoCellCadContactDecompositionSpec",
            "exportTwoCellCadContactDecompositionJson",
            "exportTwoCellCadContactDecompositionCsv",
            "rad-sim.two-cell-cad-contact-decomposition.v1",
            "twoCellExactContactHandoffPlan",
            "exportTwoCellExactContactHandoffPlanJson",
            "exportTwoCellExactContactHandoffPlanCsv",
            "rad-sim.two-cell-exact-contact-handoff-plan.v1",
            "twoCellFidelityMatrixMeasurementsFromCsv",
            "compareTwoCellFidelityMatrixMeasurements",
            "exportTwoCellFidelityMatrixMeasurementComparisonJson",
            "rad-sim.browser-two-cell-fidelity-matrix-measurement-comparison.v1",
            "calibrateTwoCellFidelityMatrixParameters",
            "exportTwoCellFidelityMatrixParameterCalibrationJson",
            "rad-sim.browser-two-cell-fidelity-matrix-parameter-calibration.v1",
            "two-cell-fidelity-matrix-calibration-applied",
            "two-cell-fidelity-matrix-measurements-load",
            "observedScalarCount",
            "cadRadCellArchiveAudit",
            "twoCellPhysicalFidelityStatus",
            "exportTwoCellPhysicalFidelityStatusJson",
            "rad-sim.cad-rad-cell-archive-audit.v1",
            "rad-sim.browser-two-cell-physical-fidelity-status.v1",
            "segmentedBodyMeshExport",
            "twoCellConnectorMeasurementTemplate",
            "exportTwoCellConnectorMeasurementTemplateCsv",
            "connectorRows",
            "highBacklashReducesAlphaResponse",
            "state-lock-alpha-only",
        ]:
            self.assertIn(symbol, two_cell_js + ui_js)
        css = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [
            "phase-diagram-grid",
            "phase-diagram-cell",
            "phase-axial-and-vertical-contact",
            "phase-free-play",
        ]:
            self.assertIn(symbol, html + ui_js + css)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in [
            "externalCaseVisible: true",
            "state.view.externalCaseVisible",
            "renderExternalTwoCellCasePreview",
            "addExternalCaseVerticalMarker",
            "externalCaseCorrected",
            "externalCaseObserved",
            "externalCaseResidual",
            "externalCaseBar",
        ]:
            self.assertIn(symbol, state_js + renderer_js + ui_js)

    def test_timeline_event_formatting_covers_optimizer(self):
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("formatEvent(event)", ui_js)
        self.assertIn("target-fit-optimized", ui_js)
        self.assertIn("inverse-plan-analyzed", ui_js)
        self.assertIn("inverse-plan-applied", ui_js)
        self.assertIn("inverse-preview-keyframe", ui_js)
        self.assertIn("keyframe:", ui_js)
        self.assertIn("score", ui_js)

    def test_sequence_export_import_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["saveSequence", "loadSequence", "analyzeSequence", "jumpBestFrame", "jumpMaxSaturationFrame", "selectWorstResidualCell", "saveSequenceCsv", "sequenceSummary", "sequenceAnalysisSummary", "sequenceFrameDetail", "sequenceChart", "sequenceFileInput", "eventList"]:
            self.assertIn(f'id="{control_id}"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["exportExperimentSequence", "importExperimentSequence", "rad-sim.sequence.v1", "frameCount", "finalCommandSummary", "initialSnapshot", "frames", "RAD sequence has no frames to import"]:
            self.assertIn(symbol, state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["RAD.exportExperimentSequence", "RAD.importExperimentSequence", "RAD.exportSequenceMetricsCsv", "rad-sim-sequence.json", "rad-sim-sequence-metrics.csv", "saveSequence()", "loadSequence()", "jumpAnalysisFrame", "selectWorstResidualCell", "sequenceSummaryText", "sequenceAnalysisText", "No sequence frames"]:
            self.assertIn(symbol, ui_js)
        self.assertIn("this.stopTimelinePlayback();\n      this.state.selection", ui_js)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [".sequence-chart", "cursor: crosshair", ".sequence-legend", ".legend-rms", ".legend-saturation", ".legend-delta"]:
            self.assertIn(symbol, styles)

    def test_browser_local_response_characterization_is_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in [
            "characterizationScope",
            "runCharacterization",
            "characterizationSummary",
            "characterizationDetail",
            "characterizationZSign",
            "characterizationZExtrema",
            "characterizationSuperposition",
            "characterizationScale",
            "characterizationRank",
            "characterizationUnderactuated",
            "frameworkLawSummary",
            "frameworkLawDetail",
            "formalizationSummary",
            "formalizationDetail",
            "saveTopologyReport",
            "saveProgrammableReport",
            "saveFormalizationTargets",
            "saveResponseAtlas",
            "runResponseAtlasSweep",
            "saveResponseAtlasSweep",
            "sweepSummary",
            "sweepTrend",
        ]:
            self.assertIn(f'id="{control_id}"', html)
        for option in ['value="single"', 'value="pair"', 'value="cluster"', 'value="lattice"']:
            self.assertIn(option, html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in [
            'characterizationScope: "single"',
            "characterization: null",
            "responseAtlasSweep: null",
            "frameworkLawCandidates: null",
            "frameworkFormalizationTargets: null",
        ]:
            self.assertIn(symbol, state_js)
        analysis_js = (WEB / "analysis.js").read_text(encoding="utf-8")
        for symbol in [
            "function characterizeLocalResponse",
            "characterizationCells",
            "scopedState",
            "superpositionError",
            "responseMatrixDiagnostic",
            "matrixRankFromColumns",
            "countReachableFromColumns",
            "activeSources",
            "responseCells",
            "alphaReachCells",
            "zReachCells",
            "positiveZReachCells",
            "negativeZReachCells",
            "maxPositiveHeightDelta",
            "maxNegativeHeightDelta",
            "responseRankAlpha",
            "responseRankHeight",
            "alphaUnderactuatedCells",
            "heightUnderactuatedCells",
            "pinHoleClearanceMm",
            "programmableDiscontinuityReport",
            "exportProgrammableDiscontinuityReport",
            "formalizationTargetManifest",
            "exportFormalizationTargetManifest",
            "frameworkOperatorLawCandidates",
            "frameworkFormalizationTargets",
            "formalizationTarget",
            "responseAtlas",
            "exportResponseAtlas",
            "responseAtlasSweep",
            "exportResponseAtlasSweep",
            "calibrationBenchNotebook",
            "exportCalibrationBenchNotebook",
            "exportCalibrationBenchNotebookCsv",
            "calibrationBenchPacket",
            "exportCalibrationBenchPacket",
            "sweepSensitivityFromTrends",
            "sweepOperatorLawCandidates",
            "endpoint finite difference over each parameter trend",
            "adjacent monotonicity over sampled parameter trend plus endpoint sensitivity",
            "operatorLawCandidates",
            "rad-sim.response-atlas.v1",
            "rad-sim.response-atlas-sweep.v1",
            "rad-sim.framework-law-candidates.v1",
            "rad-sim.formalization-targets.v1",
            "rad-sim.programmable-discontinuity-report.v1",
            "rad-sim.browser-physical-preview.v1",
            "rad-sim.vertical-load-physical-preview-report.v1",
            "rad-sim.vertical-removal-physical-comparison.v1",
            "rad-sim.mechanics-energy-certificate.v1",
            "rad-sim.vertical-load-energy-validation.v1",
            "rad-sim.vertical-load-energy-measurement-results.v1",
            "rad-sim.vertical-load-energy-experiment-protocol.v1",
            "rad-sim.vertical-load-energy-comparison-report.v1",
            "rad-sim.vertical-load-bench-packet.v1",
            "springHingePhysicalPreviewReportDescriptor",
            "build_vertical_load_physical_preview_report",
            "mechanics_energy_certificate_to_dict",
            "validate_vertical_load_energy_measurements",
            "vertical_load_energy_measurement_template",
            "vertical_load_energy_experiment_protocol",
            "vertical_load_bench_packet",
            "vertical_load_energy_measurement_results_from_json",
            "compare_vertical_load_energy_measurement_results",
            "export_vertical_load_energy_experiment_protocol_json",
            "export_vertical_load_bench_packet_json",
            "export_vertical_load_energy_measurement_template_json",
            "export_vertical_load_energy_validation_json",
            "export_vertical_load_energy_comparison_report_json",
            "export_vertical_load_physical_preview_report_csv",
            "dead_zone_zero_inside_backlash",
            "group_operator_support_decomposition",
            "LocalModeOperator.groupOperators_with_disjoint_lists_commute",
            "removed_cell_clears_group_supported_constraints",
            "CellGraph.removed_cell_constraints_clear_group_operator",
            "vertical_clearance_gates_residual_contact",
            "verticalResidualStepNat_zero_inside",
            "contactPenaltyFromClearanceNat_nonnegative",
            "fixed_cell_load_work_proxy_zero",
            "loadWorkMagnitudeNat_zero_fixed",
            "signed_vertical_load_work_residual_int",
            "signedLoadWorkInt_zero_displacement",
            "signedEnergyResidualTripleInt_zero_when_equal",
            "integer_scaled_mechanics_scaffold",
            "ScaledNatQuantity",
            "scaledSpringEnergyNat_preserves_denominator",
            "scaledSignedEnergyResidualTripleInt_zero_when_equal",
            "measurement_unit_scale_invariants",
            "MeasurementUnitScale",
            "physical_unit_scale_metadata",
            "rad-sim.physical-unit-scale-metadata.v1",
            "rad-sim.hardware-profile.v1",
            "hardware_profile_from_json",
            "measurementUnitScaleNat_preserves_denominator",
            "measurementUnitScaleResidualInt_zero_when_equal",
            "hardwareProfileCoverageMissing_zero_when_complete",
            "calibration_parameter_estimate_residual_bookkeeping",
            "rad-sim.calibration-parameter-estimates.v1",
            "calibrationFitResidualPassNat_zero",
            "calibration_model_profile_safe_update_bounds",
            "rad-sim.calibration-model-profile.v1",
            "calibrationModelProfileUpdateSafeNat_intro",
            "calibration_model_profile_selection_predicate",
            "rad-sim.calibration-model-profile-residual-comparison.v1",
            "rad-sim.calibration-model-profile-selection.v1",
            "calibrationModelProfileCandidateScoreNat_le_before",
            "calibration_model_profile_holdout_predicate",
            "rad-sim.calibration-model-profile-holdout-validation.v1",
            "calibrationModelProfileHoldoutScoreNat_le_before",
            "calibration_train_holdout_split_metadata",
            "rad-sim.calibration-train-holdout-split.v1",
            "calibrationTrainHoldoutSplitReadyNat_zero_overlap",
            "calibration_train_holdout_file_provenance",
            "rad-sim.calibration-dataset-provenance.v1",
            "calibrationTrainHoldoutProvenanceReadyNat_profile_matches",
            "calibration_bench_protocol_coverage",
            "rad-sim.calibration-bench-notebook.v1",
            "calibrationBenchProtocolCoverageReadyNat_has_two_dataset_roles",
            "export_calibration_bench_notebook_csv",
            "calibration_bench_packet_completeness",
            "rad-sim.calibration-bench-packet.v1",
            "calibrationBenchPacketCompleteNat_has_manifest",
            "write_calibration_bench_packet_artifacts",
            "calibration_bench_executed_validation_gate",
            "rad-sim.calibration-bench-execution-validation.v1",
            "calibrationBenchExecutedValidationReadyNat_zero_missing_evidence",
            "write_calibration_bench_execution_validation_artifacts",
            "mechanics_energy_certificate_nonnegative_proxy",
            "mechanicalStoredEnergyNat_nonnegative",
            "measured_vertical_load_work_validation_zero_residual",
            "measuredWorkResidualNat_zero_when_equal",
            "vertical_load_comparison_pass_predicate",
            "verticalLoadScenarioPassNat_zero_errors",
            "verticalLoadScenarioPassNat_false_when_missing",
            "physical_validation_readiness_gate",
            "rad-sim.physical-validation-readiness.v1",
            "physicalValidationReadyNat_zero_missing_evidence",
            "physical_validation_readiness_report",
            "contact_state_abstraction_gate",
            "rad-sim.contact-state-abstraction.v1",
            "contactStateAbstractionReadyNat_zero_missing_evidence",
            "contact_state_abstraction_report",
            "rad_sim.compare_vertical_load_measurements",
            "spring_hinge_removed_topology_load_comparison",
            "compare_vertical_residual_spring_hinge_3d",
            "browser_spring_preview_removed_edge_deletion",
            "physicalSkippedSpringEdges",
            "CellGraph.removal_deletes_one_step_to_removed",
            "noncommutativity_witness_from_order_error",
            "browser-cannot-inspect-path",
            "thresholded diagnostic predicates over locality, reachability, composition, and event-order metrics",
            "composition_nonadditivity",
            "event_order_noncommutativity",
            "paperSupportedAssumptions",
            "simulatorDiagnostics",
            "reportPhysicalPreview",
            "RAD.characterizeLocalResponse",
        ]:
            self.assertIn(symbol, analysis_js)
        operators_js = (WEB / "operators.js").read_text(encoding="utf-8")
        for symbol in [
            "compareVerticalResidualUnderRemoval",
            "verticalContactPenalty",
            "heightLoadContactMetrics",
            "loadWorkMagnitudeDelta",
            "topologyBlockedCells",
            "contactPenaltyDelta",
        ]:
            self.assertIn(symbol, operators_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in [
            "characterizationScope: document.getElementById(\"characterizationScope\")",
            "runCharacterization",
            "renderCharacterization",
            "this.state.experiment.characterizationScope",
            "this.state.experiment.characterization = result",
            "this.state.experiment.frameworkLawCandidates",
            "this.state.experiment.frameworkFormalizationTargets",
            'type: "characterization"',
            "characterizationSummary",
            "characterizationZSign",
            "characterizationZExtrema",
            "characterizationSuperposition",
            "characterizationRank",
            "characterizationUnderactuated",
            "saveTopologyReport",
            "saveTopologyReport()",
            "saveProgrammableReport",
            "saveFormalizationTargets",
            "saveResponseAtlas",
            "runResponseAtlasSweep",
            "saveResponseAtlasSweep",
            "saveBenchPacket",
            "saveBenchNotebook",
            "saveBenchNotebookCsv",
            "saveBenchExecutionValidation",
            "saveBenchExecutionValidationCsv",
            "this.state.experiment.responseAtlasSweep = sweep",
            "renderResponseAtlasSweep",
            "renderFrameworkLawCandidates",
            "renderFormalizationTargets",
            "currentFrameworkReportMetadata",
            "sweepSummary",
            "sweepTrend",
            "sensitivity?.dominant",
            "operatorLawCandidates?.laws",
            "law ",
            "rad-sim-response-atlas",
            "rad-sim-calibration-bench-packet",
            "rad-sim-calibration-bench-notebook",
            "rad-sim-calibration-bench-execution-validation",
            "rad-sim-topology-experiment-report",
            "rad-sim-response-atlas-sweep",
            "rad-sim-programmable-discontinuity-report",
            "rad-sim-formalization-targets",
            "characterize:",
        ]:
            self.assertIn(symbol, ui_js)
        script = (ROOT / "tests" / "validate_web_modules.js").read_text(encoding="utf-8")
        for symbol in [
            "singleCharacterization",
            "pairCharacterization",
            "pair characterization should include at least the single-cell response footprint",
            "pair characterization should report finite superposition residual",
            "characterization should carry paper-scale clearance",
            "pair characterization should report alpha response rank",
            "pair characterization should report height response rank",
            "pair characterization should report alpha underactuation",
            "characterization should count upward z response cells",
            "characterization should count downward z response cells",
            "frameworkReport",
            "framework report should count active command operators",
            "framework report should include a composition law candidate",
            "framework report should include an order law candidate",
            "framework report should expose a noncommutativity formalization witness",
            "analysis module should export formalization target manifests",
            "framework report should include browser physical preview validation",
            "browserAtlas",
            "response atlas should include residual z observation cells",
            "browserSweep",
            "response atlas sweep should show clearance-gated neighbor residual trend",
            "response atlas sweep should report negative clearance sensitivity",
            "response atlas sweep should propose a clearance residual law",
            "RAD.exportProgrammableDiscontinuityReport",
        ]:
            self.assertIn(symbol, script)

    def test_browser_obj_export_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="saveObj"', html)
        self.assertLess(html.index("./analysis.js"), html.index("./physics.js"))
        self.assertLess(html.index("./physics.js"), html.index("./two_cell_bench.js"))
        self.assertLess(html.index("./two_cell_bench.js"), html.index("./mesh_export.js"))
        self.assertLess(html.index("./mesh_export.js"), html.index("./renderer.js"))
        mesh_js = (WEB / "mesh_export.js").read_text(encoding="utf-8")
        for symbol in [
            "buildPaperRadMesh",
            "exportPaperRadMeshObj",
            "outer_plate",
            "inner_plate",
            "connector",
            "pinSegments",
            "RAD.buildPaperRadMesh",
            "RAD.exportPaperRadMeshObj",
        ]:
            self.assertIn(symbol, mesh_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["saveObj", "rad-sim-paper-rad.obj", "RAD.buildPaperRadMesh", "RAD.simulateActive", "RAD.exportPaperRadMeshObj"]:
            self.assertIn(symbol, ui_js)

    def test_browser_physical_preview_mode_is_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="simulationMode"', html)
        self.assertIn('value="springPreview"', html)
        self.assertLess(html.index("./analysis.js"), html.index("./physics.js"))
        self.assertLess(html.index("./physics.js"), html.index("./two_cell_bench.js"))
        self.assertLess(html.index("./two_cell_bench.js"), html.index("./mesh_export.js"))

        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn('simulationMode: "kinematic"', state_js)

        physics_js = (WEB / "physics.js").read_text(encoding="utf-8")
        for symbol in [
            "simulatePhysicalRelaxation",
            "simulateActive",
            "spring-preview",
            "physicalRmsHeightDelta",
            "physicalActiveSpringEdges",
            "physicalSkippedSpringEdges",
            "cellRemoved(state, r, c)",
            "recomputeLinkStrain",
            "recomputeSlope",
            "RAD.simulatePhysicalRelaxation",
            "RAD.simulateActive",
        ]:
            self.assertIn(symbol, physics_js)

        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in ["RAD.simulateActive", "Model: spring preview"]:
            self.assertIn(symbol, app_js)

    def test_sequence_analysis_module_exports_csv_metrics(self):
        analysis_js = (WEB / "analysis.js").read_text(encoding="utf-8")
        for symbol in [
            "function analyzeExperimentSequence",
            "function exportSequenceMetricsCsv",
            "frameCommandStats",
            "frameDelta",
            "worstTargetResidual",
            "worstResidual",
            "rmsTargetError",
            "command_alpha",
            "theta_deg",
            "target_error",
            "normal_tilt_deg",
            "frame_worst_row",
            "frame_worst_col",
            "frame_worst_residual",
        ]:
            self.assertIn(symbol, analysis_js)

    def test_inverse_design_module_is_loaded_and_serialized(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertLess(html.index("./math.js"), html.index("./inverse.js"))
        self.assertLess(html.index("./inverse.js"), html.index("./renderer.js"))
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["inverse:", "maxActuators", "lastScore", "plan", "preview: null", "jacobian:", "linearSolution:"]:
            self.assertIn(symbol, state_js)
        self.assertIn("cloneData(value)", state_js)
        self.assertIn("value === undefined", state_js)
        inverse_js = (WEB / "inverse.js").read_text(encoding="utf-8")
        for symbol in ["buildInverseDesignPlan", "analyzeActuatorSensitivity", "buildResponseJacobian", "solveLinearizedTargetFit", "applyLinearizedTargetFit", "applyInverseDesignPlan", "setInversePreview", "setInversePlanStepPreview", "previewFromCommands", "commandFromResidual", "objective"]:
            self.assertIn(symbol, inverse_js)
        for symbol in ["projectedScore", "projectedError", "totalImprovement", "iterative-greedy", "history", "planningStateFrom", "contribution", "errorReduction"]:
            self.assertIn(symbol, inverse_js)
        self.assertIn("state.cells.locked = locked", inverse_js)

    def test_inverse_response_jacobian_is_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["buildJacobian", "selectUnderactuatedTarget", "solveLinearFit", "applyLinearFit", "saveInverseReport", "buildPreviewPacket", "runPreviewReplay", "runPreviewPhysical", "savePreviewPacket", "savePreviewReplay", "savePreviewPhysical", "jacobianColumns", "meanReachability", "underTargetCells", "linearFitSteps", "linearFitError", "inverseReachabilitySummary", "inversePacketSummary", "inversePacketReplaySummary", "inversePacketPhysicalSummary"]:
            self.assertIn(f'id="{control_id}"', html)
        self.assertIn('value="reachability"', html)
        self.assertIn('value="underactuated"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["jacobian:", "linearSolution:", "coverageMap", "columnCount", "meanCoverage", "conditionEstimate", "targetReachability", "predictedError", "projectedError"]:
            self.assertIn(symbol, state_js)
        inverse_js = (WEB / "inverse.js").read_text(encoding="utf-8")
        for symbol in ["buildResponseJacobian", "solveLinearizedTargetFit", "applyLinearizedTargetFit", "inverseDesignReport", "exportInverseDesignReport", "rad-sim.inverse-design-report.v1", "linearized-jacobian-greedy-fit", "finite-difference-response-jacobian", "finite-response-height-reachability", "sign-compatible finite-response-height-reachability", "targetReachabilityReport", "positiveHeightReachableMap", "negativeHeightReachableMap", "positiveUnderactuatedHeightCells", "negativeUnderactuatedHeightCells", "upwardTargetHeightCells", "downwardTargetHeightCells", "underactuatedHeightCells", "worstUnderactuatedCell", "flattenDelta", "targetResidualVector", "heightDelta", "alphaDelta", "targetAlignment", "conditionEstimate", "commandZ", "commandAlpha", "z+", "z-", "alpha-", "alpha+", "RAD.buildResponseJacobian"]:
            self.assertIn(symbol, inverse_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["buildJacobian", "selectUnderactuatedTarget", "solveLinearFit", "applyLinearFit", "saveInverseReport", "buildProfileInversePreviewPacket", "runProfileInversePreviewReplay", "runProfileInversePreviewPhysical", "saveProfileInversePreviewPacket", "saveProfileInversePreviewReplay", "saveProfileInversePreviewPhysical", "profile-inverse-preview-packet-built", "profile-inverse-preview-replayed", "profile-inverse-preview-physical-checked", "jacobian-built", "linear-fit-solved", "linear-fit-applied", "jacobianColumns", "meanReachability", "underTargetCells", "linearFitSteps", "linearFitError", "positiveUnderactuatedHeightCells", "negativeUnderactuatedHeightCells", '"underactuated" : "reachability"']:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ['mode === "reachability"', 'mode === "underactuated"', 'mode === "topologyBlocked"', "reachabilityStrength", "underactuatedStrength", "topologyBlockedStrength", "state.inverse?.jacobian", "linearSolution", "coverageMap", "maxCoverage", "underactuatedHeightMap"]:
            self.assertIn(symbol, renderer_js)

    def test_inverse_sensitivity_analysis_is_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in ["analyzeSensitivity", "sensitivityCells", "meanSensitivity"]:
            self.assertIn(f'id="{control_id}"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["sensitivity:", "stepZ: 0.12", "stepAlpha: 0.12", "controllableCells", "meanGain", "maxGain"]:
            self.assertIn(symbol, state_js)
        inverse_js = (WEB / "inverse.js").read_text(encoding="utf-8")
        for symbol in ["analyzeActuatorSensitivity", "finite-difference-command-response", "rmsDelta", "countResponsiveCells", "combinedGain", "dilationGain", "RAD.analyzeActuatorSensitivity"]:
            self.assertIn(symbol, inverse_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["sensitivity-analyzed", 'this.state.view.overlayMode = "sensitivity"', "sensitivityCandidates", "meanSensitivity"]:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["sensitivityStrength", "state.inverse?.sensitivity", 'mode === "sensitivity"']:
            self.assertIn(symbol, renderer_js)

    def test_configurable_actuator_limits_and_saturation_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for control_id in [
            "zTravelLimit",
            "zTravelLimitOut",
            "alphaContractLimit",
            "alphaContractLimitOut",
            "alphaExpandLimit",
            "alphaExpandLimitOut",
            "maxSaturation",
            "saturatedActuators",
        ]:
            self.assertIn(f'id="{control_id}"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["zTravelLimit: 0.8", "alphaContractLimit: 0.55", "alphaExpandLimit: 0.35"]:
            self.assertIn(symbol, state_js)
        math_js = (WEB / "math.js").read_text(encoding="utf-8")
        for symbol in ["commandLimits", "clampCommandZ", "clampCommandAlpha", "commandSaturation", "clampAllCommands", "saturatedActuators", "maxSaturation"]:
            self.assertIn(symbol, math_js)
        inverse_js = (WEB / "inverse.js").read_text(encoding="utf-8")
        for symbol in ["RAD.clampCommandZ", "RAD.clampCommandAlpha", "RAD.commandLimits"]:
            self.assertIn(symbol, inverse_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["updateTravelLimits", "RAD.clampAllCommands", "maxSaturation", "saturatedActuators"]:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ['mode === "saturation"', "sim.saturation", "RAD.commandSaturation", "normalizedZ", "limits.z"]:
            self.assertIn(symbol, renderer_js)

    def test_ui_renders_inverse_plan_table_and_rich_measurements(self):
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["renderInversePlan", "previewInverseCandidate", "captureInversePreviewKeyframe", "formatInversePreview", "isPreviewedCandidate", "candidate-row", "history-row", "inversePlanSummary", "inversePreviewSummary", "residual"]:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["cmd a", "cmd z", "residual", "influence"]:
            self.assertIn(symbol, renderer_js)
        for symbol in ['mode === "zresidual"', "sim.zResidual", "z residual", "measurementLabelPosition", "labelLeader", "state.view.measurementLabelsVisible === true"]:
            self.assertIn(symbol, renderer_js)
        for symbol in ["halfX + margin", "halfY + margin * 0.78", "Math.max(center.z + 0.84, 0.86)"]:
            self.assertIn(symbol, renderer_js)
        self.assertIn("z residual", ui_js)

    def test_measurement_modes_are_present_end_to_end(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="measurementMode"', html)
        self.assertIn('id="showMeasurementLabels"', html)
        for option in ['value="all"', 'value="alpha"', 'value="theta"', 'value="height"', 'value="backlash"', 'value="hardware"']:
            self.assertIn(option, html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn('measurementMode: "alpha"', state_js)
        self.assertIn("measurementLabelsVisible: false", state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("this.els.measurementMode.addEventListener", ui_js)
        self.assertIn("showMeasurementLabels", ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["measurementLabel", "measurementRingRadius", "renderThetaGuide", "renderBacklashGuide", "gapLine", 'mode === "hardware"', "calibrationProfileSummary", "pinHoleClearanceMm"]:
            self.assertIn(symbol, renderer_js)

    def test_renderer_uses_persistent_surface_meshes(self):
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["syncSurfaceMeshes", "updateSurfaceMesh", "this.membraneMesh", "this.targetMesh", "surfaceContourRoot"]:
            self.assertIn(symbol, renderer_js)
        self.assertIn("position.needsUpdate = true", renderer_js)
        for symbol in ["renderSurfaceContours", "surfaceContourPoint", "addSurfaceContourLine", "targetContour"]:
            self.assertIn(symbol, renderer_js)
        for symbol in ["rows <= 1 ? 0", "cols <= 1 ? 0"]:
            self.assertIn(symbol, renderer_js)

    def test_renderer_exposes_diagnostics_and_reuses_overlay_materials(self):
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["getDiagnostics", "updateDiagnostics", "this.overlayMaterials", "overlayMaterial(mode, value)"]:
            self.assertIn(symbol, renderer_js)
        self.assertNotIn("this.materials.plate.clone();\n        let t = 0.5", renderer_js)

    def test_inverse_overlay_includes_accepted_commands(self):
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        self.assertIn("...(plan.commands || [])", renderer_js)
        self.assertIn("...(plan.candidates || [])", renderer_js)
        self.assertIn("inversePreviewContributionStrength", renderer_js)
        self.assertIn("preview?.contribution", renderer_js)

    def test_renderer_has_smooth_membrane_sampling_and_vertical_actuators(self):
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["catmullRom", "sampleScalarField", "surfaceSubdivisions", "surfaceInterpolation"]:
            self.assertIn(symbol, renderer_js)
        self.assertNotIn("grid.map((row) => row.map((p) => p.x))", renderer_js)
        self.assertIn("rail.rotation.x = Math.PI / 2", renderer_js)
        self.assertIn("sliderZ", renderer_js)
        self.assertIn("animationResponse", renderer_js)

    def test_actuator_display_modes_and_travel_gauge_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        for option in ['value="selected"', 'value="planned"', 'value="active"']:
            self.assertIn(option, html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn('actuatorDisplayMode: "selected"', state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("this.els.actuatorDisplayMode.addEventListener", ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["shouldShowActuator", "isPlannedActuator", "actuatorTravelGauge", "travelNeedle", "actuatorGaugeTick", "zTravelTicks", "updateZTravelGaugeTicks", "actuatorPositive", "actuatorNegative"]:
            self.assertIn(symbol, renderer_js)

    def test_in_plane_alpha_actuator_hardware_is_present(self):
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in [
            "alphaActuatorRail",
            "alphaActuatorSleeve",
            "alphaActuatorStop",
            "alphaRail",
            "alphaSleeve",
            "alphaNeedle",
            "alphaMinStop",
            "alphaMaxStop",
            "normalizedAlpha",
            "alphaTravel",
            "alpha travel",
            "alphaTravelTicks",
            "updateAlphaTravelGaugeTicks",
        ]:
            self.assertIn(symbol, renderer_js)

    def test_selected_measurement_travel_guides_are_present(self):
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in [
            "renderActuatorTravelGuide",
            "addMeasurementLine",
            "commandZ",
            "commandAlpha",
            "zActual",
            "alphaActual",
            "limits.alphaContract",
            "limits.alphaExpand",
        ]:
            self.assertIn(symbol, renderer_js)


if __name__ == "__main__":
    unittest.main()
