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
            "RAD.compareEventOrder",
            "RAD.paperRadCalibration",
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
            "browser entry should not require external scripts or styles",
            "scriptOrder",
            "state.grid.backlash",
            "state.grid.couplingGain",
            "state.grid.zCouplingGain",
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
        self.assertLess(html.index("./operators.js"), html.index("./inverse.js"))
        self.assertLess(html.index("./math.js"), html.index("./analysis.js"))
        self.assertLess(html.index("./analysis.js"), html.index("./physics.js"))
        self.assertLess(html.index("./physics.js"), html.index("./mesh_export.js"))
        self.assertLess(html.index("./mesh_export.js"), html.index("./renderer.js"))

    def test_browser_modules_expose_expected_api(self):
        expected = {
            "state.js": ["RAD.createState", "RAD.serialize", "RAD.updateDerivedCells", "RAD.exportExperimentSequence", "RAD.importExperimentSequence", "RAD.deserialize"],
            "operators.js": ["RAD.localActuationEvent", "RAD.lockEvent", "RAD.applyEventSequence", "RAD.compareEventOrder", "RAD.finiteDieOffRadius"],
            "analysis.js": ["RAD.analyzeExperimentSequence", "RAD.exportSequenceMetricsCsv", "RAD.sequenceFrames"],
            "physics.js": ["RAD.simulatePhysicalRelaxation", "RAD.simulateActive"],
            "mesh_export.js": ["RAD.buildPaperRadMesh", "RAD.exportPaperRadMeshObj"],
            "math.js": [
                "RAD.backlashActivation",
                "RAD.alphaToTheta",
                "RAD.compileCustomTargetExpression",
                "RAD.simulate",
                "RAD.applyPreset",
                "RAD.optimizeCommandsForTarget",
                "RAD.paperRadCalibration",
                "RAD.modelLengthToMm",
                "RAD.mmToModelLength",
            ],
            "inverse.js": ["RAD.buildInverseDesignPlan", "RAD.applyInverseDesignPlan"],
            "renderer.js": ["RAD.RadRenderer"],
            "ui.js": ["RAD.RadUI"],
            "app.js": ["RAD_APP"],
        }
        for filename, names in expected.items():
            text = (WEB / filename).read_text(encoding="utf-8")
            for name in names:
                self.assertIn(name, text, f"{filename} should expose {name}")

    def test_app_exposes_browser_validation_handles(self):
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in ["get state()", "get simulation()", "get renderer()", "update: (nextState = state) => renderAll(nextState)"]:
            self.assertIn(symbol, app_js)

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

    def test_desktop_layout_keeps_lattice_visible_while_controls_scroll(self):
        css = (WEB / "styles.css").read_text(encoding="utf-8")
        for symbol in [
            "body {\n  margin: 0;\n  min-height: 100vh;\n  overflow: hidden;",
            ".app-shell {\n  display: grid;\n  grid-template-columns: minmax(0, 1fr) 350px;\n  height: 100vh;",
            ".workspace {\n  min-width: 0;\n  height: 100vh;",
            "display: flex;\n  flex-direction: column;\n  overflow: hidden;",
            ".viewport {\n  position: relative;\n  --viewport-toolbar-height: 64px;\n  flex: 1 1 auto;",
            ".three-mount {\n  width: 100%;\n  height: 100%;",
            ".control-panel {\n  height: 100vh;",
            "overflow-y: auto;\n  overscroll-behavior: contain;",
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
        for mode in ["inspect", "edit", "analyze"]:
            self.assertIn(f'data-workspace-mode="{mode}"', html)
        app_js = (WEB / "app.js").read_text(encoding="utf-8")
        for symbol in [
            "workspaceMode",
            "setActiveWorkspaceMode",
            "applyWorkspaceMode",
            "modeInspect",
            "modeEdit",
            "modeAnalyze",
            "Workspace:",
            "targetAvailable",
            "actuatorDisplayMode: \"planned\"",
            "overlayMode: targetAvailable ? \"error\" : \"height\"",
        ]:
            self.assertIn(symbol, app_js)
        styles = (WEB / "styles.css").read_text(encoding="utf-8")
        self.assertIn("[data-workspace-mode]", styles)

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
        for symbol in ["overlayLegendText", "updateOverlayLegend", "Target Error", "Reachability", "Cell State", "Z Residual"]:
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
        self.assertIn('value="paperRad"', html)
        self.assertIn('value="mechanism"', html)
        self.assertIn('value="springPreview"', html)
        for option in ['value="zresidual"', 'value="error"', 'value="travel"', 'value="saturation"', 'value="strain"', 'value="displacement"', 'value="slope"', 'value="inverse"', 'value="sensitivity"', 'value="reachability"']:
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
            '["abstract", "paperRad", "mechanism"].includes(mode)',
            "this.state.view.explodedSelected = false",
            "stopsVisible: true",
            "pivotsVisible: true",
            "linkagesVisible: true",
            "fastenersVisible: true",
            "actuatorsVisible: true",
            "measurementMode: \"all\"",
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
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["targetErrorVectorsVisible", "renderErrorRods", "mean signed"]:
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
            "finiteDieOffRadius",
            "lockAlphaCommutes",
            "finalHeightError",
        ]:
            self.assertIn(symbol, operators_js)

        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in [
            "commitSelectedLockEvent",
            "releaseSelectedLockEvent",
            "checkSelectedEventOrder",
            "selectedEventBaseline",
            "updateOperatorInspector",
            "RAD.applyProgrammableEvent",
            "RAD.compareEventOrder",
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
            "const deadZone",
            "backlashActivation(signal, deadZone) * zCouplingGain",
            "height[r][c] = locked ? 0 : -0.65 * influence[r][c] + zResidual[r][c]",
        ]:
            self.assertIn(symbol, math_js)
        for symbol in ["zCoupling: document.getElementById(\"zCoupling\")", "pinRadius: document.getElementById(\"pinRadius\")", "holeRadius: document.getElementById(\"holeRadius\")", "paperSideLengthMm: document.getElementById(\"paperSideLengthMm\")", "paperHoleToleranceMm: document.getElementById(\"paperHoleToleranceMm\")", "this.state.grid.zCouplingGain", "this.state.grid.pinRadius", "this.state.grid.holeRadius", "this.state.grid.paperSideLengthMm", "this.state.grid.paperHoleToleranceMm", "zCouplingOut", "pinHoleClearanceOut", "paperSideLengthMmOut", "paperHoleToleranceMmOut", "backlashMmOut", "pinHoleClearanceMmOut", "holeToleranceModelOut", "RAD.paperRadCalibration"]:
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
        for symbol in ["backlashStop", "updateBacklashStops", "state.view.stopsVisible", "this.materials.stop", 'state.view.cellVisualMode || "abstract") !== "abstract"']:
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
        for symbol in [
            "abstractBody",
            "abstractInner",
            "abstractEdge",
            "abstractCorner",
            "abstractGuide",
            "abstractDatum",
            "cellVisualMode",
            "abstractMode",
            "paperRadMode",
            "record.abstract",
            "record.paperRad",
            "Cell abstraction",
            "Paper RAD cell",
            "paperOuter",
            "paperInner",
            "paperClearance",
            "pinVisualRadius",
            "holeVisualRadius",
            "state.view.cellVisualMode || \"abstract\") !== \"abstract\"",
            "state.view.cellVisualMode || \"abstract\") === \"abstract\"",
        ]:
            if symbol in {"Cell abstraction", "Paper RAD cell"}:
                self.assertIn(symbol, (WEB / "index.html").read_text(encoding="utf-8"))
            else:
                self.assertIn(symbol, renderer_js)

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

    def test_browser_obj_export_controls_are_present(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        self.assertIn('id="saveObj"', html)
        self.assertLess(html.index("./analysis.js"), html.index("./physics.js"))
        self.assertLess(html.index("./physics.js"), html.index("./mesh_export.js"))
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
        self.assertLess(html.index("./physics.js"), html.index("./mesh_export.js"))

        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn('simulationMode: "kinematic"', state_js)

        physics_js = (WEB / "physics.js").read_text(encoding="utf-8")
        for symbol in [
            "simulatePhysicalRelaxation",
            "simulateActive",
            "spring-preview",
            "physicalRmsHeightDelta",
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
        for control_id in ["buildJacobian", "solveLinearFit", "applyLinearFit", "jacobianColumns", "meanReachability", "linearFitSteps", "linearFitError"]:
            self.assertIn(f'id="{control_id}"', html)
        self.assertIn('value="reachability"', html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        for symbol in ["jacobian:", "linearSolution:", "coverageMap", "columnCount", "meanCoverage", "conditionEstimate", "predictedError", "projectedError"]:
            self.assertIn(symbol, state_js)
        inverse_js = (WEB / "inverse.js").read_text(encoding="utf-8")
        for symbol in ["buildResponseJacobian", "solveLinearizedTargetFit", "applyLinearizedTargetFit", "linearized-jacobian-greedy-fit", "finite-difference-response-jacobian", "flattenDelta", "targetResidualVector", "heightDelta", "alphaDelta", "targetAlignment", "conditionEstimate", "commandZ", "commandAlpha", "z+", "z-", "alpha-", "alpha+", "RAD.buildResponseJacobian"]:
            self.assertIn(symbol, inverse_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        for symbol in ["buildJacobian", "solveLinearFit", "applyLinearFit", "jacobian-built", "linear-fit-solved", "linear-fit-applied", "jacobianColumns", "meanReachability", "linearFitSteps", "linearFitError", 'this.state.view.overlayMode = "reachability"']:
            self.assertIn(symbol, ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ['mode === "reachability"', "reachabilityStrength", "state.inverse?.jacobian", "linearSolution", "coverageMap", "maxCoverage"]:
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
        for option in ['value="all"', 'value="alpha"', 'value="theta"', 'value="height"', 'value="backlash"']:
            self.assertIn(option, html)
        state_js = (WEB / "state.js").read_text(encoding="utf-8")
        self.assertIn('measurementMode: "alpha"', state_js)
        self.assertIn("measurementLabelsVisible: false", state_js)
        ui_js = (WEB / "ui.js").read_text(encoding="utf-8")
        self.assertIn("this.els.measurementMode.addEventListener", ui_js)
        self.assertIn("showMeasurementLabels", ui_js)
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["measurementLabel", "measurementRingRadius", "renderThetaGuide", "renderBacklashGuide", "gapLine"]:
            self.assertIn(symbol, renderer_js)

    def test_renderer_uses_persistent_surface_meshes(self):
        renderer_js = (WEB / "renderer.js").read_text(encoding="utf-8")
        for symbol in ["syncSurfaceMeshes", "updateSurfaceMesh", "this.membraneMesh", "this.targetMesh", "surfaceContourRoot"]:
            self.assertIn(symbol, renderer_js)
        self.assertIn("position.needsUpdate = true", renderer_js)
        for symbol in ["renderSurfaceContours", "surfaceContourPoint", "addSurfaceContourLine", "targetContour"]:
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
