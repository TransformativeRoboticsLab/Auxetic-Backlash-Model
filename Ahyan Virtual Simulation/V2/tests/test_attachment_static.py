import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / "web"
SIM = WEB / "two-cell-attachment"


class CurrentSimulatorStaticTests(unittest.TestCase):
    def test_current_simulator_references_existing_local_assets(self):
        html = (SIM / "index.html").read_text(encoding="utf-8")
        refs = re.findall(r'(?:src|href)="([^"]+)"', html)
        ref_paths = [ref.split("?", 1)[0] for ref in refs]
        self.assertIn("../vendor/three.min.js", refs)
        self.assertIn("./styles.css", refs)
        self.assertIn("./two_cell_attachment.js", ref_paths)
        self.assertNotIn("cdn.jsdelivr.net", html)
        self.assertNotIn("../one-cell-rotation/index.html", html)
        self.assertNotIn("../index.html", html)
        for ref in ref_paths:
            if ref.startswith("data:"):
                continue
            self.assertTrue((SIM / ref).resolve().exists(), ref)

    def test_current_simulator_has_full_functionality_controls(self):
        html = (SIM / "index.html").read_text(encoding="utf-8")
        script = (SIM / "two_cell_attachment.js").read_text(encoding="utf-8")
        for control_id in [
            "aTopAngle",
            "actuationMode",
            "animate",
            "animationHalfPeriodSeconds",
            "rowCount",
            "colCount",
            "applyGrid",
            "showPins",
            "showCenterLines",
            "pinRadiusRatio",
            "plateThickness",
            "multiSelectMode",
            "positionLockMode",
            "lockSelectedXY",
            "unlockSelectedXY",
            "clearXYLocks",
            "selectedCellsMetric",
            "selectedLockMetric",
            "activeLocksMetric",
            "maxLiftMetric",
            "maxTiltMetric",
            "maxLockDriftMetric",
            "freeAngleMetric",
            "diametralBacklashMetric",
            "computedContactCommandMetric",
        ]:
            self.assertIn(f'id="{control_id}"', html)

        for symbol in [
            "function contactDriveField",
            "function boundaryActuationScale",
            "cellDriveCommands",
            "function selectedDriveCommand",
            "function seedContactCommandsFromUniformDrive",
            "selectedDriveCommandDirty",
            "input === aTopInput",
            'lastActuationMode !== "contact"',
            "function scanContinuityScenario",
            "function scanContinuityScenarios",
            "function lockedStringBacklashScenario",
            "function feasiblePinCenter",
            "function pinFitSnapshot",
            "latestPinFitReport",
            "maxAdmissibleEdgeRise",
            "backlashLimitedLift",
            "verticalRiseLimitForSpan",
            "function constrainedLayoutSolve",
            "function lockSelectedCellXY",
            "function unlockSelectedCellXY",
            'lockZ: positionLockModeInput.value !== "xy"',
            "function driveColor",
            "function applyCellColor",
            "function selectCellFromPointer",
            "multiSelectModeInput.checked",
            'actuationModeInput.value === "contact"',
            "window.RADAttachmentDiagnostics",
        ]:
            self.assertIn(symbol, script)

        self.assertNotIn('id="lockSelectedZ"', html)
        self.assertNotIn("lockSelectedZInput", script)
        self.assertNotIn('id="contactTestCommand"', html)
        self.assertNotIn("contactTestCommandInput", script)


if __name__ == "__main__":
    unittest.main()
