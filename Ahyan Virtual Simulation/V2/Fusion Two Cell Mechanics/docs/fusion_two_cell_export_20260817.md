# Fusion Two-Cell Export - 2026-08-17

## Source

- Fusion document: `TRL RADs two cells v2`
- Export runner: Fusion Text Commands
- Command:
  `Python.RunScript "C:\Users\ahyan\AppData\Roaming\Autodesk\Autodesk Fusion 360\API\Scripts\RADTwoCellExport\RADTwoCellExport.py"`
- Project script:
  `Fusion Two Cell Mechanics\fusion_scripts\export_two_cell_design_graph.py`

## Exported Files

- `Fusion Two Cell Mechanics\exports\two_cell_design_graph.json`
- `Fusion Two Cell Mechanics\exports\two_cell.step`

## Export Contents

- Parameters: 23
- Occurrences: 6
- Joints: 7
- Rigid groups: 0
- Motion links: 0
- Contact candidates: 6
- Pose samples: 6

Fusion did not expose named user parameters such as `actuation`, `alpha`,
`theta`, or `drive`. The exporter therefore sampled the first limited revolute
joint it found.

## Driven Pose Sweep

The sampled drive coordinate is:

- Joint: `Revolute 2`
- Motion type: `adsk::fusion::RevoluteJointMotion`
- Value: `rotationValue`
- Units: radians
- Limits: 0.4363323129985824 to 0.9773843811168246 rad
- Limits in degrees: 25.0 to 56.0 deg

Pitch is measured between the two top-level `RADs unit cell` occurrence
origins.

| Sample | Rotation (rad) | Rotation (deg) | Two-cell pitch (cm) |
|---:|---:|---:|---:|
| 1 | 0.4363323129985824 | 25.000 | 3.786935596859252 |
| 2 | 0.5715953300281429 | 32.750 | 3.942255994640801 |
| 3 | 0.747000919853573 | 42.800 | 4.114586005635456 |
| 4 | 0.842121364087264 | 48.250 | 4.196093960613915 |
| 5 | 0.9773843811168246 | 56.000 | 4.29393081917228 |

Current document pose before sampling had pitch 4.114586005635456 cm. The
exporter restored the driven joint to its original value after sampling.

## Mechanism Evidence

Fusion exposed one ball joint and six revolute joints. The limited revolute
joint above is the best current CAD-derived scalar actuation coordinate. The
export also stores each pose's occurrence transforms, joint axes, joint origins,
joint limits, cylindrical face hints, and STEP geometry.

The design contains no exported `MotionLink` objects, so the coupling is not a
Fusion motion-link equation. At this stage the exact simulator relation should
be fit from the sampled transforms/pitch values and then validated against a
held-out CAD pose or bench measurement.

## Next Simulator Step

Use `two_cell_design_graph.json` as the CAD ground truth for the two-cell
no-backlash mode:

1. Load the pose samples.
2. Fit a monotone interpolation from `Revolute 2.rotationValue` to top-level
   occurrence transforms and cell pitch.
3. Render the two top-level unit-cell instances from the STEP/transform data.
4. Keep the current abstract cell model available as a separate view.
5. Do not call the resulting model physically accurate until it reproduces a
   held-out Fusion pose or retroreflector measurement.

## Implemented Calibration Artifact

The simulator now builds a browser-loadable calibration artifact with:

```powershell
python -m rad_sim.build_two_cell_fusion_calibration --input "Fusion Two Cell Mechanics\exports\two_cell_design_graph.json" --out outputs\two_cell_fusion_calibration --web-data web\data
```

Generated files:

- `outputs\two_cell_fusion_calibration\two_cell_fusion_calibration.json`
- `outputs\two_cell_fusion_calibration\two_cell_fusion_calibration_summary.csv`
- `web\data\two_cell_fusion_calibration.json`
- `web\data\two_cell_fusion_calibration.js`

Validation summary:

- `Revolute 2.rotationValue` confirmed as the active limited revolute drive in
  this export.
- Drive range: 25.0 to 56.0 degrees.
- CAD pose samples: 5 driven samples.
- Pitch range: 37.86935596859252 to 42.9393081917228 mm.
- Left top-level cell origin drift: 0.0 mm.
- Right top-level cell origin drift: 11.999345312806689 mm.
- Current-pose interpolation error: 0.0 mm.

Remaining physical blockers are unchanged: no exported Fusion motion links, no
segmented pin-hole contact geometry, no friction/contact parameters, and no
held-out bench validation yet.
