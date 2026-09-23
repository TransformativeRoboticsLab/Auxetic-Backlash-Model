# Fusion Two Cell Mechanics

This folder is the handoff packet for connecting the real Autodesk Fusion
two-cell RAD mechanism to the Reconfigurable Auxetic Devices Digital Workbench.

The immediate target is not the full backlash lattice. The target is the exact
kinematic function of two connected, uniform, non-backlash cells: how the bodies,
pins, holes, hinges, and shared connector geometry move as the actuation
parameter changes.

## Why This Exists

The current browser simulator is useful for inspection, UI development, and
normalized backlash experiments, but its two-cell expansion behavior is still a
proxy. To make it accurate, the simulator needs a CAD-derived two-cell ground
truth model before adding backlash, locks, vertical clearance, contact, or neural
surrogates.

## Folder Contents

- `raw/cad_fusion_codex_workflow_source.txt`: copied source note from the
  conversation about Fusion, Codex, MCP/API access, JSON graph exports, and STEP
  exports.
- `docs/fusion_to_simulation_workflow.md`: cleaned workflow for turning the
  Fusion model into simulator equations and validation data.
- `docs/two_cell_mechanics_data_contract.md`: the exact data we should extract
  from Fusion for a two-cell model.
- `schemas/two_cell_fusion_export.schema.json`: a structured JSON target for the
  Fusion export.
- `fusion_scripts/export_two_cell_design_graph.py`: starter Fusion Python script
  for exporting parameters, components, joints, transforms, and STEP geometry.
- `prompts/goal_fusion_two_cell_kinematics.md`: next goal prompt to run when the
  Fusion two-cell file or exported packet is available.
- `exports/`: put generated `.json`, `.step`, `.stl`, `.obj`, `.csv`, or sampled
  pose files here.

## Existing Related Assets

The broader simulator project already contains one-cell CAD reference material:

- `C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation\assets\cad`
- `C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation\cad_probe`
- `C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation\outputs\two_cell_bench_packet`

Those assets should be treated as reference material only. The two-cell
non-backlash Fusion model should become the first source of truth for correcting
the pair mechanics.

## Practical Next Step

Open the exact two-connected-cell Fusion file, run the exporter script or export
the equivalent data through Fusion MCP/API, place the outputs in `exports/`, and
then run the goal prompt in `prompts/goal_fusion_two_cell_kinematics.md`.
