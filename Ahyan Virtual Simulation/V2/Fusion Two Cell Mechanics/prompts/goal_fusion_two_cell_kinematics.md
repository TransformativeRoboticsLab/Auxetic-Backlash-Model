# Goal Prompt: Fusion-Derived Two-Cell RAD Kinematics

Use this as the next Codex goal after the exact two-connected-cell Fusion file
or a Fusion export packet has been placed in:

`C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation\Fusion Two Cell Mechanics\exports`

## Prompt

Work in `C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation`.

The project is the Reconfigurable Auxetic Devices Digital Workbench. The current
browser and Python simulator contains a normalized RAD lattice and a proxy
two-cell bench, but the two-cell expansion mechanics are not accurate enough. I
want the next phase to be math-first and CAD-grounded.

Your objective is to connect the exact function of the Autodesk Fusion model for
two connected, uniform, non-backlash RAD cells to the simulator. Treat the Fusion
model or exported Fusion graph as the source of truth. Do not infer exact
mechanics from screenshots or the current proxy model. The first milestone is
only the non-backlash two-cell mechanism. Backlash, locks, vertical pin-hole
clearance, contact, gravity, sheet-level propagation, and neural surrogates come
after the no-backlash pair is validated.

Start by reading:

- `Fusion Two Cell Mechanics\README.md`
- `Fusion Two Cell Mechanics\docs\fusion_to_simulation_workflow.md`
- `Fusion Two Cell Mechanics\docs\two_cell_mechanics_data_contract.md`
- `Fusion Two Cell Mechanics\schemas\two_cell_fusion_export.schema.json`
- every file currently in `Fusion Two Cell Mechanics\exports`
- existing two-cell simulator code, especially `rad_sim\two_cell_bench.py`,
  `web\two_cell_bench.js`, `web\math.js`, and related tests

If the Fusion export packet is incomplete, create a precise missing-data report
and update the exporter or instructions rather than guessing. If Fusion MCP/API
is available, use it to query the live Fusion session for components, joints,
parameters, transforms, and pose samples. If Fusion MCP/API is not available,
use the provided export schema and Fusion script workflow.

Derive a reduced two-cell kinematic model:

- define the actuation coordinate used by the CAD model;
- define the physically meaningful alpha and theta for each cell;
- define left and right cell centers;
- define pitch as the distance between the connected cell centers;
- define rigid body transforms for each moving body or plate;
- define closure constraints for all shared joints/connectors;
- define tracked points that can later align with retroreflector data;
- fit or derive the maps from actuation to center spacing, plate rotations, and
  cell state;
- use analytic geometry where possible and interpolation only as a temporary
  bridge;
- do not add backlash or empirical dead zones to this baseline.

Then integrate the result carefully:

- replace arbitrary two-cell expansion constants in the two-cell bench with the
  CAD-derived no-backlash functions;
- keep the current normalized/abstract cell model selectable;
- add a clear browser mode for exactly two connected cells;
- expose the CAD-derived pair metrics in the UI: actuation, alpha, theta, pitch,
  center coordinates, closure error, and fitted/held-out error;
- add tests that prove the no-backlash two-cell model has no artificial dead
  zone, preserves symmetric behavior under symmetric commands, and reproduces
  the provided CAD pose samples within reported tolerance;
- keep existing Python and browser validation passing.

Report the governing equations explicitly. Separate what is proven from what is
fit from data. Mark any simulator feature that remains a proxy. The output
should include code changes, tests, a short engineering note explaining the
two-cell equations, and a list of missing measurements needed before calling the
model physically accurate.

Run validation with the narrow relevant tests first, then the broader available
test suite:

- `python -m unittest discover -s tests`
- `node tests\validate_web_modules.js` if Node is available

Stop when the two-cell non-backlash model is integrated and validated, or when
the blocker is a missing Fusion/CAD export that cannot be recovered locally. Do
not claim the model is exact unless it has been checked against CAD pose samples
or bench retroreflector data.
