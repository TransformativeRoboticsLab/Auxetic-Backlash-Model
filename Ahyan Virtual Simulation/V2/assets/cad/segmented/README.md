# Segmented CAD Intake

Place Fusion-exported one-cell and two-cell assets here when moving from the
current reduced proxy to external rigid-body contact simulation.

Expected filenames:

- `upper_free_cell_body.step` or `.stl` or `.obj`
- `lower_cell_body.step` or `.stl` or `.obj`
- `screw_pin_body.step` or `.stl` or `.obj`
- `radial_pad_hole_surfaces.step` or `radial_pad_hole_surfaces.json`
- `joint_axes.json`
- `two_cell_connector_pairs.json`
- `mass_inertia.json`
- `contact_parameters.json`
- `lock_crown_geometry.json` or `.csv`
- `actuator_force_displacement.json` or `.csv`
- `bench_coordinate_truth.json` or `.csv`
- `mujoco/rad_two_cell.xml`
- `gazebo/rad_two_cell.sdf` or `rad_two_cell.urdf`
- `isaac/rad_two_cell.usd`

Run `python -m rad_sim.export_two_cell_bench_packet` after adding files. The
generated `two_cell_segmented_cad_readiness.json` will mark detected assets and
show detected filenames, while `two_cell_segmented_cad_intake_validation.json`
parses the filled JSON/CSV contents and blocks exact contact simulation when
required CAD-measured values are still blank.

The packet export also writes starter versions of the JSON/CSV measurement
files into `outputs/two_cell_bench_packet/segmented_cad_intake_templates/`.
Use those as fillable templates, then copy the filled files back into this
folder using the same names.

The packet also writes `two_cell_cad_contact_decomposition.json` and
`two_cell_cad_contact_decomposition.csv`. Those files define the current
CAD-to-engine contact contract: two one-cell instances, upper/lower/pin body
roles, 16 pad holes, three X-neighbor connector pairs, and segmented
convex/analytic hole-wall primitives. Use that decomposition for active
collision bodies; keep the full CAD mesh as a visual/reference mesh until its
hole walls are split into engine-safe contact primitives.

The validation gate requires:

- the three separated body assets: upper free cell, lower cell, and screw/pin;
- filled joint axes, vertical travel, eight pin/hole axes, hole surfaces, and
  two-cell connector origins;
- filled mass/inertia and contact/friction parameters;
- filled lock crown, actuator, and bench coordinate CSV measurements;
- at least one external engine handoff file, such as `mujoco/rad_two_cell.xml`.

Passing this intake gate only means an exact contact run can be attempted. A
real-life accuracy claim still requires external engine results and bench
coordinate comparison.
