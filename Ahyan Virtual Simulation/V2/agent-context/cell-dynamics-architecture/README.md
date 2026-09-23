# Cell Dynamics Architecture

This folder is a self-contained reconstruction guide for the current
Reconfigurable Auxetic Devices virtual simulation primitives. It is written for
an agent or collaborator with no prior chat context.

## Folder Map

```text
cell-dynamics-architecture/
  README.md
  handoff_prompt.md
  specs/
    geometry_and_state.md
    current_browser_primitives.md
    no_backlash_model_spec.json
    arbitrary_shape_lattice_spec.md
    generated_lattice_schema.json
  mechanics/
    single_cell_motion.md
    two_pin_neighbor_constraint.md
    three_cell_uniform_row.md
    contact_and_motion_limits.md
    shape_constrained_lattice_graph.md
  implementation/
    recreate_from_scratch.md
    current_files_map.md
    generate_lattice_from_shape.md
  validation/
    acceptance_tests.md
```

## Current Scope

The current primitive is a normalized, CAD-like, no-backlash simulation of RAD
cell motion. It is not yet a fabrication-accurate contact model. The goal of
this architecture is to preserve the mechanism logic clearly:

- A single RAD cell is represented as two identical stacked cross bodies.
- The lower cross and upper cross share the same center axis.
- The top cross rotates relative to the lower cross.
- Neighboring cells are coupled by shared pin-hole constraints.
- In the no-backlash baseline, shared pins are hard constraints with zero
  dead-zone and zero die-off.
- A row of connected cells should contract uniformly from one rotational
  actuation coordinate.
- Cell centers should move along fixed straight centerlines, not orbit the
  system center.
- Body overlap is forbidden by a primitive contact guard.

## Read Order

1. `specs/geometry_and_state.md`
2. `mechanics/single_cell_motion.md`
3. `mechanics/two_pin_neighbor_constraint.md`
4. `mechanics/three_cell_uniform_row.md`
5. `mechanics/contact_and_motion_limits.md`
6. `implementation/recreate_from_scratch.md`
7. `specs/arbitrary_shape_lattice_spec.md`
8. `mechanics/shape_constrained_lattice_graph.md`
9. `implementation/generate_lattice_from_shape.md`
10. `validation/acceptance_tests.md`

## Arbitrary Shape Generation

The architecture now includes a future-facing contract for building RAD
lattices from arbitrary shapes. A generated shape should become a reproducible
cell graph, not just a visual arrangement of crosses.

Generated outputs should preserve:

- the original shape source or procedural seed;
- the accepted, removed, boundary, and hole-boundary cells;
- all neighbor edges and pin constraints;
- the model fidelity level;
- target surface samples if the shape is 3D;
- validation metrics and warnings.

## Current Browser Entry Points

- `../../web/one-cell-rotation/index.html`
- `../../web/two-cell-attachment/index.html`
- `../../primitive-htmls/index.html`
