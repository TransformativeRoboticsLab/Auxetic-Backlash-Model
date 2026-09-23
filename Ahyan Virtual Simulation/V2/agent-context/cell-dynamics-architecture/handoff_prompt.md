# Handoff Prompt

Use this prompt for a future agent that has no chat context:

```text
You are working on the Reconfigurable Auxetic Devices Digital Workbench in:

C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation

First read:

agent-context/cell-dynamics-architecture/README.md
agent-context/cell-dynamics-architecture/specs/geometry_and_state.md
agent-context/cell-dynamics-architecture/mechanics/single_cell_motion.md
agent-context/cell-dynamics-architecture/mechanics/two_pin_neighbor_constraint.md
agent-context/cell-dynamics-architecture/mechanics/three_cell_uniform_row.md
agent-context/cell-dynamics-architecture/mechanics/contact_and_motion_limits.md
agent-context/cell-dynamics-architecture/specs/arbitrary_shape_lattice_spec.md
agent-context/cell-dynamics-architecture/mechanics/shape_constrained_lattice_graph.md
agent-context/cell-dynamics-architecture/implementation/generate_lattice_from_shape.md
agent-context/cell-dynamics-architecture/implementation/recreate_from_scratch.md
agent-context/cell-dynamics-architecture/validation/acceptance_tests.md

Your task is to preserve and extend the no-backlash RAD cell dynamics. The
current primitive models each cell as two stacked cross bodies. Neighboring
cells are connected by hard shared pin constraints. The current row primitive
has three cells, one uniform rotational actuation slider, mirrored A-B and B-C
two-pin constraints, fixed center-of-mass display framing, and a primitive
contact guard. Do not introduce backlash until the no-backlash hard-constraint
mechanism is validated.

If extending to more cells or a 2D sheet, propagate the same hard pin equality
constraints through the connected graph so all connected cells contract
uniformly. Validate pin residuals, equal neighbor pitch, collinearity or sheet
compatibility, and zero inter-cell penetration before changing the larger
workbench.

If asked to generate a lattice from an arbitrary shape, first convert the shape
to a domain representation, sample candidate cell centers, classify kept,
removed, boundary, and hole-boundary cells, then build a graph of cells and
mechanical neighbor constraints. Export enough JSON metadata that the same
lattice can be recreated without chat context. Mark the fidelity level honestly;
the current large-shape generator should be treated as a normalized
no-backlash proxy until CAD-derived 2D constraints and physical calibration are
available.
```
