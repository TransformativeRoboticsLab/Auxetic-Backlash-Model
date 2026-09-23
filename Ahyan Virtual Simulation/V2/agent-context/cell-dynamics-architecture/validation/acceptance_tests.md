# Acceptance Tests

Use these checks before claiming that the no-backlash primitive was recreated
correctly.

## Static Checks

From:

```text
C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation
```

Run:

```powershell
python -m unittest tests.test_web_static
```

Expected:

```text
OK
```

## Browser Primitive Checks

Run the one-cell check:

```powershell
$env:NODE_PATH="C:\Users\ahyan\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\node_modules"
& "C:\Users\ahyan\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe" tests\validate_one_cell_rotation_page.js
```

Run the row check:

```powershell
$env:NODE_PATH="C:\Users\ahyan\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\node_modules"
& "C:\Users\ahyan\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe" tests\validate_two_cell_attachment_page.js
```

Expected:

```text
one-cell rotation page validation passed
two-cell attachment page validation passed
```

## Visual Acceptance

Open:

```text
../../primitive-htmls/index.html
```

Confirm:

- The one-cell page shows two stacked crosses sharing one center axis.
- The row page shows exactly one rotational actuation slider.
- A, B, and C centers stay on one straight line.
- The middle cell stays at the system center in the three-cell display frame.
- A-B and B-C center spacing are equal.
- All pin residuals read `0.000 mm`.
- Max penetration reads `0.000 mm`.
- The allowed drive range is shown and the slider cannot leave it.
- The row contracts uniformly, with no backlash, no die-off, and no local delay.

## Current Expected Default Readouts

For the current default row, the exact rounded display values can vary slightly
after renderer or camera changes, but the important invariants are:

```text
allowed uniform drive: includes 26 to 64 deg
A-B pin residual: 0.000 mm
second pin residual: 0.000 mm
B lower to C upper: 0.000 mm
B upper to C lower: 0.000 mm
max penetration: 0.000 mm
A-B center pitch = B-C center pitch
system center = (0.0, 0.0)
```

## Failure Conditions

Do not accept the implementation if:

- B or C has a user-controlled independent motion slider.
- Cell A is fixed while the other cells move.
- Centers orbit the midpoint during actuation.
- A-B and B-C center distances differ in the no-backlash row.
- Any pin residual is nonzero in a nominal valid pose.
- Any inter-cell body penetration is nonzero in a valid pose.
- Backlash, dead-zone, or die-off is active in this baseline primitive.

## Shape-To-Lattice Acceptance Tests

Use these checks before claiming that an arbitrary target shape has been turned
into a valid RAD lattice.

### Rectangular Shape

Input:

```text
rows = 3
cols = 5
```

Expected:

```text
cell_count = 15
horizontal_edges = rows * (cols - 1) = 12
vertical_edges = cols * (rows - 1) = 10
connected_component_count = 1
```

### Circle Footprint

Input:

```text
shape = circle radius R
cell_pitch = p
```

Expected:

- accepted cells lie inside or near the circle;
- boundary cells are detected near the shape edge;
- no edge connects through outside removed space;
- largest component contains most cells unless the radius is too small.

### Annulus Or Hole

Input:

```text
shape = ring with inner and outer radius
```

Expected:

- cells inside the inner hole are removed;
- edges crossing the inner hole are absent;
- cells around the inner hole are marked `hole_boundary`;
- exported JSON preserves the hole or the removed-cell set.

### Generated Shape Reproducibility

Input:

```text
generator_name
parameters
seed
```

Expected:

- running the generator twice with the same seed produces identical cells;
- changing the seed changes only stochastic/generated choices;
- exported JSON contains the seed and parameters.

### Surface Target

Input:

```text
footprint D
height field h(x, y)
```

Expected:

- every kept cell has `targetZ`;
- contour view can be generated without showing cells;
- target error fields are null until a forward simulation is run;
- inverse-design code can read the target without inspecting UI state.

### Removed Cell Topology

Input:

```text
rectangular grid plus removed cell set
```

Expected:

- removed cells are not rendered;
- incident edges are removed;
- neighboring cells become boundary or hole-boundary cells;
- connected components are recomputed.
