# Arbitrary Shape To RAD Lattice Spec

This document describes how a future agent should generate a Reconfigurable
Auxetic Devices lattice from an arbitrary target shape. The input can be a
rectangle, drawn silhouette, binary image mask, signed distance field, polygon,
3D mesh footprint, target surface, or procedural generated shape.

Core rule:

```text
Do not start by drawing cells. Start by defining the target domain, then build
a cell graph whose nodes and constraints approximate that domain.
```

## Expected Output

Given a target shape, the generator should produce:

- cells with row/column indices and physical center coordinates;
- a neighbor graph with explicit mechanical edges;
- pin constraints for every connected neighbor pair;
- boundary, corner, hole-boundary, support, removed-cell, and actuator labels;
- optional height/normal/curvature samples for a 3D target surface;
- validation metrics;
- a reproducible JSON state.

## Shape Inputs

### Rectangular Grid

Use this when the user asks for a regular sheet:

```text
rows: integer
cols: integer
cell_spacing_mode: "mechanism_pitch" or "fixed_pitch"
```

The current browser primitive is this type.

### Binary Mask

Use this when the input is an image or rasterized silhouette:

```text
mask(x, y) in {0, 1}
```

Candidate cell centers are sampled on a lattice basis. Keep a cell when its
center lies inside the mask or enough of its footprint overlaps the mask.

Recommended footprint rule:

```text
keep_cell = area(cell_footprint intersection shape) / area(cell_footprint) >= tau
tau default = 0.35
```

### Signed Distance Field

Use this for mathematically generated shapes:

```text
phi(x, y) < 0  inside
phi(x, y) = 0  boundary
phi(x, y) > 0  outside
```

Examples:

```text
circle:  phi(x, y) = sqrt(x^2 + y^2) - R
ellipse: phi(x, y) = sqrt((x/a)^2 + (y/b)^2) - 1
annulus: max(R_inner - sqrt(x^2 + y^2), sqrt(x^2 + y^2) - R_outer)
```

Signed distance fields are preferred for generated shapes because they support
holes, offsets, smooth boundaries, and deterministic procedural generation.

### Vector Boundary

Use this for SVG, CAD sketches, or user-drawn outlines:

```text
outer_boundary: closed polygon or spline
holes: optional closed polygons or splines
```

Convert to point-in-polygon plus distance-to-boundary queries.

### Mesh Footprint And Surface

Use this for a 3D target:

```text
mesh vertices/faces
projection_axis: usually Z
footprint = projection(mesh onto XY)
target_height(x, y) = sampled surface height
```

The footprint determines where cells exist. The height field determines future
vertical actuation or inverse-design targets.

### Procedural Generated Shape

Use this when no file exists and the user asks for a generated type:

```text
generator_name: "circle" | "ring" | "cross" | "wave_patch" | "branching" | ...
parameters: object
seed: integer
```

Every procedural generator must be deterministic from its parameters and seed.

## Coordinate Frames

Use these frames consistently:

```text
world X: column direction
world Y: row direction
world Z: vertical actuation direction
cell local x: east-west arm axis at zero lower-cross rotation
cell local y: north-south arm axis at zero lower-cross rotation
```

The displayed lattice basis is axis aligned by default:

```text
column basis = (pitch_x, 0, 0)
row basis    = (0, pitch_y, 0)
```

If a mechanism solve gives a tilted center-to-center vector, rotate it into
world +X for display. Internal cross rotations still come from the mechanical
closure.

## Cell Classification

For candidate center `c_ij`, classify the cell as:

```text
interior:      fully inside the target domain
boundary:      near or intersecting the outer boundary
hole_boundary: near or intersecting an internal hole
removed:       intentionally absent from the lattice
outside:       rejected candidate
anchor:        supported or position-constrained cell
actuator:      driven or actuator-candidate cell
```

Removed cells are topology changes, not visual hiding. If a cell is removed:

- do not render it;
- remove all incident mechanical edges;
- remove it from collision checks;
- remove its actuation variable;
- recompute neighboring boundary labels;
- preserve the removal in exported JSON.

## Neighbor Graph

A generated lattice is a graph:

```text
G = (V, E)
V = kept cells
E = mechanical neighbor connections
```

For a rectangular candidate grid:

```text
edge between (i, j) and (i, j+1) if both cells exist
edge between (i, j) and (i+1, j) if both cells exist
```

Do not create an edge if the connection path crosses removed space, a shape
hole, or an explicit cut.

Each edge stores:

```text
from_cell
to_cell
direction
constraint_type
pin_pairs
rest_pitch
backlash_gap
clearance
status
```

## No-Backlash Baseline

For the current no-backlash baseline:

```text
all connected cells use the same uniform actuation coordinate q
shared pins are hard equality constraints
pin residual should be zero within tolerance
```

The fixed A-B seed relation is:

```text
A upper east  = B lower south
B upper south = A lower north
```

Generated larger sheets should mark their fidelity honestly:

```text
fidelity: "normalized_no_backlash_proxy"
```

until CAD-derived 2D constraints or physical calibration data are added.

## Backlash Extension

Backlash belongs on graph edges:

```text
edge.backlash_gap
edge.pin_clearance
edge.dead_zone_operator
edge.engaged_state
```

The current dead-zone form is:

```text
f(x) = max(0, x - b) + min(x + b, 0)
```

Eventually each edge can have its own calibrated `b_e`.

## Surface Targets

For a 3D target surface, separate footprint from surface:

```text
footprint D: where cells exist in XY
height h(x, y): desired Z surface
normal n(x, y): optional desired normal
curvature k(x, y): optional desired curvature
```

Each kept cell receives:

```text
target_z_ij = h(c_ij)
target_normal_ij = n(c_ij)
target_curvature_ij = k(c_ij)
```

The UI should eventually support:

- cells only;
- surface only;
- cells plus surface skin;
- contours only;
- graph only;
- target error overlay.

## Required Metrics

Every generated lattice should report:

```text
cell_count
edge_count
removed_cell_count
boundary_cell_count
connected_component_count
largest_component_fraction
footprint_coverage_ratio
mean_neighbor_pitch
max_pin_residual
estimated_min_clearance
```

For target surfaces, also report:

```text
height_rms_error
height_max_error
normal_mean_error
target_reachable_fraction
underactuated_cell_count
```

