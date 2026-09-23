# Generate A RAD Lattice From A Shape

Follow this recipe when the user asks for a RAD lattice from an arbitrary
outline, generated shape, imported file, or target surface.

## Step 1: Normalize The Shape

Convert the input into one of:

```text
rectangular_grid
binary_mask
signed_distance_field
polygon_boundary
mesh_footprint_plus_height
procedural_generator
```

Keep the original file path, units, scale, seed, and parameters.

## Step 2: Choose Fidelity

Record the model fidelity:

```text
normalized_no_backlash_proxy
calibrated_no_backlash
backlash_graph_proxy
spring_hinge_proxy
cad_contact_calibrated
```

The current configurable browser lattice is:

```text
normalized_no_backlash_proxy
```

Do not claim fabrication accuracy unless CAD-derived constraints or physical
measurements were used.

## Step 3: Determine Cell Pitch

For the current no-backlash proxy, compute reference pitch from the two-pin
mechanism at a reference angle:

```text
q_ref = midpoint of feasible drive interval
pitch_ref = center distance from no-backlash closure at q_ref
```

Then set:

```text
col_basis = (pitch_ref, 0, 0)
row_basis = (0, pitch_ref, 0)
```

Future calibrated versions may use:

```text
pitch_x(q)
pitch_y(q)
```

## Step 4: Sample Candidate Centers

Compute shape bounds plus a one-cell margin:

```text
x_min, x_max, y_min, y_max
```

Generate candidate centers:

```text
c_ij = origin + i * row_basis + j * col_basis
```

For a generated shape, center the accepted cells around the world origin after
selection.

## Step 5: Keep Or Remove Cells

For each candidate cell, compute:

```text
center_inside
footprint_overlap_ratio
distance_to_boundary
inside_hole
```

First pass:

```text
keep if center_inside and not inside_hole
```

Better pass:

```text
keep if overlap_ratio >= 0.35 and not inside_hole
```

Classify each candidate as:

```text
interior
boundary
hole_boundary
removed
outside
anchor
actuator_candidate
```

## Step 6: Build The Graph

Create one graph node per kept cell. Add neighbor edges only when both cells
exist and the connection path does not cross a removed region or explicit cut.

Grid directions:

```text
east:  (i, j+1)
west:  (i, j-1)
north: (i-1, j)
south: (i+1, j)
```

## Step 7: Assign Pin Constraints

Use the fixed no-backlash seed for the current proxy:

```text
from upper east = to lower south
to upper south  = from lower north
```

For a true 2D CAD-derived solver, replace this with direction-specific edge
templates:

```text
east_edge_pin_template
west_edge_pin_template
north_edge_pin_template
south_edge_pin_template
```

Until then, mark edge constraints:

```text
constraint_source: "repeated_current_browser_proxy"
validated: false
```

## Step 8: Add Boundary Constraints

Boundary constraints are separate from topology:

```text
free_edge
fixed_position
fixed_height
wall_contact
support_plate
actuator_mount
symmetry_plane
```

Do not encode supports by secretly moving cells.

## Step 9: Add Surface Targets

If a target surface exists, sample:

```text
target_z[cell] = h(center_x, center_y)
target_normal[cell] = normal(center_x, center_y)
target_curvature[cell] = curvature(center_x, center_y)
```

If the target is only a 2D footprint, leave surface fields null.

## Step 10: Select Actuator Candidates

Allowed actuator states:

```text
unactuated
candidate_actuator
required_actuator
forbidden_actuator
```

Heuristics:

- prefer interior cells for smooth shape control;
- include boundary cells when edge lift is needed;
- forbid removed cells;
- respect underside access;
- respect maximum actuator count.

## Step 11: Validate Before Rendering

Compute:

```text
cell_count > 0
all edges connect existing cells
no duplicate edges
no edge crosses removed space
connected_component_count
largest_component_fraction
boundary_cell_count
removed_cell_count
```

If validation fails, show warnings in the UI and export them.

## Step 12: Export Reproducibly

Export:

```text
schema
created_at
generator_type
generator_parameters
shape_source
unit_scale
cell_geometry
cells
edges
removed_cells
boundary_constraints
actuator_candidates
target_surface
validation
fidelity
notes
```

If another agent cannot recreate the same lattice from the JSON alone, the
generator is incomplete.

