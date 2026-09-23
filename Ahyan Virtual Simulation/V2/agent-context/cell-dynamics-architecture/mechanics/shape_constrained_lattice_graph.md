# Shape-Constrained Lattice Graph Mechanics

This document gives the mathematical model for building an auxetic lattice from
an arbitrary shape.

## Graph Model

A lattice is an embedded graph:

```text
G = (V, E)
```

Each vertex `v in V` is one RAD cell. Each edge `e in E` is a mechanical
connection between two neighboring cells.

Each cell has state:

```text
c_v in R^3       center position
R_v in SO(3)     cell body orientation
q_v in R         upper/lower relative rotation
z_v in R         vertical displacement or command
locked_v         boolean or constraint set
removed_v        boolean
```

For the current 2D no-backlash browser primitive:

```text
c_v = (x_v, y_v, 0)
R_v = rotation about world Z
q_v = uniform drive q
z_v = 0
```

## Cell Bodies

Each cell contains two rigid cross bodies:

```text
lower_cross_v
upper_cross_v
```

The transforms are:

```text
T_lower(v) = T(c_v) R_v
T_upper(v) = T(c_v) R_v RotZ(q_v)
```

The lower and upper crosses share the center revolute axis.

## Hole Sites

Each cross has four idealized hole sites:

```text
east  = ( r,  0, 0)
north = ( 0,  r, 0)
west  = (-r,  0, 0)
south = ( 0, -r, 0)
```

A site in world coordinates is:

```text
site_world(v, layer, site) =
  c_v + R_v * layer_rotation(q_v, layer) * site_local(site)
```

where:

```text
layer_rotation(q, lower) = I
layer_rotation(q, upper) = RotZ(q)
```

## Hard Pin Constraint

A no-backlash shared pin is a point equality:

```text
site_world(u, layer_a, site_a) = site_world(v, layer_b, site_b)
```

The pin residual is:

```text
rho_e = || site_world(u, layer_a, site_a)
          - site_world(v, layer_b, site_b) ||
```

In this simulator, pin residual is geometric constraint error, not the physical
clearance between a pin and a hole.

Acceptance:

```text
max_e rho_e <= pinResidualToleranceMm
```

## Two-Pin Neighbor Edge

A realistic neighbor pair needs two pin equalities. One pin alone leaves an
unwanted relative rotation.

Current fixed A-B no-backlash seed:

```text
A upper east  = B lower south
B upper south = A lower north
```

For arbitrary lattice generation, every edge should define:

```text
pin_1: (u, layer, site) = (v, layer, site)
pin_2: (u, layer, site) = (v, layer, site)
```

until CAD-derived direction-specific templates are available, generated 2D
edges should be marked as proxy constraints.

## Shape Domain

Let `D` be the desired planar footprint. It may come from a polygon, mask,
signed distance field, generated function, or mesh projection.

Candidate centers:

```text
c_ij = origin + i * row_basis + j * col_basis
```

Accepted cells:

```text
V = { (i, j) : c_ij passes the shape keep rule }
```

Edges:

```text
E = adjacent kept cells whose connection path is allowed
```

## Boundary Cells

A cell is a boundary cell if:

```text
distance(c_v, boundary(D)) <= cell_effective_radius
```

Boundary cells matter because they have fewer neighbors, may need supports, and
may respond differently from interior cells.

## Discontinuity Operators

Removed cells, locks, cuts, backlash changes, and supports should be modeled as
operators on the graph and feasible configuration set.

Examples:

```text
RemoveCell_k(G) = remove node k and all incident edges
CutEdge_e(G) = remove mechanical edge e
LockCell_k(Q) = add angular/state constraint to cell k
PinCellPosition_k(Q) = fix c_k
ChangeBacklashGap_e(Q) = modify dead-zone on edge e
AddActuator_k(Q) = add command variable to cell k
```

The feasible set is:

```text
Q(G) = { states satisfying pin constraints,
         contact constraints,
         locks,
         actuator bounds,
         and boundary constraints }
```

A programmable discontinuity operator maps:

```text
P: Q(G) -> Q(P(G))
```

This is the intended mathematical direction for the broader framework.

## Target Surface Problem

For a desired surface height `h(x, y)`:

```text
z_v_target = h(c_v.x, c_v.y)
```

The forward simulator predicts:

```text
z_v_actual = F_v(commands, locks, removed_cells, boundary_constraints)
```

The inverse design problem is:

```text
minimize sum_v w_v * (z_v_actual - z_v_target)^2
subject to pin constraints, contact constraints, actuator limits,
support constraints, and maximum actuator count
```

This should be solved only after local one-cell, two-cell, and string mechanics
are calibrated.

## Conformal And Metric View

For shape-conforming lattices, treat the generated sheet as a map:

```text
F: reference_grid_domain -> target_domain
```

The target metric is:

```text
g = J_F^T J_F
```

where `J_F` is the Jacobian of the map.

The auxetic mechanism provides a discrete mechanism metric depending on:

```text
q
locks
backlash gaps
removed cells
boundary constraints
```

Future inverse design should compare the target metric to the mechanism metric.
Use conformal maps when local angle preservation matters, quasi-conformal maps
when anisotropic stretch is allowed, and graph constraints when the discrete
pin topology dominates.

