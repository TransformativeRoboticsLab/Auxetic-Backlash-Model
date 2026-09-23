# Three-Cell Uniform Row

## Goal

When a new cell is added in any direction in the no-backlash model, the added
cell should become part of the same connected kinematic mechanism. One
rotational actuation should propagate through all connected cells uniformly.

The current implemented row is:

```text
A - B - C
```

## Pattern

Cells A and C use the same A-pattern:

```text
A.lowerRotation = 0
A.upperRotation = q
C.lowerRotation = 0
C.upperRotation = q
```

Cell B is the connector pattern:

```text
B.lowerRotation = beta
B.upperRotation = beta + gamma
```

`beta` and `gamma` are solved from the A-B two-pin closure.

## Center Positions Before Display Reframing

The A-B loop closure produces:

```text
A.center = (0, 0)
B.center = c_B
C.center = 2 c_B
```

This makes A-B and B-C pitches equal.

## Display Reframing

For visual inspection, the row is moved into a fixed center-of-mass frame:

```text
A.displayCenter = -v
B.displayCenter =  0
C.displayCenter = +v
```

where `v` is the displayed A-B center vector.

The centerline direction is fixed on first render for the selected attachment
configuration. Later actuation changes scale `v` along that same line instead of
rotating the row in the viewport.

## Current Mirrored B-C Constraints

The B-C constraints are the opposite-side mirror of A-B:

```text
B lower north = C upper west
B upper north = C lower south
```

For the current default A-B choices, this preserves the same no-backlash
two-pin mechanism across the next neighbor.

## Uniform Contraction Requirement

A valid no-backlash row satisfies:

```text
distance(A.center, B.center) = distance(B.center, C.center)
A.center, B.center, C.center are collinear
centerOfMass(A, B, C) = (0, 0)
all pin residuals = 0
max inter-cell penetration = 0
```

During contraction:

```text
||v(q)|| decreases
```

and all connected cells move along the same centerline.

## Extension To A 2D Sheet

For a future sheet, every neighboring pair should impose equivalent hard pin
constraints. The no-backlash sheet should be solved as one connected constraint
graph:

```text
nodes = cell rigid bodies and cross rotations
edges = shared pin equality constraints
valid state = all constraints satisfied and no inter-cell penetration
```

Backlash should only be introduced after this exact no-backlash graph is
validated.

