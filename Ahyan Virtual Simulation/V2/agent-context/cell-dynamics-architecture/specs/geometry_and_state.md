# Geometry And State

## Coordinate Convention

Use a right-handed 3D frame:

- `X` and `Y` define the sheet plane.
- `Z` is vertical.
- Rotations of cell crosses are in the `X-Y` plane about the local `Z` axis.
- Angles are stored in radians internally and displayed in degrees.
- Lengths are normalized millimeter-like values for the browser primitive.

## Current Normalized Geometry

These values match the current browser primitive constants:

```text
cellWidthMm = 55.604331
nominalHoleDiameterMm = 3.4
bodyThicknessMm = 4.0
padRadiusMm = 4.9
hubRadiusMm = 4.6
armWidthMm = 5.4
siteRadiusMm = 22.1
```

The four local arm-hole sites of a cross are:

```text
east  = ( siteRadiusMm, 0)
north = (0,  siteRadiusMm)
west  = (-siteRadiusMm, 0)
south = (0, -siteRadiusMm)
```

## Single Cell Body Layout

Each cell has two rigid cross bodies:

- lower cross: fixed relative orientation for the cell layer
- upper cross: rotates about the same center axis

The two crosses touch with zero visual gap:

```text
lower cross center z = -bodyThicknessMm / 2
upper cross center z =  bodyThicknessMm / 2
lower z band = [-bodyThicknessMm, 0]
upper z band = [0, bodyThicknessMm]
```

Within one cell, the two crosses are allowed to overlap in projection because
they are stacked layers sharing the same center axis. Between different cells,
solid body penetration is forbidden.

## Cell State Variables

For a no-backlash primitive, the minimum state is:

```text
cell.center = (x, y, z)
cell.lowerRotation = beta
cell.upperRotation = beta + gamma
cell.relativeRotation = gamma
```

For the first A-pattern cell in the current primitive:

```text
A.center = (0, 0) before display-frame recentering
A.lowerRotation = 0
A.upperRotation = q
```

For the B-pattern cell:

```text
B.center = solved from pin closure
B.lowerRotation = beta
B.upperRotation = beta + gamma
```

For the repeated C-pattern cell in the current no-backlash row:

```text
C.center = 2 * B.center before display-frame recentering
C.lowerRotation = 0
C.upperRotation = q
```

## Rendering Frame

The internal solver may produce a centerline at any world angle. The renderer
rotates and translates the whole solved pose into a fixed display frame:

- For two cells, centers should appear as `(-v/2, +v/2)`.
- For three cells, centers should appear as `(-v, 0, +v)`.
- The system center of mass stays at `(0, 0)`.
- During actuation, centers move along the same fixed line.

