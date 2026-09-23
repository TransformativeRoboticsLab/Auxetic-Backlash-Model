# Single Cell Motion

## Physical Abstraction

A RAD cell is modeled as two identical cross-shaped rigid bodies stacked on the
same central axis:

- lower cross
- upper cross

The upper and lower crosses are not separated. They occupy adjacent vertical
bands and rotate about the same center hole.

## Generalized Coordinate

Use one scalar rotational coordinate for a single cell:

```text
q = upper cross rotation relative to lower cross
```

For the A-pattern cell:

```text
lowerRotation_A = 0
upperRotation_A = q
```

The browser primitive uses `q` as the single visible actuation coordinate.

## Rigid Transform

For any local point `p = (x, y)` on a cross, its world position is:

```text
world(p; c, phi) = c + R(phi) p
```

where:

```text
c = cell center
phi = cross rotation
R(phi) = [[cos(phi), -sin(phi)],
          [sin(phi),  cos(phi)]]
```

## What Must Be Preserved

A correct single-cell recreation must preserve:

- two stacked cross bodies
- one common central rotation axis
- no top-bottom separation
- four arm-hole sites per cross
- same local geometry for upper and lower crosses
- visible rotation of the upper body relative to the lower body

## What Is Not Yet Modeled

The primitive does not yet claim:

- measured mass or inertia
- friction
- elastic deformation of the printed body
- exact Fusion body geometry
- exact pin-hole contact surfaces

