# Recreate From Scratch

This is the minimum implementation path for an agent recreating the current
cell dynamics without prior context.

## Step 1: Build Rigid Cross Geometry

Create one cross body from:

- horizontal rectangular arm or capsule
- vertical rectangular arm or capsule
- center hub disk/cylinder
- four end pads with holes

Use the geometry constants from `../specs/geometry_and_state.md`.

## Step 2: Build One Cell

Create:

```text
cell.lowerCross
cell.upperCross
```

Place them on the same center axis:

```text
lower z = -bodyThickness / 2
upper z =  bodyThickness / 2
```

No separation gap.

## Step 3: Add One Uniform Drive Coordinate

Expose only:

```text
q = uniform rotational actuation
```

For A-pattern cells:

```text
lowerRotation = 0
upperRotation = q
```

## Step 4: Add B Connector Cell

For a neighboring B-pattern cell, solve:

```text
R(q) p_Au = c_B + R(beta) p_Bl
p_Al = c_B + R(beta + gamma) p_Bu
```

Do not expose `beta` or `gamma` as user controls in the no-backlash row. They
are derived values.

Also do not expose A-B attachment-site dropdowns. The no-backlash primitive uses
one fixed A-B closure:

```text
A upper east = B lower south
B upper south = A lower north
```

## Step 5: Add C Repeat Cell

Set:

```text
C.center = 2 * c_B
C.lowerRotation = 0
C.upperRotation = q
```

Add mirrored B-C constraints:

```text
B lower north = C upper west
B upper north = C lower south
```

## Step 6: Recenter For Display

For three cells:

```text
v = displayed center step
A.center = -v
B.center = 0
C.center = +v
```

Keep the display centerline direction fixed across actuation so centers move
straight along one line.

## Step 7: Enforce Validity

For each sampled `q`:

1. Solve loop closure.
2. Reject if any pin residual is above tolerance.
3. Reject if any inter-cell body penetration is above tolerance.
4. Keep only the allowed range.
5. Clamp the UI slider and animation to the allowed range.

## Step 8: Validate

Run the checks in `../validation/acceptance_tests.md`.
