# Two-Pin Neighbor Constraint

## Purpose

Two neighboring cells should not be connected by a single loose visual overlap.
They should be connected by pin-hole equality constraints. In the current
no-backlash primitive, two neighboring cells are connected by two shared pin
relationships.

## A-B Constraint Variables

Let:

```text
q     = A upper rotation
beta  = B lower rotation
gamma = B upper rotation relative to B lower
c_B   = B center relative to A center
```

Use local sites:

```text
p_Au = A upper east site
p_Bl = B lower south site
p_Al = A lower north site
p_Bu = B upper south site
```

These are fixed mechanism definitions in the no-backlash browser primitive, not
user-selectable attachment choices.

```text
p_Au = east
p_Bl = south
p_Al = north
p_Bu = south
```

## Hard Pin Equations

The primary pin is:

```text
R(q) p_Au = c_B + R(beta) p_Bl
```

So:

```text
c_B = R(q) p_Au - R(beta) p_Bl
```

The second pin is:

```text
p_Al = c_B + R(beta + gamma) p_Bu
```

Together these two constraints close the A-B loop.

## Solving The Closure

The browser primitive solves the loop by reducing it to a two-circle closure.

Define:

```text
d = p_Al - R(q) p_Au
r_l = ||p_Bl||
r_u = ||p_Bu||
```

For the current symmetric sites:

```text
r_l = r_u = siteRadiusMm
```

If `||d|| > r_l + r_u`, the requested angle cannot close the loop.

If the loop can close, two geometric branches can exist. The solver chooses the
branch closest to the reference B orientations:

```text
beta_reference = 0 deg
gamma_reference = -25 deg
```

These reference angles are not user controls. They only select the visually
consistent branch of the no-backlash closure.

## No-Backlash Meaning

No backlash means the pin equality constraints are exact:

```text
pin_residual <= 0.01 mm
```

There is no dead-zone function, no delayed coupling, and no distance-based
die-off in this primitive.
