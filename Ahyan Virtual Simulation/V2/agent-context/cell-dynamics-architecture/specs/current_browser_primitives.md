# Current Browser Primitives

## One-Cell Primitive

Location:

```text
../../web/one-cell-rotation/index.html
```

Purpose:

- Shows a single RAD cell as two stacked identical cross bodies.
- The lower cross is fixed.
- The upper cross rotates about the shared center-hole axis.
- There is no separation between top and bottom crosses.

Use this page to verify the local cell concept before adding neighboring
constraints.

## Three-Cell No-Backlash Row Primitive

Location:

```text
../../web/two-cell-attachment/index.html
```

Despite the folder name, this page currently contains the no-backlash row
primitive with three connected cells: `A - B - C`.

Purpose:

- Demonstrates how one rotational actuation propagates uniformly through a row.
- Uses two pin constraints between A and B.
- Uses mirrored two pin constraints between B and C.
- Keeps A, B, and C centers collinear.
- Keeps A-B and B-C center pitch equal.
- Uses one visible user control: uniform rotational actuation.
- Clamps the actuation range so the hard pin constraints and contact guard stay
  feasible.

## Current Cross Colors

The current browser primitive colors the six cross bodies separately:

```text
Cell A upper cross: orange
Cell A lower cross: dark gray
Cell B upper cross: blue
Cell B lower cross: purple
Cell C upper cross: green
Cell C lower cross: rose
```

The current browser primitive now uses a neutral CAD-style material scheme
instead of exposing color-coded cross types as a feature. Earlier color naming
was only a debugging aid while the cross-to-cross relationships were being
identified.

## Fixed A-B Attachments

The A-B attachments are fixed in the simulator. They are not user-selectable
dropdown options.

Primary A-B pin:

```text
A upper east arm = B lower south arm
```

A-B loop-closure pin:

```text
B upper south arm = A lower north arm
```

In the old sidebar wording, this means:

```text
Cell A upper cross hole: east arm
Cell B lower cross hole: south / lower-facing arm
Cell B upper cross hole: south / lower-facing arm
Cell A lower cross hole: north / upper-facing arm
```

These are no-backlash hard constraints in the primitive and should not be
changed without creating a separate experimental branch/page.

Every rendered cell should also show an internal center-axis pin through the
overlapped upper and lower crosses. This center pin represents the revolute
axis of the single RAD cell itself. Neighbor-to-neighbor connections should
show two shared pins per adjacent cell pair.

The browser display frame is axis aligned: the solved A-to-B center vector is
rotated onto world +X, and repeated rows use world +Y. Internal part rotations
still come from the no-backlash closure; only the displayed lattice basis is
aligned to the viewer axes.

## Repeated Row Attachments

B-C mirrored primary pin:

```text
B lower north arm = C upper west arm
```

B-C mirrored second pin:

```text
B upper north arm = C lower south arm
```

These are derived repeat constraints used for the row/lattice proxy.
