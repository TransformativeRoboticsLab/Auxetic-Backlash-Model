# Current Files Map

## Primitive Launcher

```text
../../primitive-htmls/index.html
../../primitive-htmls/README.md
```

Purpose:

- Simple browser launcher for primitive pages.

## One-Cell Primitive

```text
../../web/one-cell-rotation/index.html
../../web/one-cell-rotation/one_cell_rotation.js
../../web/one-cell-rotation/styles.css
```

Important implementation concepts:

- one fixed lower cross
- one rotating upper cross
- shared center-hole axis
- zero top-bottom separation

## Three-Cell No-Backlash Primitive

```text
../../web/two-cell-attachment/index.html
../../web/two-cell-attachment/two_cell_attachment.js
../../web/two-cell-attachment/styles.css
```

Important implementation concepts:

- one uniform rotational actuation slider
- derived B lower and B upper rotations
- A-B two-pin loop closure
- mirrored B-C constraints
- C is an A-pattern repeat
- all three centers rendered on a fixed centerline
- primitive collision guard
- allowed drive interval display

## Validation Scripts

```text
../../tests/validate_one_cell_rotation_page.js
../../tests/validate_two_cell_attachment_page.js
../../tests/test_web_static.py
```

The two-cell attachment validation script now validates the three-cell
no-backlash row page because that page evolved from the earlier two-cell
primitive.

## Main Workbench

```text
../../web/index.html
```

This is the larger RAD digital workbench. It is not the source of truth for the
current no-backlash cell mechanism. The primitive pages should be used first
when correcting core mechanism motion.

