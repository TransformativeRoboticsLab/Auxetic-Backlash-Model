# RAD Virtual Simulation V2

The current pin-hole lattice simulator is `web/two-cell-attachment/index.html`.
Open that file directly in Chrome; its scripts, styles, and Three.js dependency
are local. Other existing V2 entry points are legacy implementations and are
not updated by this snapshot.

## Single positional anchor

One XY or XYZ position lock preserves the same actuation range and internal
geometry as the unlocked lattice. The solver computes the free configuration,
then translates all cell centers and pins together to keep the selected center
at its captured coordinates. XY does not constrain height; XYZ does. The lock
does not freeze the cell's internal rotational actuation. Pin history is kept
in the same translated frame. Two or more locks retain the constrained solver
and its range/jamming safeguards.

This snapshot includes the current boundary/clearance solver and its required
`backlash_state.js` and `pin_geometry.js` modules, not only a standalone patch
against the older V2 implementation.

## Validation

From this directory:

```text
python -m unittest discover -s tests -p test_attachment_static.py
node tests/test_backlash_state.js
node tests/test_pin_geometry.js
node tests/validate_single_anchor.js
```

The browser test requires Node.js, the Playwright package, and installed Chrome.
It compares anchored states with unlocked references, checks captured lock
coordinates and pin containment, and exercises live selection, sequential
actuation, animation, resize, and unlock/re-anchor behavior. It also checks
desktop/mobile screenshots and nonblank canvas pixels.

Coverage includes XY and XYZ locks, horizontal/vertical strings through 30 cells,
rectangular and square sheets, and 30x30 full-travel center-anchor checks. The
larger grids are computationally expensive; these finite tests are not an
exhaustive proof for every configuration.

## Known limitations

This is an uncalibrated kinematic model, not full force/contact dynamics or a
guarantee of body-body collision freedom. A separate randomized audit found a
roughly 0.038 mm quarter-turn discrepancy in one free-sheet case. It also occurs
without locks; the anchored result matches that free solution up to translation.
Pin containment and lock-coordinate checks passed across that 60-case audit,
but the complete audit's rotation-invariance assertion did not pass.
