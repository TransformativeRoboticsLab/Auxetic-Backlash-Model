# Cylindrical Wheel Tiling Simulator

Author: Ahnaf Inkiad.

This is a standalone, interactive browser simulator for the reconfigurable
soft robotic wheel: `n` RAD (Reconfigurable Auxetic Device) cells hinged
edge-to-edge into a closed ring (the wheel's cross-section), with `m` such
rings stacked along an axle direction into a full cylindrical tiling. It is
a genuine 3D structure (cells at varying x, y, and z), not a flat lattice -
open `index.html` to see it rendered and orbit the camera around it.

## Why this exists

Jacob asked for a feature to simulate non-2D (cylindrical) arrangements of
the auxetic tiles for a reconfigurable wheel. The existing flat-plane
simulators in this repo (`main_in_3D.py`'s dome approximation, and the RAD
digital workbench's rectangular lattice in `../Ahyan Virtual Simulation/`)
lay cells out on a flat or gently-curved sheet. This project instead closes
the cells into an exact ring and stacks rings into a tube, matching the
physical hardware (see the bench photos of the assembled auxetic ring in the
project chat) far more directly than a flat-sheet approximation.

## Relationship to `Ahyan Virtual Simulation/`

This folder is independent - its own top-level directory, not nested inside
`Ahyan Virtual Simulation/` - so it's clear this is a separate contribution.
It does reuse two things from that project for consistency with the rest of
the RAD tooling:

- `vendor/three.min.js`: the same vendored Three.js build (a local copy, not
  a shared reference), so this page works offline with no build step, the
  same way the other RAD pages do.
- The RAD cell's real CAD dimensions (`cellWidthMm`, `siteRadiusMm`, arm/pad/
  hole sizes) and the paper-supported `theta = 70*alpha - 60` angle law
  documented in `../Ahyan Virtual Simulation/docs/research_grounding.md`, so
  a ring built here is dimensionally comparable to the row/one-cell
  primitives there rather than using different made-up numbers.

The ring-closure math, the cylinder stacking, the backlash dead-zone
handling, collision checking, cell picking, and the UI are new work, not
adapted from that project's files.

## Open it

Static HTML/CSS/JS, no build step, no server required:

```powershell
start index.html
```

## What it models

- **Ring closure**: `n` cells are placed by walking the same single-pin
  attach relation used elsewhere in the RAD tooling (a fixed local step,
  rotated by a turn angle of `2*pi/n` at each step) all the way around a
  loop. A regular n-gon closes exactly by construction - the "ring closure
  residual" readout is a live numerical check of that, not a fixed display
  value.
- **Drive**: the paper-supported dilation factor `alpha`, converted to the
  cell rotation `theta = 70*alpha - 60`, with a bidirectional backlash
  dead-zone (`f(x) = max(0, x-b) + min(x+b, 0)`) applied around `alpha = 1`
  so the wheel doesn't respond to small commands inside its mechanical
  backlash gap.
- **Cylinder stacking**: `m` rings are repeated along the axle (`z`) axis at
  a configurable pitch. Axial row-to-row attachment is currently a rigid
  copy, not a solved joint - see Known limitations.
- **Material Restriction**: neighbor-only clearance checking (circumferential
  and axial neighbors only, not every cell pair) using the same disk/capsule
  primitive distance math as the row primitive.
- **Selection tools**: click a cell to inspect it, then Frame Cell (recenter
  the camera on it) or Isolate Cell (hide everything else) from the Selected
  Cell panel. Focus hides the control panel for a clean view.
- **Save / Load**: JSON export/import of alpha, backlash, ring/row counts,
  axial pitch, and attachment sites.

## Known limitations (intentional, not yet done)

- Only a single-pin joint per cell-to-cell connection is modeled, not the
  row primitive's mirrored double-pin closure. Generalizing double-pin
  closure to an arbitrary n-cell ring is a real nonlinear loop-closure
  problem (2 unknowns and 2 constraints per joint); it was left for later
  rather than shipped half-verified.
- Axial (row-to-row) attachment is a rigid z-offset copy, not a solved
  joint - there's no axial backlash or clearance-driven pitch yet.
- No inverse design, timeline/keyframes, or Jacobian sensitivity analysis,
  unlike the much larger RAD digital workbench in `../Ahyan Virtual
  Simulation/web/`. Those are large, separate undertakings.

## Test

```powershell
npm install playwright pngjs
npx playwright install chromium   # once
node tests/validate_cylinder_tiling_page.js
```

The Playwright script loads the page in a real Chromium browser (desktop and
mobile viewports), checks the ring closure residual, the paper angle law,
the backlash dead-zone, neighbor clearance reporting, cell picking, Focus/
Frame/Isolate behavior, and a Save/Load JSON round-trip, then screenshots
the rendered cells and checks their colors.
