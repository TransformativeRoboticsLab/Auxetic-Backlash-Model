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

Everything else - the ring-closure math, cylinder stacking, backlash
coupling/propagation, per-cell actuation, collision checking, cell picking,
timeline, target-diameter fitting, and the UI - is new work, not adapted
from that project's files.

## Open it

Static HTML/CSS/JS, no build step. Opening `index.html` directly via
`file://` works, but Chrome treats every `file://` reload as a distinct
security origin and can block a same-tab refresh outright (harmless but
confusing - closing and reopening the tab works around it). A local server
avoids that entirely and is the more reliable way to keep working in it:

```powershell
python -m http.server 8000
# then open http://localhost:8000/ in a browser
```

## What it models

- **Ring closure**: `n` cells are placed by walking the same single-pin
  attach relation used elsewhere in the RAD tooling (a fixed local step,
  rotated by a turn angle of `2*pi/n` at each step) all the way around a
  loop. When every cell shares one alpha this is an exact regular n-gon -
  the "ring closure residual" readout is a live numerical check of that,
  not a fixed display value. With per-cell actuation (below), cells
  generally differ, so the fixed-turn construction still walks all the way
  around but the *chord lengths* differ per joint - the residual then
  measures how inconsistent the per-cell commands are with a physically
  closed ring, which the coupling propagation is what keeps small in
  practice.
- **Drive**: the paper-supported dilation factor `alpha`, converted to the
  cell rotation `theta = 70*alpha - 60`, with a bidirectional backlash
  dead-zone (`f(x) = max(0, x-b) + min(x+b, 0)`) applied around `alpha = 1`.
- **Per-cell actuation and coupling propagation**: every cell is "free"
  (its alpha relaxes toward its neighbors' through the same backlash-gated
  `reluDeadzone` law Drive uses), an "actuator" (a directly commanded
  alpha), or "locked" (frozen at whatever alpha it had when locked). This
  is a discrete relaxation over the ring/cylinder's actual neighbor graph
  (circumferential + axial neighbors), not an independent slider per cell -
  actuating one cell tugs nearby free cells by a die-off distance, matching
  the paper's coupling model. With no actuators/locks this reduces exactly
  to the uniform-ring behavior. Each of the `m` rows solves its own ring
  from that row's own per-cell alphas, so rows can end up at different
  diameters (barrel/cone profiles) instead of always stacking as a true
  cylinder.
- **Target Diameter fit**: diameter isn't monotonic in alpha (it rises then
  falls as cells over-rotate past their most-open pose), so "Fit Alpha"
  exhaustively evaluates every alpha the slider can reach and keeps
  whichever produces the closest ring diameter, reporting the achieved
  diameter and residual honestly (including when a target is out of the
  achievable range) rather than solving a closed-form inverse. Fitting
  clears any actuators/locks first, since it targets one shared alpha.
- **Cylinder stacking**: `m` rings are repeated along the axle (`z`) axis at
  a configurable pitch. Axial row-to-row attachment is currently a rigid
  copy, not a solved joint - see Known limitations.
- **Material Restriction**: neighbor-only clearance checking (circumferential
  and axial neighbors only, not every cell pair) using the same disk/capsule
  primitive distance math as the row primitive.
- **Selection tools**: click a cell to inspect it (row/index, center,
  rotations, current alpha, role), then Frame Cell (recenter the camera on
  it), Isolate Cell (hide everything else), or edit its role/alpha directly
  from the Selected Cell panel. Focus hides the control panel for a clean
  view.
- **Alpha heatmap**: an optional toggle (Cross Colors panel) that colors
  each cell by its own current alpha (blue=contracted, gray=neutral,
  red=expanded) instead of the fixed row palette, so coupling propagation
  from an actuator and per-row diameter differences are visible at a
  glance.
- **Export OBJ**: serializes the actual rendered geometry of every visible
  cell (arms, hubs, pads, hole placeholders, in real millimeters) into a
  Wavefront OBJ file for reference/measurement in CAD software. This is a
  geometry snapshot, not a print-ready model - it does not perform boolean
  subtraction, so hole locations remain solid placeholder cylinders, same
  as on screen.
- **Timeline**: capture named keyframes (full state snapshots - alpha,
  backlash, ring/row counts, axial pitch, attachment sites, and every
  cell's role/alpha), jump between them, or play through them as discrete
  1.2s-interval steps (not interpolated). Save/Load Sequence JSON
  round-trips the whole keyframe list through a file.
- **Save / Load**: JSON export/import of the full configuration, including
  per-cell roles.

## Known limitations (intentional, not yet done)

- Only a single-pin joint per cell-to-cell connection is modeled, not the
  row primitive's mirrored double-pin closure. Generalizing double-pin
  closure to an arbitrary n-cell ring with independently-actuated cells is
  a real nonlinear loop-closure problem; it was left for later rather than
  shipped half-verified.
- Axial (row-to-row) attachment is a rigid z-offset copy, not a solved
  joint - there's no axial backlash or clearance-driven pitch yet.
- Per-cell actuation is a discrete backlash-gated relaxation, not a force/
  energy equilibrium solve - it's a reasonable discrete analogue of the
  paper's coupling law, not a calibrated mechanics model.
- Timeline playback jumps between keyframes rather than interpolating
  between them.
- No Jacobian sensitivity analysis or gradient-based inverse design (the
  Target Diameter fit is a direct exhaustive search over one shared alpha,
  not a general per-cell optimizer), unlike the much larger RAD digital
  workbench in `../Ahyan Virtual Simulation/web/`.

## Test

```powershell
npm install playwright pngjs
npx playwright install chromium   # once
node tests/validate_cylinder_tiling_page.js
```

The Playwright script loads the page in a real Chromium browser (desktop and
mobile viewports) and checks: the ring closure residual (both the uniform
and per-cell-actuated cases), the paper angle law, the backlash dead-zone,
neighbor clearance reporting, cell picking, Frame/Isolate/Focus behavior,
per-cell actuator/lock roles and their coupling propagation to a neighbor
(read through a minimal `window.__cylinderTilingDebug` hook rather than
guessing screen coordinates across camera angles), Target Diameter fitting
(including that the reported achieved diameter matches what's actually
rendered), Timeline capture/jump/play/save/load, and a Save/Load JSON
round-trip - then screenshots the rendered cells and checks their colors.
It also includes a regression test for a real display-scaling bug: on a
non-100%-scaled monitor, the canvas's on-screen size could drift or run
away across frames if not pinned to its container explicitly.
