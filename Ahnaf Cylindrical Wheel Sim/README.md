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

From Ahyan's V2 two-cell attachment (`V2/web/two-cell-attachment/`, read
only - nothing in that folder is modified), this project also borrows two
modeling concepts, reimplemented here in 3D:

- Sizing backlash from a physical pin in the hole (pin radius -> radial
  clearance `b` -> `b/L` -> in-plane dead zone `dphi = asin(b/L)`) and
  rendering those pins in the joints. Here the same pin also drives the
  out-of-plane joint tilt limit (Jacob's pin-in-hole contact relation).
- His "no discontinuous motion" drive envelope: step the drive in 1 deg
  increments and refuse a step where any cell jumps more than 4 mm or turns
  more than 5 deg. Here it runs as a continuous constraint solver on the
  whole cylinder, together with a new 3D collision check.

Earlier versions also adapted his flat 2D disk/capsule clearance checker and
his single-pin A-B chain walk; both have been replaced (see "Ring
construction" and "Collisions" below). Everything else - ring construction,
cylinder stacking, collisions, constraints, backlash coupling/propagation,
per-cell actuation, cell picking, timeline, target-diameter fitting, and the
UI - is new work, not adapted from that project's files.

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

- **Ring construction**: each cell's dilation twist is split symmetrically -
  lower cross at `-theta/2`, upper at `+theta/2` about the cell's normal -
  which puts its two circumferential joint pins (upper east -> next cell,
  lower west <- previous cell) level with each other, a chord
  `p = 2L*cos(theta/2)` apart (`L = 22.1 mm`). Twisting only the upper cross
  tilts that chord by `theta/2`, which on a cylinder sends the row off
  axially and left pinned pads ~10 mm apart. The pins are the vertices of a
  polygon inscribed in a circle about the axle; each cell lies flat on one
  face. The circumradius is solved so the chords close exactly
  (`R = p/(2 sin(pi/n))` for one shared alpha, bisection on the arc sum for
  per-cell alphas), so the ring always closes and the "ring closure
  residual" is ~0. Whether a pose is physically reachable is decided by the
  tilt check and the collision check instead. With opposite attachment sites
  (the default east/west) a second pad pair also meets at every joint (cell
  i's lower east with cell i+1's upper west, `2L*sin(theta/2)` along the
  axle): two pins per joint, the mirrored double-pin closure. Every ring is
  centered on the axle and rotated so cell 0's hub sits at polar angle 0, so
  rows stay coaxial and stack in columns.
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
- **Differential dilation (batch selection + shape presets)**: the whole
  structure does not have to dilate uniformly. Shift-click any number of
  cells in the 3D view to build a batch selection (shown with a green ring
  marker, independent of the single-cell inspector below), then set a role
  and alpha in the Per-Cell Actuation panel and apply it to every selected
  cell at once - this is the direct "select individual cells to expand/
  contract independently" workflow. For a faster demo, the three Shape
  Preset buttons (Barrel, Cone, Saddle) command a different alpha per row in
  one click. They pick alphas inside the current collision-free window so
  the constraint solver can actually reach them, and since the ring gets
  *smaller* as alpha rises in that window (more twist shortens the pin
  chord), the widest row gets the lowest alpha. At the defaults: Barrel
  117 / 132 / 117 mm, Cone 117 -> 126 -> 132 mm, Saddle 132 / 117 / 132 mm.
- **Target Diameter fit**: diameter isn't monotonic in alpha over the full
  slider range (the pin chord is longest at zero twist, alpha ~0.86), so
  "Fit Alpha" exhaustively evaluates every alpha the slider can reach and
  keeps whichever produces the closest ring diameter, reporting the achieved
  diameter and residual honestly rather than solving a closed-form inverse.
  With constraints on, it only considers alphas inside the collision-free
  window the pose is currently in (the pose can't cross a collision to reach
  another). Fitting clears any actuators/locks first, since it targets one
  shared alpha, and applies live as the slider is dragged.
- **Radial cell orientation**: each cell's thickness axis (top/bottom cross
  plane normal) points radially outward from the ring's own center, and its
  flat cross face lies tangent to the cylinder (spanning the circumferential
  and axial directions) - the same way the physical hardware's cells tile
  the wheel's surface, not flat-on to the axle like earlier versions. Each
  cell's frame comes straight from the ring construction: tangent along its
  polygon face, normal outward through the face, axial along the axle
  (`THREE.Matrix4.makeBasis`). A regression test checks the normal against
  the ring's actual geometry (`getRingGeometry`), not against the
  orientation code's own inputs - an earlier version passed every
  self-consistency check while rendering cells rotated like spokes.
- **Cylinder stacking**: `m` rings are repeated along the axle (`z`) axis at
  a configurable pitch. Rows are still a set spacing apart rather than
  pinned and solved - see Known limitations. A cell's north pads and the
  south pads of the cell above it are the axial joint (drawn with pins).
- **Pin alignment**: the circumferential readout is how far each joint's
  paired pad centers sit off a common pin axis (the bisector of the two
  cells' normals), read from the actual pad positions. It is ~0: their only
  separation is along the pin, by the stacked plates (`4 mm * cos(bend/2)`).
  A test checks this from the rendered meshes. The axial readout is the gap
  between a cell's north pad and the south pad of the cell above:
  `axialPitch - 2L*cos(theta/2)`-ish, 5.8 mm at the default pitch.
- **Collisions ("no fusing through other solids")**: every cross is modeled
  at its real 4 mm thickness as a hub disc, four pad discs and four
  half-arms, each extruded along the cell's own normal in its actual 3D
  pose. Penetration is estimated by sampling points across each part (rim,
  interior, three depths) and evaluating the other part's exact signed
  distance. Each cell is checked against its ring neighbors one and two
  over and the three nearest cells in the next row, culled by bounding
  spheres (~1-2 ms per check at the defaults). Joint pads - the pinned pad
  pairs at circumferential joints, and the north/south pads at axial joints
  - are exempt from each other, with the arms leading to them; how far those
  joints can bend is the backlash tilt check. Red dots mark overlap.
- **Constraints ("no discontinuous motion")**: with constraints on (the
  default), the pose never jumps to a command. It walks from the last valid
  pose toward the command in 1 deg steps of cell twist, collision-checking
  each step; a step that would overlap, or make any cell jump more than
  4 mm or turn more than 5 deg (Ahyan's V2 thresholds), is refused, and the
  pose settles at the contact point found by bisection ("held at contact").
  Leaving an already-overlapping pose is allowed as long as the overlap
  doesn't grow. Loading a file or Reset All jumps straight to the new setup.
  The panel shows the constraint state, the realized alpha against the
  command, and the collision-free range of one shared alpha. At the
  defaults (n = 10, pitch 50 mm) that is effective alpha 1.22-1.75 plus a
  separate 0.40-0.49 window across a colliding band that the pose can't
  cross; below ~1.22 neighbors' same-layer pads clash, above ~1.75 arms do.
  The default alpha is 1.4 (effective 1.30) so the page opens in a valid
  pose; the earlier 1.3 (effective 1.20) overlapped by ~0.5 mm.
- **Backlash tilt check**: with radially facing cells, each joint bolt points
  radially, so the bend between neighboring cells is a tilt of two plates on
  one bolt. The limit comes from the exact pin-in-hole contact relation
  `d + t*sin(theta) = D*cos(theta)`, i.e.
  `theta_max = acos(d / sqrt(D^2 + t^2)) - atan(t / D)` (~= `b/t` for
  `b = D - d << t`), doubled for two plates on a floating bolt. With the CAD
  hole `D = 3.4 mm`, plate `t = 4.0 mm` and the pin diameter `d` from the
  Backlash panel (default 3.0 mm, an *assumed* M3 bolt - not measured,
  confirm against the real parts), the joint limit is ~11.03 deg. Circumferential joints bend `360/n` deg (36 deg at
  `n = 10`), so a ring needs at least 33 cells to close through backlash
  alone; fewer means the hardware must take up the bend some other way
  (linkage pieces, or a dedicated hinge cell). Axial joints are checked
  against the profile slope between rows of different radius. Each is shown
  as within/exceeds backlash, alongside the minimum cell count. A thinner pin
  raises the limit (a 2.44 mm pin gives ~25 deg and 15 cells); a pin that
  fills the hole allows no tilt at all.
- **Physical pins**: a pin-diameter control (Backlash panel) sizes a pin in
  each 3.4 mm hole, following Ahyan's V2 pin-sized backlash. It reports the
  radial clearance `b`, `b/L` (L = 22.1 mm arm), the in-plane dead zone
  `dphi = asin(b/L)` and that same dead zone expressed in alpha (`dphi/70`,
  from the angle law), and renders pins through every hub and through both
  pad pairs of every circumferential and axial joint (hideable).
  The paper's design-space "normalized gap" slider is kept separate: it is
  the backlash the coupling model uses, while the pin is the hardware - the
  readouts show how far apart the two are (the default 0.10 gap is much
  larger than a 3.0 mm pin's ~0.007 in alpha).
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
- **Measurement line**: an optional CAD-style dimension line for row 0's
  diameter, offset below the assembly with extension lines back to the
  true diameter endpoints (drawing it straight through the ring's own
  center renders mostly occluded by the cell geometry, so it's offset
  clear of the part instead, the way an engineering drawing would).
- **Reset All**: restores every control to its HTML-declared default,
  clears actuators/locks, deselects the current cell, exits Isolate/Focus,
  stops any Timeline playback, and resets the camera - a clean-slate
  button distinct from the viewbar's own Reset (which only resets the
  camera).
- **Timeline**: capture named keyframes (full state snapshots - alpha,
  backlash, ring/row counts, axial pitch, attachment sites, and every
  cell's role/alpha), jump between them, or play through them as discrete
  1.2s-interval steps (not interpolated). Save/Load Sequence JSON
  round-trips the whole keyframe list through a file.
- **Save / Load**: JSON export/import of the full configuration, including
  per-cell roles.

## Known limitations (intentional, not yet done)

- The double-pin closure only arises for opposite attachment sites
  (east/west, north/south); other site pairs get a single pin per joint.
- Axial (row-to-row) spacing is a set pitch, not solved from the pin
  geometry: the north/south pads of stacked cells sit 5.8 mm apart along the
  axle at the default pitch rather than concentric, and their exemption
  from the collision check assumes they are the axial joint.
- Collisions are checked between near neighbors only (ring neighbors one
  and two over, and the three nearest cells in the next row), and
  penetration is sampled rather than solved exactly, so very shallow
  overlaps between samples can be missed (~0.05 mm tolerance).
- The constraint solver interpolates every cell toward the command together;
  when one part of the structure hits contact, the whole pose holds there,
  rather than letting unconstrained cells keep moving.
- Per-cell actuation is a discrete backlash-gated relaxation, not a force/
  energy equilibrium solve - it's a reasonable discrete analogue of the
  paper's coupling law, not a calibrated mechanics model.
- Timeline playback jumps between keyframes rather than interpolating
  between them.
- No Jacobian sensitivity analysis or gradient-based inverse design (the
  Target Diameter fit is a direct exhaustive search over one shared alpha,
  not a general per-cell optimizer), unlike the much larger RAD digital
  workbench in `../Ahyan Virtual Simulation/web/`.
- Constraints are kinematic (poses are checked, not simulated): there are
  no contact forces or dynamics, e.g. the Lankarani-Nikravesh contact model
  in Jacob's pin one-pager.

## Test

```powershell
npm install playwright pngjs
npx playwright install chromium   # once
node tests/validate_cylinder_tiling_page.js
```

The Playwright script loads the page in a real Chromium browser (desktop and
mobile viewports) and checks: exact ring closure (uniform and per-cell
actuated), the paper angle law, the backlash dead-zone, that every joint's
pinned pads share one pin axis (read from the rendered meshes), the axial
gap constant, the backlash tilt check and pin sizing, the constraints (the
default pose is collision-free, driving past either edge of the window holds
at contact without overlap, the pose can't jump to the separate clear
window, switching constraints off lets it overlap, and it can move back out),
that all three shape presets are reachable and produce their named profile,
cell picking, Frame/Isolate/Focus behavior,
per-cell actuator/lock roles and their coupling propagation to a neighbor
(read through a minimal `window.__cylinderTilingDebug` hook rather than
guessing screen coordinates across camera angles), shift-click batch
selection and applying a role/alpha to every selected cell at once, Target
Diameter fitting
(including that the reported achieved diameter matches what's actually
rendered), Timeline capture/jump/play/save/load, and a Save/Load JSON
round-trip - then screenshots the rendered cells and checks their colors.
It also includes a regression test for a real display-scaling bug: on a
non-100%-scaled monitor, the canvas's on-screen size could drift or run
away across frames if not pinned to its container explicitly. A second
regression test checks radial cell orientation against the ring's actual
center/centroid geometry (not just the orientation code's own inputs) -
the way an earlier version's bug passed every self-consistency check on
its own math while still rendering every cell rotated off-axis.
