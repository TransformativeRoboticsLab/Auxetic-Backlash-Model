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
- The paper-supported `theta = 70*alpha - 60` angle law documented in
  `../Ahyan Virtual Simulation/docs/research_grounding.md`.

The cell dimensions come from Jacob's CAD of the real part ("RADs unit
cell.stl", measured directly): two 4-arm crosses, each a 4 mm plate (0-4 and
4-8 mm), joined at the hub by a countersunk M3 bolt (3.0 mm major diameter,
confirming the default pin); hole-to-hub distance L = 22.451 mm, pad radius
4.0 mm, hole 3.4 mm, hub radius 6.375 mm, arm width 2.5 mm. (Earlier versions
used 22.1 / 4.9 / 4.6 / 5.4 mm from Ahyan's research-grounding numbers.)

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
  `p = 2L*cos(theta/2)` apart (`L = 22.451 mm`). Twisting only the upper cross
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
- **Per-cell actuation and coupling through bolted joints**: every cell is
  "free", an "actuator" (a directly commanded alpha), or "locked" (frozen at
  whatever alpha it had when locked). Both stacked pad pairs of every joint
  are bolted, and the second pin's two holes slide apart when neighbors
  twist differently: by `2L*|sin(t_i/2) - sin(t_j/2)|` around a ring and
  `L*|...|` between pinned rows. A bolt allows at most the hole clearance
  `D - d` (0.4 mm for a 3.0 mm pin, about 1 deg of twist between neighbors).
  So free cells settle as close to the shared drive as those limits allow:
  actuating one cell drags its neighbors, fading with distance, and the pin
  size sets the die-off - e.g. one cell commanded to 1.60 in a 1.30 ring
  pulls the far side to 1.53 with a 3.0 mm pin, 1.40 with a 2.3 mm pin. (The
  limits keep 10% of the clearance in hand for small geometric extras.) The
  solver also refuses any step that would push a bolted pad pair further off
  its pin than the clearance ("a bolted pin would bind"). Earlier versions
  coupled cells through the paper's abstract backlash slider (0.10 alpha,
  ~7 deg), which let neighbors differ far more than a bolted joint can and
  visibly split the second pin's holes; that slider now only sets the
  drive's dead zone. With no actuators/locks this reduces exactly to the
  uniform-ring behavior. Each of the `m` rows solves its own ring from that
  row's own per-cell alphas.
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
  chord), the widest row gets the lowest alpha. With rows bolted together,
  adjacent rows can only differ by what the pins allow, so the profile is
  scaled down to the largest the pins permit: at the snug default 3.0 mm
  pin that is under 1 mm of bulge (e.g. Barrel 128.1 / 128.9 / 128.1 mm) - a
  real consequence of double-bolted rows, not a display limit. A thinner pin
  allows more, and unpinned rows (separate rings) allow the full window if
  they have room - Barrel 115.9 / 135.7 / 115.9 mm at a 60 mm pitch; at
  50 mm the rows' pads meet partway and the move holds early.
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
- **Cylinder stacking**: `m` rings are stacked along the axle (`z`). By
  default the rows are pinned together like the cells in a ring: a cell's
  north pads meet the south pads of the cell above, so the row spacing is
  set by the pins - `2L*cos(theta/2)` for a uniform ring (43.3 mm at the
  default) - and changes as the cells twist. Unchecking "Pin rows together"
  hands the spacing to the axial pitch slider instead; the rows are then
  separate rings with nothing joining them (no axial pins, and a large pitch
  simply leaves them floating apart). Rows at different twists can't have
  their north/south pads meet exactly - each pad sits `L*sin(theta/2)` off
  its hub along the ring - so e.g. the Barrel preset leaves ~6 mm for its
  axial joints to absorb, shown in the axial pin readout.
- **Pin alignment**: the circumferential readout is how far each joint's
  paired pad centers sit off a common pin axis (the bisector of the two
  cells' normals), read from the actual pad positions. It is ~0: their only
  separation is along the pin, by the stacked plates (`4 mm * cos(bend/2)`).
  A test checks this from the rendered meshes. The axial readout is the same
  measure for a cell's north pads against the south pads of the cell above:
  ~0 for pinned rows sharing one twist, the real gap for unpinned rows.
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
  4 mm or turn more than 5 deg (Ahyan's V2 thresholds), is refused. When a
  whole-structure step is refused, each row and then (in rows whose cells
  want different things) each cell tries on its own, shrinking the step
  down to 1/8 deg, so parts that can move keep moving and blocked parts
  creep up to their contact point and hold there ("held at contact"). Only
  rows a move touches are re-checked. Leaving an already-overlapping pose is
  allowed as long as the overlap doesn't grow. Loading a file or Reset All
  jumps straight to the new setup.
- **Reachable range and the drive limit**: the collision-free range of one
  shared alpha is swept per setting. At the defaults that is effective alpha
  1.15-1.82 (theta 20.5-67 deg), plus a separate 0.40-0.56 window no
  slider command can reach; below ~1.15 - which includes every negative
  theta - neighbors' same-layer pads clash at the joints, above ~1.82 parts
  of neighboring cells do. Like Ahyan's V2
  clamping its drive to the nearest feasible value, the shared drive is
  limited to that range, 0.02 inside its edges, so a command past it doesn't
  jam the whole structure against its stops: the Drive panel says the
  command is out of reach and where the pose actually is, and a green strip
  under the slider shows the reachable commands with a hollow marker at the
  realized pose. Individual cell commands aren't limited this way; they
  hold at contact instead. The default alpha is 1.4 (effective 1.30) so the
  page opens in a valid pose.
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
  radial clearance `b`, `b/L` (L = 22.451 mm arm), the in-plane dead zone
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
- Pinned row spacing is the average over a row's cells, and rows stay
  vertical; rows at different twists leave their north/south pads offset
  (reported, not resolved by tilting the cells).
- Small twists are out of reach: below about 20.5 deg each joint has four
  pads near one spot on two layers, and neighbors' same-layer pads overlap
  (they touch exactly when 2L*sin(theta/2) = 2 x 4.0 mm). The paper's
  reference state alpha = 1 (theta = 10 deg) falls in that band - but
  Jacob's CAD shows the real part drawn with its two crosses 36.9 deg apart
  (one on the diagonals, the other at 8.1 deg), well inside the reachable
  window. So the angle law's theta is probably not the geometric twist
  between the crosses (if the CAD is the alpha = 1 pose, it is offset by
  ~27 deg); this model uses theta as the geometric twist. To confirm with
  Jacob: which alpha the CAD pose is, and how theta is measured.
- Collisions are checked between near neighbors only (ring neighbors one
  and two over, and the three nearest cells in the next row), and
  penetration is sampled rather than solved exactly, so very shallow
  overlaps between samples can be missed (~0.05 mm tolerance).
- Coupling is kinematic: free cells take the pose closest to the drive that
  the bolt clearances allow, not a force/energy equilibrium - real parts
  have friction and compliance on top of the clearance.
- Any cell commanded past the reachable range is limited to just inside it
  before solving (like the drive), so an over-commanded actuator reports
  "limited" rather than pressing against contact.
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
