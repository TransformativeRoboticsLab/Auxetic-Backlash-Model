# RAD Backlash Lattice Simulator

This repository contains a first-pass 3D demonstration simulator for Reconfigurable
Auxetic Devices (RADs) with backlash. The model is intentionally synthetic and
dimensionless: it is meant to expose the backlash-coupled lattice mechanics before
CAD calibration or target-surface inverse design.

## Run

```powershell
python -m rad_sim.demo_kinematic
python -m rad_sim.demo_spring_hinge
python -m rad_sim.interactive
```

The interactive module opens an `ipywidgets` UI when run inside Jupyter. Outside
Jupyter it falls back to a static Matplotlib demo.

## Browser UI

Open `web/index.html` directly in a browser. The page is a static
HTML/CSS/JavaScript CAD-like 3D simulator powered by a local vendored Three.js
asset at `web/vendor/three.min.js`, so the main simulator can load without a CDN
connection. It includes
orbit camera controls, click-to-select cells, visible rotating-square plates,
plate edge outlines, screw-head fasteners, hinges, corner pivot bosses, diagonal
braces, inter-cell linkage rods, backlash gaps, visible backlash limit stops,
actuator markers, membrane shaping, overlays, presets, measurement
labels, target-surface preview, a first target-fit seeding tool, an inspectable
inverse-design plan, smooth interpolated event timeline controls, JSON save/load,
and live renderer diagnostics for FPS, draw calls, and triangle count. Display
controls can switch the membrane between linear facets and smooth Catmull-Rom
interpolation, and adjust the number of surface samples used between cell
centers. The default cell visualization is now `Cell abstraction`: each cell is
drawn as an idealized outer square, rotating inner square, fixed datum bars,
rotating theta guide, corner pivots, and backlash clearance marker. Detailed
rods, stops, fasteners, and actuator hardware are suppressed in this mode so the
view reads as a CAD-style topology diagram rather than a speculative fabrication
model. The previous physical-looking linkage/fastener view is still available
as `Mechanism detail`, but it is intentionally treated as a visualization aid
rather than an exact CAD reconstruction. Switching to `Mechanism detail`
automatically enables the detailed display layers needed to inspect the
normalized mechanism: hinge pins, pivot bosses, diagonal braces, linkage rods,
fasteners, backlash stops, and active actuator markers. Switching back to
`Cell abstraction` restores the cleaner topology diagram defaults.
Saved `rad-sim.browser.v1` JSON includes both command inputs and derived cell
state. `cells.commandAlpha`, `cells.commandZ`, and `cells.locked` remain the
authoritative controls; `cells.alpha`, `cells.theta`, and `cells.z` are refreshed
from the current simulation before snapshots, event capture, and JSON export so
saved files can be inspected without rerunning the solver.
Vertical actuation is no longer treated as perfectly local. The browser model
computes a residual `zResidual` field from `cells.commandZ`: the selected cell
receives the commanded height, and neighboring cells receive a backlash-gated,
decaying fraction of that motion. The Lattice panel's `Z coupling` slider
controls this residual spillover; setting it to zero recovers the earlier
local-only vertical command behavior.
The viewport has explicit `Iso`, `Top`, `Front`, and `Side` camera presets with
`Z` treated as the vertical axis in both the scene math and the camera up vector,
plus scene `X/Y/Z` labels for orientation. The top preset uses the positive `Y`
axis as screen-up while looking down `Z`, which keeps the top-down view stable
without changing the simulator's Z-up convention.
The floor also includes a normalized-unit scale bar, so cell spacing and surface
deformation can be read against the same dimensionless units used by the model.
The lower viewport orientation HUD mirrors CAD camera conventions with colored
`X`, `Y`, and `Z` chips plus live view-axis, screen-up-axis, and projection
readouts. It updates after preset changes, orbit/pan moves, projection toggles,
and saved camera restores.
Left-drag rotates the orbit camera, Shift-drag or middle/right-drag pans the
camera target, and the mouse wheel zooms. Any manual camera movement switches
the view readout to a custom camera state so the active preset buttons do not
misrepresent the current orbit. `Reset` restores the isometric angle and clears
the panned/framed camera target back to the lattice origin. The `Focus` button
turns the browser into an inspection workspace by hiding the metric header and
right control panel; press `Controls` to restore the full editing UI.
The viewport projection toggle switches between perspective and orthographic
camera modes. Orthographic mode keeps dimensions visually parallel for CAD-style
inspection, and the selected projection is saved with the rest of the camera
pose in browser JSON.
Hovering over a cell highlights it in the 3D scene and reports the row/column in
the viewport status line; clicking still performs persistent cell selection.
`Frame Cell` recenters the orbit camera on the selected cell's current deformed
center, which is useful for inspecting a local backlash or actuator response
without losing the current viewing angle. Automatic lattice fitting preserves
panned or framed camera targets during animation; only `Reset` clears the target
back to the lattice origin.
`Isolate Cell` temporarily hides the rest of the lattice, membrane, target
surface, reference grid, inter-cell linkages, and broad vector overlays while
keeping the selected cell and its measurement guides visible. This makes the
abstract unit-cell topology easier to inspect before switching back to
`Show Lattice` for full-array behavior.
`Explode Cell` switches the selected cell into mechanism detail, isolates it,
and offsets the plates, hinges, pivots, braces, backlash stops, and active
actuator marker away from their assembled positions. The exploded view adds
leader-line callouts for the rotating plate, hinge pin, backlash gap, vertical
actuator, and alpha actuator when those parts are visible. Use `Assemble Cell`
or `Show Lattice` to restore the normal assembled view.
The top metric strip defaults to core inspection values only. Use `All Metrics`
when debugging inverse planning, actuator saturation, linkage strain, surface
slope, or renderer diagnostics; switch back to `Core Metrics` for a cleaner CAD
inspection view.
The viewport also has workspace mode buttons. `Inspect` gives a clean abstract
cell view for reading dilation and backlash; `Edit` shows selected actuators,
state coloring, pivots, and linkages for changing commands; `Analyze` enables
target/error-oriented overlays, target vectors, contours, normals, and the full
metric strip for inverse-design debugging.
On desktop, the app uses a fixed-height CAD-style workspace: the 3D lattice
stays visible in the left viewport while the right control panel scrolls inside
its own column. On narrower screens, the layout returns to normal page scrolling
so the controls remain reachable.
The right control panel uses collapsible sections. Core editing panels open by
default, while target fitting, inverse planning, timeline analysis, grid, and
persistence tools start collapsed so the first screen stays usable during CAD
inspection.
The lower-right viewport legend updates with the active overlay, naming the
displayed field and the meaning of the low/high color range so alpha, height,
target error, strain, sensitivity, and reachability views can be interpreted
without hunting through the controls.
The lower-left selected-cell HUD follows 3D cell picks and reports row/column,
alpha, rotation theta, height, current alpha/z commands, and whether the cell is
free, locked, or masked from inverse actuator placement. It is intentionally a
compact strip so it does not cover the active cell or block the lattice.
Surface contour guides can be enabled for membrane height, target height, or
target-error comparison; the guide lines sample the same interpolation field as
the rendered membrane so they stay aligned with the visible surface.
The Surface slope overlay computes finite-difference local tilt from the height
field. The top metrics report mean and maximum slope so surface smoothness and
sharp local bends can be inspected alongside height and target-error overlays.
The Show surface normals toggle draws sampled 3D normal vectors from the
membrane. The same finite-difference field reports mean and maximum normal tilt
in degrees, and the selected-cell readout includes local normal tilt so bends can
be inspected numerically as well as visually.
The display animation response slider controls how quickly the rendered lattice
eases toward a new actuation state, which makes slow expansion/contraction
inspection possible without changing the underlying kinematic solve.
The Linkage strain overlay colors the inter-cell connector rods and cells by the
current center-to-center stretch implied by the synthetic deformation field. The
top metrics report mean and maximum absolute linkage strain, and the selected
cell readout includes its worst adjacent linkage strain.
The reference lattice overlay draws the undeformed normalized cell-center grid
behind the mechanism. The Reference displacement overlay and selected-cell
readout compare the current deformed center position against that baseline, which
helps inspect expansion/contraction and out-of-plane actuation without changing
the solve. Displacement vectors can also be shown from reference centers toward
current centers; the renderer keeps the largest vectors plus the selected cell so
the direction field remains readable on dense grids.
The Show influence footprint toggle draws source-specific coupling guides around
the selected cell. Blue footprint rings/curves show the alpha/backlash influence
that would spread from the selected in-plane command; red rings/curves show the
vertical residual spillover from the selected z command. The footprint uses the
same dead-zone propagation as the solver, so increasing backlash visibly shrinks
the affected neighborhood.
The selected-cell dock also reports the same coupling numerically: `a next` and
`z next` are the immediate post-dead-zone neighbor transmissions, while `a reach`
and `z reach` count how many other cells receive nonzero selected-source
influence under the current backlash and coupling settings. Live badges below
the inspector say whether alpha and z are still in the free backlash gap or have
crossed into coupled neighbor transmission.

Actuator display can be filtered to the selected cell and accepted inverse-plan
cells, inverse-plan cells only, or all active cells. This keeps dense patterns
inspectable while preserving a full hardware view when needed.

The Target Surface panel exposes configurable actuator travel limits for vertical
motion, alpha contraction, and alpha expansion. Manual commands, seed fitting,
optimized fitting, inverse plan candidates, and sensitivity analysis all clamp to
the same limits. The Actuator saturation overlay and top-bar saturation metrics
show which cells are running near the current travel envelope; saturated slider
collars are colored with the active actuator material in the 3D view. Active
cells show both a vertical z-travel actuator and a horizontal in-plane alpha
actuator with its own rail, moving sleeve, command needle, and contraction/
expansion stops. Both actuator gauges now include physical tick marks for the
negative limit, midpoint, and positive limit. When measurements are enabled, the
selected cell also gets z and alpha travel brackets so the current command can
be compared directly against the configured travel envelope in the 3D scene.

Selected-cell alpha, vertical command, and lock controls update the 3D lattice
live. Use `Apply` when the current cell command should be recorded as a timeline
event; the live preview state can still be saved to JSON even before committing
an event. These high-use controls are docked next to the 3D viewport so cell
actuation does not require scrolling through the full sidebar. Use the dock's
`Hide` button to collapse it to a small handle when inspecting the lattice; the
selected-cell HUD remains visible while the controls are hidden. `Paint Clicks`
turns the current dock values into a direct 3D painting tool: each clicked cell
receives the active alpha command, vertical command, lock state, and inverse
actuator mask while keeping the lattice visible. The `Paint radius` slider
expands this into a local brush: radius `0` edits one cell, radius `1` edits the
clicked cell plus its four von Neumann neighbors, and radius `2` paints the next
Manhattan ring as well. While Paint Clicks is on, the 3D viewport previews the
affected brush cells under the cursor before the click is committed. The radius
is saved with browser JSON.
`Save JSON` also captures the current orbit camera pose, including preset/custom
mode, radius, spherical angles, and target point. Loading that JSON restores the
same inspection angle in addition to grid, actuation, target, overlay, and
timeline state.

Measurement labels can be switched between a compact all-values readout and
focused alpha, rotation theta, height, or backlash modes. Focused modes add
matching guide geometry, such as a theta arc or backlash clearance rays, around
the selected cell. The optional 3D text label is placed outside the lattice
footprint with a leader line, so inspection text does not sit on top of the cell
being actuated.

Current browser limitation: inverse design is an iterative greedy search, not a
full constrained optimizer. `Seed Fit` maps target height directly into commands;
`Optimize Fit` then searches cell commands to reduce target error while adding
simple penalties for actuator count and travel. `Analyze Plan` now accepts
actuators one at a time against the updated lattice state, estimates projected
score/error, records the accepted step history, and shows selectable planned
commands. `Apply Plan` commits the planned commands. The UI exposes RMS target
error, signed target-error bias, positive/negative residual range, iteration
count, actuator penalty, travel penalty, max actuators, candidate actuator count,
planner steps, projected score, recommended actuator count, and mean travel. The
target-error vectors can be toggled independently from the target surface so the
actual-to-target residual direction field can be inspected without hiding the
mechanism. The
inverse plan list can also preview a single recommended actuator: clicking a row
selects that cell, computes its standalone height contribution and target-error
change, and recolors the inverse overlay by contribution before the full plan is
applied. The plan-step preview slider can also scrub through the first N accepted
commands to inspect cumulative error reduction before applying the full plan. The
`Analyze Sensitivity` button runs a finite-difference command-response pass over
all unlocked and actuator-allowed cells. It records a per-cell sensitivity map,
ranked candidate list, controllable-cell count, and mean/max gain, then exposes a
dedicated sensitivity overlay. This is the first architecture hook for later
Jacobian-based inverse design. The `Build Jacobian` button now materializes that
next layer: it stores finite-difference height and dilation response columns for
positive/negative z commands and contraction/expansion alpha commands, records
per-cell reachability coverage, estimates a simple response conditioning ratio,
and exposes a Reachability overlay. `Solve Linear Fit` uses those columns in a
first greedy linearized residual projection and `Apply Linear Fit` commits the
resulting clamped commands. These response columns are still computed from the
current synthetic kinematic model rather than from a calibrated quasistatic
mechanism. The next solver should replace the greedy projection with a proper
nonlinear least-squares objective with actuator placement constraints, travel
limits, and continuation from previous equilibria.

The target surface panel includes a custom analytic expression mode for early
inverse-design experiments. Expressions are evaluated on normalized coordinates
`nr` and `nc` plus grid coordinates `r`, `c`, radial distance `d`, `amplitude`,
and `frequency`, with common math functions such as `sin`, `cos`, `exp`, and
`sqrt`. The browser applies a conservative identifier/character whitelist and
clamps custom target heights; this is for local exploration, not a sandboxed
multi-user expression service.

Selected cells also have an `Allow inverse actuator` mask. This is separate from
locking: a locked cell is mechanically fixed, while a disallowed cell can still
move but will not be chosen by `Seed Fit`, `Optimize Fit`, or `Analyze Plan`.
Masked cells are shown with a hatched cell-grid style. Bulk controls can allow
all actuator placements, block the rim, allow only the center region, or invert
the current placement mask.

Renderer note: the lattice cells, hinges, backlash rings, actuator hardware,
membrane, and target surface are now persistent Three.js objects. Lattice
hardware updates by transform; membrane and target surfaces update their vertex
buffers in place, and heatmap overlays reuse cached material buckets instead of
allocating per-cell materials on each update. Measurement labels, selected-cell
travel brackets, target-error rods, and sampled surface-normal vectors are still
rebuilt as lightweight overlays during state updates. Actuator
rails now render along the vertical axis with base/top travel stops and a moving
slider collar plus a side travel gauge scaled to the configured z-travel limit.
The in-plane alpha actuator renders as a normalized telescoping rail across the
cell; it indicates contraction/expansion command and saturation, but is still a
visual abstraction rather than a calibrated actuator assembly.
In `Cell abstraction` mode, the cell body is a simplified parametric diagram of
the unit cell topology, not a fabricated part model. It emphasizes dilation,
theta rotation, backlash clearance, selection, locks, and overlays. The fixed
datum bars show the cell frame; the rotating guide bars show the current theta
state. Backlash limit stops, corner pivot bosses, diagonal braces, inter-cell
linkage rods, plate edge outlines, screw-head fasteners, and actuator assemblies
are mechanism-detail cues only and are hidden from the abstract cell view unless
the display is switched to `Mechanism detail`. Linkage strain is a diagnostic
computed from neighbor center distances, not a force equilibrium or shared-node
compatibility constraint yet; these parts are not yet generated from a CAD
constraint model or shared-node constraint solve. The reference lattice is an
undeformed normalized grid, not a scanned or fabrication-calibrated datum.

Timeline note: events store full simulator snapshots. Playback now interpolates
compatible snapshots with an ease-in/ease-out transition before handing them to
the renderer, while incompatible grid-size changes still resolve as endpoint
state changes. The Timeline panel also supports named keyframe capture and
first/previous/next/last keyframe jumps, so important lattice poses can be
replayed as a compact experiment sequence while still preserving the underlying
event log. `Save Sequence` exports a separate `rad-sim.sequence.v1` JSON file
with the initial snapshot, ordered event frames, per-frame command summaries,
and a compact sequence summary. `Load Sequence` imports that same sequence
schema, rebuilds the replay timeline, and restores the simulator to the final
frame so the sequence can be inspected immediately. This is intended for
replay/analysis pipelines; sequence JSON still assumes the current synthetic
browser model and is not a fabrication recipe. `Save JSON` remains the full
editable simulator state export.
Inverse previews can be captured directly into the timeline as keyframes: a
single-candidate preview or cumulative plan-step preview is materialized as
actuator commands, recorded with preview error metadata, and can then be replayed
like any other keyframe.

Sequence analysis note: the Timeline panel can now analyze the current event
sequence, draw a compact chart of RMS target error, actuator saturation, and
command delta over time, and export a dense `rad-sim-sequence-metrics.csv` file.
Hovering over the chart previews the nearest stored frame's RMS error,
saturation, command delta, and active actuator count without changing the
simulator pose. Clicking the chart jumps the simulator to that frame using the
same smooth timeline transition as the event slider. `Best Error` jumps to the
lowest-RMS frame, while `Max Sat` jumps to the frame with the highest actuator
saturation so overdriven commands can be inspected quickly. The frame detail
also reports the worst target-residual cell, and `Worst Cell` selects that cell
in the 3D lattice with the target-error overlay enabled so local correction
candidates can be inspected directly. The analysis replays each stored snapshot
through the same browser
kinematic simulator used by the 3D view, then reports final RMS target error,
best-error frame, actuator saturation, and per-frame command delta in the panel.
The CSV emits one row per cell per frame with alpha, theta, height, target
height, residual, backlash influence, die-off, saturation, lock/actuator-mask
state, local link-strain diagnostic, slope, normal tilt, and the frame-level
worst-residual cell coordinates/value. These are
normalized simulation diagnostics; they are useful for inverse-design debugging
but are not yet measured forces, torques, calibrated actuator strokes, or
shared-node CAD constraints.

## Test

```powershell
python -m unittest discover -s tests
& "C:\Users\ahyan\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe" tests\validate_web_modules.js
```

The Node validation script executes the browser state/math/inverse modules in a
VM context. It also checks that `web/index.html` uses existing local assets with
the required script order and no external script/style dependencies, then starts
a temporary local HTTP server and fetches the entry page plus referenced assets.
It executes the vendored Three.js bundle and instantiates scene, camera, mesh,
material, raycaster, and vector primitives without creating a WebGL context. The
same script checks backlash dead-zone behavior, browser JSON roundtrip including
camera pose, preset variation, simulation metrics, and sequence export/import
without requiring Playwright or a browser binary.

## Model Layers

- `simulate_kinematic`: rotating-square auxetic cells with dead-zone backlash
  coupling, angle/dilation mapping, and synthetic out-of-plane height.
- `solve_spring_hinge`: reduced center-node spring-hinge quasistatic solver using
  SciPy optimization and penalty locks.
- `plot_lattice`: four 3D panels showing the lattice, actuated surface,
  complex-plane displacement, and rotation-angle surface.
