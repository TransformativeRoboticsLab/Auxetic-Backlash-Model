# RAD Backlash Lattice Simulator

This repository contains a first-pass 3D demonstration simulator for Reconfigurable
Auxetic Devices (RADs) with backlash. The model is intentionally synthetic and
dimensionless: it is meant to expose the backlash-coupled lattice mechanics before
CAD calibration while target-surface inverse design is still experimental.
Paper-supported formulas, prototype values, and current modeling assumptions are
tracked in `docs/research_grounding.md`; use that note before changing RAD cell
geometry or coupling laws. `model_provenance` in Python and
`RAD.modelProvenance()` in the browser expose the same paper-supported,
assumption, diagnostic, and calibration-gap classifications inside the simulator.

## Run

```powershell
python -m rad_sim.demo_kinematic
python -m rad_sim.demo_spring_hinge
python -m rad_sim.interactive
```

The interactive module opens an `ipywidgets` UI when run inside Jupyter. Outside
Jupyter it falls back to a static Matplotlib demo.

For numerical characterization of programmable discontinuities, the Python API
also exposes `characterize_single_cell`, `characterize_pair`, and
`characterize_cluster`. These helpers return response fields, vertical residuals,
die-off distances, and pairwise superposition error from the same simulator used
by the browser. `compare_physical_response`, `compare_physical_pair`, and
`compare_physical_cluster` run the same commands through the 3D spring-hinge
relaxation and report where physical center/height deltas deviate from the
kinematic backlash prediction.
For history-dependent programmable mechanics experiments, `local_actuation_event`,
`lock_event`, `release_event`, `clear_actuation_event`, `apply_event_sequence`,
`compare_event_order`, and `compare_sequence_order` model discontinuities as
state-to-state operators. A lock event commits the cell's current alpha state,
so lock/actuation ordering can be tested as a measurable noncommutative effect.
The browser mirrors this operator layer through `web/operators.js`. The selected
cell dock includes event lock, event release, and order-check controls; committed
locks store `cells.lockAlpha` so a locked cell can hold the dilation reached
before locking instead of always snapping back to the default state. Browser
order checks now include both the two-event actuation/lock commutator and a
three-event local sequence probe that reverses the sequence and swaps adjacent
events, giving a measured `maxOrderError` for path dependence.
`build_response_matrix` stacks those responses into alpha and height matrices for
controllability checks and inverse-design experiments; the matrix can also be
serialized as `rad-sim.response-matrix.v1` JSON for lab comparisons.
`sweep_response_atlas_parameters` reruns the calibration atlas across selected
backlash and pin-hole-clearance values, producing a
`rad-sim.response-atlas-sweep.v1` artifact for comparing locality, residual
vertical motion, and operator interaction before bench calibration.
`model_provenance` and `provenance_summary` return the current evidence ledger
for solver features, separating paper-supported equations from assumptions,
diagnostics, and missing calibration data.
`RADHardwareProfile`, `config_with_hardware_profile`, and
`hardware_profile_from_config` provide the first explicit bridge from normalized
simulator controls to measured real-cell dimensions. They track which key
hardware dimensions have been measured and which remain calibration gaps.
`response_decay_profile` fits the shellwise maximum response versus Manhattan
distance from active source cells, matching the browser Response Experiments
decay readout for notebook-side locality studies.
`characterize_pairwise_interactions` evaluates command pairs and builds
alpha/height residual matrices plus per-cell hotspot and degree maps that
identify which local operators are responsible for non-additive composition.
`diagnose_programmable_discontinuity` packages the same response fields into an
operator diagnostic: locality radius, reachable cells, response rank,
underactuated cells, fitted decay ratio/length, pairwise interaction graph,
superposition residual, and optional event-sequence order sensitivity. The
dead-zone law and rotating-square kinematics are paper-supported; the
superposition residual, fitted decay profile, pairwise interaction graph, and
sequence-order error are new diagnostics for detecting when composed operators
stop behaving additively or commutatively because of backlash thresholds, locks,
or saturation. `programmable_discontinuity_report` and
`export_programmable_discontinuity_report_json` serialize that diagnostic as a
`rad-sim.programmable-discontinuity-report.v1` artifact with explicit
paper-supported assumptions, simulator-introduced diagnostics, locality,
reachability, composition, response-matrix, pairwise, optional spring-hinge
physical-validation, and sequence-order fields.
`solve_inverse_design` uses that response matrix in a bounded damped least-squares
fit from target alpha/height fields to candidate actuator commands. This is the
first Python-backed inverse layer for surface-shaping studies; it is linearized
around the unactuated baseline and should be treated as a command proposal, not
as a calibrated nonlinear optimizer. Its result also carries unweighted alpha and
height residual fields, reachable/underactuated masks, and command-saturation
counts so poor target fits can be separated from mechanically unreachable or
travel-limited regions.
`validate_inverse_design_physical` runs a solved inverse command state through
the 3D spring-hinge model and reports physical target residuals plus kinematic
versus physical center/height disagreement. This is a validation layer for
candidate commands, not a physical inverse optimizer yet.

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
The `Paper RAD cell` display mode sits between those two views. It shows the
paper-grounded rotating-square structure as two concentric square parts with
four joints, pin bosses, and hole-clearance rings driven by the normalized
pin/hole controls. It is still normalized geometry, but it reflects the RAD
papers' stated unit-cell structure more directly than the general abstraction.
The Python API mirrors that with `build_paper_rad_cell_geometry`, which returns
a normalized single-cell description with two concentric square parts, four
joints per part, lock sites, alpha/z actuator axes, backlash gap, and vertical
free play. `build_paper_rad_lattice_geometry` lifts the same cell model across a
full simulated lattice and records inter-cell connector spans, neighbor alpha
jumps, height jumps, lock states, and actuator commands for pair, cluster, and
full-lattice inspection. The embedded reference values come from the RAD
preprint text: 35 mm prototype side length, normalized backlash `b = 0.1`, and
target effective Poisson ratio `-0.4`; exact CAD thicknesses and tolerances
remain configurable until measured from the physical parts.
`calibrate_paper_rad_config` converts normalized model lengths into paper
prototype millimeters. With the default reference, one model cell side maps to
35 mm, the paper `b = 0.1` backlash maps to 3.5 mm, and the extracted 0.1 mm
hole fabrication tolerance maps back into normalized model units. Current pin
and hole radii are still configured estimates, but their clearance can now be
reported in millimeters for comparison against measured parts.
The browser Lattice panel exposes the same scale bridge through `Paper side mm`
and `Hole tol mm` controls. Its live readouts report configured backlash,
pin-hole free play, and tolerance in millimeters/model units so abstract
simulator values can be compared against prototype measurements without
changing the normalized solver. Browser JSON now also carries a
`grid.hardwareProfile` record for measured dimensions such as pin radius, hole
radius, plate thickness, joint stack height, and boss radius; the Lattice panel
now includes a compact measured-profile editor, reports how many of those
measurements are present, and lists which are still missing. `Apply Measured
Dims` deliberately maps measured backlash, pin radius, and hole radius back into
the normalized backlash/free-play controls; unfilled fields remain calibration
gaps instead of silently changing the solver. The Display panel's `Calibrated
RAD cell` mode uses the profile directly for visible plate thickness, pin/hole
radii, and stack offset without requiring those values to alter the solver.
`calibration_readiness` in Python and `RAD.calibrationReadiness()` in the browser
classify the current profile as paper-scale, partial-measured,
visual-calibrated, or mesh-calibrated, while separately listing the physical
solver gaps that still require stiffness, actuator, friction/contact, and
response measurements. `calibration_measurement_plan` and the browser `Save
Calibration Plan` action turn those gaps into an ordered JSON checklist of the
next geometry and solver measurements needed before treating a RAD cell as
physically calibrated. `build_calibration_experiment_protocol` and the browser
Response Experiments `Save Protocol` action generate a repeatable single-cell,
pair-cell, cluster, and locked-cell test protocol for collecting those response
measurements against the same simulator coordinates. The companion results
template and comparison helpers accept optional measured alpha, height, center
displacement, slip, and actuator-force values, then report error against the
current simulator without treating the measurement set as a calibrated solver.
`build_response_atlas` and `export_response_atlas_json` run the same protocol
through the simulator as a compact `rad-sim.response-atlas.v1` artifact, with
optional spring-hinge physical-preview metrics, so backlash and clearance
settings can be compared before bench data exists. The browser `Save Atlas`
action exports the same type of response atlas from the current 3D simulator
state.
The browser Response Experiments panel can export a fillable results template,
load the completed JSON, and display compact measured-cell, missing-observation,
alpha RMSE, height RMSE, maximum combined per-cell error, and worst-step readouts.
Loading completed results also switches the Display panel to the `Calibration
error` overlay, which colors the 3D lattice by the averaged per-cell alpha/height
measurement disagreement. `Select Cal Error` moves the 3D selection to the cell
with the largest imported measurement disagreement; when the `Calibration
residual` overlay is active, the same button selects the post-fit residual
hotspot instead. `Next Cal Error` cycles through the ranked raw or residual
mismatch list, depending on the active overlay. This keeps mismatch inspection
in the viewport without placing labels over the lattice. `Save Comparison` exports a
JSON report containing the imported comparison, per-cell error field, worst-cell
record, signed alpha/height bias, measured-versus-simulated gain/bias fit, and
hardware/profile metadata for lab notes or GitHub review. The fit is diagnostic:
it estimates whether measured alpha or height response needs a scale or offset
correction, but it does not automatically mutate the solver. The Display panel's
`Calibration residual` overlay colors the remaining per-cell error after that
global gain/bias fit, making local mechanism mismatch easier to distinguish from
uniform scale or offset error.
For external inspection, `build_paper_rad_lattice_mesh` converts the normalized
paper RAD lattice into extruded plate, pin, and connector mesh components, and
`export_paper_rad_mesh_obj` serializes that mesh to OBJ text. This is meant for
CAD-style review and downstream tooling; it is not yet a fabrication export.
The browser OBJ export now records the active calibration profile in the header
and uses measured profile defaults for plate thickness, pin radius, stack height,
and connector width when available.
The browser Persistence panel also has `Save OBJ`, which exports the current
visible simulation state as a normalized paper-RAD OBJ mesh.
The Display panel can switch the visible solver between `Kinematic` and
`3D spring preview`. Kinematic mode shows the direct backlash propagation model;
spring preview runs a lightweight browser relaxation over the cell-center graph,
then recomputes height, slope, and link-strain diagnostics for a more physical
inspection view without replacing the Python spring-hinge solver. The
`Model disagreement` overlay uses that same spring preview to color cells by
per-cell height disagreement from the kinematic prediction, so suspicious
regions can be inspected directly in the 3D lattice.
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
The vertical dead-zone is now tied to normalized pin-hole free play: `hole radius
- pin radius`. This follows the RAD papers' treatment of backlash as joint
clearance and keeps the user-specified vertical residual assumption explicit
until calibrated physical pin and hole dimensions are measured.
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

The Response Experiments panel characterizes the current selected cell, selected
pair, local cluster, or active lattice command set. It reports active command
sources, total responding cells, alpha and z reach, alpha/z die-off, paper-scale
pin-hole clearance, local response rank, underactuated cells, and a
superposition residual that compares the combined response against the sum of
isolated source responses. It also reports a bounded pairwise interaction scan
that identifies how many source pairs create non-additive alpha/height residuals
and the largest pairwise interaction error. Running the characterization switches
the viewport to the Operator interaction overlay, which colors command-source
cells by the strongest non-additive pair residual found in that scan; `Select
Hotspot` moves selection to the strongest cell without adding labels over the
lattice. The same readout reports max non-additive degree, counting how many
non-additive pair edges touch the strongest cells. It also compares the same
response against the browser 3D spring preview and reports height/center
disagreement from the kinematic backlash prediction. The reported decay ratio
and decay length come from a simulator-side log-linear fit to shellwise maximum
response versus Manhattan distance from the active sources; this is a locality
diagnostic, not a paper-derived material law. This is an early numerical probe for
programmable-discontinuity behavior: nonzero residual indicates that backlash,
locks, or saturation are making the local operators interact non-additively.

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
and exposes a Reachability overlay. It also compares the current target residual
against finite response columns, counts target cells with requested height motion
outside the reachable set, and exposes an `Underactuated target` overlay plus
`Select Under Target` for the worst unreachable target cell. `Solve Linear Fit`
uses those columns in a first greedy linearized residual projection and `Apply
Linear Fit` commits the resulting clamped commands. These response columns are
still computed from the current synthetic kinematic model rather than from a
calibrated quasistatic mechanism. `Validate Physical` checks the current linear fit or analyzed plan
against the browser spring-preview relaxation and reports physical target error
plus kinematic/physical center-height disagreement. `Save Inverse Report`
exports a `rad-sim.inverse-design-report.v1` JSON artifact with the current
target, commands, plan, Jacobian/reachability diagnostics, linear fit, and any
physical validation. The next solver should
replace the greedy projection with a proper nonlinear least-squares objective
with actuator placement constraints, travel limits, and continuation from
previous equilibria.

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
- `build_paper_rad_cell_geometry`: normalized single-cell RAD geometry with
  two concentric square parts, four joints per part, pin/hole clearance,
  lock sites, actuator axes, and paper reference metadata.
- `calibrate_paper_rad_config`: paper-derived prototype scale conversion for
  side length, normalized backlash, pin/hole clearance, and hole fabrication
  tolerance in millimeters and normalized model units.
- `build_paper_rad_lattice_geometry`: full-lattice paper RAD geometry records
  derived from a `LatticeState` or `SimulationResult`, including inter-cell
  connector spans and neighbor alpha/height jumps.
- `build_paper_rad_lattice_mesh` / `export_paper_rad_mesh_obj`: normalized
  CAD-style triangle mesh and OBJ export for paper RAD plates, pins, and
  connectors.
- `solve_spring_hinge`: reduced center-node spring-hinge quasistatic solver using
  SciPy optimization and penalty locks.
- `solve_spring_hinge_3d`: out-of-plane spring-hinge relaxation over 3D center
  nodes, so vertical actuation and pin-clearance residual height can participate
  in the same axial/hinge/penalty energy model.
- `local_actuation_event` / `lock_event` / `release_event` /
  `compare_event_order` / `compare_sequence_order`: executable programmable-
  discontinuity operators for sequence studies, lock invariance, locality
  checks, adjacent-swap sensitivity, and event-order noncommutativity.
- `diagnose_programmable_discontinuity`: Python operator diagnostic for locality,
  reachable sets, response rank, underactuated regions, shellwise response decay,
  pairwise interaction matrices/hotspot/degree maps, and additive versus
  non-additive operator composition. When given an `event_sequence`, it also
  reports reversal/adjacent-swap order sensitivity.
- `programmable_discontinuity_report` /
  `export_programmable_discontinuity_report_json`: versioned JSON-ready research
  artifact for the programmable-discontinuity diagnostic, including the
  paper-supported dead-zone and alpha-theta assumptions plus simulator
  diagnostics for locality, reachability, nonadditivity, pairwise interactions,
  event-order sensitivity, and optional `include_physical=True` spring-hinge
  model-disagreement evidence.
- `characterize_pairwise_interactions`: Python pairwise command interaction graph
  for identifying which actuation operators create non-additive residuals and
  mapping each command cell's strongest pair residual and non-additive degree.
- `response_decay_profile`: Python log-linear shell-max locality diagnostic that
  mirrors the browser Response Experiments decay ratio/length readout.
- `solve_inverse_design`: bounded damped least-squares inverse design using
  finite response columns for alpha/height target fields, with residual fields,
  reachable/underactuated masks, and command-saturation diagnostics.
- `inverse_design_report` / `export_inverse_design_report_json`: Python inverse
  design report artifact containing response-matrix diagnostics, residual fields,
  reachability masks, commands, and optional physical validation.
- `validate_inverse_design_physical`: spring-hinge validation of a solved inverse
  command set, reporting physical target error and kinematic/physical model
  disagreement before a command set is treated as mechanically credible.
- `web/operators.js`: browser-side version of the event operators with committed
  `lockAlpha` state, selected-cell order diagnostics, and
  `compareSequenceOrder` adjacent-swap/reversal sensitivity checks in the
  viewport dock.
- `web/math.js` `paperRadCalibration`: browser-side paper scale conversion for
  backlash, pin-hole clearance, and fabrication tolerance readouts.
- `model_provenance` / `web/provenance.js`: shared evidence ledger that marks
  formulas and simulator layers as paper-supported, assumptions, diagnostics, or
  calibration gaps.
- `RADHardwareProfile` / `web/math.js` `calibrationProfileSummary`: measured
  hardware profile scaffolding that tracks pin, hole, plate, stack, and boss
  dimensions before they are allowed to replace normalized model assumptions.
  The browser editor saves those values in `grid.hardwareProfile` and only
  applies measured backlash/radii to normalized controls on explicit command.
- `calibration_readiness` / `web/math.js` `calibrationReadiness`: conservative
  readiness gate for calibrated geometry versus still-uncalibrated physical
  solver behavior.
- `calibration_measurement_plan` / `web/math.js` `calibrationMeasurementPlan`:
  ordered geometry and solver measurement tasks that show what must be measured
  next before the visual CAD profile can become a calibrated physical model.
- `build_calibration_experiment_protocol` / `web/analysis.js`
  `calibrationExperimentProtocol`: repeatable single, pair, cluster, and lock
  response protocol exports for matching simulator response fields against
  physical bench measurements.
- `build_response_atlas` / `export_response_atlas_json`: Python
  `rad-sim.response-atlas.v1` artifact that summarizes simulator responses for
  the calibration protocol's single, pair, cluster, and lock cases, including
  observation-cell deltas, die-off, superposition residuals, and optional
  spring-hinge model-disagreement metrics.
- `sweep_response_atlas_parameters` / `export_response_atlas_sweep_json`:
  Python `rad-sim.response-atlas-sweep.v1` artifact that compares atlas
  summaries over backlash and pin-hole-clearance settings for locality,
  residual vertical motion, and non-additive operator interaction studies.
- `web/analysis.js` `responseAtlas` / `exportResponseAtlas`: browser-side
  `rad-sim.response-atlas.v1` export behind `Save Atlas`, using the current
  protocol, kinematic response fields, and browser spring-preview validation.
- `web/analysis.js` `responseAtlasSweep` / `exportResponseAtlasSweep`:
  browser-side `rad-sim.response-atlas-sweep.v1` export behind `Save Sweep`.
  `Run Sweep` computes the same artifact into the Response Experiments readout
  so backlash/clearance trends can be inspected without leaving the 3D UI.
- `calibration_experiment_results_template` /
  `compare_calibration_experiment_measurements` and the browser matching
  helpers: fillable bench-results schema plus simulator-vs-measurement error
  report for the calibration protocol. The browser can import a completed
  results JSON and summarize the comparison in the Response Experiments panel.
  `calibration_experiment_comparison_report` and the browser `Save Comparison`
  action add per-cell raw/residual fields plus measured-vs-simulated gain and
  bias diagnostics without mutating solver parameters. The report includes
  ranked raw-error and post-fit residual cell lists so repeated local mismatch
  candidates can be inspected after a bench run.
- `web/renderer.js` `calibratedRad` mode and `web/mesh_export.js`
  `calibratedMeshDimensions`: profile-aware visual/export layer that uses
  measured dimensions for CAD-like inspection while keeping solver assumptions
  explicit.
- `web/analysis.js` `characterizeLocalResponse`: browser-side single, pair,
  cluster, and active-lattice response experiment metrics with superposition
  residuals, pairwise operator-interaction counts and hotspot maps, local
  response rank, underactuated-cell counts, browser spring-preview disagreement
  metrics, and log-linear locality/decay estimates for
  programmable-discontinuity studies.
- `web/analysis.js` `buildResponseMatrix` / `exportResponseMatrix`: browser-side
  alpha and height response-matrix artifact for the selected single, pair,
  cluster, or active-lattice actuator set.
- `web/analysis.js` `programmableDiscontinuityReport` /
  `exportProgrammableDiscontinuityReport`: browser-side
  `rad-sim.programmable-discontinuity-report.v1` export behind `Save Framework
  Report`, carrying the same paper-supported assumptions and simulator
  diagnostics used by the Python framework report, plus the browser spring-preview
  physical-validation summary.
- `web/inverse.js` `validateInversePlanPhysical`: browser-side physical
  validation of analyzed inverse plans and linear fits against spring-preview
  relaxation, including physical target residual and model-disagreement metrics.
- `web/inverse.js` `buildResponseJacobian`: browser-side target reachability and
  underactuated-height diagnostics for finite response columns.
- `web/inverse.js` `inverseDesignReport` / `exportInverseDesignReport`: browser
  inverse-design report artifact for target, plan, Jacobian, linear fit, current
  commands, reachability limits, and physical validation state.
- `web/physics.js` `simulatePhysicalRelaxation`: browser spring-preview
  relaxation that now also exposes per-cell `modelErrorHeight` and
  `modelErrorCenter` fields against the kinematic state for model-disagreement
  overlays.
- `compare_physical_response` / `compare_physical_pair` /
  `compare_physical_cluster`: single, pair, and cluster diagnostics that quantify
  the deviation between kinematic backlash propagation and the 3D spring-hinge
  physical response.
- `plot_lattice`: four 3D panels showing the lattice, actuated surface,
  complex-plane displacement, and rotation-angle surface.
