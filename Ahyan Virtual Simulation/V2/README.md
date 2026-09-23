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
python -m rad_sim.export_calibration_bench_packet --out outputs/calibration_bench_packet
python -m rad_sim.compare_calibration_bench_packet --fit outputs/calibration_bench_packet/calibration_fit_results_template.json --holdout outputs/calibration_bench_packet/calibration_holdout_results_template.json --out outputs/calibration_bench_execution_validation
python -m rad_sim.export_vertical_load_bench_packet --out outputs/vertical_load_bench_packet
python -m rad_sim.compare_vertical_load_measurements --input outputs/vertical_load_bench_packet/vertical_load_energy_measurement_template.json --out outputs/vertical_load_measurement_comparison
python -m rad_sim.export_two_cell_bench_packet --out outputs/two_cell_bench_packet
python -m rad_sim.run_two_cell_external_fidelity_mjcf --index outputs/two_cell_bench_packet/two_cell_external_fidelity_matrix_mjcf_index.json --out outputs/two_cell_external_fidelity_mjcf_run
python -m rad_sim.compare_two_cell_measurements --input outputs/two_cell_bench_packet/two_cell_measurement_template.csv --out outputs/two_cell_measurement_comparison
python -m rad_sim.compare_two_cell_fidelity_matrix_measurements --input outputs/two_cell_bench_packet/two_cell_fidelity_matrix_measurement_template.csv --out outputs/two_cell_fidelity_matrix_measurement_comparison
python -m rad_sim.compare_two_cell_external_results --input outputs/two_cell_bench_packet/two_cell_external_results_template.csv --out outputs/two_cell_external_comparison
python -m rad_sim.cad_mesh_intake assets/cad/RADs_unit_cell.step --out outputs/cad_mesh_audit.json --pin-radius-mm 1.36 --hole-radius-mm 1.70 --boss-radius-mm 3.10 --backlash-mm 0.34
```

Optional external-engine execution uses the `physics` extra:

```powershell
python -m pip install -e ".[physics]"
python -m rad_sim.run_two_cell_external_fidelity_mjcf --index outputs/two_cell_bench_packet/two_cell_external_fidelity_matrix_mjcf_index.json --out outputs/two_cell_external_fidelity_mjcf_stratified_real --steps 40 --case-ids FM_000,FM_001,FM_004,FM_005,FM_006,FM_007,FM_014,FM_021,FM_028,FM_056,FM_084,FM_168 --compare
python -m rad_sim.fit_two_cell_external_fidelity_correction --input outputs/two_cell_external_fidelity_mjcf_full_real/two_cell_external_fidelity_matrix_mjcf_measurements.csv --out outputs/two_cell_external_fidelity_mjcf_full_real_correction --holdout-stride 5
python -m rad_sim.apply_two_cell_external_fidelity_correction --profile outputs/two_cell_external_fidelity_mjcf_full_real_correction/two_cell_external_fidelity_correction_profile.json --measurements outputs/two_cell_external_fidelity_mjcf_full_real/two_cell_external_fidelity_matrix_mjcf_measurements.csv --out outputs/two_cell_external_fidelity_mjcf_full_real_corrected_preview
python -m rad_sim.export_two_cell_external_fidelity_web_data
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
`group_actuation_event`, `lock_event`, `release_event`, `clear_actuation_event`,
`remove_cell_event`, `restore_cell_event`, `apply_event_sequence`,
`compare_event_order`, and `compare_sequence_order` model discontinuities as
state-to-state operators. A lock event commits the cell's current alpha state,
so lock/actuation ordering can be tested as a measurable noncommutative effect;
remove-cell events delete a cell from the propagation graph until restored.
The browser mirrors this operator layer through `web/operators.js`. The selected
cell dock includes event lock, event release, event remove, event restore, and
order-check controls; committed locks store `cells.lockAlpha` so a locked cell
can hold the dilation reached before locking instead of always snapping back to
the default state. Browser order checks now include both the two-event
actuation/lock commutator and a three-event local sequence probe that reverses
the sequence and swaps adjacent events, giving a measured `maxOrderError` for
path dependence.
`build_response_matrix` stacks those responses into alpha and height matrices for
controllability checks and inverse-design experiments; the matrix can also be
serialized as `rad-sim.response-matrix.v1` JSON for lab comparisons.
`lattice_topology_diagnostic` reports present/removed masks, connected
component labels and sizes, active edges, and deleted edges after cell removal;
the browser mirrors this as `RAD.topologyDiagnostics` and exposes component
counts through simulation metrics.
`compare_removed_topology_reachability` runs a matched intact-versus-removed
response-matrix experiment for a selected removal set, reporting reach loss,
rank loss, component-count change, and deleted-edge change without conflating
topology loss with actuators placed on removed cells.
`build_topology_experiment_report` generalizes that into named topology
scenarios with locks, removals, group commands, signed vertical reach, component
blocking, component-wise response rank, die-off, and optional target residuals. The browser
mirrors this as `RAD.topologyExperimentReport` /
`RAD.exportTopologyExperimentReport` for JSON experiment snapshots.
`sweep_response_atlas_parameters` reruns the calibration atlas across selected
backlash and pin-hole-clearance values, producing a
`rad-sim.response-atlas-sweep.v1` artifact for comparing locality, residual
vertical motion, operator interaction, and finite-difference parameter
sensitivity, plus falsifiable monotonic operator-law candidates before bench
calibration.
`model_provenance` and `provenance_summary` return the current evidence ledger
for solver features, separating paper-supported equations from assumptions,
diagnostics, and missing calibration data.
`RADHardwareProfile`, `config_with_hardware_profile`, and
`hardware_profile_from_config` provide the first explicit bridge from normalized
simulator controls to measured real-cell dimensions. They track which key
hardware dimensions have been measured and which remain calibration gaps.
`export_hardware_profile_json` / `hardware_profile_from_json` serialize the same
record as `rad-sim.hardware-profile.v1`, so a caliper-measured profile can be
used by command-line bench exports without hard-coding prototype dimensions.
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
superposition residual, fitted decay profile, pairwise interaction graph,
sequence-order error, and framework law candidates are new diagnostics for
detecting when composed operators stop behaving additively or commutatively
because of backlash thresholds, locks, or saturation.
`programmable_discontinuity_report` and
`export_programmable_discontinuity_report_json` serialize that diagnostic as a
`rad-sim.programmable-discontinuity-report.v1` artifact with explicit
paper-supported assumptions, simulator-introduced diagnostics, locality,
reachability, composition, response-matrix, pairwise, optional spring-hinge
physical-validation, sequence-order fields, and a
`rad-sim.formalization-targets.v1` theorem-target manifest for Lean-sized
claims such as backlash dead-zone lemmas, lock idempotence, response-rank
bounds, and noncommutativity witnesses. These are proof targets, not completed
Lean proofs unless Lean/Lake are available and the premises are formalized.
`formalization_target_manifest` and
`export_formalization_target_manifest_json` export that proof-target manifest
without the heavier response fields.
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
`reachable_equilibrium_profile_inverse_report` adds a read-only calibration
bridge: it scores an existing inverse solve against a reachable-equilibrium
empirical profile's alpha/height response scales, leakage tolerance,
group-sequence tolerance, and uncertainty budget. The report is intentionally a
diagnostic certificate only; it does not mutate `LatticeConfig`, rerun the
optimizer, or claim that the empirical profile is a physical law.
`reachable_equilibrium_profile_inverse_acceptance_report` adds the next guard:
it classifies that certificate as accepted for preview, review-required, or
rejected for missing evidence using explicit residual-score, band-failure,
actuator-budget, underactuation, and optional physical-validation thresholds.
Acceptance is still only a preview decision, not hardware execution approval.
`reachable_equilibrium_profile_inverse_preview_packet` packages an accepted
profile-aware inverse plan into a read-only preview/lab handoff artifact with
the exact actuator commands, target summary, residual summary, preview events,
and review protocol. The packet is intended for browser inspection, notebooks,
or supplemental tables; it never applies commands by itself.
`reachable_equilibrium_profile_inverse_preview_replay_report` adds a
deterministic read-only replay certificate for that packet. It rebuilds the
command state, reruns the kinematic simulator, checks command/event coverage,
and compares replayed residuals against packet metadata when target arrays are
supplied. This is still only a reproducibility gate for preview artifacts, not
a hardware execution, spring-hinge, or contact-validity claim.
`reachable_equilibrium_profile_inverse_preview_physical_report` then runs that
same accepted packet through the 3D spring-hinge preview and records solver
success, spring/hinge energy terms, target residuals, and kinematic-versus-
physical model disagreement. It is useful as a physical-preview audit before
lab work, but remains normalized and uncalibrated until contact, friction,
gravity, stiffness, and real geometry are measured.
The browser inverse panel exposes this same chain with `Build Packet`, `Replay
Packet`, `Physical Packet`, and matching save buttons, with readiness summaries
shown directly in the inverse-plan readout.

For the real one-cell CAD path, the A360 share `https://a360.co/4bMlzip` is
tracked as `assets/cad/RADs_unit_cell.f3d` plus a preview and reference summary.
`cad_rad_cell_archive_audit` / `cad_rad_cell_archive_audit.json` inspect that
local Fusion archive directly: ZIP readability, SHA-256, manifest entries,
preview entries, BREP blob entries, linked display labels, and the remaining
segmentation/contact evidence required before exact rigid-body simulation.
If the A360 viewer is exported as STEP/STP, OBJ, or STL, `python -m
rad_sim.cad_mesh_intake` now audits that geometry directly and writes a
`rad-sim.cad-mesh-audit.v1` JSON artifact. STEP/STP support extracts
`CARTESIAN_POINT` coordinates for the CAD envelope, while OBJ/STL support uses
mesh vertices. Bounds provide the side length and thickness fields for a
`rad-sim.hardware-profile.v1`; pin radius, hole radius, backlash, and boss radius
can be passed from calipers or segmented contact surfaces. A bounds-only mesh is
kept as partial measured evidence, while a mesh plus pin/hole/boss values can
reach `mesh-calibrated` geometry. This still does not validate solver accuracy:
exact contact requires segmented pin-hole surfaces and bench-coordinate holdouts.
The two-cell packet also runs `one_cell_cad_export_audit`; it searches for
`assets/cad/RADs_unit_cell.step`, `.stp`, `.obj`, or `.stl`, prefers STEP when
multiple exports are present, and records the derived bounds/profile evidence in
`one_cell_cad_export_audit.json`.
The two-cell bench packet writes fillable segmented-CAD intake templates under
`outputs/two_cell_bench_packet/segmented_cad_intake_templates/`. Copy filled
versions into `assets/cad/segmented/`, add the separated body meshes and at least
one engine handoff file, then rerun the packet export. The generated
`two_cell_segmented_cad_intake_validation.json` is the hard gate for attempting
exact rigid-body contact; passing that gate still does not validate physical
accuracy until external engine results are compared against bench coordinates.
The packet also writes 252 generated MuJoCo proxy files under
`outputs/two_cell_bench_packet/external_fidelity_matrix_mjcf/`, an index at
`two_cell_external_fidelity_matrix_mjcf_index.json`, and a conservative
`two_cell_external_fidelity_matrix_mjcf_run.json` status artifact. Packet export
does not auto-run the external engine; use
`python -m rad_sim.run_two_cell_external_fidelity_mjcf` to attempt the batch
with MuJoCo installed. If MuJoCo is unavailable, the runner still writes an
explicit not-run audit record instead of implying physics results exist. A
`two_cell_physical_response_atlas.json` file is also written from the same
252-case reduced matrix. It groups response by lock mode, hole radius, backlash,
and actuation case, ranks high-slip/contact-penalty cases, and identifies the
first free, lock, and clearance cases to measure on the bench. Treat it as an
experiment-selection artifact, not a physical-accuracy proof. The packet also
materializes those first atlas-ranked cases under
`response_atlas_priority_mjcf/` with
`two_cell_response_atlas_priority_mjcf_index.json` and
`two_cell_response_atlas_priority_mjcf_run.json`, so MuJoCo or bench tracking can
start from the highest-value subset before running all 252 cases.
A completed run also writes
`two_cell_external_fidelity_matrix_mjcf_measurements.csv`, whose connector
marker rows are aligned to the dense fidelity-matrix comparison path. When run
with `--compare`, it also writes
`two_cell_external_fidelity_benchmark_summary.json`, a compact report grouped by
hole radius, backlash, lock mode, and actuation case. That summary is the
recommended first file to inspect after a full 252-case external-engine run.
`python -m rad_sim.fit_two_cell_external_fidelity_correction` then fits an
affine correction profile from reduced predictions to MuJoCo proxy observations
and validates it on held-out case IDs. This is useful for previewing a more
accurate simulator response, but it remains an external-proxy surrogate and not
a real physical calibration until bench holdouts and segmented CAD contact pass.
`python -m rad_sim.apply_two_cell_external_fidelity_correction` applies that
profile back to the dense 756-row table and writes corrected preview predictions
plus baseline/corrected residuals. `python -m
rad_sim.export_two_cell_external_fidelity_web_data` compacts those artifacts into
`web/data/two_cell_external_fidelity_summary.json` and `.js`; `web/index.html`
loads the JS file locally and the Two-Cell Bench panel shows the MuJoCo matrix
coverage, contact/lock agreement, baseline-to-corrected RMS, and the remaining
segmented-CAD/bench-holdout blockers. The same panel's `External proxy case`
selector groups representative upper/middle/lower connector rows by caseId, and
`Apply External Case` transfers that case's backlash, pin/hole radius, lock
flags, and actuation commands into the 1x2 CAD preview while showing corrected
versus observed-proxy right-cell height and vertical slip. This improves the
browser preview, but it still does not claim exact physical fidelity.

## Formal Lean Layer

The `formal/` directory contains a Lean 4 scaffold for the math-first
programmable-discontinuities layer. Its theorem files target initial active-set
and dead-zone claims such as lock/unlock idempotence, admissible-set monotonicity
under constraint changes, disjoint-support operator commutation, finite mode
graphs represented by finite lists, and the three scalar branches of the RAD
backlash law. The first backlash proof is the ordered integer version so the
scaffold typechecks without pulling in the full mathlib real-analysis stack. A
successful `lake build` inside `formal/` is the authority for proof checking.
This layer is deliberately separate from the Python and browser simulator: Lean
currently formalizes the discrete operator skeleton and scalar backlash map, not
CAD-accurate geometry, spring-hinge convergence, contact/friction, vertical
residual mechanics, or inverse-design correctness.

Two framework drafts connect that proof layer to the broader research program:
`docs/programmable_discontinuities_framework.md` is a research-paper-style
draft, and `docs/programmable_discontinuities_monograph.md` is a technical
monograph note with definitions, theorem targets, simulator links, and the next
cell-removal/group-actuation/sheet-contour roadmap.

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

The browser model now treats backlash and vertical lifting as separate stages.
Backlash still uses the dead-zone law to transmit in-plane alpha commands
through the lattice graph, but alpha influence no longer creates height by
itself. `web/constraints.js` first computes the free preferred XY realization,
then applies position locks, fitted walls, or fixed-width channels. Only the
blocked XY residual is mapped into proxy vertical motion through
`z = k * max(0, residual - threshold)^p`. This residual-to-Z law is explicitly
uncalibrated until retroreflector wall/channel data is fitted.

Open `web/one-cell-rotation/index.html` for the isolated one-cell primitive:
two identical CAD-informed cross bodies, a fixed lower body, and a driven upper
body rotating about the shared center-hole axis. This page is intended to
validate the single-cell rotation concept before rebuilding the two-cell
expansion relation.
Open `web/two-cell-attachment/index.html` for the three-cell no-backlash row
primitive. It extends the corrected two-cell attachment by repeating an
A-pattern cell on the far side of Cell B, so the chain is A-B-C. A-B keeps the
same two shared pin relationships, while B-C uses the mirrored opposite-side
pin relationships. The page exposes one uniform rotational actuation slider;
Cell B and Cell C orientations are solved from the no-backlash loop constraints
instead of being independently user-driven. The drive slider is clamped to the
sampled angle range that preserves all four coincident pins and zero detected
inter-cell penetration.
The rendered row is free-floating in the center-of-mass frame: Cell A, Cell B,
and Cell C remain on one fixed centerline with equal A-B and B-C spacing, so
contraction/expansion appears as uniform centerline motion rather than a fixed
cell pulling a moving neighbor.
The `Topology graph` display mode hides the CAD cell geometry and renders the
active cell graph directly: present cells are clickable nodes, present
neighbor couplings are edges, removed cells become red deleted nodes, and
incident deleted edges are shown faintly. This is the closest browser view to
the formal graph-deletion model used in the Lean and Python topology
diagnostics.
The `Paper RAD cell` display mode sits between those two views. It shows the
paper-grounded rotating-square structure as two concentric square parts with
four joints, pin bosses, and hole-clearance rings driven by the normalized
pin/hole controls. It is still normalized geometry, but it reflects the RAD
papers' stated unit-cell structure more directly than the general abstraction.
The `CAD RAD cell` display mode is the first one-cell CAD-derived abstraction
from the Autodesk A360 unit-cell share. It renders the central hub/screw,
radial arms, eight outer pads/holes, and pin-hole clearance rings in the 3D
scene. The downloaded Fusion archive is kept at `assets/cad/RADs_unit_cell.f3d`
with preview `assets/cad/RADs_unit_cell_preview.png` and reference metadata
`assets/cad/RADs_unit_cell_reference.json`; the manifest labels the linked cell
as `RADs free cell 4mm tall 3.4mm hole`. The Lattice panel's `Apply A360 CAD
Profile` button loads the same dimensions as a
`rad-sim.cad-rad-cell-reference-profile.v1` hardware-profile wrapper so the
browser scale, CAD cell display, and exported physical packets use the same
one-cell source. The Two-Cell Bench panel
can run a reduced two-body contact proxy,
clearance sweep, and a one-click `Show 1x2 Bench` fixture with optional
position locks. The same panel can run and save the reduced 252-row two-cell
fidelity matrix over hole radius, backlash, lock mode, and alpha/z actuation.
It can also run and save a contact phase map that classifies every matrix row
as free play, axial contact, vertical contact, state-lock vertical-free, or
position-locked so two-cell experiments can be chosen by mechanism phase rather
than only by response magnitude.
`Run Phase Diagram` builds a denser hole-radius-versus-backlash grid for the
selected two-cell actuation and reports where the reduced model changes phase.
`Run Transition Report` converts those reduced phase changes into transition
brackets and a short lab-fillable measurement list for the adjacent two-cell
cases. It is meant to decide which hole/backlash cases to measure next, not to
claim the real transition boundary before bench data exists.
The panel also shows an exact-CAD physical status line and exports
`rad-sim-two-cell-physical-fidelity-status.json`; that status joins the CAD
archive audit, suite/matrix readiness, expected 756 measurement rows, and the
remaining segmentation/contact/bench-evidence blockers.
The browser two-cell panel is back to the reduced non-Fusion proxy. The sampled
Fusion center-path calibration has been removed because it did not reproduce the
actual two-cell mechanism. The next CAD step should use a Fusion API export of
segmented bodies, joint axes, limits, transforms, contact/gap geometry, and
driven pose sweeps before any CAD-derived motion is reintroduced here.
The Python companion `python -m rad_sim.export_two_cell_bench_packet` writes the
same bench as `cad_rad_cell_layout.json`, `two_cell_bench.json`,
`cad_rad_cell_archive_audit.json`,
`cad_rad_cell_reference_profile.json`,
`one_cell_cad_export_audit.json`,
`two_cell_segmented_cad_readiness.json`,
`two_cell_segmented_cad_readiness.csv`,
`two_cell_segmented_cad_intake_templates.json`,
`two_cell_segmented_cad_intake_templates.csv`,
`two_cell_backlash_sweep.csv`,
`two_cell_actuation_sweep.csv`, `two_cell_actuation_sweep.json`,
`two_cell_quasistatic.json`, `two_cell_quasistatic_sweep.csv`,
`two_cell_quasistatic_sweep.json`, `two_cell_connector_contact.json`,
`two_cell_connector_contact.csv`, `two_cell_connector_contact_sweep.json`,
`two_cell_connector_contact_sweep.csv`,
`two_cell_connector_measurement_template.json`,
`two_cell_connector_measurement_template.csv`, `two_cell_cad_proxy.xml`,
`two_cell_cad_proxy.json`, `two_cell_mujoco_proxy_run.json`,
`two_cell_physical_simulation_suite.json`,
`two_cell_physical_simulation_suite.csv`,
`two_cell_physical_fidelity_matrix.json`,
`two_cell_physical_fidelity_matrix.csv`,
`two_cell_contact_phase_map.json`,
`two_cell_contact_phase_map.csv`,
`two_cell_radius_backlash_phase_diagram.json`,
`two_cell_radius_backlash_phase_diagram.csv`,
`two_cell_radius_backlash_transition_report.json`,
`two_cell_radius_backlash_transition_report.csv`,
`two_cell_cad_contact_decomposition.json`,
`two_cell_cad_contact_decomposition.csv`,
`two_cell_exact_contact_handoff_plan.json`,
`two_cell_exact_contact_handoff_plan.csv`,
`two_cell_fidelity_matrix_measurement_template.json`,
`two_cell_fidelity_matrix_measurement_template.csv`,
`two_cell_external_fidelity_matrix_manifest.json`,
`two_cell_external_fidelity_matrix_manifest.csv`,
`two_cell_external_fidelity_matrix_mjcf_index.json`,
`external_fidelity_matrix_mjcf/FM_###.xml`,
`two_cell_response_atlas_priority_mjcf_index.json`,
`response_atlas_priority_mjcf/FM_###.xml`,
`two_cell_external_results_template.csv`,
`two_cell_external_results_template.json`, `two_cell_physical_test_packet.json`,
`two_cell_physics_validation_report.json`,
`two_cell_physics_validation_report.csv`, and
`two_cell_measurement_template.csv`.
Those artifacts are intended for the first real two-cell experiments: free
actuation, state locking, position locking, left-cell residual response, and
pin/hole clearance sweeps. After filling measured positions/angles into the CSV,
`python -m rad_sim.compare_two_cell_measurements --input ...` exports a
measured-vs-model residual report and a reduced-parameter calibration grid
search.
`two_cell_physical_simulation_suite.json` consolidates the named reduced-physics
tests for two connected cells: nominal contraction/lift, tight and loose hole
clearance, opposite actuation polarity, alpha-only and z-only commands, state
lock, position lock, left-free residual motion, zero/high backlash die-off, and
loose-clearance gravity sag. The CSV companion exposes the same cases as a
flat table for measured coordinates, connector slip, and lock-held observations.
Passing this suite means the internal reduced model is self-consistent; it does
not mean the CAD proxy is physically exact.
`two_cell_physical_fidelity_matrix.json` is the denser reduced-physics ladder:
it sweeps hole radius, backlash, lock mode, and actuation case together, then
summarizes clearance die-off, backlash die-off, lock behavior, and actuation
polarity. It is the main generated table for deciding which two-cell physical
tests to run next.
`two_cell_contact_phase_map.json/csv` derives from the same 252-row matrix and
labels each `FM_###` case by contact/lock phase. It also reports phase counts,
active-contact cases, free-play cases, lock failures, and a measurement-priority
subset for bench or external-engine runs. It is a planning artifact for phase
coverage, not a proof of exact real-cell phase boundaries.
`two_cell_radius_backlash_phase_diagram.json/csv` is the denser two-cell sweep
for the specific radius/backlash question. By default it evaluates a 9-by-9 grid
of hole radius and backlash for the `contract_lift` actuation with free lock
mode, records the phase at each grid point, and lists phase transitions. This is
the artifact to compare against physical tests where holes or pins are changed
to create different backlash distances.
`two_cell_radius_backlash_transition_report.json/csv` derives from that phase
diagram and records the cases on both sides of each predicted phase change. The
CSV includes blank observed phase, z, alpha, connector-slip, and lock-held
fields so it can be filled directly from a two-cell bench sweep.
After filling it, run
`python -m rad_sim.compare_two_cell_transition_measurements --input
outputs/two_cell_bench_packet/two_cell_radius_backlash_transition_report.csv
--out outputs/two_cell_radius_backlash_transition_comparison` to produce a
measured-vs-predicted boundary report. The browser has the same path through
`Load Transition Results` and `Save Transition Comparison`; `Apply Transition
Calibration` can then apply the proposed effective hole-radius/backlash and
vertical-coupling update to the reduced proxy. This update is deliberately
marked as a proxy calibration, not as a verified physical CAD parameter.
The same comparison writer now also creates
`two_cell_radius_backlash_transition_rerun.json/csv`. That rerun keeps the
nominal tested hole/backlash grid fixed while shifting the effective
hole-radius/backlash used by the reduced contact model, then reports how the
predicted transition brackets moved. Use this as a before/after calibration
diagnostic, not as physical validation; independent holdout sweeps and segmented
CAD contact still gate real claims.
`two_cell_fidelity_matrix_measurement_template.csv` expands that dense matrix
into every upper/middle/lower connector row. It includes predicted normalized
cell-center coordinates, predicted connector marker coordinates in millimeters,
predicted connector slip/contact mode, and blank observed fields for the
physical experiment. The companion comparison helper aligns by `caseId` and
`connector`, then reports cell, connector, contact-mode, and lock-held residuals.
After filling this dense template, run
`python -m rad_sim.compare_two_cell_fidelity_matrix_measurements --input
outputs/two_cell_bench_packet/two_cell_fidelity_matrix_measurement_template.csv
--out outputs/two_cell_fidelity_matrix_measurement_comparison`. This writes
`two_cell_fidelity_matrix_measurement_comparison.json/csv` plus
`two_cell_fidelity_matrix_parameter_calibration.json`; the calibration file only
proposes reduced-proxy updates for coupling, vertical coupling, and effective
pin-hole clearance, and does not validate exact physical accuracy. In the browser,
the same fillable CSV is available from the two-cell bench panel as
`Save Matrix Template CSV`; after recording bench or external-engine observed
fields, use `Load Matrix Results` to compute observed-vs-predicted residuals
inside the UI and estimate reduced-proxy parameter scales, then `Save Matrix
Comparison`, `Save Matrix Calibration`, and `Save Physical Status` to archive
the residuals, proposed proxy updates, and exact-CAD readiness blockers. Use
`Apply Matrix Calibration` only as a manual reduced-proxy tuning step; it updates
coupling, vertical coupling, and effective hole radius but still does not make an
exact physical-fidelity claim.
`two_cell_external_fidelity_matrix_manifest.json/csv` is the dense external
engine handoff aligned to the same `FM_###` cases. It lists all 252 radius,
backlash, lock, and actuation cases, expected body states, lock flags, per-case
MJCF/result paths, and the 756 compatible connector measurement rows. Use it
when running MuJoCo/Gazebo/Isaac or bench tracking across the full matrix; it is
a run manifest and still reports segmented-CAD, engine-result, and bench-holdout
evidence as missing until those rows are filled and compared.
The packet writer also materializes the referenced proxy MJCF files under
`external_fidelity_matrix_mjcf/` and records them in
`two_cell_external_fidelity_matrix_mjcf_index.json`, so an external solver can
run the cases directly instead of reconstructing controls from the CSV.
In the browser two-cell panel, `Save Engine Manifest` exports the same
case-aligned external-run manifest from the current UI settings.
`two_cell_connector_measurement_template.csv` expands the same named suite into
upper/middle/lower connector rows. Each row contains predicted connector marker
positions, predicted lateral/vertical slip, predicted contact mode, and blank
observed fields for measured connector positions, slip, observed contact mode,
and lock-held notes.
After filling those observed fields, run
`python -m rad_sim.compare_two_cell_connector_measurements --input
outputs/two_cell_bench_packet/two_cell_connector_measurement_template.csv --out
outputs/two_cell_connector_measurement_comparison` to generate per-connector
marker-position, slip, contact-mode, and lock-held residuals. This is the first
calibration bridge between the A360 one-cell CAD reference and measured
two-cell physical behavior.
The physics validation report is a consolidated gate: it checks internal
CAD-layout envelope consistency, clearance die-off, quasistatic convergence,
state-lock versus position-lock behavior, gravity sag, actuation polarity, and
MJCF proxy readiness, while explicitly leaving physical accuracy incomplete
until external solver or bench rows are loaded.
`cad_rad_cell_layout.json` is the shared one-cell geometry abstraction used by
Python and the browser `CAD RAD cell` mode: it uses the A360 bounding box,
3.4 mm nominal hole label, eight radial pad sites, and three horizontal
two-cell connector pairs. It is dimension-consistent, not an exact BREP mesh.
`two_cell_segmented_cad_readiness.json` is the exact-CAD handoff manifest. It
records which A360 assets are present locally and lists the still-missing
segmented bodies, joint axes, pin/hole contact surfaces, material/contact
parameters, lock-crown geometry, actuator curves, and engine handoff files
needed before the simulator can run exact rigid-body contact instead of a
reduced proxy. Add Fusion exports under `assets/cad/segmented/` using the
filenames in `assets/cad/segmented/README.md`; the manifest will detect them
automatically.
`two_cell_segmented_cad_intake_templates.json` and the generated
`segmented_cad_intake_templates/` folder provide starter files for the missing
exact-contact evidence: joint axes, hole surfaces, mass/inertia, contact
parameters, actuator force/displacement curves, lock-crown geometry, and bench
coordinate truth.
`two_cell_cad_contact_decomposition.json/csv` is the explicit CAD-to-engine
contract for the A360 one-cell geometry assembled into two connected cells. It
names the left/right upper, lower, and pin bodies; the 16 pad holes; the three
X-neighbor connector pairs; and the segmented convex/analytic hole-wall
collision primitives a contact engine should use instead of a single concave
hole mesh.
`two_cell_exact_contact_handoff_plan.json/csv` is the next experimental ladder:
it combines transition-boundary cases and high-response atlas cases into a
ranked two-cell exact-contact/bench checklist. The one-cell A360 CAD archive is
treated as the visual and dimensional source for the real cell, but the plan
keeps `physicalAccuracyValidated=false` until segmented bodies, pin-hole
contact surfaces, joint axes, contact parameters, engine results, and bench
holdouts are supplied.
The connector-contact reports evaluate those three CAD-derived pin-hole pairs
under the current two-cell pose and sweeps, recording lateral slip, vertical
slip, clearance excess, and contact penalty as hole radius, lock mode, and
actuation change.
The quasistatic two-cell solver minimizes an internal clearance/contact energy
with hard coordinate constraints for position locks, a separate state-lock
constraint on the cell angle, axial backlash slack, vertical pin-hole clearance,
ground contact, actuation targets, and gravity potential. This fixes the prior
bad abstraction where a state lock could be treated as forcing the cell back to
ground height: in this layer, a state lock can freeze alpha while still allowing
vertical pin-hole motion unless a position lock is also active.
When SciPy is installed, this solve uses `scipy.optimize.minimize`; otherwise
the two-cell packet path falls back to a small deterministic bounded pattern
search and labels that engine in the JSON. Higher-dimensional spring-hinge and
inverse-design mechanics still require SciPy and will report that dependency
when called without the full Python environment.
The MJCF proxy uses the A360 cell bounding box and 3.4 mm hole label to scale a
two-body model with radial arms, hub, screw, eight pads per cell, connector
pin-hole contact pairs, slide/yaw actuation joints, and fixed bodies for state
or position locks. It is ready for an external MuJoCo run, but it is still a
proxy until Fusion exports separated STEP/mesh bodies and bench friction/contact
data are available.
`two_cell_mujoco_proxy_run.json` is the optional external-engine execution
artifact for that proxy. If the Python `mujoco` package is installed, it loads
and steps the generated MJCF model and records final body-center coordinates.
On machines without MuJoCo, it deliberately reports `mujocoPythonPackage` and
`externalMuJoCoRun` as missing evidence rather than claiming validation.
`two_cell_external_results_template.csv` is the fillable handoff for external
MuJoCo/Gazebo/Isaac runs or tracked bench coordinates for the smaller named
suite. Fill the `external*` columns, then run
`python -m rad_sim.compare_two_cell_external_results --input ...` to compare
those independent body coordinates against the quasistatic solver, including
lock-constraint violations and tolerance failures. For the full 252-case ladder,
use `two_cell_external_fidelity_matrix_manifest.json/csv` plus the dense
fidelity-matrix measurement template instead.
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
radii, and stack offset without requiring those values to alter the solver. It
also switches the measurement readout to `Hardware profile`, which labels the
selected cell with measured pin, hole, clearance, plate, and stack dimensions
when available. `Save Profile` exports the active profile as
`rad-sim.hardware-profile.v1` JSON, and `Load Profile` imports that same schema
or the browser's flat profile shape, applies measured backlash/radii to the
effective normalized grid, and records a `hardware-profile-import` event for
traceability.
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
The cleaned 12-cell lock logs can now be converted into per-cell 3D coordinate
targets with `load_lock_dataset`, `lock_coordinate_dataset`, and
`train_lock_coordinate_mlp`. The browser loads the processed
`web/data/lock_dataset.js` artifact plus `web/lock_surrogate.js`; when the grid
has 12 columns and each row has the same lock mask, the displayed centers are
postprocessed by the measured lock-coordinate surrogate. Regenerate that browser
artifact from the local logs with `python -m rad_sim.export_lock_dataset_web`.
This is an empirical dataset/MLP scaffold, not a validated mechanics law.
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
correction, but it does not automatically mutate the solver. The report now also
includes `rad-sim.calibration-parameter-estimates.v1`, which names provisional
alpha response gain, height response gain, neighbor/direct z-coupling ratio,
measured slip/free-play summary, and force-height contact proxy while marking
backlash and contact law calibration as not identifiable from the current data
unless the required sweep or force/contact measurements exist. The Display panel's
`Calibration residual` overlay colors the remaining per-cell error after that
global gain/bias fit, making local mechanism mismatch easier to distinguish from
uniform scale or offset error. After a results JSON is imported, the calibration
readout also shows the fitted z-coupling estimate, measured slip/free-play
summary, and contact-proxy status directly in the panel; these are inspection
labels, not automatic solver mutations. `Save Model Profile` exports
`rad-sim.calibration-model-profile.v1`, a claim-labeled fitted profile with
safe bounded simulator updates separated from diagnostic-only quantities.
`Apply Model Profile` currently mutates only direct, bounded grid fields such as
`zCouplingGain`; alpha/height response gains, backlash, clearance, and contact
terms remain recorded evidence until a later calibrated solver has explicit
fields for them.
`calibration_bench_notebook` and the browser `Save Bench Notebook` action now
package the response protocol into `rad-sim.calibration-bench-notebook.v1`: a
bench-facing plan with instrument requirements, fit and holdout dataset IDs,
profile-freeze evidence, scenario rows, expected outputs, and explicit
pass/fail criteria. `Save Bench CSV` / `export_calibration_bench_notebook_csv`
writes one row per protocol step for a lab notebook or spreadsheet. This is a
protocol-coverage artifact; it helps enforce train/holdout separation and
measurement columns, but it does not prove that the physical run was collected
independently or that contact/friction/stiffness laws are calibrated.
`calibration_bench_packet` and the browser `Save Bench Packet` action bundle
that notebook with the machine-readable protocol, fit template, holdout
template, artifact manifest, and validation instructions as
`rad-sim.calibration-bench-packet.v1`. For a folder that can be sent to a TRL
collaborator, run
`python -m rad_sim.export_calibration_bench_packet --out outputs/calibration_bench_packet`.
The writer creates separate JSON/CSV/README files so the fit and holdout
templates can be filled as immutable raw bench artifacts.
After those templates are filled with measurements, `python -m
rad_sim.compare_calibration_bench_packet --fit <fit.json> --holdout
<holdout.json> --out outputs/calibration_bench_execution_validation` writes
`rad-sim.calibration-bench-execution-validation.v1`, a one-row execution CSV,
the detailed holdout replay, the model profile, and the fit comparison report.
The browser mirrors this with `Check Holdout`, `Save Execution Gate`, and `Save
Gate CSV`. This gate proves only evidence bookkeeping: fit/holdout samples are
present, a bounded simulator update was replayed, residual and provenance flags
passed, and no required metadata is missing. It still does not prove contact,
friction, stiffness, gravity, or material-law accuracy.
`physical_validation_readiness_report` and the browser
`physicalValidationReadinessReport` compose that execution gate with a filled
vertical-load comparison report into `rad-sim.physical-validation-readiness.v1`.
The report is a conservative physical-claim review gate: it requires passing
fit/holdout calibration evidence, passing vertical-load scenarios, load-work
and contact-proxy fields, configured pin-hole clearance, and zero missing
evidence before a result is treated as ready for human physical-validation
review. It still labels physical accuracy as experimentally unvalidated.
`contact_state_abstraction_report` and browser
`contactStateAbstractionReport` add the first explicit rigid-contact
abstraction layer as `rad-sim.contact-state-abstraction.v1`: every active cell
gets pin, hole, clearance, unilateral contact mode, penetration, and normalized
contact-penalty records, while removed cells are marked as deleted topology.
This is an inspectable contact-state inventory, not collision detection or
calibrated friction/contact stiffness.
`contact_graph_consistency_report` and browser
`contactGraphConsistencyReport` then check the graph/contact boundary as
`rad-sim.contact-graph-consistency.v1`: edges incident to removed cells must be
deleted, removed cells must not carry active contact, active bodies must have
contact records, and group-support cells must appear in the contact inventory.
The report also records lost support when an actuator/support references a
removed cell. This is a finite consistency gate for removable-cell mechanics,
not a proof of physical contact or rigid-body collision behavior.
`physical_realization_map_report` and browser
`physicalRealizationMapReport` add a first explicit realization layer as
`rad-sim.physical-realization-map.v1`: abstract events such as actuation,
locking, cell removal, and group commands are mapped to simulator state effects,
support cells, mechanism-channel labels, contact-graph evidence, and claim
labels. This is the bridge from the operator calculus to the RAD simulator, but
it is still an abstract-to-simulator map until physical actuator, contact, and
motion measurements identify the real hardware channel.
`external_physics_engine_audit_report` and browser
`externalPhysicsEngineAuditReport` add an external-tool readiness layer as
`rad-sim.external-physics-engine-audit.v1`: candidate independent engines such
as MuJoCo, PyBullet, and Project Chrono are checked for rigid-body contact,
gravity, joint-constraint, and headless-run support, then linked to the current
contact abstraction and realization scenario. This audit records whether the
project is ready to run an independent solver; it does not itself run that
solver or validate physical accuracy.
`mujoco_model_export_report` and browser `mujocoModelExportReport` add the
first concrete external-engine artifact as `rad-sim.mujoco-model-export.v1`.
The export emits coarse MJCF XML where each active cell is a box body, fixed
cells omit free joints, removed cells are skipped, gravity is configured, and
external force records can be attached. This is intentionally a normalized
proxy model, not fabrication-accurate CAD.
`mujoco_pin_hole_contact_geometry_report` and browser
`mujocoPinHoleContactGeometryReport` add
`rad-sim.mujoco-pin-hole-contact-geometry.v1`: each active cell receives four
pin records, four hole records, configured clearance records, contact-pair
records, and MJCF proxy geometry fragments. These records move the external
handoff closer to real RAD pin-hole mechanics, but the geometry is still
unmeasured and friction/compliance/contact stiffness remain uncalibrated.
`mujoco_contact_parameter_report` and browser
`mujocoContactParameterReport` add
`rad-sim.mujoco-contact-parameter-profile.v1`: each exported contact pair gets
proxy stiffness, damping, friction, solver-parameter, and XML-attribute records.
This makes MuJoCo contact assumptions explicit and exportable, but the values
remain experimental assumptions unless they are backed by bench measurements.
`contact_parameter_calibration_packet` and browser
`contactParameterCalibrationPacket` add
`rad-sim.contact-parameter-calibration-packet.v1`: a bench handoff that bundles
the pin-hole contact geometry, the contact-parameter profile, fit and holdout
measurement rows, and required columns for slip, force, friction, rebound, and
fitted solver parameters. This is a calibration plan, not calibration evidence.
`contact_parameter_calibration_results_template`,
`compare_contact_parameter_calibration_results`, and browser
`compareContactParameterCalibrationResults` add
`rad-sim.contact-parameter-bench-validation.v1`: blank packet rows can be filled
from bench or external contact-fitting data, then checked for missing fields,
fit residuals, holdout residuals, and independent fit-vs-holdout agreement.
Passing this gate is empirical calibration evidence for the chosen proxy
parameters, not a proof of the contact law.
`contact_parameter_interval_calibration_report` and browser
`contactParameterIntervalCalibrationReport` add
`rad-sim.contact-parameter-interval-calibration.v1`: a bounded-parameter layer
that turns filled fit/holdout rows into intervals for pin radius, hole radius,
clearance, friction, contact stiffness/damping, and MuJoCo solver parameters.
It checks whether the frozen simulator assumptions fall inside the empirical
bounds and records uncertainty and holdout-agreement counts.
`mujoco_external_run_report` runs that MJCF model only when the optional
`mujoco` Python package is actually available. If MuJoCo is unavailable, the
report remains a structured missing-evidence result rather than pretending
validation happened. The browser companion `mujocoExternalRunReport` is a
non-running placeholder because the static UI cannot execute MuJoCo.
`mujoco_external_comparison_report` and browser
`mujocoExternalComparisonReport` compare imported or executed MuJoCo body-center
positions against the simulator state as
`rad-sim.mujoco-external-comparison.v1`. A passing comparison is still only an
external-engine agreement diagnostic; bench measurements are required before
physical accuracy claims.
`equilibrium_relation_report` and browser `equilibriumRelationReport` add a
finite variational-mechanics gate as `rad-sim.equilibrium-relation.v1`: a
realized operator state is checked against spring-hinge solver success,
contact-state evidence, nonnegative spring/hinge/lock/contact/load proxy energy
terms, objective-energy balance, and a residual tolerance. This is a rigorous
bookkeeping relation between operator realization and the current numerical
solver output; it is not yet a proof of continuous minimization, friction,
material laws, or calibrated contact mechanics.
`reachable_equilibrium_controllability_report` and browser
`reachableEquilibriumControllabilityReport` add
`rad-sim.reachable-equilibrium-controllability.v1`: a selected actuator basis is
evaluated through finite response columns on the equilibrium state, then
compared against target cells and removed-cell topology components. The report
separates evidence readiness from full target reachability, so underactuated or
topology-blocked targets are recorded as design information unless strict full
target reachability is requested. This is a linearized finite reachability
proxy, not a nonlinear controllability theorem.
`reachable_equilibrium_bench_protocol` and browser
`reachableEquilibriumBenchProtocol` add
`rad-sim.reachable-equilibrium-bench-protocol.v1`: the finite controllability
report is converted into physical bench trials for baseline equilibrium capture,
per-actuator alpha/z response columns, group target probing, and explicit
topology-blocked target controls. The protocol is a measurement design and
evidence gate; it does not certify that measurements were collected or that the
hardware is physically controllable.
`reachable_equilibrium_bench_results_template` and
`compare_reachable_equilibrium_bench_results` close that loop with
`rad-sim.reachable-equilibrium-bench-results.v1` and
`rad-sim.reachable-equilibrium-bench-comparison.v1`. The template creates
fillable target-level measurement rows; the comparison checks qualitative
reachability flags, topology-blocked target leakage, component-label
consistency, and simultaneous-versus-sequenced group response mismatch. This is
still a bench-result validation scaffold, not a calibrated amplitude law or
proof of nonlinear controllability.
`reachable_equilibrium_amplitude_calibration_report` and browser
`reachableEquilibriumAmplitudeCalibrationReport` add
`rad-sim.reachable-equilibrium-amplitude-calibration.v1`: repeated bench rows
are aggregated into alpha/height amplitude estimates, sample standard
deviations, simple uncertainty half-widths, qualitative residual fields,
topology-response bands, and simultaneous-versus-sequenced group residuals.
The claim is empirical and finite-sample: it summarizes measured response
amplitudes, but it is not yet a constitutive model, a calibrated statistical
law, or a nonlinear controllability proof.
`reachable_equilibrium_empirical_profile_from_amplitude` and browser
`reachableEquilibriumEmpiricalProfileFromAmplitude` add
`rad-sim.reachable-equilibrium-empirical-profile.v1`: the measured amplitude
report is converted into bounded empirical profile proposals for alpha response
scale, height response scale, topology-leakage tolerance, group-sequence
tolerance, and uncertainty budget. These are safe-to-review update proposals
with optional holdout hooks; no `LatticeConfig` field is mutated in v1.
For external inspection, `build_paper_rad_lattice_mesh` converts the normalized
paper RAD lattice into extruded plate, pin, and connector mesh components, and
`export_paper_rad_mesh_obj` serializes that mesh to OBJ text. This is meant for
CAD-style review and downstream tooling; it is not yet a fabrication export.
The browser OBJ export now records the active calibration profile in the header
and uses measured profile defaults for plate thickness, pin radius, stack height,
and connector width when available.
The browser Persistence panel also has `Save OBJ`, which exports the current
visible simulation state as a normalized paper-RAD OBJ mesh.
The Display panel now includes `Sheet contour` as a cell visual mode. It hides
cell geometry while preserving the membrane, target surface, contour lines, and
target/topology overlays for surface-first inspection. A dedicated
`Topology blocked` overlay separates targets isolated by removed-cell component
topology from generic finite-response underactuation.
The Display panel can switch the visible solver between `Kinematic` and
`3D spring preview`. Kinematic mode shows the direct backlash propagation model;
spring preview runs a lightweight browser relaxation over the cell-center graph,
then recomputes height, slope, and link-strain diagnostics for a more physical
inspection view without replacing the Python spring-hinge solver. The
`Model disagreement` overlay uses that same spring preview to color cells by
per-cell height disagreement from the kinematic prediction, so suspicious
regions can be inspected directly in the 3D lattice.
Removed cells are now true graph deletions in this browser preview as well:
they are skipped during neighbor relaxation, excluded from slope/error averages,
and any link strain crossing the deleted cell is marked as an absent spring edge.
The preview reports active and skipped spring-edge counts so intact and
removed-topology views can be compared without mistaking a hidden cell for a
physical coupler.
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
`compare_vertical_residual_under_removal` and browser
`RAD.compareVerticalResidualUnderRemoval` compare that residual field before and
after cell removal, including topology-blocked cells and an uncalibrated
clearance-excess contact penalty. The same diagnostic accepts fixed cells and a
signed external z-load to report normalized load-work and height-contact proxy
terms for gravity/load thought experiments.
`compare_vertical_residual_spring_hinge_3d` runs the same intact/removed z-load
case through `solve_spring_hinge_3d`, whose spring and hinge terms now skip
removed cells. Components with no fixed cell are reported as unanchored and are
excluded from quasistatic external forcing.
`vertical_removal_physical_comparison_to_dict` /
`export_vertical_removal_physical_comparison_json` serialize this as
`rad-sim.vertical-removal-physical-comparison.v1`, while
`build_vertical_load_physical_preview_report` /
`export_vertical_load_physical_preview_report_json` produce a three-scenario
`rad-sim.vertical-load-physical-preview-report.v1` artifact for one-cell,
two-cell, and removed-middle z-load studies.
`export_vertical_load_physical_preview_report_csv` converts that artifact into a
bench-planning table with one row per scenario and columns for fixed cells,
removed cells, command/load values, solver convergence, unanchored components,
spring/hinge deletion counts, model-error metrics, contact/load proxy deltas,
and claim labels.
The spring-hinge solver metadata now separates nonnegative stored energy
terms from signed external-work potential, and each comparison embeds a
`rad-sim.mechanics-energy-certificate.v1` certificate. That certificate
summarizes stored spring/hinge/lock energy, contact/load-magnitude proxy
terms, signed load work, Lean theorem references, and the remaining calibration
limit.
`validate_vertical_load_energy_measurements` compares bench-measured load-cell
heights against the spring-hinge predicted signed load work, load-work
magnitude, and height-contact proxy, exporting
`rad-sim.vertical-load-energy-validation.v1` for later gravity/load-cell tests.
`vertical_load_energy_measurement_template` creates a fillable
`rad-sim.vertical-load-energy-measurement-results.v1` file for the one-cell,
two-cell, and removed-middle load cases, and
`compare_vertical_load_energy_measurement_results` turns filled bench rows into a
`rad-sim.vertical-load-energy-comparison-report.v1` summary.
`vertical_load_energy_experiment_protocol` exports the companion
`rad-sim.vertical-load-energy-experiment-protocol.v1` procedure with fixture
setup, instruments, repeat count, measurement order, safety notes, uncertainty
notes, and expected result fields.
`vertical_load_bench_packet` assembles the preview report, mechanics-certificate
summaries, blank measurement template, experiment protocol, and comparison
instructions into one `rad-sim.vertical-load-bench-packet.v1` handoff artifact.
For a folder you can send to a bench-testing collaborator, run
`python -m rad_sim.export_vertical_load_bench_packet --out outputs/vertical_load_bench_packet`.
It writes the packet JSON, preview JSON, preview CSV, blank measurement template,
experiment protocol, `physical_unit_scale_metadata.json`, `hardware_profile.json`,
and a short README with measurement cautions. Pass
`--hardware-profile path/to/hardware_profile.json` to apply measured backlash,
pin radius, and hole radius to the effective normalized simulation before the
packet is generated.
After the bench team fills `measuredHeight` in the measurement template, run
`python -m rad_sim.compare_vertical_load_measurements --input outputs/vertical_load_bench_packet/vertical_load_energy_measurement_template.json --out outputs/vertical_load_measurement_comparison`.
That writes a claim-labeled comparison report plus a one-row-per-scenario CSV
summary, the same physical unit-scale metadata, and the hardware profile while
preserving the warning that simulator agreement is not calibrated proof of
gravity, friction, fixture compliance, or rigid-body contact. The comparison CLI
also accepts `--hardware-profile` so filled measurements can be evaluated against
the same measured geometry assumptions as the packet export.
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
diagnostic, not a paper-derived material law. The panel also shows how many
framework law candidates are currently supported by the diagnostic predicates
and names the primary supported locality, reachability, composition, or
event-order candidate. A companion formalization readout counts Lean-sized proof
targets and shows whether the primary target is ready, blocked by missing Lean
tooling, or waiting for calibrated premises. This is an early numerical probe for
programmable-discontinuity behavior: nonzero residual indicates that backlash,
locks, or saturation are making the local operators interact non-additively.
The signed z readout separates upward and downward response cells and extrema so
vertical push/pull commands can be checked independently.

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
`Select Under Target` for the worst unreachable target cell. The reachability
diagnostic keeps separate upward and downward target counts so vertical
push/pull requests are not hidden inside a single absolute-height score. The
Python inverse report computes those signed masks from explicit positive and
negative finite probes. Removed-cell topology is reported separately: Python
inverse results include topology-reachable masks plus topology-blocked
alpha/height target masks, and the browser Jacobian report highlights height
targets isolated from every allowed actuator by deleted cells. `Solve Linear
Fit` uses those columns in a first greedy
linearized residual projection and `Apply Linear Fit` commits the resulting
clamped commands. These response columns are still computed from the current
synthetic kinematic model rather than from a calibrated quasistatic mechanism.
`Validate Physical` checks the current linear fit or analyzed plan
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
  SciPy optimization and penalty locks; removed cells delete incident spring and
  hinge terms.
- `solve_spring_hinge_3d`: out-of-plane spring-hinge relaxation over 3D center
  nodes, so vertical actuation and pin-clearance residual height can participate
  in the same axial/hinge/penalty energy model. Removed cells delete incident
  spring and hinge terms in this physical-preview graph as well.
- `local_actuation_event` / `group_actuation_event` / `lock_event` /
  `release_event` / `remove_cell_event` / `restore_cell_event` /
  `compare_event_order` / `compare_sequence_order` /
  `compare_group_actuation_decomposition` /
  `compare_group_actuation_under_removal` /
  `compare_vertical_residual_under_removal` /
  `compare_vertical_residual_spring_hinge_3d`: executable programmable-
  discontinuity operators for sequence studies, lock invariance, topology
  deletion/restoration, group-support decomposition, removed-support loss,
  clearance-gated vertical residual loss, spring-hinge physical-preview
  comparison, locality checks, adjacent-swap sensitivity, and event-order
  noncommutativity.
- `vertical_removal_physical_comparison_to_dict` /
  `export_vertical_removal_physical_comparison_json`: JSON-ready compact or
  field-rich export for a single intact-versus-removed z-load physical preview.
- `build_vertical_load_physical_preview_report` /
  `export_vertical_load_physical_preview_report_json`: named report generator for
  one-cell, two-cell, and removed-middle spring-hinge z-load scenarios.
- `export_vertical_load_physical_preview_report_csv`: table export for the same
  physical-preview report, intended for bench-test planning and paper tables.
- `mechanics_energy_certificate_to_dict`: claim-labeled energy certificate that
  separates nonnegative stored/proxy terms from signed external-work terms.
- `validate_vertical_load_energy_measurements` /
  `export_vertical_load_energy_validation_json`: bench-facing comparison from
  measured load-cell heights to simulated signed work and contact proxies.
- `vertical_load_energy_measurement_template` /
  `compare_vertical_load_energy_measurement_results`: generate fillable
  vertical-load bench rows and summarize completed measurement files.
- `vertical_load_energy_experiment_protocol`: step-by-step physical protocol for
  collecting the vertical-load energy measurement rows.
- `vertical_load_bench_packet`: one-file vertical-load handoff packet combining
  preview, protocol, template, certificate summaries, and comparison guidance.
- `write_vertical_load_bench_packet_artifacts`: Python helper behind
  `python -m rad_sim.export_vertical_load_bench_packet`; writes the packet,
  preview report, CSV bench table, measurement template, experiment protocol,
  physical unit-scale metadata, hardware profile, and README into a chosen
  output folder.
- `write_vertical_load_energy_comparison_artifacts`: Python helper behind
  `python -m rad_sim.compare_vertical_load_measurements`; reads filled
  vertical-load measurement JSON and writes the comparison report, scenario
  summary CSV, physical unit-scale metadata, hardware profile, and README into a
  chosen output folder.
- `physical_unit_scale_metadata` / `export_physical_unit_scale_metadata_json`:
  claim-labeled `rad-sim.physical-unit-scale-metadata.v1` export linking
  normalized-to-millimeter scale metadata to the Lean
  `measurement_unit_scale_invariants` target. These functions accept an optional
  measured `RADHardwareProfile` and report both input and effective normalized
  geometry.
- `lattice_topology_diagnostic`: graph-level component/deleted-edge diagnostic
  for removed-cell topology experiments.
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
  model-disagreement evidence. Reports include diagnostic law candidates for
  locality, rank-limited reachability, non-additive composition, and event-order
  noncommutativity, plus formalization targets that separate paper-supported
  dead-zone lemmas from simulator-diagnostic witnesses.
- `formalization_target_manifest` /
  `export_formalization_target_manifest_json`: standalone Python export for the
  `rad-sim.formalization-targets.v1` proof-target manifest embedded in the
  framework report.
- `characterize_pairwise_interactions`: Python pairwise command interaction graph
  for identifying which actuation operators create non-additive residuals and
  mapping each command cell's strongest pair residual and non-additive degree.
- `response_decay_profile`: Python log-linear shell-max locality diagnostic that
  mirrors the browser Response Experiments decay ratio/length readout.
- `compare_removed_topology_reachability`: matched intact-versus-removed
  response-matrix experiment for topology deletion, reporting alpha/height
  reach loss, rank loss, component-count change, and deleted-edge change.
- `build_topology_experiment_report` / `export_topology_experiment_report_json`:
  named topology scenario report for intact, locked, removed, and group-actuated
  experiments, including component blocking, signed height reach, component-wise
  response rank, scenario deltas, die-off, and target residual diagnostics.
- `solve_inverse_design`: bounded damped least-squares inverse design using
  finite response columns for alpha/height target fields, with residual fields,
  reachable/underactuated masks, topology-blocked target masks for removed-cell
  components, and command-saturation diagnostics.
- `inverse_design_report` / `export_inverse_design_report_json`: Python inverse
  design report artifact containing response-matrix diagnostics, residual fields,
  reachability masks, commands, and optional physical validation.
- `reachable_equilibrium_profile_inverse_report` /
  `export_reachable_equilibrium_profile_inverse_json` /
  `export_reachable_equilibrium_profile_inverse_csv`: calibration-aware inverse
  residual certificate using a reachable-equilibrium empirical profile without
  mutating simulator physics or re-optimizing commands.
- `reachable_equilibrium_profile_inverse_acceptance_report` /
  `export_reachable_equilibrium_profile_inverse_acceptance_json` /
  `export_reachable_equilibrium_profile_inverse_acceptance_csv`: thresholded
  acceptance/review/rejection gate for profile-aware inverse certificates.
- `reachable_equilibrium_profile_inverse_preview_packet` /
  `export_reachable_equilibrium_profile_inverse_preview_packet_json` /
  `export_reachable_equilibrium_profile_inverse_preview_packet_csv`: read-only
  packet of accepted inverse commands, preview events, residual metadata, and
  human-review protocol.
- `validate_inverse_design_physical`: spring-hinge validation of a solved inverse
  command set, reporting physical target error and kinematic/physical model
  disagreement before a command set is treated as mechanically credible.
- `web/operators.js`: browser-side version of the event operators with committed
  `lockAlpha` state, selected-cell order diagnostics,
  `compareGroupActuationDecomposition`, `compareGroupActuationUnderRemoval`,
  `compareVerticalResidualUnderRemoval`, and `compareSequenceOrder`
  adjacent-swap/reversal sensitivity checks in the viewport dock.
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
- `web/math.js` `exportHardwareProfileJson` / `importHardwareProfileJson`:
  browser roundtrip for `rad-sim.hardware-profile.v1`, used by the Lattice
  panel's `Save Profile` and `Load Profile` actions.
- `export_hardware_profile_json` / `hardware_profile_from_json`: Python
  roundtrip for the same `rad-sim.hardware-profile.v1` schema used by the
  bench-packet and measurement-comparison CLIs.
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
  residual vertical motion, non-additive operator interaction, and endpoint
  finite-difference sensitivity studies, then records simulator-diagnostic
  monotonic operator-law candidates.
- `web/analysis.js` `responseAtlas` / `exportResponseAtlas`: browser-side
  `rad-sim.response-atlas.v1` export behind `Save Atlas`, using the current
  protocol, kinematic response fields, and browser spring-preview validation.
- `web/analysis.js` `responseAtlasSweep` / `exportResponseAtlasSweep`:
  browser-side `rad-sim.response-atlas-sweep.v1` export behind `Save Sweep`.
  `Run Sweep` computes the same artifact into the Response Experiments readout
  so backlash/clearance trends, dominant sensitivity, and the clearance law
  candidate can be inspected without leaving the 3D UI.
- `calibration_experiment_results_template` /
  `compare_calibration_experiment_measurements` and the browser matching
  helpers: fillable bench-results schema plus simulator-vs-measurement error
  report for the calibration protocol. The browser can import a completed
  results JSON and summarize the comparison in the Response Experiments panel.
  `calibration_experiment_comparison_report` and the browser `Save Comparison`
  action add per-cell raw/residual fields plus measured-vs-simulated gain and
  bias diagnostics without mutating solver parameters. The report includes
  ranked raw-error and post-fit residual cell lists plus
  `rad-sim.calibration-parameter-estimates.v1` parameter bookkeeping so repeated
  local mismatch and underdetermined physical parameters can be inspected after
  a bench run. `calibration_model_profile_from_report`,
  `apply_calibration_model_profile`, browser `calibrationModelProfile`, and
  browser `Apply Model Profile` convert those estimates into
  `rad-sim.calibration-model-profile.v1` with safe bounded simulator updates
  separated from diagnostic-only fitted quantities.
  `calibration_model_profile_residual_comparison` and browser `Check Model
  Profile` replay the same measurement file before and after a cloned profile
  application, export
  `rad-sim.calibration-model-profile-residual-comparison.v1`, and record an
  audit history. `select_calibration_model_profile` / browser
  `selectCalibrationModelProfile` then selects only candidates with an applied
  update, non-increased missing observations, and non-increased residual score.
  `calibration_model_profile_holdout_validation` and browser `Check Holdout`
  replay the training results and a separately loaded holdout results JSON into
  `rad-sim.calibration-model-profile-holdout-validation.v1`, with separate fit,
  holdout, and independent-validation flags. `Save Holdout CSV` /
  `export_calibration_model_profile_holdout_validation_csv` writes the same
  validation into a two-row fit/holdout bench-summary table with
  `rad-sim.calibration-train-holdout-split.v1` metadata. Result JSON templates
  now include `rad-sim.calibration-dataset-provenance.v1` fields for dataset ID,
  fit/holdout role, source-file ID, collection time, operator, frozen profile
  ID, and profile-freeze timestamp. This is still residual bookkeeping unless
  the holdout file was collected independently from the fit file; it is not a
  proof of contact, friction, or material physics.
- `calibration_bench_notebook` / `export_calibration_bench_notebook_json` /
  `export_calibration_bench_notebook_csv` and browser
  `calibrationBenchNotebook` / `Save Bench Notebook` / `Save Bench CSV`:
  bench-facing calibration notebook exports with instruments, fit/holdout
  dataset-plan metadata, profile-freeze evidence, scenario rows, output
  artifacts, and pass/fail criteria. The corresponding Lean target
  `CalibrationBenchProtocolCoverageNat` proves finite coverage of scenarios,
  instruments, outputs, measurement columns, and two dataset roles; it does not
  prove real physical calibration.
- `calibration_bench_packet` / `export_calibration_bench_packet_json` /
  `write_calibration_bench_packet_artifacts` and browser
  `calibrationBenchPacket` / `Save Bench Packet`: complete calibration bench
  handoff bundle containing the notebook, notebook CSV, protocol, fit template,
  holdout template, artifact manifest, and validation instructions. The CLI
  module `python -m rad_sim.export_calibration_bench_packet` writes those pieces
  as separate lab-facing files. The corresponding Lean target
  `CalibrationBenchPacketCompletenessNat` proves only finite bundle
  completeness.
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
  diagnostics used by the Python framework report, plus diagnostic law
  candidates, formalization-target metadata, and the browser spring-preview
  physical-validation summary. `Save Formal Targets` exports only the
  `rad-sim.formalization-targets.v1` manifest for proof planning.
- `web/inverse.js` `validateInversePlanPhysical`: browser-side physical
  validation of analyzed inverse plans and linear fits against spring-preview
  relaxation, including physical target residual and model-disagreement metrics.
- `web/inverse.js` `buildResponseJacobian`: browser-side target reachability and
  underactuated-height diagnostics for finite response columns.
- `web/analysis.js` `topologyExperimentReport` /
  `exportTopologyExperimentReport`: browser-side topology scenario report for
  component blocking, component response rank, signed vertical reach, and target
  residual export.
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
