# Programmable Discontinuities as a Framework for Reconfigurable Mechanics

Status: research-paper draft v0.1  
Concrete platform: Reconfigurable Auxetic Devices (RADs)  
Scope: general programmable mechanics across auxetics, origami, kirigami,
metamaterials, soft robots, deployable mechanisms, and mechanical computation

## Abstract

Programmable mechanical systems are usually described by fixed geometry with
tunable parameters. RADs suggest a broader viewpoint: the active constraint
architecture itself can be deliberately changed. A lock can engage, a contact
can close, a cell can be removed, an actuator group can rewrite a local
deformation field, or a pin-hole clearance can make neighboring cells respond
only after a threshold. This draft frames those events as programmable
discontinuities: controlled, discrete changes to active constraints,
connectivity, admissible configurations, or constitutive/contact regime.

RADs provide the motivating example, but the same language applies to origami
fold locking, kirigami cut activation, multistable metamaterial switching,
contact-rich soft robots, modular robotic sheets, and mechanical memory. The
working thesis is that programmable mechanics can be studied as coupled
continuous-discrete dynamics on a changing graph of constraints, with
operator-level laws tested by simulation, formalized in Lean where possible,
and validated by physical experiments.

## Claim Classes

Every claim in this framework should carry one of four labels:

- Lean-proven theorem: checked by `lake build` in `formal/`.
- Mathematically plausible conjecture: stated precisely but not yet proved.
- Simulator-derived empirical law: observed in Python/Three.js diagnostics.
- Experimentally unvalidated physical assumption: plausible hardware claim
  requiring measurement.

## Core Definitions

Definition 1: Mechanical system.  
A mechanical system is a tuple

```text
S = (G, Q, q, C, sigma, Pi, A, B, P)
```

where `G` is a cell/body graph, `Q` is a configuration space, `q in Q` is the
current configuration, `C` is a constraint family, `sigma` is the active mode,
`Pi` is a mode-dependent energy or residual, `A` is an actuator set, `B`
contains backlash/contact thresholds, and `P` contains physical realization
parameters such as pin radius, hole radius, friction, stiffness, mass, and
calibration data.

Claim class: mathematically plausible conjecture until a richer Lean structure
is formalized.

Definition 2: Active mode.  
An active mode is a predicate over constraints:

```text
sigma : C -> {inactive, active}
```

The current Lean scaffold represents this as `Mode Constraint := Constraint ->
Bool`.

Claim class: Lean-proven theorem infrastructure.

Definition 3: Admissible configuration set.  
Given a system and active mode, the admissible configurations are

```text
Q_sigma = { q in Q | c(q) holds for every active constraint c in sigma }.
```

In Lean, `MechanicalSystem.admissible` formalizes this for abstract predicates.

Claim class: Lean-proven theorem infrastructure.

Definition 4: Programmable discontinuity operator.  
A programmable discontinuity is a controlled operator

```text
D : S^- -> S^+
```

or at the active-mode abstraction level,

```text
D_sigma : Mode -> Mode.
```

Examples include lock, unlock, actuation command, cell removal, contact
activation, graph deletion, topology rewrite, snap-through branch change, and
clearance-threshold crossing.

Claim class: Lean-proven at the active-mode abstraction; broader physical
interpretation remains a conjectural framework.

Definition 5: Support and locality.  
The support of an operator is the set of constraints, cells, contacts, or graph
edges it can directly rewrite. A locality radius or die-off distance measures
how far nonzero response propagates after backlash, clearance, stiffness, or
membrane coupling.

Claim class: Lean-proven for abstract operator support; simulator-derived for
RAD die-off diagnostics.

Definition 6: Cell removal.  
Cell removal is graph deletion plus constraint deletion:

```text
G' = G \ cell
sigma'(c) = inactive if c touches removed cell, otherwise sigma(c)
```

Claim class: Lean-proven for abstract graph/constraint deletion; physical
realization unvalidated until removed-cell experiments are run.

Definition 7: Group actuation.  
Group actuation is simultaneous or sequenced application of several local
actuator operators. If supports are disjoint and no shared contact/clearance
threshold is crossed, a commutative approximation may hold. Otherwise order
effects and non-additivity are expected.

Claim class: Lean-proven for disjoint abstract support commutation;
simulator-derived for RAD response matrices; physical assumption until tested.

Definition 8: Rigid-body contact state.  
A contact-rich realization has bodies, joints, possible contacts, active-contact
flags, and clearance parameters. For RADs, pin-hole clearance contributes to
backlash and vertical residual motion.

Claim class: typed Lean target plus experimentally unvalidated physical
assumption.

Definition 9: Energy landscape and equilibrium map.  
For each mode `sigma`, an energy or residual `Pi_sigma : Q -> R` defines a
landscape over admissible configurations. The equilibrium map is the relation

```text
Eq(sigma, load) = { q in Q_sigma | q is stationary or minimizes Pi_sigma under load }.
```

This relation may be set-valued because contact, backlash, buckling, and
multistability can produce multiple branches.

Claim class: typed Lean target and variational-mechanics roadmap; not yet
proved for the RAD solver.

Definition 10: Physical realization map.  
The physical realization map connects an abstract operator to hardware:

```text
rho(D_abstract, P) = D_physical
```

where `P` contains geometry, stiffness, mass, clearance, friction, actuator
limits, membrane parameters, and calibration data. Two abstractly identical
operators may realize differently when printed tolerances, friction, or gravity
change.

Claim class: mathematical framework definition; experimentally unvalidated
until calibration data identifies `P`.

## Mathematical Fields Integrated

Graph theory:  
Cells, hinges, locks, removed modules, support, adjacency, and reachability are
graphs or graph rewrites. Cell removal is vertex deletion plus incident-edge
deletion.

Rigidity theory:  
Locks and removed cells change constraints, ranks, mechanisms, and degrees of
freedom. Candidate theorem: a local lock can increase rigidity rank by at most
the number of independent constraints it activates.

Hybrid systems:  
Within each mode the system evolves continuously. Discontinuity operators jump
between modes. Contact activation, lock events, snap-through, and cell removal
are mode transitions.

Variational mechanics:  
Each mode can carry an energy `Pi_sigma(q)`. Equilibria are minimizers or
stationary points subject to active constraints. Current Lean code only proves
toy natural-number nonnegativity facts and records real-valued theorem targets.

Contact mechanics:  
Pins, holes, unilateral contact, clearance, friction, and collision state define
when constraints become active. The current RAD vertical residual model should
be treated as a clearance/contact hypothesis, not a paper-proven law.

Topology and geometry:  
Each active mode defines an admissible set. Mode transitions move between
strata of a stratified configuration space. Origami, kirigami, and auxetic
systems differ in physical realization but share this changing-strata language.

Algebra:  
Events compose. Identity events, associativity, idempotent locks/unlocks, and
commuting disjoint-support operators are the first algebraic laws. Noncommuting
event pairs become candidates for mechanical memory and computation.

Optimization and control:  
Inverse design becomes a search over actuator commands, group operators, lock
sets, removed cells, and target surfaces. Reachability and underactuation must
be reported before claiming shape synthesis.

Numerical analysis:  
Solver outputs require residual norms, constraint violation, conditioning,
continuation robustness, and comparison with external physics engines or
physical experiments.

## RAD Simulator Mapping

| Framework object | Current simulator object | Claim class |
| --- | --- | --- |
| Backlash dead-zone | `rad_sim.coupling.backlash_activation`, `RAD.backlashActivation` | paper-supported formula; Lean-proven integer branches |
| Active lock mode | `LatticeState.locked_mask`, browser `cells.locked` | Lean-proven abstraction; simulator implementation |
| Actuation event | `local_actuation_event`, browser `RAD.localActuationEvent` | simulator operator |
| Event sequence | `apply_event_sequence`, browser `RAD.applyEventSequence` | Lean-proven sequence algebra; simulator diagnostic |
| Disjoint-support commutation | `compare_event_order`, `compare_sequence_order` | Lean-proven abstract theorem; simulator-derived violations when supports interact |
| Die-off radius | `finite_die_off_radius`, `response_decay_profile` | simulator-derived empirical law |
| Cell removal | `LatticeState.removed_mask`, browser `cells.removed`, remove/restore events | Lean-proven graph deletion; simulator topology implementation |
| Group actuation | `group_actuation_event`, browser `RAD.groupActuationEvent`, brush/paint event capture | simulator operator; Lean-proven group-support list facts |
| Group decomposition diagnostic | `compare_group_actuation_decomposition`, browser `RAD.compareGroupActuationDecomposition` | simulator-derived equality check for current event semantics |
| Group under removal diagnostic | `compare_group_actuation_under_removal`, browser `RAD.compareGroupActuationUnderRemoval` | simulator-derived support/topology-loss diagnostic |
| Graph topology diagnostic | `lattice_topology_diagnostic`, browser `RAD.topologyDiagnostics` | simulator graph abstraction |
| Graph topology view | browser `Topology graph` cell visual mode | simulator UI diagnostic for graph deletion/reachability |
| Topology-blocked inverse targets | Python `topology_blocked_*_mask`, browser `topologyBlockedHeightMap` | simulator-derived graph reachability diagnostic |
| Vertical residual | `vertical_clearance_operator`, `computeVerticalResidual` | experimentally unvalidated contact/clearance hypothesis |
| Vertical residual under removal | `compare_vertical_residual_under_removal`, browser `RAD.compareVerticalResidualUnderRemoval` | simulator-derived clearance/topology/contact diagnostic |
| Fixed-load proxy | `fixed_cells` and `external_z_load` options in vertical-removal diagnostics | Lean-proven finite magnitude scaffold; experimentally unvalidated load model |
| Signed load-work residual | `python -m rad_sim.compare_vertical_load_measurements` and signed work fields in comparison reports | Lean-proven finite `Int` sign-convention scaffold; external-load calibration remains experimental |
| Integer-scaled mechanics scaffold | `ScaledNatQuantity`, `ScaledIntQuantity`, and scaled numerator theorems in `formal/` | Lean-proven denominator-carrying bridge toward rational mechanics; not yet a `Rat`/`Real` field proof |
| Measurement unit scaling | `MeasurementUnitScale`, `physical_unit_scale_metadata`, `rad-sim.physical-unit-scale-metadata.v1`, `rad-sim.hardware-profile.v1`, and `calibrate_paper_rad_config` | Lean-proven finite zero/denominator invariants plus hardware-profile coverage metadata; calibration remains experimental |
| Spring-hinge load comparison | `compare_vertical_residual_spring_hinge_3d`, `solve_spring_hinge_3d` | simulator-derived physical preview; uncalibrated |
| Spring-hinge load report | `build_vertical_load_physical_preview_report`, `export_vertical_load_physical_preview_report_json` | simulator-derived research artifact with explicit claim labels |
| Spring-hinge load CSV | `export_vertical_load_physical_preview_report_csv` | simulator-derived bench-planning table with scenario commands, solver status, topology deltas, model errors, proxy load/contact terms, and claim labels |
| Mechanics energy certificate | `mechanics_energy_certificate_to_dict`, `rad-sim.mechanics-energy-certificate.v1` | Lean-proven finite nonnegative proxy scaffold; signed work/contact calibration remains unvalidated |
| Calibration parameter estimates | `calibration_experiment_comparison_report`, browser `calibrationParameterEstimates`, and `rad-sim.calibration-parameter-estimates.v1` | Lean-proven finite residual bookkeeping; fitted parameters remain empirical |
| Calibration model profile | `calibration_model_profile_from_report`, browser `calibrationModelProfile`, and `rad-sim.calibration-model-profile.v1` | claim-labeled safe-update gate for bounded simulator parameters; physical-law status remains unproved |
| Calibration profile residual selection | `calibration_model_profile_residual_comparison`, `select_calibration_model_profile`, browser `Check Model Profile`, and `rad-sim.calibration-model-profile-selection.v1` | Lean-proven finite candidate-selection predicate; residual improvement is not independent physical validation |
| Calibration profile holdout validation | `calibration_model_profile_holdout_validation`, browser `Check Holdout`, and `rad-sim.calibration-model-profile-holdout-validation.v1` | Lean-proven finite train/holdout gate; independence depends on bench-data collection protocol |
| Calibration holdout CSV and split metadata | `export_calibration_model_profile_holdout_validation_csv`, browser `Save Holdout CSV`, and `rad-sim.calibration-train-holdout-split.v1` | bench-summary table plus Lean-proven finite split metadata predicate; file provenance remains experimental |
| Calibration dataset provenance | `rad-sim.calibration-dataset-provenance.v1`, result-template provenance fields, browser fit/holdout import roles | Lean-proven finite provenance predicate for known IDs, distinct files, roles, frozen profile, and matching holdout profile ID |
| Calibration bench notebook | `calibration_bench_notebook`, browser `Save Bench Notebook` / `Save Bench CSV`, and `rad-sim.calibration-bench-notebook.v1` | Lean-proven finite protocol-coverage predicate; physical independence and calibration remain lab obligations |
| Calibration bench packet | `calibration_bench_packet`, `python -m rad_sim.export_calibration_bench_packet`, browser `Save Bench Packet`, and `rad-sim.calibration-bench-packet.v1` | Lean-proven finite packet-completeness predicate; raw-file independence and physical calibration remain external evidence |
| Executed calibration bench validation | `calibration_bench_execution_validation`, `python -m rad_sim.compare_calibration_bench_packet`, browser `Save Execution Gate` / `Save Gate CSV`, and `rad-sim.calibration-bench-execution-validation.v1` | Lean-proven finite evidence gate for positive fit/holdout samples, applied update, residual plus independent-validation flags, and zero missing evidence; physical-law calibration remains external |
| Measured load-work validation | `validate_vertical_load_energy_measurements`, `rad-sim.vertical-load-energy-validation.v1` | Lean-proven zero-residual scaffold plus simulator-vs-bench comparison hook |
| Filled measurement comparison | `python -m rad_sim.compare_vertical_load_measurements`, `vertical_load_energy_comparison_summary.csv` | Lean-proven finite pass/fail predicate for zero-error complete scenarios; physical calibration remains experimental |
| Physical validation readiness | `physical_validation_readiness_report`, browser `physicalValidationReadinessReport`, and `rad-sim.physical-validation-readiness.v1` | Lean-proven finite evidence-composition gate linking executed calibration, vertical-load pass, load/contact proxy terms, clearance configuration, and zero missing evidence; not a contact/gravity proof |
| Contact-state abstraction | `contact_state_abstraction_report`, browser `contactStateAbstractionReport`, and `rad-sim.contact-state-abstraction.v1` | Lean-proven finite inventory gate for active bodies, pins, holes, clearance pairs, contact-state records, penalty terms, and zero missing evidence; physical contact remains uncalibrated |
| Contact-graph consistency | `contact_graph_consistency_report`, browser `contactGraphConsistencyReport`, and `rad-sim.contact-graph-consistency.v1` | Lean-proven finite graph/contact/support gate for deleted removed-cell edges, inactive removed contacts, active-body contact records, represented group supports, and zero missing evidence; physical contact remains uncalibrated |
| Physical realization map | `physical_realization_map_report`, browser `physicalRealizationMapReport`, and `rad-sim.physical-realization-map.v1` | Lean-proven finite map from abstract operators to simulator state effects, support records, claim labels, and contact-graph evidence; hardware realization remains experimental |
| External physics engine audit | `external_physics_engine_audit_report`, browser `externalPhysicsEngineAuditReport`, and `rad-sim.external-physics-engine-audit.v1` | Lean-proven finite readiness gate for independent rigid-body/contact engine candidates, available tool records, required features, validation scenarios, contact-model records, and zero missing evidence; external solver execution remains future evidence |
| MuJoCo model export | `mujoco_model_export_report`, browser `mujocoModelExportReport`, `export_mujoco_model_xml`, and `rad-sim.mujoco-model-export.v1` | Lean-proven finite MJCF export gate for model records, active bodies, fixed bodies, gravity, XML bytes, and zero missing evidence; geometry is a normalized proxy |
| MuJoCo pin-hole contact geometry | `mujoco_pin_hole_contact_geometry_report`, browser `mujocoPinHoleContactGeometryReport`, `export_mujoco_pin_hole_contact_geometry_csv`, and `rad-sim.mujoco-pin-hole-contact-geometry.v1` | Lean-proven finite geometry-inventory gate for pin, hole, clearance, contact-pair, active-pair, MJCF-fragment, and zero missing evidence records; proxy contact geometry remains uncalibrated |
| MuJoCo contact parameter profile | `mujoco_contact_parameter_report`, browser `mujocoContactParameterReport`, `export_mujoco_contact_parameter_csv`, and `rad-sim.mujoco-contact-parameter-profile.v1` | Lean-proven finite profile gate for contact-pair, stiffness, damping, friction, solver-parameter, XML-attribute, and zero missing evidence records; parameter values remain uncalibrated until measured |
| Contact-parameter calibration packet | `contact_parameter_calibration_packet`, browser `contactParameterCalibrationPacket`, `export_contact_parameter_calibration_packet_csv`, and `rad-sim.contact-parameter-calibration-packet.v1` | Lean-proven finite packet gate linking contact geometry/profile evidence to fit and holdout measurement rows for slip, force, friction, rebound, and fitted solver parameters; blank rows are not calibration evidence |
| Contact-parameter bench validation | `compare_contact_parameter_calibration_results`, browser `compareContactParameterCalibrationResults`, `export_contact_parameter_bench_validation_csv`, and `rad-sim.contact-parameter-bench-validation.v1` | Lean-proven finite filled-result gate for ready packet evidence, fit rows, holdout rows, completed measurements, parameter residuals, fit pass, holdout pass, independent holdout agreement, and zero missing evidence; calibrated contact mechanics remains empirical |
| Contact-parameter interval calibration | `contact_parameter_interval_calibration_report`, browser `contactParameterIntervalCalibrationReport`, `export_contact_parameter_interval_calibration_csv`, and `rad-sim.contact-parameter-interval-calibration.v1` | Lean-proven finite interval gate requiring a passing bench validation, parameter intervals, accepted interval records, uncertainty records, holdout agreement records, simulator parameters inside bounds, and zero missing evidence; intervals are empirical bounds, not constitutive laws |
| MuJoCo external run | `mujoco_external_run_report`, browser `mujocoExternalRunReport`, and `rad-sim.mujoco-external-run.v1` | Lean-proven finite run-completeness gate for engine availability, ready export, body-result coverage, completed steps, and zero missing evidence; browser cannot execute MuJoCo |
| MuJoCo external comparison | `mujoco_external_comparison_report`, browser `mujocoExternalComparisonReport`, `export_mujoco_external_comparison_csv`, and `rad-sim.mujoco-external-comparison.v1` | Lean-proven finite simulator-vs-external comparison gate for matched body records and tolerance pass; bench validation remains external evidence |
| Equilibrium relation | `equilibrium_relation_report`, browser `equilibriumRelationReport`, and `rad-sim.equilibrium-relation.v1` | Lean-proven finite variational evidence gate tying realization readiness, solver success, contact-state readiness, nonnegative energy terms, residual tolerance, and zero missing evidence; continuous mechanics remains unproved |
| Reachable equilibrium controllability | `reachable_equilibrium_controllability_report`, browser `reachableEquilibriumControllabilityReport`, and `rad-sim.reachable-equilibrium-controllability.v1` | Lean-proven finite response-matrix reachability gate over equilibrium evidence, actuator basis, target cells, topology components, optional full-target reachability policy, and zero missing evidence; nonlinear controllability remains unproved |
| Reachable equilibrium bench protocol | `reachable_equilibrium_bench_protocol`, browser `reachableEquilibriumBenchProtocol`, and `rad-sim.reachable-equilibrium-bench-protocol.v1` | Lean-proven finite protocol gate turning controllability evidence into baseline, response-column, group-probe, and topology-blocked target measurement trials; bench execution and physical controllability remain external evidence |
| Reachable equilibrium bench validation | `reachable_equilibrium_bench_results_template`, `compare_reachable_equilibrium_bench_results`, browser `compareReachableEquilibriumBenchResults`, and `rad-sim.reachable-equilibrium-bench-comparison.v1` | Lean-proven finite filled-result validation gate for completed target measurements, qualitative reachability checks, topology-leakage checks, group simultaneous/sequenced checks, and zero missing evidence; calibrated physical laws remain external evidence |
| Reachable equilibrium amplitude calibration | `reachable_equilibrium_amplitude_calibration_report`, browser `reachableEquilibriumAmplitudeCalibrationReport`, and `rad-sim.reachable-equilibrium-amplitude-calibration.v1` | Lean-proven finite repeated-trial amplitude gate for measured alpha/height estimates, residual fields, uncertainty bands, topology-response bands, and group-sequence residuals; statistical and constitutive mechanics laws remain external evidence |
| Reachable equilibrium empirical profile | `reachable_equilibrium_empirical_profile_from_amplitude`, browser `reachableEquilibriumEmpiricalProfileFromAmplitude`, and `rad-sim.reachable-equilibrium-empirical-profile.v1` | Lean-proven finite bounded-profile gate for alpha/height response scales, leakage and group-sequence tolerances, uncertainty metadata, safe proposals, and holdout hooks; no simulator mutation or physical-law proof in v1 |
| Profile-aware inverse diagnostic | `reachable_equilibrium_profile_inverse_report`, browser `reachableEquilibriumProfileInverseReport`, and `rad-sim.reachable-equilibrium-profile-inverse.v1` | Lean-proven finite diagnostic gate connecting a ready empirical profile to an existing inverse solve through target cells, residual records, normalized score terms, and read-only profile use; it is not a nonlinear reachability proof or calibrated physical inverse solver |
| Profile-aware inverse acceptance | `reachable_equilibrium_profile_inverse_acceptance_report`, browser `reachableEquilibriumProfileInverseAcceptanceReport`, and `rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1` | Lean-proven finite threshold gate for preview acceptance using a ready profile-inverse report, scaled residual-score limit, band-failure limit, actuator budget, and zero missing evidence; it is not hardware execution authorization |
| Profile-aware inverse preview packet | `reachable_equilibrium_profile_inverse_preview_packet`, browser `reachableEquilibriumProfileInversePreviewPacket`, and `rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1` | Lean-proven finite handoff gate for accepted inverse commands, preview event records, target/residual records, read-only status, and zero missing evidence; it supports inspection and lab review without applying commands |
| Profile-aware inverse preview replay | `reachable_equilibrium_profile_inverse_preview_replay_report`, browser `reachableEquilibriumProfileInversePreviewReplayReport`, and `rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1` | Lean-proven finite replay gate for accepted packet readiness, exact command replay, event coverage, simulator/target residual records, residual agreement, read-only replay, and zero missing evidence; it checks reproducibility of a preview packet, not physical execution |
| Profile-aware inverse physical preview | `reachable_equilibrium_profile_inverse_preview_physical_report`, browser `reachableEquilibriumProfileInversePreviewPhysicalReport`, and `rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1` | Lean-proven finite evidence gate for replay readiness, physical-preview solver success, command and target-residual records, model-comparison records, energy/proxy records, read-only status, and zero missing evidence; the physical model remains uncalibrated |
| Vertical-load bench template | `vertical_load_energy_measurement_template`, `compare_vertical_load_energy_measurement_results` | fillable physical experiment bridge for one-cell, two-cell, and removed-middle load tests |
| Vertical-load protocol | `vertical_load_energy_experiment_protocol`, `rad-sim.vertical-load-energy-experiment-protocol.v1` | step-by-step experimental procedure with instruments, repeats, safety, and uncertainty notes |
| Vertical-load bench packet | `vertical_load_bench_packet`, `rad-sim.vertical-load-bench-packet.v1` | publication/bench handoff bundle for preview, protocol, template, certificates, and comparison guidance |
| Browser removed-edge spring preview | `RAD.simulatePhysicalRelaxation`, `physicalActiveSpringEdges`, `physicalSkippedSpringEdges` | simulator diagnostic tying UI spring relaxation to graph deletion semantics |
| Rigid contact | `physics.js` preview and hardware profile fields | roadmap; not calibrated |
| Sheet-only contour view | browser `Sheet contour` cell visual mode plus membrane/contour controls | simulator UI diagnostic |

## Theorem and Conjecture Table

| Statement | Field | Current evidence | Lean status | Simulator experiment | Physical experiment |
| --- | --- | --- | --- | --- | --- |
| Locking twice is equivalent to locking once. | algebra | active-mode logic | proved | none | lock repeatability |
| Unlocking twice is equivalent to unlocking once. | algebra | active-mode logic | proved | none | release repeatability |
| Adding active constraints shrinks admissible configurations. | geometry/constraints | predicate logic | proved | none | compare lock/free motion |
| Disjoint-support operators commute. | algebra/locality | support logic | proved | run event pairs | physically actuate disjoint regions |
| Group supports include every listed cell and empty groups act as neutral mode operators. | algebra/support | list-supported operator logic | proved | group-decomposition diagnostic | group actuation repeatability |
| Disjoint list-supported group operators commute in the abstract mode model. | algebra/locality | list disjointness plus support commutation | proved | compare disjoint group events | physically actuate separated groups |
| Removed-cell constraint deletion clears constraints touching the removed cell even after a group operator is applied. | graph/support | constraint-touch deletion plus group operator abstraction | proved | group-under-removal diagnostic | remove a grouped cell and measure lost support |
| Clearance excess is zero inside the pin-hole gap, vertical residual transmission is zero inside that gap, and the corresponding clearance-excess contact penalty is nonnegative. | contact/algebra | discrete clearance-excess model | proved for `Nat` | vertical-removal residual diagnostic | one-cell and pair z-clearance sweep |
| Load-work magnitude is nonnegative and fixed zero-displacement cells contribute zero load work. | load/variational mechanics | finite load-work proxy | proved for `Nat` | fixed-cell z-load diagnostic | clamp cells and measure gravity sag |
| Signed finite load work is zero for zero load or zero displacement, and signed residual triples vanish when signed measured and simulated quantities agree. | load/validation mechanics | signed integer work scaffold | proved for `Int` | `compare_vertical_load_measurements` report | load-cell sign-convention check |
| Integer-scaled spring/contact energies and signed residuals preserve denominators and retain zero/nonnegative numerator facts. | numerical mechanics | numerator/denominator records | proved for `Nat`/`Int` records | scaled artifact export target | compare normalized units to calibrated units |
| Measurement-unit scaling preserves zero commands, zero equality residuals, and denominator metadata for unsigned and signed finite quantities. | calibration/numerical analysis | finite unit-scale records | proved for `Nat`/`Int` records | paper-scale calibration export | verify normalized-to-mm conversion against bench measurements |
| A mechanics certificate can separate nonnegative stored/proxy terms from signed external-work terms. | variational mechanics | finite stored-energy sum | proved for `Nat` | spring-hinge energy-certificate report | compare load-cell work and measured deflection |
| A calibration fit with at least one finite sample and zero fitted residual passes the finite residual predicate. | calibration/numerical analysis | finite residual bookkeeping | proved for `Nat` record | calibration parameter-estimate report | one-cell and two-cell response replay |
| A calibration model-profile update is safe only when it has a positive sample count and the proposed bounded parameter lies inside declared bounds. | calibration/numerical analysis | finite bounded-update gate | proved for `Nat` record | model-profile apply audit | replay fitted z-coupling against pair measurements |
| A calibration model-profile candidate is selectable only when it has an applied update, non-increased missing observations, and non-increased finite residual score. | calibration/numerical analysis | finite candidate-selection predicate | proved for `Nat` record | model-profile before/after residual comparison | hold out a second bench file for independent validation |
| A calibration model-profile holdout validation passes only when both fit and holdout residual scores and missing-observation counts do not increase after an applied update. | calibration/numerical analysis | finite train/holdout predicate | proved for `Nat` record | load separate fit and holdout JSON files | collect holdout measurements after freezing the fitted profile |
| A finite train/holdout split is ready only when fit samples and holdout samples are positive, overlap is zero, and the profile is frozen. | calibration/provenance | finite split metadata predicate | proved for `Nat` record | holdout CSV export and split schema | enforce file IDs and profile-freeze timestamp in bench notebook |
| A finite train/holdout provenance record is ready only when fit and holdout dataset IDs are known, dataset and source-file IDs are distinct, roles are marked, the profile is frozen, and the holdout names the frozen profile. | calibration/provenance | finite provenance predicate | proved for `Nat` record | provenance warnings in holdout validation | maintain lab notebook IDs and immutable raw files |
| A calibration bench notebook is coverage-ready only when it has scenarios, instruments, output artifacts, measurement columns, and fit plus holdout dataset roles. | calibration/protocol design | finite protocol-coverage predicate | proved for `Nat` record | bench notebook JSON/CSV export | execute the recorded fit and holdout protocols with immutable raw files |
| A calibration bench packet is complete only when it contains notebook artifacts, protocol artifacts, fit and holdout templates, validation instructions, and a sufficiently populated artifact manifest. | calibration/protocol design | finite packet-completeness predicate | proved for `Nat` record | packet JSON plus folder writer | send artifact folder to a collaborator and preserve returned raw files |
| A filled calibration bench execution is validation-ready only when fit and holdout samples are positive, an update was applied, residual and independent-validation flags pass, and missing evidence is zero. | calibration/provenance | finite executed-evidence predicate | proved for `Nat` record | execution-validation JSON/CSV writer and browser export | collect real independent fit/holdout files and audit immutable provenance |
| A measured-versus-simulated finite work residual is zero when the measured and simulated quantities are equal. | validation/numerical analysis | absolute finite residual | proved for `Nat` | vertical-load energy validation | load-cell and motion-capture replay |
| A complete vertical-load comparison row with zero finite signed-work, magnitude, and contact-proxy errors passes the finite predicate, while a row with missing measurements fails. | validation/numerical analysis | finite pass/fail predicate | proved for `Nat` | `compare_vertical_load_measurements` CLI and summary CSV | fill template with load-cell/motion-capture heights |
| A physical-validation readiness packet passes only when executed calibration is ready, vertical-load comparison passes, load/contact proxy terms are present, clearance is configured, and no evidence is missing. | calibration/contact validation | finite evidence-composition predicate | proved for `Nat` record | readiness JSON/CSV export | human audit before contact/gravity/material claims |
| A contact-state abstraction passes only when active bodies have matching pin/hole/clearance inventories, contact-state and penalty records cover them, clearance is configured, and evidence is complete. | contact mechanics/graph inventory | finite contact-state predicate | proved for `Nat` record | per-cell contact-state JSON/CSV export | compare inferred contact modes to motion/contact sensing |
| A contact graph is internally consistent only when removed-cell incident edges are deleted, removed cells have no active contact, active bodies have contact records, group-support cells are represented, and evidence is complete. | graph/contact mechanics | finite contact-graph predicate | proved for `Nat` record | contact-graph JSON/CSV export | compare removed-cell contacts and support loss to motion/contact sensing |
| A physical realization map is ready only when abstract operators have simulator state effects, support records, claim labels, contact-graph evidence, and no missing realization evidence. | hybrid systems/operator theory | finite abstract-to-simulator predicate | proved for `Nat` record | realization-map JSON/CSV export | identify real actuator/contact channels for each operator |
| A MuJoCo contact-parameter profile is ready only when each contact pair has parameter, friction, solver, stiffness, damping, XML-attribute, and zero missing-evidence records. | contact mechanics/external validation | finite contact-parameter predicate | proved for `Nat` record | contact-parameter JSON/CSV export | fit stiffness, damping, friction, and solver settings against pin-hole contact measurements |
| A contact-parameter calibration packet is complete only when profile evidence is attached and each contact-pair inventory has fit and holdout measurement rows, measurement columns, manifest entries, and zero missing evidence. | contact calibration/provenance | finite packet-completeness predicate | proved for `Nat` record | contact calibration packet JSON/CSV export | collect independent fit/holdout pin-hole force, slip, friction, rebound, and solver-parameter measurements |
| A finite equilibrium relation is ready only when realization evidence, solver success, contact evidence, nonnegative energy terms, residual tolerance, and zero missing evidence all hold. | variational mechanics/numerical analysis | finite equilibrium evidence predicate | proved for `Nat` record | equilibrium-relation JSON/CSV export | compare solver residuals to calibrated force/displacement data |
| A reachable-equilibrium controllability report is ready only when equilibrium evidence, an actuator basis, response columns, target cells, reachable responses, topology evidence, optional target policy, and zero missing evidence all hold. | control/reachability/graph theory | finite response-matrix predicate | proved for `Nat` record | reachable-equilibrium JSON/CSV export | compare predicted reachable targets to physical actuation trials |
| Removed-topology vertical load cases can be compared against a 3D spring-hinge physical preview on the same deleted spring graph. | numerical mechanics | center-node spring-hinge solver | simulator diagnostic | spring-hinge load comparison | gravity/load bench test |
| Browser spring-preview relaxation treats removed cells as deleted nodes and marks incident link-strain entries as absent edges. | graph/numerical mechanics | UI spring-preview validator plus graph deletion theorems | simulator diagnostic using proved graph-deletion theorem references | remove middle cell in browser spring preview | compare intact/removed physical sheet response |
| Cell removal deletes incident edges. | graph theory | graph deletion | proved | removed-cell mask blocks propagation | remove/disable printed cell |
| Removed cells prevent one-step reachability through the deleted node while preserving nonincident one-step edges. | graph theory | graph deletion plus one-step reachability | proved | one-by-three removed-middle topology test | removable-cell lattice test |
| Removed cells prevent any positive-length abstract path starting or ending at the removed node. | graph theory | endpoint-present path abstraction | proved | component-blocked scenario report | removable-cell lattice test |
| Removed cells prevent present-cell reflexive reachability from or to the deleted node. | graph theory | endpoint-present reachability closure | proved | component-blocked scenario report | removable-cell lattice test |
| Backlash is zero inside the dead zone and linear outside. | nonsmooth analysis | RAD formula | proved for `Int`; real target pending | sweep commands across gap | single-joint calibration |
| Larger backlash decreases locality radius. | graph/numerics | response atlas | empirical law | sweep `backlash` | print/mechanically vary clearance |
| Vertical residual die-off follows pin-hole clearance. | contact mechanics | current model intuition | unvalidated assumption | sweep clearance and z commands | motion tracking under vertical actuation |
| Removed cells can re-route deformation fields. | graph/rigidity | framework conjecture | target pending | compare intact/removed rank and reach | removable-cell lattice test |
| Event noncommutativity can encode mechanical memory. | algebra/hybrid systems | simulator order diagnostics | conjecture | sequence-order atlas | repeat event-order experiments |

## Simulator and UI Roadmap

Priority A: research primitives

- Use the implemented component-wise response rank and topology report to
  compare locked, removed, and group-actuated scenarios.
- Extend beyond finite-response diagnostics into nonlinear component-wise
  controllability and die-off radius after deletion.
- Compare simultaneous group actuation events against sequenced local events
  over disjoint and overlapping supports.

Priority B: physical abstraction

- Add rigid-body contact representation: bodies, pin joints, hole radii,
  clearance, contact active flags, and collision/penetration penalty preview.
- Keep it separate from the calibrated solver until pin/hole/friction
  measurements exist.

Priority C: visualization

- Use the implemented sheet-only contour mode for target residual, topology
  blocked, removed-cell, and signed vertical reach inspection.
- Use the implemented topology-graph display mode to inspect removed-cell
  deletion, active couplings, and selected supports beside the CAD and sheet
  views.
- Add removed-cell visual state and deletion/restoration controls.

Priority D: inverse design

- Extend response matrices to include explicit group-actuation columns and
  topology-change candidate columns.
- Upgrade the implemented alpha/height/topology reachability report into a
  constrained nonlinear inverse solver with topology-aware actuator placement.

## Next Concrete Experiment

The next experiment should be a two-cell and three-cell topology test:

1. Simulate baseline actuation response.
2. Remove the middle or neighboring cell in the graph abstraction.
3. Recompute alpha/height reachability and die-off radius. The current Python
   entrypoints are `compare_removed_topology_reachability` for a matched
   two-state deletion comparison and `build_topology_experiment_report` for
   intact/locked/removed/group scenario sets.
4. Compare whether cell removal re-routes, isolates, or amplifies response.
5. Repeat physically by disabling or removing one printed cell if the hardware
   can be modified without damage.

This directly tests whether topology rewriting deserves equal status with
locks and actuators as a programmable discontinuity.
