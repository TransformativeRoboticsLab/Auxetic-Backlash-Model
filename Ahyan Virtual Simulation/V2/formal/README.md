# RAD programmable discontinuities formal layer

This directory is a Lean 4 scaffold for the math-first layer of the RAD
simulator. It formalizes statements that are independent of the current
visualization and numerical solver:

- active mechanical modes as Boolean constraint maps, including finite-support
  witnesses and `FiniteMode n` for explicitly finite constraint universes;
- abstract constraints and admissible configuration sets;
- lock and unlock events as mode operations;
- local discontinuity operators with explicit support and finite-support
  witnesses;
- list-supported group operators for simultaneous lock/unlock abstractions;
- reachability in a finite mode graph;
- the scalar backlash dead-zone map used by the RAD preprint and simulator.

The Lean package is intentionally separate from the Python and Three.js code.
The simulator remains responsible for numerical geometry, physical previews,
visualization, inverse design, and calibration workflows. Lean is used only for
small first-principles claims whose premises are explicit enough to prove.

## Build

Install Lean 4 with Lake and run:

```powershell
cd "C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation\formal"
lake build
```

This scaffold pins Lean to `v4.21.0` in `lean-toolchain` and intentionally avoids
mathlib for the first pass so the proof layer stays lightweight.

## Proof targets

The current Lean files define proof terms for:

- `MechanicalSystem.lock_idempotent`
- `MechanicalSystem.unlock_idempotent`
- `MechanicalSystem.lock_shrinks_admissible`
- `MechanicalSystem.unlock_enlarges_admissible`
- `finite_active_constraints`
- `LocalModeOperator.finite_operator_support`
- `LocalModeOperator.disjoint_support_operators_commute`
- `LocalModeOperator.groupOperator_apply_inside`
- `LocalModeOperator.groupOperator_apply_outside`
- `LocalModeOperator.groupOperator_support_inside`
- `LocalModeOperator.groupOperator_support_outside`
- `LocalModeOperator.groupOperator_empty_apply`
- `LocalModeOperator.groupOperator_append_support_left`
- `LocalModeOperator.groupOperator_append_support_right`
- `LocalModeOperator.groupOperators_with_disjoint_lists_commute`
- `identity_compose_left`
- `identity_compose_right`
- `composeOperators_assoc`
- `applyEventSequence_append`
- `CellGraph.removed_cell_not_present`
- `CellGraph.removal_deletes_incident_edges_left`
- `CellGraph.removal_deletes_incident_edges_right`
- `CellGraph.removal_deletes_one_step_from_removed`
- `CellGraph.removal_deletes_one_step_to_removed`
- `CellGraph.removal_preserves_nonincident_one_step`
- `CellGraph.positive_path_source_present`
- `CellGraph.positive_path_target_present`
- `CellGraph.removal_deletes_positive_path_from_removed`
- `CellGraph.removal_deletes_positive_path_to_removed`
- `CellGraph.reachable_source_present`
- `CellGraph.reachable_target_present`
- `CellGraph.removal_deletes_reachable_from_removed`
- `CellGraph.removal_deletes_reachable_to_removed`
- `CellGraph.removed_cell_constraints_inactive`
- `CellGraph.removed_cell_constraints_clear_group_operator`
- `CellGraph.removed_cell_constraints_not_active_after_group_operator`
- `reachable_trans`

These graph-deletion theorems are also the formal reference for the browser
spring-preview diagnostic: removed cells should not act as hidden averaging
nodes, and incident link-strain entries should be treated as absent edges.
- `springEnergyNat_nonnegative`
- `hingeEnergyNat_nonnegative`
- `contactPenaltyNat_nonnegative`
- `clearanceExcessNat_zero_inside`
- `verticalResidualStepNat_zero_inside`
- `contactPenaltyFromClearanceNat_nonnegative`
- `loadWorkMagnitudeNat_nonnegative`
- `loadWorkMagnitudeNat_zero_fixed`
- `signedLoadWorkInt_zero_displacement`
- `signedLoadWorkInt_zero_load`
- `scaledSpringEnergyNat_numerator_nonnegative`
- `scaledSpringEnergyNat_zero_displacement`
- `scaledSpringEnergyNat_preserves_denominator`
- `scaledContactPenaltyNat_numerator_nonnegative`
- `scaledContactPenaltyNat_zero_penetration`
- `scaledSignedLoadWorkInt_zero_displacement`
- `scaledSignedLoadWorkInt_zero_load`
- `scaledSignedLoadWorkInt_preserves_denominator`
- `mechanicalStoredEnergyNat_nonnegative`
- `mechanicalStoredEnergyNat_zero_components`
- `absoluteErrorNat_self`
- `absoluteErrorInt_self`
- `measuredWorkResidualNat_zero_when_equal`
- `signedWorkResidualInt_zero_when_equal`
- `measuredEnergyResidualTripleNat_nonnegative`
- `measuredEnergyResidualTripleNat_zero_when_equal`
- `signedEnergyResidualTripleInt_nonnegative`
- `signedEnergyResidualTripleInt_zero_when_equal`
- `scaledSignedWorkResidualInt_numerator_nonnegative`
- `scaledSignedWorkResidualInt_zero_when_equal`
- `scaledSignedEnergyResidualTripleInt_numerator_nonnegative`
- `scaledSignedEnergyResidualTripleInt_zero_when_equal`
- `measurementUnitScaleNat_zero`
- `measurementUnitScaleNat_preserves_denominator`
- `measurementUnitScaleInt_zero`
- `measurementUnitScaleInt_preserves_denominator`
- `measurementUnitScaleResidualNat_zero_when_equal`
- `measurementUnitScaleResidualNat_preserves_denominator`
- `measurementUnitScaleResidualInt_zero_when_equal`
- `measurementUnitScaleResidualInt_preserves_denominator`
- `hardwareProfileCoverageComplete_true_when_equal`
- `hardwareProfileCoverageMissing_zero_when_complete`
- `hardwareProfileCoverageMissing_preserves_total`
- `calibrationFitResidualPassNat_zero`
- `calibrationFitResidualNat_preserves_sample_count`
- `calibrationModelProfileUpdateSafeNat_intro`
- `calibrationModelProfileUpdateSafeNat_bounds`
- `calibrationModelProfileCandidateImprovesNat_intro`
- `calibrationModelProfileCandidateImprovesNat_applied_updates`
- `calibrationModelProfileCandidateImprovesNat_residual_nonincrease`
- `calibrationModelProfileCandidateImprovesNat_missing_nonincrease`
- `calibrationModelProfileCandidateScoreNat_le_before`
- `calibrationModelProfileHoldoutPassNat_intro`
- `calibrationModelProfileHoldoutPassNat_applied_updates`
- `calibrationModelProfileHoldoutPassNat_fit_residual_nonincrease`
- `calibrationModelProfileHoldoutPassNat_holdout_residual_nonincrease`
- `calibrationModelProfileHoldoutPassNat_missing_nonincrease`
- `calibrationModelProfileHoldoutScoreNat_le_before`
- `calibrationTrainHoldoutSplitReadyNat_intro`
- `calibrationTrainHoldoutSplitReadyNat_has_fit_samples`
- `calibrationTrainHoldoutSplitReadyNat_has_holdout_samples`
- `calibrationTrainHoldoutSplitReadyNat_zero_overlap`
- `calibrationTrainHoldoutSplitReadyNat_profile_frozen`
- `calibrationTrainHoldoutProvenanceReadyNat_intro`
- `calibrationTrainHoldoutProvenanceReadyNat_distinct_files`
- `calibrationTrainHoldoutProvenanceReadyNat_profile_frozen`
- `calibrationTrainHoldoutProvenanceReadyNat_profile_matches`
- `calibrationBenchProtocolCoverageReadyNat_intro`
- `calibrationBenchProtocolCoverageReadyNat_has_scenarios`
- `calibrationBenchProtocolCoverageReadyNat_has_outputs`
- `calibrationBenchProtocolCoverageReadyNat_has_two_dataset_roles`
- `calibrationBenchPacketCompleteNat_intro`
- `calibrationBenchPacketCompleteNat_has_notebook`
- `calibrationBenchPacketCompleteNat_has_fit_template`
- `calibrationBenchPacketCompleteNat_has_holdout_template`
- `calibrationBenchPacketCompleteNat_has_manifest`
- `calibrationBenchExecutedValidationReadyNat_intro`
- `calibrationBenchExecutedValidationReadyNat_has_fit_samples`
- `calibrationBenchExecutedValidationReadyNat_has_holdout_samples`
- `calibrationBenchExecutedValidationReadyNat_has_applied_updates`
- `calibrationBenchExecutedValidationReadyNat_has_independent_validation`
- `calibrationBenchExecutedValidationReadyNat_zero_missing_evidence`
- `physicalValidationReadyNat_intro`
- `physicalValidationReadyNat_has_calibration`
- `physicalValidationReadyNat_has_vertical_load`
- `physicalValidationReadyNat_has_load_proxy`
- `physicalValidationReadyNat_has_contact_proxy`
- `physicalValidationReadyNat_zero_missing_evidence`
- `contactStateAbstractionReadyNat_intro`
- `contactStateAbstractionReadyNat_has_bodies`
- `contactStateAbstractionReadyNat_pins_match_bodies`
- `contactStateAbstractionReadyNat_holes_match_bodies`
- `contactStateAbstractionReadyNat_has_contact_records`
- `contactStateAbstractionReadyNat_has_penalty_terms`
- `contactStateAbstractionReadyNat_zero_missing_evidence`
- `contactParameterCalibrationPacketCompleteNat_intro`
- `contactParameterCalibrationPacketCompleteNat_zero_missing_evidence`
- `contactParameterBenchValidationReadyNat_intro`
- `contactParameterBenchValidationReadyNat_zero_missing_evidence`
- `contactParameterIntervalCalibrationReadyNat_intro`
- `contactParameterIntervalCalibrationReadyNat_zero_missing_evidence`
- `verticalLoadScenarioPassNat_zero_errors`
- `verticalLoadScenarioPassNat_false_when_missing`
- `ModeGraph.finite_mode_graph_vertices`
- `ModeGraph.finite_mode_graph_edges`
- `backlash_zero_inside`
- `backlash_positive_branch`
- `backlash_negative_branch`

A successful `lake build` is the authority that these proof terms typecheck
against the pinned Lean version.

The backlash lemmas currently use `Int` rather than `Real` to avoid a heavy
mathlib dependency. They prove the same max/min dead-zone branches over an
ordered signed scalar type available in Lean's standard library. A real-valued
version is represented by `RealValuedBacklashTargets` and is the next proof step
once mathlib cache/tooling is available.
The mechanics file follows the same discipline: simple `Nat` energy
nonnegativity lemmas are proved now. `mechanicalStoredEnergyNat` records the
finite nonnegative stored/proxy sum used by the simulator energy certificate,
and `measuredWorkResidualNat_zero_when_equal` records the finite validation
identity behind measured-versus-simulated load-work comparisons.
`signedLoadWorkInt_zero_displacement` and `signedLoadWorkInt_zero_load` add a
signed integer work scaffold for the simulator convention `-load *
displacement`, and `signedWorkResidualInt_zero_when_equal` proves the matching
signed residual identity using integer absolute value.
`ScaledNatQuantity` and `ScaledIntQuantity` add explicit denominator-carrying
records. The scaled theorems prove numerator nonnegativity, zero numerators for
zero displacement/penetration/load cases, denominator preservation, and
zero signed residual numerators under equality. This is an integer-scaled
bridge toward rational mechanics, not a proof over Lean's `Rat` or `Real`.
`MeasurementUnitScale` records the next unit-conversion abstraction: normalized
simulator values can be mapped into denominator-carrying physical-unit records
while preserving zero commands, equality residual zeros, and denominator
metadata for both unsigned and signed quantities. These are unit-scale metadata
invariants, not proof that a measured hardware profile is calibrated.
`HardwareProfileCoverage` adds the finite bookkeeping layer used by the
`rad-sim.hardware-profile.v1` artifact: complete measured-field coverage is
recognized when measured and total field counts agree, and missing-field count
is zero in that complete case. This proves metadata consistency, not physical
accuracy of the measured dimensions.
`CalibrationFitResidualNat` adds the finite bookkeeping layer used by
`rad-sim.calibration-parameter-estimates.v1`: a fit with at least one sample and
zero fitted residual passes any finite tolerance, and the sample count is
preserved as metadata. This is a residual-accounting theorem, not proof that
the fitted gain, z-coupling estimate, or contact proxy is a physical law.
`CalibrationModelProfileUpdateNat` adds the finite safety gate used by
`rad-sim.calibration-model-profile.v1`: a profile update is safe when it has a
positive sample count and the proposed bounded parameter lies between declared
lower and upper bounds. This is an update-admissibility theorem, not proof that
the fitted update is a calibrated physical law.
`CalibrationModelProfileCandidateNat` adds the finite selection predicate used
by `rad-sim.calibration-model-profile-selection.v1`: a selectable profile
candidate must have at least one applied update, non-increased residual score,
and non-increased missing-observation count. This proves the bookkeeping rule
used to rank profile candidates, not that the fitted profile is independently
validated hardware physics.
`CalibrationModelProfileHoldoutNat` adds the finite train/holdout predicate
used by `rad-sim.calibration-model-profile-holdout-validation.v1`: both the
fit replay and holdout replay must have non-increased residual scores and
non-increased missing-observation counts after an applied update. This proves
the holdout gate logic, not that a holdout file was experimentally independent
or that contact/friction/material parameters are correct.
`CalibrationTrainHoldoutSplitNat` adds the finite metadata predicate for a
proper split: positive fit samples, positive holdout samples, zero overlap, and
a frozen profile. This proves the split bookkeeping rule once those premises
are supplied; it does not itself verify file provenance or lab procedure.
`CalibrationTrainHoldoutProvenanceNat` adds the finite provenance predicate:
fit and holdout dataset IDs are known, dataset IDs and source-file IDs are
distinct, roles are marked, the profile is frozen, and the holdout file names
the frozen profile. This is still metadata logic; laboratory notebook
discipline must provide the premises.
`CalibrationBenchProtocolCoverageNat` adds the finite protocol-coverage
predicate used by `rad-sim.calibration-bench-notebook.v1`: a notebook is ready
only when it has at least one scenario, instrument, output artifact,
measurement column, and both fit and holdout dataset roles. This proves that
the handoff artifact is structurally populated; it does not prove that the
bench run was performed independently or that any contact/friction law is
calibrated.
`CalibrationBenchPacketCompletenessNat` adds the finite completeness predicate
used by `rad-sim.calibration-bench-packet.v1`: the packet must contain notebook
artifacts, protocol artifacts, fit and holdout templates, validation
instructions, and a populated artifact manifest. This proves bundle
completeness for the handoff folder/JSON; it does not prove that the raw files
were collected correctly or that simulator residual agreement is physical
truth.
`CalibrationBenchExecutedValidationNat` adds the finite evidence gate used by
`rad-sim.calibration-bench-execution-validation.v1`: fit and holdout samples
must be positive, at least one bounded update must be applied, residual and
independent-validation flags must be positive, and missing evidence must be
zero. This proves only the finite bookkeeping predicate once those fields are
supplied; it does not prove laboratory procedure truth, measurement accuracy,
contact mechanics, gravity, friction, stiffness, or material-law calibration.
`PhysicalValidationReadinessNat` composes that calibration gate with the
vertical-load comparison gate used by
`rad-sim.physical-validation-readiness.v1`: executed calibration must be ready,
at least one vertical-load scenario must pass, load-work and contact-proxy terms
must be present, clearance must be configured, and missing evidence must be
zero. This is a physical-claim review gate, not a proof of calibrated rigid-body
contact or gravity.
`ContactStateAbstractionNat` records the finite bookkeeping behind
`rad-sim.contact-state-abstraction.v1`: active bodies must have matching pin,
hole, and clearance records, contact-state and penalty records must cover those
bodies, clearance must be configured, and missing evidence must be zero. This
formalizes the abstract contact-state inventory, not collision geometry,
friction, or calibrated contact stiffness.
`ContactGraphConsistencyNat` records the finite graph/contact bookkeeping behind
`rad-sim.contact-graph-consistency.v1`: active bodies must exist, contact
records must cover active bodies, active graph edges must not touch removed
cells, removed cells must not carry active contact, group-support records must
cover group-support cells, and missing evidence must be zero. The compiled
projection theorems are `contactGraphConsistentNat_has_active_bodies`,
`contactGraphConsistentNat_contact_records_cover_bodies`,
`contactGraphConsistentNat_removed_edges_clear`,
`contactGraphConsistentNat_removed_contacts_clear`,
`contactGraphConsistentNat_support_records_cover_support`, and
`contactGraphConsistentNat_zero_missing_evidence`.
`PhysicalRealizationMapNat` records the finite bridge behind
`rad-sim.physical-realization-map.v1`: abstract operators must be present,
realized simulator state effects, support records, and claim labels must cover
those operators, contact-graph evidence must be ready, and missing evidence must
be zero. The compiled projection theorems are
`physicalRealizationMapReadyNat_has_operators`,
`physicalRealizationMapReadyNat_realizes_all`,
`physicalRealizationMapReadyNat_has_support_records`,
`physicalRealizationMapReadyNat_has_state_effects`,
`physicalRealizationMapReadyNat_has_claim_labels`,
`physicalRealizationMapReadyNat_has_contact_graph`, and
`physicalRealizationMapReadyNat_zero_missing_evidence`.
`ExternalPhysicsEngineAuditNat` records the finite readiness gate behind
`rad-sim.external-physics-engine-audit.v1`: external rigid-body/contact engine
candidates must be present, at least one independent candidate must be
available, required feature and validation-scenario records must exist,
contact-model and independent-tool records must be present, and missing evidence
must be zero. The compiled projection theorems are
`externalPhysicsEngineAuditReadyNat_has_engines`,
`externalPhysicsEngineAuditReadyNat_has_available`,
`externalPhysicsEngineAuditReadyNat_has_features`,
`externalPhysicsEngineAuditReadyNat_has_scenarios`,
`externalPhysicsEngineAuditReadyNat_has_contact_model`,
`externalPhysicsEngineAuditReadyNat_has_independent_tool`, and
`externalPhysicsEngineAuditReadyNat_zero_missing_evidence`.
`ExternalPhysicsModelExportNat`, `ExternalContactGeometryNat`,
`ExternalContactParameterProfileNat`, `ContactParameterCalibrationPacketNat`,
`ContactParameterBenchValidationNat`, `ContactParameterIntervalCalibrationNat`,
`ExternalPhysicsRunNat`, and `ExternalPhysicsComparisonNat` record the finite
MuJoCo handoff gates behind
`rad-sim.mujoco-model-export.v1`,
`rad-sim.mujoco-pin-hole-contact-geometry.v1`,
`rad-sim.mujoco-contact-parameter-profile.v1`,
`rad-sim.contact-parameter-calibration-packet.v1`,
`rad-sim.contact-parameter-bench-validation.v1`,
`rad-sim.contact-parameter-interval-calibration.v1`,
`rad-sim.mujoco-external-run.v1`, and
`rad-sim.mujoco-external-comparison.v1`. The export gate requires model, body,
fixed-body, gravity, and XML-byte records with zero missing evidence. The
contact-geometry gate requires pin, hole, clearance, contact-pair, active-pair,
and XML-fragment records with zero missing evidence. The contact-parameter gate
requires each contact pair to have parameter, friction, solver, stiffness,
damping, and XML-attribute records with zero missing evidence. The calibration
packet gate requires contact-pair records, parameter coverage, measurement
columns, fit and holdout template rows, artifact-manifest entries, attached
profile evidence, and zero missing evidence. The contact bench-validation gate
requires a ready packet, fit rows, holdout rows, completed measurements covering
both, parameter residuals, fit pass, holdout pass, independent holdout pass, and
zero missing evidence. The interval-calibration gate requires a passing bench
validation, parameter intervals, accepted interval records, uncertainty records,
holdout agreement records, simulator parameters inside bounds, and zero missing
evidence. The run gate requires an available engine, a ready model export,
body-result coverage, positive requested steps, completed steps, and zero
missing evidence. The comparison gate requires a
complete run, simulator records, external and matched records covering simulator
records, tolerance evidence, within-tolerance status, and zero missing evidence.
The compiled projection theorem families end in
`externalPhysicsModelExportReadyNat_zero_missing_evidence`,
`externalContactGeometryReadyNat_zero_missing_evidence`,
`externalContactParameterProfileReadyNat_zero_missing_evidence`,
`contactParameterCalibrationPacketCompleteNat_zero_missing_evidence`,
`contactParameterBenchValidationReadyNat_zero_missing_evidence`,
`contactParameterIntervalCalibrationReadyNat_zero_missing_evidence`,
`externalPhysicsRunReadyNat_zero_missing_evidence`, and
`externalPhysicsComparisonReadyNat_zero_missing_evidence`.
`EquilibriumRelationNat` records the finite variational evidence gate behind
`rad-sim.equilibrium-relation.v1`: physical realization evidence, solver
success, and contact-state evidence must be present, energy terms must be
covered by nonnegative terms, the residual must be within tolerance, and missing
evidence must be zero. The compiled projection theorems are
`equilibriumRelationReadyNat_has_realization`,
`equilibriumRelationReadyNat_has_solver`,
`equilibriumRelationReadyNat_has_contact`,
`equilibriumRelationReadyNat_has_energy_terms`,
`equilibriumRelationReadyNat_energy_terms_nonnegative`,
`equilibriumRelationReadyNat_residual_within_tolerance`, and
`equilibriumRelationReadyNat_zero_missing_evidence`.
`ReachableEquilibriumControllabilityNat` records the finite response-basis gate
behind `rad-sim.reachable-equilibrium-controllability.v1`: equilibrium evidence
must be ready, an actuator basis and response columns must exist, target cells
and reachable responses must be recorded, topology evidence must be present,
the optional full-target policy must be satisfied, and missing evidence must be
zero. The compiled projection theorems are
`reachableEquilibriumControllabilityReadyNat_has_equilibrium`,
`reachableEquilibriumControllabilityReadyNat_has_actuator_basis`,
`reachableEquilibriumControllabilityReadyNat_has_response_columns`,
`reachableEquilibriumControllabilityReadyNat_has_target_cells`,
`reachableEquilibriumControllabilityReadyNat_has_reachable_responses`,
`reachableEquilibriumControllabilityReadyNat_has_topology`,
`reachableEquilibriumControllabilityReadyNat_target_policy`, and
`reachableEquilibriumControllabilityReadyNat_zero_missing_evidence`.
`ReachableEquilibriumBenchProtocolNat` records the finite bench-protocol gate
behind `rad-sim.reachable-equilibrium-bench-protocol.v1`: reachability evidence
must be attached, scenarios and actuator-column trials must exist, target
observations and measurement columns must be recorded, the topology-control
policy must be represented, pass/fail criteria must be present, and missing
evidence must be zero. The compiled projection theorems are
`reachableEquilibriumBenchProtocolReadyNat_has_reachability`,
`reachableEquilibriumBenchProtocolReadyNat_has_steps`,
`reachableEquilibriumBenchProtocolReadyNat_has_actuator_trials`,
`reachableEquilibriumBenchProtocolReadyNat_has_targets`,
`reachableEquilibriumBenchProtocolReadyNat_has_measurement_columns`,
`reachableEquilibriumBenchProtocolReadyNat_has_topology_policy`,
`reachableEquilibriumBenchProtocolReadyNat_has_pass_fail`, and
`reachableEquilibriumBenchProtocolReadyNat_zero_missing_evidence`.
`ReachableEquilibriumBenchValidationNat` records the finite filled-result gate
behind `rad-sim.reachable-equilibrium-bench-comparison.v1`: the protocol must
be ready, result rows and target measurements must exist, completed
measurements must cover the target measurements, reachability, topology, and
group-sequence checks must be present, the comparison must pass, and missing
evidence must be zero. The compiled projection theorems are
`reachableEquilibriumBenchValidationReadyNat_has_protocol`,
`reachableEquilibriumBenchValidationReadyNat_has_rows`,
`reachableEquilibriumBenchValidationReadyNat_has_targets`,
`reachableEquilibriumBenchValidationReadyNat_measurements_complete`,
`reachableEquilibriumBenchValidationReadyNat_has_reachability_checks`,
`reachableEquilibriumBenchValidationReadyNat_has_topology_checks`,
`reachableEquilibriumBenchValidationReadyNat_has_group_checks`,
`reachableEquilibriumBenchValidationReadyNat_has_pass`, and
`reachableEquilibriumBenchValidationReadyNat_zero_missing_evidence`.
`ReachableEquilibriumAmplitudeCalibrationNat` records the finite repeated-trial
amplitude gate behind
`rad-sim.reachable-equilibrium-amplitude-calibration.v1`: bench validation must
be ready, repeated-trial groups, amplitude estimates, residual-field cells,
uncertainty bands, topology-leakage bands, and group-sequence residuals must be
present, and missing evidence must be zero. The compiled projection theorems
are `reachableEquilibriumAmplitudeCalibrationReadyNat_has_validation`,
`reachableEquilibriumAmplitudeCalibrationReadyNat_has_repeated_trials`,
`reachableEquilibriumAmplitudeCalibrationReadyNat_has_amplitudes`,
`reachableEquilibriumAmplitudeCalibrationReadyNat_has_residual_field`,
`reachableEquilibriumAmplitudeCalibrationReadyNat_has_uncertainty`,
`reachableEquilibriumAmplitudeCalibrationReadyNat_has_topology_bands`,
`reachableEquilibriumAmplitudeCalibrationReadyNat_has_group_residuals`, and
`reachableEquilibriumAmplitudeCalibrationReadyNat_zero_missing_evidence`.
`ReachableEquilibriumEmpiricalProfileNat` records the finite bounded-profile
gate behind `rad-sim.reachable-equilibrium-empirical-profile.v1`: amplitude
calibration must be ready, bounded and safe proposals must exist, tolerance and
uncertainty records must be present, the holdout hook must be represented, and
missing evidence must be zero. The compiled projection theorems are
`reachableEquilibriumEmpiricalProfileReadyNat_has_amplitude`,
`reachableEquilibriumEmpiricalProfileReadyNat_has_bounded_proposals`,
`reachableEquilibriumEmpiricalProfileReadyNat_has_safe_proposals`,
`reachableEquilibriumEmpiricalProfileReadyNat_has_tolerances`,
`reachableEquilibriumEmpiricalProfileReadyNat_has_uncertainty`,
`reachableEquilibriumEmpiricalProfileReadyNat_has_holdout_hooks`, and
`reachableEquilibriumEmpiricalProfileReadyNat_zero_missing_evidence`.
`ReachableEquilibriumProfileInverseNat` records the finite profile-aware inverse
diagnostic behind `rad-sim.reachable-equilibrium-profile-inverse.v1`: an
empirical profile must be ready, safe proposals and target cells must exist, the
inverse solve must be represented, residual and score records must be present,
profile use must be read-only, and missing evidence must be zero. The compiled
projection theorems are
`reachableEquilibriumProfileInverseReadyNat_has_profile`,
`reachableEquilibriumProfileInverseReadyNat_has_safe_proposals`,
`reachableEquilibriumProfileInverseReadyNat_has_targets`,
`reachableEquilibriumProfileInverseReadyNat_has_solve`,
`reachableEquilibriumProfileInverseReadyNat_has_residuals`,
`reachableEquilibriumProfileInverseReadyNat_has_scores`,
`reachableEquilibriumProfileInverseReadyNat_read_only_profile_use`, and
`reachableEquilibriumProfileInverseReadyNat_zero_missing_evidence`.
`ReachableEquilibriumProfileInverseAcceptanceNat` records the finite acceptance
gate behind `rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1`: the
profile-inverse report must be ready, score records must exist, scaled residual
score and band failures must fit declared limits, actuator count must stay
within budget, the plan must be accepted for preview, and missing evidence must
be zero. The compiled projection theorems are
`reachableEquilibriumProfileInverseAcceptanceReadyNat_has_profile_inverse`,
`reachableEquilibriumProfileInverseAcceptanceReadyNat_has_score_records`,
`reachableEquilibriumProfileInverseAcceptanceReadyNat_score_within_limit`,
`reachableEquilibriumProfileInverseAcceptanceReadyNat_band_failures_within_limit`,
`reachableEquilibriumProfileInverseAcceptanceReadyNat_actuators_within_limit`,
`reachableEquilibriumProfileInverseAcceptanceReadyNat_accepted`, and
`reachableEquilibriumProfileInverseAcceptanceReadyNat_zero_missing_evidence`.
`ReachableEquilibriumProfileInversePreviewPacketNat` records the finite handoff
gate behind
`rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1`: acceptance
must be ready, command records must exist, event records must cover commands,
target and residual records must exist, the packet must be marked read-only, and
missing evidence must be zero. The compiled projection theorems are
`reachableEquilibriumProfileInversePreviewPacketReadyNat_has_acceptance`,
`reachableEquilibriumProfileInversePreviewPacketReadyNat_has_commands`,
`reachableEquilibriumProfileInversePreviewPacketReadyNat_events_cover_commands`,
`reachableEquilibriumProfileInversePreviewPacketReadyNat_has_targets`,
`reachableEquilibriumProfileInversePreviewPacketReadyNat_has_residuals`,
`reachableEquilibriumProfileInversePreviewPacketReadyNat_read_only`, and
`reachableEquilibriumProfileInversePreviewPacketReadyNat_zero_missing_evidence`.
`ReachableEquilibriumProfileInversePreviewReplayNat` records the finite
read-only replay gate behind
`rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1`: the source
packet must be ready, command records must exist and replay exactly, event
coverage must be represented, simulation and target-residual records must
exist, residual agreement must hold, replay must be read-only, and missing
evidence must be zero. The compiled projection theorems are
`reachableEquilibriumProfileInversePreviewReplayReadyNat_has_packet`,
`reachableEquilibriumProfileInversePreviewReplayReadyNat_has_commands`,
`reachableEquilibriumProfileInversePreviewReplayReadyNat_replays_all_commands`,
`reachableEquilibriumProfileInversePreviewReplayReadyNat_has_event_coverage`,
`reachableEquilibriumProfileInversePreviewReplayReadyNat_has_simulation`,
`reachableEquilibriumProfileInversePreviewReplayReadyNat_has_target_residuals`,
`reachableEquilibriumProfileInversePreviewReplayReadyNat_has_residual_agreement`,
`reachableEquilibriumProfileInversePreviewReplayReadyNat_read_only`, and
`reachableEquilibriumProfileInversePreviewReplayReadyNat_zero_missing_evidence`.
`ReachableEquilibriumProfileInversePreviewPhysicalNat` records the finite
physical-preview gate behind
`rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1`: replay
must be ready, the physical-preview solver must succeed, command and target
residual records must exist, model-comparison and energy records must exist,
the preview must be read-only, and missing evidence must be zero. The compiled
projection theorems are
`reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_replay`,
`reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_solver`,
`reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_commands`,
`reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_targets`,
`reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_model_comparison`,
`reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_energy`,
`reachableEquilibriumProfileInversePreviewPhysicalReadyNat_read_only`, and
`reachableEquilibriumProfileInversePreviewPhysicalReadyNat_zero_missing_evidence`.
`measuredEnergyResidualTripleNat_zero_when_equal` extends that identity to the
signed-work, work-magnitude, and contact-proxy error triple used by the bench
comparison report, while `signedEnergyResidualTripleInt_zero_when_equal` records
the analogous signed integer triple. `verticalLoadScenarioPassNat_zero_errors` proves that a
scenario with no missing measurements and zero finite errors passes any
tolerance, while `verticalLoadScenarioPassNat_false_when_missing` proves that
missing measurements fail the finite predicate. These are validation logic
theorems, not calibration theorems for gravity or contact. The
`BenchValidationTargets` structure records the premises needed for a future
real-valued bench-validation model,
while `RealValuedEnergyTargets` records the premises needed for future
real-valued spring, hinge, and contact proofs.

## Not yet formalized

The following simulator features are not Lean-proven physics claims yet:

- CAD-accurate RAD cell geometry;
- spring-hinge equilibrium and solver convergence;
- gravity, contact, friction, pin-hole clearance, and vertical residual motion;
- calibrated stiffness, actuator force, or material behavior;
- inverse design correctness for target surfaces;
- conformal mapping beyond the existing research-grounding notes.

Those require measured premises or a more detailed mechanics model before Lean
can prove statements about them.
