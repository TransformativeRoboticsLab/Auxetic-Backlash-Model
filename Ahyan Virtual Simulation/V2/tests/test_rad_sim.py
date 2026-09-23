import csv
import io
import json
import math
import subprocess
import sys
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rad_sim import (
    DeadZonePropagationOperator,
    CAD_MESH_AUDIT_SCHEMA,
    CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA,
    CAD_RAD_CELL_LAYOUT_SCHEMA,
    CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA,
    ONE_CELL_CAD_EXPORT_AUDIT_SCHEMA,
    HARDWARE_PROFILE_DIMENSIONS,
    HARDWARE_PROFILE_SCHEMA,
    LatticeConfig,
    LatticeState,
    LoadCase,
    calibrate_repeated_row_locks_from_lock_dataset,
    calibrate_single_strand_from_lock_dataset,
    CALIBRATION_SOLVER_GAPS,
    CALIBRATION_SOLVER_TASKS,
    PAPER_RAD_REFERENCE,
    RADHardwareProfile,
    SourceCommand,
    TWO_CELL_ACTUATION_SWEEP_SCHEMA,
    TWO_CELL_BENCH_SCHEMA,
    TWO_CELL_CONNECTOR_MEASUREMENT_COMPARISON_SCHEMA,
    TWO_CELL_CONNECTOR_MEASUREMENT_TEMPLATE_SCHEMA,
    TWO_CELL_CONNECTOR_CONTACT_SCHEMA,
    TWO_CELL_CONNECTOR_CONTACT_SWEEP_SCHEMA,
    TWO_CELL_EXTERNAL_COMPARISON_SCHEMA,
    TWO_CELL_EXTERNAL_FIDELITY_BENCHMARK_SCHEMA,
    TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_APPLICATION_SCHEMA,
    TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_SCHEMA,
    TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA,
    TWO_CELL_EXTERNAL_FIDELITY_MJCF_RUN_SCHEMA,
    TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA,
    TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA,
    TWO_CELL_EXTERNAL_TEMPLATE_SCHEMA,
    WEB_TWO_CELL_EXTERNAL_FIDELITY_SCHEMA,
    TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA,
    TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA,
    TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA,
    TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA,
    TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA,
    TWO_CELL_CONTACT_PHASE_MAP_SCHEMA,
    TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA,
    TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA,
    TWO_CELL_MEASUREMENT_COMPARISON_SCHEMA,
    TWO_CELL_MJCF_SCHEMA,
    TWO_CELL_MUJOCO_RUN_SCHEMA,
    TWO_CELL_PACKET_SCHEMA,
    TWO_CELL_PARAMETER_CALIBRATION_SCHEMA,
    TWO_CELL_PHYSICS_VALIDATION_SCHEMA,
    TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA,
    TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA,
    TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA,
    TWO_CELL_QUASISTATIC_SCHEMA,
    TWO_CELL_QUASISTATIC_SWEEP_SCHEMA,
    TWO_CELL_SEGMENTED_CAD_INTAKE_TEMPLATES_SCHEMA,
    TWO_CELL_SEGMENTED_CAD_INTAKE_VALIDATION_SCHEMA,
    TWO_CELL_SEGMENTED_CAD_READINESS_SCHEMA,
    TWO_CELL_SWEEP_SCHEMA,
    TopologyExperimentScenario,
    TwoCellBenchControls,
    TwoCellConnectorMeasurement,
    TwoCellMeasurement,
    apply_event_sequence,
    backlash_activation,
    normalized_backlash_to_alpha_dead_zone,
    normalized_backlash_to_theta_dead_zone,
    active_neighbor_edges,
    build_calibration_experiment_protocol,
    build_response_atlas,
    build_response_matrix,
    build_paper_rad_cell_geometry,
    build_paper_rad_lattice_geometry,
    build_paper_rad_lattice_mesh,
    build_topology_experiment_report,
    build_vertical_load_physical_preview_report,
    cad_mesh_audit,
    cad_rad_cell_archive_audit,
    cad_rad_cell_layout,
    cad_rad_cell_reference_profile,
    one_cell_cad_export_audit,
    calibrate_paper_rad_config,
    calibration_experiment_measurements_from_json,
    calibration_experiment_comparison_report,
    calibration_model_profile_from_report,
    calibration_model_profile_holdout_validation,
    calibration_model_profile_residual_comparison,
    calibration_experiment_results_template,
    apply_calibration_model_profile,
    select_calibration_model_profile,
    physical_validation_readiness_report,
    contact_state_abstraction_report,
    contact_graph_consistency_report,
    physical_realization_map_report,
    external_physics_engine_audit_report,
    mujoco_model_export_report,
    mujoco_pin_hole_contact_geometry_report,
    mujoco_contact_parameter_report,
    contact_parameter_calibration_packet,
    contact_parameter_calibration_results_template,
    contact_parameter_calibration_results_from_json,
    compare_contact_parameter_calibration_results,
    contact_parameter_interval_calibration_report,
    mujoco_external_run_report,
    mujoco_external_comparison_report,
    equilibrium_relation_report,
    reachable_equilibrium_controllability_report,
    reachable_equilibrium_bench_protocol,
    reachable_equilibrium_bench_results_template,
    reachable_equilibrium_bench_results_from_json,
    compare_reachable_equilibrium_bench_results,
    reachable_equilibrium_amplitude_calibration_report,
    reachable_equilibrium_empirical_profile_from_amplitude,
    calibration_measurement_plan,
    calibration_readiness,
    characterize_cluster,
    characterize_pair,
    characterize_pairwise_interactions,
    characterize_response,
    characterize_single_cell,
    compare_physical_cluster,
    compare_physical_pair,
    compare_physical_response,
    compare_removed_topology_reachability,
    compare_calibration_experiment_measurements,
    compare_event_order,
    compare_group_actuation_decomposition,
    compare_group_actuation_under_removal,
    compare_vertical_residual_spring_hinge_3d,
    compare_vertical_residual_under_removal,
    compare_vertical_load_energy_measurement_results,
    compare_sequence_order,
    compare_two_cell_external_results,
    compare_two_cell_fidelity_matrix_measurements,
    compare_two_cell_radius_backlash_transition_measurements,
    calibrate_two_cell_fidelity_matrix_parameters,
    compare_two_cell_connector_measurements,
    compare_two_cell_measurements,
    calibration_bench_notebook,
    calibration_bench_packet,
    calibration_bench_execution_validation,
    config_with_hardware_profile,
    clear_actuation_event,
    diagnose_programmable_discontinuity,
    evaluate_programmable_operators,
    export_vertical_load_bench_packet_json,
    export_cad_mesh_audit_json,
    export_cad_rad_cell_archive_audit_json,
    export_cad_rad_cell_layout_json,
    export_cad_rad_cell_reference_profile_json,
    export_one_cell_cad_export_audit_json,
    write_vertical_load_bench_packet_artifacts,
    write_calibration_bench_packet_artifacts,
    write_calibration_bench_execution_validation_artifacts,
    write_two_cell_bench_packet_artifacts,
    write_two_cell_external_comparison_artifacts,
    write_two_cell_measurement_comparison_artifacts,
    write_vertical_load_energy_comparison_artifacts,
    finite_die_off_radius,
    export_calibration_experiment_protocol_json,
    export_calibration_bench_notebook_json,
    export_calibration_bench_notebook_csv,
    export_calibration_bench_packet_json,
    export_calibration_bench_execution_validation_csv,
    export_calibration_bench_execution_validation_json,
    export_calibration_experiment_comparison_json,
    export_calibration_experiment_comparison_report_json,
    export_calibration_model_profile_json,
    export_calibration_model_profile_holdout_validation_csv,
    export_calibration_model_profile_holdout_validation_json,
    export_calibration_model_profile_residual_comparison_json,
    export_calibration_model_profile_selection_json,
    export_calibration_experiment_results_template_json,
    export_inverse_design_report_json,
    export_reachable_equilibrium_profile_inverse_acceptance_csv,
    export_reachable_equilibrium_profile_inverse_acceptance_json,
    export_reachable_equilibrium_profile_inverse_csv,
    export_reachable_equilibrium_profile_inverse_json,
    export_reachable_equilibrium_profile_inverse_preview_packet_csv,
    export_reachable_equilibrium_profile_inverse_preview_packet_json,
    export_reachable_equilibrium_profile_inverse_preview_physical_csv,
    export_reachable_equilibrium_profile_inverse_preview_physical_json,
    export_reachable_equilibrium_profile_inverse_preview_replay_csv,
    export_reachable_equilibrium_profile_inverse_preview_replay_json,
    export_programmable_discontinuity_report_json,
    export_formalization_target_manifest_json,
    export_hardware_profile_json,
    export_lock_coordinate_dataset_json,
    export_lock_coordinate_mlp_json,
    export_lock_dataset_summary_json,
    export_response_atlas_json,
    export_response_atlas_sweep_json,
    export_response_matrix_json,
    export_removed_topology_reachability_comparison_json,
    export_topology_experiment_report_json,
    export_two_cell_actuation_sweep_csv,
    export_two_cell_actuation_sweep_json,
    export_two_cell_backlash_sweep_csv,
    export_two_cell_bench_json,
    export_two_cell_cad_contact_decomposition_csv,
    export_two_cell_cad_contact_decomposition_json,
    export_two_cell_connector_contact_csv,
    export_two_cell_connector_contact_json,
    export_two_cell_connector_contact_sweep_csv,
    export_two_cell_connector_contact_sweep_json,
    export_two_cell_connector_measurement_comparison_csv,
    export_two_cell_connector_measurement_comparison_json,
    export_two_cell_connector_measurement_template_csv,
    export_two_cell_connector_measurement_template_json,
    export_two_cell_external_comparison_json,
    export_two_cell_external_fidelity_matrix_manifest_csv,
    export_two_cell_external_fidelity_matrix_manifest_json,
    export_two_cell_external_fidelity_mjcf_measurements_csv,
    export_two_cell_external_fidelity_mjcf_run_json,
    export_two_cell_exact_contact_handoff_plan_csv,
    export_two_cell_exact_contact_handoff_plan_json,
    export_two_cell_external_results_template_csv,
    export_two_cell_external_results_template_json,
    export_two_cell_fidelity_matrix_measurement_comparison_csv,
    export_two_cell_fidelity_matrix_measurement_comparison_json,
    export_two_cell_fidelity_matrix_parameter_calibration_json,
    export_two_cell_radius_backlash_transition_comparison_csv,
    export_two_cell_radius_backlash_transition_comparison_json,
    export_two_cell_radius_backlash_transition_rerun_csv,
    export_two_cell_radius_backlash_transition_rerun_json,
    export_two_cell_external_fidelity_correction_application_csv,
    export_two_cell_external_fidelity_correction_application_json,
    export_two_cell_external_fidelity_correction_profile_json,
    export_two_cell_fidelity_matrix_measurement_template_csv,
    export_two_cell_fidelity_matrix_measurement_template_json,
    export_two_cell_measurement_template_csv,
    export_two_cell_measurement_comparison_json,
    export_two_cell_mjcf_proxy_json,
    export_two_cell_mjcf_xml,
    export_two_cell_mujoco_proxy_run_json,
    export_two_cell_physical_fidelity_matrix_csv,
    export_two_cell_physical_fidelity_matrix_json,
    export_two_cell_contact_phase_map_csv,
    export_two_cell_contact_phase_map_json,
    export_two_cell_radius_backlash_phase_diagram_csv,
    export_two_cell_radius_backlash_phase_diagram_json,
    export_two_cell_radius_backlash_transition_report_csv,
    export_two_cell_radius_backlash_transition_report_json,
    export_two_cell_physical_response_atlas_csv,
    export_two_cell_physical_response_atlas_json,
    export_two_cell_parameter_calibration_json,
    export_two_cell_physical_test_packet_json,
    export_two_cell_physical_simulation_suite_csv,
    export_two_cell_physical_simulation_suite_json,
    export_two_cell_physics_validation_report_csv,
    export_two_cell_physics_validation_report_json,
    export_two_cell_quasistatic_json,
    export_two_cell_quasistatic_sweep_csv,
    export_two_cell_quasistatic_sweep_json,
    export_two_cell_segmented_cad_intake_templates_csv,
    export_two_cell_segmented_cad_intake_templates_json,
    export_two_cell_segmented_cad_intake_validation_csv,
    export_two_cell_segmented_cad_intake_validation_json,
    export_two_cell_segmented_cad_readiness_csv,
    export_two_cell_segmented_cad_readiness_json,
    export_physical_unit_scale_metadata_json,
    export_physical_validation_readiness_csv,
    export_physical_validation_readiness_json,
    export_contact_state_abstraction_csv,
    export_contact_state_abstraction_json,
    export_contact_graph_consistency_csv,
    export_contact_graph_consistency_json,
    export_physical_realization_map_csv,
    export_physical_realization_map_json,
    export_external_physics_engine_audit_csv,
    export_external_physics_engine_audit_json,
    export_mujoco_model_xml,
    export_mujoco_model_report_json,
    export_mujoco_pin_hole_contact_geometry_csv,
    export_mujoco_pin_hole_contact_geometry_json,
    export_mujoco_contact_parameter_csv,
    export_mujoco_contact_parameter_json,
    export_contact_parameter_calibration_packet_csv,
    export_contact_parameter_calibration_packet_json,
    export_contact_parameter_calibration_results_template_json,
    export_contact_parameter_bench_validation_csv,
    export_contact_parameter_bench_validation_json,
    export_contact_parameter_interval_calibration_csv,
    export_contact_parameter_interval_calibration_json,
    export_mujoco_external_run_json,
    export_mujoco_external_comparison_csv,
    export_mujoco_external_comparison_json,
    export_equilibrium_relation_csv,
    export_equilibrium_relation_json,
    export_reachable_equilibrium_controllability_csv,
    export_reachable_equilibrium_controllability_json,
    export_reachable_equilibrium_bench_protocol_csv,
    export_reachable_equilibrium_bench_protocol_json,
    export_reachable_equilibrium_bench_results_template_json,
    export_reachable_equilibrium_bench_comparison_csv,
    export_reachable_equilibrium_bench_comparison_json,
    export_reachable_equilibrium_amplitude_calibration_csv,
    export_reachable_equilibrium_amplitude_calibration_json,
    export_reachable_equilibrium_empirical_profile_csv,
    export_reachable_equilibrium_empirical_profile_json,
    export_vertical_load_energy_comparison_report_json,
    export_vertical_load_energy_experiment_protocol_json,
    export_vertical_load_energy_measurement_template_json,
    export_vertical_load_energy_validation_json,
    export_vertical_load_physical_preview_report_json,
    export_vertical_load_physical_preview_report_csv,
    export_vertical_removal_physical_comparison_json,
    hardware_profile_from_config,
    hardware_profile_from_json,
    group_actuation_event,
    inverse_design_report,
    reachable_equilibrium_profile_inverse_acceptance_report,
    reachable_equilibrium_profile_inverse_preview_packet,
    reachable_equilibrium_profile_inverse_preview_physical_report,
    reachable_equilibrium_profile_inverse_preview_replay_report,
    reachable_equilibrium_profile_inverse_report,
    lattice_topology_diagnostic,
    local_actuation_event,
    lock_event,
    lock_coordinate_dataset,
    lock_coordinate_training_arrays,
    lock_dataset_summary,
    load_lock_dataset,
    parse_lock_filename,
    predict_lock_coordinate_mlp,
    predict_lock_strand_shape,
    lock_projection,
    mechanics_energy_certificate_to_dict,
    model_provenance,
    physical_unit_scale_metadata,
    provenance_by_status,
    provenance_summary,
    programmable_discontinuity_report,
    formalization_target_manifest,
    export_paper_rad_mesh_obj,
    iter_obj_vertices,
    mesh_bounds_from_file,
    release_event,
    realized_lock_cells,
    remove_cell_event,
    restore_cell_event,
    run_calibration_experiment_protocol,
    response_decay_profile,
    simulate_kinematic,
    simulate_two_cell_bench,
    solve_two_cell_quasistatic,
    solve_inverse_design,
    solve_spring_hinge,
    solve_spring_hinge_3d,
    sweep_two_cell_backlash,
    sweep_two_cell_connector_contact,
    sweep_two_cell_quasistatic,
    sweep_response_atlas_parameters,
    two_cell_actuation_sweep,
    two_cell_cad_contact_decomposition_spec,
    two_cell_connector_contact_report,
    two_cell_connector_measurement_template,
    two_cell_connector_measurement_template_rows,
    two_cell_connector_measurements_from_csv,
    two_cell_connector_measurements_from_json,
    two_cell_external_results_from_csv,
    two_cell_external_results_from_json,
    two_cell_external_fidelity_matrix_manifest,
    two_cell_external_fidelity_matrix_manifest_rows,
    two_cell_external_fidelity_benchmark_summary,
    apply_two_cell_external_fidelity_correction_profile,
    build_two_cell_external_fidelity_web_summary,
    fit_two_cell_external_fidelity_correction_profile,
    two_cell_external_fidelity_mjcf_measurement_rows,
    two_cell_external_fidelity_mjcf_measurements_csv,
    two_cell_external_fidelity_mjcf_run_report,
    two_cell_exact_contact_handoff_plan,
    two_cell_external_results_template,
    two_cell_external_results_template_rows,
    two_cell_fidelity_matrix_measurement_template,
    two_cell_fidelity_matrix_measurement_template_rows,
    two_cell_fidelity_matrix_measurements_from_csv,
    two_cell_fidelity_matrix_measurements_from_json,
    two_cell_radius_backlash_transition_measurements_from_csv,
    two_cell_radius_backlash_transition_rerun_report,
    two_cell_radius_backlash_transition_measurements_from_json,
    two_cell_mjcf_proxy_report,
    two_cell_mujoco_proxy_run_report,
    two_cell_physical_fidelity_matrix,
    two_cell_contact_phase_map,
    two_cell_radius_backlash_phase_diagram,
    two_cell_radius_backlash_transition_report,
    two_cell_physical_response_atlas,
    two_cell_physical_test_packet,
    two_cell_physical_simulation_suite,
    two_cell_physics_validation_report,
    two_cell_segmented_cad_intake_templates,
    two_cell_segmented_cad_intake_validation_report,
    two_cell_segmented_cad_readiness_report,
    two_cell_measurement_template_rows,
    two_cell_measurements_from_csv,
    two_cell_measurements_from_json,
    write_two_cell_connector_measurement_comparison_artifacts,
    write_two_cell_external_fidelity_benchmark_summary_artifact,
    write_two_cell_external_fidelity_correction_application_artifacts,
    write_two_cell_external_fidelity_correction_profile_artifact,
    write_two_cell_external_fidelity_web_summary,
    write_two_cell_fidelity_matrix_measurement_comparison_artifacts,
    write_two_cell_radius_backlash_transition_comparison_artifacts,
    calibrate_two_cell_parameters,
    train_lock_coordinate_mlp,
    validate_inverse_design_physical,
    vertical_load_bench_packet,
    validate_vertical_load_energy_measurements,
    vertical_load_energy_experiment_protocol,
    vertical_load_energy_measurement_results_from_json,
    vertical_load_energy_measurement_template,
    vertical_removal_physical_comparison_to_dict,
    vertical_clearance_operator,
    pair_distance_from_theta,
    solve_constraint_kinematic,
)
from rad_sim.coupling import alpha_to_theta, theta_to_alpha
from rad_sim.constraint_kinematic_reports import constraint_solver_scenario_reports
from rad_sim.spring_hinge import hinge_energy, spring_energy


class RadSimTests(unittest.TestCase):
    def test_backlash_activation_dead_zone(self):
        x = np.array([-0.2, -0.05, 0.0, 0.05, 0.2])
        y = backlash_activation(x, 0.1)
        np.testing.assert_allclose(y, [-0.1, 0.0, 0.0, 0.0, 0.1])

    def test_angle_alpha_mapping(self):
        alpha = np.array([0.5, 1.0, 1.5])
        theta = alpha_to_theta(alpha)
        np.testing.assert_allclose(theta, [-25.0, 10.0, 45.0])
        np.testing.assert_allclose(theta_to_alpha(theta), alpha)

    def test_normalized_backlash_uses_research_angle_conversion(self):
        theta_dead_zone = normalized_backlash_to_theta_dead_zone(0.1)
        self.assertAlmostEqual(theta_dead_zone, math.degrees(math.asin(0.1)))
        self.assertAlmostEqual(
            normalized_backlash_to_alpha_dead_zone(0.1),
            theta_dead_zone / 70.0,
        )

    def test_constraint_solver_satisfies_single_edge_distance_equation(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.0, coupling_gain=0.0)
        state = LatticeState.uniform(config)
        result = solve_constraint_kinematic(config, state)
        distance = np.linalg.norm(
            result.deformed_centers_3d[0, 1] - result.deformed_centers_3d[0, 0]
        )
        desired = pair_distance_from_theta(
            config, result.theta_degrees[0, 0], result.theta_degrees[0, 1]
        )
        self.assertTrue(result.metadata["success"])
        self.assertAlmostEqual(distance, desired, places=4)

    def test_constraint_solver_deletes_edges_incident_to_removed_cells(self):
        config = LatticeConfig(rows=1, cols=3)
        state = LatticeState.uniform(config)
        state.removed_mask[0, 1] = True
        self.assertEqual(active_neighbor_edges(config, state.normalized(config)), [])
        result = solve_constraint_kinematic(config, state)
        self.assertEqual(result.metadata["edge_count"], 0)

    def test_constraint_solver_position_locked_cells_hold_xy(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.0, coupling_gain=0.0)
        state = LatticeState.uniform(config)
        state.position_locked_mask[0, 0] = True
        result = solve_constraint_kinematic(config, state)
        np.testing.assert_allclose(
            result.deformed_centers_3d[0, 0, :2],
            result.original_centers_3d[0, 0, :2],
            atol=1e-5,
        )

    def test_constraint_solver_fixed_end_strand_lifts_when_xy_compressed(self):
        config = LatticeConfig(rows=1, cols=15, backlash=0.0, coupling_gain=0.0)
        state = LatticeState.uniform(config)
        state.position_locked_mask[0, 0] = True
        state.position_locked_mask[0, -1] = True
        result = solve_constraint_kinematic(config, state)
        self.assertTrue(result.metadata["success"])
        self.assertGreater(np.max(np.abs(result.deformed_centers_3d[..., 2])), 0.01)

    def test_constraint_solver_scenario_reports_cover_required_cases(self):
        report = constraint_solver_scenario_reports()
        self.assertEqual(report["schema"], "rad-sim.constraint-kinematic-scenario-report.v1")
        sizes = {(scenario["rows"], scenario["cols"]) for scenario in report["scenarios"]}
        for size in [(1, 2), (1, 15), (15, 1), (3, 3), (15, 15)]:
            self.assertIn(size, sizes)
        for scenario in report["scenarios"]:
            self.assertTrue(np.isfinite(scenario["maxEdgeError"]))
            self.assertTrue(np.isfinite(scenario["maxAbsHeight"]))

    def test_neighbor_die_off_increases_with_backlash(self):
        state_low = LatticeState.uniform(LatticeConfig(rows=5, cols=5, backlash=0.02))
        state_high = LatticeState.uniform(LatticeConfig(rows=5, cols=5, backlash=0.2))
        state_low.actuator_grid[2, 2] = -0.3
        state_high.actuator_grid[2, 2] = -0.3
        low = simulate_kinematic(LatticeConfig(rows=5, cols=5, backlash=0.02), state_low)
        high = simulate_kinematic(LatticeConfig(rows=5, cols=5, backlash=0.2), state_high)
        low_far = abs(low.alpha[0, 0] - 1.0)
        high_far = abs(high.alpha[0, 0] - 1.0)
        self.assertGreaterEqual(low_far, high_far)

    def test_uniform_alpha_preserves_symmetry(self):
        config = LatticeConfig(rows=4, cols=4)
        state = LatticeState.uniform(config, alpha=1.2)
        result = simulate_kinematic(config, state)
        xs = result.deformed_centers[..., 0]
        ys = result.deformed_centers[..., 1]
        self.assertAlmostEqual(abs(xs.min()), abs(xs.max()))
        self.assertAlmostEqual(abs(ys.min()), abs(ys.max()))

    def test_locked_cell_remains_fixed(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        state = LatticeState.uniform(config)
        state.locked_mask[1, 1] = True
        state.actuator_grid[1, 1] = -0.5
        result = simulate_kinematic(config, state)
        self.assertAlmostEqual(result.alpha[1, 1], 1.0)

    def test_lock_dataset_calibrates_single_strand_and_repeated_rows(self):
        logs = Path(r"C:\Users\ahyan\OneDrive\Desktop\logs")
        if not (logs / "marker_xyz_means_all_37_files.csv").exists():
            self.skipTest("local lock dataset is not available")

        parsed = parse_lock_filename("AOM_L3_1_2_3_303030_2.csv")
        self.assertEqual(parsed, (3, (1, 2, 3), (30, 30, 30), 2))
        self.assertEqual(realized_lock_cells((1, 2, 3), (2,)), (1, 3))

        dataset = load_lock_dataset(logs)
        summary = lock_dataset_summary(dataset)
        self.assertEqual(summary["schema"], "rad-sim.single-strand-lock-dataset-summary.v1")
        self.assertEqual(dataset.configuration_count, 37)
        self.assertEqual(summary["stableStateHistogram"], {3: 3, 2: 32, 1: 2})
        self.assertGreater(len(summary["failedConsecutiveLockRecords"]), 0)

        prediction = predict_lock_strand_shape(
            dataset,
            (1, 2, 3),
            (30, 30, 30),
            state_index=2,
            method="nearest",
        )
        self.assertEqual(prediction["schema"], "rad-sim.single-strand-lock-prediction.v1")
        self.assertEqual(prediction["failedLockCells"], (2,))
        self.assertEqual(prediction["realizedLockCells"], (1, 3))
        self.assertEqual(np.asarray(prediction["predictedMarkers"]).shape[1], 3)
        edge_prediction = predict_lock_strand_shape(
            dataset,
            (10, 11, 12),
            (30, 30, 30),
            state_index=2,
            method="nearest",
        )
        self.assertEqual(edge_prediction["failedLockCells"], (10, 11))
        self.assertEqual(edge_prediction["realizedLockCells"], (12,))

        config = LatticeConfig(rows=1, cols=12)
        state = LatticeState.uniform(config)
        state.locked_mask[0, [0, 1, 2]] = True
        state.position_locked_mask[0, [0, 11]] = True
        result = simulate_kinematic(config, state)
        overlay = calibrate_single_strand_from_lock_dataset(
            result,
            dataset,
            state_index=2,
            method="nearest",
        )
        self.assertEqual(
            overlay["schema"],
            "rad-sim.single-strand-lock-calibration-overlay.v1",
        )
        self.assertEqual(np.asarray(overlay["measuredMarkers"]).shape, (12, 3))
        self.assertTrue(np.isfinite(overlay["metrics"]["rmsResidual"]))
        self.assertTrue(overlay["metrics"]["endpointAnchored"])
        np.testing.assert_allclose(
            np.asarray(overlay["measuredMarkers"])[[0, -1]],
            np.asarray(overlay["simulatedCenters"])[[0, -1]],
            atol=1e-9,
        )

        sheet_config = LatticeConfig(rows=3, cols=12)
        sheet_state = LatticeState.uniform(sheet_config)
        sheet_state.locked_mask[:, [0, 1, 2]] = True
        sheet_state.position_locked_mask[:, [0, 11]] = True
        sheet_result = simulate_kinematic(sheet_config, sheet_state)
        sheet_overlay = calibrate_repeated_row_locks_from_lock_dataset(
            sheet_result,
            dataset,
            state_index=2,
            method="nearest",
        )
        calibrated = np.asarray(sheet_overlay["calibratedSheetCenters"])
        self.assertEqual(
            sheet_overlay["schema"],
            "rad-sim.repeated-row-lock-calibration-overlay.v1",
        )
        np.testing.assert_allclose(calibrated[0, :, 0], calibrated[1, :, 0])
        np.testing.assert_allclose(calibrated[0, :, 2], calibrated[2, :, 2])
        self.assertEqual(sheet_overlay["metrics"]["endpointAnchoredRows"], 3)
        np.testing.assert_allclose(
            calibrated[:, [0, -1], :],
            sheet_result.deformed_centers_3d[:, [0, -1], :],
            atol=1e-9,
        )
        self.assertIn("markerCount", export_lock_dataset_summary_json(dataset))

        features, targets = lock_coordinate_training_arrays(dataset)
        self.assertEqual(features.shape, (37, 49))
        self.assertEqual(targets.shape, (37, 36))
        coordinate_payload = lock_coordinate_dataset(dataset)
        self.assertEqual(coordinate_payload["schema"], "rad-sim.lock-coordinate-dataset.v1")
        self.assertEqual(coordinate_payload["outputDimension"], 36)
        self.assertIn("cellCoordinates", coordinate_payload["records"][0])
        self.assertIn("cellCoordinates", export_lock_coordinate_dataset_json(dataset))

        model = train_lock_coordinate_mlp(
            dataset,
            hidden_units=8,
            epochs=40,
            learning_rate=0.01,
            seed=11,
        )
        self.assertEqual(model["schema"], "rad-sim.lock-coordinate-mlp.v1")
        self.assertEqual(model["outputDimension"], 36)
        self.assertTrue(np.isfinite(model["metrics"]["trainMseStandardized"]))
        mlp_prediction = predict_lock_coordinate_mlp(
            model,
            (1, 2, 3),
            (30, 30, 30),
            state_index=2,
        )
        self.assertEqual(
            mlp_prediction["schema"],
            "rad-sim.lock-coordinate-mlp-prediction.v1",
        )
        self.assertEqual(np.asarray(mlp_prediction["cellCoordinates"]).shape, (12, 3))
        self.assertIn("lock-coordinate-mlp", export_lock_coordinate_mlp_json(model))

    def test_vertical_actuation_leaves_neighbor_residual(self):
        config = LatticeConfig(rows=5, cols=5, backlash=0.02, z_coupling_gain=0.35)
        state = LatticeState.uniform(config)
        state.z_actuator_grid[2, 2] = 0.4
        result = simulate_kinematic(config, state)
        height = result.metadata["height"]
        residual = result.metadata["z_residual"]
        self.assertAlmostEqual(height[2, 2], 0.4)
        self.assertGreater(height[2, 3], 0.0)
        self.assertLess(abs(height[2, 3]), abs(height[2, 2]))
        self.assertAlmostEqual(height[2, 3], residual[2, 3])
        self.assertAlmostEqual(result.metadata["z_dead_zone"], config.pin_hole_clearance)
        upward = characterize_response(config, (SourceCommand((2, 2), z=0.4),))
        self.assertGreater(upward.positive_z_reach, 0)
        self.assertEqual(upward.negative_z_reach, 0)
        self.assertGreater(upward.max_positive_height_delta, 0.0)
        self.assertAlmostEqual(upward.max_negative_height_delta, 0.0)
        downward = characterize_response(config, (SourceCommand((2, 2), z=-0.4),))
        self.assertGreater(downward.negative_z_reach, 0)
        self.assertEqual(downward.positive_z_reach, 0)
        self.assertLess(downward.max_negative_height_delta, 0.0)
        self.assertAlmostEqual(downward.max_positive_height_delta, 0.0)

    def test_pin_hole_clearance_controls_vertical_residual_dieoff(self):
        tight = LatticeConfig(
            rows=5,
            cols=5,
            z_coupling_gain=0.35,
            pin_radius=0.18,
            hole_radius=0.20,
        )
        loose = LatticeConfig(
            rows=5,
            cols=5,
            z_coupling_gain=0.35,
            pin_radius=0.18,
            hole_radius=0.32,
        )
        tight_state = LatticeState.uniform(tight)
        loose_state = LatticeState.uniform(loose)
        tight_state.z_actuator_grid[2, 2] = 0.4
        loose_state.z_actuator_grid[2, 2] = 0.4
        tight_result = simulate_kinematic(tight, tight_state)
        loose_result = simulate_kinematic(loose, loose_state)
        self.assertGreater(
            tight_result.metadata["z_residual"][2, 3],
            loose_result.metadata["z_residual"][2, 3],
        )

    def test_pin_hole_radius_validation(self):
        with self.assertRaises(ValueError):
            LatticeConfig(pin_radius=0.2, hole_radius=0.1)

    def test_cad_rad_cell_layout_matches_a360_envelope_and_connector_sites(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        layout = cad_rad_cell_layout(config)
        self.assertEqual(layout["schema"], CAD_RAD_CELL_LAYOUT_SCHEMA)
        self.assertEqual(layout["cadReference"]["nominalHoleDiameterMm"], 3.4)
        self.assertIn("STEP", layout["cadReference"]["viewerAccess"]["downloadFormatsObserved"])
        self.assertIn("outer circular pin-hole pads", layout["cadReference"]["visibleFeatures"])
        self.assertAlmostEqual(layout["visualModel"]["holeToPadRadiusRatio"], 1.7 / 3.145)
        self.assertGreater(layout["visualModel"]["pinToHoleRadiusRatio"], 0.0)
        self.assertEqual(len(layout["padSites"]), 8)
        self.assertEqual(
            [site["name"] for site in layout["padSites"]],
            ["n", "ne", "e", "se", "s", "sw", "w", "nw"],
        )
        self.assertTrue(layout["envelopeCheck"]["matchesBoundingBox"])
        self.assertAlmostEqual(
            layout["envelopeCheck"]["maxOuterEnvelopeMm"],
            layout["envelopeCheck"]["halfWidthMm"],
        )
        self.assertEqual(
            [pair["name"] for pair in layout["twoCellConnectorPairs"]],
            ["upper", "middle", "lower"],
        )
        self.assertEqual(layout["twoCellConnectorPairs"][1]["leftSite"], "e")
        self.assertEqual(layout["twoCellConnectorPairs"][1]["rightSite"], "w")
        self.assertAlmostEqual(layout["dimensionsMm"]["holeRadiusMm"], 1.7)
        self.assertAlmostEqual(layout["dimensionsMm"]["pinRadiusMm"], 1.36)
        self.assertGreater(layout["visualModel"]["siteRadiusToWidth"], 0.4)

        wider_hole = cad_rad_cell_layout(config, hole_radius=0.5)
        self.assertGreater(wider_hole["dimensionsMm"]["holeRadiusMm"], layout["dimensionsMm"]["holeRadiusMm"])
        self.assertGreater(
            wider_hole["dimensionsMm"]["pinHoleClearanceMm"],
            layout["dimensionsMm"]["pinHoleClearanceMm"],
        )
        exported = json.loads(export_cad_rad_cell_layout_json(config))
        self.assertEqual(exported["schema"], CAD_RAD_CELL_LAYOUT_SCHEMA)

    def test_cad_rad_cell_reference_profile_packages_a360_dimensions_for_ui_and_packets(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        profile = cad_rad_cell_reference_profile(config)
        self.assertEqual(profile["schema"], CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA)
        self.assertEqual(profile["cadLayout"]["schema"], CAD_RAD_CELL_LAYOUT_SCHEMA)
        hardware = profile["hardwareProfile"]
        self.assertEqual(hardware["schema"], HARDWARE_PROFILE_SCHEMA)
        self.assertEqual(hardware["name"], "a360-rads-unit-cell-visual-profile")
        self.assertAlmostEqual(hardware["dimensionsMm"]["sideLengthMm"], 55.604331129396634)
        self.assertAlmostEqual(hardware["dimensionsMm"]["holeRadiusMm"], 1.7)
        self.assertAlmostEqual(hardware["dimensionsMm"]["pinRadiusMm"], 1.36)
        self.assertAlmostEqual(hardware["dimensionsMm"]["plateThicknessMm"], 4.0)
        self.assertAlmostEqual(hardware["dimensionsMm"]["jointStackHeightMm"], 19.99999621152464)
        self.assertIn("backlashMm", profile["assumptionLabels"])
        self.assertIn("exact real-cell dynamics", profile["claimBoundary"]["blockedClaim"])
        exported = json.loads(export_cad_rad_cell_reference_profile_json(config))
        self.assertEqual(exported["schema"], CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA)

    def test_cad_rad_cell_archive_audit_reads_fusion_archive_evidence(self):
        audit = cad_rad_cell_archive_audit()
        self.assertEqual(audit["schema"], CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA)
        self.assertTrue(audit["summary"]["archivePresent"])
        self.assertTrue(audit["summary"]["zipReadable"])
        self.assertEqual(audit["summary"]["brepEntryCount"], 2)
        self.assertEqual(audit["summary"]["previewEntryCount"], 1)
        self.assertGreaterEqual(audit["summary"]["manifestEntryCount"], 3)
        self.assertTrue(audit["summary"]["canDeriveVisualReference"])
        self.assertFalse(audit["summary"]["canAttemptExactRigidBodyContact"])
        self.assertFalse(audit["summary"]["physicalAccuracyValidated"])
        self.assertIn("segmentedBodyMeshExport", audit["summary"]["missingEvidence"])
        self.assertIn("FusionAssetName[Active]/Previews/small.png", {entry["name"] for entry in audit["previewEntries"]})
        self.assertIn(
            "FusionAssetName[Active]/Breps.BlobParts/BREP.7d2e0b8b-fbfb-4c1a-a973-6c500ff05860.smb",
            {entry["name"] for entry in audit["brepEntries"]},
        )
        self.assertEqual(len(audit["archive"]["sha256"]), 64)
        self.assertTrue(any("RADs free cell" in label for label in audit["linkedLabels"]))

        exported = json.loads(export_cad_rad_cell_archive_audit_json())
        self.assertEqual(exported["schema"], CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA)

    def test_cad_mesh_audit_derives_hardware_profile_from_exported_obj(self):
        obj_text = "\n".join(
            [
                "o rad_cell_export",
                "v 0 0 0",
                "v 55.604331 0 0",
                "v 55.604331 55.604331 0",
                "v 0 55.604331 0",
                "v 0 0 20",
                "v 55.604331 0 20",
                "v 55.604331 55.604331 20",
                "v 0 55.604331 20",
                "f 1 2 3",
                "f 1 3 4",
                "f 5 7 6",
                "f 5 8 7",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "one_cell.obj"
            path.write_text(obj_text, encoding="utf-8")
            bounds = mesh_bounds_from_file(path)
            self.assertEqual(bounds.vertex_count, 8)
            self.assertEqual(bounds.face_count, 4)
            self.assertAlmostEqual(bounds.planar_span_mm, 55.604331)
            self.assertAlmostEqual(bounds.thickness_mm, 20.0)

            partial = cad_mesh_audit(path)
            self.assertEqual(partial["schema"], CAD_MESH_AUDIT_SCHEMA)
            self.assertEqual(partial["mesh"]["vertexCount"], 8)
            self.assertEqual(partial["mesh"]["faceCount"], 4)
            self.assertEqual(partial["calibrationReadiness"]["level"], "partial-measured")
            self.assertFalse(partial["claimBoundary"]["physicalAccuracyValidated"])
            self.assertIn("pin_radius_mm", partial["calibrationReadiness"]["meshMissingFields"])

            ready = json.loads(
                export_cad_mesh_audit_json(
                    path,
                    pin_radius_mm=1.36,
                    hole_radius_mm=1.7,
                    boss_radius_mm=3.1,
                    backlash_mm=0.34,
                )
            )
            self.assertEqual(ready["schema"], CAD_MESH_AUDIT_SCHEMA)
            self.assertEqual(ready["hardwareProfile"]["dimensionsMm"]["pinRadiusMm"], 1.36)
            self.assertEqual(ready["hardwareProfile"]["dimensionsMm"]["holeRadiusMm"], 1.7)
            self.assertEqual(ready["calibrationReadiness"]["level"], "mesh-calibrated")
            self.assertFalse(ready["calibrationReadiness"]["solverReady"])

    def test_cad_mesh_audit_derives_bounds_from_exported_step(self):
        step_text = "\n".join(
            [
                "ISO-10303-21;",
                "DATA;",
                "#1=CARTESIAN_POINT('',(0.,0.,0.));",
                "#2=CARTESIAN_POINT('',(55.604331,0.,0.));",
                "#3=CARTESIAN_POINT('',(55.604331,55.604331,0.));",
                "#4=CARTESIAN_POINT('',(0.,55.604331,0.));",
                "#5=CARTESIAN_POINT('',(0.,0.,20.));",
                "#6=CARTESIAN_POINT('',(55.604331,55.604331,20.));",
                "#10=ADVANCED_FACE('',(),#20,.T.);",
                "#11=ADVANCED_FACE('',(),#21,.T.);",
                "ENDSEC;",
                "END-ISO-10303-21;",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "one_cell.step"
            path.write_text(step_text, encoding="utf-8")
            bounds = mesh_bounds_from_file(path)
            self.assertEqual(bounds.vertex_count, 6)
            self.assertEqual(bounds.face_count, 2)
            self.assertAlmostEqual(bounds.planar_span_mm, 55.604331)
            self.assertAlmostEqual(bounds.thickness_mm, 20.0)

            audit = cad_mesh_audit(
                path,
                pin_radius_mm=1.36,
                hole_radius_mm=1.7,
                boss_radius_mm=3.1,
                backlash_mm=0.34,
            )
            self.assertEqual(audit["schema"], CAD_MESH_AUDIT_SCHEMA)
            self.assertIn(".step", audit["supportedFormats"])
            self.assertEqual(audit["mesh"]["format"], ".step")
            self.assertEqual(audit["mesh"]["vertexCount"], 6)
            self.assertEqual(audit["calibrationReadiness"]["level"], "mesh-calibrated")
            self.assertFalse(audit["claimBoundary"]["physicalAccuracyValidated"])

    def test_one_cell_cad_export_audit_detects_step_export_for_packet_readiness(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        with tempfile.TemporaryDirectory() as tmp:
            missing = one_cell_cad_export_audit(config, cad_dir=tmp)
            self.assertEqual(missing["schema"], ONE_CELL_CAD_EXPORT_AUDIT_SCHEMA)
            self.assertFalse(missing["summary"]["exportDetected"])
            self.assertIn("oneCellStepOrMeshExport", missing["summary"]["missingEvidence"])

            step_path = Path(tmp) / "RADs_unit_cell.step"
            step_path.write_text(
                "\n".join(
                    [
                        "ISO-10303-21;",
                        "DATA;",
                        "#1=CARTESIAN_POINT('',(0.,0.,0.));",
                        "#2=CARTESIAN_POINT('',(55.604331,0.,0.));",
                        "#3=CARTESIAN_POINT('',(55.604331,55.604331,20.));",
                        "#10=ADVANCED_FACE('',(),#20,.T.);",
                        "ENDSEC;",
                        "END-ISO-10303-21;",
                    ]
                ),
                encoding="utf-8",
            )
            detected = one_cell_cad_export_audit(config, cad_dir=tmp)
            self.assertTrue(detected["summary"]["exportDetected"])
            self.assertEqual(detected["summary"]["status"], "one-cell-cad-export-detected-needs-segmentation")
            self.assertEqual(detected["cadMeshAudit"]["schema"], CAD_MESH_AUDIT_SCHEMA)
            self.assertEqual(detected["cadMeshAudit"]["mesh"]["format"], ".step")
            self.assertAlmostEqual(
                detected["cadMeshAudit"]["hardwareProfile"]["dimensionsMm"]["holeRadiusMm"],
                1.7,
            )
            exported = json.loads(export_one_cell_cad_export_audit_json(config, cad_dir=tmp))
            self.assertEqual(exported["schema"], ONE_CELL_CAD_EXPORT_AUDIT_SCHEMA)
            self.assertTrue(exported["summary"]["exportDetected"])

    def test_two_cell_segmented_cad_readiness_lists_missing_exact_contact_assets(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        report = two_cell_segmented_cad_readiness_report(config)
        self.assertEqual(report["schema"], TWO_CELL_SEGMENTED_CAD_READINESS_SCHEMA)
        self.assertTrue(report["localAssets"]["fusionArchive"]["exists"])
        self.assertTrue(report["localAssets"]["viewerSummary"]["exists"])
        self.assertEqual(len(report["fusionFragments"]), 3)
        self.assertEqual(len(report["bodySegmentationTasks"]), 4)
        self.assertEqual(len(report["jointAndContactTasks"]), 4)
        self.assertEqual(len(report["physicsParameterTasks"]), 5)
        self.assertEqual(len(report["engineHandoffTasks"]), 3)
        self.assertFalse(report["summary"]["segmentedCadReady"])
        self.assertFalse(report["summary"]["canRunExactRigidBodyContact"])
        self.assertTrue(report["summary"]["canRunReducedProxy"])
        self.assertFalse(report["summary"]["oneCellCadExportDetected"])
        self.assertEqual(report["oneCellCadExportAudit"]["schema"], ONE_CELL_CAD_EXPORT_AUDIT_SCHEMA)
        self.assertEqual(report["summary"]["detectedSegmentedCadAssetCount"], 0)
        self.assertIn("Export separated STEP/STL/OBJ bodies", report["summary"]["nextAction"])

        exported = json.loads(export_two_cell_segmented_cad_readiness_json(config))
        rows = list(csv.DictReader(io.StringIO(export_two_cell_segmented_cad_readiness_csv(config))))
        self.assertEqual(exported["schema"], TWO_CELL_SEGMENTED_CAD_READINESS_SCHEMA)
        self.assertEqual(len(rows), report["summary"]["requiredTaskCount"])
        self.assertIn("screw_pin_body", {row["id"] for row in rows})
        self.assertIn("assets/cad/segmented/screw_pin_body.step", next(row for row in rows if row["id"] == "screw_pin_body")["expectedAssets"])

        with tempfile.TemporaryDirectory() as tmp:
            intake = Path(tmp)
            for relative in [
                "upper_free_cell_body.step",
                "lower_cell_body.step",
                "screw_pin_body.step",
                "radial_pad_hole_surfaces.json",
                "joint_axes.json",
                "two_cell_connector_pairs.json",
                "mass_inertia.json",
                "contact_parameters.json",
                "lock_crown_geometry.json",
                "actuator_force_displacement.csv",
                "bench_coordinate_truth.csv",
                "mujoco/rad_two_cell.xml",
                "gazebo/rad_two_cell.sdf",
                "isaac/rad_two_cell.usd",
            ]:
                path = intake / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("placeholder", encoding="utf-8")
            detected = two_cell_segmented_cad_readiness_report(config, intake_dir=str(intake))
            self.assertTrue(detected["summary"]["segmentedCadReady"])
            self.assertEqual(detected["summary"]["missingRequiredTaskCount"], 0)
            self.assertGreater(detected["summary"]["detectedSegmentedCadAssetCount"], 0)
            self.assertEqual(detected["summary"]["status"], "segmented-cad-assets-detected-needs-validation")

    def test_two_cell_segmented_cad_intake_templates_generate_fillable_files(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        templates = two_cell_segmented_cad_intake_templates(config)
        self.assertEqual(templates["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_TEMPLATES_SCHEMA)
        self.assertEqual(templates["templateCount"], 8)
        by_path = {item["path"]: item for item in templates["templates"]}
        self.assertIn("assets/cad/segmented/joint_axes.json", by_path)
        joint_axes = json.loads(by_path["assets/cad/segmented/joint_axes.json"]["content"])
        self.assertEqual(len(joint_axes["pinHoleAxes"]), 8)
        self.assertEqual(joint_axes["pinHoleAxes"][1]["site"], "ne")
        connector_pairs = json.loads(by_path["assets/cad/segmented/two_cell_connector_pairs.json"]["content"])
        self.assertEqual([pair["name"] for pair in connector_pairs["pairs"]], ["upper", "middle", "lower"])
        actuator_rows = list(
            csv.DictReader(io.StringIO(by_path["assets/cad/segmented/actuator_force_displacement.csv"]["content"]))
        )
        bench_rows = list(
            csv.DictReader(io.StringIO(by_path["assets/cad/segmented/bench_coordinate_truth.csv"]["content"]))
        )
        self.assertEqual(len(actuator_rows), templates["summary"]["actuatorRows"])
        self.assertEqual(len(bench_rows), templates["summary"]["benchCoordinateRows"])
        self.assertGreater(len(bench_rows), 0)

        exported = json.loads(export_two_cell_segmented_cad_intake_templates_json(config))
        summary_rows = list(csv.DictReader(io.StringIO(export_two_cell_segmented_cad_intake_templates_csv(config))))
        self.assertEqual(exported["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_TEMPLATES_SCHEMA)
        self.assertEqual(len(summary_rows), templates["templateCount"])
        self.assertIn("contact_parameters.json", {Path(row["path"]).name for row in summary_rows})

    def test_two_cell_segmented_cad_intake_validation_flags_missing_default_intake(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        with tempfile.TemporaryDirectory() as tmp:
            report = two_cell_segmented_cad_intake_validation_report(config, intake_dir=tmp)
        self.assertEqual(report["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_VALIDATION_SCHEMA)
        self.assertEqual(report["summary"]["existingTemplateCount"], 0)
        self.assertEqual(report["summary"]["completeTemplateCount"], 0)
        self.assertFalse(report["summary"]["intakeValidationReady"])
        self.assertFalse(report["summary"]["canAttemptExactRigidBodyContact"])
        self.assertIn("joint_axes.json", {Path(row["path"]).name for row in report["templateValidationRows"]})
        self.assertIn("upper_free_cell_body", report["summary"]["missingEvidence"])

        exported = json.loads(export_two_cell_segmented_cad_intake_validation_json(config, intake_dir=tmp))
        rows = list(csv.DictReader(io.StringIO(export_two_cell_segmented_cad_intake_validation_csv(config, intake_dir=tmp))))
        self.assertEqual(exported["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_VALIDATION_SCHEMA)
        self.assertEqual(len(rows), 15)
        self.assertIn("joint_axes.json", {row["id"] for row in rows})

    def test_two_cell_segmented_cad_intake_validation_accepts_filled_temp_intake(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        with tempfile.TemporaryDirectory() as tmp:
            intake = Path(tmp)
            templates = two_cell_segmented_cad_intake_templates(config, intake_dir=str(intake))
            by_name = {Path(item["path"]).name: item for item in templates["templates"]}

            for relative in ("upper_free_cell_body.step", "lower_cell_body.step", "screw_pin_body.step"):
                (intake / relative).write_text("solid segmented-body", encoding="utf-8")
            (intake / "mujoco").mkdir(parents=True, exist_ok=True)
            (intake / "mujoco" / "rad_two_cell.xml").write_text("<mujoco/>", encoding="utf-8")

            joint_axes = json.loads(by_name["joint_axes.json"]["content"])
            for axis in joint_axes["axes"]:
                if axis["id"] == "vertical_slide_axis":
                    axis["travelMinMm"] = -0.7
                    axis["travelMaxMm"] = 0.7
                    axis["source"] = "bench-filled"
            (intake / "joint_axes.json").write_text(json.dumps(joint_axes), encoding="utf-8")

            surfaces = json.loads(by_name["radial_pad_hole_surfaces.json"]["content"])
            for surface in surfaces["surfaces"]:
                surface["fusionBodyOrFaceId"] = f"face_{surface['site']}"
            (intake / "radial_pad_hole_surfaces.json").write_text(json.dumps(surfaces), encoding="utf-8")

            connectors = json.loads(by_name["two_cell_connector_pairs.json"]["content"])
            for pair in connectors["pairs"]:
                pair["measuredLeftOriginMm"] = [0.0, 0.0, 0.0]
                pair["measuredRightOriginMm"] = [55.6, 0.0, 0.0]
                pair["source"] = "bench-filled"
            (intake / "two_cell_connector_pairs.json").write_text(json.dumps(connectors), encoding="utf-8")

            mass_inertia = json.loads(by_name["mass_inertia.json"]["content"])
            for body in mass_inertia["bodies"]:
                body["massKg"] = 0.015
                body["centerOfMassMm"] = [0.0, 0.0, 3.0]
                body["inertiaTensorKgMm2"] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            (intake / "mass_inertia.json").write_text(json.dumps(mass_inertia), encoding="utf-8")

            contacts = json.loads(by_name["contact_parameters.json"]["content"])
            for contact in contacts["contacts"]:
                contact["normalStiffness"] = 1200
                contact["normalDamping"] = 4.5
                contact["staticFriction"] = 0.45
                contact["dynamicFriction"] = 0.35
            (intake / "contact_parameters.json").write_text(json.dumps(contacts), encoding="utf-8")

            lock_rows = list(csv.DictReader(io.StringIO(by_name["lock_crown_geometry.csv"]["content"])))
            for row in lock_rows:
                row["effectiveLatticeThetaDeg"] = "24.0"
                row["armWidthCorrectionDeg"] = "6.0"
                row["falloutObserved"] = "no"
            with (intake / "lock_crown_geometry.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(lock_rows[0]))
                writer.writeheader()
                writer.writerows(lock_rows)

            actuator_rows = list(csv.DictReader(io.StringIO(by_name["actuator_force_displacement.csv"]["content"])))
            for row in actuator_rows:
                row["measuredDisplacementMm"] = str(float(row["command"]) * 2.0)
                row["measuredForceN"] = "1.2"
            with (intake / "actuator_force_displacement.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(actuator_rows[0]))
                writer.writeheader()
                writer.writerows(actuator_rows)

            bench_rows = list(csv.DictReader(io.StringIO(by_name["bench_coordinate_truth.csv"]["content"])))
            for row in bench_rows:
                row["xMm"] = "0.0"
                row["yMm"] = "0.0"
                row["zMm"] = "0.0"
                row["measurementSource"] = "fixture-tracker"
            with (intake / "bench_coordinate_truth.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(bench_rows[0]))
                writer.writeheader()
                writer.writerows(bench_rows)

            report = two_cell_segmented_cad_intake_validation_report(config, intake_dir=str(intake))

        self.assertEqual(report["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_VALIDATION_SCHEMA)
        self.assertEqual(report["summary"]["completeTemplateCount"], report["summary"]["templateCount"])
        self.assertEqual(report["summary"]["completeBodyAssetCount"], report["summary"]["requiredBodyAssetCount"])
        self.assertEqual(report["summary"]["engineHandoffAssetCount"], 1)
        self.assertTrue(report["summary"]["intakeValidationReady"])
        self.assertTrue(report["summary"]["canAttemptExactRigidBodyContact"])

    def test_two_cell_bench_reports_cad_contact_proxy(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.08, pin_radius=0.18, hole_radius=0.225)
        controls = TwoCellBenchControls(alpha_command=-0.42, z_command=0.4)
        result = simulate_two_cell_bench(config, controls)
        self.assertEqual(result["schema"], TWO_CELL_BENCH_SCHEMA)
        self.assertEqual(result["model"], "cad-derived-two-cell-rigid-contact-proxy")
        self.assertEqual(result["cadReference"]["boundingBoxMm"]["heightZ"], 19.99999621152464)
        self.assertEqual(len(result["cells"]), 2)
        self.assertEqual(result["connector"]["contactMode"], "clearance-taken-up")
        self.assertGreater(abs(result["cells"][1]["residualAlpha"]), 0.0)
        self.assertGreater(result["cells"][1]["center"]["z"], 0.0)
        self.assertAlmostEqual(result["lockChecks"]["leftFixtureError"], 0.0)

    def test_two_cell_bench_lock_modes_suppress_motion(self):
        config = LatticeConfig(rows=1, cols=2)
        state_locked = simulate_two_cell_bench(
            config,
            TwoCellBenchControls(alpha_command=-0.5, z_command=0.5, right_locked=True),
        )
        self.assertEqual(state_locked["connector"]["contactMode"], "blocked-by-lock")
        self.assertTrue(state_locked["lockChecks"]["rightStateMotionSuppressed"])
        self.assertAlmostEqual(state_locked["cells"][1]["residualAlpha"], 0.0)
        self.assertAlmostEqual(state_locked["cells"][1]["center"]["z"], 0.0)
        position_locked = simulate_two_cell_bench(
            config,
            TwoCellBenchControls(alpha_command=-0.5, z_command=0.5, right_position_locked=True),
        )
        self.assertEqual(position_locked["connector"]["contactMode"], "blocked-by-lock")
        self.assertAlmostEqual(position_locked["cells"][1]["center"]["x"], config.cell_size)
        self.assertAlmostEqual(position_locked["lockChecks"]["rightFixtureError"], 0.0)

    def test_two_cell_clearance_sweep_tracks_backlash_and_dieoff(self):
        config = LatticeConfig(
            rows=1,
            cols=2,
            backlash=0.02,
            z_coupling_gain=0.35,
            pin_radius=0.18,
            hole_radius=0.2,
        )
        controls = TwoCellBenchControls(
            alpha_command=-0.45,
            z_command=0.45,
            left_position_locked=False,
            hole_sweep_max=0.5,
            hole_sweep_steps=7,
        )
        sweep = sweep_two_cell_backlash(config, controls)
        self.assertEqual(sweep["schema"], TWO_CELL_SWEEP_SCHEMA)
        self.assertEqual(len(sweep["rows"]), 7)
        self.assertTrue(sweep["trend"]["clearanceIncreasesBacklash"])
        self.assertTrue(sweep["trend"]["neighborResponseDropsWithClearance"])
        self.assertGreater(abs(sweep["trend"]["neighborZStart"]), abs(sweep["trend"]["neighborZEnd"]))

    def test_two_cell_actuation_sweep_covers_lock_modes_and_hole_radii(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.03, pin_radius=0.18, hole_radius=0.225)
        sweep = two_cell_actuation_sweep(
            config,
            alpha_commands=[-0.4, 0.4],
            z_commands=[-0.3, 0.3],
            hole_radii=[0.18, 0.225, 0.5],
        )
        self.assertEqual(sweep["schema"], TWO_CELL_ACTUATION_SWEEP_SCHEMA)
        self.assertEqual(sweep["summary"]["rowCount"], 48)
        self.assertEqual(len(sweep["rows"]), 48)
        self.assertEqual(
            set(sweep["summary"]["lockModes"]),
            {"free", "right_state_locked", "right_position_locked", "left_free"},
        )
        locked = next(row for row in sweep["rows"] if row["lockMode"] == "right_state_locked")
        self.assertTrue(locked["rightStateMotionSuppressed"])
        self.assertAlmostEqual(locked["rightZ"], 0.0)
        position_locked = next(row for row in sweep["rows"] if row["lockMode"] == "right_position_locked")
        self.assertAlmostEqual(position_locked["rightFixtureError"], 0.0)
        clearances = sorted({row["clearance"] for row in sweep["rows"]})
        self.assertEqual(clearances, sorted(clearances))
        self.assertGreater(clearances[-1], clearances[0])
        sweep_json = json.loads(
            export_two_cell_actuation_sweep_json(
                config,
                alpha_commands=[-0.4, 0.4],
                z_commands=[-0.3, 0.3],
                hole_radii=[0.18, 0.225, 0.5],
            )
        )
        sweep_rows = list(
            csv.DictReader(
                io.StringIO(
                    export_two_cell_actuation_sweep_csv(
                        config,
                        alpha_commands=[-0.4, 0.4],
                        z_commands=[-0.3, 0.3],
                        hole_radii=[0.18, 0.225, 0.5],
                    )
                )
            )
        )
        self.assertEqual(sweep_json["schema"], TWO_CELL_ACTUATION_SWEEP_SCHEMA)
        self.assertEqual(len(sweep_rows), 48)

    def test_two_cell_connector_contact_reports_three_cad_pin_hole_pairs(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        controls = TwoCellBenchControls(alpha_command=-0.35, z_command=0.35)
        report = two_cell_connector_contact_report(config, controls, gravity_force=0.0)
        self.assertEqual(report["schema"], TWO_CELL_CONNECTOR_CONTACT_SCHEMA)
        self.assertTrue(report["summary"]["connectorContactReady"])
        self.assertEqual(report["summary"]["connectorCount"], 3)
        self.assertEqual(
            [connector["connector"] for connector in report["connectors"]],
            ["upper", "middle", "lower"],
        )
        self.assertEqual(report["connectors"][1]["leftSite"], "e")
        self.assertEqual(report["connectors"][1]["rightSite"], "w")
        self.assertAlmostEqual(report["dimensions"]["holeRadiusMm"], 1.7)
        self.assertGreater(report["summary"]["maxVerticalSlipMm"], 0.0)

        state_locked = two_cell_connector_contact_report(
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35, right_locked=True),
            gravity_force=0.0,
        )
        self.assertTrue(state_locked["summary"]["stateLocksHold"])
        self.assertTrue(state_locked["summary"]["rightLockedVerticalMotionAllowed"])
        position_locked = two_cell_connector_contact_report(
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35, right_position_locked=True),
            gravity_force=0.0,
        )
        self.assertTrue(position_locked["summary"]["positionLocksHold"])
        self.assertAlmostEqual(position_locked["quasistatic"]["cells"][1]["center"]["z"], 0.0)

        tight = two_cell_connector_contact_report(config, controls, hole_radius=0.18, gravity_force=0.0)
        loose = two_cell_connector_contact_report(config, controls, hole_radius=0.5, gravity_force=0.0)
        self.assertGreater(loose["dimensions"]["pinHoleClearanceMm"], tight["dimensions"]["pinHoleClearanceMm"])
        self.assertLessEqual(
            loose["summary"]["maxVerticalExcessMm"],
            tight["summary"]["maxVerticalExcessMm"],
        )

        sweep = sweep_two_cell_connector_contact(
            config,
            alpha_commands=[-0.4],
            z_commands=[0.3],
            hole_radii=[0.18, 0.225, 0.5],
            lock_modes=["free", "right_position_locked"],
            gravity_force=0.0,
        )
        self.assertEqual(sweep["schema"], TWO_CELL_CONNECTOR_CONTACT_SWEEP_SCHEMA)
        self.assertEqual(sweep["summary"]["rowCount"], 6)
        self.assertEqual(sweep["summary"]["connectorRowCount"], 18)
        self.assertEqual(sweep["summary"]["solverFailureCount"], 0)

        report_json = json.loads(export_two_cell_connector_contact_json(config, controls, gravity_force=0.0))
        report_rows = list(csv.DictReader(io.StringIO(export_two_cell_connector_contact_csv(config, controls, gravity_force=0.0))))
        sweep_json = json.loads(
            export_two_cell_connector_contact_sweep_json(
                config,
                alpha_commands=[-0.4],
                z_commands=[0.3],
                hole_radii=[0.18, 0.225],
                gravity_force=0.0,
            )
        )
        sweep_rows = list(
            csv.DictReader(
                io.StringIO(
                    export_two_cell_connector_contact_sweep_csv(
                        config,
                        alpha_commands=[-0.4],
                        z_commands=[0.3],
                        hole_radii=[0.18, 0.225],
                        gravity_force=0.0,
                    )
                )
            )
        )
        self.assertEqual(report_json["schema"], TWO_CELL_CONNECTOR_CONTACT_SCHEMA)
        self.assertEqual(len(report_rows), 3)
        self.assertEqual(sweep_json["schema"], TWO_CELL_CONNECTOR_CONTACT_SWEEP_SCHEMA)
        self.assertEqual(len(sweep_rows), sweep_json["summary"]["rowCount"])

    def test_two_cell_mjcf_proxy_exports_cad_dimensioned_locks_and_clearance(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        free = two_cell_mjcf_proxy_report(
            config,
            TwoCellBenchControls(alpha_command=-0.4, z_command=0.35),
        )
        self.assertEqual(free["schema"], TWO_CELL_MJCF_SCHEMA)
        self.assertEqual(free["summary"]["bodyCount"], 2)
        self.assertEqual(free["summary"]["fixedBodyCount"], 1)
        self.assertEqual(free["summary"]["connectorContactPairCount"], 3)
        self.assertGreater(free["summary"]["actuatorCount"], 0)
        self.assertIn("<mujoco", free["xml"])
        self.assertIn("left_cell_e_pin", free["xml"])
        self.assertIn("left_cell_e_marker", free["xml"])
        self.assertIn("right_cell_w_hole_clearance", free["xml"])
        self.assertIn('joint="right_cell_slide_z"', free["xml"])
        self.assertEqual(ET.fromstring(free["xml"]).tag, "mujoco")
        self.assertAlmostEqual(free["scale"]["holeRadiusMm"], 1.7)
        self.assertGreater(free["scale"]["clearanceMm"], 0.0)

        wide_hole = two_cell_mjcf_proxy_report(
            config,
            TwoCellBenchControls(alpha_command=-0.4, z_command=0.35),
            hole_radius=0.5,
        )
        self.assertGreater(wide_hole["summary"]["clearanceMm"], free["summary"]["clearanceMm"])

        locked = two_cell_mjcf_proxy_report(
            config,
            TwoCellBenchControls(alpha_command=-0.4, z_command=0.35, right_locked=True),
        )
        self.assertEqual(locked["summary"]["fixedBodyCount"], 2)
        self.assertEqual(locked["summary"]["actuatorCount"], 0)
        self.assertTrue(locked["summary"]["rightStateMotionSuppressed"])
        self.assertTrue(
            next(record for record in locked["locks"] if record["body"] == "right_cell")["fixedInMjcf"]
        )

        position_locked = two_cell_mjcf_proxy_report(
            config,
            TwoCellBenchControls(alpha_command=-0.4, z_command=0.35, right_position_locked=True),
        )
        self.assertEqual(position_locked["summary"]["fixedBodyCount"], 2)
        self.assertEqual(position_locked["summary"]["actuatorCount"], 0)
        self.assertAlmostEqual(position_locked["bench"]["lockChecks"]["rightFixtureError"], 0.0)

        exported_json = json.loads(export_two_cell_mjcf_proxy_json(config))
        exported_xml = export_two_cell_mjcf_xml(config)
        self.assertEqual(exported_json["schema"], TWO_CELL_MJCF_SCHEMA)
        self.assertIn("rad_two_cell_cad_proxy", exported_xml)

    def test_two_cell_mujoco_proxy_run_reports_optional_engine_status(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        unavailable = two_cell_mujoco_proxy_run_report(
            config,
            TwoCellBenchControls(alpha_command=-0.4, z_command=0.35),
            engine_availability={"mujoco": False},
            steps=4,
        )
        self.assertEqual(unavailable["schema"], TWO_CELL_MUJOCO_RUN_SCHEMA)
        self.assertFalse(unavailable["summary"]["twoCellMujocoProxyRunComplete"])
        self.assertFalse(unavailable["summary"]["mujocoAvailable"])
        self.assertIn("mujocoPythonPackage", unavailable["summary"]["missingEvidence"])
        self.assertEqual(unavailable["summary"]["expectedBodyResultCount"], 2)
        self.assertEqual(unavailable["mjcfProxy"]["schema"], TWO_CELL_MJCF_SCHEMA)

        current = json.loads(export_two_cell_mujoco_proxy_run_json(config, steps=4))
        self.assertEqual(current["schema"], TWO_CELL_MUJOCO_RUN_SCHEMA)
        self.assertEqual(current["summary"]["expectedBodyResultCount"], 2)
        if current["summary"]["mujocoAvailable"]:
            self.assertEqual(current["solver"]["stepsRequested"], 4)
            self.assertIn("bodies", current["results"])
        else:
            self.assertFalse(current["summary"]["twoCellMujocoProxyRunComplete"])
            self.assertIn("mujocoPythonPackage", current["summary"]["missingEvidence"])

    def test_two_cell_quasistatic_solver_distinguishes_state_and_position_locks(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        state_locked = solve_two_cell_quasistatic(
            config,
            TwoCellBenchControls(alpha_command=-0.5, z_command=0.45, right_locked=True),
            gravity_force=0.0,
        )
        self.assertEqual(state_locked["schema"], TWO_CELL_QUASISTATIC_SCHEMA)
        self.assertTrue(state_locked["solver"]["success"])
        self.assertTrue(state_locked["constraints"]["stateLocksHold"])
        self.assertAlmostEqual(state_locked["cells"][1]["alpha"], config.initial_alpha)
        self.assertGreater(state_locked["cells"][1]["center"]["z"], 0.0)
        self.assertTrue(state_locked["constraints"]["rightLockedVerticalMotionAllowed"])

        position_locked = solve_two_cell_quasistatic(
            config,
            TwoCellBenchControls(alpha_command=-0.5, z_command=0.45, right_position_locked=True),
            gravity_force=0.0,
        )
        self.assertTrue(position_locked["constraints"]["positionLocksHold"])
        self.assertAlmostEqual(position_locked["cells"][1]["center"]["x"], config.cell_size)
        self.assertAlmostEqual(position_locked["cells"][1]["center"]["z"], 0.0)
        self.assertAlmostEqual(position_locked["constraints"]["rightFixtureError"], 0.0)

    def test_two_cell_quasistatic_gravity_sag_increases_with_hole_clearance(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        controls = TwoCellBenchControls(alpha_command=0.0, z_command=0.0)
        tight = solve_two_cell_quasistatic(
            config,
            controls,
            hole_radius=0.18,
            gravity_force=0.5,
        )
        loose = solve_two_cell_quasistatic(
            config,
            controls,
            hole_radius=0.5,
            gravity_force=0.5,
        )
        self.assertTrue(tight["solver"]["success"])
        self.assertTrue(loose["solver"]["success"])
        self.assertLess(loose["cells"][1]["center"]["z"], tight["cells"][1]["center"]["z"])
        self.assertLessEqual(
            abs(tight["contact"]["verticalShear"]),
            tight["contact"]["verticalClearance"] + 0.01,
        )
        self.assertGreater(loose["contact"]["verticalClearance"], tight["contact"]["verticalClearance"])

    def test_two_cell_quasistatic_sweep_exports_solver_rows(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        sweep = sweep_two_cell_quasistatic(
            config,
            alpha_commands=[-0.3, 0.3],
            z_commands=[-0.2, 0.2],
            hole_radii=[0.18, 0.225],
            lock_modes=["free", "right_state_locked", "right_position_locked"],
            gravity_force=0.0,
        )
        self.assertEqual(sweep["schema"], TWO_CELL_QUASISTATIC_SWEEP_SCHEMA)
        self.assertEqual(sweep["summary"]["rowCount"], 24)
        self.assertEqual(sweep["summary"]["solverFailureCount"], 0)
        self.assertTrue(all(row["solverSuccess"] for row in sweep["rows"]))
        locked = next(row for row in sweep["rows"] if row["lockMode"] == "right_state_locked")
        self.assertAlmostEqual(locked["rightStateAlphaError"], 0.0)
        position_locked = next(row for row in sweep["rows"] if row["lockMode"] == "right_position_locked")
        self.assertAlmostEqual(position_locked["rightFixtureError"], 0.0)
        sweep_json = json.loads(
            export_two_cell_quasistatic_sweep_json(
                config,
                alpha_commands=[-0.3, 0.3],
                z_commands=[-0.2, 0.2],
                hole_radii=[0.18, 0.225],
                lock_modes=["free", "right_state_locked", "right_position_locked"],
                gravity_force=0.0,
            )
        )
        sweep_rows = list(
            csv.DictReader(
                io.StringIO(
                    export_two_cell_quasistatic_sweep_csv(
                        config,
                        alpha_commands=[-0.3, 0.3],
                        z_commands=[-0.2, 0.2],
                        hole_radii=[0.18, 0.225],
                        lock_modes=["free", "right_state_locked", "right_position_locked"],
                        gravity_force=0.0,
                    )
                )
            )
        )
        solved = json.loads(export_two_cell_quasistatic_json(config, gravity_force=0.0))
        self.assertEqual(sweep_json["schema"], TWO_CELL_QUASISTATIC_SWEEP_SCHEMA)
        self.assertEqual(len(sweep_rows), 24)
        self.assertEqual(solved["schema"], TWO_CELL_QUASISTATIC_SCHEMA)

    def test_two_cell_physical_simulation_suite_covers_locks_clearance_and_backlash(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        suite = two_cell_physical_simulation_suite(
            config,
            TwoCellBenchControls(alpha_command=-0.4, z_command=0.35),
            gravity_force=0.0,
        )
        self.assertEqual(suite["schema"], TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA)
        self.assertEqual(suite["summary"]["caseCount"], 12)
        self.assertEqual(suite["summary"]["solverFailureCount"], 0)
        self.assertTrue(suite["summary"]["internalSuiteReady"])
        self.assertFalse(suite["summary"]["physicalAccuracyValidated"])
        self.assertTrue(suite["summary"]["clearanceTrend"]["looseHasLargerClearance"])
        self.assertTrue(suite["summary"]["backlashTrend"]["highBacklashReducesAlphaResponse"])
        self.assertIn("externalMuJoCoRun", suite["summary"]["missingEvidence"])

        by_id = {row["caseId"]: row for row in suite["rows"]}
        self.assertGreater(by_id["loose_clearance_contract_lift"]["clearance"], by_id["tight_clearance_contract_lift"]["clearance"])
        self.assertAlmostEqual(by_id["right_state_locked_lift"]["rightAlpha"], config.initial_alpha)
        self.assertTrue(by_id["right_state_locked_lift"]["rightLockedVerticalMotionAllowed"])
        self.assertAlmostEqual(by_id["right_position_locked_lift"]["rightX"], config.cell_size)
        self.assertAlmostEqual(by_id["right_position_locked_lift"]["rightZ"], 0.0)
        self.assertLessEqual(
            abs(by_id["high_backlash_left_free_lift"]["leftAlphaDeltaFromInitial"]),
            abs(by_id["zero_backlash_left_free_lift"]["leftAlphaDeltaFromInitial"]) + 1e-7,
        )

        suite_json = json.loads(export_two_cell_physical_simulation_suite_json(config, gravity_force=0.0))
        suite_rows = list(
            csv.DictReader(io.StringIO(export_two_cell_physical_simulation_suite_csv(config, gravity_force=0.0)))
        )
        self.assertEqual(suite_json["schema"], TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA)
        self.assertEqual(len(suite_rows), suite["summary"]["caseCount"])
        self.assertIn("connectorMaxVerticalSlipMm", suite_rows[0])
        self.assertIn("leftAlphaDeltaFromInitial", suite_rows[0])

    def test_two_cell_physical_fidelity_matrix_sweeps_radius_backlash_locks_and_actuation(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        matrix = two_cell_physical_fidelity_matrix(
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35),
            tolerance=1e-7,
        )
        self.assertEqual(matrix["schema"], TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA)
        self.assertEqual(matrix["summary"]["rowCount"], 252)
        self.assertEqual(matrix["summary"]["connectorRowCount"], 756)
        self.assertEqual(set(matrix["axes"]["lockModes"]), {"free", "right_state_locked", "right_position_locked", "left_free"})
        self.assertIn("contract_lift", matrix["axes"]["actuationCases"])
        self.assertTrue(matrix["summary"]["internalMatrixReady"])
        self.assertEqual(matrix["summary"]["status"], "reduced-physics-fidelity-matrix-ready-needs-external-data")
        self.assertFalse(matrix["summary"]["physicalAccuracyValidated"])
        self.assertTrue(matrix["summary"]["clearanceTrend"]["clearanceNondecreasing"])
        self.assertFalse(matrix["summary"]["clearanceTrend"]["rawVerticalSlipNonincreasing"])
        self.assertFalse(matrix["summary"]["clearanceTrend"]["verticalExcessNonincreasing"])
        self.assertTrue(matrix["summary"]["clearanceTrend"]["penaltyNonincreasing"])
        self.assertEqual(len(matrix["summary"]["clearanceTrend"]["verticalFreePlayUtilization"]), 3)
        self.assertTrue(matrix["summary"]["backlashTrend"]["alphaResponseNonincreasing"])
        self.assertFalse(matrix["summary"]["backlashTrend"]["verticalResponseNonincreasing"])
        self.assertTrue(matrix["summary"]["backlashTrend"]["verticalResponseMateriallyFlat"])
        self.assertTrue(matrix["summary"]["lockChecks"]["stateLockAlphaHold"])
        self.assertGreater(matrix["summary"]["lockChecks"]["stateLockVerticalMotionSamples"], 0)
        self.assertTrue(matrix["summary"]["lockChecks"]["positionLockFixtureHold"])
        self.assertTrue(all(matrix["summary"]["actuationPolarity"].values()))
        self.assertIn("filledConnectorMeasurements", matrix["summary"]["missingEvidence"])

        exported = json.loads(export_two_cell_physical_fidelity_matrix_json(config))
        exported_rows = list(csv.DictReader(io.StringIO(export_two_cell_physical_fidelity_matrix_csv(config))))
        self.assertEqual(exported["schema"], TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA)
        self.assertEqual(len(exported_rows), matrix["summary"]["rowCount"])
        self.assertIn("connectorMaxVerticalSlipMm", exported_rows[0])

    def test_two_cell_contact_phase_map_classifies_matrix_rows_for_bench_planning(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        phase_map = two_cell_contact_phase_map(
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35),
            tolerance=1e-7,
        )
        self.assertEqual(phase_map["schema"], TWO_CELL_CONTACT_PHASE_MAP_SCHEMA)
        self.assertEqual(phase_map["summary"]["rowCount"], 252)
        self.assertEqual(len(phase_map["rows"]), 252)
        self.assertFalse(phase_map["summary"]["physicalAccuracyValidated"])
        self.assertEqual(
            phase_map["summary"]["status"],
            "contact-phase-map-ready-needs-bench-or-external-validation",
        )
        self.assertGreater(phase_map["summary"]["activeContactCaseCount"], 0)
        self.assertEqual(sum(phase_map["summary"]["phaseCounts"].values()), phase_map["summary"]["rowCount"])
        self.assertGreaterEqual(phase_map["summary"]["freePlayCaseCount"], 0)
        self.assertGreater(phase_map["summary"]["lockCaseCount"], 0)
        self.assertIn("axial-contact", phase_map["summary"]["phaseCounts"])
        self.assertIn("position-locked", phase_map["summary"]["phaseCounts"])
        self.assertIn("state-locked-vertical-free", phase_map["summary"]["phaseCounts"])
        self.assertTrue(phase_map["summary"]["measurementPriority"])
        self.assertIn("exact real-cell phase boundaries", phase_map["claimBoundary"]["blockedClaim"])
        self.assertTrue(all("caseId" in row and "phase" in row for row in phase_map["rows"]))

        exported = json.loads(export_two_cell_contact_phase_map_json(config))
        exported_rows = list(csv.DictReader(io.StringIO(export_two_cell_contact_phase_map_csv(config))))
        self.assertEqual(exported["schema"], TWO_CELL_CONTACT_PHASE_MAP_SCHEMA)
        self.assertEqual(len(exported_rows), phase_map["summary"]["rowCount"])
        self.assertIn("phase", exported_rows[0])
        self.assertIn("contactMode", exported_rows[0])

    def test_two_cell_radius_backlash_phase_diagram_sweeps_dense_hole_backlash_grid(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        diagram = two_cell_radius_backlash_phase_diagram(
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35, hole_sweep_max=0.5),
            hole_steps=5,
            backlash_steps=6,
            actuation_case="contract_lift",
            lock_mode="free",
        )
        self.assertEqual(diagram["schema"], TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA)
        self.assertEqual(diagram["summary"]["holeStepCount"], 5)
        self.assertEqual(diagram["summary"]["backlashStepCount"], 6)
        self.assertEqual(diagram["summary"]["rowCount"], 30)
        self.assertEqual(len(diagram["phaseGrid"]), 6)
        self.assertEqual(len(diagram["phaseGrid"][0]), 5)
        self.assertEqual(sum(diagram["summary"]["phaseCounts"].values()), 30)
        self.assertGreater(diagram["summary"]["activeContactCaseCount"], 0)
        self.assertFalse(diagram["summary"]["physicalAccuracyValidated"])
        self.assertIn("real transition boundaries", diagram["claimBoundary"]["blockedClaim"])
        self.assertTrue(all(row["actuationCase"] == "contract_lift" for row in diagram["rows"]))
        self.assertTrue(all(row["lockMode"] == "free" for row in diagram["rows"]))

        exported = json.loads(export_two_cell_radius_backlash_phase_diagram_json(config, hole_steps=4, backlash_steps=4))
        exported_rows = list(csv.DictReader(io.StringIO(export_two_cell_radius_backlash_phase_diagram_csv(config, hole_steps=4, backlash_steps=4))))
        self.assertEqual(exported["schema"], TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA)
        self.assertEqual(len(exported_rows), 16)
        self.assertIn("phase", exported_rows[0])
        self.assertIn("holeRadius", exported_rows[0])

    def test_two_cell_radius_backlash_transition_report_prioritizes_boundary_measurements(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        report = two_cell_radius_backlash_transition_report(
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35, hole_sweep_max=0.5),
            hole_steps=5,
            backlash_steps=6,
            actuation_case="contract_lift",
            lock_mode="free",
        )
        self.assertEqual(report["schema"], TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA)
        self.assertFalse(report["summary"]["physicalAccuracyValidated"])
        self.assertGreater(report["summary"]["transitionBracketCount"], 0)
        self.assertGreater(report["summary"]["measurementCaseCount"], 0)
        self.assertIn("observed contact/free-play phase", report["measurementProtocol"]["requiredObservables"])
        self.assertTrue(all(item["physicalAccuracyValidated"] is False for item in report["transitionBrackets"]))
        self.assertTrue(all("observedPhase" in item and "phaseDiagramCaseId" in item for item in report["measurementCases"]))
        self.assertIn("real RAD phase-transition law", report["claimBoundary"]["blockedClaim"])

        exported = json.loads(
            export_two_cell_radius_backlash_transition_report_json(config, hole_steps=4, backlash_steps=4)
        )
        exported_rows = list(
            csv.DictReader(
                io.StringIO(export_two_cell_radius_backlash_transition_report_csv(config, hole_steps=4, backlash_steps=4))
            )
        )
        self.assertEqual(exported["schema"], TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA)
        self.assertGreaterEqual(len(exported_rows), 1)
        self.assertIn("observedConnectorMaxVerticalSlipMm", exported_rows[0])
        self.assertIn("measurementTarget", exported_rows[0])

    def test_two_cell_radius_backlash_transition_measurements_compare_phase_boundary_shift(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        report = two_cell_radius_backlash_transition_report(
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35, hole_sweep_max=0.5),
            hole_steps=5,
            backlash_steps=6,
            actuation_case="contract_lift",
            lock_mode="free",
        )
        exact = dict(report["measurementCases"][0])
        exact.update(
            {
                "observedPhase": exact["phase"],
                "observedRightZ": exact["rightZ"],
                "observedRightAlpha": exact["rightAlpha"],
                "observedConnectorMaxVerticalSlipMm": exact["connectorMaxVerticalSlipMm"],
                "observedConnectorMaxVerticalExcessMm": exact["connectorMaxVerticalExcessMm"],
            }
        )
        shifted = dict(report["measurementCases"][0])
        shifted["observedPhase"] = "axial-contact"
        comparison = compare_two_cell_radius_backlash_transition_measurements(
            [exact, shifted],
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35, hole_sweep_max=0.5),
            hole_steps=5,
            backlash_steps=6,
            tolerance=1e-9,
        )
        self.assertEqual(comparison["schema"], TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA)
        self.assertEqual(comparison["summary"]["providedRowCount"], 2)
        self.assertEqual(comparison["summary"]["matchedRowCount"], 2)
        self.assertEqual(comparison["summary"]["phaseObservationCount"], 2)
        self.assertLess(comparison["summary"]["phaseAccuracy"], 1.0)
        self.assertEqual(
            comparison["summary"]["transitionShiftDirection"],
            "effective-hole-transition-lower-than-reduced-model",
        )
        self.assertFalse(comparison["acceptance"]["physicalAccuracyValidated"])
        self.assertIn("rightZScale", comparison["calibrationEstimate"])
        self.assertGreater(
            comparison["calibrationEstimate"]["proposedReducedProxyUpdates"]["holeRadius"],
            config.hole_radius,
        )
        self.assertEqual(
            comparison["calibrationEstimate"]["proposedReducedProxyUpdates"]["transitionBoundaryDirection"],
            "effective-hole-transition-lower-than-reduced-model",
        )
        rerun = two_cell_radius_backlash_transition_rerun_report(
            [exact, shifted],
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35, hole_sweep_max=0.5),
            hole_steps=5,
            backlash_steps=6,
            tolerance=1e-9,
        )
        self.assertEqual(rerun["schema"], TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA)
        self.assertGreater(rerun["axisBiases"]["effectiveHoleRadiusBias"], 0)
        self.assertEqual(
            rerun["summary"]["transitionShiftDirection"],
            "effective-hole-transition-lower-than-reduced-model",
        )
        self.assertGreaterEqual(rerun["summary"]["transitionMovementCount"], 1)
        self.assertFalse(rerun["claimBoundary"]["physicalAccuracyValidated"])
        self.assertGreater(rerun["calibrated"]["diagram"]["axes"]["holeRadiusBias"], 0)

        csv_text = export_two_cell_radius_backlash_transition_report_csv(
            config,
            TwoCellBenchControls(alpha_command=-0.35, z_command=0.35, hole_sweep_max=0.5),
            hole_steps=5,
            backlash_steps=6,
        )
        rows = list(csv.DictReader(io.StringIO(csv_text)))
        rows[0]["observedPhase"] = rows[0]["phase"]
        rows[0]["observedRightZ"] = rows[0]["rightZ"]
        parsed = two_cell_radius_backlash_transition_measurements_from_csv(
            "\n".join(
                [
                    ",".join(rows[0].keys()),
                    ",".join(str(value) for value in rows[0].values()),
                ]
            )
        )
        exported = json.loads(
            export_two_cell_radius_backlash_transition_comparison_json(
                parsed,
                config,
                hole_steps=5,
                backlash_steps=6,
            )
        )
        exported_rows = list(
            csv.DictReader(
                io.StringIO(
                    export_two_cell_radius_backlash_transition_comparison_csv(
                        parsed,
                        config,
                        hole_steps=5,
                        backlash_steps=6,
                    )
                )
            )
        )
        self.assertEqual(exported["schema"], TWO_CELL_RADIUS_BACKLASH_TRANSITION_COMPARISON_SCHEMA)
        self.assertGreaterEqual(len(exported_rows), 1)
        self.assertIn("transitionShiftVote", exported_rows[0])
        rerun_exported = json.loads(
            export_two_cell_radius_backlash_transition_rerun_json(
                [exact, shifted],
                config,
                hole_steps=5,
                backlash_steps=6,
            )
        )
        rerun_exported_rows = list(
            csv.DictReader(
                io.StringIO(
                    export_two_cell_radius_backlash_transition_rerun_csv(
                        [exact, shifted],
                        config,
                        hole_steps=5,
                        backlash_steps=6,
                    )
                )
            )
        )
        self.assertEqual(rerun_exported["schema"], TWO_CELL_RADIUS_BACKLASH_TRANSITION_RERUN_SCHEMA)
        self.assertIn("midpointDelta", rerun_exported_rows[0])
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "transition_measurements.csv"
            out = Path(tmp) / "out"
            source.write_text(
                "\n".join(
                    [
                        ",".join(rows[0].keys()),
                        ",".join(str(value) for value in rows[0].values()),
                    ]
                ),
                encoding="utf-8",
            )
            files = write_two_cell_radius_backlash_transition_comparison_artifacts(
                source,
                out,
                hole_steps=5,
                backlash_steps=6,
            )
            self.assertTrue(files["comparison"].exists())
            self.assertTrue(files["comparison_csv"].exists())
            self.assertTrue(files["rerun"].exists())
            self.assertTrue(files["rerun_csv"].exists())

    def test_two_cell_physical_response_atlas_prioritizes_lock_clearance_and_slip_cases(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        atlas = two_cell_physical_response_atlas(config, rank_limit=6)
        self.assertEqual(atlas["schema"], TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA)
        self.assertEqual(atlas["matrixSummary"]["rowCount"], 252)
        self.assertEqual(atlas["connectorSummary"]["connectorRowCount"], 756)
        self.assertFalse(atlas["claimBoundary"]["physicalAccuracyValidated"])
        self.assertEqual(atlas["claimBoundary"]["allowedClaim"], "internal reduced two-cell response atlas for selecting physical tests")
        self.assertTrue(atlas["invariants"]["solverAllConverged"])
        self.assertTrue(atlas["invariants"]["lockAllPass"])
        self.assertTrue(atlas["invariants"]["stateLockAlphaHolds"])
        self.assertGreater(atlas["invariants"]["stateLockAllowsVerticalSamples"], 0)
        self.assertTrue(atlas["invariants"]["positionLockFixtureHolds"])
        self.assertTrue(atlas["invariants"]["clearanceNondecreasing"])
        self.assertTrue(atlas["invariants"]["clearancePenaltyNonincreasing"])
        self.assertTrue(atlas["invariants"]["zLiftPositive"])
        self.assertTrue(atlas["invariants"]["zPushdownNegative"])
        self.assertEqual(len(atlas["groupSummaries"]["byLockMode"]), 4)
        self.assertEqual(len(atlas["groupSummaries"]["byHoleRadius"]), 3)
        self.assertGreater(atlas["connectorSummary"]["maxVerticalSlipMm"], 0.0)
        self.assertTrue(atlas["rankedCases"]["maxVerticalSlip"])
        self.assertTrue(atlas["benchPriority"]["firstCasesToMeasure"])
        self.assertTrue(atlas["benchPriority"]["lockCasesToMeasure"])
        self.assertIn("segmentedCADBodies", atlas["claimBoundary"]["remainingEvidence"])

        exported = json.loads(export_two_cell_physical_response_atlas_json(config))
        rows = list(csv.DictReader(io.StringIO(export_two_cell_physical_response_atlas_csv(config))))
        self.assertEqual(exported["schema"], TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA)
        self.assertGreaterEqual(len(rows), len(atlas["benchPriority"]["firstCasesToMeasure"]))
        self.assertIn("connectorMaxVerticalSlipMm", rows[0])
        self.assertIn(rows[0]["category"], {"firstCasesToMeasure", "lockCasesToMeasure", "clearanceCasesToMeasure"})

    def test_two_cell_fidelity_matrix_measurement_template_covers_all_matrix_connectors(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        template = two_cell_fidelity_matrix_measurement_template(config)
        rows = two_cell_fidelity_matrix_measurement_template_rows(config)
        self.assertEqual(template["schema"], TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA)
        self.assertEqual(template["summary"]["matrixRowCount"], 252)
        self.assertEqual(template["summary"]["connectorRowCount"], 756)
        self.assertEqual(len(rows), 756)
        self.assertEqual({row["connector"] for row in rows}, {"upper", "middle", "lower"})
        self.assertEqual(rows[0]["caseId"], "FM_000")
        self.assertIn("predictedRightCellZ", rows[0])
        self.assertIn("predictedVerticalSlipMm", rows[0])
        self.assertEqual(rows[0]["observedRightCellZ"], "")
        self.assertEqual(rows[0]["observedVerticalSlipMm"], "")

        exported = json.loads(export_two_cell_fidelity_matrix_measurement_template_json(config))
        exported_rows = list(
            csv.DictReader(io.StringIO(export_two_cell_fidelity_matrix_measurement_template_csv(config)))
        )
        self.assertEqual(exported["schema"], TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA)
        self.assertEqual(len(exported_rows), 756)
        self.assertIn("observedRightCellZ", exported_rows[0])

    def test_two_cell_external_fidelity_matrix_manifest_covers_engine_handoff(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        manifest = two_cell_external_fidelity_matrix_manifest(config, engine="mujoco")
        rows = two_cell_external_fidelity_matrix_manifest_rows(config)
        csv_rows = list(csv.DictReader(io.StringIO(export_two_cell_external_fidelity_matrix_manifest_csv(config))))
        exported = json.loads(export_two_cell_external_fidelity_matrix_manifest_json(config))
        self.assertEqual(manifest["schema"], TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA)
        self.assertEqual(manifest["summary"]["caseCount"], 252)
        self.assertEqual(manifest["summary"]["connectorMeasurementRowCount"], 756)
        self.assertTrue(manifest["summary"]["proxyEngineReady"])
        self.assertFalse(manifest["summary"]["exactSegmentedCadReady"])
        self.assertFalse(manifest["summary"]["physicalAccuracyValidated"])
        self.assertIn("externalEngineRunResults", manifest["summary"]["missingEvidence"])
        self.assertIn("filledFidelityMatrixMeasurements", manifest["summary"]["missingEvidence"])
        self.assertEqual(rows[0]["caseId"], "FM_000")
        self.assertEqual(rows[0]["mjcfProxyPath"], "external_fidelity_matrix_mjcf/FM_000.xml")
        self.assertEqual(rows[0]["measurementRows"], 3)
        self.assertEqual(len(csv_rows), manifest["summary"]["caseCount"])
        self.assertEqual(exported["schema"], TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA)
        self.assertEqual(exported["resultSchema"]["fillableCsv"], "two_cell_fidelity_matrix_measurement_template.csv")

    def test_two_cell_external_fidelity_mjcf_run_report_tracks_not_run_status(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        proxy = two_cell_mjcf_proxy_report(config, TwoCellBenchControls(alpha_command=-0.35, z_command=0.35))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xml_dir = root / "external_fidelity_matrix_mjcf"
            xml_dir.mkdir()
            (xml_dir / "FM_000.xml").write_text(str(proxy["xml"]), encoding="utf-8")
            index_path = root / "two_cell_external_fidelity_matrix_mjcf_index.json"
            index_path.write_text(
                json.dumps(
                    {
                        "schema": "rad-sim.two-cell-external-fidelity-matrix-mjcf-index.v1",
                        "caseCount": 1,
                        "directory": "external_fidelity_matrix_mjcf",
                        "files": [
                            {
                                "caseId": "FM_000",
                                "matrixIndex": 0,
                                "path": "external_fidelity_matrix_mjcf/FM_000.xml",
                                "lockMode": "free",
                                "backlash": 0.1,
                                "pinRadius": 0.18,
                                "holeRadius": 0.225,
                                "pinHoleClearanceMm": proxy["scale"]["clearanceMm"],
                                "pinRadiusMm": proxy["scale"]["pinRadiusMm"],
                                "holeRadiusMm": proxy["scale"]["holeRadiusMm"],
                                "xyMmPerModelUnit": proxy["scale"]["xyMmPerModelUnit"],
                                "zMmPerModelUnit": proxy["scale"]["zMmPerModelUnit"],
                                "alphaCommand": -0.35,
                                "zCommand": 0.35,
                                "bodyCount": 2,
                                "actuatorCount": proxy["summary"]["actuatorCount"],
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            report = two_cell_external_fidelity_mjcf_run_report(
                index_path,
                output_dir=root / "run",
                steps=2,
                engine_availability={"mujoco": False},
                case_ids=["FM_000"],
            )
            exported = json.loads(
                export_two_cell_external_fidelity_mjcf_run_json(
                    index_path,
                    output_dir=root / "run-export",
                    steps=2,
                    engine_availability={"mujoco": False},
                )
            )
            self.assertEqual(report["schema"], TWO_CELL_EXTERNAL_FIDELITY_MJCF_RUN_SCHEMA)
            self.assertFalse(report["summary"]["externalFidelityMjcfRunComplete"])
            self.assertFalse(report["summary"]["mujocoAvailable"])
            self.assertEqual(report["summary"]["caseCount"], 1)
            self.assertEqual(report["solver"]["caseIdsRequested"], ["FM_000"])
            self.assertEqual(report["summary"]["caseRunCount"], 0)
            self.assertEqual(report["summary"]["expectedBodyResultCount"], 2)
            self.assertEqual(report["summary"]["expectedConnectorMeasurementRowCount"], 3)
            self.assertEqual(report["summary"]["connectorMeasurementRowCount"], 0)
            self.assertIn("mujocoPythonPackage", report["summary"]["missingEvidence"])
            self.assertIn("externalFidelityMjcfRun", report["summary"]["missingEvidence"])
            self.assertIn("mujocoConnectorMeasurements", report["summary"]["missingEvidence"])
            self.assertEqual(exported["schema"], TWO_CELL_EXTERNAL_FIDELITY_MJCF_RUN_SCHEMA)
            self.assertEqual(two_cell_external_fidelity_mjcf_measurement_rows(report), [])
            report_csv = two_cell_external_fidelity_mjcf_measurements_csv(report)
            exported_csv = export_two_cell_external_fidelity_mjcf_measurements_csv(
                index_path,
                output_dir=root / "csv-export",
                steps=2,
                engine_availability={"mujoco": False},
            )
            self.assertIn("observedLeftXmm", report_csv.splitlines()[0])
            self.assertIn("observedLeftXmm", exported_csv.splitlines()[0])

            cli_out = root / "cli-run"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "rad_sim.run_two_cell_external_fidelity_mjcf",
                    "--index",
                    str(index_path),
                    "--out",
                    str(cli_out),
                    "--steps",
                    "2",
                    "--case-id",
                    "FM_000",
                    "--compare",
                    "--assume-mujoco-unavailable",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("run:", completed.stdout)
            self.assertIn("measurements:", completed.stdout)
            self.assertIn("comparison:", completed.stdout)
            self.assertIn("benchmark:", completed.stdout)
            cli_report = json.loads(
                (cli_out / "two_cell_external_fidelity_matrix_mjcf_run.json").read_text(encoding="utf-8")
            )
            self.assertEqual(cli_report["schema"], TWO_CELL_EXTERNAL_FIDELITY_MJCF_RUN_SCHEMA)
            self.assertTrue((cli_out / "two_cell_external_fidelity_matrix_mjcf_measurements.csv").exists())
            self.assertTrue(
                (
                    root
                    / "cli-run_comparison"
                    / "two_cell_fidelity_matrix_measurement_comparison.json"
                ).exists()
            )
            cli_benchmark = json.loads(
                (
                    root
                    / "cli-run_comparison"
                    / "two_cell_external_fidelity_benchmark_summary.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(cli_benchmark["schema"], TWO_CELL_EXTERNAL_FIDELITY_BENCHMARK_SCHEMA)
            self.assertFalse(cli_benchmark["summary"]["externalRunComplete"])

            direct_benchmark = two_cell_external_fidelity_benchmark_summary(
                report,
                json.loads(
                    (
                        root
                        / "cli-run_comparison"
                        / "two_cell_fidelity_matrix_measurement_comparison.json"
                    ).read_text(encoding="utf-8")
                ),
                {},
            )
            self.assertEqual(direct_benchmark["schema"], TWO_CELL_EXTERNAL_FIDELITY_BENCHMARK_SCHEMA)
            written_benchmark = write_two_cell_external_fidelity_benchmark_summary_artifact(
                cli_out / "two_cell_external_fidelity_matrix_mjcf_run.json",
                root / "cli-run_comparison" / "two_cell_fidelity_matrix_measurement_comparison.json",
                root / "cli-run_comparison" / "two_cell_fidelity_matrix_parameter_calibration.json",
                root / "written-benchmark",
            )
            self.assertTrue(written_benchmark["benchmark"].exists())

    def test_two_cell_fidelity_matrix_measurement_comparison_requires_and_matches_data(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        blank_rows = two_cell_fidelity_matrix_measurement_template_rows(config)
        blank_comparison = compare_two_cell_fidelity_matrix_measurements(blank_rows, config)
        self.assertEqual(blank_comparison["schema"], TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA)
        self.assertEqual(blank_comparison["sampleCount"], 756)
        self.assertEqual(blank_comparison["observedScalarCount"], 0)
        self.assertFalse(blank_comparison["acceptance"]["passesTolerance"])
        self.assertIn("observedNumericFields", blank_comparison["missingEvidence"])

        filled_rows = deepcopy(blank_rows)
        for row in filled_rows:
            for observed, predicted in {
                "observedLeftCellX": "predictedLeftCellX",
                "observedLeftCellY": "predictedLeftCellY",
                "observedLeftCellZ": "predictedLeftCellZ",
                "observedRightCellX": "predictedRightCellX",
                "observedRightCellY": "predictedRightCellY",
                "observedRightCellZ": "predictedRightCellZ",
                "observedLeftAlpha": "predictedLeftAlpha",
                "observedRightAlpha": "predictedRightAlpha",
                "observedLeftTheta": "predictedLeftTheta",
                "observedRightTheta": "predictedRightTheta",
                "observedLeftXmm": "predictedLeftXmm",
                "observedLeftYmm": "predictedLeftYmm",
                "observedLeftZmm": "predictedLeftZmm",
                "observedRightXmm": "predictedRightXmm",
                "observedRightYmm": "predictedRightYmm",
                "observedRightZmm": "predictedRightZmm",
                "observedLateralSlipMm": "predictedLateralSlipMm",
                "observedVerticalSlipMm": "predictedVerticalSlipMm",
                "observedTotalSlipMm": "predictedTotalSlipMm",
            }.items():
                row[observed] = row[predicted]
            row["observedContactMode"] = row["predictedContactMode"]
            if row["lockMode"] in {"right_state_locked", "right_position_locked"}:
                row["lockHeldObserved"] = "held"
            row["measurementSource"] = "synthetic-roundtrip"
        parsed = two_cell_fidelity_matrix_measurements_from_json(json.dumps({"rows": filled_rows}))
        comparison = compare_two_cell_fidelity_matrix_measurements(parsed, config, tolerance=1e-9)
        self.assertEqual(comparison["matchedRowCount"], 756)
        self.assertGreater(comparison["observedScalarCount"], 10000)
        self.assertAlmostEqual(comparison["rmsError"], 0.0)
        self.assertAlmostEqual(comparison["cellRmsError"], 0.0)
        self.assertAlmostEqual(comparison["connectorRmsError"], 0.0)
        self.assertAlmostEqual(comparison["contactModeAccuracy"], 1.0)
        self.assertAlmostEqual(comparison["lockHeldAccuracy"], 1.0)
        self.assertFalse(comparison["acceptance"]["requiresMoreData"])
        self.assertTrue(comparison["acceptance"]["passesTolerance"])

        csv_text = io.StringIO()
        writer = csv.DictWriter(csv_text, fieldnames=list(filled_rows[0]))
        writer.writeheader()
        writer.writerows(filled_rows)
        parsed_csv = two_cell_fidelity_matrix_measurements_from_csv(csv_text.getvalue())
        self.assertEqual(len(parsed_csv), 756)
        comparison_json = json.loads(export_two_cell_fidelity_matrix_measurement_comparison_json(parsed_csv, config))
        comparison_rows = list(
            csv.DictReader(
                io.StringIO(export_two_cell_fidelity_matrix_measurement_comparison_csv(parsed_csv, config))
            )
        )
        self.assertEqual(comparison_json["schema"], TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA)
        self.assertEqual(len(comparison_rows), 756)
        calibration = calibrate_two_cell_fidelity_matrix_parameters(parsed_csv, config, tolerance=1e-9)
        self.assertEqual(calibration["schema"], TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA)
        self.assertTrue(calibration["acceptance"]["readyForReducedProxyCalibration"])
        self.assertFalse(calibration["acceptance"]["physicalAccuracyValidated"])
        self.assertAlmostEqual(calibration["estimates"]["alphaResponseScale"]["estimate"], 1.0)
        self.assertAlmostEqual(calibration["estimates"]["zResponseScale"]["estimate"], 1.0)
        self.assertAlmostEqual(calibration["estimates"]["verticalSlipScale"]["estimate"], 1.0)
        self.assertAlmostEqual(
            calibration["proposedReducedProxyUpdates"]["zCouplingGain"],
            config.z_coupling_gain,
        )
        self.assertAlmostEqual(
            calibration["proposedReducedProxyUpdates"]["holeRadius"],
            config.hole_radius,
        )
        exported_calibration = json.loads(
            export_two_cell_fidelity_matrix_parameter_calibration_json(parsed_csv, config, tolerance=1e-9)
        )
        self.assertEqual(exported_calibration["schema"], TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA)
        self.assertGreaterEqual(len(exported_calibration["axisDiagnostics"]["byHoleRadius"]), 3)

    def test_two_cell_fidelity_matrix_parameter_calibration_tracks_vertical_scale(self):
        config = LatticeConfig(
            rows=1,
            cols=2,
            backlash=0.1,
            pin_radius=0.18,
            hole_radius=0.225,
            z_coupling_gain=0.3,
        )
        rows = deepcopy(two_cell_fidelity_matrix_measurement_template_rows(config))
        copied_fields = {
            "observedLeftCellX": "predictedLeftCellX",
            "observedLeftCellY": "predictedLeftCellY",
            "observedLeftCellZ": "predictedLeftCellZ",
            "observedRightCellX": "predictedRightCellX",
            "observedRightCellY": "predictedRightCellY",
            "observedRightCellZ": "predictedRightCellZ",
            "observedLeftAlpha": "predictedLeftAlpha",
            "observedRightAlpha": "predictedRightAlpha",
            "observedLeftTheta": "predictedLeftTheta",
            "observedRightTheta": "predictedRightTheta",
            "observedLeftXmm": "predictedLeftXmm",
            "observedLeftYmm": "predictedLeftYmm",
            "observedLeftZmm": "predictedLeftZmm",
            "observedRightXmm": "predictedRightXmm",
            "observedRightYmm": "predictedRightYmm",
            "observedRightZmm": "predictedRightZmm",
            "observedLateralSlipMm": "predictedLateralSlipMm",
            "observedVerticalSlipMm": "predictedVerticalSlipMm",
            "observedTotalSlipMm": "predictedTotalSlipMm",
        }
        scaled_observations = {
            "observedLeftCellZ",
            "observedRightCellZ",
            "observedLeftZmm",
            "observedRightZmm",
            "observedVerticalSlipMm",
            "observedTotalSlipMm",
        }
        for row in rows:
            for observed, predicted in copied_fields.items():
                value = row[predicted]
                row[observed] = float(value) * 1.5 if observed in scaled_observations else value
            row["observedContactMode"] = row["predictedContactMode"]
            if row["lockMode"] in {"right_state_locked", "right_position_locked"}:
                row["lockHeldObserved"] = "held"
            row["measurementSource"] = "synthetic-vertical-scale"

        calibration = calibrate_two_cell_fidelity_matrix_parameters(rows, config)
        self.assertGreater(calibration["estimates"]["zResponseScale"]["estimate"], 1.2)
        self.assertGreater(calibration["estimates"]["verticalSlipScale"]["estimate"], 1.2)
        self.assertGreater(
            calibration["proposedReducedProxyUpdates"]["zCouplingGain"],
            config.z_coupling_gain,
        )
        self.assertGreater(
            calibration["proposedReducedProxyUpdates"]["holeRadius"],
            config.hole_radius,
        )

    def test_two_cell_external_fidelity_correction_profile_improves_holdout(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        rows = deepcopy(two_cell_fidelity_matrix_measurement_template_rows(config))
        correction_pairs = {
            "observedLeftCellX": "predictedLeftCellX",
            "observedLeftCellY": "predictedLeftCellY",
            "observedLeftCellZ": "predictedLeftCellZ",
            "observedRightCellX": "predictedRightCellX",
            "observedRightCellY": "predictedRightCellY",
            "observedRightCellZ": "predictedRightCellZ",
            "observedLeftXmm": "predictedLeftXmm",
            "observedLeftYmm": "predictedLeftYmm",
            "observedLeftZmm": "predictedLeftZmm",
            "observedRightXmm": "predictedRightXmm",
            "observedRightYmm": "predictedRightYmm",
            "observedRightZmm": "predictedRightZmm",
            "observedLateralSlipMm": "predictedLateralSlipMm",
            "observedVerticalSlipMm": "predictedVerticalSlipMm",
            "observedTotalSlipMm": "predictedTotalSlipMm",
        }
        for row in rows:
            for observed, predicted in correction_pairs.items():
                row[observed] = 2.0 * float(row[predicted]) + 0.25
            row["observedContactMode"] = row["predictedContactMode"]
            if row["lockMode"] in {"right_state_locked", "right_position_locked"}:
                row["lockHeldObserved"] = "held"
            row["measurementSource"] = "synthetic-affine-external-proxy"

        profile = fit_two_cell_external_fidelity_correction_profile(rows, config, holdout_stride=4)
        self.assertEqual(profile["schema"], TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_SCHEMA)
        self.assertTrue(profile["summary"]["readyForProxyPreview"])
        self.assertTrue(profile["summary"]["coversFullFidelityMatrix"])
        self.assertGreater(profile["holdout"]["baselineRms"], profile["holdout"]["correctedRms"])
        self.assertLess(profile["holdout"]["correctedRms"], profile["holdout"]["baselineRms"] * 0.05)
        self.assertAlmostEqual(profile["fieldCorrections"]["observedRightCellZ"]["slope"], 2.0)
        self.assertAlmostEqual(profile["fieldCorrections"]["observedRightCellZ"]["intercept"], 0.25)
        exported = json.loads(export_two_cell_external_fidelity_correction_profile_json(rows, config, holdout_stride=4))
        self.assertEqual(exported["schema"], TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_SCHEMA)

        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "matrix_measurements.csv"
            out = Path(tmp) / "correction"
            csv_text = io.StringIO()
            writer = csv.DictWriter(csv_text, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
            source.write_text(csv_text.getvalue(), encoding="utf-8")
            files = write_two_cell_external_fidelity_correction_profile_artifact(
                source,
                out,
                holdout_stride=4,
            )
            written = json.loads(files["profile"].read_text(encoding="utf-8"))
            self.assertEqual(written["schema"], TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_SCHEMA)
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "rad_sim.fit_two_cell_external_fidelity_correction",
                    "--input",
                    str(source),
                    "--out",
                    str(out / "cli"),
                    "--holdout-stride",
                    "4",
                ],
                cwd=Path(__file__).resolve().parents[1],
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("profile:", completed.stdout)
            self.assertTrue((out / "cli" / "two_cell_external_fidelity_correction_profile.json").exists())

    def test_two_cell_external_fidelity_correction_application_writes_corrected_preview(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        rows = deepcopy(two_cell_fidelity_matrix_measurement_template_rows(config))
        for row in rows:
            for observed, predicted in {
                "observedLeftCellX": "predictedLeftCellX",
                "observedLeftCellY": "predictedLeftCellY",
                "observedLeftCellZ": "predictedLeftCellZ",
                "observedRightCellX": "predictedRightCellX",
                "observedRightCellY": "predictedRightCellY",
                "observedRightCellZ": "predictedRightCellZ",
                "observedLeftXmm": "predictedLeftXmm",
                "observedLeftYmm": "predictedLeftYmm",
                "observedLeftZmm": "predictedLeftZmm",
                "observedRightXmm": "predictedRightXmm",
                "observedRightYmm": "predictedRightYmm",
                "observedRightZmm": "predictedRightZmm",
                "observedLateralSlipMm": "predictedLateralSlipMm",
                "observedVerticalSlipMm": "predictedVerticalSlipMm",
                "observedTotalSlipMm": "predictedTotalSlipMm",
            }.items():
                row[observed] = 1.25 * float(row[predicted]) - 0.1
            row["observedContactMode"] = row["predictedContactMode"]
            if row["lockMode"] in {"right_state_locked", "right_position_locked"}:
                row["lockHeldObserved"] = "held"
            row["measurementSource"] = "synthetic-application-test"
        profile = fit_two_cell_external_fidelity_correction_profile(rows, config, holdout_stride=3)
        application = apply_two_cell_external_fidelity_correction_profile(profile, rows, config)
        self.assertEqual(application["schema"], TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_APPLICATION_SCHEMA)
        self.assertEqual(application["summary"]["rowCount"], 756)
        self.assertTrue(application["summary"]["readyForProxyPreview"])
        self.assertGreater(application["evaluation"]["baselineRms"], application["evaluation"]["correctedRms"])
        self.assertIn("correctedRightCellZ", application["rows"][0])
        self.assertIn("correctedResidualRightCellZ", application["rows"][0])
        exported = json.loads(export_two_cell_external_fidelity_correction_application_json(profile, rows, config))
        self.assertEqual(exported["schema"], TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_APPLICATION_SCHEMA)
        exported_csv = export_two_cell_external_fidelity_correction_application_csv(application)
        self.assertIn("correctedRightCellZ", exported_csv.splitlines()[0])

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "measurements.csv"
            profile_path = root / "profile.json"
            out = root / "applied"
            csv_text = io.StringIO()
            writer = csv.DictWriter(csv_text, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
            source.write_text(csv_text.getvalue(), encoding="utf-8")
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            files = write_two_cell_external_fidelity_correction_application_artifacts(
                profile_path,
                out,
                measurements_path=source,
            )
            written = json.loads(files["application"].read_text(encoding="utf-8"))
            self.assertEqual(written["schema"], TWO_CELL_EXTERNAL_FIDELITY_CORRECTION_APPLICATION_SCHEMA)
            self.assertTrue(files["application_csv"].exists())
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "rad_sim.apply_two_cell_external_fidelity_correction",
                    "--profile",
                    str(profile_path),
                    "--measurements",
                    str(source),
                    "--out",
                    str(out / "cli"),
                ],
                cwd=Path(__file__).resolve().parents[1],
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("application:", completed.stdout)
            self.assertTrue((out / "cli" / "two_cell_external_fidelity_correction_application.json").exists())

    def test_two_cell_external_fidelity_web_summary_exports_compact_browser_data(self):
        benchmark = {
            "summary": {
                "externalRunComplete": True,
                "coversFullFidelityMatrix": True,
                "passesReducedModelTolerance": False,
                "remainingEvidence": ["contactModeAgreement"],
            },
            "runSummary": {
                "caseCount": 252,
                "expectedConnectorMeasurementRowCount": 756,
            },
            "comparisonSummary": {
                "observedScalarCount": 12,
                "contactModeAccuracy": 0.5,
                "lockHeldAccuracy": 1.0,
            },
            "proposedReducedProxyUpdates": {"zCouplingGain": 0.22},
            "axisSummaries": {"byHoleRadius": [{"holeRadius": 0.225, "residualRms": 2.0}]},
            "worstAxisResiduals": [],
        }
        profile = {
            "cadReference": {"source": "Autodesk A360 public share https://a360.co/4bMlzip"},
            "summary": {
                "readyForProxyPreview": True,
                "holdoutImprovement": 0.07,
                "remainingEvidence": ["benchCoordinateHoldout"],
            },
            "fieldCorrections": {
                "observedRightCellZ": {"predictedField": "predictedRightCellZ", "slope": 1.1, "intercept": 0.02},
                "observedVerticalSlipMm": {"predictedField": "predictedVerticalSlipMm", "slope": 1.9, "intercept": 1.5},
            },
        }
        application = {
            "summary": {
                "rowCount": 3,
                "measurementRowCount": 3,
                "readyForProxyPreview": True,
                "remainingEvidence": ["segmentedCadContactModel"],
            },
            "evaluation": {
                "caseCount": 1,
                "observedScalarCount": 12,
                "baselineRms": 2.0,
                "correctedRms": 1.5,
                "improvement": 0.25,
                "contactModeAccuracy": 0.5,
                "lockHeldAccuracy": 1.0,
            },
            "rows": [
                {
                    "caseId": "FM_000",
                    "matrixIndex": 0,
                    "connector": "upper",
                    "predictedRightCellZ": 0.1,
                    "correctedRightCellZ": 0.13,
                    "observedRightCellZ": 0.12,
                    "predictedVerticalSlipMm": 1.0,
                    "correctedVerticalSlipMm": 3.4,
                    "observedVerticalSlipMm": 3.1,
                }
            ],
        }
        summary = build_two_cell_external_fidelity_web_summary(benchmark, profile, application)
        self.assertEqual(summary["schema"], WEB_TWO_CELL_EXTERNAL_FIDELITY_SCHEMA)
        self.assertTrue(summary["summary"]["externalRunComplete"])
        self.assertEqual(summary["summary"]["caseCount"], 252)
        self.assertEqual(summary["summary"]["baselineRms"], 2.0)
        self.assertEqual(summary["summary"]["correctedRms"], 1.5)
        self.assertFalse(summary["summary"]["physicalAccuracyValidated"])
        self.assertEqual(len(summary["representativeRows"]), 1)
        self.assertIn("segmentedCadContactModel", summary["summary"]["remainingEvidence"])

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            benchmark_path = root / "benchmark.json"
            profile_path = root / "profile.json"
            application_path = root / "application.json"
            out_path = root / "web" / "two_cell_external_fidelity_summary.json"
            benchmark_path.write_text(json.dumps(benchmark), encoding="utf-8")
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            application_path.write_text(json.dumps(application), encoding="utf-8")
            files = write_two_cell_external_fidelity_web_summary(
                benchmark_path,
                profile_path,
                application_path,
                out_path,
            )
            written = json.loads(files["json"].read_text(encoding="utf-8"))
            self.assertEqual(written["schema"], WEB_TWO_CELL_EXTERNAL_FIDELITY_SCHEMA)
            self.assertIn("window.RAD_EXTERNAL_FIDELITY_SUMMARY", files["js"].read_text(encoding="utf-8"))

    def test_two_cell_fidelity_matrix_measurement_comparison_writer_creates_reports(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "matrix_measurements.csv"
            out = Path(tmp) / "comparison"
            source.write_text(export_two_cell_fidelity_matrix_measurement_template_csv(config), encoding="utf-8")
            files = write_two_cell_fidelity_matrix_measurement_comparison_artifacts(source, out)
            comparison = json.loads(files["comparison"].read_text(encoding="utf-8"))
            calibration = json.loads(files["calibration"].read_text(encoding="utf-8"))
            rows = list(csv.DictReader(io.StringIO(files["comparison_csv"].read_text(encoding="utf-8"))))
            self.assertEqual(comparison["schema"], TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA)
            self.assertEqual(calibration["schema"], TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA)
            self.assertEqual(len(rows), 756)

    def test_two_cell_fidelity_matrix_measurement_comparison_cli_writes_artifacts(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "matrix_measurements.csv"
            out = Path(tmp) / "comparison"
            source.write_text(export_two_cell_fidelity_matrix_measurement_template_csv(config), encoding="utf-8")
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "rad_sim.compare_two_cell_fidelity_matrix_measurements",
                    "--input",
                    str(source),
                    "--out",
                    str(out),
                    "--backlash",
                    "0.1",
                    "--pin-radius",
                    "0.18",
                    "--hole-radius",
                    "0.225",
                ],
                cwd=Path(__file__).resolve().parents[1],
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("comparison:", completed.stdout)
            self.assertIn("calibration:", completed.stdout)
            comparison = json.loads((out / "two_cell_fidelity_matrix_measurement_comparison.json").read_text(encoding="utf-8"))
            calibration = json.loads((out / "two_cell_fidelity_matrix_parameter_calibration.json").read_text(encoding="utf-8"))
            rows = list(
                csv.DictReader(
                    io.StringIO((out / "two_cell_fidelity_matrix_measurement_comparison.csv").read_text(encoding="utf-8"))
                )
            )
            self.assertEqual(comparison["schema"], TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_COMPARISON_SCHEMA)
            self.assertEqual(calibration["schema"], TWO_CELL_FIDELITY_MATRIX_PARAMETER_CALIBRATION_SCHEMA)
            self.assertEqual(len(rows), 756)

    def test_two_cell_connector_measurement_template_targets_suite_connectors(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        template = two_cell_connector_measurement_template(config)
        rows = two_cell_connector_measurement_template_rows(config)
        self.assertEqual(template["schema"], TWO_CELL_CONNECTOR_MEASUREMENT_TEMPLATE_SCHEMA)
        self.assertEqual(template["summary"]["caseCount"], 12)
        self.assertEqual(template["summary"]["connectorRowCount"], 36)
        self.assertEqual(len(rows), 36)
        self.assertEqual({row["connector"] for row in rows}, {"upper", "middle", "lower"})
        first = rows[0]
        self.assertIn("predictedVerticalSlipMm", first)
        self.assertIn("predictedContactMode", first)
        self.assertEqual(first["observedVerticalSlipMm"], "")
        self.assertEqual(first["lockHeldObserved"], "")

        exported = json.loads(export_two_cell_connector_measurement_template_json(config))
        exported_rows = list(csv.DictReader(io.StringIO(export_two_cell_connector_measurement_template_csv(config))))
        self.assertEqual(exported["schema"], TWO_CELL_CONNECTOR_MEASUREMENT_TEMPLATE_SCHEMA)
        self.assertEqual(len(exported_rows), 36)
        self.assertIn("observedRightZmm", exported_rows[0])

    def test_two_cell_connector_measurement_comparison_requires_observed_connector_data(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        blank_csv = export_two_cell_connector_measurement_template_csv(config)
        measurements = two_cell_connector_measurements_from_csv(blank_csv)
        comparison = compare_two_cell_connector_measurements(measurements, config)
        self.assertEqual(comparison["schema"], TWO_CELL_CONNECTOR_MEASUREMENT_COMPARISON_SCHEMA)
        self.assertEqual(comparison["sampleCount"], 36)
        self.assertEqual(comparison["observedConnectorRowCount"], 0)
        self.assertEqual(comparison["observedScalarCount"], 0)
        self.assertTrue(comparison["acceptance"]["requiresMoreData"])
        self.assertFalse(comparison["acceptance"]["readyForConnectorCalibration"])
        self.assertIn("observedConnectorRows", comparison["missingEvidence"])

    def test_two_cell_connector_measurement_comparison_is_zero_for_model_generated_rows(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        rows = deepcopy(two_cell_connector_measurement_template_rows(config))
        for row in rows:
            row["observedLeftXmm"] = row["predictedLeftXmm"]
            row["observedLeftYmm"] = row["predictedLeftYmm"]
            row["observedLeftZmm"] = row["predictedLeftZmm"]
            row["observedRightXmm"] = row["predictedRightXmm"]
            row["observedRightYmm"] = row["predictedRightYmm"]
            row["observedRightZmm"] = row["predictedRightZmm"]
            row["observedLateralSlipMm"] = row["predictedLateralSlipMm"]
            row["observedVerticalSlipMm"] = row["predictedVerticalSlipMm"]
            row["observedTotalSlipMm"] = row["predictedTotalSlipMm"]
            row["observedContactMode"] = row["predictedContactMode"]
            if row["lockMode"] in {"right_state_locked", "right_position_locked"}:
                row["lockHeldObserved"] = "held"
        measurements = two_cell_connector_measurements_from_json(json.dumps({"rows": rows}))
        comparison = compare_two_cell_connector_measurements(measurements, config, tolerance=1e-9)
        self.assertEqual(comparison["schema"], TWO_CELL_CONNECTOR_MEASUREMENT_COMPARISON_SCHEMA)
        self.assertEqual(comparison["matchedRowCount"], 36)
        self.assertEqual(comparison["observedScalarCount"], 324)
        self.assertAlmostEqual(comparison["markerRmsErrorMm"], 0.0)
        self.assertAlmostEqual(comparison["slipRmsErrorMm"], 0.0)
        self.assertAlmostEqual(comparison["contactModeAccuracy"], 1.0)
        self.assertAlmostEqual(comparison["lockHeldAccuracy"], 1.0)
        self.assertFalse(comparison["acceptance"]["requiresMoreData"])
        self.assertTrue(comparison["acceptance"]["passesTolerance"])

        comparison_json = json.loads(export_two_cell_connector_measurement_comparison_json(measurements, config))
        comparison_rows = list(
            csv.DictReader(io.StringIO(export_two_cell_connector_measurement_comparison_csv(measurements, config)))
        )
        self.assertEqual(comparison_json["schema"], TWO_CELL_CONNECTOR_MEASUREMENT_COMPARISON_SCHEMA)
        self.assertEqual(len(comparison_rows), 36)
        self.assertEqual(comparison_rows[0]["status"], "matched")

    def test_two_cell_connector_measurement_comparison_detects_contact_and_slip_mismatch(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        first = two_cell_connector_measurement_template_rows(config)[0]
        measurement = TwoCellConnectorMeasurement(
            case_id=first["caseId"],
            connector=first["connector"],
            vertical_slip_mm=float(first["predictedVerticalSlipMm"]) + 0.25,
            contact_mode_observed="wrong-contact",
        )
        comparison = compare_two_cell_connector_measurements([measurement], config, tolerance=1e-9)
        self.assertGreater(comparison["slipMaxAbsErrorMm"], 0.2)
        self.assertLess(comparison["contactModeAccuracy"], 1.0)
        self.assertFalse(comparison["acceptance"]["passesTolerance"])
        self.assertIn("connectorToleranceAgreement", comparison["missingEvidence"])
        self.assertIn("contactModeAgreement", comparison["missingEvidence"])

    def test_two_cell_connector_measurement_comparison_cli_writes_artifacts(self):
        config = LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225)
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "connector_measurements.csv"
            out = Path(tmp) / "comparison"
            source.write_text(export_two_cell_connector_measurement_template_csv(config), encoding="utf-8")
            files = write_two_cell_connector_measurement_comparison_artifacts(source, out)
            comparison = json.loads(files["comparison"].read_text(encoding="utf-8"))
            rows = list(csv.DictReader(io.StringIO(files["comparison_csv"].read_text(encoding="utf-8"))))
            self.assertEqual(comparison["schema"], TWO_CELL_CONNECTOR_MEASUREMENT_COMPARISON_SCHEMA)
            self.assertEqual(len(rows), 36)

    def test_two_cell_external_template_and_comparison_roundtrip(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        template = two_cell_external_results_template(config)
        self.assertEqual(template["schema"], TWO_CELL_EXTERNAL_TEMPLATE_SCHEMA)
        self.assertEqual(template["summary"]["rowCount"], 10)
        self.assertEqual(len(two_cell_external_results_template_rows(config)), 10)
        csv_template = export_two_cell_external_results_template_csv(config)
        json_template = json.loads(export_two_cell_external_results_template_json(config))
        self.assertIn("externalZ", csv_template)
        self.assertEqual(json_template["schema"], TWO_CELL_EXTERNAL_TEMPLATE_SCHEMA)

        rows = deepcopy(template["rows"])
        for row in rows:
            row["externalX"] = row["expectedX"]
            row["externalY"] = row["expectedY"]
            row["externalZ"] = row["expectedZ"]
            row["externalAlpha"] = row["expectedAlpha"]
            row["externalTheta"] = row["expectedTheta"]
            row["sourceEngine"] = "synthetic-external"
        parsed_json = two_cell_external_results_from_json(json.dumps({"rows": rows}))
        comparison = compare_two_cell_external_results(parsed_json, config, tolerance=1e-9)
        self.assertEqual(comparison["schema"], TWO_CELL_EXTERNAL_COMPARISON_SCHEMA)
        self.assertTrue(comparison["summary"]["externalComparisonReady"])
        self.assertEqual(comparison["summary"]["missingExternalResultCount"], 0)
        self.assertAlmostEqual(comparison["summary"]["rmsScalarError"], 0.0)

        comparison_json = json.loads(export_two_cell_external_comparison_json(parsed_json, config, tolerance=1e-9))
        self.assertEqual(comparison_json["summary"], comparison["summary"])

        csv_text = io.StringIO()
        writer = csv.DictWriter(csv_text, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        parsed_csv = two_cell_external_results_from_csv(csv_text.getvalue())
        self.assertEqual(len(parsed_csv), 10)

    def test_two_cell_external_comparison_detects_lock_and_tolerance_mismatch(self):
        config = LatticeConfig(rows=1, cols=2, pin_radius=0.18, hole_radius=0.225)
        rows = deepcopy(two_cell_external_results_template(config)["rows"])
        for row in rows:
            row["externalX"] = row["expectedX"]
            row["externalY"] = row["expectedY"]
            row["externalZ"] = row["expectedZ"]
            row["externalAlpha"] = row["expectedAlpha"]
            row["externalTheta"] = row["expectedTheta"]
        target = next(row for row in rows if row["caseId"] == "right_position_locked" and row["bodyId"] == "right")
        target["externalZ"] = float(target["expectedZ"]) + 0.08
        comparison = compare_two_cell_external_results(rows, config, tolerance=1e-6)
        self.assertFalse(comparison["summary"]["externalComparisonReady"])
        self.assertEqual(comparison["summary"]["lockViolationCount"], 1)
        self.assertGreater(comparison["summary"]["maxPositionError"], 0.01)
        self.assertIn("externalTolerance", comparison["summary"]["missingEvidence"])
        self.assertIn("lockConstraintAgreement", comparison["summary"]["missingEvidence"])

        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "external.csv"
            with input_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            files = write_two_cell_external_comparison_artifacts(input_path, Path(tmp) / "reports")
            saved = json.loads(files["comparison"].read_text(encoding="utf-8"))
            self.assertEqual(saved["schema"], TWO_CELL_EXTERNAL_COMPARISON_SCHEMA)
            self.assertFalse(saved["summary"]["externalComparisonReady"])

    def test_two_cell_physics_validation_report_collects_internal_and_external_gates(self):
        report = two_cell_physics_validation_report(LatticeConfig(rows=1, cols=2))
        self.assertEqual(report["schema"], TWO_CELL_PHYSICS_VALIDATION_SCHEMA)
        self.assertEqual(report["summary"]["status"], "needs-external-physical-validation")
        self.assertTrue(report["summary"]["internalValidationPass"])
        self.assertFalse(report["summary"]["externalValidationPass"])
        test_ids = {test["testId"] for test in report["tests"]}
        self.assertEqual(
            test_ids,
            {
                "cad_layout_envelope",
                "hole_clearance_monotonicity",
                "quasistatic_sweep_success",
                "connector_clearance_contact_dieoff",
                "state_lock_alpha_hold_vertical_free",
                "position_lock_fixture_hold",
                "gravity_clearance_sag",
                "actuation_polarity",
                "mjcf_proxy_ready",
                "mujoco_proxy_run_gate",
                "external_result_gate",
            },
        )
        internal_tests = [
            test
            for test in report["tests"]
            if test["testId"] not in {"external_result_gate", "mujoco_proxy_run_gate"}
        ]
        self.assertTrue(all(test["pass"] for test in internal_tests))
        mujoco_gate = next(test for test in report["tests"] if test["testId"] == "mujoco_proxy_run_gate")
        if mujoco_gate["metrics"]["mujocoAvailable"]:
            self.assertIn(mujoco_gate["status"], {"pass", "incomplete"})
        else:
            self.assertEqual(mujoco_gate["status"], "incomplete")
            self.assertIn("mujocoPythonPackage", mujoco_gate["missingEvidence"])
        external_gate = next(test for test in report["tests"] if test["testId"] == "external_result_gate")
        self.assertEqual(external_gate["status"], "incomplete")
        self.assertIn("externalResultRows", external_gate["missingEvidence"])
        self.assertIn("externalMuJoCoRun", report["summary"]["missingEvidence"])

        report_json = json.loads(export_two_cell_physics_validation_report_json())
        report_csv = list(csv.DictReader(io.StringIO(export_two_cell_physics_validation_report_csv())))
        self.assertEqual(report_json["schema"], TWO_CELL_PHYSICS_VALIDATION_SCHEMA)
        self.assertEqual(len(report_csv), report["summary"]["testCount"])
        self.assertIn("mjcf_proxy_ready", {row["testId"] for row in report_csv})

    def test_two_cell_cad_contact_decomposition_spec_defines_engine_contract(self):
        spec = two_cell_cad_contact_decomposition_spec(LatticeConfig(rows=1, cols=2), ring_segments=8)
        self.assertEqual(spec["schema"], TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA)
        self.assertEqual(spec["summary"]["holeCount"], 16)
        self.assertEqual(spec["summary"]["twoCellConnectorCount"], 3)
        self.assertEqual(spec["summary"]["ringSegments"], 8)
        self.assertGreaterEqual(spec["summary"]["collisionPrimitiveCount"], 16 * 8)
        self.assertFalse(spec["summary"]["physicalAccuracyValidated"])
        self.assertFalse(spec["summary"]["canBuildExactContactModel"])
        self.assertIn("convex-or-analytic-hole-wall-contact-primitives", spec["summary"]["missingEvidence"])
        self.assertEqual(len(spec["holeDecomposition"][0]["ringSegments"]), 8)
        self.assertEqual(
            spec["holeDecomposition"][0]["recommendedCollisionPrimitive"],
            "convex_ring_sector_or_capsule_wall",
        )
        self.assertEqual(spec["twoCellContactPairs"][0]["leftCell"], "left")

        exported = json.loads(export_two_cell_cad_contact_decomposition_json(ring_segments=8))
        rows = list(csv.DictReader(io.StringIO(export_two_cell_cad_contact_decomposition_csv(ring_segments=8))))
        self.assertEqual(exported["schema"], TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA)
        self.assertGreaterEqual(len(rows), 16 * 8)
        self.assertIn("convex_ring_sector_or_capsule_wall", {row["primitive"] for row in rows})

    def test_two_cell_exact_contact_handoff_plan_prioritizes_real_contact_cases(self):
        plan = two_cell_exact_contact_handoff_plan(LatticeConfig(rows=1, cols=2), priority_limit=8)
        self.assertEqual(plan["schema"], TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA)
        self.assertEqual(plan["summary"]["status"], "needs-segmented-cad-contact-inputs")
        self.assertFalse(plan["summary"]["physicalAccuracyValidated"])
        self.assertFalse(plan["summary"]["canRunExactContact"])
        self.assertGreater(plan["summary"]["caseCount"], 0)
        self.assertGreater(plan["summary"]["transitionCaseCount"], 0)
        self.assertIn("segmentedCADBodies", plan["summary"]["missingEvidence"])
        self.assertIn("externalRigidBodyContactResults", plan["summary"]["missingEvidence"])
        self.assertIn("upper_free_cell_body", plan["caseRows"][0]["requiredCadAssets"])
        self.assertIn("bench_coordinate_truth", plan["caseRows"][0]["requiredMeasuredInputs"])
        self.assertEqual(plan["caseRows"][0]["expectedConnectorCount"], 3)

        exported = json.loads(export_two_cell_exact_contact_handoff_plan_json(priority_limit=8))
        rows = list(csv.DictReader(io.StringIO(export_two_cell_exact_contact_handoff_plan_csv(priority_limit=8))))
        self.assertEqual(exported["schema"], TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA)
        self.assertEqual(len(rows), plan["summary"]["caseCount"])
        self.assertIn("requiredCadAssets", rows[0])
        self.assertIn("canRunExactContact", rows[0])

    def test_two_cell_physical_packet_exports_lab_measurement_template(self):
        packet = two_cell_physical_test_packet(LatticeConfig(rows=1, cols=2))
        self.assertEqual(packet["schema"], TWO_CELL_PACKET_SCHEMA)
        self.assertEqual(packet["cadLayout"]["schema"], CAD_RAD_CELL_LAYOUT_SCHEMA)
        self.assertEqual(packet["clearanceSweep"]["schema"], TWO_CELL_SWEEP_SCHEMA)
        self.assertEqual(packet["actuationSweep"]["schema"], TWO_CELL_ACTUATION_SWEEP_SCHEMA)
        self.assertEqual(packet["quasistaticPhysics"]["schema"], TWO_CELL_QUASISTATIC_SCHEMA)
        self.assertEqual(packet["quasistaticSweep"]["schema"], TWO_CELL_QUASISTATIC_SWEEP_SCHEMA)
        self.assertEqual(packet["physicalSimulationSuite"]["schema"], TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA)
        self.assertEqual(packet["physicalResponseAtlas"]["schema"], TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA)
        self.assertEqual(packet["connectorMeasurementTemplate"]["schema"], TWO_CELL_CONNECTOR_MEASUREMENT_TEMPLATE_SCHEMA)
        self.assertEqual(
            packet["fidelityMatrixMeasurementTemplate"]["schema"],
            TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA,
        )
        self.assertEqual(
            packet["externalFidelityMatrixManifest"]["schema"],
            TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA,
        )
        self.assertEqual(packet["cadContactDecomposition"]["schema"], TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA)
        self.assertEqual(packet["exactContactHandoffPlan"]["schema"], TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA)
        self.assertEqual(packet["connectorContact"]["schema"], TWO_CELL_CONNECTOR_CONTACT_SCHEMA)
        self.assertEqual(packet["connectorContactSweep"]["schema"], TWO_CELL_CONNECTOR_CONTACT_SWEEP_SCHEMA)
        self.assertEqual(packet["segmentedCadReadiness"]["schema"], TWO_CELL_SEGMENTED_CAD_READINESS_SCHEMA)
        self.assertEqual(packet["segmentedCadIntakeTemplates"]["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_TEMPLATES_SCHEMA)
        self.assertEqual(packet["segmentedCadIntakeValidation"]["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_VALIDATION_SCHEMA)
        self.assertEqual(packet["mjcfProxy"]["schema"], TWO_CELL_MJCF_SCHEMA)
        self.assertEqual(packet["mujocoProxyRun"]["schema"], TWO_CELL_MUJOCO_RUN_SCHEMA)
        self.assertEqual(packet["externalResultsTemplate"]["schema"], TWO_CELL_EXTERNAL_TEMPLATE_SCHEMA)
        self.assertEqual(packet["physicsValidationReport"]["schema"], TWO_CELL_PHYSICS_VALIDATION_SCHEMA)
        self.assertIn("rightTheta", packet["measurementTemplateColumns"])
        self.assertIn("predictedRightZ", packet["measurementTemplatePredictedColumns"])
        self.assertEqual(
            {case["caseId"] for case in packet["cases"]},
            {
                "free_contract_lift",
                "free_expand_pushdown",
                "right_state_locked",
                "right_position_locked",
                "left_free_neighbor_residual",
            },
        )
        bench_json = json.loads(export_two_cell_bench_json())
        packet_json = json.loads(export_two_cell_physical_test_packet_json())
        sweep_rows = list(csv.DictReader(io.StringIO(export_two_cell_backlash_sweep_csv())))
        template_rows = list(csv.DictReader(io.StringIO(export_two_cell_measurement_template_csv())))
        self.assertEqual(bench_json["schema"], TWO_CELL_BENCH_SCHEMA)
        self.assertEqual(packet_json["schema"], TWO_CELL_PACKET_SCHEMA)
        self.assertEqual(len(sweep_rows), 9)
        self.assertEqual(len(template_rows), 5)
        self.assertIn("predictedRightZ", template_rows[0])
        self.assertEqual(len(two_cell_measurement_template_rows()), 5)

    def test_two_cell_bench_packet_writer_creates_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = write_two_cell_bench_packet_artifacts(
                tmp,
                backlash=0.03,
                hole_sweep_steps=5,
            )
            self.assertEqual(
                set(files),
                {
                    "cad_archive_audit",
                    "cad_layout",
                    "cad_reference_profile",
                    "one_cell_cad_export_audit",
                    "segmented_cad_readiness",
                    "segmented_cad_readiness_csv",
                    "segmented_cad_intake_templates",
                    "segmented_cad_intake_templates_csv",
                    "segmented_cad_intake_validation",
                    "segmented_cad_intake_validation_csv",
                    "segmented_cad_template_dir",
                    "bench",
                    "sweep",
                    "actuation_sweep",
                    "actuation_report",
                    "mjcf_xml",
                    "mjcf_report",
                    "mujoco_run",
                    "physical_suite",
                    "physical_suite_csv",
                    "physical_fidelity_matrix",
                    "physical_fidelity_matrix_csv",
                    "physical_response_atlas",
                    "physical_response_atlas_csv",
                    "contact_phase_map",
                    "contact_phase_map_csv",
                    "radius_backlash_phase_diagram",
                    "radius_backlash_phase_diagram_csv",
                    "radius_backlash_transition_report",
                    "radius_backlash_transition_report_csv",
                    "cad_contact_decomposition",
                    "cad_contact_decomposition_csv",
                    "exact_contact_handoff_plan",
                    "exact_contact_handoff_plan_csv",
                    "fidelity_matrix_measurement_template",
                    "fidelity_matrix_measurement_template_csv",
                    "external_fidelity_matrix_manifest",
                    "external_fidelity_matrix_manifest_csv",
                    "external_fidelity_matrix_mjcf_dir",
                    "external_fidelity_matrix_mjcf_index",
                    "external_fidelity_matrix_mjcf_run",
                    "response_atlas_priority_mjcf_dir",
                    "response_atlas_priority_mjcf_index",
                    "response_atlas_priority_mjcf_run",
                    "quasistatic",
                    "quasistatic_sweep",
                    "quasistatic_report",
                    "connector_contact",
                    "connector_contact_csv",
                    "connector_contact_sweep",
                    "connector_contact_sweep_csv",
                    "connector_measurement_template",
                    "connector_measurement_template_csv",
                    "external_template",
                    "external_template_json",
                    "physics_validation",
                    "physics_validation_csv",
                    "packet",
                    "template",
                },
            )
            cad_archive_audit_report = json.loads(files["cad_archive_audit"].read_text(encoding="utf-8"))
            cad_layout = json.loads(files["cad_layout"].read_text(encoding="utf-8"))
            cad_reference_profile = json.loads(files["cad_reference_profile"].read_text(encoding="utf-8"))
            one_cell_cad_export = json.loads(files["one_cell_cad_export_audit"].read_text(encoding="utf-8"))
            segmented_cad_readiness = json.loads(files["segmented_cad_readiness"].read_text(encoding="utf-8"))
            segmented_cad_templates = json.loads(files["segmented_cad_intake_templates"].read_text(encoding="utf-8"))
            segmented_cad_validation = json.loads(files["segmented_cad_intake_validation"].read_text(encoding="utf-8"))
            bench = json.loads(files["bench"].read_text(encoding="utf-8"))
            actuation_report = json.loads(files["actuation_report"].read_text(encoding="utf-8"))
            mjcf_report = json.loads(files["mjcf_report"].read_text(encoding="utf-8"))
            mujoco_run = json.loads(files["mujoco_run"].read_text(encoding="utf-8"))
            physical_suite = json.loads(files["physical_suite"].read_text(encoding="utf-8"))
            physical_fidelity_matrix = json.loads(files["physical_fidelity_matrix"].read_text(encoding="utf-8"))
            physical_response_atlas = json.loads(files["physical_response_atlas"].read_text(encoding="utf-8"))
            contact_phase_map = json.loads(files["contact_phase_map"].read_text(encoding="utf-8"))
            radius_backlash_phase_diagram = json.loads(files["radius_backlash_phase_diagram"].read_text(encoding="utf-8"))
            radius_backlash_transition_report = json.loads(
                files["radius_backlash_transition_report"].read_text(encoding="utf-8")
            )
            cad_contact_decomposition = json.loads(files["cad_contact_decomposition"].read_text(encoding="utf-8"))
            exact_contact_handoff_plan = json.loads(files["exact_contact_handoff_plan"].read_text(encoding="utf-8"))
            fidelity_matrix_measurement_template = json.loads(
                files["fidelity_matrix_measurement_template"].read_text(encoding="utf-8")
            )
            external_fidelity_matrix_manifest = json.loads(
                files["external_fidelity_matrix_manifest"].read_text(encoding="utf-8")
            )
            external_fidelity_matrix_mjcf_index = json.loads(
                files["external_fidelity_matrix_mjcf_index"].read_text(encoding="utf-8")
            )
            external_fidelity_matrix_mjcf_run = json.loads(
                files["external_fidelity_matrix_mjcf_run"].read_text(encoding="utf-8")
            )
            response_atlas_priority_mjcf_index = json.loads(
                files["response_atlas_priority_mjcf_index"].read_text(encoding="utf-8")
            )
            response_atlas_priority_mjcf_run = json.loads(
                files["response_atlas_priority_mjcf_run"].read_text(encoding="utf-8")
            )
            quasistatic = json.loads(files["quasistatic"].read_text(encoding="utf-8"))
            quasistatic_report = json.loads(files["quasistatic_report"].read_text(encoding="utf-8"))
            connector_contact = json.loads(files["connector_contact"].read_text(encoding="utf-8"))
            connector_contact_sweep = json.loads(files["connector_contact_sweep"].read_text(encoding="utf-8"))
            connector_measurement_template = json.loads(files["connector_measurement_template"].read_text(encoding="utf-8"))
            external_template = json.loads(files["external_template_json"].read_text(encoding="utf-8"))
            physics_validation = json.loads(files["physics_validation"].read_text(encoding="utf-8"))
            packet = json.loads(files["packet"].read_text(encoding="utf-8"))
            sweep_rows = list(csv.DictReader(io.StringIO(files["sweep"].read_text(encoding="utf-8"))))
            segmented_cad_rows = list(
                csv.DictReader(io.StringIO(files["segmented_cad_readiness_csv"].read_text(encoding="utf-8")))
            )
            segmented_template_rows = list(
                csv.DictReader(io.StringIO(files["segmented_cad_intake_templates_csv"].read_text(encoding="utf-8")))
            )
            segmented_validation_rows = list(
                csv.DictReader(io.StringIO(files["segmented_cad_intake_validation_csv"].read_text(encoding="utf-8")))
            )
            actuation_rows = list(csv.DictReader(io.StringIO(files["actuation_sweep"].read_text(encoding="utf-8"))))
            quasistatic_rows = list(csv.DictReader(io.StringIO(files["quasistatic_sweep"].read_text(encoding="utf-8"))))
            connector_contact_rows = list(
                csv.DictReader(io.StringIO(files["connector_contact_csv"].read_text(encoding="utf-8")))
            )
            connector_contact_sweep_rows = list(
                csv.DictReader(io.StringIO(files["connector_contact_sweep_csv"].read_text(encoding="utf-8")))
            )
            connector_measurement_rows = list(
                csv.DictReader(io.StringIO(files["connector_measurement_template_csv"].read_text(encoding="utf-8")))
            )
            physical_suite_rows = list(
                csv.DictReader(io.StringIO(files["physical_suite_csv"].read_text(encoding="utf-8")))
            )
            physical_fidelity_rows = list(
                csv.DictReader(io.StringIO(files["physical_fidelity_matrix_csv"].read_text(encoding="utf-8")))
            )
            physical_response_rows = list(
                csv.DictReader(io.StringIO(files["physical_response_atlas_csv"].read_text(encoding="utf-8")))
            )
            contact_phase_rows = list(
                csv.DictReader(io.StringIO(files["contact_phase_map_csv"].read_text(encoding="utf-8")))
            )
            radius_backlash_phase_rows = list(
                csv.DictReader(io.StringIO(files["radius_backlash_phase_diagram_csv"].read_text(encoding="utf-8")))
            )
            radius_backlash_transition_rows = list(
                csv.DictReader(io.StringIO(files["radius_backlash_transition_report_csv"].read_text(encoding="utf-8")))
            )
            cad_contact_decomposition_rows = list(
                csv.DictReader(io.StringIO(files["cad_contact_decomposition_csv"].read_text(encoding="utf-8")))
            )
            exact_contact_handoff_rows = list(
                csv.DictReader(io.StringIO(files["exact_contact_handoff_plan_csv"].read_text(encoding="utf-8")))
            )
            fidelity_matrix_measurement_rows = list(
                csv.DictReader(io.StringIO(files["fidelity_matrix_measurement_template_csv"].read_text(encoding="utf-8")))
            )
            external_fidelity_manifest_rows = list(
                csv.DictReader(io.StringIO(files["external_fidelity_matrix_manifest_csv"].read_text(encoding="utf-8")))
            )
            external_rows = list(csv.DictReader(io.StringIO(files["external_template"].read_text(encoding="utf-8"))))
            physics_validation_rows = list(
                csv.DictReader(io.StringIO(files["physics_validation_csv"].read_text(encoding="utf-8")))
            )
            template_rows = list(csv.DictReader(io.StringIO(files["template"].read_text(encoding="utf-8"))))
            self.assertEqual(cad_archive_audit_report["schema"], CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA)
            self.assertEqual(cad_archive_audit_report["summary"]["brepEntryCount"], 2)
            self.assertEqual(cad_layout["schema"], CAD_RAD_CELL_LAYOUT_SCHEMA)
            self.assertEqual(cad_reference_profile["schema"], CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA)
            self.assertEqual(cad_reference_profile["hardwareProfile"]["schema"], HARDWARE_PROFILE_SCHEMA)
            self.assertEqual(one_cell_cad_export["schema"], ONE_CELL_CAD_EXPORT_AUDIT_SCHEMA)
            self.assertFalse(one_cell_cad_export["summary"]["physicalAccuracyValidated"])
            self.assertEqual(segmented_cad_readiness["schema"], TWO_CELL_SEGMENTED_CAD_READINESS_SCHEMA)
            self.assertEqual(segmented_cad_templates["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_TEMPLATES_SCHEMA)
            self.assertEqual(segmented_cad_validation["schema"], TWO_CELL_SEGMENTED_CAD_INTAKE_VALIDATION_SCHEMA)
            self.assertEqual(bench["schema"], TWO_CELL_BENCH_SCHEMA)
            self.assertEqual(actuation_report["schema"], TWO_CELL_ACTUATION_SWEEP_SCHEMA)
            self.assertEqual(mjcf_report["schema"], TWO_CELL_MJCF_SCHEMA)
            self.assertEqual(mujoco_run["schema"], TWO_CELL_MUJOCO_RUN_SCHEMA)
            self.assertEqual(physical_suite["schema"], TWO_CELL_PHYSICAL_SIMULATION_SUITE_SCHEMA)
            self.assertEqual(physical_fidelity_matrix["schema"], TWO_CELL_PHYSICAL_FIDELITY_MATRIX_SCHEMA)
            self.assertEqual(physical_response_atlas["schema"], TWO_CELL_PHYSICAL_RESPONSE_ATLAS_SCHEMA)
            self.assertEqual(contact_phase_map["schema"], TWO_CELL_CONTACT_PHASE_MAP_SCHEMA)
            self.assertEqual(radius_backlash_phase_diagram["schema"], TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA)
            self.assertEqual(
                radius_backlash_transition_report["schema"],
                TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA,
            )
            self.assertEqual(cad_contact_decomposition["schema"], TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA)
            self.assertEqual(exact_contact_handoff_plan["schema"], TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA)
            self.assertEqual(
                fidelity_matrix_measurement_template["schema"],
                TWO_CELL_FIDELITY_MATRIX_MEASUREMENT_TEMPLATE_SCHEMA,
            )
            self.assertEqual(
                external_fidelity_matrix_manifest["schema"],
                TWO_CELL_EXTERNAL_FIDELITY_MANIFEST_SCHEMA,
            )
            self.assertEqual(
                external_fidelity_matrix_mjcf_index["schema"],
                "rad-sim.two-cell-external-fidelity-matrix-mjcf-index.v1",
            )
            self.assertEqual(
                external_fidelity_matrix_mjcf_run["schema"],
                TWO_CELL_EXTERNAL_FIDELITY_MJCF_RUN_SCHEMA,
            )
            self.assertEqual(
                response_atlas_priority_mjcf_index["schema"],
                "rad-sim.two-cell-response-atlas-priority-mjcf-index.v1",
            )
            self.assertEqual(
                response_atlas_priority_mjcf_run["schema"],
                TWO_CELL_EXTERNAL_FIDELITY_MJCF_RUN_SCHEMA,
            )
            self.assertEqual(quasistatic["schema"], TWO_CELL_QUASISTATIC_SCHEMA)
            self.assertEqual(quasistatic_report["schema"], TWO_CELL_QUASISTATIC_SWEEP_SCHEMA)
            self.assertEqual(connector_contact["schema"], TWO_CELL_CONNECTOR_CONTACT_SCHEMA)
            self.assertEqual(connector_contact_sweep["schema"], TWO_CELL_CONNECTOR_CONTACT_SWEEP_SCHEMA)
            self.assertEqual(connector_measurement_template["schema"], TWO_CELL_CONNECTOR_MEASUREMENT_TEMPLATE_SCHEMA)
            self.assertEqual(external_template["schema"], TWO_CELL_EXTERNAL_TEMPLATE_SCHEMA)
            self.assertEqual(physics_validation["schema"], TWO_CELL_PHYSICS_VALIDATION_SCHEMA)
            self.assertIn("<mujoco", files["mjcf_xml"].read_text(encoding="utf-8"))
            self.assertEqual(packet["schema"], TWO_CELL_PACKET_SCHEMA)
            self.assertEqual(packet["cadArchiveAudit"]["schema"], CAD_RAD_CELL_ARCHIVE_AUDIT_SCHEMA)
            self.assertEqual(packet["cadReferenceProfile"]["schema"], CAD_RAD_CELL_REFERENCE_PROFILE_SCHEMA)
            self.assertEqual(packet["oneCellCadExportAudit"]["schema"], ONE_CELL_CAD_EXPORT_AUDIT_SCHEMA)
            self.assertEqual(packet["contactPhaseMap"]["schema"], TWO_CELL_CONTACT_PHASE_MAP_SCHEMA)
            self.assertEqual(packet["radiusBacklashPhaseDiagram"]["schema"], TWO_CELL_RADIUS_BACKLASH_PHASE_DIAGRAM_SCHEMA)
            self.assertEqual(
                packet["radiusBacklashTransitionReport"]["schema"],
                TWO_CELL_RADIUS_BACKLASH_TRANSITION_REPORT_SCHEMA,
            )
            self.assertEqual(packet["exactContactHandoffPlan"]["schema"], TWO_CELL_EXACT_CONTACT_HANDOFF_PLAN_SCHEMA)
            self.assertEqual(len(sweep_rows), 5)
            self.assertEqual(len(segmented_cad_rows), segmented_cad_readiness["summary"]["requiredTaskCount"])
            self.assertEqual(len(segmented_template_rows), segmented_cad_templates["templateCount"])
            self.assertEqual(len(segmented_validation_rows), 15)
            self.assertTrue((files["segmented_cad_template_dir"] / "joint_axes.json").exists())
            self.assertTrue((files["segmented_cad_template_dir"] / "bench_coordinate_truth.csv").exists())
            self.assertEqual(len(actuation_rows), actuation_report["summary"]["rowCount"])
            self.assertEqual(len(quasistatic_rows), quasistatic_report["summary"]["rowCount"])
            self.assertEqual(len(connector_contact_rows), connector_contact["summary"]["connectorCount"])
            self.assertEqual(len(connector_contact_sweep_rows), connector_contact_sweep["summary"]["rowCount"])
            self.assertEqual(len(connector_measurement_rows), connector_measurement_template["summary"]["connectorRowCount"])
            self.assertEqual(len(physical_suite_rows), physical_suite["summary"]["caseCount"])
            self.assertEqual(len(physical_fidelity_rows), physical_fidelity_matrix["summary"]["rowCount"])
            self.assertGreaterEqual(len(physical_response_rows), 5)
            self.assertIn("caseId", physical_response_rows[0])
            self.assertEqual(len(contact_phase_rows), contact_phase_map["summary"]["rowCount"])
            self.assertIn("phase", contact_phase_rows[0])
            self.assertEqual(
                len(radius_backlash_phase_rows),
                radius_backlash_phase_diagram["summary"]["rowCount"],
            )
            self.assertIn("holeRadius", radius_backlash_phase_rows[0])
            self.assertEqual(
                len(radius_backlash_transition_rows),
                radius_backlash_transition_report["summary"]["measurementCaseCount"],
            )
            self.assertIn("observedPhase", radius_backlash_transition_rows[0])
            self.assertGreaterEqual(len(cad_contact_decomposition_rows), 16 * 12)
            self.assertIn("convex_ring_sector_or_capsule_wall", {row["primitive"] for row in cad_contact_decomposition_rows})
            self.assertEqual(packet["cadContactDecomposition"]["schema"], TWO_CELL_CAD_CONTACT_DECOMPOSITION_SCHEMA)
            self.assertEqual(len(exact_contact_handoff_rows), exact_contact_handoff_plan["summary"]["caseCount"])
            self.assertIn("requiredCadAssets", exact_contact_handoff_rows[0])
            self.assertFalse(exact_contact_handoff_plan["summary"]["canRunExactContact"])
            self.assertEqual(len(fidelity_matrix_measurement_rows), 756)
            self.assertEqual(
                len(external_fidelity_manifest_rows),
                external_fidelity_matrix_manifest["summary"]["caseCount"],
            )
            self.assertEqual(
                external_fidelity_matrix_mjcf_index["caseCount"],
                external_fidelity_matrix_manifest["summary"]["caseCount"],
            )
            self.assertIn("xyMmPerModelUnit", external_fidelity_matrix_mjcf_index["files"][0])
            self.assertIn("pinHoleClearanceMm", external_fidelity_matrix_mjcf_index["files"][0])
            self.assertEqual(
                external_fidelity_matrix_mjcf_run["summary"]["caseCount"],
                external_fidelity_matrix_manifest["summary"]["caseCount"],
            )
            self.assertEqual(
                external_fidelity_matrix_mjcf_run["summary"]["expectedConnectorMeasurementRowCount"],
                3 * external_fidelity_matrix_manifest["summary"]["caseCount"],
            )
            self.assertFalse(external_fidelity_matrix_mjcf_run["summary"]["externalFidelityMjcfRunComplete"])
            self.assertIn("externalFidelityMjcfRun", external_fidelity_matrix_mjcf_run["summary"]["missingEvidence"])
            self.assertIn("mujocoConnectorMeasurements", external_fidelity_matrix_mjcf_run["summary"]["missingEvidence"])
            self.assertEqual(
                len(list(files["external_fidelity_matrix_mjcf_dir"].glob("FM_*.xml"))),
                external_fidelity_matrix_manifest["summary"]["caseCount"],
            )
            self.assertEqual(
                response_atlas_priority_mjcf_index["caseCount"],
                len(physical_response_atlas["benchPriority"]["firstCasesToMeasure"]),
            )
            self.assertEqual(
                response_atlas_priority_mjcf_run["summary"]["caseCount"],
                response_atlas_priority_mjcf_index["caseCount"],
            )
            self.assertEqual(
                len(list(files["response_atlas_priority_mjcf_dir"].glob("FM_*.xml"))),
                response_atlas_priority_mjcf_index["caseCount"],
            )
            self.assertEqual(
                {case["caseId"] for case in response_atlas_priority_mjcf_index["files"]},
                {case["caseId"] for case in physical_response_atlas["benchPriority"]["firstCasesToMeasure"]},
            )
            fm0_xml = (files["external_fidelity_matrix_mjcf_dir"] / "FM_000.xml").read_text(encoding="utf-8")
            self.assertIn("<mujoco", fm0_xml)
            self.assertIn("left_cell_ne_marker", fm0_xml)
            self.assertEqual(len(external_rows), external_template["summary"]["rowCount"])
            self.assertEqual(len(physics_validation_rows), physics_validation["summary"]["testCount"])
            self.assertEqual(len(template_rows), 5)

    def test_two_cell_measurement_comparison_is_zero_for_model_generated_rows(self):
        packet = two_cell_physical_test_packet(LatticeConfig(rows=1, cols=2))
        rows = []
        for case in packet["cases"]:
            bench = case["bench"]
            rows.append(
                {
                    "caseId": case["caseId"],
                    "leftX": bench["cells"][0]["center"]["x"],
                    "leftY": bench["cells"][0]["center"]["y"],
                    "leftZ": bench["cells"][0]["center"]["z"],
                    "rightX": bench["cells"][1]["center"]["x"],
                    "rightY": bench["cells"][1]["center"]["y"],
                    "rightZ": bench["cells"][1]["center"]["z"],
                    "rightAlpha": bench["cells"][1]["alpha"],
                    "rightTheta": bench["cells"][1]["theta"],
                    "contactModeObserved": bench["connector"]["contactMode"],
                }
            )
        comparison = compare_two_cell_measurements(two_cell_measurements_from_json(json.dumps(rows)))
        self.assertEqual(comparison["schema"], TWO_CELL_MEASUREMENT_COMPARISON_SCHEMA)
        self.assertEqual(comparison["sampleCount"], 5)
        self.assertAlmostEqual(comparison["rmsError"], 0.0)
        self.assertAlmostEqual(comparison["contactModeAccuracy"], 1.0)
        self.assertFalse(comparison["acceptance"]["requiresMoreData"])

    def test_two_cell_measurement_comparison_detects_biased_data(self):
        csv_text = "\n".join(
            [
                "caseId,rightZ,rightAlpha,contactModeObserved",
                "free_contract_lift,0.99,0.50,wrong-contact",
                "right_state_locked,0.25,1.00,blocked-by-lock",
            ]
        )
        measurements = two_cell_measurements_from_csv(csv_text)
        comparison = compare_two_cell_measurements(measurements)
        self.assertGreater(comparison["rmsError"], 0.1)
        self.assertLess(comparison["contactModeAccuracy"], 1.0)
        self.assertTrue(comparison["acceptance"]["requiresMoreData"])
        payload = json.loads(export_two_cell_measurement_comparison_json(measurements))
        self.assertEqual(payload["schema"], TWO_CELL_MEASUREMENT_COMPARISON_SCHEMA)

    def test_two_cell_parameter_calibration_selects_best_candidate(self):
        true_config = LatticeConfig(rows=1, cols=2, backlash=0.06, pin_radius=0.18, hole_radius=0.24)
        protocol_controls = TwoCellBenchControls(alpha_command=-0.5, z_command=0.45)
        packet = two_cell_physical_test_packet(
            true_config,
            protocol_controls,
        )
        measurements = []
        for case in packet["cases"]:
            bench = case["bench"]
            measurements.append(
                TwoCellMeasurement(
                    case_id=case["caseId"],
                    left_z=bench["cells"][0]["center"]["z"],
                    right_x=bench["cells"][1]["center"]["x"],
                    right_z=bench["cells"][1]["center"]["z"],
                    right_alpha=bench["cells"][1]["alpha"],
                    right_theta=bench["cells"][1]["theta"],
                    contact_mode_observed=bench["connector"]["contactMode"],
                )
            )
        calibration = calibrate_two_cell_parameters(
            measurements,
            LatticeConfig(rows=1, cols=2, backlash=0.1, pin_radius=0.18, hole_radius=0.225),
            protocol_controls,
            backlash_candidates=[0.04, 0.06, 0.1],
            hole_radius_candidates=[0.225, 0.24, 0.30],
            z_coupling_gain_candidates=[0.32],
        )
        self.assertEqual(calibration["schema"], TWO_CELL_PARAMETER_CALIBRATION_SCHEMA)
        self.assertEqual(calibration["best"]["backlash"], 0.06)
        self.assertEqual(calibration["best"]["holeRadius"], 0.24)
        calibration_json = json.loads(export_two_cell_parameter_calibration_json(measurements))
        self.assertEqual(calibration_json["schema"], TWO_CELL_PARAMETER_CALIBRATION_SCHEMA)

    def test_two_cell_measurement_comparison_writer_creates_reports(self):
        bench = simulate_two_cell_bench()
        csv_text = "\n".join(
            [
                "caseId,rightZ,rightAlpha,rightTheta,contactModeObserved",
                f"free_contract_lift,{bench['cells'][1]['center']['z']},{bench['cells'][1]['alpha']},{bench['cells'][1]['theta']},{bench['connector']['contactMode']}",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "measurements.csv"
            input_path.write_text(csv_text, encoding="utf-8")
            files = write_two_cell_measurement_comparison_artifacts(input_path, Path(tmp) / "reports")
            comparison = json.loads(files["comparison"].read_text(encoding="utf-8"))
            calibration = json.loads(files["calibration"].read_text(encoding="utf-8"))
            self.assertEqual(comparison["schema"], TWO_CELL_MEASUREMENT_COMPARISON_SCHEMA)
            self.assertEqual(calibration["schema"], TWO_CELL_PARAMETER_CALIBRATION_SCHEMA)

    def test_paper_rad_reference_contains_extracted_prototype_values(self):
        self.assertEqual(PAPER_RAD_REFERENCE.concentric_parts, 2)
        self.assertEqual(PAPER_RAD_REFERENCE.joints_per_part, 4)
        self.assertAlmostEqual(PAPER_RAD_REFERENCE.normalized_backlash, 0.1)
        self.assertAlmostEqual(PAPER_RAD_REFERENCE.poisson_ratio, -0.4)
        self.assertAlmostEqual(PAPER_RAD_REFERENCE.side_length_mm, 35.0)
        self.assertAlmostEqual(PAPER_RAD_REFERENCE.fabrication_hole_tolerance_mm, 0.1)
        self.assertAlmostEqual(PAPER_RAD_REFERENCE.reference_backlash_mm, 3.5)

    def test_model_provenance_classifies_evidence_and_assumptions(self):
        items = model_provenance()
        ids = {item.id for item in items}
        self.assertIn("backlash_dead_zone", ids)
        self.assertIn("vertical_residual_coupling", ids)
        self.assertGreaterEqual(len(provenance_by_status("paper-supported")), 3)
        summary = provenance_summary()
        self.assertGreaterEqual(summary["paper-supported"], 3)
        self.assertGreaterEqual(summary["implementation-assumption"], 1)
        self.assertGreaterEqual(summary["simulator-diagnostic"], 1)
        self.assertGreaterEqual(summary["calibration-gap"], 1)
        backlash = next(item for item in items if item.id == "backlash_dead_zone")
        self.assertEqual(backlash.status, "paper-supported")
        self.assertTrue(any("Eq. 2" in ref for ref in backlash.page_refs))

    def test_paper_rad_calibration_converts_model_units_to_mm(self):
        config = LatticeConfig(
            backlash=0.1,
            cell_size=1.0,
            pin_radius=0.18,
            hole_radius=0.225,
        )
        calibration = calibrate_paper_rad_config(config)
        self.assertAlmostEqual(calibration.mm_per_model_unit, 35.0)
        self.assertAlmostEqual(calibration.configured_backlash_mm, 3.5)
        self.assertAlmostEqual(calibration.reference_backlash_mm, 3.5)
        self.assertAlmostEqual(calibration.pin_radius_mm, 6.3)
        self.assertAlmostEqual(calibration.hole_radius_mm, 7.875)
        self.assertAlmostEqual(calibration.pin_hole_clearance_mm, 1.575)
        self.assertAlmostEqual(calibration.fabrication_hole_tolerance_model, 0.1 / 35.0)
        self.assertAlmostEqual(calibration.model_length_to_mm(0.25), 8.75)
        self.assertAlmostEqual(calibration.mm_to_model_length(8.75), 0.25)

    def test_paper_rad_calibration_respects_nonunit_model_cell_size(self):
        config = LatticeConfig(cell_size=2.0, backlash=0.1, pin_radius=0.2, hole_radius=0.3)
        calibration = calibrate_paper_rad_config(config)
        self.assertAlmostEqual(calibration.mm_per_model_unit, 17.5)
        self.assertAlmostEqual(calibration.configured_backlash_mm, 3.5)
        self.assertAlmostEqual(calibration.pin_hole_clearance_mm, 1.75)
        self.assertAlmostEqual(calibration.fabrication_hole_tolerance_model, 0.1 * 2.0 / 35.0)

    def test_physical_unit_scale_metadata_links_calibration_to_lean_target(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            cell_size=2.0,
            backlash=0.1,
            pin_radius=0.2,
            hole_radius=0.3,
        )
        metadata = physical_unit_scale_metadata(config)
        self.assertEqual(metadata["schema"], "rad-sim.physical-unit-scale-metadata.v1")
        self.assertEqual(
            metadata["formalizationTarget"]["id"],
            "measurement_unit_scale_invariants",
        )
        self.assertIn(
            "measurementUnitScaleResidualInt_zero_when_equal",
            metadata["formalizationTarget"]["leanTheorems"],
        )
        self.assertEqual(
            metadata["scale"]["modelUnitToMillimeter"]["numerator"],
            35,
        )
        self.assertEqual(
            metadata["scale"]["modelUnitToMillimeter"]["denominator"],
            2,
        )
        self.assertAlmostEqual(
            metadata["configuredPhysicalQuantities"]["pinHoleClearanceMm"]["value"],
            1.75,
        )
        self.assertEqual(
            metadata["zeroAndResidualExamples"]["zeroModelLengthToMm"]["scaledNumerator"],
            0,
        )
        exported = json.loads(export_physical_unit_scale_metadata_json(config))
        self.assertEqual(exported["schema"], metadata["schema"])

    def test_hardware_profile_tracks_measured_dimension_coverage(self):
        profile = RADHardwareProfile(
            name="bench-measurement-v1",
            source="caliper",
            pin_radius_mm=6.3,
            hole_radius_mm=7.875,
            plate_thickness_mm=2.1,
        )
        self.assertEqual(
            profile.measured_fields,
            ("pin_radius_mm", "hole_radius_mm", "plate_thickness_mm"),
        )
        self.assertEqual(len(profile.missing_fields), len(HARDWARE_PROFILE_DIMENSIONS) - 3)
        self.assertAlmostEqual(profile.coverage_ratio, 3 / len(HARDWARE_PROFILE_DIMENSIONS))
        self.assertAlmostEqual(profile.pin_hole_clearance_mm, 1.575)
        reference = profile.to_reference()
        self.assertAlmostEqual(reference.side_length_mm, 35.0)
        self.assertAlmostEqual(reference.fabrication_hole_tolerance_mm, 0.1)

    def test_hardware_profile_json_roundtrip_preserves_measured_dimensions(self):
        profile = RADHardwareProfile(
            name="bench-measured-v2",
            source="digital calipers",
            side_length_mm=40.0,
            fabrication_hole_tolerance_mm=0.08,
            backlash_mm=4.0,
            pin_radius_mm=5.0,
            hole_radius_mm=5.8,
            plate_thickness_mm=2.4,
            joint_stack_height_mm=4.5,
            boss_radius_mm=2.0,
            notes="prototype A",
        )
        payload = profile.to_dict()
        self.assertEqual(payload["schema"], HARDWARE_PROFILE_SCHEMA)
        self.assertEqual(payload["dimensionsMm"]["pinRadiusMm"], 5.0)
        self.assertAlmostEqual(payload["derived"]["pinHoleClearanceMm"], 0.8)
        restored = hardware_profile_from_json(export_hardware_profile_json(profile))
        self.assertEqual(restored.name, profile.name)
        self.assertEqual(restored.source, profile.source)
        self.assertAlmostEqual(restored.side_length_mm, 40.0)
        self.assertAlmostEqual(restored.pin_hole_clearance_mm, 0.8)
        self.assertEqual(restored.measured_fields, HARDWARE_PROFILE_DIMENSIONS)
        bom_restored = hardware_profile_from_json(
            "\ufeff" + export_hardware_profile_json(profile)
        )
        self.assertEqual(bom_restored.name, profile.name)

    def test_physical_unit_scale_metadata_accepts_measured_hardware_profile(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            cell_size=2.0,
            backlash=0.05,
            pin_radius=0.1,
            hole_radius=0.12,
        )
        profile = RADHardwareProfile(
            name="bench-measured-v2",
            side_length_mm=40.0,
            backlash_mm=4.0,
            pin_radius_mm=5.0,
            hole_radius_mm=5.8,
            plate_thickness_mm=2.4,
            joint_stack_height_mm=4.5,
            boss_radius_mm=2.0,
        )
        metadata = physical_unit_scale_metadata(config, hardware_profile=profile)
        self.assertEqual(metadata["hardwareProfile"]["schema"], HARDWARE_PROFILE_SCHEMA)
        self.assertEqual(metadata["hardwareProfile"]["name"], "bench-measured-v2")
        self.assertTrue(metadata["profileApplication"]["applied"])
        self.assertAlmostEqual(metadata["grid"]["backlash"], 0.1)
        self.assertAlmostEqual(metadata["grid"]["pinRadius"], 0.25)
        self.assertAlmostEqual(metadata["grid"]["holeRadius"], 0.29)
        self.assertAlmostEqual(
            metadata["configuredPhysicalQuantities"]["pinHoleClearanceMm"]["value"],
            0.8,
        )
        self.assertEqual(metadata["calibrationReadiness"]["level"], "mesh-calibrated")
        self.assertIn(
            "hardwareProfileCoverageMissing_zero_when_complete",
            metadata["formalizationTarget"]["hardwareProfileLeanTheorems"],
        )

    def test_calibration_readiness_reports_partial_and_solver_gaps(self):
        profile = RADHardwareProfile(
            pin_radius_mm=6.3,
            hole_radius_mm=7.875,
            plate_thickness_mm=2.1,
        )
        readiness = calibration_readiness(profile)
        self.assertEqual(readiness.level, "partial-measured")
        self.assertFalse(readiness.visual_ready)
        self.assertFalse(readiness.mesh_ready)
        self.assertFalse(readiness.solver_ready)
        self.assertIn("joint_stack_height_mm", readiness.visual_missing_fields)
        self.assertEqual(readiness.solver_gaps, CALIBRATION_SOLVER_GAPS)
        self.assertIn("solver still", readiness.summary)

    def test_calibration_readiness_reports_mesh_calibrated_geometry(self):
        profile = RADHardwareProfile(
            pin_radius_mm=6.3,
            hole_radius_mm=7.875,
            plate_thickness_mm=2.1,
            joint_stack_height_mm=4.0,
            boss_radius_mm=2.2,
        )
        readiness = calibration_readiness(profile)
        self.assertEqual(readiness.level, "mesh-calibrated")
        self.assertTrue(readiness.visual_ready)
        self.assertTrue(readiness.mesh_ready)
        self.assertFalse(readiness.solver_ready)
        self.assertEqual(readiness.mesh_missing_fields, ())
        self.assertIn("solver still", readiness.summary)

    def test_calibration_measurement_plan_lists_missing_geometry_before_solver(self):
        profile = RADHardwareProfile(
            pin_radius_mm=6.3,
            hole_radius_mm=7.875,
            plate_thickness_mm=2.1,
        )
        plan = calibration_measurement_plan(profile)
        missing = [task for task in plan if task.status == "missing"]
        done = [task for task in plan if task.status == "done"]
        self.assertEqual(missing[0].id, "geometry_joint_stack_height_mm")
        self.assertEqual(missing[1].id, "geometry_boss_radius_mm")
        self.assertEqual(missing[0].category, "geometry")
        self.assertEqual(missing[0].evidence_field, "joint_stack_height_mm")
        self.assertIn("Required", missing[0].notes)
        self.assertTrue(any(task.id == "solver_response_data" for task in missing))
        self.assertEqual(len([task for task in plan if task.category == "solver"]), len(CALIBRATION_SOLVER_TASKS))
        self.assertTrue(all(task.status == "missing" for task in plan if task.category == "solver"))
        self.assertEqual(done[0].id, "geometry_pin_radius_mm")

    def test_calibration_measurement_plan_keeps_solver_tasks_missing_after_mesh_calibration(self):
        profile = RADHardwareProfile(
            pin_radius_mm=6.3,
            hole_radius_mm=7.875,
            plate_thickness_mm=2.1,
            joint_stack_height_mm=4.0,
            boss_radius_mm=2.2,
        )
        plan = calibration_measurement_plan(profile)
        geometry = [task for task in plan if task.category == "geometry"]
        solver = [task for task in plan if task.category == "solver"]
        self.assertTrue(all(task.status == "done" for task in geometry))
        self.assertTrue(all(task.status == "missing" for task in solver))
        self.assertEqual(solver[0].id, "solver_axial_hinge_stiffness")

    def test_hardware_profile_updates_config_from_measured_radii(self):
        config = LatticeConfig(cell_size=1.0, pin_radius=0.1, hole_radius=0.12, backlash=0.05)
        profile = RADHardwareProfile(
            side_length_mm=35.0,
            backlash_mm=3.5,
            pin_radius_mm=6.3,
            hole_radius_mm=7.875,
        )
        updated = config_with_hardware_profile(config, profile)
        self.assertAlmostEqual(updated.backlash, 0.1)
        self.assertAlmostEqual(updated.pin_radius, 0.18)
        self.assertAlmostEqual(updated.hole_radius, 0.225)
        self.assertAlmostEqual(updated.pin_hole_clearance, 0.045)

    def test_partial_hardware_profile_keeps_config_radii_valid(self):
        config = LatticeConfig(cell_size=1.0, pin_radius=0.1, hole_radius=0.12)
        profile = RADHardwareProfile(pin_radius_mm=7.0)
        updated = config_with_hardware_profile(config, profile)
        self.assertAlmostEqual(updated.pin_radius, 0.2)
        self.assertGreaterEqual(updated.hole_radius, updated.pin_radius)

    def test_hardware_profile_from_config_marks_current_values_as_estimates(self):
        config = LatticeConfig(cell_size=2.0, pin_radius=0.2, hole_radius=0.3, backlash=0.1)
        profile = hardware_profile_from_config(config)
        self.assertEqual(profile.name, "current-config-estimate")
        self.assertEqual(profile.source, "normalized simulator controls")
        self.assertAlmostEqual(profile.pin_radius_mm, 3.5)
        self.assertAlmostEqual(profile.hole_radius_mm, 5.25)
        self.assertAlmostEqual(profile.pin_hole_clearance_mm, 1.75)
        self.assertIn("plate_thickness_mm", profile.missing_fields)

    def test_hardware_profile_rejects_inverted_pin_hole_radii(self):
        with self.assertRaises(ValueError):
            RADHardwareProfile(pin_radius_mm=2.0, hole_radius_mm=1.0)

    def test_paper_rad_cell_has_two_four_joint_parts(self):
        config = LatticeConfig(pin_radius=0.12, hole_radius=0.19)
        cell = build_paper_rad_cell_geometry(config, center=(0.2, -0.1), z=0.3)
        self.assertEqual(cell.outer_part.shape, (4, 3))
        self.assertEqual(cell.inner_part.shape, (4, 3))
        self.assertEqual(len(cell.outer_joints), 4)
        self.assertEqual(len(cell.inner_joints), 4)
        self.assertEqual(cell.joint_count, 8)
        self.assertAlmostEqual(cell.center[2], 0.3)
        self.assertAlmostEqual(cell.backlash_gap, config.backlash * config.cell_size)
        self.assertAlmostEqual(cell.vertical_free_play, config.pin_hole_clearance)
        self.assertAlmostEqual(cell.backlash_gap_mm, 3.5)
        self.assertAlmostEqual(cell.vertical_free_play_mm, 2.45)
        self.assertTrue(all(joint.clearance == config.pin_hole_clearance for joint in cell.outer_joints))

    def test_paper_rad_inner_part_rotates_with_alpha(self):
        config = LatticeConfig(backlash=0.1)
        contracted = build_paper_rad_cell_geometry(config, alpha=0.8)
        expanded = build_paper_rad_cell_geometry(config, alpha=1.2)
        self.assertAlmostEqual(contracted.theta_degrees, alpha_to_theta(0.8))
        self.assertAlmostEqual(expanded.theta_degrees, alpha_to_theta(1.2))
        self.assertFalse(np.allclose(contracted.inner_part, expanded.inner_part))
        contracted_edges = np.linalg.norm(
            np.roll(contracted.inner_part[:, :2], -1, axis=0) - contracted.inner_part[:, :2],
            axis=1,
        )
        np.testing.assert_allclose(contracted_edges, contracted_edges[0])

    def test_paper_rad_cell_exposes_lock_and_actuator_interfaces(self):
        config = LatticeConfig(pin_radius=0.18, hole_radius=0.225)
        cell = build_paper_rad_cell_geometry(config)
        self.assertEqual(cell.lock_sites.shape, (4, 3))
        self.assertEqual(cell.alpha_actuator_axis.shape, (2, 3))
        self.assertEqual(cell.z_actuator_axis.shape, (2, 3))
        self.assertAlmostEqual(
            np.linalg.norm(cell.z_actuator_axis[1] - cell.z_actuator_axis[0]),
            config.pin_hole_clearance,
        )

    def test_paper_rad_lattice_geometry_builds_cells_and_connectors(self):
        config = LatticeConfig(rows=2, cols=3)
        lattice = build_paper_rad_lattice_geometry(config)
        self.assertEqual(lattice.shape, (2, 3))
        self.assertEqual(lattice.cell_count, 6)
        self.assertEqual(lattice.connector_count, 7)
        self.assertEqual(lattice.centers.shape, (2, 3, 3))
        self.assertEqual(lattice.cell_at(1, 2).index, (1, 2))
        self.assertTrue(all(connector.length > 0.0 for connector in lattice.connectors))

    def test_paper_rad_lattice_geometry_preserves_state_flags(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.2)
        state = LatticeState.uniform(config)
        state.locked_mask[1, 1] = True
        state.actuator_grid[0, 1] = -0.25
        state.z_actuator_grid[2, 1] = 0.3
        lattice = build_paper_rad_lattice_geometry(config, state)
        self.assertTrue(lattice.cell_at(1, 1).locked)
        self.assertEqual(lattice.active_actuator_count, 2)
        self.assertAlmostEqual(lattice.cell_at(0, 1).command_alpha, -0.25)
        self.assertAlmostEqual(lattice.cell_at(2, 1).command_z, 0.3)
        np.testing.assert_array_equal(lattice.locked_mask, state.locked_mask)

    def test_paper_rad_lattice_geometry_accepts_physical_result(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.0)
        state = LatticeState.uniform(config)
        state.z_actuator_grid[1, 1] = 0.3
        physical = solve_spring_hinge_3d(
            config,
            state,
            LoadCase(lock_stiffness=500.0, maxiter=400),
        )
        lattice = build_paper_rad_lattice_geometry(config, physical)
        np.testing.assert_allclose(lattice.centers, physical.deformed_centers_3d)
        np.testing.assert_allclose(lattice.height, physical.deformed_centers_3d[..., 2])

    def test_paper_rad_lattice_connectors_track_neighbor_jumps(self):
        config = LatticeConfig(rows=1, cols=2, z_coupling_gain=0.0, backlash=0.0)
        state = LatticeState.uniform(config)
        state.z_actuator_grid[0, 1] = 0.2
        lattice = build_paper_rad_lattice_geometry(config, state)
        connector = lattice.connectors[0]
        self.assertEqual(connector.first, (0, 0))
        self.assertEqual(connector.second, (0, 1))
        self.assertEqual(connector.axis, "x")
        self.assertAlmostEqual(connector.height_delta, 0.2)
        self.assertAlmostEqual(connector.backlash_gap, 0.0)

    def test_paper_rad_lattice_mesh_contains_plates_pins_and_connectors(self):
        config = LatticeConfig(rows=2, cols=2)
        mesh = build_paper_rad_lattice_mesh(config, pin_segments=8)
        counts = mesh.kind_counts()
        self.assertEqual(counts["outer_plate"], 4)
        self.assertEqual(counts["inner_plate"], 4)
        self.assertEqual(counts["pin"], 32)
        self.assertEqual(counts["connector"], 4)
        self.assertGreater(mesh.vertex_count, 0)
        self.assertGreater(mesh.face_count, 0)

    def test_paper_rad_lattice_mesh_tracks_vertical_deformation(self):
        config = LatticeConfig(rows=1, cols=2, z_coupling_gain=0.0)
        state = LatticeState.uniform(config)
        state.z_actuator_grid[0, 1] = 0.3
        mesh = build_paper_rad_lattice_mesh(config, state, include_pins=False)
        lower, upper = mesh.bounds
        self.assertGreater(upper[2] - lower[2], 0.3)
        self.assertEqual(mesh.kind_counts()["connector"], 1)

    def test_paper_rad_lattice_mesh_uses_hardware_profile_dimensions(self):
        config = LatticeConfig(rows=1, cols=1, cell_size=1.0)
        profile = RADHardwareProfile(
            pin_radius_mm=7.0,
            hole_radius_mm=8.0,
            plate_thickness_mm=3.5,
            joint_stack_height_mm=4.2,
            boss_radius_mm=1.4,
        )
        mesh = build_paper_rad_lattice_mesh(
            config,
            hardware_profile=profile,
            pin_segments=8,
            include_connectors=False,
        )
        outer_plate = next(component for component in mesh.components if component.kind == "outer_plate")
        pin = next(component for component in mesh.components if component.kind == "pin")
        self.assertAlmostEqual(np.ptp(outer_plate.vertices[:, 2]), 0.1)
        self.assertAlmostEqual(np.ptp(pin.vertices[:, 2]), 0.12)
        center = pin.vertices[-2]
        radial_distance = np.linalg.norm(pin.vertices[0, :2] - center[:2])
        self.assertAlmostEqual(radial_distance, 0.2)

    def test_paper_rad_mesh_obj_export_has_valid_indices(self):
        config = LatticeConfig(rows=1, cols=1)
        mesh = build_paper_rad_lattice_mesh(
            config,
            include_pins=False,
            include_connectors=False,
        )
        obj = export_paper_rad_mesh_obj(mesh)
        self.assertIn("o cell_0_0_outer_plate", obj)
        self.assertIn("o cell_0_0_inner_plate", obj)
        vertices = list(iter_obj_vertices(obj))
        faces = [line for line in obj.splitlines() if line.startswith("f ")]
        self.assertEqual(len(vertices), mesh.vertex_count)
        self.assertEqual(len(faces), mesh.face_count)
        max_face_index = max(int(index) for face in faces for index in face.split()[1:])
        self.assertLessEqual(max_face_index, mesh.vertex_count)

    def test_paper_rad_mesh_rejects_invalid_dimensions(self):
        with self.assertRaises(ValueError):
            build_paper_rad_lattice_mesh(
                LatticeConfig(rows=1, cols=1),
                plate_thickness=0.0,
            )

    def test_dead_zone_operator_has_no_neighbor_transmission_inside_gap(self):
        op = DeadZonePropagationOperator("test", dead_zone=0.1, gain=0.5)
        result = op.propagate(3, 3, [(1, 1, 0.05)])
        self.assertAlmostEqual(result.field[1, 1], 0.05)
        self.assertAlmostEqual(result.field[1, 2], 0.0)
        self.assertTrue(np.isinf(result.die_off[1, 2]))

    def test_dead_zone_operator_is_monotone_outside_gap(self):
        op = DeadZonePropagationOperator("test", dead_zone=0.1, gain=0.5)
        self.assertAlmostEqual(op.transmit(0.2), 0.05)
        self.assertAlmostEqual(op.transmit(0.3), 0.1)
        self.assertGreater(op.transmit(0.3), op.transmit(0.2))
        self.assertLess(op.transmit(-0.3), op.transmit(-0.2))

    def test_lock_projection_keeps_locked_cell_invariant(self):
        reference = np.array([[1.0, 1.0], [1.0, 1.0]])
        proposed = np.array([[0.7, 1.2], [1.4, 0.8]])
        locked = np.array([[True, False], [False, True]])
        projected = lock_projection(reference, proposed, locked)
        np.testing.assert_allclose(projected, [[1.0, 1.2], [1.4, 1.0]])

    def test_operator_stack_matches_kinematic_fields(self):
        config = LatticeConfig(rows=4, cols=4, backlash=0.04, z_coupling_gain=0.35)
        state = LatticeState.uniform(config)
        state.actuator_grid[1, 1] = -0.3
        state.z_actuator_grid[1, 1] = 0.35
        state.locked_mask[0, 0] = True
        operators = evaluate_programmable_operators(config, state)
        result = simulate_kinematic(config, state)
        np.testing.assert_allclose(operators["alpha"], result.alpha)
        np.testing.assert_allclose(operators["height"], result.metadata["height"])
        np.testing.assert_allclose(
            operators["z_residual"], result.metadata["z_residual"]
        )

    def test_vertical_clearance_operator_uses_pin_hole_gap(self):
        config = LatticeConfig(pin_radius=0.12, hole_radius=0.19)
        op = vertical_clearance_operator(config)
        self.assertAlmostEqual(op.dead_zone, 0.07)

    def test_contact_state_abstraction_reports_pin_hole_modes_and_penalty(self):
        config = LatticeConfig(rows=1, cols=3, pin_radius=0.10, hole_radius=0.12)
        state = LatticeState.uniform(config)
        state.z_actuator_grid[0, 0] = 0.01
        state.z_actuator_grid[0, 1] = 0.08
        state.removed_mask[0, 2] = True
        report = contact_state_abstraction_report(
            config,
            state,
            contact_stiffness=2.0,
        )
        self.assertEqual(report["schema"], "rad-sim.contact-state-abstraction.v1")
        self.assertTrue(report["summary"]["contactStateAbstractionReady"])
        self.assertEqual(report["topology"]["activeBodyCount"], 2)
        self.assertEqual(report["topology"]["removedBodyCount"], 1)
        self.assertEqual(report["topology"]["pinCount"], 2)
        self.assertEqual(report["topology"]["holeCount"], 2)
        self.assertEqual(report["topology"]["penaltyTermCount"], 2)
        modes = {(cell["row"], cell["col"]): cell["mode"] for cell in report["cells"]}
        self.assertEqual(modes[(0, 0)], "free-clearance")
        self.assertEqual(modes[(0, 1)], "engaged")
        self.assertEqual(modes[(0, 2)], "removed")
        engaged = next(cell for cell in report["cells"] if cell["mode"] == "engaged")
        self.assertAlmostEqual(engaged["penetration"], 0.06)
        self.assertAlmostEqual(engaged["contactPenalty"], 0.0036)
        self.assertEqual(
            report["formalization"]["targetId"],
            "contact_state_abstraction_gate",
        )
        csv_text = export_contact_state_abstraction_csv(report)
        self.assertIn("contact_penalty", csv_text)
        self.assertIn("engaged", csv_text)
        exported = json.loads(
            export_contact_state_abstraction_json(
                config,
                state,
                contact_stiffness=2.0,
            )
        )
        self.assertEqual(exported["summary"], report["summary"])
        not_ready = contact_state_abstraction_report(
            config,
            state,
            contact_stiffness=0.0,
        )
        self.assertFalse(not_ready["summary"]["contactStateAbstractionReady"])
        self.assertIn("positiveContactStiffness", not_ready["summary"]["missingEvidence"])

    def test_contact_graph_consistency_deletes_removed_edges_and_tracks_support(self):
        config = LatticeConfig(rows=2, cols=2, pin_radius=0.10, hole_radius=0.12)
        state = LatticeState.uniform(config)
        state.removed_mask[0, 1] = True
        state.z_actuator_grid[1, 1] = 0.08
        state.actuator_grid[0, 0] = 0.2
        contact_report = contact_state_abstraction_report(config, state)
        group_support = ((0, 1), (1, 1))
        report = contact_graph_consistency_report(
            config,
            state,
            contact_report=contact_report,
            group_support=group_support,
        )
        self.assertEqual(
            report["schema"],
            "rad-sim.contact-graph-consistency.v1",
        )
        self.assertTrue(report["summary"]["contactGraphConsistent"])
        self.assertEqual(report["formalization"]["targetId"], "contact_graph_consistency_gate")
        self.assertEqual(report["graph"]["activeBodyCount"], 3)
        self.assertEqual(report["graph"]["removedBodyCount"], 1)
        self.assertEqual(report["graph"]["totalEdgeCount"], 4)
        self.assertEqual(report["graph"]["activeEdgeCount"], 2)
        self.assertEqual(report["graph"]["deletedEdgeCount"], 2)
        self.assertEqual(report["graph"]["removedIncidentActiveEdgeCount"], 0)
        self.assertEqual(report["contact"]["removedActiveContactCount"], 0)
        self.assertEqual(report["support"]["supportCellCount"], 2)
        self.assertEqual(report["support"]["removedSupportCellCount"], 1)
        self.assertEqual(report["support"]["supportContactRecordCount"], 2)
        csv_text = export_contact_graph_consistency_csv(report)
        self.assertIn("contact_graph_consistent", csv_text)
        exported = json.loads(
            export_contact_graph_consistency_json(
                config,
                state,
                contact_report=contact_report,
                group_support=group_support,
            )
        )
        self.assertEqual(exported["summary"], report["summary"])

        inconsistent_contact = json.loads(json.dumps(contact_report))
        removed_cell = next(
            cell
            for cell in inconsistent_contact["cells"]
            if cell["row"] == 0 and cell["col"] == 1
        )
        removed_cell["mode"] = "engaged"
        removed_cell["contactState"]["unilateralContactActive"] = True
        inconsistent = contact_graph_consistency_report(
            config,
            state,
            contact_report=inconsistent_contact,
            group_support=group_support,
        )
        self.assertFalse(inconsistent["summary"]["contactGraphConsistent"])
        self.assertIn(
            "removedActiveContacts",
            inconsistent["summary"]["missingEvidence"],
        )

    def test_physical_realization_map_tracks_event_support_and_effects(self):
        config = LatticeConfig(rows=2, cols=2, pin_radius=0.10, hole_radius=0.12)
        state = LatticeState.uniform(config)
        events = (
            group_actuation_event(((0, 0), (0, 1)), alpha=0.2, z=0.05),
            lock_event((1, 1)),
            remove_cell_event((0, 1)),
        )
        report = physical_realization_map_report(
            config,
            state,
            event_sequence=events,
        )
        self.assertEqual(report["schema"], "rad-sim.physical-realization-map.v1")
        self.assertTrue(report["summary"]["physicalRealizationMapReady"])
        self.assertEqual(report["summary"]["abstractOperatorCount"], 3)
        self.assertEqual(report["summary"]["realizedOperatorCount"], 3)
        self.assertEqual(report["summary"]["supportRecordCount"], 3)
        self.assertEqual(report["support"]["supportCellCount"], 3)
        self.assertEqual(
            report["formalization"]["targetId"],
            "physical_realization_map_gate",
        )
        self.assertEqual(report["operators"][0]["operatorClass"], "group actuation operator")
        self.assertEqual(report["operators"][2]["operatorClass"], "graph deletion operator")
        self.assertEqual(report["operators"][2]["stateEffect"]["removedDeltaCount"], 1)
        self.assertTrue(report["contactGraph"]["contactGraphConsistent"])
        csv_text = export_physical_realization_map_csv(report)
        self.assertIn("hardware_channel", csv_text)
        self.assertIn("graph deletion operator", csv_text)
        exported = json.loads(
            export_physical_realization_map_json(
                config,
                state,
                event_sequence=events,
            )
        )
        self.assertEqual(exported["summary"], report["summary"])

        no_effect = physical_realization_map_report(
            config,
            state,
            event_sequence=(release_event((0, 0)),),
        )
        self.assertFalse(no_effect["summary"]["physicalRealizationMapReady"])
        self.assertIn("stateEffectRecords", no_effect["summary"]["missingEvidence"])

    def test_external_physics_engine_audit_tracks_independent_tool_readiness(self):
        config = LatticeConfig(rows=2, cols=2, pin_radius=0.10, hole_radius=0.12)
        state = LatticeState.uniform(config)
        state.z_actuator_grid[0, 0] = 0.2
        events = (group_actuation_event(((0, 0), (0, 1)), alpha=0.2, z=0.05),)
        report = external_physics_engine_audit_report(
            config,
            state,
            event_sequence=events,
            engine_availability={"mujoco": True},
        )
        self.assertEqual(report["schema"], "rad-sim.external-physics-engine-audit.v1")
        self.assertTrue(report["summary"]["externalPhysicsEngineAuditReady"])
        self.assertEqual(report["summary"]["availableEngineCount"], 1)
        self.assertGreaterEqual(report["summary"]["feasibleEngineCount"], 1)
        self.assertGreater(report["summary"]["contactModelRecordCount"], 0)
        self.assertEqual(
            report["formalization"]["targetId"],
            "external_physics_engine_audit_gate",
        )
        csv_text = export_external_physics_engine_audit_csv(report)
        self.assertIn("MuJoCo", csv_text)
        exported = json.loads(
            export_external_physics_engine_audit_json(
                config,
                state,
                event_sequence=events,
                engine_availability={"mujoco": True},
            )
        )
        self.assertEqual(exported["summary"], report["summary"])

        missing = external_physics_engine_audit_report(
            config,
            state,
            event_sequence=events,
            engine_availability={
                "mujoco": False,
                "pybullet": False,
                "pychrono": False,
            },
        )
        self.assertFalse(missing["summary"]["externalPhysicsEngineAuditReady"])
        self.assertIn("availableExternalEngine", missing["summary"]["missingEvidence"])
        self.assertIn("independentToolRecords", missing["summary"]["missingEvidence"])

    def test_mujoco_model_export_and_unavailable_run_are_claim_labeled(self):
        config = LatticeConfig(rows=2, cols=2, pin_radius=0.10, hole_radius=0.12)
        state = LatticeState.uniform(config)
        state.removed_mask[0, 1] = True
        load = LoadCase(
            fixed_cells=((0, 0),),
            external_forces={(1, 1): (0.0, 0.0, -0.25)},
        )
        report = mujoco_model_export_report(config, state, load_case=load)
        self.assertEqual(report["schema"], "rad-sim.mujoco-model-export.v1")
        self.assertTrue(report["summary"]["mujocoModelExportReady"])
        self.assertEqual(report["summary"]["bodyRecordCount"], 3)
        self.assertEqual(report["summary"]["fixedBodyCount"], 1)
        self.assertEqual(report["summary"]["removedBodyCount"], 1)
        self.assertEqual(report["summary"]["loadRecordCount"], 1)
        self.assertEqual(report["summary"]["pinRecordCount"], 12)
        self.assertEqual(report["summary"]["holeRecordCount"], 12)
        self.assertEqual(report["summary"]["clearanceRecordCount"], 12)
        self.assertEqual(report["summary"]["contactPairRecordCount"], 12)
        self.assertEqual(report["formalization"]["targetId"], "mujoco_model_export_gate")
        self.assertIn("<mujoco", report["xml"])
        self.assertIn("rad_cell_0_0", report["xml"])
        self.assertIn("rad_cell_0_0_nw_pin", report["xml"])
        self.assertIn("rad_cell_0_0_nw_hole_clearance", report["xml"])
        self.assertIn('friction="0.4 0.02 0.001"', report["xml"])
        self.assertIn('condim="3"', report["xml"])
        self.assertNotIn("rad_cell_0_1_plate", report["xml"])
        self.assertIn("normalized proxy geometry", report["claimLabels"]["geometry"])
        self.assertIn("rad_cell_1_1", export_mujoco_model_xml(config, state, load_case=load))
        exported = json.loads(export_mujoco_model_report_json(config, state, load_case=load))
        self.assertEqual(exported["summary"], report["summary"])
        contact_geometry = mujoco_pin_hole_contact_geometry_report(
            config,
            state,
            export_report=report,
        )
        self.assertEqual(
            contact_geometry["schema"],
            "rad-sim.mujoco-pin-hole-contact-geometry.v1",
        )
        self.assertTrue(
            contact_geometry["summary"]["mujocoPinHoleContactGeometryReady"]
        )
        self.assertEqual(contact_geometry["summary"]["pinRecordCount"], 12)
        self.assertEqual(contact_geometry["summary"]["holeRecordCount"], 12)
        self.assertEqual(contact_geometry["summary"]["activeContactPairCount"], 12)
        self.assertEqual(contact_geometry["summary"]["contactParameterRecordCount"], 12)
        self.assertEqual(contact_geometry["summary"]["frictionRecordCount"], 12)
        self.assertEqual(contact_geometry["summary"]["solverParameterRecordCount"], 12)
        self.assertAlmostEqual(contact_geometry["summary"]["minClearance"], 0.02)
        self.assertEqual(
            contact_geometry["formalization"]["targetId"],
            "mujoco_pin_hole_contact_geometry_gate",
        )
        self.assertTrue(
            all(
                not (record["row"] == 0 and record["col"] == 1)
                for record in contact_geometry["pins"]
            )
        )
        self.assertIn("clearance", export_mujoco_pin_hole_contact_geometry_csv(contact_geometry))
        exported_contact = json.loads(
            export_mujoco_pin_hole_contact_geometry_json(
                config,
                state,
                export_report=report,
            )
        )
        self.assertEqual(exported_contact["summary"], contact_geometry["summary"])
        parameters = mujoco_contact_parameter_report(
            config,
            state,
            export_report=report,
            contact_stiffness=1500.0,
            contact_damping=3.0,
            friction=(0.6, 0.03, 0.002),
            calibrated_contact=True,
        )
        self.assertEqual(
            parameters["schema"],
            "rad-sim.mujoco-contact-parameter-profile.v1",
        )
        self.assertTrue(parameters["summary"]["mujocoContactParameterProfileReady"])
        self.assertEqual(parameters["summary"]["parameterRecordCount"], 12)
        self.assertEqual(parameters["summary"]["calibratedRecordCount"], 12)
        self.assertEqual(
            parameters["formalization"]["targetId"],
            "mujoco_contact_parameter_profile_gate",
        )
        self.assertEqual(parameters["parameters"][0]["contactStiffness"], 1500.0)
        self.assertEqual(parameters["parameters"][0]["friction"], [0.6, 0.03, 0.002])
        self.assertIn("contact_stiffness", export_mujoco_contact_parameter_csv(parameters))
        exported_parameters = json.loads(
            export_mujoco_contact_parameter_json(
                config,
                state,
                export_report=report,
                contact_stiffness=1500.0,
                contact_damping=3.0,
                friction=(0.6, 0.03, 0.002),
                calibrated_contact=True,
            )
        )
        self.assertEqual(exported_parameters["summary"], parameters["summary"])
        calibration_packet = contact_parameter_calibration_packet(
            config,
            state,
            export_report=report,
            contact_geometry_report=contact_geometry,
            contact_parameter_report=parameters,
            repeat_count=2,
            fit_dataset_id="fit-contact",
            holdout_dataset_id="holdout-contact",
        )
        self.assertEqual(
            calibration_packet["schema"],
            "rad-sim.contact-parameter-calibration-packet.v1",
        )
        self.assertTrue(
            calibration_packet["summary"]["contactParameterCalibrationPacketReady"]
        )
        self.assertEqual(calibration_packet["summary"]["contactPairRecordCount"], 12)
        self.assertEqual(calibration_packet["summary"]["parameterRecordCount"], 12)
        self.assertEqual(calibration_packet["summary"]["fitTemplateRowCount"], 24)
        self.assertEqual(calibration_packet["summary"]["holdoutTemplateRowCount"], 24)
        self.assertIn(
            "measured_pin_hole_slip_mm",
            calibration_packet["measurementColumns"],
        )
        self.assertEqual(
            calibration_packet["formalization"]["targetId"],
            "contact_parameter_calibration_packet_completeness",
        )
        self.assertIn(
            "fitted_contact_stiffness",
            export_contact_parameter_calibration_packet_csv(calibration_packet),
        )
        exported_calibration_packet = json.loads(
            export_contact_parameter_calibration_packet_json(
                config,
                state,
                export_report=report,
                contact_parameter_report=parameters,
                repeat_count=2,
            )
        )
        self.assertEqual(
            exported_calibration_packet["summary"]["fitTemplateRowCount"],
            24,
        )
        blank_contact_results = contact_parameter_calibration_results_template(
            calibration_packet
        )
        self.assertEqual(
            blank_contact_results["schema"],
            "rad-sim.contact-parameter-calibration-results.v1",
        )
        self.assertEqual(len(blank_contact_results["measurements"]), 48)
        self.assertEqual(
            contact_parameter_calibration_results_from_json(
                json.dumps(blank_contact_results)
            )["schema"],
            "rad-sim.contact-parameter-calibration-results.v1",
        )
        self.assertEqual(
            json.loads(
                export_contact_parameter_calibration_results_template_json(
                    calibration_packet
                )
            )["sourcePacketSchema"],
            "rad-sim.contact-parameter-calibration-packet.v1",
        )
        blank_contact_validation = compare_contact_parameter_calibration_results(
            calibration_packet,
            blank_contact_results,
        )
        self.assertFalse(
            blank_contact_validation["summary"]["contactParameterBenchValidationPass"]
        )
        self.assertGreater(
            blank_contact_validation["metrics"]["missingMeasurementCount"],
            0,
        )
        blank_interval_calibration = contact_parameter_interval_calibration_report(
            calibration_packet,
            blank_contact_results,
            bench_validation=blank_contact_validation,
        )
        self.assertFalse(
            blank_interval_calibration["summary"][
                "contactParameterIntervalCalibrationReady"
            ]
        )
        self.assertIn(
            "contactParameterBenchValidationPass",
            blank_interval_calibration["summary"]["missingEvidence"],
        )
        filled_contact_results = deepcopy(blank_contact_results)
        for result_row in filled_contact_results["measurements"]:
            measured = result_row["requiredMeasurements"]
            measured["pinRadiusMm"] = calibration_packet["grid"]["pinRadius"]
            measured["holeRadiusMm"] = calibration_packet["grid"]["holeRadius"]
            measured["clearanceMm"] = calibration_packet["grid"]["pinHoleClearance"]
            measured["normalLoadN"] = 1.0
            measured["tangentialLoadN"] = 0.2
            measured["imposedZMm"] = 0.1
            measured["measuredPinHoleSlipMm"] = 0.02
            measured["measuredNormalForceN"] = 1.0
            measured["measuredTangentForceN"] = 0.6
            measured["measuredReboundRatio"] = 0.5
            measured["measuredContactDurationS"] = 0.1
            measured["measuredStaticFrictionCoeff"] = result_row["assumedFriction"][0]
            measured["measuredDynamicFrictionCoeff"] = result_row["assumedFriction"][1]
            measured["fittedContactStiffness"] = result_row[
                "assumedContactStiffness"
            ]
            measured["fittedContactDamping"] = result_row["assumedContactDamping"]
            measured["fittedSolrefTimeconst"] = result_row["assumedSolref"][0]
            measured["fittedSolrefDampingRatio"] = result_row["assumedSolref"][1]
            measured["fittedSolimpWidth"] = result_row["assumedSolimp"][0]
            measured["fixtureNotes"] = "synthetic validation row"
        contact_validation = compare_contact_parameter_calibration_results(
            calibration_packet,
            filled_contact_results,
            tolerance=1e-9,
            holdout_tolerance=1e-9,
        )
        self.assertEqual(
            contact_validation["schema"],
            "rad-sim.contact-parameter-bench-validation.v1",
        )
        self.assertTrue(
            contact_validation["summary"]["contactParameterBenchValidationPass"]
        )
        self.assertTrue(contact_validation["summary"]["fitParameterPass"])
        self.assertTrue(contact_validation["summary"]["holdoutParameterPass"])
        self.assertTrue(contact_validation["summary"]["independentHoldoutPass"])
        self.assertEqual(contact_validation["metrics"]["missingMeasurementCount"], 0)
        self.assertEqual(contact_validation["metrics"]["maxFitParameterResidual"], 0.0)
        self.assertEqual(
            contact_validation["formalization"]["targetId"],
            "contact_parameter_bench_validation_gate",
        )
        self.assertIn(
            "max_parameter_residual",
            export_contact_parameter_bench_validation_csv(contact_validation),
        )
        self.assertTrue(
            json.loads(
                export_contact_parameter_bench_validation_json(
                    calibration_packet,
                    filled_contact_results,
                    tolerance=1e-9,
                    holdout_tolerance=1e-9,
                )
            )["summary"]["contactParameterBenchValidationPass"]
        )
        interval_calibration = contact_parameter_interval_calibration_report(
            calibration_packet,
            filled_contact_results,
            bench_validation=contact_validation,
        )
        self.assertEqual(
            interval_calibration["schema"],
            "rad-sim.contact-parameter-interval-calibration.v1",
        )
        self.assertTrue(
            interval_calibration["summary"][
                "contactParameterIntervalCalibrationReady"
            ]
        )
        self.assertEqual(
            interval_calibration["summary"]["parameterIntervalCount"],
            10,
        )
        self.assertEqual(
            interval_calibration["summary"]["acceptedParameterIntervalCount"],
            10,
        )
        self.assertTrue(
            interval_calibration["summary"]["simulatorParametersInsideBounds"]
        )
        self.assertEqual(
            interval_calibration["formalization"]["targetId"],
            "contact_parameter_interval_calibration_gate",
        )
        self.assertIn(
            "lower_bound",
            export_contact_parameter_interval_calibration_csv(interval_calibration),
        )
        self.assertTrue(
            json.loads(
                export_contact_parameter_interval_calibration_json(
                    calibration_packet,
                    filled_contact_results,
                    bench_validation=contact_validation,
                )
            )["summary"]["contactParameterIntervalCalibrationReady"]
        )

        unavailable = mujoco_external_run_report(
            config,
            state,
            load_case=load,
            export_report=report,
            engine_availability={"mujoco": False},
        )
        self.assertEqual(unavailable["schema"], "rad-sim.mujoco-external-run.v1")
        self.assertFalse(unavailable["summary"]["mujocoExternalRunComplete"])
        self.assertIn("mujocoPackage", unavailable["summary"]["missingEvidence"])
        self.assertEqual(unavailable["formalization"]["targetId"], "mujoco_external_run_gate")
        exported_run = json.loads(
            export_mujoco_external_run_json(
                config,
                state,
                load_case=load,
                export_report=report,
                engine_availability={"mujoco": False},
            )
        )
        self.assertEqual(exported_run["summary"], unavailable["summary"])

    def test_mujoco_external_comparison_accepts_complete_imported_results(self):
        config = LatticeConfig(rows=2, cols=2, pin_radius=0.10, hole_radius=0.12)
        state = LatticeState.uniform(config)
        events = (group_actuation_event(((0, 0), (0, 1)), alpha=0.2, z=0.05),)
        export_report = mujoco_model_export_report(
            config,
            state,
            event_sequence=events,
            load_case=LoadCase(fixed_cells=((0, 0),)),
        )
        final_state = apply_event_sequence(config, state, events)
        sim = simulate_kinematic(config, final_state)
        bodies = []
        for body in export_report["bodies"]:
            row = int(body["row"])
            col = int(body["col"])
            position = [float(value) for value in sim.deformed_centers_3d[row, col]]
            bodies.append(
                {
                    "row": row,
                    "col": col,
                    "name": body["name"],
                    "initialPosition": body["position"],
                    "finalPosition": position,
                    "displacement": [0.0, 0.0, 0.0],
                    "fixed": body["fixed"],
                }
            )
        run_report = {
            "schema": "rad-sim.mujoco-external-run.v1",
            "summary": {
                "mujocoExternalRunComplete": True,
                "bodyResultCount": len(bodies),
                "expectedBodyResultCount": len(bodies),
                "missingEvidence": [],
            },
            "results": {"bodies": bodies},
            "modelExport": export_report,
        }
        comparison = mujoco_external_comparison_report(
            config,
            state,
            run_report,
            event_sequence=events,
            tolerance=1e-9,
        )
        self.assertEqual(comparison["schema"], "rad-sim.mujoco-external-comparison.v1")
        self.assertTrue(comparison["summary"]["mujocoExternalComparisonReady"])
        self.assertEqual(comparison["summary"]["comparisonRecordCount"], len(bodies))
        self.assertEqual(comparison["summary"]["maxPositionError"], 0.0)
        self.assertEqual(
            comparison["formalization"]["targetId"],
            "mujoco_external_comparison_gate",
        )
        csv_text = export_mujoco_external_comparison_csv(comparison)
        self.assertIn("position_error", csv_text)
        exported = json.loads(
            export_mujoco_external_comparison_json(
                config,
                state,
                run_report,
                event_sequence=events,
                tolerance=1e-9,
            )
        )
        self.assertEqual(exported["summary"], comparison["summary"])

        incomplete = dict(run_report)
        incomplete["summary"] = {"mujocoExternalRunComplete": False}
        missing = mujoco_external_comparison_report(config, state, incomplete, event_sequence=events)
        self.assertFalse(missing["summary"]["mujocoExternalComparisonReady"])
        self.assertIn("mujocoExternalRun", missing["summary"]["missingEvidence"])

    def test_equilibrium_relation_report_links_realization_energy_and_solver(self):
        config = LatticeConfig(rows=2, cols=2, pin_radius=0.10, hole_radius=0.12)
        state = LatticeState.uniform(config)
        events = (
            group_actuation_event(((0, 0), (0, 1)), alpha=0.2, z=0.05),
            lock_event((1, 1)),
        )
        report = equilibrium_relation_report(
            config,
            state,
            event_sequence=events,
            load_case=LoadCase(lock_stiffness=400.0, maxiter=250),
            tolerance=1e-8,
            residual_tolerance=1.0,
        )
        self.assertEqual(report["schema"], "rad-sim.equilibrium-relation.v1")
        self.assertTrue(report["summary"]["equilibriumRelationReady"])
        self.assertTrue(report["solver"]["success"])
        self.assertTrue(report["residual"]["passesTolerance"])
        self.assertEqual(
            report["energy"]["nonnegativeTermCount"],
            report["energy"]["requiredNonnegativeTermCount"],
        )
        self.assertEqual(
            report["formalization"]["targetId"],
            "equilibrium_relation_gate",
        )
        csv_text = export_equilibrium_relation_csv(report)
        self.assertIn("equilibrium_relation_ready", csv_text)
        self.assertIn("stored_energy", csv_text)
        exported = json.loads(
            export_equilibrium_relation_json(
                config,
                state,
                event_sequence=events,
                load_case=LoadCase(lock_stiffness=400.0, maxiter=250),
                tolerance=1e-8,
                residual_tolerance=1.0,
            )
        )
        self.assertEqual(exported["summary"], report["summary"])

        missing_realization = equilibrium_relation_report(
            config,
            state,
            event_sequence=events,
            realization_report={
                "schema": "rad-sim.physical-realization-map.v1",
                "summary": {"physicalRealizationMapReady": False},
            },
            residual_tolerance=1.0,
        )
        self.assertFalse(missing_realization["summary"]["equilibriumRelationReady"])
        self.assertIn(
            "physicalRealizationMap",
            missing_realization["summary"]["missingEvidence"],
        )

    def test_reachable_equilibrium_controllability_reports_topology_blocked_targets(self):
        config = LatticeConfig(
            rows=1,
            cols=3,
            backlash=0.0,
            coupling_gain=1.0,
            z_coupling_gain=1.0,
            pin_radius=0.10,
            hole_radius=0.12,
        )
        state = LatticeState.uniform(config)
        events = (
            group_actuation_event(((0, 0),), alpha=0.2, z=0.05),
            remove_cell_event((0, 1)),
        )
        report = reachable_equilibrium_controllability_report(
            config,
            state,
            event_sequence=events,
            target_cells=((0, 2),),
            residual_tolerance=1.0,
        )
        self.assertEqual(
            report["schema"],
            "rad-sim.reachable-equilibrium-controllability.v1",
        )
        self.assertTrue(report["summary"]["reachableEquilibriumControllabilityReady"])
        self.assertFalse(report["summary"]["targetFullyReachable"])
        self.assertEqual(report["actuatorBasis"]["cellCount"], 1)
        self.assertEqual(report["targets"]["targetCellCount"], 1)
        self.assertEqual(report["topology"]["topologyBlockedTargetCells"], 1)
        self.assertEqual(
            report["formalization"]["targetId"],
            "reachable_equilibrium_controllability_gate",
        )
        csv_text = export_reachable_equilibrium_controllability_csv(report)
        self.assertIn("reachable_equilibrium_controllability_ready", csv_text)
        self.assertIn("topology_blocked_target_cells", csv_text)
        exported = json.loads(
            export_reachable_equilibrium_controllability_json(
                config,
                state,
                event_sequence=events,
                target_cells=((0, 2),),
                residual_tolerance=1.0,
            )
        )
        self.assertEqual(exported["summary"], report["summary"])
        self.assertTrue(report["response"]["alphaReachableMap"][0][0])
        self.assertTrue(report["response"]["heightReachableMap"][0][0])

        protocol = reachable_equilibrium_bench_protocol(
            config,
            state,
            controllability_report=report,
            repeat_count=2,
        )
        self.assertEqual(
            protocol["schema"],
            "rad-sim.reachable-equilibrium-bench-protocol.v1",
        )
        self.assertTrue(protocol["summary"]["benchProtocolReady"])
        self.assertEqual(protocol["summary"]["topologyBlockedControlCount"], 1)
        self.assertTrue(protocol["topology"]["topologyPolicySatisfied"])
        self.assertTrue(protocol["targetSummaries"][0]["blockedByTopology"])
        self.assertIn(
            "topology_blocked_target_control",
            {step["id"] for step in protocol["steps"]},
        )
        self.assertEqual(
            protocol["formalization"]["targetId"],
            "reachable_equilibrium_bench_protocol_gate",
        )
        protocol_csv = export_reachable_equilibrium_bench_protocol_csv(protocol)
        self.assertIn("protocol_ready", protocol_csv)
        self.assertIn("topology_blocked_target_control", protocol_csv)
        exported_protocol = json.loads(
            export_reachable_equilibrium_bench_protocol_json(
                config,
                state,
                controllability_report=report,
                repeat_count=2,
            )
        )
        self.assertEqual(exported_protocol["summary"], protocol["summary"])
        template = reachable_equilibrium_bench_results_template(
            protocol,
            dataset_id="blocked-target-run",
        )
        self.assertEqual(
            template["schema"],
            "rad-sim.reachable-equilibrium-bench-results.v1",
        )
        self.assertGreater(len(template["measurements"]), 0)
        for row in template["measurements"]:
            modes = set(row["measurementMode"])
            for measurement in row["targetMeasurements"]:
                measurement["measuredTopologyComponentLabel"] = measurement[
                    "expectedTopologyComponentLabel"
                ]
                measurement["measuredAlphaDelta"] = (
                    0.02
                    if "alpha" in modes and measurement["predictedAlphaReachable"]
                    else 0.0
                )
                measurement["measuredHeightDelta"] = (
                    0.02
                    if "height" in modes and measurement["predictedHeightReachable"]
                    else 0.0
                )
        parsed_results = reachable_equilibrium_bench_results_from_json(
            json.dumps(template)
        )
        comparison = compare_reachable_equilibrium_bench_results(
            protocol,
            parsed_results,
            response_tolerance=1e-6,
            blocked_tolerance=1e-6,
            group_sequence_tolerance=1e-6,
        )
        self.assertEqual(
            comparison["schema"],
            "rad-sim.reachable-equilibrium-bench-comparison.v1",
        )
        self.assertTrue(comparison["summary"]["benchComparisonPass"])
        self.assertTrue(comparison["summary"]["topologyBlockedLeakagePass"])
        self.assertTrue(comparison["summary"]["groupSequencePass"])
        self.assertEqual(comparison["metrics"]["topologyLeakageCount"], 0)
        self.assertEqual(
            comparison["formalization"]["targetId"],
            "reachable_equilibrium_bench_validation_gate",
        )
        template_json = json.loads(
            export_reachable_equilibrium_bench_results_template_json(
                protocol,
                dataset_id="blocked-target-run",
            )
        )
        self.assertEqual(template_json["schema"], template["schema"])
        comparison_csv = export_reachable_equilibrium_bench_comparison_csv(
            comparison
        )
        self.assertIn("comparison_pass", comparison_csv)
        self.assertIn("group_target_reachability_probe", comparison_csv)
        comparison_json = json.loads(
            export_reachable_equilibrium_bench_comparison_json(
                protocol,
                parsed_results,
                response_tolerance=1e-6,
                blocked_tolerance=1e-6,
                group_sequence_tolerance=1e-6,
            )
        )
        self.assertEqual(comparison_json["summary"], comparison["summary"])
        amplitude = reachable_equilibrium_amplitude_calibration_report(
            protocol,
            parsed_results,
            comparison_report=comparison,
            response_tolerance=1e-6,
            blocked_tolerance=1e-6,
            group_sequence_tolerance=1e-6,
            confidence_sigma=2.0,
        )
        self.assertEqual(
            amplitude["schema"],
            "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
        )
        self.assertTrue(amplitude["summary"]["amplitudeCalibrationReady"])
        self.assertTrue(amplitude["summary"]["topologyLeakageBandsPass"])
        self.assertTrue(amplitude["summary"]["groupSequenceResidualPass"])
        self.assertGreater(amplitude["metrics"]["amplitudeEstimateCount"], 0)
        self.assertGreater(amplitude["metrics"]["repeatedTrialGroupCount"], 0)
        self.assertEqual(amplitude["metrics"]["failedTopologyBandCount"], 0)
        self.assertEqual(
            amplitude["formalization"]["targetId"],
            "reachable_equilibrium_amplitude_calibration_gate",
        )
        amplitude_csv = export_reachable_equilibrium_amplitude_calibration_csv(
            amplitude
        )
        self.assertIn("calibration_ready", amplitude_csv)
        self.assertIn("alpha_uncertainty", amplitude_csv)
        amplitude_json = json.loads(
            export_reachable_equilibrium_amplitude_calibration_json(
                protocol,
                parsed_results,
                comparison_report=comparison,
                response_tolerance=1e-6,
                blocked_tolerance=1e-6,
                group_sequence_tolerance=1e-6,
                confidence_sigma=2.0,
            )
        )
        self.assertEqual(amplitude_json["summary"], amplitude["summary"])
        profile_report = reachable_equilibrium_controllability_report(
            config,
            state,
            event_sequence=(group_actuation_event(((0, 0),), alpha=0.2, z=0.05),),
            target_cells=((0, 0),),
            residual_tolerance=1.0,
        )
        profile_protocol = reachable_equilibrium_bench_protocol(
            config,
            state,
            controllability_report=profile_report,
            repeat_count=2,
        )
        profile_template = reachable_equilibrium_bench_results_template(
            profile_protocol,
            dataset_id="reachable-target-run",
        )
        for row in profile_template["measurements"]:
            modes = set(row["measurementMode"])
            for measurement in row["targetMeasurements"]:
                measurement["measuredTopologyComponentLabel"] = measurement[
                    "expectedTopologyComponentLabel"
                ]
                measurement["measuredAlphaDelta"] = (
                    0.02
                    if "alpha" in modes and measurement["predictedAlphaReachable"]
                    else 0.0
                )
                measurement["measuredHeightDelta"] = (
                    0.02
                    if "height" in modes and measurement["predictedHeightReachable"]
                    else 0.0
                )
        profile_comparison = compare_reachable_equilibrium_bench_results(
            profile_protocol,
            profile_template,
            response_tolerance=1e-6,
            blocked_tolerance=1e-6,
            group_sequence_tolerance=1e-6,
        )
        profile_amplitude = reachable_equilibrium_amplitude_calibration_report(
            profile_protocol,
            profile_template,
            comparison_report=profile_comparison,
            response_tolerance=1e-6,
            blocked_tolerance=1e-6,
            group_sequence_tolerance=1e-6,
            confidence_sigma=2.0,
        )
        profile = reachable_equilibrium_empirical_profile_from_amplitude(
            profile_amplitude,
            safety_factor=1.5,
            min_samples=1,
        )
        self.assertEqual(
            profile["schema"],
            "rad-sim.reachable-equilibrium-empirical-profile.v1",
        )
        self.assertTrue(profile["summary"]["empiricalProfileReady"])
        self.assertGreater(profile["metrics"]["safeProposalCount"], 0)
        self.assertEqual(profile["holdoutValidation"]["holdoutPass"], True)
        self.assertIn(
            "alphaResponseScale",
            {update["name"] for update in profile["recommendedUpdates"]},
        )
        self.assertEqual(
            profile["formalization"]["targetId"],
            "reachable_equilibrium_empirical_profile_gate",
        )
        profile_csv = export_reachable_equilibrium_empirical_profile_csv(profile)
        self.assertIn("profile_ready", profile_csv)
        self.assertIn("alphaResponseScale", profile_csv)
        exported_profile = json.loads(
            export_reachable_equilibrium_empirical_profile_json(
                profile_amplitude,
                safety_factor=1.5,
                min_samples=1,
            )
        )
        self.assertEqual(exported_profile["summary"], profile["summary"])
        target_height = np.zeros((config.rows, config.cols), dtype=float)
        target_height[0, 0] = 0.03
        inverse = solve_inverse_design(
            config,
            target_height=target_height,
            actuator_cells=((0, 0),),
            include_alpha=False,
            include_z=True,
        )
        before_z_gain = config.z_coupling_gain
        profile_inverse = reachable_equilibrium_profile_inverse_report(
            inverse,
            profile,
        )
        self.assertEqual(
            profile_inverse["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse.v1",
        )
        self.assertTrue(profile_inverse["summary"]["profileInverseReady"])
        self.assertEqual(profile_inverse["target"]["targetCellCount"], 1)
        self.assertGreaterEqual(
            profile_inverse["profile"]["safeBoundedProposalCount"], 1
        )
        self.assertEqual(
            profile_inverse["formalization"]["targetId"],
            "reachable_equilibrium_profile_inverse_gate",
        )
        self.assertEqual(config.z_coupling_gain, before_z_gain)
        profile_inverse_csv = export_reachable_equilibrium_profile_inverse_csv(
            profile_inverse
        )
        self.assertIn("profile_inverse_ready", profile_inverse_csv)
        self.assertIn(
            "rad-sim.reachable-equilibrium-profile-inverse.v1",
            profile_inverse_csv,
        )
        exported_profile_inverse = json.loads(
            export_reachable_equilibrium_profile_inverse_json(inverse, profile)
        )
        self.assertEqual(
            exported_profile_inverse["summary"],
            profile_inverse["summary"],
        )
        acceptance = reachable_equilibrium_profile_inverse_acceptance_report(
            profile_inverse,
            max_weighted_residual_score=1.0,
            max_band_failures=1,
            max_active_actuators=1,
        )
        self.assertEqual(
            acceptance["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
        )
        self.assertTrue(acceptance["summary"]["profileInverseAcceptanceReady"])
        self.assertEqual(acceptance["decision"]["decision"], "accept-for-preview")
        self.assertEqual(
            acceptance["formalization"]["targetId"],
            "reachable_equilibrium_profile_inverse_acceptance_gate",
        )
        acceptance_csv = export_reachable_equilibrium_profile_inverse_acceptance_csv(
            acceptance
        )
        self.assertIn("acceptance_ready", acceptance_csv)
        self.assertIn("accept-for-preview", acceptance_csv)
        exported_acceptance = json.loads(
            export_reachable_equilibrium_profile_inverse_acceptance_json(
                profile_inverse,
                max_weighted_residual_score=1.0,
                max_band_failures=1,
                max_active_actuators=1,
            )
        )
        self.assertEqual(exported_acceptance["summary"], acceptance["summary"])
        review_candidate = json.loads(json.dumps(profile_inverse))
        review_candidate["summary"]["profileWeightedResidualScore"] = 5.0
        review = reachable_equilibrium_profile_inverse_acceptance_report(
            review_candidate,
            max_weighted_residual_score=1.0,
        )
        self.assertFalse(review["summary"]["profileInverseAcceptanceReady"])
        self.assertEqual(review["decision"]["decision"], "review-required")
        self.assertIn(
            "profileWeightedResidualScore",
            review["decision"]["failedCriteria"],
        )
        packet = reachable_equilibrium_profile_inverse_preview_packet(
            inverse,
            profile_inverse,
            acceptance,
            packet_id="unit-test-profile-inverse-preview",
            notes="unit test packet",
        )
        self.assertEqual(
            packet["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
        )
        self.assertTrue(packet["summary"]["profileInversePreviewPacketReady"])
        self.assertEqual(packet["summary"]["commandCount"], 1)
        self.assertEqual(packet["summary"]["eventCount"], 1)
        self.assertEqual(packet["commands"][0]["row"], 0)
        self.assertEqual(packet["commands"][0]["col"], 0)
        self.assertTrue(packet["reviewProtocol"]["requiresHumanReviewBeforeHardware"])
        self.assertEqual(
            packet["formalization"]["targetId"],
            "reachable_equilibrium_profile_inverse_preview_packet_gate",
        )
        packet_csv = export_reachable_equilibrium_profile_inverse_preview_packet_csv(
            packet
        )
        self.assertIn("packet_ready", packet_csv)
        self.assertIn("unit-test-profile-inverse-preview", packet["packetId"])
        exported_packet = json.loads(
            export_reachable_equilibrium_profile_inverse_preview_packet_json(
                inverse,
                profile_inverse,
                acceptance,
                packet_id="unit-test-profile-inverse-preview",
            )
        )
        self.assertEqual(exported_packet["summary"], packet["summary"])
        rejected_packet = reachable_equilibrium_profile_inverse_preview_packet(
            inverse,
            profile_inverse,
            review,
        )
        self.assertFalse(
            rejected_packet["summary"]["profileInversePreviewPacketReady"]
        )
        self.assertIn(
            "profileInverseAcceptanceReady",
            rejected_packet["summary"]["missingEvidence"],
        )
        replay = reachable_equilibrium_profile_inverse_preview_replay_report(
            config,
            packet,
            target_height=target_height,
            residual_agreement_tolerance=1e-8,
        )
        self.assertEqual(
            replay["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
        )
        self.assertTrue(replay["summary"]["profileInversePreviewReplayReady"])
        self.assertEqual(replay["commands"]["replayedCommandCount"], 1)
        self.assertEqual(replay["commands"]["invalidCommandCount"], 0)
        self.assertTrue(replay["residuals"]["residualAgreementPass"])
        self.assertEqual(
            replay["formalization"]["targetId"],
            "reachable_equilibrium_profile_inverse_preview_replay_gate",
        )
        replay_csv = export_reachable_equilibrium_profile_inverse_preview_replay_csv(
            replay
        )
        self.assertIn("replay_ready", replay_csv)
        exported_replay = json.loads(
            export_reachable_equilibrium_profile_inverse_preview_replay_json(
                config,
                packet,
                target_height=target_height,
                residual_agreement_tolerance=1e-8,
            )
        )
        self.assertEqual(exported_replay["summary"], replay["summary"])
        physical = reachable_equilibrium_profile_inverse_preview_physical_report(
            config,
            packet,
            target_height=target_height,
            residual_agreement_tolerance=1e-8,
        )
        self.assertEqual(
            physical["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
        )
        self.assertTrue(
            physical["summary"]["profileInversePreviewPhysicalReady"]
        )
        self.assertTrue(physical["physical"]["physicalSuccess"])
        self.assertEqual(physical["commands"]["replayedCommandCount"], 1)
        self.assertGreaterEqual(physical["comparison"]["heightRmsModelError"], 0.0)
        self.assertGreaterEqual(physical["comparison"]["centerRmsModelError"], 0.0)
        self.assertGreaterEqual(physical["physical"]["physicalEnergy"], 0.0)
        self.assertEqual(
            physical["formalization"]["targetId"],
            "reachable_equilibrium_profile_inverse_preview_physical_gate",
        )
        physical_csv = export_reachable_equilibrium_profile_inverse_preview_physical_csv(
            physical
        )
        self.assertIn("physical_ready", physical_csv)
        exported_physical = json.loads(
            export_reachable_equilibrium_profile_inverse_preview_physical_json(
                config,
                packet,
                target_height=target_height,
                residual_agreement_tolerance=1e-8,
            )
        )
        self.assertEqual(exported_physical["summary"], physical["summary"])
        failed_physical = reachable_equilibrium_profile_inverse_preview_physical_report(
            config,
            packet,
            target_height=target_height,
            max_height_model_error=0.0,
        )
        if failed_physical["comparison"]["maxAbsHeightModelError"] > 1e-9:
            self.assertFalse(
                failed_physical["summary"]["profileInversePreviewPhysicalReady"]
            )
            self.assertIn(
                "modelErrorThreshold",
                failed_physical["summary"]["missingEvidence"],
            )
        invalid_packet = json.loads(json.dumps(packet))
        invalid_packet["commands"][0]["row"] = config.rows + 10
        invalid_replay = reachable_equilibrium_profile_inverse_preview_replay_report(
            config,
            invalid_packet,
            target_height=target_height,
        )
        self.assertFalse(
            invalid_replay["summary"]["profileInversePreviewReplayReady"]
        )
        self.assertIn(
            "validCommandRecords",
            invalid_replay["summary"]["missingEvidence"],
        )
        needs_holdout = reachable_equilibrium_empirical_profile_from_amplitude(
            profile_amplitude,
            require_holdout=True,
        )
        self.assertFalse(needs_holdout["summary"]["empiricalProfileReady"])
        self.assertIn("holdoutValidation", needs_holdout["summary"]["missingEvidence"])

        leaked = json.loads(json.dumps(template))
        for row in leaked["measurements"]:
            if row["stepId"] == "topology_blocked_target_control":
                row["targetMeasurements"][0]["measuredHeightDelta"] = 0.01
                break
        leaked_comparison = compare_reachable_equilibrium_bench_results(
            protocol,
            leaked,
            blocked_tolerance=1e-6,
        )
        self.assertFalse(leaked_comparison["summary"]["benchComparisonPass"])
        self.assertFalse(
            leaked_comparison["summary"]["topologyBlockedLeakagePass"]
        )
        self.assertGreater(leaked_comparison["metrics"]["topologyLeakageCount"], 0)
        leaked_amplitude = reachable_equilibrium_amplitude_calibration_report(
            protocol,
            leaked,
            comparison_report=leaked_comparison,
            blocked_tolerance=1e-6,
        )
        self.assertFalse(
            leaked_amplitude["summary"]["amplitudeCalibrationReady"]
        )
        self.assertIn(
            "topologyLeakageBand",
            leaked_amplitude["summary"]["missingEvidence"],
        )

        strict = reachable_equilibrium_controllability_report(
            config,
            state,
            event_sequence=events,
            target_cells=((0, 2),),
            require_full_target_reachability=True,
            residual_tolerance=1.0,
        )
        self.assertFalse(
            strict["summary"]["reachableEquilibriumControllabilityReady"]
        )
        self.assertIn("targetReachability", strict["summary"]["missingEvidence"])

    def test_finite_die_off_radius_respects_max_coupling_steps(self):
        config = LatticeConfig(rows=5, cols=5, backlash=0.0, max_coupling_steps=1)
        op = DeadZonePropagationOperator(
            "test",
            dead_zone=config.backlash,
            gain=config.coupling_gain,
            max_steps=config.max_coupling_steps,
        )
        result = op.propagate(config.rows, config.cols, [(2, 2, 0.3)])
        self.assertEqual(finite_die_off_radius(result.die_off), 1)
        self.assertTrue(np.isinf(result.die_off[0, 0]))

    def test_removed_cell_blocks_dead_zone_propagation_path(self):
        op = DeadZonePropagationOperator("test", dead_zone=0.0, gain=1.0)
        active = np.array([[True, False, True]], dtype=bool)
        result = op.propagate(1, 3, [(0, 0, 0.3)], active_mask=active)
        self.assertAlmostEqual(result.field[0, 0], 0.3)
        self.assertAlmostEqual(result.field[0, 1], 0.0)
        self.assertAlmostEqual(result.field[0, 2], 0.0)
        self.assertTrue(np.isinf(result.die_off[0, 2]))

    def test_lock_event_commits_current_alpha_and_height_state(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        state = LatticeState.uniform(config)
        actuated = apply_event_sequence(
            config,
            state,
            (local_actuation_event((1, 1), alpha=-0.3, z=0.2),),
        )
        actuated_fields = evaluate_programmable_operators(config, actuated)
        final = apply_event_sequence(
            config,
            actuated,
            (
                lock_event((1, 1)),
                clear_actuation_event(),
            ),
        )
        self.assertTrue(final.locked_mask[1, 1])
        self.assertAlmostEqual(final.alpha_grid[1, 1], 0.7)
        self.assertAlmostEqual(final.lock_z_grid[1, 1], actuated_fields["height"][1, 1])
        fields = evaluate_programmable_operators(config, final)
        self.assertAlmostEqual(fields["alpha"][1, 1], 0.7)
        self.assertAlmostEqual(fields["height"][1, 1], actuated_fields["height"][1, 1])
        self.assertNotAlmostEqual(fields["height"][1, 1], 0.0)

    def test_lock_and_actuation_event_order_is_noncommutative(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        state = LatticeState.uniform(config)
        diagnostic = compare_event_order(
            config,
            state,
            local_actuation_event((1, 1), alpha=-0.3),
            lock_event((1, 1)),
        )
        self.assertTrue(diagnostic.mode_commutes)
        self.assertFalse(diagnostic.command_commutes)
        self.assertFalse(diagnostic.alpha_grid_commutes)
        self.assertGreater(diagnostic.final_alpha_error, 0.1)
        self.assertGreater(diagnostic.final_height_error, 0.1)

    def test_sequence_order_diagnostic_reports_adjacent_swap_sensitivity(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        state = LatticeState.uniform(config)
        diagnostic = compare_sequence_order(
            config,
            state,
            (
                local_actuation_event((1, 1), alpha=-0.3, z=0.2),
                lock_event((1, 1)),
                clear_actuation_event((1, 1)),
            ),
        )
        self.assertEqual(diagnostic.event_count, 3)
        self.assertEqual(diagnostic.adjacent_pair_count, 2)
        self.assertGreater(diagnostic.noncommuting_adjacent_pairs, 0)
        self.assertTrue(diagnostic.order_sensitive)
        self.assertGreater(diagnostic.max_order_error, 0.1)
        self.assertGreater(diagnostic.reverse.final_alpha_error, 0.1)
        self.assertGreater(diagnostic.reverse.final_height_error, 0.1)
        self.assertEqual(len(diagnostic.adjacent), 2)
        self.assertTrue(any(item.sensitive for item in diagnostic.adjacent))

    def test_release_event_allows_later_actuation(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        state = LatticeState.uniform(config)
        final = apply_event_sequence(
            config,
            state,
            (
                lock_event((1, 1)),
                release_event((1, 1)),
                local_actuation_event((1, 1), alpha=-0.25),
            ),
        )
        self.assertFalse(final.locked_mask[1, 1])
        fields = evaluate_programmable_operators(config, final)
        self.assertLess(fields["alpha"][1, 1], config.initial_alpha)

    def test_clear_actuation_event_can_clear_one_cell_or_all_cells(self):
        config = LatticeConfig(rows=2, cols=2)
        state = apply_event_sequence(
            config,
            LatticeState.uniform(config),
            (
                local_actuation_event((0, 0), alpha=0.2, z=0.1),
                local_actuation_event((1, 1), alpha=-0.3, z=-0.2),
                clear_actuation_event((0, 0)),
            ),
        )
        self.assertAlmostEqual(state.actuator_grid[0, 0], 0.0)
        self.assertAlmostEqual(state.z_actuator_grid[0, 0], 0.0)
        self.assertNotAlmostEqual(state.actuator_grid[1, 1], 0.0)
        final = apply_event_sequence(config, state, (clear_actuation_event(),))
        np.testing.assert_allclose(final.actuator_grid, 0.0)
        np.testing.assert_allclose(final.z_actuator_grid, 0.0)

    def test_remove_and_restore_cell_events_change_topology_state(self):
        config = LatticeConfig(rows=1, cols=3, backlash=0.0, coupling_gain=1.0)
        final = apply_event_sequence(
            config,
            LatticeState.uniform(config),
            (
                local_actuation_event((0, 0), alpha=0.3, z=0.2),
                remove_cell_event((0, 1)),
                local_actuation_event((0, 1), alpha=0.3, z=0.2),
            ),
        )
        self.assertTrue(final.removed_mask[0, 1])
        self.assertAlmostEqual(final.actuator_grid[0, 1], 0.0)
        self.assertAlmostEqual(final.z_actuator_grid[0, 1], 0.0)
        fields = evaluate_programmable_operators(config, final)
        self.assertAlmostEqual(fields["actuator_influence"][0, 2], 0.0)
        self.assertTrue(np.isinf(fields["die_off"][0, 2]))

        restored = apply_event_sequence(config, final, (restore_cell_event((0, 1)),))
        self.assertFalse(restored.removed_mask[0, 1])

    def test_topology_diagnostic_labels_components_after_removal(self):
        config = LatticeConfig(rows=1, cols=3)
        state = apply_event_sequence(
            config,
            LatticeState.uniform(config),
            (remove_cell_event((0, 1)),),
        )
        topology = lattice_topology_diagnostic(config, state)
        self.assertEqual(topology["component_count"], 2)
        self.assertEqual(topology["component_sizes"], (1, 1))
        self.assertEqual(topology["removed_cell_count"], 1)
        self.assertEqual(topology["deleted_edge_count"], 2)
        self.assertEqual(topology["component_labels"][0, 1], -1)

    def test_group_actuation_event_applies_to_multiple_present_cells(self):
        config = LatticeConfig(rows=2, cols=2)
        state = apply_event_sequence(
            config,
            LatticeState.uniform(config),
            (
                remove_cell_event((1, 1)),
                group_actuation_event(((0, 0), (1, 1), (1, 0)), alpha=-0.2, z=0.15),
            ),
        )
        self.assertAlmostEqual(state.actuator_grid[0, 0], -0.2)
        self.assertAlmostEqual(state.z_actuator_grid[0, 0], 0.15)
        self.assertAlmostEqual(state.actuator_grid[1, 0], -0.2)
        self.assertAlmostEqual(state.z_actuator_grid[1, 0], 0.15)
        self.assertAlmostEqual(state.actuator_grid[1, 1], 0.0)
        self.assertAlmostEqual(state.z_actuator_grid[1, 1], 0.0)

    def test_group_actuation_decomposes_into_local_sequence(self):
        config = LatticeConfig(rows=2, cols=2)
        state = apply_event_sequence(
            config,
            LatticeState.uniform(config),
            (remove_cell_event((1, 1)),),
        )
        diagnostic = compare_group_actuation_decomposition(
            config,
            state,
            ((0, 0), (1, 0), (1, 1)),
            alpha=-0.2,
            z=0.15,
        )
        self.assertTrue(diagnostic.decomposes)
        self.assertEqual(diagnostic.cells, ((0, 0), (1, 0), (1, 1)))
        self.assertEqual(diagnostic.applied_cells, ((0, 0), (1, 0)))
        self.assertEqual(diagnostic.skipped_removed_cells, ((1, 1),))
        self.assertEqual(len(diagnostic.local_events), 3)
        self.assertEqual(diagnostic.distance.mode_changes, 0)
        self.assertLessEqual(diagnostic.distance.command_error, 1e-9)
        np.testing.assert_allclose(
            diagnostic.simultaneous_state.actuator_grid,
            diagnostic.sequenced_state.actuator_grid,
        )

    def test_group_actuation_under_removal_reports_support_and_topology_loss(self):
        config = LatticeConfig(rows=1, cols=3, backlash=0.0, coupling_gain=0.8)
        diagnostic = compare_group_actuation_under_removal(
            config,
            LatticeState.uniform(config),
            ((0, 0), (0, 1), (0, 2)),
            ((0, 1),),
            alpha=-0.2,
            z=0.15,
        )
        self.assertTrue(diagnostic.intact.decomposes)
        self.assertTrue(diagnostic.removed.decomposes)
        self.assertEqual(diagnostic.lost_applied_cells, ((0, 1),))
        self.assertEqual(diagnostic.newly_skipped_cells, ((0, 1),))
        self.assertEqual(diagnostic.removed.skipped_removed_cells, ((0, 1),))
        self.assertEqual(diagnostic.component_count_delta, 1)
        self.assertEqual(diagnostic.deleted_edge_delta, 2)
        self.assertGreater(diagnostic.distance.command_error, 0)
        self.assertGreater(diagnostic.distance.final_alpha_error, 0)
        self.assertEqual(diagnostic.removed_topology["component_count"], 2)

    def test_vertical_residual_under_removal_reports_clearance_topology_and_contact(self):
        config = LatticeConfig(
            rows=1,
            cols=5,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        diagnostic = compare_vertical_residual_under_removal(
            config,
            LatticeState.uniform(config),
            ((0, 0), (0, 2)),
            ((0, 2),),
            z=0.4,
            fixed_cells=((0, 0),),
            external_z_load=-0.2,
        )
        self.assertEqual(diagnostic.fixed_cells, ((0, 0),))
        self.assertAlmostEqual(diagnostic.external_z_load, -0.2)
        self.assertEqual(diagnostic.removed_source_cells, ((0, 2),))
        self.assertEqual(diagnostic.affected_neighbor_cells, ((0, 1),))
        self.assertIn((0, 3), diagnostic.topology_blocked_cells)
        self.assertIn((0, 4), diagnostic.topology_blocked_cells)
        self.assertEqual(diagnostic.component_count_delta, 1)
        self.assertEqual(diagnostic.deleted_edge_delta, 2)
        self.assertGreater(diagnostic.intact_die_off_radius, diagnostic.removed_die_off_radius)
        self.assertGreater(diagnostic.max_abs_residual_delta, 0.0)
        self.assertAlmostEqual(diagnostic.clearance, 0.02)
        self.assertEqual(
            diagnostic.intact_contact_engaged_cells,
            ((0, 0), (0, 2)),
        )
        self.assertEqual(diagnostic.removed_contact_engaged_cells, ((0, 0),))
        self.assertAlmostEqual(
            diagnostic.intact_contact_penalty_energy,
            2 * 0.5 * (0.4 - config.pin_hole_clearance) ** 2,
        )
        self.assertLess(diagnostic.contact_penalty_delta, 0.0)
        self.assertNotIn((0, 0), diagnostic.intact_load_active_cells)
        self.assertIn((0, 1), diagnostic.intact_load_active_cells)
        self.assertNotIn((0, 2), diagnostic.removed_load_active_cells)
        self.assertGreater(
            diagnostic.intact_load_work_magnitude,
            diagnostic.removed_load_work_magnitude,
        )
        self.assertLess(diagnostic.load_work_magnitude_delta, 0.0)
        self.assertNotIn((0, 0), diagnostic.intact_height_contact_engaged_cells)
        self.assertGreater(
            diagnostic.intact_height_contact_penalty_energy,
            diagnostic.removed_height_contact_penalty_energy,
        )
        self.assertLess(diagnostic.height_contact_penalty_delta, 0.0)

    def test_vertical_residual_spring_hinge_comparison_uses_removed_topology(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        comparison = compare_vertical_residual_spring_hinge_3d(
            config,
            LatticeState.uniform(config),
            ((0, 0), (0, 2)),
            ((0, 2),),
            z=0.25,
            fixed_cells=((0, 0),),
            external_z_load=-0.05,
            load_case=LoadCase(
                fixed_cells=((0, 0),),
                axial_stiffness=20.0,
                hinge_stiffness=1.0,
                lock_stiffness=250.0,
                maxiter=250,
            ),
        )
        self.assertTrue(comparison.physical_success)
        self.assertEqual(comparison.vertical_diagnostic.removed_source_cells, ((0, 2),))
        self.assertEqual(comparison.intact_solver_removed_cells, 0)
        self.assertEqual(comparison.removed_solver_removed_cells, 1)
        self.assertGreater(comparison.intact_spring_edges, comparison.removed_spring_edges)
        self.assertEqual(comparison.spring_edge_delta, -2)
        self.assertEqual(comparison.intact_solver_load_cells, ((0, 1), (0, 2), (0, 3)))
        self.assertEqual(comparison.removed_solver_load_cells, ((0, 1),))
        self.assertEqual(comparison.removed_unanchored_load_cells, ((0, 3),))
        self.assertEqual(comparison.intact_physical_height.shape, (1, 4))
        self.assertEqual(comparison.removed_physical_height.shape, (1, 4))
        self.assertGreaterEqual(comparison.intact_height_rms_model_error, 0.0)
        self.assertGreaterEqual(comparison.removed_height_rms_model_error, 0.0)
        self.assertGreaterEqual(comparison.delta_rms_model_error, 0.0)
        self.assertGreaterEqual(comparison.max_abs_delta_model_error, 0.0)
        self.assertGreaterEqual(comparison.intact_energy_breakdown.stored_energy, 0.0)
        self.assertGreaterEqual(comparison.removed_energy_breakdown.stored_energy, 0.0)
        self.assertAlmostEqual(
            comparison.intact_energy_breakdown.objective_energy,
            comparison.intact_energy_breakdown.stored_energy
            + comparison.intact_energy_breakdown.external_potential_energy,
        )
        self.assertIn("unanchored load components", comparison.solver_topology_note)
        payload = vertical_removal_physical_comparison_to_dict(
            comparison, include_fields=True
        )
        self.assertEqual(
            payload["schema"],
            "rad-sim.vertical-removal-physical-comparison.v1",
        )
        self.assertEqual(
            payload["claimLabels"]["gravityContact"],
            "experimentally unvalidated physical assumption",
        )
        self.assertEqual(payload["solver"]["springEdgeDelta"], -2)
        self.assertIn("intactEnergyBreakdown", payload["solver"])
        self.assertIn("mechanicsCertificate", payload)
        self.assertEqual(
            payload["mechanicsCertificate"]["schema"],
            "rad-sim.mechanics-energy-certificate.v1",
        )
        self.assertTrue(payload["mechanicsCertificate"]["passesNonnegativeCheck"])
        self.assertIn(
            "mechanicalStoredEnergyNat_nonnegative",
            payload["mechanicsCertificate"]["leanTheorems"],
        )
        certificate = mechanics_energy_certificate_to_dict(comparison)
        self.assertEqual(certificate["schema"], payload["mechanicsCertificate"]["schema"])
        self.assertIn("signedTerms", certificate["intact"])
        self.assertIn("externalPotentialEnergy", certificate["intact"]["signedTerms"])
        self.assertAlmostEqual(
            certificate["intact"]["nonnegativeProxyTotal"],
            sum(certificate["intact"]["nonnegativeTerms"].values()),
        )
        self.assertEqual(
            payload["solver"]["removedUnanchoredLoadCells"],
            [{"row": 0, "col": 3}],
        )
        self.assertIn("fields", payload)
        self.assertIn("physicalVsKinematicDeltaError", payload["fields"])
        exported = json.loads(
            export_vertical_removal_physical_comparison_json(
                comparison, include_fields=False
            )
        )
        self.assertEqual(exported["schema"], payload["schema"])
        self.assertNotIn("fields", exported)

    def test_vertical_load_physical_preview_report_exports_named_scenarios(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        report = build_vertical_load_physical_preview_report(
            config,
            load_case=LoadCase(
                fixed_cells=((0, 0),),
                axial_stiffness=20.0,
                hinge_stiffness=1.0,
                lock_stiffness=250.0,
                maxiter=250,
            ),
            include_fields=False,
        )
        self.assertEqual(
            report["schema"],
            "rad-sim.vertical-load-physical-preview-report.v1",
        )
        self.assertEqual(report["summary"]["scenarioCount"], 3)
        self.assertTrue(report["summary"]["allPhysicalSuccess"])
        names = [scenario["name"] for scenario in report["scenarios"]]
        self.assertEqual(
            names,
            ["one-cell-z-load", "two-cell-z-load", "removed-middle-z-load"],
        )
        removed_middle = report["scenarios"][2]["comparison"]
        self.assertEqual(removed_middle["solver"]["springEdgeDelta"], -2)
        self.assertGreaterEqual(
            report["summary"]["totalUnanchoredLoadCells"],
            1,
        )
        self.assertNotIn("fields", removed_middle)
        exported = json.loads(export_vertical_load_physical_preview_report_json(report))
        self.assertEqual(exported["schema"], report["schema"])
        csv_rows = list(
            csv.DictReader(
                io.StringIO(export_vertical_load_physical_preview_report_csv(report))
            )
        )
        self.assertEqual(len(csv_rows), 3)
        self.assertEqual(csv_rows[0]["scenario"], "one-cell-z-load")
        self.assertEqual(csv_rows[0]["fixed_cells"], "0:0")
        self.assertEqual(csv_rows[0]["z_command"], "0.25")
        self.assertEqual(csv_rows[0]["external_z_load"], "-0.05")
        self.assertEqual(csv_rows[2]["scenario"], "removed-middle-z-load")
        self.assertEqual(csv_rows[2]["removed_cells"], "0:1")
        self.assertEqual(csv_rows[2]["spring_edge_delta"], "-2")
        self.assertIn(
            "gravityContact=experimentally unvalidated physical assumption",
            csv_rows[2]["claim_labels"],
        )
        self.assertEqual(
            removed_middle["mechanicsCertificate"]["schema"],
            "rad-sim.mechanics-energy-certificate.v1",
        )

    def test_vertical_load_energy_validation_compares_bench_heights_to_preview(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        comparison = compare_vertical_residual_spring_hinge_3d(
            config,
            LatticeState.uniform(config),
            ((0, 0), (0, 2)),
            ((0, 2),),
            z=0.25,
            fixed_cells=((0, 0),),
            external_z_load=-0.05,
            load_case=LoadCase(
                fixed_cells=((0, 0),),
                axial_stiffness=20.0,
                hinge_stiffness=1.0,
                lock_stiffness=250.0,
                maxiter=250,
            ),
        )
        intact_heights = {
            cell: float(comparison.intact_physical_height[cell])
            for cell in comparison.intact_solver_load_cells
        }
        removed_heights = {
            cell: float(comparison.removed_physical_height[cell])
            for cell in comparison.removed_solver_load_cells
        }
        validation = validate_vertical_load_energy_measurements(
            comparison,
            intact_heights,
            removed_heights,
            tolerance=1e-9,
        )
        self.assertEqual(
            validation["schema"],
            "rad-sim.vertical-load-energy-validation.v1",
        )
        self.assertTrue(validation["summary"]["passesTolerance"])
        self.assertEqual(validation["summary"]["missingMeasurementCount"], 0)
        self.assertAlmostEqual(
            validation["summary"]["maxAbsWorkOrContactError"],
            0.0,
        )
        self.assertIn(
            "measuredWorkResidualNat_zero_when_equal",
            validation["leanTheorems"],
        )
        exported = json.loads(export_vertical_load_energy_validation_json(validation))
        self.assertEqual(exported["schema"], validation["schema"])

        perturbed = dict(intact_heights)
        first_cell = comparison.intact_solver_load_cells[0]
        perturbed[first_cell] += 0.02
        perturbed_validation = validate_vertical_load_energy_measurements(
            comparison,
            perturbed,
            removed_heights,
            tolerance=1e-9,
        )
        self.assertFalse(perturbed_validation["summary"]["passesTolerance"])
        self.assertGreater(
            perturbed_validation["summary"]["maxAbsWorkOrContactError"],
            0.0,
        )
        self.assertGreater(
            abs(perturbed_validation["intact"]["errors"]["signedLoadWork"]),
            0.0,
        )

    def test_vertical_load_energy_measurement_template_roundtrips_to_report(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        load_case = LoadCase(
            fixed_cells=((0, 0),),
            axial_stiffness=20.0,
            hinge_stiffness=1.0,
            lock_stiffness=250.0,
            maxiter=250,
        )
        template = vertical_load_energy_measurement_template(
            config,
            load_case=load_case,
            tolerance=1e-9,
        )
        self.assertEqual(
            template["schema"],
            "rad-sim.vertical-load-energy-measurement-results.v1",
        )
        self.assertEqual(len(template["scenarios"]), 3)
        self.assertIn("predictedHeight", template["scenarios"][0]["measurements"]["intact"][0])
        exported_template = json.loads(
            export_vertical_load_energy_measurement_template_json(
                config,
                load_case=load_case,
                tolerance=1e-9,
            )
        )
        self.assertEqual(exported_template["schema"], template["schema"])

        filled = json.loads(json.dumps(template))
        for scenario in filled["scenarios"]:
            for side in ("intact", "removed"):
                for measurement in scenario["measurements"][side]:
                    measurement["measuredHeight"] = measurement["predictedHeight"]
        parsed = vertical_load_energy_measurement_results_from_json(json.dumps(filled))
        report = compare_vertical_load_energy_measurement_results(
            config,
            parsed,
            load_case=load_case,
        )
        self.assertEqual(
            report["schema"],
            "rad-sim.vertical-load-energy-comparison-report.v1",
        )
        self.assertTrue(report["summary"]["allScenariosPassTolerance"])
        self.assertEqual(report["summary"]["missingMeasurementCount"], 0)
        self.assertAlmostEqual(report["summary"]["maxAbsWorkOrContactError"], 0.0)
        exported_report = json.loads(
            export_vertical_load_energy_comparison_report_json(report)
        )
        self.assertEqual(exported_report["schema"], report["schema"])
        calibration_execution = {
            "schema": "rad-sim.calibration-bench-execution-validation.v1",
            "summary": {
                "executionValidationPass": True,
                "fitStepCount": 7,
                "holdoutStepCount": 7,
            },
        }
        readiness = physical_validation_readiness_report(
            config,
            calibration_execution,
            report,
        )
        self.assertEqual(
            readiness["schema"],
            "rad-sim.physical-validation-readiness.v1",
        )
        self.assertTrue(readiness["summary"]["physicalValidationReady"])
        self.assertEqual(readiness["summary"]["missingEvidence"], [])
        self.assertGreater(readiness["evidence"]["loadProxyTermCount"], 0)
        self.assertGreater(readiness["evidence"]["contactProxyTermCount"], 0)
        self.assertEqual(
            readiness["formalization"]["targetId"],
            "physical_validation_readiness_gate",
        )
        readiness_csv = export_physical_validation_readiness_csv(readiness)
        self.assertIn("physical_validation_ready", readiness_csv)
        self.assertIn("ready-for-physical-claim-review", readiness_csv)
        exported_readiness = json.loads(
            export_physical_validation_readiness_json(
                config,
                calibration_execution,
                report,
            )
        )
        self.assertEqual(exported_readiness["summary"], readiness["summary"])
        missing_readiness = physical_validation_readiness_report(config)
        self.assertFalse(missing_readiness["summary"]["physicalValidationReady"])
        self.assertIn(
            "verticalLoadComparisonReport",
            missing_readiness["summary"]["missingEvidence"],
        )

        filled["scenarios"][0]["measurements"]["intact"][0]["measuredHeight"] += 0.02
        perturbed = compare_vertical_load_energy_measurement_results(
            config,
            filled,
            load_case=load_case,
        )
        self.assertFalse(perturbed["summary"]["allScenariosPassTolerance"])
        self.assertGreater(perturbed["summary"]["maxAbsWorkOrContactError"], 0.0)

    def test_vertical_load_energy_experiment_protocol_exports_bench_steps(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        load_case = LoadCase(
            fixed_cells=((0, 0),),
            axial_stiffness=20.0,
            hinge_stiffness=1.0,
            lock_stiffness=250.0,
            maxiter=250,
        )
        protocol = vertical_load_energy_experiment_protocol(
            config,
            load_case=load_case,
            repeat_count=4,
            tolerance=1e-9,
        )
        self.assertEqual(
            protocol["schema"],
            "rad-sim.vertical-load-energy-experiment-protocol.v1",
        )
        self.assertEqual(
            protocol["measurementSchema"],
            "rad-sim.vertical-load-energy-measurement-results.v1",
        )
        self.assertEqual(protocol["repeatCount"], 4)
        self.assertEqual(len(protocol["steps"]), 3)
        first = protocol["steps"][0]
        self.assertEqual(first["fixture"]["fixedCells"], [{"row": 0, "col": 0}])
        self.assertIn("measuredHeight", first["requiredMeasurements"]["fields"])
        self.assertIn("height gauge", protocol["requiredInstruments"][0])
        self.assertIn("compare_vertical_load_energy_measurement_results", protocol["outputs"]["comparisonFunction"])
        self.assertEqual(
            protocol["claimLabels"]["benchMeasurement"],
            "experimentally unvalidated physical assumption",
        )
        exported = json.loads(
            export_vertical_load_energy_experiment_protocol_json(
                config,
                load_case=load_case,
                repeat_count=4,
                tolerance=1e-9,
            )
        )
        self.assertEqual(exported["schema"], protocol["schema"])
        self.assertEqual(exported["steps"][0]["id"], first["id"])

    def test_vertical_load_bench_packet_bundles_preview_protocol_and_template(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        load_case = LoadCase(
            fixed_cells=((0, 0),),
            axial_stiffness=20.0,
            hinge_stiffness=1.0,
            lock_stiffness=250.0,
            maxiter=250,
        )
        packet = vertical_load_bench_packet(
            config,
            load_case=load_case,
            repeat_count=2,
            tolerance=1e-9,
        )
        self.assertEqual(packet["schema"], "rad-sim.vertical-load-bench-packet.v1")
        self.assertEqual(
            packet["schemas"]["protocol"],
            "rad-sim.vertical-load-energy-experiment-protocol.v1",
        )
        self.assertEqual(
            packet["previewReport"]["schema"],
            "rad-sim.vertical-load-physical-preview-report.v1",
        )
        self.assertEqual(
            packet["measurementTemplate"]["schema"],
            "rad-sim.vertical-load-energy-measurement-results.v1",
        )
        self.assertEqual(
            packet["experimentProtocol"]["schema"],
            "rad-sim.vertical-load-energy-experiment-protocol.v1",
        )
        self.assertEqual(packet["parameters"]["repeatCount"], 2)
        self.assertEqual(len(packet["mechanicsCertificateSummaries"]), 3)
        self.assertIn(
            "compare_vertical_load_energy_measurement_results",
            packet["comparisonInstructions"]["compareFunction"],
        )
        self.assertEqual(
            packet["claimLabels"]["benchMeasurement"],
            "experimentally unvalidated physical assumption",
        )
        exported = json.loads(
            export_vertical_load_bench_packet_json(
                config,
                load_case=load_case,
                repeat_count=2,
                tolerance=1e-9,
            )
        )
        self.assertEqual(exported["schema"], packet["schema"])
        self.assertEqual(exported["previewReport"]["summary"]["scenarioCount"], 3)

    def test_vertical_load_bench_packet_applies_measured_hardware_profile(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            cell_size=2.0,
            backlash=0.05,
            pin_radius=0.1,
            hole_radius=0.12,
        )
        profile = RADHardwareProfile(
            name="bench-measured-v2",
            side_length_mm=40.0,
            backlash_mm=4.0,
            pin_radius_mm=5.0,
            hole_radius_mm=5.8,
        )
        packet = vertical_load_bench_packet(config, hardware_profile=profile)
        self.assertAlmostEqual(packet["grid"]["pinHoleClearance"], 0.04)
        self.assertEqual(
            packet["unitScaleMetadata"]["hardwareProfile"]["name"],
            "bench-measured-v2",
        )
        self.assertTrue(packet["unitScaleMetadata"]["profileApplication"]["applied"])

    def test_vertical_load_bench_packet_writer_creates_lab_artifacts(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        load_case = LoadCase(
            fixed_cells=((0, 0),),
            axial_stiffness=20.0,
            hinge_stiffness=1.0,
            lock_stiffness=250.0,
            maxiter=250,
        )
        with tempfile.TemporaryDirectory() as tmp:
            written = write_vertical_load_bench_packet_artifacts(
                tmp,
                config=config,
                load_case=load_case,
                repeat_count=2,
                tolerance=1e-9,
            )
            expected = {
                "packet",
                "preview_json",
                "preview_csv",
                "measurement_template",
                "experiment_protocol",
                "unit_scale_metadata",
                "hardware_profile",
                "readme",
            }
            self.assertEqual(set(written), expected)
            for path in written.values():
                self.assertTrue(Path(path).exists(), path)

            packet = json.loads(Path(written["packet"]).read_text(encoding="utf-8"))
            preview = json.loads(Path(written["preview_json"]).read_text(encoding="utf-8"))
            template = json.loads(
                Path(written["measurement_template"]).read_text(encoding="utf-8")
            )
            protocol = json.loads(
                Path(written["experiment_protocol"]).read_text(encoding="utf-8")
            )
            unit_scale = json.loads(
                Path(written["unit_scale_metadata"]).read_text(encoding="utf-8")
            )
            hardware_profile = json.loads(
                Path(written["hardware_profile"]).read_text(encoding="utf-8")
            )
            csv_text = Path(written["preview_csv"]).read_text(encoding="utf-8")
            readme = Path(written["readme"]).read_text(encoding="utf-8")

            self.assertEqual(packet["schema"], "rad-sim.vertical-load-bench-packet.v1")
            self.assertEqual(
                preview["schema"], "rad-sim.vertical-load-physical-preview-report.v1"
            )
            self.assertEqual(
                template["schema"], "rad-sim.vertical-load-energy-measurement-results.v1"
            )
            self.assertEqual(
                protocol["schema"], "rad-sim.vertical-load-energy-experiment-protocol.v1"
            )
            self.assertEqual(
                unit_scale["schema"], "rad-sim.physical-unit-scale-metadata.v1"
            )
            self.assertEqual(hardware_profile["schema"], HARDWARE_PROFILE_SCHEMA)
            self.assertEqual(
                packet["unitScaleMetadata"]["formalizationTarget"]["id"],
                "measurement_unit_scale_invariants",
            )
            self.assertIn("scenario", csv_text.splitlines()[0])
            self.assertIn("claim_labels", csv_text.splitlines()[0])
            self.assertIn("unit scale metadata", readme)
            self.assertIn("Do not copy predictedHeight", readme)

    def test_vertical_load_measurement_comparison_writer_creates_report_artifacts(self):
        config = LatticeConfig(
            rows=1,
            cols=4,
            pin_radius=0.10,
            hole_radius=0.12,
            z_coupling_gain=0.5,
        )
        load_case = LoadCase(
            fixed_cells=((0, 0),),
            axial_stiffness=20.0,
            hinge_stiffness=1.0,
            lock_stiffness=250.0,
            maxiter=250,
        )
        filled = vertical_load_energy_measurement_template(
            config,
            load_case=load_case,
            tolerance=1e-9,
        )
        for scenario in filled["scenarios"]:
            for side in ("intact", "removed"):
                for measurement in scenario["measurements"][side]:
                    measurement["measuredHeight"] = measurement["predictedHeight"]

        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "filled_measurements.json"
            input_path.write_text(json.dumps(filled, indent=2), encoding="utf-8")
            out_dir = Path(tmp) / "comparison"
            written = write_vertical_load_energy_comparison_artifacts(
                input_path,
                out_dir,
                config=config,
                load_case=load_case,
            )
            self.assertEqual(
                set(written),
                {"report", "summary_csv", "unit_scale_metadata", "hardware_profile", "readme"},
            )
            for path in written.values():
                self.assertTrue(Path(path).exists(), path)

            report = json.loads(Path(written["report"]).read_text(encoding="utf-8"))
            unit_scale = json.loads(
                Path(written["unit_scale_metadata"]).read_text(encoding="utf-8")
            )
            hardware_profile = json.loads(
                Path(written["hardware_profile"]).read_text(encoding="utf-8")
            )
            summary_csv = Path(written["summary_csv"]).read_text(encoding="utf-8")
            readme = Path(written["readme"]).read_text(encoding="utf-8")
            self.assertEqual(
                report["schema"],
                "rad-sim.vertical-load-energy-comparison-report.v1",
            )
            self.assertTrue(report["summary"]["allScenariosPassTolerance"])
            self.assertEqual(report["summary"]["missingMeasurementCount"], 0)
            self.assertEqual(
                report["unitScaleMetadata"]["schema"],
                "rad-sim.physical-unit-scale-metadata.v1",
            )
            self.assertEqual(
                unit_scale["formalizationTarget"]["id"],
                "measurement_unit_scale_invariants",
            )
            self.assertEqual(hardware_profile["schema"], HARDWARE_PROFILE_SCHEMA)
            rows = list(csv.DictReader(io.StringIO(summary_csv)))
            self.assertEqual(len(rows), 3)
            self.assertIn("passes_tolerance", rows[0])
            self.assertEqual(rows[0]["passes_tolerance"], "True")
            self.assertIn("unit scale metadata", readme)
            self.assertIn("does not prove calibrated gravity", readme)

    def test_single_cell_characterization_reports_reach_and_dieoff(self):
        config = LatticeConfig(rows=5, cols=5, backlash=0.02, z_coupling_gain=0.35)
        response = characterize_single_cell(config, (2, 2), alpha=-0.3, z=0.4)
        self.assertGreater(response.alpha_reach, 1)
        self.assertGreater(response.z_reach, 1)
        self.assertGreaterEqual(response.effective_alpha_die_off, 1)
        self.assertGreaterEqual(response.effective_z_die_off, 1)
        self.assertAlmostEqual(response.height_delta[2, 2], 0.4 - 0.65 * -0.3)

    def test_pair_characterization_reports_superposition_error(self):
        config = LatticeConfig(rows=5, cols=5, backlash=0.15, z_coupling_gain=0.35)
        pair = characterize_pair(
            config,
            SourceCommand((2, 1), alpha=-0.25, z=0.2),
            SourceCommand((2, 3), alpha=-0.25, z=0.2),
        )
        self.assertLess(pair.alpha_superposition_error, 1e-9)
        self.assertLess(pair.height_superposition_error, 1e-9)
        self.assertGreater(pair.combined.alpha_reach, pair.first.alpha_reach)

    def test_cluster_characterization_combines_multiple_sources(self):
        config = LatticeConfig(rows=5, cols=5, backlash=0.02, z_coupling_gain=0.35)
        cluster = characterize_cluster(
            config,
            [
                SourceCommand((1, 1), alpha=-0.25),
                SourceCommand((3, 3), z=0.35),
                SourceCommand((1, 3), alpha=0.2, z=-0.2),
            ],
        )
        self.assertEqual(len(cluster.commands), 3)
        self.assertGreater(cluster.alpha_reach, 1)
        self.assertGreater(cluster.z_reach, 1)
        self.assertGreater(cluster.max_abs_height_delta, 0)

    def test_calibration_experiment_protocol_covers_single_pair_cluster_and_lock(self):
        config = LatticeConfig(rows=5, cols=5, z_coupling_gain=0.35)
        profile = RADHardwareProfile(name="bench-v1", pin_radius_mm=6.3, hole_radius_mm=7.875)
        protocol = build_calibration_experiment_protocol(
            config,
            center_cell=(2, 2),
            hardware_profile=profile,
        )
        self.assertEqual(protocol.schema, "rad-sim.calibration-experiment-protocol.v1")
        self.assertEqual(protocol.hardware_profile_name, "bench-v1")
        self.assertEqual(protocol.center_cell, (2, 2))
        self.assertEqual(len(protocol.steps), 7)
        self.assertEqual(
            {step.scope for step in protocol.steps},
            {"single", "pair", "cluster", "lock"},
        )
        residual = next(step for step in protocol.steps if step.id == "pair_z_residual")
        self.assertEqual(residual.scope, "pair")
        self.assertEqual(len(residual.commands), 1)
        self.assertGreaterEqual(len(residual.observation_cells), 2)
        locked = next(step for step in protocol.steps if step.id == "locked_cell_control")
        self.assertEqual(locked.locked_cells, ((2, 2),))
        exported = json.loads(export_calibration_experiment_protocol_json(protocol))
        self.assertEqual(exported["schema"], protocol.schema)
        self.assertEqual(exported["hardwareProfile"], "bench-v1")
        self.assertEqual(exported["steps"][0]["id"], "single_alpha_contract")
        self.assertIn("pin_hole_slip_mm", exported["measurementFields"])

    def test_calibration_experiment_protocol_runs_kinematic_expectations(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.35)
        protocol = build_calibration_experiment_protocol(config, center_cell=(1, 1))
        simulations = run_calibration_experiment_protocol(config, protocol)
        self.assertEqual(len(simulations), len(protocol.steps))
        single_z = next(result for result in simulations if result.step_id == "single_z_lift")
        self.assertGreater(single_z.z_reach, 1)
        self.assertGreater(single_z.max_abs_height_delta, 0.0)
        locked = next(result for result in simulations if result.step_id == "locked_cell_control")
        self.assertGreaterEqual(locked.max_abs_height_delta, 0.0)
        locked_step = next(step for step in protocol.steps if step.id == "locked_cell_control")
        locked_response = characterize_response(
            config,
            locked_step.commands,
            locked_cells=locked_step.locked_cells,
        )
        row, col = locked_step.locked_cells[0]
        self.assertEqual(locked_response.alpha_delta[row, col], 0.0)
        self.assertEqual(locked_response.height_delta[row, col], 0.0)

    def test_response_atlas_exports_single_pair_cluster_and_lock_metrics(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.35)
        atlas = build_response_atlas(config, center_cell=(1, 1))
        payload = atlas.to_dict()

        self.assertEqual(payload["schema"], "rad-sim.response-atlas.v1")
        self.assertEqual(payload["summary"]["entryCount"], 7)
        self.assertEqual(payload["summary"]["scopeCounts"]["single"], 3)
        self.assertEqual(payload["summary"]["scopeCounts"]["pair"], 2)
        self.assertEqual(payload["summary"]["scopeCounts"]["cluster"], 1)
        self.assertEqual(payload["summary"]["scopeCounts"]["lock"], 1)
        single_z = next(entry for entry in payload["entries"] if entry["stepId"] == "single_z_lift")
        self.assertEqual(len(single_z["observationCells"]), 2)
        self.assertGreater(abs(single_z["observationCells"][1]["heightDelta"]), 0.0)
        pair = next(entry for entry in payload["entries"] if entry["stepId"] == "pair_superposition")
        self.assertEqual(len(pair["commands"]), 2)
        self.assertGreaterEqual(pair["heightSuperpositionError"], 0.0)
        locked = next(entry for entry in payload["entries"] if entry["stepId"] == "locked_cell_control")
        self.assertAlmostEqual(locked["observationCells"][0]["alphaDelta"], 0.0)
        self.assertAlmostEqual(locked["observationCells"][0]["heightDelta"], 0.0)
        exported = json.loads(export_response_atlas_json(atlas))
        self.assertEqual(exported["summary"], payload["summary"])

    def test_response_atlas_can_attach_physical_metrics(self):
        config = LatticeConfig(rows=2, cols=2, z_coupling_gain=0.0)
        atlas = build_response_atlas(
            config,
            center_cell=(0, 0),
            physical=True,
            load_case=LoadCase(lock_stiffness=350.0, maxiter=250),
        )
        payload = atlas.to_dict()

        self.assertTrue(payload["physical"])
        self.assertEqual(payload["summary"]["physicalEntryCount"], len(payload["entries"]))
        self.assertGreater(payload["summary"]["physicalSuccessCount"], 0)
        physical_entry = next(entry for entry in payload["entries"] if entry["physicalSuccess"])
        self.assertIsNotNone(physical_entry["physicalHeightRmsError"])
        self.assertIsNotNone(physical_entry["physicalCenterRmsError"])
        self.assertGreaterEqual(physical_entry["physicalEnergy"], 0.0)

    def test_response_atlas_sweep_exports_backlash_and_clearance_trends(self):
        config = LatticeConfig(
            rows=3,
            cols=3,
            backlash=0.02,
            z_coupling_gain=0.35,
            pin_radius=0.18,
            hole_radius=0.20,
        )
        sweep = sweep_response_atlas_parameters(
            config,
            backlash_values=(0.02, 0.14),
            clearance_values=(0.02, 0.12),
            center_cell=(1, 1),
        )
        payload = sweep.to_dict()

        self.assertEqual(payload["schema"], "rad-sim.response-atlas-sweep.v1")
        self.assertEqual(payload["summary"]["sampleCount"], 4)
        self.assertEqual(payload["parameters"]["backlashValues"], [0.02, 0.14])
        self.assertEqual(payload["parameters"]["pinHoleClearanceValues"], [0.02, 0.12])
        self.assertEqual(len(payload["samples"]), 4)
        self.assertGreater(payload["summary"]["maxObservedNeighborZResidual"], 0.0)
        by_clearance = payload["trends"]["byPinHoleClearance"]
        self.assertGreater(
            by_clearance[0]["maxObservedNeighborZResidual"],
            by_clearance[1]["maxObservedNeighborZResidual"],
        )
        sensitivity = payload["sensitivity"]
        self.assertEqual(
            sensitivity["method"],
            "endpoint finite difference over each parameter trend",
        )
        self.assertLess(
            sensitivity["metrics"]["pinHoleClearance"]["maxObservedNeighborZResidual"],
            0.0,
        )
        self.assertIsNotNone(sensitivity["dominant"])
        laws = payload["operatorLawCandidates"]
        self.assertEqual(
            laws["method"],
            "adjacent monotonicity over sampled parameter trend plus endpoint sensitivity",
        )
        clearance_law = next(
            law
            for law in laws["laws"]
            if law["parameter"] == "pinHoleClearance"
            and law["metric"] == "maxObservedNeighborZResidual"
        )
        self.assertEqual(clearance_law["monotonicity"], "decreasing")
        self.assertLess(clearance_law["slope"], 0.0)
        self.assertTrue(clearance_law["supportedBySweep"])
        self.assertEqual(clearance_law["status"], "simulator-diagnostic")
        first_sample = payload["samples"][0]
        self.assertEqual(first_sample["atlas"]["schema"], "rad-sim.response-atlas.v1")
        self.assertAlmostEqual(first_sample["settings"]["holeRadius"], 0.20)

        exported = json.loads(export_response_atlas_sweep_json(sweep))
        self.assertEqual(exported["summary"], payload["summary"])
        self.assertEqual(exported["sensitivity"], payload["sensitivity"])
        self.assertEqual(exported["operatorLawCandidates"], payload["operatorLawCandidates"])

    def test_calibration_results_template_roundtrips_and_compares_to_simulation(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.35)
        protocol = build_calibration_experiment_protocol(
            config,
            center_cell=(1, 1),
            repeat_count=2,
        )
        notebook = calibration_bench_notebook(
            config,
            protocol,
            fit_dataset_id="fit-test",
            holdout_dataset_id="holdout-test",
            profile_id="profile-test-v1",
        )
        self.assertEqual(notebook["schema"], "rad-sim.calibration-bench-notebook.v1")
        self.assertEqual(notebook["datasetPlan"]["fit"]["role"], "fit")
        self.assertEqual(notebook["datasetPlan"]["holdout"]["role"], "holdout")
        self.assertEqual(notebook["datasetPlan"]["holdout"]["profileId"], "profile-test-v1")
        self.assertEqual(notebook["protocol"]["stepCount"], len(protocol.steps))
        self.assertEqual(len(notebook["scenarios"]), len(protocol.steps))
        self.assertEqual(len(notebook["instruments"]), 6)
        self.assertIn("profileFrozenAt", notebook["measurementColumns"])
        self.assertTrue(
            any(output["id"] == "holdout-validation-csv" for output in notebook["outputs"])
        )
        self.assertTrue(
            any(output["id"] == "bench-notebook-csv" for output in notebook["outputs"])
        )
        self.assertTrue(
            any(scenario["stepId"] == "pair_z_residual" for scenario in notebook["scenarios"])
        )
        exported_notebook = json.loads(
            export_calibration_bench_notebook_json(
                config,
                protocol,
                fit_dataset_id="fit-test",
                holdout_dataset_id="holdout-test",
                profile_id="profile-test-v1",
            )
        )
        self.assertEqual(exported_notebook["datasetPlan"], notebook["datasetPlan"])
        notebook_rows = list(
            csv.DictReader(io.StringIO(export_calibration_bench_notebook_csv(notebook)))
        )
        self.assertEqual(len(notebook_rows), len(protocol.steps))
        self.assertEqual(notebook_rows[3]["step_id"], "pair_z_residual")
        self.assertEqual(notebook_rows[3]["fit_dataset_id"], "fit-test")
        self.assertEqual(notebook_rows[3]["holdout_dataset_id"], "holdout-test")
        self.assertIn("pin_hole_slip_mm", notebook_rows[3]["measurement_columns"])
        packet = calibration_bench_packet(
            config,
            protocol,
            fit_dataset_id="fit-test",
            holdout_dataset_id="holdout-test",
            profile_id="profile-test-v1",
            profile_frozen_at="2026-08-09T10:30:00Z",
        )
        self.assertEqual(packet["schema"], "rad-sim.calibration-bench-packet.v1")
        self.assertEqual(
            packet["schemas"]["notebook"],
            "rad-sim.calibration-bench-notebook.v1",
        )
        self.assertEqual(
            packet["benchNotebook"]["datasetPlan"]["fit"]["datasetId"],
            notebook["datasetPlan"]["fit"]["datasetId"],
        )
        self.assertEqual(packet["benchNotebook"]["datasetPlan"]["holdout"]["role"], "holdout")
        self.assertEqual(
            packet["fitResultsTemplate"]["provenance"]["datasetRole"],
            "fit",
        )
        self.assertEqual(
            packet["holdoutResultsTemplate"]["provenance"]["datasetRole"],
            "holdout",
        )
        self.assertEqual(
            packet["holdoutResultsTemplate"]["provenance"]["profileId"],
            "profile-test-v1",
        )
        self.assertIn("pair_z_residual", packet["benchNotebookCsv"])
        self.assertTrue(
            any(item["id"] == "holdout-template" for item in packet["artifactManifest"])
        )
        exported_packet = json.loads(
            export_calibration_bench_packet_json(
                config,
                protocol,
                fit_dataset_id="fit-test",
                holdout_dataset_id="holdout-test",
                profile_id="profile-test-v1",
                profile_frozen_at="2026-08-09T10:30:00Z",
            )
        )
        self.assertEqual(exported_packet["filenames"], packet["filenames"])
        template = calibration_experiment_results_template(protocol)
        self.assertEqual(
            len(template.steps),
            sum(step.repeat_count for step in protocol.steps),
        )
        template_json = json.loads(export_calibration_experiment_results_template_json(protocol))
        self.assertEqual(template_json["schema"], "rad-sim.calibration-experiment-results.v1")
        step_result = next(
            step
            for step in template_json["steps"]
            if step["stepId"] == "single_z_lift" and step["repeatIndex"] == 1
        )
        protocol_step = next(step for step in protocol.steps if step.id == "single_z_lift")
        simulated = characterize_response(
            config,
            protocol_step.commands,
            locked_cells=protocol_step.locked_cells,
        )
        for cell in step_result["cells"]:
            row = cell["row"]
            col = cell["col"]
            cell["alphaDelta"] = float(simulated.alpha_delta[row, col])
            cell["heightDelta"] = float(simulated.height_delta[row, col])
            cell["pinHoleSlipMm"] = 0.12
            cell["actuatorForceN"] = 3.4
        measurements = calibration_experiment_measurements_from_json(json.dumps(template_json))
        comparisons = compare_calibration_experiment_measurements(config, protocol, measurements)
        comparison = next(item for item in comparisons if item.step_id == "single_z_lift" and item.repeat_index == 1)
        self.assertEqual(comparison.missing_observation_count, 0)
        self.assertAlmostEqual(comparison.alpha_rmse or 0.0, 0.0)
        self.assertAlmostEqual(comparison.height_rmse or 0.0, 0.0)
        self.assertAlmostEqual(comparison.mean_signed_height_error or 0.0, 0.0)
        self.assertAlmostEqual(comparison.mean_abs_height_error or 0.0, 0.0)
        self.assertAlmostEqual(comparison.mean_pin_hole_slip_mm or 0.0, 0.12)
        self.assertAlmostEqual(comparison.mean_actuator_force_n or 0.0, 3.4)
        exported = json.loads(export_calibration_experiment_comparison_json(comparisons))
        self.assertEqual(exported["schema"], "rad-sim.calibration-experiment-comparison.v1")
        report = calibration_experiment_comparison_report(config, protocol, measurements)
        self.assertEqual(report["schema"], "rad-sim.calibration-comparison-report.v1")
        self.assertEqual(
            report["parameterEstimates"]["schema"],
            "rad-sim.calibration-parameter-estimates.v1",
        )
        self.assertEqual(report["comparison"]["fit"]["height"]["sampleCount"], len(step_result["cells"]))
        self.assertAlmostEqual(report["comparison"]["fit"]["height"]["suggestedGain"], 1.0)
        self.assertAlmostEqual(report["comparison"]["fit"]["height"]["suggestedBias"], 0.0)
        self.assertEqual(
            report["parameterEstimates"]["estimates"]["alphaResponse"]["status"],
            "estimated-from-results",
        )
        self.assertEqual(
            report["parameterEstimates"]["estimates"]["zCouplingGain"]["sampleCount"],
            1,
        )
        self.assertIsNotNone(
            report["parameterEstimates"]["estimates"]["zCouplingGain"]["estimate"]
        )
        self.assertEqual(
            report["parameterEstimates"]["estimates"]["verticalFreePlay"]["sampleCount"],
            len(step_result["cells"]),
        )
        self.assertEqual(
            report["parameterEstimates"]["estimates"]["backlash"]["status"],
            "not-identifiable-from-current-results",
        )
        self.assertGreater(
            report["parameterEstimates"]["estimates"]["contactProxy"]["forceHeightSlopeNPerModelUnit"],
            0.0,
        )
        profile = report["modelProfile"]
        self.assertEqual(profile["schema"], "rad-sim.calibration-model-profile.v1")
        self.assertEqual(
            profile["sourceReportSchema"],
            "rad-sim.calibration-comparison-report.v1",
        )
        z_update = next(
            update
            for update in profile["recommendedUpdates"]
            if update["name"] == "zCouplingGain"
        )
        self.assertTrue(z_update["safeToApply"])
        self.assertEqual(z_update["configField"], "z_coupling_gain")
        self.assertEqual(z_update["sampleCount"], 1)
        alpha_update = next(
            update
            for update in profile["recommendedUpdates"]
            if update["name"] == "alphaResponseGain"
        )
        self.assertFalse(alpha_update["safeToApply"])
        fitted_config, application = apply_calibration_model_profile(config, profile)
        self.assertEqual(
            application["schema"],
            "rad-sim.calibration-model-profile-application.v1",
        )
        self.assertEqual(
            application["appliedUpdates"][0]["configField"],
            "z_coupling_gain",
        )
        self.assertAlmostEqual(
            fitted_config.z_coupling_gain,
            z_update["proposed"],
        )
        exported_profile = json.loads(export_calibration_model_profile_json(config, report))
        self.assertEqual(exported_profile, calibration_model_profile_from_report(config, report))
        self.assertAlmostEqual(report["comparison"]["field"]["maxCombinedError"], 0.0)
        self.assertAlmostEqual(report["comparison"]["fitResidualField"]["maxCombinedError"], 0.0)
        self.assertIsNotNone(report["summary"]["worstCell"])
        self.assertGreater(len(report["summary"]["topCells"]), 0)
        self.assertEqual(report["summary"]["topCells"][0], report["summary"]["worstCell"])
        perturbed_json = json.loads(json.dumps(template_json))
        perturbed_step = next(
            step
            for step in perturbed_json["steps"]
            if step["stepId"] == "single_z_lift" and step["repeatIndex"] == 1
        )
        perturbed_step["cells"][0]["heightDelta"] += 0.05
        perturbed = calibration_experiment_measurements_from_json(json.dumps(perturbed_json))
        perturbed_report = calibration_experiment_comparison_report(config, protocol, perturbed)
        self.assertGreater(perturbed_report["comparison"]["field"]["maxCombinedError"], 0.0)
        self.assertGreater(perturbed_report["comparison"]["fit"]["height"]["rmsRawError"], 0.0)
        self.assertGreater(len(perturbed_report["summary"]["topCells"]), 0)
        self.assertGreater(len(perturbed_report["summary"]["fitResidualTopCells"]), 0)
        self.assertTrue(
            any(
                any(value > 0 for value in row)
                for row in perturbed_report["comparison"]["fitResidualField"]["sampleCount"]
            )
        )
        exported_report = json.loads(
            export_calibration_experiment_comparison_report_json(config, protocol, perturbed)
        )
        self.assertEqual(exported_report["schema"], report["schema"])

        low_z_config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.05)
        high_z_config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.5)
        profile_protocol = build_calibration_experiment_protocol(
            low_z_config,
            center_cell=(1, 1),
            repeat_count=1,
        )
        profile_template = calibration_experiment_results_template(profile_protocol).to_dict()
        self.assertEqual(
            profile_template["provenance"]["schema"],
            "rad-sim.calibration-dataset-provenance.v1",
        )
        self.assertEqual(profile_template["provenance"]["datasetRole"], "unassigned")
        profile_step = next(
            step
            for step in profile_template["steps"]
            if step["stepId"] == "single_z_lift"
        )
        profile_protocol_step = next(
            step for step in profile_protocol.steps if step.id == "single_z_lift"
        )
        high_z_response = characterize_response(
            high_z_config,
            profile_protocol_step.commands,
            locked_cells=profile_protocol_step.locked_cells,
        )
        for cell in profile_step["cells"]:
            row = cell["row"]
            col = cell["col"]
            cell["alphaDelta"] = float(high_z_response.alpha_delta[row, col])
            cell["heightDelta"] = float(high_z_response.height_delta[row, col])
        profile_measurements = calibration_experiment_measurements_from_json(
            json.dumps(profile_template)
        )
        residual_comparison = calibration_model_profile_residual_comparison(
            low_z_config,
            profile_protocol,
            profile_measurements,
        )
        self.assertEqual(
            residual_comparison["schema"],
            "rad-sim.calibration-model-profile-residual-comparison.v1",
        )
        self.assertTrue(residual_comparison["application"]["appliedUpdates"])
        self.assertLess(
            residual_comparison["after"]["metrics"]["residualScore"],
            residual_comparison["before"]["metrics"]["residualScore"],
        )
        exported_residual = json.loads(
            export_calibration_model_profile_residual_comparison_json(
                low_z_config,
                profile_protocol,
                profile_measurements,
            )
        )
        self.assertEqual(exported_residual["delta"], residual_comparison["delta"])
        rejected = json.loads(json.dumps(residual_comparison))
        rejected["application"]["appliedUpdates"] = []
        rejected["delta"]["residualScore"] = 0.0
        selection = select_calibration_model_profile([rejected, residual_comparison])
        self.assertEqual(
            selection["schema"],
            "rad-sim.calibration-model-profile-selection.v1",
        )
        self.assertEqual(selection["candidateCount"], 2)
        self.assertEqual(selection["eligibleCount"], 1)
        self.assertEqual(selection["selectedIndex"], 1)
        exported_selection = json.loads(
            export_calibration_model_profile_selection_json([rejected, residual_comparison])
        )
        self.assertEqual(exported_selection["selectedIndex"], 1)
        holdout_config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.45)
        holdout_template = calibration_experiment_results_template(profile_protocol).to_dict()
        protocol_steps = {step.id: step for step in profile_protocol.steps}
        for measured_step in holdout_template["steps"]:
            protocol_step = protocol_steps[measured_step["stepId"]]
            holdout_response = characterize_response(
                holdout_config,
                protocol_step.commands,
                locked_cells=protocol_step.locked_cells,
            )
            for cell in measured_step["cells"]:
                row = cell["row"]
                col = cell["col"]
                cell["alphaDelta"] = float(holdout_response.alpha_delta[row, col])
                cell["heightDelta"] = float(holdout_response.height_delta[row, col])
        holdout_measurements = calibration_experiment_measurements_from_json(
            json.dumps(holdout_template)
        )
        holdout_validation = calibration_model_profile_holdout_validation(
            low_z_config,
            profile_protocol,
            profile_measurements,
            holdout_measurements,
        )
        self.assertEqual(
            holdout_validation["schema"],
            "rad-sim.calibration-model-profile-holdout-validation.v1",
        )
        self.assertTrue(holdout_validation["fitPass"])
        self.assertTrue(holdout_validation["holdoutPass"])
        self.assertTrue(holdout_validation["residualValidationPass"])
        self.assertFalse(holdout_validation["independentValidationPass"])
        self.assertEqual(
            holdout_validation["splitMetadata"]["schema"],
            "rad-sim.calibration-train-holdout-split.v1",
        )
        self.assertEqual(holdout_validation["splitMetadata"]["fitStepCount"], 7)
        self.assertEqual(holdout_validation["splitMetadata"]["holdoutStepCount"], 7)
        self.assertEqual(
            holdout_validation["splitMetadata"]["independenceStatus"],
            "requires-external-bench-protocol",
        )
        self.assertIn("fitDatasetId", holdout_validation["splitMetadata"]["missingEvidence"])
        self.assertGreater(len(holdout_validation["provenanceWarnings"]), 0)
        self.assertLess(
            holdout_validation["holdout"]["after"]["metrics"]["residualScore"],
            holdout_validation["holdout"]["before"]["metrics"]["residualScore"],
        )
        holdout_csv_rows = list(
            csv.DictReader(
                io.StringIO(
                    export_calibration_model_profile_holdout_validation_csv(
                        holdout_validation
                    )
                )
            )
        )
        self.assertEqual(len(holdout_csv_rows), 2)
        self.assertEqual(holdout_csv_rows[0]["dataset"], "fit")
        self.assertEqual(holdout_csv_rows[1]["dataset"], "holdout")
        self.assertEqual(holdout_csv_rows[1]["pass"], "true")
        self.assertEqual(holdout_csv_rows[1]["residual_validation_pass"], "true")
        self.assertEqual(holdout_csv_rows[1]["independent_validation_pass"], "false")
        self.assertEqual(
            holdout_csv_rows[1]["independence_status"],
            "requires-external-bench-protocol",
        )
        profile_with_id = json.loads(json.dumps(holdout_validation["modelProfile"]))
        profile_with_id["profileId"] = "synthetic-fit-profile-v1"
        fit_template_with_provenance = json.loads(json.dumps(profile_template))
        fit_template_with_provenance["provenance"] = {
            "schema": "rad-sim.calibration-dataset-provenance.v1",
            "datasetId": "fit-synthetic-z-050",
            "datasetRole": "fit",
            "sourceFileId": "fit-synthetic-z-050.json",
            "collectedAt": "2026-08-09T10:00:00Z",
            "operator": "unit-test",
        }
        holdout_template_with_provenance = json.loads(json.dumps(holdout_template))
        holdout_template_with_provenance["provenance"] = {
            "schema": "rad-sim.calibration-dataset-provenance.v1",
            "datasetId": "holdout-synthetic-z-045",
            "datasetRole": "holdout",
            "sourceFileId": "holdout-synthetic-z-045.json",
            "collectedAt": "2026-08-09T11:00:00Z",
            "operator": "unit-test",
            "profileId": "synthetic-fit-profile-v1",
            "profileFrozenAt": "2026-08-09T10:30:00Z",
        }
        documented_fit = calibration_experiment_measurements_from_json(
            json.dumps(fit_template_with_provenance)
        )
        documented_holdout = calibration_experiment_measurements_from_json(
            json.dumps(holdout_template_with_provenance)
        )
        documented_validation = calibration_model_profile_holdout_validation(
            low_z_config,
            profile_protocol,
            documented_fit,
            documented_holdout,
            profile_with_id,
        )
        self.assertTrue(documented_validation["residualValidationPass"])
        self.assertTrue(documented_validation["independentValidationPass"])
        self.assertEqual(documented_validation["provenanceWarnings"], [])
        self.assertEqual(
            documented_validation["splitMetadata"]["independenceStatus"],
            "documented-independent-holdout",
        )
        self.assertEqual(
            documented_validation["splitMetadata"]["fitDataset"]["datasetId"],
            "fit-synthetic-z-050",
        )
        documented_csv_rows = list(
            csv.DictReader(
                io.StringIO(
                    export_calibration_model_profile_holdout_validation_csv(
                        documented_validation
                    )
                )
            )
        )
        self.assertEqual(documented_csv_rows[1]["independent_validation_pass"], "true")
        self.assertEqual(
            documented_csv_rows[1]["independence_status"],
            "documented-independent-holdout",
        )
        execution_validation = calibration_bench_execution_validation(
            low_z_config,
            profile_protocol,
            documented_fit,
            documented_holdout,
            profile_with_id,
        )
        self.assertEqual(
            execution_validation["schema"],
            "rad-sim.calibration-bench-execution-validation.v1",
        )
        self.assertEqual(
            execution_validation["summary"]["status"],
            "documented-independent-validation",
        )
        self.assertTrue(execution_validation["summary"]["residualValidationPass"])
        self.assertTrue(execution_validation["summary"]["independentValidationPass"])
        self.assertTrue(execution_validation["summary"]["executionValidationPass"])
        self.assertEqual(execution_validation["summary"]["missingEvidenceCount"], 0)
        self.assertEqual(
            execution_validation["formalization"]["targetId"],
            "calibration_bench_executed_validation_gate",
        )
        self.assertEqual(
            execution_validation["modelProfile"]["profileId"],
            "synthetic-fit-profile-v1",
        )
        execution_csv_rows = list(
            csv.DictReader(
                io.StringIO(
                    export_calibration_bench_execution_validation_csv(
                        execution_validation
                    )
                )
            )
        )
        self.assertEqual(execution_csv_rows[0]["execution_validation_pass"], "true")
        self.assertEqual(execution_csv_rows[0]["fit_dataset_id"], "fit-synthetic-z-050")
        self.assertEqual(
            execution_csv_rows[0]["holdout_dataset_id"],
            "holdout-synthetic-z-045",
        )
        exported_execution = json.loads(
            export_calibration_bench_execution_validation_json(
                low_z_config,
                profile_protocol,
                documented_fit,
                documented_holdout,
                profile_with_id,
            )
        )
        self.assertEqual(exported_execution["summary"], execution_validation["summary"])
        with tempfile.TemporaryDirectory() as tmp:
            written_execution = write_calibration_bench_execution_validation_artifacts(
                tmp,
                low_z_config,
                profile_protocol,
                json.dumps(fit_template_with_provenance),
                json.dumps(holdout_template_with_provenance),
                profile=profile_with_id,
            )
            self.assertEqual(
                set(written_execution),
                {
                    "report",
                    "summary_csv",
                    "holdout_validation",
                    "holdout_csv",
                    "model_profile",
                    "fit_comparison_report",
                    "readme",
                },
            )
            written_report = json.loads(
                Path(written_execution["report"]).read_text(encoding="utf-8")
            )
            self.assertTrue(written_report["summary"]["executionValidationPass"])
            self.assertIn(
                "execution_validation_pass",
                Path(written_execution["summary_csv"]).read_text(encoding="utf-8"),
            )
        exported_holdout = json.loads(
            export_calibration_model_profile_holdout_validation_json(
                low_z_config,
                profile_protocol,
                profile_measurements,
                holdout_measurements,
            )
        )
        self.assertEqual(exported_holdout["holdoutPass"], holdout_validation["holdoutPass"])

    def test_calibration_bench_packet_writer_creates_lab_artifacts(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.35)
        with tempfile.TemporaryDirectory() as tmp:
            written = write_calibration_bench_packet_artifacts(
                tmp,
                config=config,
                center_cell=(1, 1),
                repeat_count=1,
                fit_dataset_id="fit-writer",
                holdout_dataset_id="holdout-writer",
                profile_id="profile-writer-v1",
                profile_frozen_at="2026-08-09T12:00:00Z",
            )
            self.assertEqual(
                set(written),
                {
                    "packet",
                    "notebook",
                    "notebook_csv",
                    "protocol",
                    "fit_template",
                    "holdout_template",
                    "readme",
                    "hardware_profile",
                },
            )
            for path in written.values():
                self.assertTrue(Path(path).exists())
            packet = json.loads(Path(written["packet"]).read_text(encoding="utf-8"))
            self.assertEqual(packet["schema"], "rad-sim.calibration-bench-packet.v1")
            self.assertEqual(packet["fitResultsTemplate"]["provenance"]["datasetId"], "fit-writer")
            self.assertEqual(
                packet["holdoutResultsTemplate"]["provenance"]["profileFrozenAt"],
                "2026-08-09T12:00:00Z",
            )
            notebook_csv = Path(written["notebook_csv"]).read_text(encoding="utf-8")
            self.assertIn("pair_z_residual", notebook_csv)
            self.assertIn("pin_hole_slip_mm", notebook_csv)
            readme = Path(written["readme"]).read_text(encoding="utf-8")
            self.assertIn("Keep the fit and holdout raw files separate", readme)
            profile = json.loads(Path(written["hardware_profile"]).read_text(encoding="utf-8"))
            self.assertEqual(profile["schema"], HARDWARE_PROFILE_SCHEMA)

    def test_characterization_respects_locked_cells(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        response = characterize_single_cell(
            config,
            (1, 1),
            alpha=-0.4,
            z=0.3,
            locked_cells=((1, 1),),
        )
        self.assertAlmostEqual(response.alpha_delta[1, 1], 0.0)
        self.assertAlmostEqual(response.height_delta[1, 1], 0.0)

    def test_response_matrix_builds_alpha_and_height_columns(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.35)
        matrix = build_response_matrix(config, actuator_cells=[(1, 1), (1, 2)])
        self.assertEqual(matrix.alpha.shape, (9, 4))
        self.assertEqual(matrix.height.shape, (9, 4))
        self.assertEqual(len(matrix.commands), 4)
        self.assertGreater(matrix.alpha_rank, 0)
        self.assertGreater(matrix.height_rank, 0)
        self.assertGreater(matrix.reachable_alpha_cells(), 1)
        self.assertGreater(matrix.reachable_height_cells(), 1)
        payload = matrix.to_dict()
        self.assertEqual(payload["schema"], "rad-sim.response-matrix.v1")
        self.assertEqual(payload["grid"], {"rows": 3, "cols": 3})
        self.assertEqual(len(payload["cellOrder"]), 9)
        self.assertEqual(payload["diagnostics"]["columnCount"], 4)
        self.assertEqual(payload["diagnostics"]["alphaRank"], matrix.alpha_rank)
        self.assertEqual(payload["diagnostics"]["reachableHeightCells"], matrix.reachable_height_cells())
        self.assertEqual(json.loads(export_response_matrix_json(matrix))["schema"], payload["schema"])

    def test_response_matrix_can_select_command_family(self):
        config = LatticeConfig(rows=3, cols=3)
        matrix = build_response_matrix(
            config,
            actuator_cells=[(1, 1)],
            include_alpha=False,
            include_z=True,
        )
        self.assertEqual(matrix.alpha.shape, (9, 1))
        self.assertEqual(matrix.commands[0].alpha, 0.0)
        self.assertGreater(matrix.commands[0].z, 0.0)
        self.assertAlmostEqual(matrix.alpha[4, 0], 0.0)
        self.assertGreater(matrix.height[4, 0], 0.0)

    def test_removed_topology_reachability_comparison_reports_losses(self):
        config = LatticeConfig(
            rows=1,
            cols=3,
            backlash=0.0,
            coupling_gain=1.0,
            z_coupling_gain=1.0,
        )
        comparison = compare_removed_topology_reachability(
            config,
            removed_cells=[(0, 1)],
            actuator_cells=[(0, 0)],
        )
        self.assertEqual(comparison.intact_topology["component_count"], 1)
        self.assertEqual(comparison.removed_topology["component_count"], 2)
        self.assertEqual(comparison.component_count_delta, 1)
        self.assertEqual(comparison.deleted_edge_delta, 2)
        self.assertGreater(comparison.alpha_reachable_cell_loss, 0)
        self.assertGreater(comparison.height_reachable_cell_loss, 0)
        payload = comparison.to_dict()
        self.assertEqual(
            payload["schema"],
            "rad-sim.removed-topology-reachability-comparison.v1",
        )
        self.assertEqual(payload["removed"]["topology"]["component_count"], 2)
        self.assertEqual(
            payload["summary"]["heightReachableCellLoss"],
            comparison.height_reachable_cell_loss,
        )
        exported = json.loads(
            export_removed_topology_reachability_comparison_json(comparison)
        )
        self.assertEqual(exported["schema"], payload["schema"])

    def test_topology_experiment_report_compares_discontinuity_scenarios(self):
        config = LatticeConfig(
            rows=1,
            cols=3,
            backlash=0.0,
            coupling_gain=1.0,
            z_coupling_gain=1.0,
            pin_radius=0.0,
            hole_radius=0.0,
        )
        target = np.zeros((1, 3), dtype=float)
        target[0, 2] = 0.2
        scenarios = (
            TopologyExperimentScenario(
                name="intact-left-source",
                commands=(SourceCommand((0, 0), alpha=-0.2, z=0.2),),
                actuator_cells=((0, 0),),
            ),
            TopologyExperimentScenario(
                name="locked-left-source",
                commands=(SourceCommand((0, 0), alpha=-0.2, z=0.2),),
                locked_cells=((0, 0),),
                actuator_cells=((0, 0),),
            ),
            TopologyExperimentScenario(
                name="removed-middle-left-source",
                commands=(SourceCommand((0, 0), alpha=-0.2, z=0.2),),
                removed_cells=((0, 1),),
                actuator_cells=((0, 0),),
            ),
        )
        report = build_topology_experiment_report(
            config,
            scenarios,
            target_height=target,
        )
        self.assertEqual(len(report.scenarios), 3)
        intact, locked, removed = report.scenarios
        self.assertGreaterEqual(intact.component_reachable_cells, 3)
        self.assertEqual(locked.response_matrix.reachable_height_cells(), 0)
        self.assertEqual(removed.topology["component_count"], 2)
        self.assertEqual(removed.component_reachable_cells, 1)
        self.assertEqual(removed.component_blocked_cells, 1)
        self.assertEqual(len(removed.component_summaries), 2)
        self.assertEqual(removed.component_summaries[0].cell_count, 1)
        self.assertEqual(removed.component_summaries[0].actuator_cell_count, 1)
        self.assertEqual(removed.component_summaries[1].actuator_cell_count, 0)
        self.assertEqual(removed.component_summaries[1].height_rank, 0)
        self.assertEqual(removed.component_summaries[1].blocked_height_cells, 1)
        self.assertGreater(removed.target_height_rms, 0.0)
        self.assertGreaterEqual(intact.positive_height_reachable_cells, 1)
        payload = report.to_dict()
        self.assertEqual(payload["schema"], "rad-sim.topology-experiment-report.v1")
        self.assertEqual(payload["summary"]["scenarioCount"], 3)
        self.assertGreaterEqual(payload["summary"]["maxComponentHeightRank"], 1)
        self.assertGreaterEqual(payload["summary"]["maxComponentBlockedHeightCells"], 1)
        self.assertEqual(payload["scenarioComparisons"][1]["componentCountDelta"], 1)
        self.assertLess(payload["scenarioComparisons"][1]["reachableHeightCellDelta"], 0)
        self.assertGreaterEqual(payload["summary"]["maxComponentBlockedCells"], 1)
        self.assertEqual(
            payload["scenarios"][2]["metrics"]["componentBlockedCells"],
            1,
        )
        self.assertEqual(
            payload["scenarios"][2]["metrics"]["componentResponseRank"][1]["blockedHeightCells"],
            1,
        )
        exported = json.loads(export_topology_experiment_report_json(report))
        self.assertEqual(exported["schema"], payload["schema"])

    def test_response_matrix_respects_locked_cells(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        matrix = build_response_matrix(
            config,
            actuator_cells=[(1, 1)],
            locked_cells=((1, 1),),
        )
        center_row = 1 * config.cols + 1
        np.testing.assert_allclose(matrix.alpha[center_row, :], 0.0)
        np.testing.assert_allclose(matrix.height[center_row, :], 0.0)

    def test_response_matrix_requires_at_least_one_command_family(self):
        with self.assertRaises(ValueError):
            build_response_matrix(
                LatticeConfig(rows=2, cols=2),
                include_alpha=False,
                include_z=False,
            )

    def test_response_decay_profile_fits_shell_decay(self):
        config = LatticeConfig(rows=5, cols=5, backlash=0.02, z_coupling_gain=0.35)
        response = characterize_single_cell(config, (2, 2), alpha=-0.3, z=0.4)
        decay = response_decay_profile(response)
        self.assertEqual(decay.model, "log-linear shell max")
        self.assertGreater(decay.alpha_shells, 1)
        self.assertGreater(decay.z_shells, 1)
        self.assertGreaterEqual(decay.alpha_reach, 1)
        self.assertGreaterEqual(decay.z_reach, 1)
        self.assertTrue(np.isfinite(decay.alpha_ratio))
        self.assertTrue(np.isfinite(decay.z_ratio))
        self.assertGreaterEqual(decay.alpha_length, 0.0)
        self.assertGreaterEqual(decay.z_length, 0.0)

    def test_pairwise_interaction_graph_detects_nonadditive_pair(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.1, coupling_gain=0.5)
        graph = characterize_pairwise_interactions(
            config,
            (
                SourceCommand((1, 1), alpha=0.08),
                SourceCommand((1, 1), alpha=0.08),
            ),
        )
        self.assertEqual(graph.total_pair_count, 1)
        self.assertEqual(graph.evaluated_pair_count, 1)
        self.assertFalse(graph.truncated)
        self.assertEqual(graph.nonadditive_pair_count, 1)
        self.assertGreater(graph.max_alpha_error, 0.0)
        self.assertEqual(graph.interaction_hotspot_map.shape, (3, 3))
        self.assertGreater(graph.interaction_hotspot_map[1, 1], 0.0)
        self.assertEqual(graph.max_hotspot_error, graph.interaction_hotspot_map[1, 1])
        self.assertEqual(graph.interaction_degree_map.shape, (3, 3))
        self.assertEqual(graph.interaction_degree_map[1, 1], 1)
        self.assertEqual(graph.max_interaction_degree, 1)
        self.assertEqual(graph.interaction_density, 1.0)
        self.assertEqual(graph.interactions[0].manhattan_distance, 0)
        self.assertTrue(graph.interactions[0].nonadditive)
        np.testing.assert_allclose(graph.alpha_error_matrix, graph.alpha_error_matrix.T)

    def test_programmable_discontinuity_diagnostic_reports_locality_and_rank(self):
        config = LatticeConfig(rows=5, cols=5, backlash=0.02, z_coupling_gain=0.35)
        diagnostic = diagnose_programmable_discontinuity(
            config,
            (SourceCommand((2, 2), alpha=-0.3, z=0.4),),
        )
        self.assertEqual(diagnostic.active_operator_count, 1)
        self.assertGreaterEqual(diagnostic.alpha_locality_radius, 1)
        self.assertGreaterEqual(diagnostic.z_locality_radius, 1)
        self.assertGreater(diagnostic.reachable_alpha_cells, 1)
        self.assertGreater(diagnostic.reachable_height_cells, 1)
        self.assertGreater(diagnostic.alpha_rank, 0)
        self.assertGreater(diagnostic.height_rank, 0)
        self.assertIsNotNone(diagnostic.decay_profile)
        self.assertTrue(np.isfinite(diagnostic.alpha_decay_ratio))
        self.assertTrue(np.isfinite(diagnostic.z_decay_ratio))
        self.assertGreaterEqual(diagnostic.alpha_decay_length, 0.0)
        self.assertGreaterEqual(diagnostic.z_decay_length, 0.0)
        self.assertFalse(diagnostic.nonadditive)

    def test_programmable_discontinuity_diagnostic_reports_sequence_order(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        diagnostic = diagnose_programmable_discontinuity(
            config,
            (SourceCommand((1, 1), alpha=-0.3, z=0.2),),
            event_sequence=(
                local_actuation_event((1, 1), alpha=-0.3, z=0.2),
                lock_event((1, 1)),
                clear_actuation_event((1, 1)),
            ),
        )
        self.assertEqual(len(diagnostic.event_sequence), 3)
        self.assertIsNotNone(diagnostic.sequence_order)
        self.assertTrue(diagnostic.order_sensitive)
        self.assertGreater(diagnostic.noncommuting_adjacent_pairs, 0)
        self.assertGreater(diagnostic.max_order_error, 0.1)

    def test_programmable_discontinuity_diagnostic_counts_underactuated_regions(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.2, z_coupling_gain=0.0)
        diagnostic = diagnose_programmable_discontinuity(
            config,
            (SourceCommand((1, 1), alpha=-0.12),),
            include_z=False,
        )
        self.assertLess(diagnostic.reachable_alpha_cells, diagnostic.total_cells)
        self.assertGreater(diagnostic.alpha_underactuated_cells, 0)
        self.assertLess(diagnostic.reachable_height_cells, diagnostic.total_cells)
        self.assertGreater(diagnostic.height_underactuated_cells, 0)

    def test_programmable_discontinuity_diagnostic_detects_dead_zone_nonadditivity(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.1, coupling_gain=0.5)
        diagnostic = diagnose_programmable_discontinuity(
            config,
            (
                SourceCommand((1, 1), alpha=0.08),
                SourceCommand((1, 1), alpha=0.08),
            ),
        )
        self.assertTrue(diagnostic.nonadditive)
        self.assertGreater(diagnostic.alpha_superposition_error, 0.0)
        self.assertGreater(diagnostic.combined.alpha_reach, 1)
        self.assertIsNotNone(diagnostic.pairwise_interactions)
        self.assertEqual(diagnostic.nonadditive_pair_count, 1)
        self.assertGreater(diagnostic.max_pairwise_interaction_error, 0.0)
        self.assertEqual(diagnostic.pairwise_interaction_hotspot_map.shape, (3, 3))
        self.assertGreater(diagnostic.max_pairwise_hotspot_error, 0.0)
        self.assertEqual(diagnostic.pairwise_interaction_degree_map.shape, (3, 3))
        self.assertEqual(diagnostic.max_pairwise_interaction_degree, 1)
        self.assertEqual(diagnostic.pairwise_interaction_density, 1.0)

    def test_programmable_discontinuity_report_exports_research_artifact(self):
        config = LatticeConfig(
            rows=3,
            cols=3,
            backlash=0.1,
            coupling_gain=0.5,
            z_coupling_gain=0.25,
        )
        diagnostic = diagnose_programmable_discontinuity(
            config,
            (
                SourceCommand((1, 1), alpha=0.08, z=0.05),
                SourceCommand((1, 1), alpha=0.08, z=0.05),
            ),
            locked_cells=((0, 0),),
            event_sequence=(
                local_actuation_event((1, 1), alpha=0.08, z=0.05),
                lock_event((1, 1)),
                clear_actuation_event((1, 1)),
            ),
        )
        payload = programmable_discontinuity_report(
            diagnostic,
            config=config,
            include_response_matrix=False,
            include_fields=False,
        )

        self.assertEqual(payload["schema"], "rad-sim.programmable-discontinuity-report.v1")
        self.assertEqual(payload["config"]["pinHoleClearance"], config.pin_hole_clearance)
        self.assertEqual(payload["operators"]["activeOperatorCount"], 2)
        self.assertEqual(payload["operators"]["lockedCells"], [{"row": 0, "col": 0}])
        self.assertEqual(
            payload["paperSupportedAssumptions"][0]["formula"],
            "f(x)=max(0,x-b)+min(x+b,0)",
        )
        self.assertEqual(
            payload["simulatorDiagnostics"][1]["name"],
            "superposition residual",
        )
        laws = payload["operatorLawCandidates"]
        self.assertEqual(laws["schema"], "rad-sim.framework-law-candidates.v1")
        self.assertEqual(
            laws["method"],
            "thresholded diagnostic predicates over locality, reachability, composition, and event-order metrics",
        )
        law_by_id = {law["id"]: law for law in laws["laws"]}
        self.assertTrue(law_by_id["composition_nonadditivity"]["supportedByDiagnostic"])
        self.assertTrue(law_by_id["event_order_noncommutativity"]["supportedByDiagnostic"])
        self.assertEqual(
            law_by_id["composition_nonadditivity"]["evidence"]["maxSuperpositionError"],
            max(
                diagnostic.alpha_superposition_error,
                diagnostic.height_superposition_error,
            ),
        )
        formalization = payload["formalizationTargets"]
        self.assertEqual(formalization["schema"], "rad-sim.formalization-targets.v1")
        self.assertEqual(formalization["tooling"]["engine"], "Lean")
        self.assertIn("available", formalization["tooling"])
        target_by_id = {target["id"]: target for target in formalization["targets"]}
        self.assertEqual(
            target_by_id["dead_zone_zero_inside_backlash"]["source"],
            "paper-supported",
        )
        self.assertIn(
            "real max/min lemmas",
            target_by_id["dead_zone_zero_inside_backlash"]["dependencies"],
        )
        group_target = target_by_id["group_operator_support_decomposition"]
        self.assertEqual(group_target["status"], "lean-proved-discrete")
        self.assertIn(
            "LocalModeOperator.groupOperators_with_disjoint_lists_commute",
            group_target["evidence"]["leanTheorems"],
        )
        removal_group_target = target_by_id[
            "removed_cell_clears_group_supported_constraints"
        ]
        self.assertEqual(removal_group_target["status"], "lean-proved-discrete")
        self.assertIn(
            "CellGraph.removed_cell_constraints_clear_group_operator",
            removal_group_target["evidence"]["leanTheorems"],
        )
        clearance_target = target_by_id["vertical_clearance_gates_residual_contact"]
        self.assertEqual(clearance_target["status"], "lean-proved-discrete")
        self.assertIn(
            "verticalResidualStepNat_zero_inside",
            clearance_target["evidence"]["leanTheorems"],
        )
        self.assertIn(
            "contactPenaltyFromClearanceNat_nonnegative",
            clearance_target["evidence"]["leanTheorems"],
        )
        load_target = target_by_id["fixed_cell_load_work_proxy_zero"]
        self.assertEqual(load_target["status"], "lean-proved-discrete")
        self.assertIn(
            "loadWorkMagnitudeNat_zero_fixed",
            load_target["evidence"]["leanTheorems"],
        )
        signed_load_target = target_by_id["signed_vertical_load_work_residual_int"]
        self.assertEqual(signed_load_target["status"], "lean-proved-discrete")
        self.assertEqual(
            signed_load_target["evidence"]["cliModule"],
            "rad_sim.compare_vertical_load_measurements",
        )
        self.assertIn(
            "signedLoadWorkInt_zero_displacement",
            signed_load_target["evidence"]["leanTheorems"],
        )
        self.assertIn(
            "signedEnergyResidualTripleInt_zero_when_equal",
            signed_load_target["evidence"]["leanTheorems"],
        )
        scaled_target = target_by_id["integer_scaled_mechanics_scaffold"]
        self.assertEqual(scaled_target["status"], "lean-proved-discrete")
        self.assertIn(
            "scaledSpringEnergyNat_preserves_denominator",
            scaled_target["evidence"]["leanTheorems"],
        )
        self.assertIn(
            "scaledSignedEnergyResidualTripleInt_zero_when_equal",
            scaled_target["evidence"]["leanTheorems"],
        )
        self.assertIn("Rat/Real", scaled_target["evidence"]["nextFormalStep"])
        unit_scale_target = target_by_id["measurement_unit_scale_invariants"]
        self.assertEqual(unit_scale_target["status"], "lean-proved-discrete")
        self.assertEqual(
            unit_scale_target["evidence"]["pythonFunction"],
            "calibrate_paper_rad_config",
        )
        self.assertEqual(
            unit_scale_target["evidence"]["unitScaleSchema"],
            "rad-sim.physical-unit-scale-metadata.v1",
        )
        self.assertEqual(
            unit_scale_target["evidence"]["unitScaleFunction"],
            "physical_unit_scale_metadata",
        )
        self.assertEqual(
            unit_scale_target["evidence"]["hardwareProfileSchema"],
            HARDWARE_PROFILE_SCHEMA,
        )
        self.assertIn(
            "hardware_profile_from_json",
            unit_scale_target["evidence"]["hardwareProfileFunctions"],
        )
        self.assertIn(
            "measurementUnitScaleNat_preserves_denominator",
            unit_scale_target["evidence"]["leanTheorems"],
        )
        self.assertIn(
            "measurementUnitScaleResidualInt_zero_when_equal",
            unit_scale_target["evidence"]["leanTheorems"],
        )
        self.assertIn(
            "hardwareProfileCoverageMissing_zero_when_complete",
            unit_scale_target["evidence"]["leanTheorems"],
        )
        self.assertIn("physical-units", unit_scale_target["evidence"]["claimLimit"])
        calibration_estimate_target = target_by_id[
            "calibration_parameter_estimate_residual_bookkeeping"
        ]
        self.assertEqual(
            calibration_estimate_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            calibration_estimate_target["evidence"]["parameterEstimateSchema"],
            "rad-sim.calibration-parameter-estimates.v1",
        )
        self.assertIn(
            "calibrationFitResidualPassNat_zero",
            calibration_estimate_target["evidence"]["leanTheorems"],
        )
        calibration_profile_target = target_by_id[
            "calibration_model_profile_safe_update_bounds"
        ]
        self.assertEqual(
            calibration_profile_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            calibration_profile_target["evidence"]["profileSchema"],
            "rad-sim.calibration-model-profile.v1",
        )
        self.assertEqual(
            calibration_profile_target["evidence"]["applicationSchema"],
            "rad-sim.calibration-model-profile-application.v1",
        )
        self.assertIn(
            "calibrationModelProfileUpdateSafeNat_intro",
            calibration_profile_target["evidence"]["leanTheorems"],
        )
        calibration_profile_selection_target = target_by_id[
            "calibration_model_profile_selection_predicate"
        ]
        self.assertEqual(
            calibration_profile_selection_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            calibration_profile_selection_target["evidence"]["comparisonSchema"],
            "rad-sim.calibration-model-profile-residual-comparison.v1",
        )
        self.assertEqual(
            calibration_profile_selection_target["evidence"]["selectionSchema"],
            "rad-sim.calibration-model-profile-selection.v1",
        )
        self.assertIn(
            "calibrationModelProfileCandidateScoreNat_le_before",
            calibration_profile_selection_target["evidence"]["leanTheorems"],
        )
        calibration_profile_holdout_target = target_by_id[
            "calibration_model_profile_holdout_predicate"
        ]
        self.assertEqual(
            calibration_profile_holdout_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            calibration_profile_holdout_target["evidence"]["schema"],
            "rad-sim.calibration-model-profile-holdout-validation.v1",
        )
        self.assertEqual(
            calibration_profile_holdout_target["evidence"]["csvExportFunction"],
            "export_calibration_model_profile_holdout_validation_csv",
        )
        self.assertIn(
            "calibrationModelProfileHoldoutScoreNat_le_before",
            calibration_profile_holdout_target["evidence"]["leanTheorems"],
        )
        calibration_split_target = target_by_id[
            "calibration_train_holdout_split_metadata"
        ]
        self.assertEqual(calibration_split_target["status"], "lean-proved-discrete")
        self.assertEqual(
            calibration_split_target["evidence"]["schema"],
            "rad-sim.calibration-train-holdout-split.v1",
        )
        self.assertEqual(
            calibration_split_target["evidence"]["csvExportFunction"],
            "export_calibration_model_profile_holdout_validation_csv",
        )
        self.assertIn(
            "calibrationTrainHoldoutSplitReadyNat_zero_overlap",
            calibration_split_target["evidence"]["leanTheorems"],
        )
        calibration_provenance_target = target_by_id[
            "calibration_train_holdout_file_provenance"
        ]
        self.assertEqual(calibration_provenance_target["status"], "lean-proved-discrete")
        self.assertEqual(
            calibration_provenance_target["evidence"]["datasetProvenanceSchema"],
            "rad-sim.calibration-dataset-provenance.v1",
        )
        self.assertIn(
            "calibrationTrainHoldoutProvenanceReadyNat_profile_matches",
            calibration_provenance_target["evidence"]["leanTheorems"],
        )
        calibration_bench_target = target_by_id["calibration_bench_protocol_coverage"]
        self.assertEqual(calibration_bench_target["status"], "lean-proved-discrete")
        self.assertEqual(
            calibration_bench_target["evidence"]["schema"],
            "rad-sim.calibration-bench-notebook.v1",
        )
        self.assertEqual(
            calibration_bench_target["evidence"]["csvExportFunction"],
            "export_calibration_bench_notebook_csv",
        )
        self.assertIn(
            "calibrationBenchProtocolCoverageReadyNat_has_two_dataset_roles",
            calibration_bench_target["evidence"]["leanTheorems"],
        )
        calibration_packet_target = target_by_id["calibration_bench_packet_completeness"]
        self.assertEqual(calibration_packet_target["status"], "lean-proved-discrete")
        self.assertEqual(
            calibration_packet_target["evidence"]["schema"],
            "rad-sim.calibration-bench-packet.v1",
        )
        self.assertEqual(
            calibration_packet_target["evidence"]["writerFunction"],
            "write_calibration_bench_packet_artifacts",
        )
        self.assertIn(
            "calibrationBenchPacketCompleteNat_has_manifest",
            calibration_packet_target["evidence"]["leanTheorems"],
        )
        calibration_execution_target = target_by_id[
            "calibration_bench_executed_validation_gate"
        ]
        self.assertEqual(calibration_execution_target["status"], "lean-proved-discrete")
        self.assertEqual(
            calibration_execution_target["evidence"]["schema"],
            "rad-sim.calibration-bench-execution-validation.v1",
        )
        self.assertEqual(
            calibration_execution_target["evidence"]["writerFunction"],
            "write_calibration_bench_execution_validation_artifacts",
        )
        self.assertEqual(
            calibration_execution_target["evidence"]["browserFunction"],
            "calibrationBenchExecutionValidation",
        )
        self.assertIn(
            "calibrationBenchExecutedValidationReadyNat_zero_missing_evidence",
            calibration_execution_target["evidence"]["leanTheorems"],
        )
        energy_target = target_by_id["mechanics_energy_certificate_nonnegative_proxy"]
        self.assertEqual(energy_target["status"], "lean-proved-discrete")
        self.assertEqual(
            energy_target["evidence"]["schema"],
            "rad-sim.mechanics-energy-certificate.v1",
        )
        self.assertIn(
            "mechanicalStoredEnergyNat_nonnegative",
            energy_target["evidence"]["leanTheorems"],
        )
        measured_target = target_by_id[
            "measured_vertical_load_work_validation_zero_residual"
        ]
        self.assertEqual(measured_target["status"], "lean-proved-discrete")
        self.assertEqual(
            measured_target["evidence"]["schema"],
            "rad-sim.vertical-load-energy-validation.v1",
        )
        self.assertEqual(
            measured_target["evidence"]["protocolSchema"],
            "rad-sim.vertical-load-energy-experiment-protocol.v1",
        )
        self.assertEqual(
            measured_target["evidence"]["protocolFunction"],
            "vertical_load_energy_experiment_protocol",
        )
        self.assertEqual(
            measured_target["evidence"]["packetSchema"],
            "rad-sim.vertical-load-bench-packet.v1",
        )
        self.assertEqual(
            measured_target["evidence"]["packetFunction"],
            "vertical_load_bench_packet",
        )
        self.assertIn(
            "measuredWorkResidualNat_zero_when_equal",
            measured_target["evidence"]["leanTheorems"],
        )
        comparison_pass_target = target_by_id[
            "vertical_load_comparison_pass_predicate"
        ]
        self.assertEqual(comparison_pass_target["status"], "lean-proved-discrete")
        self.assertEqual(
            comparison_pass_target["evidence"]["cliModule"],
            "rad_sim.compare_vertical_load_measurements",
        )
        self.assertEqual(
            comparison_pass_target["evidence"]["reportSchema"],
            "rad-sim.vertical-load-energy-comparison-report.v1",
        )
        self.assertIn(
            "verticalLoadScenarioPassNat_zero_errors",
            comparison_pass_target["evidence"]["leanTheorems"],
        )
        self.assertIn(
            "verticalLoadScenarioPassNat_false_when_missing",
            comparison_pass_target["evidence"]["leanTheorems"],
        )
        physical_readiness_target = target_by_id["physical_validation_readiness_gate"]
        self.assertEqual(physical_readiness_target["status"], "lean-proved-discrete")
        self.assertEqual(
            physical_readiness_target["evidence"]["schema"],
            "rad-sim.physical-validation-readiness.v1",
        )
        self.assertEqual(
            physical_readiness_target["evidence"]["pythonFunction"],
            "physical_validation_readiness_report",
        )
        self.assertIn(
            "physicalValidationReadyNat_zero_missing_evidence",
            physical_readiness_target["evidence"]["leanTheorems"],
        )
        contact_state_target = target_by_id["contact_state_abstraction_gate"]
        self.assertEqual(contact_state_target["status"], "lean-proved-discrete")
        self.assertEqual(
            contact_state_target["evidence"]["schema"],
            "rad-sim.contact-state-abstraction.v1",
        )
        self.assertEqual(
            contact_state_target["evidence"]["pythonFunction"],
            "contact_state_abstraction_report",
        )
        self.assertIn(
            "contactStateAbstractionReadyNat_zero_missing_evidence",
            contact_state_target["evidence"]["leanTheorems"],
        )
        contact_graph_target = target_by_id["contact_graph_consistency_gate"]
        self.assertEqual(contact_graph_target["status"], "lean-proved-discrete")
        self.assertEqual(
            contact_graph_target["evidence"]["schema"],
            "rad-sim.contact-graph-consistency.v1",
        )
        self.assertEqual(
            contact_graph_target["evidence"]["pythonFunction"],
            "contact_graph_consistency_report",
        )
        self.assertIn(
            "contactGraphConsistentNat_zero_missing_evidence",
            contact_graph_target["evidence"]["leanTheorems"],
        )
        realization_target = target_by_id["physical_realization_map_gate"]
        self.assertEqual(realization_target["status"], "lean-proved-discrete")
        self.assertEqual(
            realization_target["evidence"]["schema"],
            "rad-sim.physical-realization-map.v1",
        )
        self.assertEqual(
            realization_target["evidence"]["pythonFunction"],
            "physical_realization_map_report",
        )
        self.assertIn(
            "physicalRealizationMapReadyNat_zero_missing_evidence",
            realization_target["evidence"]["leanTheorems"],
        )
        external_audit_target = target_by_id["external_physics_engine_audit_gate"]
        self.assertEqual(external_audit_target["status"], "lean-proved-discrete")
        self.assertEqual(
            external_audit_target["evidence"]["schema"],
            "rad-sim.external-physics-engine-audit.v1",
        )
        self.assertEqual(
            external_audit_target["evidence"]["pythonFunction"],
            "external_physics_engine_audit_report",
        )
        self.assertIn(
            "externalPhysicsEngineAuditReadyNat_zero_missing_evidence",
            external_audit_target["evidence"]["leanTheorems"],
        )
        mujoco_export_target = target_by_id["mujoco_model_export_gate"]
        self.assertEqual(mujoco_export_target["status"], "lean-proved-discrete")
        self.assertEqual(
            mujoco_export_target["evidence"]["schema"],
            "rad-sim.mujoco-model-export.v1",
        )
        self.assertEqual(
            mujoco_export_target["evidence"]["pythonFunction"],
            "mujoco_model_export_report",
        )
        self.assertIn(
            "externalPhysicsModelExportReadyNat_zero_missing_evidence",
            mujoco_export_target["evidence"]["leanTheorems"],
        )
        mujoco_contact_target = target_by_id["mujoco_pin_hole_contact_geometry_gate"]
        self.assertEqual(mujoco_contact_target["status"], "lean-proved-discrete")
        self.assertEqual(
            mujoco_contact_target["evidence"]["schema"],
            "rad-sim.mujoco-pin-hole-contact-geometry.v1",
        )
        self.assertEqual(
            mujoco_contact_target["evidence"]["pythonFunction"],
            "mujoco_pin_hole_contact_geometry_report",
        )
        self.assertIn(
            "externalContactGeometryReadyNat_zero_missing_evidence",
            mujoco_contact_target["evidence"]["leanTheorems"],
        )
        mujoco_parameter_target = target_by_id["mujoco_contact_parameter_profile_gate"]
        self.assertEqual(mujoco_parameter_target["status"], "lean-proved-discrete")
        self.assertEqual(
            mujoco_parameter_target["evidence"]["schema"],
            "rad-sim.mujoco-contact-parameter-profile.v1",
        )
        self.assertEqual(
            mujoco_parameter_target["evidence"]["pythonFunction"],
            "mujoco_contact_parameter_report",
        )
        self.assertIn(
            "externalContactParameterProfileReadyNat_zero_missing_evidence",
            mujoco_parameter_target["evidence"]["leanTheorems"],
        )
        contact_calibration_target = target_by_id[
            "contact_parameter_calibration_packet_completeness"
        ]
        self.assertEqual(contact_calibration_target["status"], "lean-proved-discrete")
        self.assertEqual(
            contact_calibration_target["evidence"]["schema"],
            "rad-sim.contact-parameter-calibration-packet.v1",
        )
        self.assertEqual(
            contact_calibration_target["evidence"]["pythonFunction"],
            "contact_parameter_calibration_packet",
        )
        self.assertIn(
            "contactParameterCalibrationPacketCompleteNat_zero_missing_evidence",
            contact_calibration_target["evidence"]["leanTheorems"],
        )
        contact_bench_target = target_by_id["contact_parameter_bench_validation_gate"]
        self.assertEqual(contact_bench_target["status"], "lean-proved-discrete")
        self.assertEqual(
            contact_bench_target["evidence"]["schema"],
            "rad-sim.contact-parameter-bench-validation.v1",
        )
        self.assertEqual(
            contact_bench_target["evidence"]["pythonFunction"],
            "compare_contact_parameter_calibration_results",
        )
        self.assertIn(
            "contactParameterBenchValidationReadyNat_zero_missing_evidence",
            contact_bench_target["evidence"]["leanTheorems"],
        )
        contact_interval_target = target_by_id[
            "contact_parameter_interval_calibration_gate"
        ]
        self.assertEqual(contact_interval_target["status"], "lean-proved-discrete")
        self.assertEqual(
            contact_interval_target["evidence"]["schema"],
            "rad-sim.contact-parameter-interval-calibration.v1",
        )
        self.assertEqual(
            contact_interval_target["evidence"]["pythonFunction"],
            "contact_parameter_interval_calibration_report",
        )
        self.assertIn(
            "contactParameterIntervalCalibrationReadyNat_zero_missing_evidence",
            contact_interval_target["evidence"]["leanTheorems"],
        )
        mujoco_run_target = target_by_id["mujoco_external_run_gate"]
        self.assertEqual(mujoco_run_target["status"], "lean-proved-discrete")
        self.assertEqual(
            mujoco_run_target["evidence"]["schema"],
            "rad-sim.mujoco-external-run.v1",
        )
        self.assertEqual(
            mujoco_run_target["evidence"]["pythonFunction"],
            "mujoco_external_run_report",
        )
        self.assertIn(
            "externalPhysicsRunReadyNat_zero_missing_evidence",
            mujoco_run_target["evidence"]["leanTheorems"],
        )
        mujoco_comparison_target = target_by_id["mujoco_external_comparison_gate"]
        self.assertEqual(mujoco_comparison_target["status"], "lean-proved-discrete")
        self.assertEqual(
            mujoco_comparison_target["evidence"]["schema"],
            "rad-sim.mujoco-external-comparison.v1",
        )
        self.assertEqual(
            mujoco_comparison_target["evidence"]["pythonFunction"],
            "mujoco_external_comparison_report",
        )
        self.assertIn(
            "externalPhysicsComparisonReadyNat_zero_missing_evidence",
            mujoco_comparison_target["evidence"]["leanTheorems"],
        )
        equilibrium_target = target_by_id["equilibrium_relation_gate"]
        self.assertEqual(equilibrium_target["status"], "lean-proved-discrete")
        self.assertEqual(
            equilibrium_target["evidence"]["schema"],
            "rad-sim.equilibrium-relation.v1",
        )
        self.assertEqual(
            equilibrium_target["evidence"]["pythonFunction"],
            "equilibrium_relation_report",
        )
        self.assertIn(
            "equilibriumRelationReadyNat_zero_missing_evidence",
            equilibrium_target["evidence"]["leanTheorems"],
        )
        reachable_equilibrium_target = target_by_id[
            "reachable_equilibrium_controllability_gate"
        ]
        self.assertEqual(
            reachable_equilibrium_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            reachable_equilibrium_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-controllability.v1",
        )
        self.assertEqual(
            reachable_equilibrium_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_controllability_report",
        )
        self.assertIn(
            "reachableEquilibriumControllabilityReadyNat_zero_missing_evidence",
            reachable_equilibrium_target["evidence"]["leanTheorems"],
        )
        reachable_bench_target = target_by_id[
            "reachable_equilibrium_bench_protocol_gate"
        ]
        self.assertEqual(reachable_bench_target["status"], "lean-proved-discrete")
        self.assertEqual(
            reachable_bench_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-bench-protocol.v1",
        )
        self.assertEqual(
            reachable_bench_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_bench_protocol",
        )
        self.assertIn(
            "reachableEquilibriumBenchProtocolReadyNat_zero_missing_evidence",
            reachable_bench_target["evidence"]["leanTheorems"],
        )
        reachable_bench_validation_target = target_by_id[
            "reachable_equilibrium_bench_validation_gate"
        ]
        self.assertEqual(
            reachable_bench_validation_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            reachable_bench_validation_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-bench-comparison.v1",
        )
        self.assertEqual(
            reachable_bench_validation_target["evidence"]["pythonFunction"],
            "compare_reachable_equilibrium_bench_results",
        )
        self.assertIn(
            "reachableEquilibriumBenchValidationReadyNat_zero_missing_evidence",
            reachable_bench_validation_target["evidence"]["leanTheorems"],
        )
        amplitude_target = target_by_id[
            "reachable_equilibrium_amplitude_calibration_gate"
        ]
        self.assertEqual(amplitude_target["status"], "lean-proved-discrete")
        self.assertEqual(
            amplitude_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
        )
        self.assertEqual(
            amplitude_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_amplitude_calibration_report",
        )
        self.assertIn(
            "reachableEquilibriumAmplitudeCalibrationReadyNat_zero_missing_evidence",
            amplitude_target["evidence"]["leanTheorems"],
        )
        empirical_profile_target = target_by_id[
            "reachable_equilibrium_empirical_profile_gate"
        ]
        self.assertEqual(
            empirical_profile_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            empirical_profile_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-empirical-profile.v1",
        )
        self.assertEqual(
            empirical_profile_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_empirical_profile_from_amplitude",
        )
        self.assertIn(
            "reachableEquilibriumEmpiricalProfileReadyNat_zero_missing_evidence",
            empirical_profile_target["evidence"]["leanTheorems"],
        )
        profile_inverse_target = target_by_id[
            "reachable_equilibrium_profile_inverse_gate"
        ]
        self.assertEqual(
            profile_inverse_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            profile_inverse_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse.v1",
        )
        self.assertEqual(
            profile_inverse_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_profile_inverse_report",
        )
        self.assertIn(
            "reachableEquilibriumProfileInverseReadyNat_zero_missing_evidence",
            profile_inverse_target["evidence"]["leanTheorems"],
        )
        profile_inverse_acceptance_target = target_by_id[
            "reachable_equilibrium_profile_inverse_acceptance_gate"
        ]
        self.assertEqual(
            profile_inverse_acceptance_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            profile_inverse_acceptance_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
        )
        self.assertEqual(
            profile_inverse_acceptance_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_profile_inverse_acceptance_report",
        )
        self.assertIn(
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_zero_missing_evidence",
            profile_inverse_acceptance_target["evidence"]["leanTheorems"],
        )
        profile_inverse_packet_target = target_by_id[
            "reachable_equilibrium_profile_inverse_preview_packet_gate"
        ]
        self.assertEqual(
            profile_inverse_packet_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            profile_inverse_packet_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
        )
        self.assertEqual(
            profile_inverse_packet_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_profile_inverse_preview_packet",
        )
        self.assertIn(
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_zero_missing_evidence",
            profile_inverse_packet_target["evidence"]["leanTheorems"],
        )
        profile_inverse_replay_target = target_by_id[
            "reachable_equilibrium_profile_inverse_preview_replay_gate"
        ]
        self.assertEqual(
            profile_inverse_replay_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            profile_inverse_replay_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
        )
        self.assertEqual(
            profile_inverse_replay_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_profile_inverse_preview_replay_report",
        )
        self.assertIn(
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_zero_missing_evidence",
            profile_inverse_replay_target["evidence"]["leanTheorems"],
        )
        profile_inverse_physical_target = target_by_id[
            "reachable_equilibrium_profile_inverse_preview_physical_gate"
        ]
        self.assertEqual(
            profile_inverse_physical_target["status"],
            "lean-proved-discrete",
        )
        self.assertEqual(
            profile_inverse_physical_target["evidence"]["schema"],
            "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
        )
        self.assertEqual(
            profile_inverse_physical_target["evidence"]["pythonFunction"],
            "reachable_equilibrium_profile_inverse_preview_physical_report",
        )
        self.assertIn(
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_zero_missing_evidence",
            profile_inverse_physical_target["evidence"]["leanTheorems"],
        )
        spring_target = target_by_id["spring_hinge_removed_topology_load_comparison"]
        self.assertEqual(spring_target["status"], "simulator-diagnostic")
        self.assertEqual(
            spring_target["evidence"]["pythonFunction"],
            "compare_vertical_residual_spring_hinge_3d",
        )
        self.assertEqual(
            spring_target["evidence"]["reportSchema"],
            "rad-sim.vertical-load-physical-preview-report.v1",
        )
        self.assertEqual(
            spring_target["evidence"]["csvExportFunction"],
            "export_vertical_load_physical_preview_report_csv",
        )
        browser_preview_target = target_by_id[
            "browser_spring_preview_removed_edge_deletion"
        ]
        self.assertEqual(browser_preview_target["status"], "simulator-diagnostic")
        self.assertEqual(
            browser_preview_target["evidence"]["browserFunction"],
            "RAD.simulatePhysicalRelaxation",
        )
        self.assertIn(
            "physicalSkippedSpringEdges",
            browser_preview_target["evidence"]["browserMetrics"],
        )
        self.assertIn(
            "CellGraph.removal_deletes_one_step_to_removed",
            browser_preview_target["evidence"]["leanTheorems"],
        )
        self.assertTrue(
            target_by_id["noncommutativity_witness_from_order_error"]["evidence"][
                "orderSensitive"
            ]
        )
        standalone_formalization = formalization_target_manifest(diagnostic)
        self.assertEqual(standalone_formalization, formalization)
        self.assertTrue(payload["composition"]["nonadditive"])
        self.assertTrue(payload["composition"]["orderSensitive"])
        self.assertIsNotNone(payload["sequenceOrder"])
        self.assertEqual(
            payload["combinedResponse"]["positiveZReachCells"],
            diagnostic.combined.positive_z_reach,
        )
        self.assertEqual(
            payload["combinedResponse"]["negativeZReachCells"],
            diagnostic.combined.negative_z_reach,
        )
        self.assertEqual(
            payload["combinedResponse"]["maxPositiveHeightDelta"],
            diagnostic.combined.max_positive_height_delta,
        )
        self.assertEqual(
            payload["combinedResponse"]["maxNegativeHeightDelta"],
            diagnostic.combined.max_negative_height_delta,
        )
        self.assertNotIn("fields", payload["combinedResponse"])
        self.assertNotIn("alpha", payload["responseMatrix"])
        self.assertEqual(
            payload["responseMatrix"]["diagnostics"]["reachableHeightCells"],
            diagnostic.reachable_height_cells,
        )

        exported = json.loads(
            export_programmable_discontinuity_report_json(
                diagnostic,
                config=config,
                include_response_matrix=False,
                include_fields=False,
            )
        )
        self.assertEqual(exported["schema"], payload["schema"])
        self.assertEqual(exported["composition"], payload["composition"])
        self.assertEqual(exported["operatorLawCandidates"], payload["operatorLawCandidates"])
        self.assertEqual(exported["formalizationTargets"], payload["formalizationTargets"])
        exported_formalization = json.loads(
            export_formalization_target_manifest_json(diagnostic)
        )
        self.assertEqual(exported_formalization, payload["formalizationTargets"])

    def test_programmable_discontinuity_report_can_include_physical_validation(self):
        config = LatticeConfig(rows=2, cols=2, z_coupling_gain=0.0)
        diagnostic = diagnose_programmable_discontinuity(
            config,
            (
                SourceCommand((0, 1), z=0.18),
                SourceCommand((1, 0), alpha=-0.12),
            ),
            include_physical=True,
            load_case=LoadCase(lock_stiffness=350.0, maxiter=250),
        )
        self.assertIsNotNone(diagnostic.physical_validation)
        self.assertTrue(diagnostic.physical_validation.physical_success)

        payload = programmable_discontinuity_report(
            diagnostic,
            config=config,
            include_response_matrix=False,
            include_fields=False,
        )
        physical = payload["physicalValidation"]
        self.assertEqual(physical["schema"], "rad-sim.physical-response-comparison.v1")
        self.assertEqual(physical["model"], "spring_hinge_3d")
        self.assertEqual(len(physical["commands"]), 2)
        self.assertTrue(physical["physicalSuccess"])
        self.assertGreaterEqual(physical["heightRmsModelError"], 0.0)
        self.assertGreaterEqual(physical["centerRmsModelError"], 0.0)
        self.assertGreater(physical["iterations"], 0)
        self.assertNotIn("fields", physical)

        exported = json.loads(
            export_programmable_discontinuity_report_json(
                diagnostic,
                config=config,
                include_response_matrix=False,
                include_fields=False,
            )
        )
        self.assertEqual(exported["physicalValidation"]["schema"], physical["schema"])

    def test_physical_response_comparison_reports_spring_hinge_deviation(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.0)
        comparison = compare_physical_response(
            config,
            (SourceCommand((1, 1), z=0.35),),
            load_case=LoadCase(lock_stiffness=500.0, maxiter=400),
        )
        self.assertTrue(comparison.physical_success)
        self.assertEqual(comparison.physical_height_delta.shape, (3, 3))
        self.assertEqual(comparison.physical_center_delta.shape, (3, 3, 3))
        self.assertGreater(comparison.physical_height_delta[1, 1], 0.1)
        self.assertGreaterEqual(comparison.height_rms_error, 0.0)
        self.assertGreaterEqual(comparison.center_rms_error, 0.0)
        self.assertGreaterEqual(comparison.physical_energy, 0.0)

    def test_physical_pair_comparison_reports_pairwise_nonlocality(self):
        config = LatticeConfig(rows=2, cols=2, z_coupling_gain=0.0)
        pair = compare_physical_pair(
            config,
            SourceCommand((0, 1), z=0.2),
            SourceCommand((1, 0), z=0.2),
            load_case=LoadCase(lock_stiffness=400.0, maxiter=300),
        )
        self.assertTrue(pair.combined.physical_success)
        self.assertEqual(len(pair.combined.commands), 2)
        self.assertGreaterEqual(pair.physical_height_superposition_error, 0.0)
        self.assertGreaterEqual(pair.physical_center_superposition_error, 0.0)
        self.assertTrue(np.isfinite(pair.combined.center_rms_error))

    def test_physical_cluster_comparison_accepts_multiple_commands(self):
        config = LatticeConfig(rows=2, cols=2, z_coupling_gain=0.0)
        cluster = compare_physical_cluster(
            config,
            (
                SourceCommand((0, 1), z=0.18),
                SourceCommand((1, 0), alpha=-0.12),
                SourceCommand((1, 1), z=-0.08),
            ),
            load_case=LoadCase(lock_stiffness=350.0, maxiter=250),
        )
        self.assertTrue(cluster.physical_success)
        self.assertEqual(len(cluster.commands), 3)
        self.assertEqual(cluster.physical_center_delta.shape, (2, 2, 3))
        self.assertTrue(np.isfinite(cluster.max_abs_center_error))

    def test_inverse_design_reduces_height_target_error(self):
        config = LatticeConfig(rows=5, cols=5, z_coupling_gain=0.25)
        target = np.zeros((5, 5), dtype=float)
        target[2, 2] = 0.36
        target[1:4, 2] = [0.1, 0.36, 0.1]
        target[2, 1:4] = [0.1, 0.36, 0.1]
        solution = solve_inverse_design(
            config,
            target_height=target,
            include_alpha=False,
            include_z=True,
        )
        self.assertTrue(solution.success)
        self.assertGreater(solution.active_actuator_count, 0)
        self.assertLess(solution.rms_error_after, solution.rms_error_before)
        self.assertGreater(solution.state.z_actuator_grid[2, 2], 0.0)
        self.assertEqual(solution.target_height.shape, target.shape)
        self.assertEqual(solution.height_residual.shape, target.shape)
        self.assertEqual(solution.reachable_height_mask.shape, target.shape)
        self.assertEqual(solution.underactuated_height_mask.shape, target.shape)
        self.assertGreaterEqual(solution.height_rms_residual, 0.0)
        self.assertGreaterEqual(solution.max_abs_height_residual, 0.0)
        self.assertGreaterEqual(solution.saturated_column_fraction, 0.0)
        self.assertLessEqual(solution.saturated_column_fraction, 1.0)
        payload = solution.to_dict(include_response_matrix=False)
        self.assertEqual(payload["schema"], "rad-sim.inverse-design-result.v1")
        self.assertEqual(payload["metrics"]["activeActuatorCount"], solution.active_actuator_count)
        self.assertEqual(payload["metrics"]["heightUnderactuatedCells"], solution.height_underactuated_cells)
        self.assertEqual(payload["target"]["height"][2][2], target[2, 2])
        report = inverse_design_report(solution, include_response_matrix=False)
        self.assertEqual(report["schema"], "rad-sim.inverse-design-report.v1")
        self.assertEqual(report["inverse"]["responseMatrix"]["schema"], "rad-sim.response-matrix.v1")
        self.assertEqual(json.loads(export_inverse_design_report_json(solution, include_response_matrix=False))["schema"], report["schema"])

    def test_inverse_design_fits_alpha_contraction(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        target_alpha = np.full((3, 3), config.initial_alpha)
        target_alpha[1, 1] = 0.75
        solution = solve_inverse_design(
            config,
            target_alpha=target_alpha,
            actuator_cells=[(1, 1)],
            include_alpha=True,
            include_z=False,
            alpha_weight=1.0,
        )
        self.assertTrue(solution.success)
        self.assertLess(solution.rms_error_after, solution.rms_error_before)
        self.assertLess(solution.state.actuator_grid[1, 1], 0.0)
        self.assertLess(solution.result.alpha[1, 1], config.initial_alpha)

    def test_inverse_design_excludes_locked_actuators(self):
        config = LatticeConfig(rows=3, cols=3)
        target = np.zeros((3, 3), dtype=float)
        target[1, 1] = 0.3
        solution = solve_inverse_design(
            config,
            target_height=target,
            actuator_cells=[(1, 1), (0, 1)],
            locked_cells=[(1, 1)],
            include_alpha=False,
            include_z=True,
        )
        self.assertTrue(solution.state.locked_mask[1, 1])
        self.assertAlmostEqual(solution.state.z_actuator_grid[1, 1], 0.0)
        self.assertTrue(all(command.cell != (1, 1) for command in solution.commands))

    def test_inverse_design_marks_underactuated_height_targets(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.0)
        target = np.zeros((3, 3), dtype=float)
        target[2, 2] = 0.25
        solution = solve_inverse_design(
            config,
            target_height=target,
            actuator_cells=[(0, 0)],
            include_alpha=False,
            include_z=True,
        )
        self.assertTrue(solution.reachable_height_mask[0, 0])
        self.assertFalse(solution.reachable_height_mask[2, 2])
        self.assertTrue(solution.underactuated_height_mask[2, 2])
        self.assertEqual(solution.height_underactuated_cells, 1)
        self.assertEqual(solution.positive_height_underactuated_cells, 1)
        self.assertEqual(solution.negative_height_underactuated_cells, 0)
        self.assertEqual(solution.alpha_underactuated_cells, 0)
        self.assertGreater(solution.height_rms_residual, 0.0)
        payload = solution.to_dict(include_response_matrix=False)
        self.assertEqual(payload["metrics"]["positiveHeightUnderactuatedCells"], 1)
        self.assertEqual(payload["metrics"]["negativeHeightUnderactuatedCells"], 0)
        self.assertFalse(payload["reachability"]["positiveHeight"][2][2])

    def test_inverse_design_marks_topology_blocked_removed_component_targets(self):
        config = LatticeConfig(rows=1, cols=3, backlash=0.0, coupling_gain=1.0)
        target = np.zeros((1, 3), dtype=float)
        target[0, 2] = 0.3
        solution = solve_inverse_design(
            config,
            target_height=target,
            actuator_cells=[(0, 0)],
            removed_cells=[(0, 1)],
            alpha_weight=0.0,
            height_weight=1.0,
        )
        self.assertTrue(solution.state.removed_mask[0, 1])
        self.assertFalse(solution.topology_reachable_mask[0, 2])
        self.assertTrue(solution.topology_blocked_height_mask[0, 2])
        self.assertEqual(solution.topology_blocked_height_cells, 1)
        payload = solution.to_dict(include_response_matrix=False)
        self.assertTrue(payload["reachability"]["topologyBlockedHeight"][0][2])
        self.assertEqual(payload["metrics"]["topologyBlockedHeightCells"], 1)

    def test_inverse_design_tracks_downward_underactuated_height_targets(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.0)
        target = np.zeros((3, 3), dtype=float)
        target[2, 2] = -0.25
        solution = solve_inverse_design(
            config,
            target_height=target,
            actuator_cells=[(0, 0)],
            include_alpha=False,
            include_z=True,
        )
        self.assertTrue(solution.underactuated_height_mask[2, 2])
        self.assertEqual(solution.positive_height_underactuated_cells, 0)
        self.assertEqual(solution.negative_height_underactuated_cells, 1)
        payload = solution.to_dict(include_response_matrix=False)
        self.assertEqual(payload["metrics"]["positiveHeightUnderactuatedCells"], 0)
        self.assertEqual(payload["metrics"]["negativeHeightUnderactuatedCells"], 1)

    def test_inverse_design_signed_height_reachability_uses_directional_probes(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.3)
        target = np.zeros((3, 3), dtype=float)
        target[1, 1] = 0.2
        solution = solve_inverse_design(
            config,
            target_height=target,
            actuator_cells=[(1, 1)],
            include_alpha=False,
            include_z=True,
            z_step=0.12,
            tolerance=1e-8,
        )
        upward = characterize_response(
            config,
            (SourceCommand((1, 1), z=0.12),),
            tolerance=1e-8,
        )
        downward = characterize_response(
            config,
            (SourceCommand((1, 1), z=-0.12),),
            tolerance=1e-8,
        )
        np.testing.assert_array_equal(
            solution.positive_height_reachable_mask,
            upward.height_delta > 1e-8,
        )
        np.testing.assert_array_equal(
            solution.negative_height_reachable_mask,
            downward.height_delta < -1e-8,
        )

    def test_inverse_design_physical_validation_reports_spring_hinge_fit(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.0)
        target = np.zeros((3, 3), dtype=float)
        target[1, 1] = 0.24
        solution = solve_inverse_design(
            config,
            target_height=target,
            actuator_cells=[(1, 1)],
            include_alpha=False,
            include_z=True,
        )
        validation = validate_inverse_design_physical(
            solution,
            LoadCase(lock_stiffness=600.0, maxiter=350),
        )
        self.assertTrue(validation.physical_success)
        self.assertEqual(validation.height_residual_after.shape, target.shape)
        self.assertEqual(validation.height_model_error.shape, target.shape)
        self.assertEqual(validation.center_model_error.shape, (3, 3, 3))
        self.assertGreater(validation.physical_rms_height_error_before, 0.0)
        self.assertGreaterEqual(validation.physical_rms_height_error_after, 0.0)
        self.assertTrue(np.isfinite(validation.physical_height_error_improvement))
        self.assertGreaterEqual(validation.model_agreement_score, 0.0)
        self.assertLessEqual(validation.model_agreement_score, 1.0)
        payload = validation.to_dict()
        self.assertEqual(payload["schema"], "rad-sim.inverse-physical-validation.v1")
        self.assertAlmostEqual(
            payload["fields"]["centerModelErrorNorm"][1][1],
            np.linalg.norm(validation.center_model_error[1, 1]),
        )
        report = inverse_design_report(solution, validation, include_response_matrix=False)
        self.assertEqual(report["physicalValidation"]["schema"], payload["schema"])
        self.assertTrue(np.isfinite(validation.center_rms_model_error))
        self.assertTrue(np.isfinite(validation.physical_energy))

    def test_inverse_design_requires_target(self):
        with self.assertRaises(ValueError):
            solve_inverse_design(LatticeConfig(rows=2, cols=2))

    def test_zero_vertical_coupling_recovers_local_z_motion(self):
        config = LatticeConfig(rows=5, cols=5, backlash=0.02, z_coupling_gain=0.0)
        state = LatticeState.uniform(config)
        state.z_actuator_grid[2, 2] = 0.4
        result = simulate_kinematic(config, state)
        height = result.metadata["height"]
        self.assertAlmostEqual(height[2, 2], 0.4)
        self.assertAlmostEqual(height[2, 3], 0.0)

    def test_spring_energy_formula(self):
        points = np.array([[0.0, 0.0], [1.2, 0.0]])
        energy = spring_energy(points, [(0, 1)], np.array([1.0]), 10.0)
        self.assertAlmostEqual(energy, 0.5 * 10.0 * 0.2**2)

    def test_hinge_energy_formula(self):
        angle = math.pi / 2
        points = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 1.0]])
        energy = hinge_energy(points, [(0, 1, 2)], np.array([math.pi]), 2.0)
        self.assertAlmostEqual(energy, 0.5 * 2.0 * (angle - math.pi) ** 2)

    def test_spring_hinge_3d_relaxes_vertical_actuation(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.0)
        state = LatticeState.uniform(config)
        state.z_actuator_grid[1, 1] = 0.35
        result = solve_spring_hinge_3d(
            config,
            state,
            LoadCase(lock_stiffness=500.0, maxiter=400),
        )
        self.assertEqual(result.metadata["model"], "spring_hinge_3d")
        self.assertTrue(result.metadata["success"])
        self.assertGreater(result.deformed_centers_3d[1, 1, 2], 0.1)
        self.assertAlmostEqual(result.deformed_centers_3d[0, 0, 2], 0.0)
        self.assertGreaterEqual(result.metadata["stored_energy"], 0.0)
        self.assertAlmostEqual(
            result.metadata["objective_energy"],
            result.metadata["stored_energy"]
            + result.metadata["external_potential_energy"],
            places=7,
        )
        np.testing.assert_allclose(
            result.metadata["height"], result.deformed_centers_3d[..., 2]
        )

    def test_spring_hinge_3d_respects_locked_height_target(self):
        config = LatticeConfig(rows=3, cols=3, z_coupling_gain=0.0)
        state = LatticeState.uniform(config)
        state.locked_mask[1, 1] = True
        state.z_actuator_grid[1, 1] = 0.4
        result = solve_spring_hinge_3d(
            config,
            state,
            LoadCase(lock_stiffness=800.0, maxiter=400),
        )
        self.assertTrue(result.metadata["success"])
        self.assertLess(abs(result.deformed_centers_3d[1, 1, 2]), 0.02)

    def test_position_locked_cell_keeps_physical_coordinates(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.4)
        state = LatticeState.uniform(config)
        state.position_locked_mask[1, 1] = True
        state.actuator_grid[1, 1] = -0.4
        state.z_actuator_grid[1, 1] = 0.5

        kinematic = simulate_kinematic(config, state)
        np.testing.assert_allclose(
            kinematic.deformed_centers_3d[1, 1],
            kinematic.original_centers_3d[1, 1],
        )

        physical_2d = solve_spring_hinge(
            config,
            state,
            LoadCase(fixed_cells=(), lock_stiffness=800.0, maxiter=400),
        )
        np.testing.assert_allclose(
            physical_2d.deformed_centers[1, 1],
            physical_2d.original_centers[1, 1],
            atol=1e-9,
        )
        self.assertEqual(physical_2d.metadata["position_locked_cells"], 1)

        physical_3d = solve_spring_hinge_3d(
            config,
            state,
            LoadCase(fixed_cells=(), lock_stiffness=800.0, maxiter=400),
        )
        np.testing.assert_allclose(
            physical_3d.deformed_centers_3d[1, 1],
            physical_3d.original_centers_3d[1, 1],
            atol=1e-9,
        )
        self.assertEqual(physical_3d.metadata["position_locked_cells"], 1)

    def test_auxetic_poisson_trend(self):
        c1 = LatticeConfig(rows=4, cols=4)
        expanded = simulate_kinematic(c1, LatticeState.uniform(c1, alpha=1.2))
        contracted = simulate_kinematic(c1, LatticeState.uniform(c1, alpha=0.8))
        width_e = np.ptp(expanded.deformed_centers[..., 0])
        height_e = np.ptp(expanded.deformed_centers[..., 1])
        width_c = np.ptp(contracted.deformed_centers[..., 0])
        height_c = np.ptp(contracted.deformed_centers[..., 1])
        axial_strain = (width_e - width_c) / width_c
        transverse_strain = (height_e - height_c) / height_c
        poisson = -transverse_strain / axial_strain
        self.assertLess(poisson, 0.0)


if __name__ == "__main__":
    unittest.main()
