import json
import math
import unittest

import numpy as np

from rad_sim import (
    DeadZonePropagationOperator,
    HARDWARE_PROFILE_DIMENSIONS,
    LatticeConfig,
    LatticeState,
    LoadCase,
    CALIBRATION_SOLVER_GAPS,
    CALIBRATION_SOLVER_TASKS,
    PAPER_RAD_REFERENCE,
    RADHardwareProfile,
    SourceCommand,
    apply_event_sequence,
    backlash_activation,
    build_calibration_experiment_protocol,
    build_response_matrix,
    build_paper_rad_cell_geometry,
    build_paper_rad_lattice_geometry,
    build_paper_rad_lattice_mesh,
    calibrate_paper_rad_config,
    calibration_experiment_measurements_from_json,
    calibration_experiment_results_template,
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
    compare_calibration_experiment_measurements,
    compare_event_order,
    compare_sequence_order,
    config_with_hardware_profile,
    clear_actuation_event,
    diagnose_programmable_discontinuity,
    evaluate_programmable_operators,
    finite_die_off_radius,
    export_calibration_experiment_protocol_json,
    export_calibration_experiment_comparison_json,
    export_calibration_experiment_results_template_json,
    hardware_profile_from_config,
    local_actuation_event,
    lock_event,
    lock_projection,
    model_provenance,
    provenance_by_status,
    provenance_summary,
    export_paper_rad_mesh_obj,
    iter_obj_vertices,
    release_event,
    run_calibration_experiment_protocol,
    response_decay_profile,
    simulate_kinematic,
    solve_inverse_design,
    solve_spring_hinge_3d,
    validate_inverse_design_physical,
    vertical_clearance_operator,
)
from rad_sim.coupling import alpha_to_theta, theta_to_alpha
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

    def test_lock_event_commits_current_alpha_state(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.0)
        state = LatticeState.uniform(config)
        final = apply_event_sequence(
            config,
            state,
            (
                local_actuation_event((1, 1), alpha=-0.3),
                lock_event((1, 1)),
                clear_actuation_event(),
            ),
        )
        self.assertTrue(final.locked_mask[1, 1])
        self.assertAlmostEqual(final.alpha_grid[1, 1], 0.7)
        fields = evaluate_programmable_operators(config, final)
        self.assertAlmostEqual(fields["alpha"][1, 1], 0.7)
        self.assertAlmostEqual(fields["height"][1, 1], 0.0)

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

    def test_calibration_results_template_roundtrips_and_compares_to_simulation(self):
        config = LatticeConfig(rows=3, cols=3, backlash=0.02, z_coupling_gain=0.35)
        protocol = build_calibration_experiment_protocol(
            config,
            center_cell=(1, 1),
            repeat_count=2,
        )
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
        self.assertAlmostEqual(comparison.mean_pin_hole_slip_mm or 0.0, 0.12)
        self.assertAlmostEqual(comparison.mean_actuator_force_n or 0.0, 3.4)
        exported = json.loads(export_calibration_experiment_comparison_json(comparisons))
        self.assertEqual(exported["schema"], "rad-sim.calibration-experiment-comparison.v1")

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
        self.assertEqual(solution.alpha_underactuated_cells, 0)
        self.assertGreater(solution.height_rms_residual, 0.0)

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
