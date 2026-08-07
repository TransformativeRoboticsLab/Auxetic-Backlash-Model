import math
import unittest

import numpy as np

from rad_sim import (
    DeadZonePropagationOperator,
    LatticeConfig,
    LatticeState,
    LoadCase,
    PAPER_RAD_REFERENCE,
    SourceCommand,
    backlash_activation,
    build_response_matrix,
    build_paper_rad_cell_geometry,
    build_paper_rad_lattice_geometry,
    characterize_cluster,
    characterize_pair,
    characterize_single_cell,
    evaluate_programmable_operators,
    lock_projection,
    simulate_kinematic,
    solve_inverse_design,
    solve_spring_hinge_3d,
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
