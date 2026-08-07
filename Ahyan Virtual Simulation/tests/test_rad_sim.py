import math
import unittest

import numpy as np

from rad_sim import (
    DeadZonePropagationOperator,
    LatticeConfig,
    LatticeState,
    SourceCommand,
    backlash_activation,
    build_response_matrix,
    characterize_cluster,
    characterize_pair,
    characterize_single_cell,
    evaluate_programmable_operators,
    lock_projection,
    simulate_kinematic,
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
