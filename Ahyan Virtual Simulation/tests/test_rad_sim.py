import math
import unittest

import numpy as np

from rad_sim import LatticeConfig, LatticeState, backlash_activation, simulate_kinematic
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
