import unittest

import numpy as np

from simulator import MonteCarloSIRSimulator, infection_probability, periodic_distance


class InfectionProbabilityTests(unittest.TestCase):
    def test_normal_formula_and_cutoff(self) -> None:
        distances = np.array([0.0, 0.5, 1.0, 1.01])
        probabilities = infection_probability(
            distances,
            source_is_superspreader=False,
            model="normal",
            r0=1.0,
            w0=1.0,
        )
        np.testing.assert_allclose(probabilities, [1.0, 0.25, 0.0, 0.0])

    def test_strong_superspreader_formula(self) -> None:
        distances = np.array([0.0, 0.5, 1.0, 1.01])
        probabilities = infection_probability(
            distances,
            source_is_superspreader=True,
            model="strong",
            r0=1.0,
            w0=0.8,
        )
        np.testing.assert_allclose(probabilities, [0.8, 0.8, 0.8, 0.0])

    def test_hub_superspreader_uses_sqrt_six_cutoff(self) -> None:
        cutoff = np.sqrt(6.0)
        distances = np.array([0.0, cutoff / 2.0, cutoff, cutoff + 0.01])
        probabilities = infection_probability(
            distances,
            source_is_superspreader=True,
            model="hub",
            r0=1.0,
            w0=1.0,
        )
        np.testing.assert_allclose(probabilities, [1.0, 0.25, 0.0, 0.0])

    def test_vertical_wrapping_is_explicitly_optional(self) -> None:
        source = np.array([5.0, 0.1])
        target = np.array([[5.0, 9.9]])
        self.assertAlmostEqual(periodic_distance(source, target, 10.0)[0], 9.8)
        self.assertAlmostEqual(
            periodic_distance(source, target, 10.0, periodic_y=True)[0],
            0.2,
        )

    def test_simulator_infects_across_vertical_boundary(self) -> None:
        simulator = MonteCarloSIRSimulator(
            model="strong",
            lambda_ss=0.5,
            positions=np.array([[1.0, 1.0], [5.0, 9.5]]),
            is_superspreader=np.array([True, False]),
            L=10.0,
            r0=1.0,
            w0=1.0,
            gamma=1.0,
            max_steps=2,
            seed=0,
        )

        result = simulator.run()

        self.assertEqual(result.total_infected, 2)
        self.assertAlmostEqual(result.unwrapped_positions[1, 1], -0.5)


if __name__ == "__main__":
    unittest.main()
