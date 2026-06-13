import unittest
from types import SimpleNamespace

import numpy as np

from experiments import (
    ExperimentConfig,
    PERCOLATION_N_SWEEP,
    RunMetrics,
    critical_density_rows,
    has_percolated_to_top,
    metric_rows,
    route_configs_from_rows,
    sars_epidemic_curve_rows,
    sars_secondary_patient_rows,
)


class ExperimentDataTests(unittest.TestCase):
    def test_sars_secondary_distribution_accounts_for_all_cases(self) -> None:
        rows = sars_secondary_patient_rows()
        frequencies = {int(row["secondary_infections"]): int(row["frequency"]) for row in rows}

        self.assertEqual(
            frequencies,
            {0: 162, 1: 19, 2: 8, 3: 6, 7: 1, 12: 1, 21: 1, 23: 2, 40: 1},
        )
        self.assertEqual(sum(frequencies.values()), 201)
        self.assertAlmostEqual(sum(float(row["probability"]) for row in rows), 1.0, places=5)
        self.assertTrue(all(int(row["denominator"]) == 201 for row in rows))

    def test_critical_density_interpolates_half_probability_crossing(self) -> None:
        rows = [
            {
                "model": model,
                "lambda_ss": "1.0000",
                "density": density,
                "percolation_probability": probability,
            }
            for model in ("strong", "hub")
            for density, probability in ((0.75, 0.2), (1.0, 0.4), (1.25, 0.8))
        ]

        critical_rows = critical_density_rows(rows)

        self.assertEqual(len(critical_rows), 2)
        for row in critical_rows:
            self.assertAlmostEqual(float(row["critical_density"]), 1.0625)

    def test_percolation_sweep_covers_hub_reference_density(self) -> None:
        self.assertEqual(PERCOLATION_N_SWEEP[:4], (50, 75, 100, 125))
        normalized_densities = [n / 100.0 * 3.141592653589793 for n in PERCOLATION_N_SWEEP]
        self.assertLess(min(normalized_densities), 3.2)
        self.assertGreater(max(normalized_densities), 3.2)

    def test_top_band_percolation_rule(self) -> None:
        below = SimpleNamespace(
            infected_time=np.array([0.0, 1.0]),
            positions=np.array([[5.0, 0.0], [5.0, 8.999]]),
        )
        reaches_band = SimpleNamespace(
            infected_time=np.array([0.0, 1.0]),
            positions=np.array([[5.0, 0.0], [5.0, 9.0]]),
        )

        self.assertFalse(has_percolated_to_top(below, 10.0))
        self.assertTrue(has_percolated_to_top(reaches_band, 10.0))

    def test_route_selection_rows_restore_exact_configs(self) -> None:
        rows = [
            {
                "experiment": "paper_fig_9_13_networks",
                "model": model,
                "lambda_ss": "0.0000" if model == "normal" else "0.2000",
                "N": "477",
                "density": "4.770000",
                "L": "10.000000",
                "seed": str(seed),
                "max_steps": "200",
                "total_infected": "100",
                "percolated": "1",
            }
            for model, seed in (("normal", 100001), ("strong", 100002), ("hub", 100003))
        ]

        configs = route_configs_from_rows(rows)

        self.assertEqual(
            configs,
            [
                ExperimentConfig(
                    experiment="paper_fig_9_13_networks",
                    model=model,
                    lambda_ss=0.0 if model == "normal" else 0.2,
                    N=477,
                    density=4.77,
                    seed=seed,
                    max_steps=200,
                )
                for model, seed in (("normal", 100001), ("strong", 100002), ("hub", 100003))
            ],
        )

    def test_output_row_schemas_are_stable(self) -> None:
        metric = RunMetrics(
            experiment="test",
            model="normal",
            lambda_ss=0.0,
            N=10,
            density=0.1,
            L=10.0,
            seed=1,
            total_infected=1,
            attack_rate=0.1,
            duration=1,
            peak_new_infections=0,
            peak_active_infected=1,
            propagation_speed=0.0,
            mean_secondary_infections=0.0,
            max_secondary_infections=0,
            percolated=False,
        )

        self.assertEqual(
            tuple(metric_rows([metric])[0]),
            (
                "experiment",
                "model",
                "lambda_ss",
                "N",
                "density",
                "L",
                "seed",
                "total_infected",
                "attack_rate",
                "duration",
                "peak_new_infections",
                "peak_active_infected",
                "propagation_speed",
                "mean_secondary_infections",
                "max_secondary_infections",
                "percolated",
            ),
        )
        self.assertEqual(
            tuple(sars_secondary_patient_rows()[0]),
            ("secondary_infections", "frequency", "probability", "denominator", "source"),
        )
        self.assertEqual(
            tuple(sars_epidemic_curve_rows()[0]),
            (
                "time_step",
                "days_since_start",
                "window_start",
                "window_end",
                "new_cases",
                "total_cases",
                "source",
            ),
        )


if __name__ == "__main__":
    unittest.main()
