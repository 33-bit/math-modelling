import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from epidemic_model.experiment_pipeline import (
    SimulationConfig,
    DENSITY_SWEEP_POPULATIONS,
    SimulationMetrics,
    critical_density_theory_rows,
    critical_density_estimate_rows,
    top_band_reached,
    simulation_summary_rows,
    route_configs_from_saved_rows,
    sars_case_time_series_rows,
    sars_secondary_case_rows,
)
from epidemic_model.experiment_pipeline.files import normalize_csv_row
from epidemic_model.experiment_pipeline.cli import (
    default_output_dir,
    normalize_experiment_selection,
)


class ExperimentDataTests(unittest.TestCase):
    def test_sars_secondary_distribution_accounts_for_all_cases(self) -> None:
        rows = sars_secondary_case_rows()
        frequencies = {int(row["secondary_cases"]): int(row["frequency"]) for row in rows}

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
                "transmission_model": transmission_model,
                "superspreader_fraction": "1.0000",
                "population_density": population_density,
                "percolation_probability": probability,
            }
            for transmission_model in ("strong", "hub")
            for population_density, probability in ((0.75, 0.2), (1.0, 0.4), (1.25, 0.8))
        ]

        critical_rows = critical_density_estimate_rows(rows)

        self.assertEqual(len(critical_rows), 2)
        for row in critical_rows:
            self.assertAlmostEqual(float(row["critical_density"]), 1.0625)

    def test_critical_density_reference_equations(self) -> None:
        rows = critical_density_theory_rows(points=6)
        indexed = {(row["transmission_model"], float(row["superspreader_fraction"])): row for row in rows}

        self.assertEqual(len(rows), 12)
        for transmission_model, critical_r0 in (("strong", 4.5), ("hub", 3.2)):
            at_zero = indexed[(transmission_model, 0.0)]
            at_one = indexed[(transmission_model, 1.0)]
            at_point_four = indexed[(transmission_model, 0.4)]

            self.assertAlmostEqual(float(at_zero["reproduction_factor"]), 1.0 / 6.0)
            self.assertAlmostEqual(float(at_one["reproduction_factor"]), 1.0)
            self.assertAlmostEqual(float(at_point_four["reproduction_factor"]), 0.5)
            self.assertAlmostEqual(float(at_zero["critical_reproduction_number"]), critical_r0)
            self.assertAlmostEqual(
                float(at_zero["normalized_critical_density"]),
                6.0 * critical_r0,
                places=5,
            )
            self.assertAlmostEqual(
                float(at_point_four["normalized_critical_density"]),
                2.0 * critical_r0,
                places=5,
            )
            self.assertAlmostEqual(float(at_one["normalized_critical_density"]), critical_r0)

    def test_percolation_sweep_covers_hub_reference_density(self) -> None:
        self.assertEqual(DENSITY_SWEEP_POPULATIONS[:4], (50, 75, 100, 125))
        normalized_densities = [n / 100.0 * 3.141592653589793 for n in DENSITY_SWEEP_POPULATIONS]
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

        self.assertFalse(top_band_reached(below, 10.0))
        self.assertTrue(top_band_reached(reaches_band, 10.0))

    def test_route_selection_output_rows_restore_exact_configs(self) -> None:
        rows = [
            {
                "scenario_name": "route_and_secondary_distribution",
                "transmission_model": transmission_model,
                "superspreader_fraction": "0.0000" if transmission_model == "normal" else "0.2000",
                "population_size": "477",
                "population_density": "4.770000",
                "domain_size": "10.000000",
                "random_seed": str(random_seed),
                "max_time_steps": "200",
                "total_cases": "100",
                "reached_top_boundary": "1",
            }
            for transmission_model, random_seed in (("normal", 100001), ("strong", 100002), ("hub", 100003))
        ]

        configs = route_configs_from_saved_rows(rows)

        self.assertEqual(
            configs,
            [
                SimulationConfig(
                    scenario_name="route_and_secondary_distribution",
                    transmission_model=transmission_model,
                    superspreader_fraction=0.0 if transmission_model == "normal" else 0.2,
                    population_size=477,
                    population_density=4.77,
                    random_seed=random_seed,
                    max_time_steps=200,
                )
                for transmission_model, random_seed in (("normal", 100001), ("strong", 100002), ("hub", 100003))
            ],
        )

    def test_legacy_csv_headers_are_normalized(self) -> None:
        row = normalize_csv_row(
            {
                "experiment": "route_and_secondary_distribution",
                "model": "hub",
                "lambda_ss": "0.2000",
                "N": "477",
                "density": "4.770000",
                "L": "10.000000",
                "seed": "100003",
                "max_steps": "200",
                "mean_attack_rate": "0.500000",
                "mean_propagation_speed": "1.250000",
                "secondary_infections": "3",
                "critical_R0": "3.200000",
                "percolated": "1",
            }
        )

        self.assertEqual(row["scenario_name"], "route_and_secondary_distribution")
        self.assertEqual(row["transmission_model"], "hub")
        self.assertEqual(row["superspreader_fraction"], "0.2000")
        self.assertEqual(row["population_size"], "477")
        self.assertEqual(row["population_density"], "4.770000")
        self.assertEqual(row["domain_size"], "10.000000")
        self.assertEqual(row["random_seed"], "100003")
        self.assertEqual(row["max_time_steps"], "200")
        self.assertEqual(row["mean_case_fraction"], "0.500000")
        self.assertEqual(row["mean_front_speed"], "1.250000")
        self.assertEqual(row["secondary_cases"], "3")
        self.assertEqual(row["critical_reproduction_number"], "3.200000")
        self.assertEqual(row["reached_top_boundary"], "1")

    def test_experiment_selection_defaults_to_all_groups(self) -> None:
        self.assertEqual(
            normalize_experiment_selection(None),
            ("baseline", "percolation", "epidemic-curves", "routes", "sars", "sensitivity"),
        )
        self.assertEqual(normalize_experiment_selection(["all"]), normalize_experiment_selection(None))

    def test_experiment_selection_preserves_requested_order_without_duplicates(self) -> None:
        self.assertEqual(
            normalize_experiment_selection(["routes", "percolation", "routes"]),
            ("routes", "percolation"),
        )

    def test_partial_experiment_default_output_does_not_target_full_results(self) -> None:
        self.assertEqual(default_output_dir(("percolation",)), Path("results/percolation"))
        self.assertEqual(
            default_output_dir(("routes", "sars")),
            Path("results/selected_routes_sars"),
        )

    def test_output_row_schemas_are_stable(self) -> None:
        metric = SimulationMetrics(
            scenario_name="test",
            transmission_model="normal",
            superspreader_fraction=0.0,
            population_size=10,
            population_density=0.1,
            domain_size=10.0,
            random_seed=1,
            total_cases=1,
            case_fraction=0.1,
            outbreak_duration=1,
            peak_new_cases=0,
            peak_active_cases=1,
            front_speed=0.0,
            mean_secondary_cases=0.0,
            max_secondary_cases=0,
            reached_top_boundary=False,
        )

        self.assertEqual(
            tuple(simulation_summary_rows([metric])[0]),
            (
                "scenario_name",
                "transmission_model",
                "superspreader_fraction",
                "population_size",
                "population_density",
                "domain_size",
                "random_seed",
                "total_cases",
                "case_fraction",
                "outbreak_duration",
                "peak_new_cases",
                "peak_active_cases",
                "front_speed",
                "mean_secondary_cases",
                "max_secondary_cases",
                "reached_top_boundary",
            ),
        )
        self.assertEqual(
            tuple(sars_secondary_case_rows()[0]),
            ("secondary_cases", "frequency", "probability", "denominator", "source"),
        )
        self.assertEqual(
            tuple(critical_density_theory_rows(points=2)[0]),
            (
                "transmission_model",
                "superspreader_fraction",
                "reproduction_factor",
                "critical_reproduction_number",
                "normalized_critical_density",
            ),
        )
        self.assertEqual(
            tuple(sars_case_time_series_rows()[0]),
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
