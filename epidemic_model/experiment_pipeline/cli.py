"""Command-line orchestration for the experiment pipeline."""

import argparse
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from epidemic_model.config import DEFAULT_MAX_STEPS
from epidemic_model.simulator import SimulationResult

from .files import prepare_output_dir, write_csv
from .metrics import run_configured_simulation, run_simulation
from .plotting import (
    regenerate_plots_from_outputs,
    route_selection_output_rows,
    select_representative_routes,
    write_selected_plots,
)
from .records import SimulationConfig, SimulationMetrics
from .scenarios import (
    density_sweep_configs,
    epidemic_curve_configs,
    fixed_density_baseline_configs,
    route_and_secondary_configs,
    sars_curve_comparison_configs,
    superspreader_fraction_sweep_configs,
)
from .settings import (
    DEFAULT_SEED_INDICES,
    SPREAD_ANALYSIS_DENSITY,
    SPREAD_ANALYSIS_POPULATION,
)
from .tables import (
    critical_density_estimate_rows,
    critical_density_theory_rows,
    grouped_metric_rows,
    infection_curve_rows,
    infection_front_rows,
    infection_kernel_rows,
    sars_case_time_series_rows,
    sars_secondary_case_rows,
    secondary_case_distribution_rows,
    simulation_summary_rows,
)


EXPERIMENT_BUILDERS = {
    "baseline": fixed_density_baseline_configs,
    "percolation": density_sweep_configs,
    "epidemic-curves": epidemic_curve_configs,
    "routes": route_and_secondary_configs,
    "sars": sars_curve_comparison_configs,
    "sensitivity": superspreader_fraction_sweep_configs,
}
EXPERIMENT_CHOICES = ("all", *EXPERIMENT_BUILDERS)


def normalize_experiment_selection(selected: list[str] | None) -> tuple[str, ...]:
    if not selected or "all" in selected:
        return tuple(EXPERIMENT_BUILDERS)

    normalized: list[str] = []
    for name in selected:
        if name not in EXPERIMENT_BUILDERS:
            raise ValueError(f"unknown experiment group: {name}")
        if name not in normalized:
            normalized.append(name)
    return tuple(normalized)


def default_output_dir(selected_experiments: tuple[str, ...]) -> Path:
    if selected_experiments == tuple(EXPERIMENT_BUILDERS):
        return Path("results")
    if len(selected_experiments) == 1:
        return Path("results") / selected_experiments[0]
    return Path("results") / ("selected_" + "_".join(selected_experiments))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for generated CSV, plot, and report files.",
    )
    parser.add_argument(
        "--experiment",
        action="append",
        choices=EXPERIMENT_CHOICES,
        help=(
            "Experiment group to run. Use multiple --experiment flags to combine groups. "
            "Choices: all, baseline, percolation, epidemic-curves, routes, sars, sensitivity."
        ),
    )
    parser.add_argument(
        "--list-experiments",
        action="store_true",
        help="Print available experiment groups and exit.",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        default=len(DEFAULT_SEED_INDICES),
        help="Number of random seeds per parameter setting.",
    )
    parser.add_argument(
        "--seed-offset",
        dest="seed_offset",
        type=int,
        default=0,
        help="First random seed to use; increase this to generate a fresh reproducible result set.",
    )
    parser.add_argument(
        "--max-steps",
        dest="max_time_steps",
        type=int,
        default=DEFAULT_MAX_STEPS,
        help="Maximum timesteps per simulation.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 1),
        help="Number of worker processes for simulations.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=20,
        help="Number of simulations sent to each worker chunk when --jobs > 1.",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Regenerate plots from existing CSV files without rerunning the full Monte Carlo batch.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.list_experiments:
        print("available experiments:")
        for name in EXPERIMENT_CHOICES:
            print(f"  {name}")
        return

    if args.seeds <= 0:
        raise ValueError("--seeds must be positive")
    if args.seed_offset < 0:
        raise ValueError("--seed-offset must be non-negative")
    if args.max_time_steps <= 0:
        raise ValueError("--max-steps must be positive")
    if args.jobs <= 0:
        raise ValueError("--jobs must be positive")
    if args.chunksize <= 0:
        raise ValueError("--chunksize must be positive")

    selected_experiments = normalize_experiment_selection(args.experiment)
    selected_set = set(selected_experiments)
    random_seeds = tuple(range(args.seed_offset, args.seed_offset + args.seeds))
    output_dir = args.output_dir or default_output_dir(selected_experiments)
    plot_experiments = selected_set - {"baseline"}
    if args.plot_only:
        if plot_experiments:
            regenerate_plots_from_outputs(output_dir, plot_experiments)
            print(f"regenerated plots in {output_dir / 'plots'}")
        else:
            print("baseline has no plots to regenerate")
        return

    plots_dir = prepare_output_dir(output_dir)

    epidemic_curve_results: list[tuple[SimulationConfig, SimulationResult]] = []
    secondary_route_results: list[tuple[SimulationConfig, SimulationResult]] = []
    sars_comparison_results: list[tuple[SimulationConfig, SimulationResult]] = []
    front_results: list[tuple[SimulationConfig, SimulationResult]] = []
    route_results: list[tuple[SimulationConfig, SimulationResult]] = []
    all_metrics: list[SimulationMetrics] = []

    all_configs: list[SimulationConfig] = []
    for experiment_name in selected_experiments:
        all_configs.extend(EXPERIMENT_BUILDERS[experiment_name](random_seeds, args.max_time_steps))

    def collect_result(index: int, config: SimulationConfig, result: SimulationResult, metrics: SimulationMetrics) -> None:
        all_metrics.append(metrics)
        if config.scenario_name == "epidemic_curve_comparison":
            epidemic_curve_results.append((config, result))
        if config.scenario_name == "route_and_secondary_distribution":
            secondary_route_results.append((config, result))
            route_results.append((config, result))
        if config.scenario_name == "sars_epidemic_comparison":
            sars_comparison_results.append((config, result))
        if (
            config.scenario_name == "superspreader_fraction_sweep"
            and config.population_size == SPREAD_ANALYSIS_POPULATION
            and abs(config.population_density - SPREAD_ANALYSIS_DENSITY) < 1e-9
        ):
            front_results.append((config, result))
        if index % 500 == 0 or index == total_configs:
            print(f"completed {index}/{total_configs} simulations", flush=True)

    total_configs = len(all_configs)
    if args.jobs == 1:
        for index, config in enumerate(all_configs, start=1):
            result, metrics = run_simulation(config)
            collect_result(index, config, result, metrics)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            completed_runs = executor.map(run_configured_simulation, all_configs, chunksize=args.chunksize)
            for index, (config, result, metrics) in enumerate(completed_runs, start=1):
                collect_result(index, config, result, metrics)

    summary_rows = simulation_summary_rows(all_metrics)
    write_csv(output_dir / "summary_metrics.csv", summary_rows)

    baseline_rows: list[dict[str, object]] | None = None
    percolation_rows: list[dict[str, object]] | None = None
    critical_rows: list[dict[str, object]] | None = None
    critical_reference_rows: list[dict[str, object]] | None = None
    front_rows: list[dict[str, object]] | None = None
    curves_rows: list[dict[str, object]] | None = None
    sars_model_curve_rows: list[dict[str, object]] | None = None
    secondary_rows: list[dict[str, object]] | None = None
    infection_rows: list[dict[str, object]] | None = None
    sars_secondary_rows: list[dict[str, object]] | None = None
    sars_epicurve_rows: list[dict[str, object]] | None = None
    sensitivity_rows: list[dict[str, object]] | None = None
    selected_route_results: list[tuple[SimulationConfig, SimulationResult]] = []

    if "baseline" in selected_set:
        baseline_metrics = [metric for metric in all_metrics if metric.scenario_name == "fixed_density_baseline"]
        baseline_rows = grouped_metric_rows(
            baseline_metrics,
            ("transmission_model", "superspreader_fraction", "population_size", "population_density"),
        )
        write_csv(output_dir / "baseline_summary.csv", baseline_rows)

    if "percolation" in selected_set:
        percolation_metrics = [metric for metric in all_metrics if metric.scenario_name == "percolation_density_sweep"]
        percolation_rows = grouped_metric_rows(
            percolation_metrics,
            ("transmission_model", "superspreader_fraction", "population_size", "population_density"),
        )
        write_csv(output_dir / "percolation_probability.csv", percolation_rows)

        critical_rows = critical_density_estimate_rows(percolation_rows)
        write_csv(output_dir / "critical_density.csv", critical_rows)
        critical_reference_rows = critical_density_theory_rows()
        write_csv(output_dir / "critical_density_reference_curves.csv", critical_reference_rows)

        infection_rows = infection_kernel_rows()
        write_csv(output_dir / "infection_probability_functions.csv", infection_rows)

    if "sensitivity" in selected_set:
        propagation_rows = grouped_metric_rows(
            all_metrics,
            ("scenario_name", "transmission_model", "superspreader_fraction", "population_size", "population_density"),
        )
        write_csv(output_dir / "front_speed.csv", propagation_rows)

        front_rows = infection_front_rows(front_results, args.max_time_steps)
        write_csv(output_dir / "front_distance.csv", front_rows)

        sensitivity_metrics = [metric for metric in all_metrics if metric.scenario_name == "superspreader_fraction_sweep"]
        sensitivity_rows = grouped_metric_rows(
            sensitivity_metrics,
            ("scenario_name", "transmission_model", "superspreader_fraction", "population_size", "population_density"),
        )
        write_csv(output_dir / "sensitivity_summary.csv", sensitivity_rows)

    if "epidemic-curves" in selected_set:
        curves_rows = infection_curve_rows(epidemic_curve_results, args.max_time_steps)
        write_csv(output_dir / "epidemic_curves.csv", curves_rows)

    if "routes" in selected_set:
        secondary_rows = secondary_case_distribution_rows(secondary_route_results)
        write_csv(output_dir / "secondary_distribution.csv", secondary_rows)

        selected_route_results = select_representative_routes(route_results)
        write_csv(output_dir / "route_plot_selections.csv", route_selection_output_rows(selected_route_results))

    if "sars" in selected_set:
        sars_model_curve_rows = infection_curve_rows(sars_comparison_results, args.max_time_steps)
        write_csv(output_dir / "sars_epidemic_model_curves.csv", sars_model_curve_rows)

        sars_secondary_rows = sars_secondary_case_rows()
        write_csv(output_dir / "sars_singapore_secondary_patients.csv", sars_secondary_rows)

        sars_epicurve_rows = sars_case_time_series_rows()
        write_csv(output_dir / "sars_singapore_epidemic_curve.csv", sars_epicurve_rows)

    write_selected_plots(
        plots_dir,
        selected_experiments=plot_experiments,
        infection_rows=infection_rows,
        percolation_rows=percolation_rows,
        critical_rows=critical_rows,
        critical_reference_rows=critical_reference_rows,
        front_rows=front_rows,
        sensitivity_rows=sensitivity_rows,
        curves_rows=curves_rows,
        secondary_rows=secondary_rows,
        sars_secondary_rows=sars_secondary_rows,
        sars_epicurve_rows=sars_epicurve_rows,
        sars_model_curve_rows=sars_model_curve_rows,
        route_results=selected_route_results,
    )

    print(f"wrote clean results to {output_dir}")
