"""Paper-style experiments and sensitivity analysis.

This script runs reproducible Monte Carlo batches for the paper's superspreader SIR
setup, writes clean CSV files, generates plots, and creates a concise report.
"""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path
from statistics import mean
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from config import DEFAULT_MAX_STEPS, DEFAULT_N, DEFAULT_R0, DEFAULT_W0
from simulator import MonteCarloSIRSimulator, SimulationResult, infection_probability, periodic_distance


DEFAULT_SEEDS = tuple(range(1000))
PERCOLATION_TOP_MARGIN = DEFAULT_R0
LAMBDA_SWEEP = (0.00, 0.20, 0.40, 0.60, 0.80, 1.00)
PERCOLATION_N_SWEEP = (50, 75, 100, 125, *range(150, 901, 50))
PLOT_DPI = 300
PAPER_FIGSIZE = (4.8, 3.35)
PAPER_ROUTE_FIGSIZE = (3.8, 3.8)
PAPER_LINEWIDTH = 1.0
PAPER_MARKERSIZE = 4.2
MODEL_COLORS = {
    "normal": "#6B7280",
    "strong": "#2563EB",
    "hub": "#E11D48",
}
INFECTION_SOURCE_COLORS = {
    "normal": "#6B7280",
    "superspreader": "#E11D48",
}
LAMBDA_COLORS = {
    0.0: "#6B7280",
    0.2: "#2563EB",
    0.4: "#059669",
    0.6: "#D97706",
    0.8: "#7C3AED",
    1.0: "#DC2626",
}
ROUTE_COLORS = {
    "edge": "#94A3B8",
    "susceptible_normal": "#CBD5E1",
    "susceptible_superspreader": "#F59E0B",
    "infected_normal": "#2563EB",
    "infected_superspreader": "#E11D48",
    "initial": "#111827",
}
SARS_BAR_COLOR = "#FDE68A"
LAMBDA_MARKERS = {
    0.0: "o",
    0.2: "s",
    0.4: "^",
    0.6: "D",
    0.8: "v",
    1.0: "*",
}
PAPER_BOX_L = 10.0 * DEFAULT_R0
DENSITY_SCALE = np.pi * DEFAULT_R0**2
PAPER_FIG_6_8_NORMALIZED_DENSITY = 20.0
PAPER_FIG_6_8_N = int(round(PAPER_FIG_6_8_NORMALIZED_DENSITY * PAPER_BOX_L**2 / DENSITY_SCALE))
PAPER_FIG_6_8_DENSITY = PAPER_FIG_6_8_N / PAPER_BOX_L**2
PAPER_FIG_9_13_NORMALIZED_DENSITY = 15.0
PAPER_FIG_9_13_N = int(round(PAPER_FIG_9_13_NORMALIZED_DENSITY * PAPER_BOX_L**2 / DENSITY_SCALE))
PAPER_FIG_9_13_DENSITY = PAPER_FIG_9_13_N / PAPER_BOX_L**2
SARS_COMPARISON_NORMALIZED_DENSITY = 15.0
SARS_COMPARISON_N = DEFAULT_N
SARS_COMPARISON_DENSITY = SARS_COMPARISON_N / PAPER_BOX_L**2
PAPER_SPREAD_LAMBDA = 0.20
SARS_COMPARISON_LAMBDA = 0.40
SARS_TIMESTEP_DAYS = 6
STRONG_CRITICAL_R0_REFERENCE = 4.5
HUB_CRITICAL_R0_REFERENCE = 3.2
SARS_SECONDARY_TOTAL_PROBABLE_CASES = 201
PLOT_FILENAMES = {
    "infection_probability_strong": "fig01_infection_probability_strong.png",
    "infection_probability_hub": "fig02_infection_probability_hub.png",
    "percolation_probability_strong": "fig03_percolation_probability_strong.png",
    "percolation_probability_hub": "fig04_percolation_probability_hub.png",
    "critical_density": "fig05_critical_density.png",
    "front_distance_strong": "fig06_front_distance_strong.png",
    "velocity_vs_lambda": "fig07_velocity_vs_lambda.png",
    "epidemic_curves": "fig08_epidemic_curves.png",
    "infection_route_strong": "fig09_infection_route_strong.png",
    "infection_route_hub": "fig10_infection_route_hub.png",
    "infection_route_normal": "fig11_infection_route_normal.png",
    "secondary_distribution_normal": "fig12_secondary_distribution_normal.png",
    "secondary_distribution_superspreaders": "fig13_secondary_distribution_superspreaders.png",
    "sars_secondary_patients": "fig14_sars_secondary_patients.png",
    "sars_epidemic_curve_comparison": "fig15_sars_epidemic_curve_comparison.png",
    "percolation_probability": "fig_extra_percolation_probability_comparison.png",
    "front_distance_hub": "fig_extra_front_distance_hub.png",
    "secondary_distribution": "fig_extra_secondary_distribution_comparison.png",
    "sensitivity_lambda_attack_rate": "fig_extra_sensitivity_lambda_attack_rate.png",
}
SARS_SECONDARY_PATIENT_FREQUENCIES = (
    (0, 162, "CDC MMWR: 162 probable SARS cases had no evidence of transmission"),
    (1, 19, "Digitized from CDC MMWR Figure 3"),
    (2, 8, "Digitized from CDC MMWR Figure 3"),
    (3, 6, "Digitized from CDC MMWR Figure 3"),
    (7, 1, "Digitized from CDC MMWR Figure 3"),
    (12, 1, "CDC MMWR Case 5 direct probable SARS links"),
    (21, 1, "CDC MMWR Case 1 direct probable SARS links"),
    (23, 2, "CDC MMWR Cases 2 and 3 direct probable SARS links"),
    (40, 1, "CDC MMWR Case 4 direct probable SARS links"),
)
SARS_EPICURVE_6_DAY_COUNTS = (
    0,
    0,
    1,
    2,
    17,
    42,
    36,
    20,
    44,
    29,
    20,
    15,
    8,
    4,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Parameters that define one simulation replicate."""

    experiment: str
    model: str
    lambda_ss: float
    N: int
    density: float
    seed: int
    max_steps: int = DEFAULT_MAX_STEPS

    @property
    def L(self) -> float:
        return float(np.sqrt(self.N / self.density))


@dataclass(frozen=True)
class RunMetrics:
    """Compact metrics collected from one simulation replicate."""

    experiment: str
    model: str
    lambda_ss: float
    N: int
    density: float
    L: float
    seed: int
    total_infected: int
    attack_rate: float
    duration: int
    peak_new_infections: int
    peak_active_infected: int
    propagation_speed: float
    mean_secondary_infections: float
    max_secondary_infections: int
    percolated: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory for generated CSV, plot, and report files.",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        default=len(DEFAULT_SEEDS),
        help="Number of random seeds per parameter setting.",
    )
    parser.add_argument(
        "--seed-offset",
        type=int,
        default=0,
        help="First random seed to use; increase this to generate a fresh reproducible result set.",
    )
    parser.add_argument(
        "--max-steps",
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


def active_infected_curve(result: SimulationResult, max_steps: int) -> list[int]:
    """Return active infected counts for t = 0..max_steps."""
    curve: list[int] = []
    recovered_time = result.recovered_time
    infected_time = result.infected_time

    for time in range(max_steps + 1):
        infected = np.isfinite(infected_time) & (infected_time <= time)
        not_recovered = np.isnan(recovered_time) | (recovered_time > time)
        curve.append(int(np.count_nonzero(infected & not_recovered)))
    return curve


def cumulative_infected_curve(result: SimulationResult, max_steps: int) -> list[int]:
    """Return cumulative infected counts for t = 0..max_steps."""
    curve: list[int] = []
    infected_time = result.infected_time
    for time in range(max_steps + 1):
        curve.append(int(np.count_nonzero(np.isfinite(infected_time) & (infected_time <= time))))
    return curve


def new_infections_curve(result: SimulationResult, max_steps: int) -> list[int]:
    """Return new infections for t = 0..max_steps."""
    curve = [1]
    curve.extend(int(value) for value in result.new_infections_per_step)
    if len(curve) < max_steps + 1:
        curve.extend([0] * (max_steps + 1 - len(curve)))
    return curve[: max_steps + 1]


def front_distance_curve(result: SimulationResult, L: float, max_steps: int) -> list[float]:
    """Return the furthest infected distance from the initial case over time."""
    curve: list[float] = []
    source_position = result.positions[0]
    infected_time = result.infected_time
    for time in range(max_steps + 1):
        infected_ids = np.flatnonzero(np.isfinite(infected_time) & (infected_time <= time))
        if infected_ids.size == 0:
            curve.append(0.0)
            continue
        distances = periodic_distance(source_position, result.positions[infected_ids], L)
        curve.append(float(np.max(distances)))
    return curve


def has_percolated_to_top(result: SimulationResult, L: float) -> bool:
    """Return whether infection reaches the top boundary band."""
    infected_ids = np.flatnonzero(np.isfinite(result.infected_time))
    if infected_ids.size == 0:
        return False
    top_threshold = max(0.0, L - PERCOLATION_TOP_MARGIN)
    return bool(np.any(result.positions[infected_ids, 1] >= top_threshold))


def estimate_propagation_speed(result: SimulationResult, L: float) -> float:
    """Estimate front speed from advancing front points before the plateau."""
    infected_ids = np.flatnonzero(np.isfinite(result.infected_time))
    if infected_ids.size < 3 or result.duration < 2:
        return 0.0

    source_position = result.positions[0]
    distances = periodic_distance(source_position, result.positions[infected_ids], L)
    times = result.infected_time[infected_ids]

    front_times: list[float] = []
    front_radii: list[float] = []
    previous_front_radius = 0.0
    for time in range(1, result.duration + 1):
        reached = times <= time
        if np.any(reached):
            front_radius = float(np.max(distances[reached]))
            if front_radius > previous_front_radius + 1e-9:
                front_times.append(float(time))
                front_radii.append(front_radius)
                previous_front_radius = front_radius

    if len(front_times) < 2:
        return 0.0

    slope, _intercept = np.polyfit(front_times, front_radii, deg=1)
    return max(0.0, float(slope))


def run_simulation(config: ExperimentConfig) -> tuple[SimulationResult, RunMetrics]:
    simulator = MonteCarloSIRSimulator(
        N=config.N,
        model=config.model,
        lambda_ss=config.lambda_ss,
        L=config.L,
        max_steps=config.max_steps,
        seed=config.seed,
    )
    result = simulator.run()
    active_curve = active_infected_curve(result, config.max_steps)
    secondary_counts = result.secondary_counts[np.isfinite(result.infected_time)]
    attack_rate = result.total_infected / config.N

    metrics = RunMetrics(
        experiment=config.experiment,
        model=config.model,
        lambda_ss=config.lambda_ss,
        N=config.N,
        density=config.density,
        L=config.L,
        seed=config.seed,
        total_infected=result.total_infected,
        attack_rate=attack_rate,
        duration=result.duration,
        peak_new_infections=max(result.new_infections_per_step, default=0),
        peak_active_infected=max(active_curve, default=0),
        propagation_speed=estimate_propagation_speed(result, config.L),
        mean_secondary_infections=float(np.mean(secondary_counts)) if secondary_counts.size else 0.0,
        max_secondary_infections=int(np.max(secondary_counts)) if secondary_counts.size else 0,
        percolated=has_percolated_to_top(result, config.L),
    )
    return result, metrics


def run_config(config: ExperimentConfig) -> tuple[ExperimentConfig, SimulationResult, RunMetrics]:
    result, metrics = run_simulation(config)
    return config, result, metrics


def baseline_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    for model, lambda_ss in (("normal", 0.0), ("strong", PAPER_SPREAD_LAMBDA), ("hub", PAPER_SPREAD_LAMBDA)):
        for seed in seeds:
            configs.append(
                ExperimentConfig(
                    experiment="baseline_curves",
                    model=model,
                    lambda_ss=lambda_ss,
                    N=PAPER_FIG_6_8_N,
                    density=PAPER_FIG_6_8_DENSITY,
                    seed=seed,
                    max_steps=max_steps,
                )
            )
    return configs


def percolation_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    for N in PERCOLATION_N_SWEEP:
        density = N / (10.0 * DEFAULT_R0) ** 2
        for model in ("strong", "hub"):
            for lambda_ss in LAMBDA_SWEEP:
                for seed in seeds:
                    configs.append(
                        ExperimentConfig(
                            experiment="percolation_density_sweep",
                            model=model,
                            lambda_ss=lambda_ss,
                            N=N,
                            density=density,
                            seed=seed,
                            max_steps=max_steps,
                        )
                    )

        for seed in seeds:
            configs.append(
                ExperimentConfig(
                    experiment="percolation_density_sweep",
                    model="normal",
                    lambda_ss=0.0,
                    N=N,
                    density=density,
                    seed=seed,
                    max_steps=max_steps,
                )
            )
    return configs


def paper_curve_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    for model, lambda_ss in (("normal", 0.0), ("strong", PAPER_SPREAD_LAMBDA), ("hub", PAPER_SPREAD_LAMBDA)):
        for seed in seeds:
            configs.append(
                ExperimentConfig(
                    experiment="paper_fig_8_curves",
                    model=model,
                    lambda_ss=lambda_ss,
                    N=PAPER_FIG_6_8_N,
                    density=PAPER_FIG_6_8_DENSITY,
                    seed=seed,
                    max_steps=max_steps,
                )
            )
    return configs


def secondary_route_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    for model, lambda_ss in (("normal", 0.0), ("strong", PAPER_SPREAD_LAMBDA), ("hub", PAPER_SPREAD_LAMBDA)):
        for seed in seeds:
            configs.append(
                ExperimentConfig(
                    experiment="paper_fig_9_13_networks",
                    model=model,
                    lambda_ss=lambda_ss,
                    N=PAPER_FIG_9_13_N,
                    density=PAPER_FIG_9_13_DENSITY,
                    seed=seed,
                    max_steps=max_steps,
                )
            )
    return configs


def sars_comparison_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    for model, lambda_ss in (("normal", 0.0), ("strong", SARS_COMPARISON_LAMBDA), ("hub", SARS_COMPARISON_LAMBDA)):
        for seed in seeds:
            configs.append(
                ExperimentConfig(
                    experiment="sars_epidemic_comparison",
                    model=model,
                    lambda_ss=lambda_ss,
                    N=SARS_COMPARISON_N,
                    density=SARS_COMPARISON_DENSITY,
                    seed=seed,
                    max_steps=max_steps,
                )
            )
    return configs


def sensitivity_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    for model in ("strong", "hub"):
        for lambda_ss in LAMBDA_SWEEP:
            for seed in seeds:
                configs.append(
                    ExperimentConfig(
                        experiment="paper_lambda_sweep",
                        model=model,
                        lambda_ss=lambda_ss,
                        N=PAPER_FIG_6_8_N,
                        density=PAPER_FIG_6_8_DENSITY,
                        seed=seed,
                        max_steps=max_steps,
                    )
                )
    return configs


def metric_rows(metrics: Iterable[RunMetrics]) -> list[dict[str, object]]:
    return [
        {
            "experiment": metric.experiment,
            "model": metric.model,
            "lambda_ss": f"{metric.lambda_ss:.4f}",
            "N": metric.N,
            "density": f"{metric.density:.4f}",
            "L": f"{metric.L:.4f}",
            "seed": metric.seed,
            "total_infected": metric.total_infected,
            "attack_rate": f"{metric.attack_rate:.6f}",
            "duration": metric.duration,
            "peak_new_infections": metric.peak_new_infections,
            "peak_active_infected": metric.peak_active_infected,
            "propagation_speed": f"{metric.propagation_speed:.6f}",
            "mean_secondary_infections": f"{metric.mean_secondary_infections:.6f}",
            "max_secondary_infections": metric.max_secondary_infections,
            "percolated": int(metric.percolated),
        }
        for metric in metrics
    ]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write for {path}")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        raise FileNotFoundError(f"missing required CSV: {path}")
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


def infection_probability_rows(points: int = 401) -> list[dict[str, object]]:
    max_radius = np.sqrt(6.0) * DEFAULT_R0
    distances = np.linspace(0.0, max_radius, points)
    rows: list[dict[str, object]] = []
    for model in ("strong", "hub"):
        for source_label, source_is_superspreader in (("normal", False), ("superspreader", True)):
            probabilities = infection_probability(
                distances,
                source_is_superspreader=source_is_superspreader,
                model=model,
                r0=DEFAULT_R0,
                w0=DEFAULT_W0,
            )
            for distance, probability in zip(distances, probabilities):
                rows.append(
                    {
                        "model": model,
                        "source_type": source_label,
                        "distance_over_r0": f"{distance / DEFAULT_R0:.6f}",
                        "probability_over_w0": f"{probability / DEFAULT_W0:.6f}",
                    }
                )
    return rows


def sars_secondary_patient_rows() -> list[dict[str, object]]:
    total_frequency = sum(frequency for _count, frequency, _source in SARS_SECONDARY_PATIENT_FREQUENCIES)
    if total_frequency != SARS_SECONDARY_TOTAL_PROBABLE_CASES:
        raise ValueError(
            "SARS secondary-patient frequencies must account for all "
            f"{SARS_SECONDARY_TOTAL_PROBABLE_CASES} probable cases; got {total_frequency}"
        )

    rows: list[dict[str, object]] = []
    for secondary_count, frequency, source in SARS_SECONDARY_PATIENT_FREQUENCIES:
        rows.append(
            {
                "secondary_infections": secondary_count,
                "frequency": frequency,
                "probability": f"{frequency / SARS_SECONDARY_TOTAL_PROBABLE_CASES:.12f}",
                "denominator": SARS_SECONDARY_TOTAL_PROBABLE_CASES,
                "source": source,
            }
        )
    return rows


def sars_epidemic_curve_rows() -> list[dict[str, object]]:
    start = date(2003, 2, 13)
    total_cases = sum(SARS_EPICURVE_6_DAY_COUNTS)
    rows: list[dict[str, object]] = []
    for time_step, cases in enumerate(SARS_EPICURVE_6_DAY_COUNTS):
        window_start = start + timedelta(days=SARS_TIMESTEP_DAYS * time_step)
        window_end = window_start + timedelta(days=SARS_TIMESTEP_DAYS - 1)
        rows.append(
            {
                "time_step": time_step,
                "days_since_start": SARS_TIMESTEP_DAYS * time_step,
                "window_start": window_start.isoformat(),
                "window_end": window_end.isoformat(),
                "new_cases": cases,
                "total_cases": total_cases,
                "source": "Approximate 6-day bins digitized from the published Singapore SARS epidemic curve",
            }
        )
    return rows


def group_metrics(metrics: Iterable[RunMetrics], keys: tuple[str, ...]) -> list[dict[str, object]]:
    buckets: dict[tuple[object, ...], list[RunMetrics]] = {}
    for metric in metrics:
        key = tuple(getattr(metric, name) for name in keys)
        buckets.setdefault(key, []).append(metric)

    rows: list[dict[str, object]] = []
    for key, values in sorted(buckets.items()):
        row = {name: value for name, value in zip(keys, key)}
        row.update(
            {
                "runs": len(values),
                "mean_attack_rate": f"{mean(item.attack_rate for item in values):.6f}",
                "std_attack_rate": f"{np.std([item.attack_rate for item in values], ddof=0):.6f}",
                "percolation_probability": f"{mean(float(item.percolated) for item in values):.6f}",
                "mean_duration": f"{mean(item.duration for item in values):.6f}",
                "mean_peak_active_infected": f"{mean(item.peak_active_infected for item in values):.6f}",
                "mean_propagation_speed": f"{mean(item.propagation_speed for item in values):.6f}",
                "mean_secondary_infections": f"{mean(item.mean_secondary_infections for item in values):.6f}",
                "mean_max_secondary_infections": f"{mean(item.max_secondary_infections for item in values):.6f}",
            }
        )
        rows.append(row)
    return rows


def epidemic_curve_rows(
    results: list[tuple[ExperimentConfig, SimulationResult]],
    max_steps: int,
) -> list[dict[str, object]]:
    buckets: dict[tuple[str, float], list[SimulationResult]] = {}
    for config, result in results:
        buckets.setdefault((config.model, config.lambda_ss), []).append(result)

    rows: list[dict[str, object]] = []
    for (model, lambda_ss), result_group in sorted(buckets.items()):
        new_curves = np.array([new_infections_curve(result, max_steps) for result in result_group])
        active_curves = np.array([active_infected_curve(result, max_steps) for result in result_group])
        cumulative_curves = np.array([cumulative_infected_curve(result, max_steps) for result in result_group])

        for time in range(max_steps + 1):
            rows.append(
                {
                    "model": model,
                    "lambda_ss": f"{lambda_ss:.4f}",
                    "time": time,
                    "mean_new_infections": f"{np.mean(new_curves[:, time]):.6f}",
                    "mean_active_infected": f"{np.mean(active_curves[:, time]):.6f}",
                    "mean_cumulative_infected": f"{np.mean(cumulative_curves[:, time]):.6f}",
                }
            )
    return rows


def secondary_distribution_rows(
    results: list[tuple[ExperimentConfig, SimulationResult]],
) -> list[dict[str, object]]:
    buckets: dict[tuple[str, float], list[int]] = {}
    for config, result in results:
        infected_counts = result.secondary_counts[np.isfinite(result.infected_time)]
        buckets.setdefault((config.model, config.lambda_ss), []).extend(int(value) for value in infected_counts)

    rows: list[dict[str, object]] = []
    for (model, lambda_ss), counts in sorted(buckets.items()):
        total = len(counts)
        for secondary_count in range(max(counts) + 1 if counts else 1):
            frequency = counts.count(secondary_count)
            rows.append(
                {
                    "model": model,
                    "lambda_ss": f"{lambda_ss:.4f}",
                    "secondary_infections": secondary_count,
                    "frequency": frequency,
                    "probability": f"{frequency / total:.6f}" if total else "0.000000",
                }
            )
    return rows


def front_distance_rows(
    results: list[tuple[ExperimentConfig, SimulationResult]],
    max_steps: int,
) -> list[dict[str, object]]:
    buckets: dict[tuple[str, float], list[list[float]]] = {}
    for config, result in results:
        curve = front_distance_curve(result, config.L, max_steps)
        buckets.setdefault((config.model, config.lambda_ss), []).append(curve)

    rows: list[dict[str, object]] = []
    for (model, lambda_ss), curves in sorted(buckets.items()):
        curve_array = np.array(curves, dtype=float)
        for time in range(max_steps + 1):
            rows.append(
                {
                    "model": model,
                    "lambda_ss": f"{lambda_ss:.4f}",
                    "time": time,
                    "mean_front_distance": f"{np.mean(curve_array[:, time]):.6f}",
                }
            )
    return rows


def critical_density_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    critical_rows: list[dict[str, object]] = []
    for model in ("strong", "hub"):
        model_rows = [row for row in rows if row["model"] == model]
        for lambda_ss in sorted({float(row["lambda_ss"]) for row in model_rows}):
            lambda_rows = [row for row in model_rows if abs(float(row["lambda_ss"]) - lambda_ss) < 1e-9]
            ordered = sorted(lambda_rows, key=lambda row: float(row["density"]))
            critical_density: float | None = None

            previous_density = float(ordered[0]["density"])
            previous_probability = float(ordered[0]["percolation_probability"])
            for row in ordered:
                density = float(row["density"])
                probability = float(row["percolation_probability"])
                if probability >= 0.5:
                    if probability == previous_probability:
                        critical_density = density
                    else:
                        fraction = (0.5 - previous_probability) / (probability - previous_probability)
                        critical_density = previous_density + fraction * (density - previous_density)
                    break
                previous_density = density
                previous_probability = probability

            critical_rows.append(
                {
                    "model": model,
                    "lambda_ss": f"{lambda_ss:.4f}",
                    "critical_density": (
                        f"{critical_density:.6f}" if critical_density is not None else ""
                    ),
                }
            )
    return critical_rows


def critical_density_reference_rows(points: int = 201) -> list[dict[str, object]]:
    """Return the analytical Eq. 3-5 reference curves used in Fig. 5."""
    if points < 2:
        raise ValueError("points must be at least 2")

    rows: list[dict[str, object]] = []
    for lambda_ss in np.linspace(0.0, 1.0, points):
        reproduction_factor = lambda_ss + (1.0 - lambda_ss) / 6.0
        for model, critical_r0 in (
            ("strong", STRONG_CRITICAL_R0_REFERENCE),
            ("hub", HUB_CRITICAL_R0_REFERENCE),
        ):
            rows.append(
                {
                    "model": model,
                    "lambda_ss": f"{lambda_ss:.6f}",
                    "reproduction_factor": f"{reproduction_factor:.12f}",
                    "critical_R0": f"{critical_r0:.6f}",
                    "normalized_critical_density": f"{critical_r0 / reproduction_factor:.6f}",
                }
            )
    return rows


def prepare_output_dir(output_dir: Path) -> Path:
    plots_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    for pattern in ("*.csv", "*.md"):
        for path in output_dir.glob(pattern):
            path.unlink()
    for path in plots_dir.glob("*.png"):
        path.unlink()
    return plots_dir


def apply_paper_axes() -> None:
    axis = plt.gca()
    axis.grid(False)
    axis.tick_params(direction="in", top=False, right=False, length=4, width=0.8)
    for spine in axis.spines.values():
        spine.set_linewidth(0.8)


def finish_plot(
    path: Path,
    *,
    legend: bool = True,
    legend_loc: str = "best",
    legend_ncol: int = 1,
    legend_fontsize: float = 7.6,
    legend_bbox: tuple[float, float] | None = None,
) -> None:
    apply_paper_axes()
    if legend:
        plt.legend(
            loc=legend_loc,
            bbox_to_anchor=legend_bbox,
            ncol=legend_ncol,
            frameon=True,
            framealpha=0.9,
            facecolor="white",
            edgecolor="0.86",
            fontsize=legend_fontsize,
            handlelength=2.0,
            borderpad=0.45,
            labelspacing=0.35,
            columnspacing=0.9,
        )
    plt.tight_layout()
    plt.savefig(path, dpi=PLOT_DPI, bbox_inches="tight", pad_inches=0.03)
    plt.close()


def plot_infection_probability(rows: list[dict[str, object]], path: Path, model: str) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    model_rows = [row for row in rows if row["model"] == model]
    for source_type, style in (
        (
            "normal",
            {
                "linestyle": "--",
                "color": INFECTION_SOURCE_COLORS["normal"],
                "label": "Normal",
                "linewidth": PAPER_LINEWIDTH + 0.2,
            },
        ),
        (
            "superspreader",
            {
                "linestyle": "-",
                "color": INFECTION_SOURCE_COLORS["superspreader"],
                "label": "Superspreader",
                "linewidth": PAPER_LINEWIDTH + 0.2,
            },
        ),
    ):
        source_rows = [row for row in model_rows if row["source_type"] == source_type]
        distances = np.array([float(row["distance_over_r0"]) for row in source_rows])
        probabilities = np.array([float(row["probability_over_w0"]) for row in source_rows])
        plt.plot(distances, probabilities, **style)

    plt.xlabel(r"$r / r_0$")
    plt.ylabel(r"$w(r) / w_0$")
    plt.xlim(0.0, 1.0 if model == "strong" else np.sqrt(6.0))
    plt.ylim(-0.05, 1.05)
    finish_plot(path, legend_loc="upper right")


def plot_percolation(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    relevant_rows = [
        row
        for row in rows
        if row["model"] == "normal"
        or abs(float(row["lambda_ss"]) - PAPER_SPREAD_LAMBDA) < 1e-9
    ]
    style_by_model = {
        "normal": {"marker": "^", "linestyle": ":", "label": "No superspreaders"},
        "strong": {"marker": "o", "linestyle": "-", "label": "Strong"},
        "hub": {"marker": "s", "linestyle": "--", "label": "Hub"},
    }
    for model in ("normal", "strong", "hub"):
        model_rows = [row for row in relevant_rows if row["model"] == model]
        densities = np.array([float(row["density"]) for row in model_rows]) * DENSITY_SCALE
        probabilities = np.array([float(row["percolation_probability"]) for row in model_rows])
        order = np.argsort(densities)
        style = style_by_model[model]
        plt.plot(
            densities[order],
            probabilities[order],
            color=MODEL_COLORS[model],
            marker=style["marker"],
            linestyle=style["linestyle"],
            markerfacecolor="white",
            markeredgecolor=MODEL_COLORS[model],
            markersize=PAPER_MARKERSIZE,
            linewidth=PAPER_LINEWIDTH,
            label=style["label"],
        )

    plt.xlabel(r"$\rho \pi r_0^2$")
    plt.ylabel("Percolation probability")
    plt.xlim(0.0, 25.0)
    plt.ylim(-0.05, 1.05)
    finish_plot(path, legend_loc="lower right")


def plot_percolation_model(rows: list[dict[str, object]], path: Path, model: str) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    model_rows = [row for row in rows if row["model"] == model]
    for lambda_ss in sorted({float(row["lambda_ss"]) for row in model_rows}):
        lambda_rows = [row for row in model_rows if abs(float(row["lambda_ss"]) - lambda_ss) < 1e-9]
        densities = np.array([float(row["density"]) for row in lambda_rows]) * DENSITY_SCALE
        probabilities = np.array([float(row["percolation_probability"]) for row in lambda_rows])
        order = np.argsort(densities)
        marker = LAMBDA_MARKERS.get(round(lambda_ss, 1), "o")
        color = LAMBDA_COLORS.get(round(lambda_ss, 1), "#111827")
        plt.plot(
            densities[order],
            probabilities[order],
            linestyle="-",
            marker=marker,
            color=color,
            linewidth=PAPER_LINEWIDTH,
            markersize=PAPER_MARKERSIZE,
            markerfacecolor="white",
            markeredgecolor=color,
            label=rf"$\lambda$={lambda_ss:.1f}",
        )

    plt.xlabel(r"$\rho \pi r_0^2$")
    plt.ylabel("Percolation probability")
    plt.ylim(-0.05, 1.05)
    plt.xlim(0.0, 25.0)
    finish_plot(
        path,
        legend_loc="upper center",
        legend_ncol=3,
        legend_fontsize=7.0,
        legend_bbox=(0.5, -0.16),
    )


def plot_critical_density(
    rows: list[dict[str, object]],
    reference_rows: list[dict[str, object]],
    path: Path,
) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    style_by_model = {
        "strong": {"marker": "o", "label": "Strong sim."},
        "hub": {"marker": "s", "label": "Hub sim."},
    }
    for model in ("strong", "hub"):
        model_rows = [row for row in rows if row["model"] == model and str(row["critical_density"]).strip()]
        lambdas = np.array([float(row["lambda_ss"]) for row in model_rows])
        densities = np.array([float(row["critical_density"]) for row in model_rows]) * DENSITY_SCALE
        order = np.argsort(lambdas)
        style = style_by_model[model]
        plt.plot(
            lambdas[order],
            densities[order],
            linestyle="None",
            marker=style["marker"],
            color=MODEL_COLORS[model],
            markerfacecolor=MODEL_COLORS[model],
            markeredgecolor=MODEL_COLORS[model],
            markersize=PAPER_MARKERSIZE + 0.5,
            label=style["label"],
        )

    for model, linestyle in (("strong", "-"), ("hub", "--")):
        model_reference_rows = [row for row in reference_rows if row["model"] == model]
        lambdas = np.array([float(row["lambda_ss"]) for row in model_reference_rows])
        densities = np.array(
            [float(row["normalized_critical_density"]) for row in model_reference_rows]
        )
        order = np.argsort(lambdas)
        plt.plot(
            lambdas[order],
            densities[order],
            color=MODEL_COLORS[model],
            linestyle=linestyle,
            linewidth=PAPER_LINEWIDTH,
            alpha=0.85,
            label=rf"{model.title()} $R_0=R_c$",
        )

    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$\rho_c \pi r_0^2$")
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 25.0)
    finish_plot(path, legend_loc="upper right", legend_fontsize=7.0)


def plot_epidemic_curves(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    style_by_model = {
        "strong": {"marker": "o", "label": rf"Strong ($\lambda$={PAPER_SPREAD_LAMBDA:.1f})"},
        "hub": {"marker": "s", "label": rf"Hub ($\lambda$={PAPER_SPREAD_LAMBDA:.1f})"},
        "normal": {"marker": "^", "label": "No superspreaders"},
    }
    for model in ("strong", "hub", "normal"):
        model_rows = [row for row in rows if row["model"] == model]
        times = np.array([int(row["time"]) for row in model_rows])
        new_infections = np.array([float(row["mean_new_infections"]) for row in model_rows])
        style = style_by_model[model]
        plt.plot(
            times,
            new_infections,
            color=MODEL_COLORS[model],
            linewidth=PAPER_LINEWIDTH,
            marker=style["marker"],
            markersize=PAPER_MARKERSIZE,
            markerfacecolor="white",
            markeredgecolor=MODEL_COLORS[model],
            markevery=2,
            label=style["label"],
        )

    plt.xlabel("time step")
    plt.ylabel("the number of infected")
    plt.xlim(0.0, 40.0)
    finish_plot(path, legend_loc="upper right", legend_fontsize=7.4)


def plot_secondary_distribution(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    style_by_model = {
        "normal": {"marker": "^", "linestyle": ":", "label": "No superspreaders"},
        "strong": {"marker": "o", "linestyle": "-", "label": "Strong"},
        "hub": {"marker": "s", "linestyle": "--", "label": "Hub"},
    }
    for model in ("normal", "strong", "hub"):
        model_rows = [row for row in rows if row["model"] == model]
        counts = np.array([int(row["secondary_infections"]) for row in model_rows])
        probabilities = np.array([float(row["probability"]) for row in model_rows])
        mask = probabilities > 0
        style = style_by_model[model]
        plt.semilogy(
            counts[mask],
            probabilities[mask],
            color=MODEL_COLORS[model],
            linestyle=style["linestyle"],
            marker=style["marker"],
            markersize=PAPER_MARKERSIZE,
            markerfacecolor="white",
            markeredgecolor=MODEL_COLORS[model],
            linewidth=PAPER_LINEWIDTH,
            label=style["label"],
        )

    plt.xlabel("Secondary infections caused by one individual")
    plt.ylabel("Probability")
    finish_plot(path, legend_loc="upper right")


def plot_secondary_distribution_normal(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    model_rows = [row for row in rows if row["model"] == "normal"]
    counts = np.array([int(row["secondary_infections"]) for row in model_rows])
    probabilities = np.array([float(row["probability"]) for row in model_rows])
    mask = probabilities > 0
    plt.plot(
        counts[mask],
        probabilities[mask],
        marker="^",
        color=MODEL_COLORS["normal"],
        markerfacecolor="white",
        markeredgecolor=MODEL_COLORS["normal"],
        markersize=PAPER_MARKERSIZE,
        linewidth=PAPER_LINEWIDTH,
        label="No superspreaders",
    )

    plt.xlabel("the number of links")
    plt.ylabel("Probability")
    plt.xlim(0.0, 20.0)
    plt.ylim(0.0, 0.85)
    finish_plot(path, legend_loc="upper right")


def plot_secondary_distribution_superspreaders(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    for model, marker, linestyle in (("strong", "o", "-"), ("hub", "s", "--")):
        model_rows = [row for row in rows if row["model"] == model]
        counts = np.array([int(row["secondary_infections"]) for row in model_rows])
        probabilities = np.array([float(row["probability"]) for row in model_rows])
        mask = probabilities > 0
        plt.plot(
            counts[mask],
            probabilities[mask],
            color=MODEL_COLORS[model],
            linestyle=linestyle,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=MODEL_COLORS[model],
            markersize=PAPER_MARKERSIZE,
            linewidth=PAPER_LINEWIDTH,
            label="Strong" if model == "strong" else "Hub",
        )

    plt.xlabel("the number of links")
    plt.ylabel("Probability")
    plt.xlim(0.0, 20.0)
    plt.ylim(0.0, 0.85)
    finish_plot(path, legend_loc="upper right")


def plot_sars_secondary_distribution(
    empirical_rows: list[dict[str, object]],
    simulation_rows: list[dict[str, object]],
    path: Path,
) -> None:
    _ = simulation_rows
    plt.figure(figsize=PAPER_FIGSIZE)
    empirical_counts = np.array([int(row["secondary_infections"]) for row in empirical_rows])
    empirical_frequencies = np.array([int(row["frequency"]) for row in empirical_rows])
    plt.bar(
        empirical_counts,
        empirical_frequencies,
        width=0.8,
        color=SARS_BAR_COLOR,
        edgecolor="#92400E",
        hatch="//",
        linewidth=0.8,
    )

    plt.xlabel("number of direct secondary cases")
    plt.ylabel("number")
    plt.xlim(0, 40)
    plt.ylim(0, 180)
    finish_plot(path, legend=False)


def plot_sars_epidemic_comparison(
    empirical_rows: list[dict[str, object]],
    model_rows: list[dict[str, object]],
    path: Path,
) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    time_steps = np.array([int(row["time_step"]) for row in empirical_rows])
    cases = np.array([int(row["new_cases"]) for row in empirical_rows])
    plt.bar(
        time_steps,
        cases,
        width=0.85,
        color=SARS_BAR_COLOR,
        edgecolor="#92400E",
        hatch="//",
        linewidth=0.8,
        label="SARS data (digitized, approximate)",
    )

    for model, marker, linestyle in (("normal", "^", ":"), ("strong", "o", "-"), ("hub", "s", "--")):
        rows = [row for row in model_rows if row["model"] == model]
        times = np.array([int(row["time"]) for row in rows])
        new_infections = np.array([float(row["mean_new_infections"]) for row in rows])
        if model == "strong":
            label = rf"Strong ($\lambda$={SARS_COMPARISON_LAMBDA:.1f})"
        elif model == "hub":
            label = rf"Hub ($\lambda$={SARS_COMPARISON_LAMBDA:.1f})"
        else:
            label = "No superspreaders"
        plt.plot(
            times,
            new_infections,
            color=MODEL_COLORS[model],
            linestyle=linestyle,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=MODEL_COLORS[model],
            markersize=PAPER_MARKERSIZE,
            markevery=3,
            linewidth=PAPER_LINEWIDTH,
            label=label,
        )

    plt.xlabel("time step")
    plt.ylabel("number of patients")
    plt.xlim(0, 25)
    plt.ylim(0, 80)
    finish_plot(path, legend_loc="upper right", legend_fontsize=7.2)


def baseline_sensitivity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        row
        for row in rows
        if row["experiment"] == "paper_lambda_sweep"
        and int(row["N"]) == PAPER_FIG_6_8_N
        and abs(float(row["density"]) - PAPER_FIG_6_8_DENSITY) < 1e-9
    ]


def plot_sensitivity(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    relevant_rows = baseline_sensitivity_rows(rows)
    for model in ("strong", "hub"):
        model_rows = [row for row in relevant_rows if row["model"] == model]
        lambdas = np.array([float(row["lambda_ss"]) for row in model_rows])
        attack_rates = np.array([float(row["mean_attack_rate"]) for row in model_rows])
        order = np.argsort(lambdas)
        marker = "o" if model == "strong" else "s"
        linestyle = "-" if model == "strong" else "--"
        plt.plot(
            lambdas[order],
            attack_rates[order],
            color=MODEL_COLORS[model],
            linestyle=linestyle,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=MODEL_COLORS[model],
            markersize=PAPER_MARKERSIZE,
            linewidth=PAPER_LINEWIDTH,
            label="Strong" if model == "strong" else "Hub",
        )

    plt.xlabel(r"$\lambda$")
    plt.ylabel("Mean attack rate")
    plt.ylim(-0.05, 1.05)
    finish_plot(path, legend_loc="lower right")


def plot_velocity(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    relevant_rows = baseline_sensitivity_rows(rows)
    for model in ("strong", "hub"):
        model_rows = [row for row in relevant_rows if row["model"] == model]
        lambdas = np.array([float(row["lambda_ss"]) for row in model_rows])
        speeds = np.array([float(row["mean_propagation_speed"]) for row in model_rows])
        order = np.argsort(lambdas)
        marker = "o" if model == "strong" else "s"
        linestyle = "-" if model == "strong" else "--"
        plt.plot(
            lambdas[order],
            speeds[order],
            color=MODEL_COLORS[model],
            linestyle=linestyle,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=MODEL_COLORS[model],
            markersize=PAPER_MARKERSIZE,
            linewidth=PAPER_LINEWIDTH,
            label="Strong" if model == "strong" else "Hub",
        )

    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"velocity (/ $r_0 \cdot s$)")
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.6)
    finish_plot(path, legend_loc="upper left")


def plot_front_distance(rows: list[dict[str, object]], path: Path, model: str) -> None:
    plt.figure(figsize=PAPER_FIGSIZE)
    model_rows = [row for row in rows if row["model"] == model]
    for lambda_ss in sorted({float(row["lambda_ss"]) for row in model_rows}):
        lambda_rows = [row for row in model_rows if abs(float(row["lambda_ss"]) - lambda_ss) < 1e-9]
        times = np.array([int(row["time"]) for row in lambda_rows])
        front_distance = np.array([float(row["mean_front_distance"]) for row in lambda_rows])
        marker = LAMBDA_MARKERS.get(round(lambda_ss, 1), "o")
        color = LAMBDA_COLORS.get(round(lambda_ss, 1), "#111827")
        plt.plot(
            times,
            front_distance,
            color=color,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=color,
            markersize=PAPER_MARKERSIZE,
            markevery=4,
            linewidth=PAPER_LINEWIDTH,
            label=rf"$\lambda$={lambda_ss:.1f}",
        )

    plt.xlabel("time step")
    plt.ylabel(r"$r_f / r_0$")
    plt.xlim(0.0, 40.0)
    plt.ylim(0.0, 12.0)
    finish_plot(
        path,
        legend_loc="upper center",
        legend_ncol=3,
        legend_fontsize=7.0,
        legend_bbox=(0.5, -0.16),
    )


def plot_infection_route(
    results: list[tuple[ExperimentConfig, SimulationResult]],
    path: Path,
    model: str,
) -> None:
    selection = select_route_result(results, model)
    if selection is None:
        return
    config, result = selection

    positions = result.positions
    infected_ids = np.flatnonzero(np.isfinite(result.infected_time))
    superspreader_ids = np.flatnonzero(result.is_superspreader)
    infected_mask = np.isfinite(result.infected_time)
    normal_mask = ~result.is_superspreader
    susceptible_normal = normal_mask & ~infected_mask
    susceptible_superspreader = result.is_superspreader & ~infected_mask
    infected_normal = normal_mask & infected_mask
    infected_superspreader = result.is_superspreader & infected_mask

    plt.figure(figsize=PAPER_ROUTE_FIGSIZE)
    for target_id, source_id in enumerate(result.infection_source):
        if source_id < 0:
            continue
        source_x, source_y = positions[source_id]
        target_x, target_y = positions[target_id]
        delta_x = target_x - source_x
        if abs(delta_x) <= config.L / 2.0:
            plt.plot(
                [source_x, target_x],
                [source_y, target_y],
                color=ROUTE_COLORS["edge"],
                alpha=0.55,
                linewidth=0.45,
            )
        else:
            if delta_x > 0:
                target_x_unwrapped = target_x - config.L
                boundary_x = 0.0
                wrapped_boundary_x = config.L
            else:
                target_x_unwrapped = target_x + config.L
                boundary_x = config.L
                wrapped_boundary_x = 0.0
            fraction = (boundary_x - source_x) / (target_x_unwrapped - source_x)
            boundary_y = source_y + fraction * (target_y - source_y)
            plt.plot(
                [source_x, boundary_x],
                [source_y, boundary_y],
                color=ROUTE_COLORS["edge"],
                alpha=0.55,
                linewidth=0.45,
            )
            plt.plot(
                [wrapped_boundary_x, target_x],
                [boundary_y, target_y],
                color=ROUTE_COLORS["edge"],
                alpha=0.55,
                linewidth=0.45,
            )

    plt.scatter(
        positions[susceptible_normal, 0],
        positions[susceptible_normal, 1],
        s=7,
        facecolors="white",
        edgecolors=ROUTE_COLORS["susceptible_normal"],
        linewidths=0.45,
        label="S (normal)",
        zorder=2,
    )
    if superspreader_ids.size:
        plt.scatter(
            positions[susceptible_superspreader, 0],
            positions[susceptible_superspreader, 1],
            s=24,
            facecolors="white",
            edgecolors=ROUTE_COLORS["susceptible_superspreader"],
            linewidths=0.8,
            label="S (superspreader)",
            zorder=3,
        )
    plt.scatter(
        positions[infected_normal, 0],
        positions[infected_normal, 1],
        s=11,
        color=ROUTE_COLORS["infected_normal"],
        linewidths=0,
        label="I (normal)",
        zorder=3,
    )
    if superspreader_ids.size:
        plt.scatter(
            positions[infected_superspreader, 0],
            positions[infected_superspreader, 1],
            s=28,
            facecolors=ROUTE_COLORS["infected_superspreader"],
            edgecolors=ROUTE_COLORS["infected_superspreader"],
            linewidths=0.8,
            label="I (superspreader)",
            zorder=4,
        )
    plt.scatter(positions[0, 0], positions[0, 1], s=44, marker="*", color=ROUTE_COLORS["initial"], zorder=5)
    plt.title("route of infection", fontsize=10)
    plt.xlabel("")
    plt.ylabel("")
    plt.xlim(0, config.L)
    plt.ylim(0, config.L)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.legend(
        frameon=True,
        framealpha=0.88,
        facecolor="white",
        edgecolor="0.86",
        fontsize=6.8,
        handlelength=1.0,
        loc="lower left",
        borderpad=0.35,
        labelspacing=0.25,
    )
    finish_plot(path, legend=False)


def select_route_result(
    results: list[tuple[ExperimentConfig, SimulationResult]],
    model: str,
) -> tuple[ExperimentConfig, SimulationResult] | None:
    """Select the deterministic representative used for one infection-route figure."""
    candidates = [(config, result) for config, result in results if config.model == model]
    if not candidates:
        return None
    if model == "normal":
        non_percolating = [
            (config, result)
            for config, result in candidates
            if not has_percolated_to_top(result, config.L)
        ]
        return max(non_percolating or candidates, key=lambda item: item[1].total_infected)
    return max(candidates, key=lambda item: item[1].total_infected)


def select_route_results(
    results: list[tuple[ExperimentConfig, SimulationResult]],
) -> list[tuple[ExperimentConfig, SimulationResult]]:
    selected = [select_route_result(results, model) for model in ("normal", "strong", "hub")]
    return [item for item in selected if item is not None]


def route_selection_rows(
    results: list[tuple[ExperimentConfig, SimulationResult]],
) -> list[dict[str, object]]:
    return [
        {
            "experiment": config.experiment,
            "model": config.model,
            "lambda_ss": f"{config.lambda_ss:.4f}",
            "N": config.N,
            "density": f"{config.density:.6f}",
            "L": f"{config.L:.6f}",
            "seed": config.seed,
            "max_steps": config.max_steps,
            "total_infected": result.total_infected,
            "percolated": int(has_percolated_to_top(result, config.L)),
        }
        for config, result in results
    ]


def route_configs_from_rows(rows: list[dict[str, object]]) -> list[ExperimentConfig]:
    configs = [
        ExperimentConfig(
            experiment=str(row["experiment"]),
            model=str(row["model"]),
            lambda_ss=float(row["lambda_ss"]),
            N=int(row["N"]),
            density=float(row["density"]),
            seed=int(row["seed"]),
            max_steps=int(row["max_steps"]),
        )
        for row in rows
    ]
    if {config.model for config in configs} != {"normal", "strong", "hub"}:
        raise ValueError("route_plot_selections.csv must contain normal, strong, and hub rows")
    return configs


def write_plots(
    plots_dir: Path,
    *,
    infection_rows: list[dict[str, object]],
    percolation_rows: list[dict[str, object]],
    critical_rows: list[dict[str, object]],
    critical_reference_rows: list[dict[str, object]],
    front_rows: list[dict[str, object]],
    sensitivity_rows: list[dict[str, object]],
    curves_rows: list[dict[str, object]],
    secondary_rows: list[dict[str, object]],
    sars_secondary_rows: list[dict[str, object]],
    sars_epicurve_rows: list[dict[str, object]],
    sars_model_curve_rows: list[dict[str, object]],
    route_results: list[tuple[ExperimentConfig, SimulationResult]],
) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    for path in plots_dir.glob("*.png"):
        path.unlink()

    plot_infection_probability(infection_rows, plots_dir / PLOT_FILENAMES["infection_probability_strong"], "strong")
    plot_infection_probability(infection_rows, plots_dir / PLOT_FILENAMES["infection_probability_hub"], "hub")
    plot_percolation(percolation_rows, plots_dir / PLOT_FILENAMES["percolation_probability"])
    plot_percolation_model(percolation_rows, plots_dir / PLOT_FILENAMES["percolation_probability_strong"], "strong")
    plot_percolation_model(percolation_rows, plots_dir / PLOT_FILENAMES["percolation_probability_hub"], "hub")
    plot_critical_density(
        critical_rows,
        critical_reference_rows,
        plots_dir / PLOT_FILENAMES["critical_density"],
    )
    plot_front_distance(front_rows, plots_dir / PLOT_FILENAMES["front_distance_strong"], "strong")
    plot_front_distance(front_rows, plots_dir / PLOT_FILENAMES["front_distance_hub"], "hub")
    plot_velocity(sensitivity_rows, plots_dir / PLOT_FILENAMES["velocity_vs_lambda"])
    plot_epidemic_curves(curves_rows, plots_dir / PLOT_FILENAMES["epidemic_curves"])
    plot_secondary_distribution(secondary_rows, plots_dir / PLOT_FILENAMES["secondary_distribution"])
    plot_secondary_distribution_normal(secondary_rows, plots_dir / PLOT_FILENAMES["secondary_distribution_normal"])
    plot_secondary_distribution_superspreaders(
        secondary_rows,
        plots_dir / PLOT_FILENAMES["secondary_distribution_superspreaders"],
    )
    plot_sars_secondary_distribution(
        sars_secondary_rows,
        secondary_rows,
        plots_dir / PLOT_FILENAMES["sars_secondary_patients"],
    )
    plot_sars_epidemic_comparison(
        sars_epicurve_rows,
        sars_model_curve_rows,
        plots_dir / PLOT_FILENAMES["sars_epidemic_curve_comparison"],
    )
    plot_sensitivity(sensitivity_rows, plots_dir / PLOT_FILENAMES["sensitivity_lambda_attack_rate"])
    plot_infection_route(route_results, plots_dir / PLOT_FILENAMES["infection_route_normal"], "normal")
    plot_infection_route(route_results, plots_dir / PLOT_FILENAMES["infection_route_strong"], "strong")
    plot_infection_route(route_results, plots_dir / PLOT_FILENAMES["infection_route_hub"], "hub")


def regenerate_plots_from_csv(output_dir: Path) -> None:
    route_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    selection_rows = read_csv(output_dir / "route_plot_selections.csv")
    for config in route_configs_from_rows(selection_rows):
        result, _metrics = run_simulation(config)
        route_results.append((config, result))

    write_plots(
        output_dir / "plots",
        infection_rows=read_csv(output_dir / "infection_probability_functions.csv"),
        percolation_rows=read_csv(output_dir / "percolation_probability.csv"),
        critical_rows=read_csv(output_dir / "critical_density.csv"),
        critical_reference_rows=read_csv(output_dir / "critical_density_reference_curves.csv"),
        front_rows=read_csv(output_dir / "front_distance.csv"),
        sensitivity_rows=read_csv(output_dir / "sensitivity_summary.csv"),
        curves_rows=read_csv(output_dir / "epidemic_curves.csv"),
        secondary_rows=read_csv(output_dir / "secondary_distribution.csv"),
        sars_secondary_rows=read_csv(output_dir / "sars_singapore_secondary_patients.csv"),
        sars_epicurve_rows=read_csv(output_dir / "sars_singapore_epidemic_curve.csv"),
        sars_model_curve_rows=read_csv(output_dir / "sars_epidemic_model_curves.csv"),
        route_results=route_results,
    )


def main() -> None:
    args = parse_args()
    if args.seeds <= 0:
        raise ValueError("--seeds must be positive")
    if args.seed_offset < 0:
        raise ValueError("--seed-offset must be non-negative")
    if args.max_steps <= 0:
        raise ValueError("--max-steps must be positive")
    if args.jobs <= 0:
        raise ValueError("--jobs must be positive")
    if args.chunksize <= 0:
        raise ValueError("--chunksize must be positive")

    seeds = tuple(range(args.seed_offset, args.seed_offset + args.seeds))
    output_dir = args.output_dir
    if args.plot_only:
        regenerate_plots_from_csv(output_dir)
        print(f"regenerated paper-style plots in {output_dir / 'plots'}")
        return

    plots_dir = prepare_output_dir(output_dir)

    baseline_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    paper_curve_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    secondary_route_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    sars_comparison_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    front_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    route_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    all_metrics: list[RunMetrics] = []

    all_configs = (
        baseline_configs(seeds, args.max_steps)
        + percolation_configs(seeds, args.max_steps)
        + paper_curve_configs(seeds, args.max_steps)
        + secondary_route_configs(seeds, args.max_steps)
        + sars_comparison_configs(seeds, args.max_steps)
        + sensitivity_configs(seeds, args.max_steps)
    )

    def collect_result(index: int, config: ExperimentConfig, result: SimulationResult, metrics: RunMetrics) -> None:
        all_metrics.append(metrics)
        if config.experiment == "baseline_curves":
            baseline_results.append((config, result))
        if config.experiment == "paper_fig_8_curves":
            paper_curve_results.append((config, result))
        if config.experiment == "paper_fig_9_13_networks":
            secondary_route_results.append((config, result))
            route_results.append((config, result))
        if config.experiment == "sars_epidemic_comparison":
            sars_comparison_results.append((config, result))
        if (
            config.experiment == "paper_lambda_sweep"
            and config.N == PAPER_FIG_6_8_N
            and abs(config.density - PAPER_FIG_6_8_DENSITY) < 1e-9
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
            completed_runs = executor.map(run_config, all_configs, chunksize=args.chunksize)
            for index, (config, result, metrics) in enumerate(completed_runs, start=1):
                collect_result(index, config, result, metrics)

    summary_rows = metric_rows(all_metrics)
    write_csv(output_dir / "summary_metrics.csv", summary_rows)

    baseline_metrics = [metric for metric in all_metrics if metric.experiment == "baseline_curves"]
    baseline_rows = group_metrics(baseline_metrics, ("model", "lambda_ss", "N", "density"))
    write_csv(output_dir / "baseline_summary.csv", baseline_rows)

    percolation_metrics = [metric for metric in all_metrics if metric.experiment == "percolation_density_sweep"]
    percolation_rows = group_metrics(percolation_metrics, ("model", "lambda_ss", "N", "density"))
    write_csv(output_dir / "percolation_probability.csv", percolation_rows)

    critical_rows = critical_density_rows(percolation_rows)
    write_csv(output_dir / "critical_density.csv", critical_rows)
    critical_reference_rows = critical_density_reference_rows()
    write_csv(output_dir / "critical_density_reference_curves.csv", critical_reference_rows)

    propagation_rows = group_metrics(all_metrics, ("experiment", "model", "lambda_ss", "N", "density"))
    write_csv(output_dir / "propagation_speed.csv", propagation_rows)

    front_rows = front_distance_rows(front_results, args.max_steps)
    write_csv(output_dir / "front_distance.csv", front_rows)

    curves_rows = epidemic_curve_rows(paper_curve_results, args.max_steps)
    write_csv(output_dir / "epidemic_curves.csv", curves_rows)

    sars_model_curve_rows = epidemic_curve_rows(sars_comparison_results, args.max_steps)
    write_csv(output_dir / "sars_epidemic_model_curves.csv", sars_model_curve_rows)

    secondary_rows = secondary_distribution_rows(secondary_route_results)
    write_csv(output_dir / "secondary_distribution.csv", secondary_rows)

    infection_rows = infection_probability_rows()
    write_csv(output_dir / "infection_probability_functions.csv", infection_rows)

    sars_secondary_rows = sars_secondary_patient_rows()
    write_csv(output_dir / "sars_singapore_secondary_patients.csv", sars_secondary_rows)

    sars_epicurve_rows = sars_epidemic_curve_rows()
    write_csv(output_dir / "sars_singapore_epidemic_curve.csv", sars_epicurve_rows)

    sensitivity_metrics = [metric for metric in all_metrics if metric.experiment == "paper_lambda_sweep"]
    sensitivity_rows = group_metrics(sensitivity_metrics, ("experiment", "model", "lambda_ss", "N", "density"))
    write_csv(output_dir / "sensitivity_summary.csv", sensitivity_rows)

    selected_route_results = select_route_results(route_results)
    write_csv(output_dir / "route_plot_selections.csv", route_selection_rows(selected_route_results))

    plot_infection_probability(infection_rows, plots_dir / PLOT_FILENAMES["infection_probability_strong"], "strong")
    plot_infection_probability(infection_rows, plots_dir / PLOT_FILENAMES["infection_probability_hub"], "hub")
    plot_percolation(percolation_rows, plots_dir / PLOT_FILENAMES["percolation_probability"])
    plot_percolation_model(percolation_rows, plots_dir / PLOT_FILENAMES["percolation_probability_strong"], "strong")
    plot_percolation_model(percolation_rows, plots_dir / PLOT_FILENAMES["percolation_probability_hub"], "hub")
    plot_critical_density(
        critical_rows,
        critical_reference_rows,
        plots_dir / PLOT_FILENAMES["critical_density"],
    )
    plot_front_distance(front_rows, plots_dir / PLOT_FILENAMES["front_distance_strong"], "strong")
    plot_front_distance(front_rows, plots_dir / PLOT_FILENAMES["front_distance_hub"], "hub")
    plot_velocity(sensitivity_rows, plots_dir / PLOT_FILENAMES["velocity_vs_lambda"])
    plot_epidemic_curves(curves_rows, plots_dir / PLOT_FILENAMES["epidemic_curves"])
    plot_secondary_distribution(secondary_rows, plots_dir / PLOT_FILENAMES["secondary_distribution"])
    plot_secondary_distribution_normal(secondary_rows, plots_dir / PLOT_FILENAMES["secondary_distribution_normal"])
    plot_secondary_distribution_superspreaders(
        secondary_rows,
        plots_dir / PLOT_FILENAMES["secondary_distribution_superspreaders"],
    )
    plot_sars_secondary_distribution(
        sars_secondary_rows,
        secondary_rows,
        plots_dir / PLOT_FILENAMES["sars_secondary_patients"],
    )
    plot_sars_epidemic_comparison(
        sars_epicurve_rows,
        sars_model_curve_rows,
        plots_dir / PLOT_FILENAMES["sars_epidemic_curve_comparison"],
    )
    plot_sensitivity(sensitivity_rows, plots_dir / PLOT_FILENAMES["sensitivity_lambda_attack_rate"])
    plot_infection_route(selected_route_results, plots_dir / PLOT_FILENAMES["infection_route_normal"], "normal")
    plot_infection_route(selected_route_results, plots_dir / PLOT_FILENAMES["infection_route_strong"], "strong")
    plot_infection_route(selected_route_results, plots_dir / PLOT_FILENAMES["infection_route_hub"], "hub")

    print(f"wrote clean results to {output_dir}")


if __name__ == "__main__":
    main()
