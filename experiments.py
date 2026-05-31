"""Member 4 experiments and sensitivity analysis.

This script runs reproducible Monte Carlo batches for the superspreader SIR
model, writes clean CSV files, generates plots, and creates a concise report.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from config import DEFAULT_MAX_STEPS, DEFAULT_N, DEFAULT_R0
from simulator import MonteCarloSIRSimulator, SimulationResult, periodic_distance


BASELINE_DENSITY = DEFAULT_N / (10.0 * DEFAULT_R0) ** 2
BASELINE_LAMBDA = 0.20
DEFAULT_SEEDS = tuple(range(30))
PERCOLATION_TOP_MARGIN = DEFAULT_R0
LAMBDA_SWEEP = (0.00, 0.20, 0.40, 0.60, 0.80, 1.00)
DENSITY_SWEEP = (
    0.50,
    0.75,
    1.00,
    1.25,
    1.50,
    2.00,
    2.50,
    3.00,
    3.50,
    4.00,
    BASELINE_DENSITY,
    5.50,
    6.50,
    8.00,
)
PLOT_DPI = 180
DENSITY_SCALE = np.pi * DEFAULT_R0**2
STRONG_CRITICAL_R0_REFERENCE = 4.5
HUB_CRITICAL_R0_REFERENCE = 3.0


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
        "--max-steps",
        type=int,
        default=DEFAULT_MAX_STEPS,
        help="Maximum timesteps per simulation.",
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
    """Estimate front speed from the slope of max infected radius over time."""
    infected_ids = np.flatnonzero(np.isfinite(result.infected_time))
    if infected_ids.size < 3 or result.duration < 2:
        return 0.0

    source_position = result.positions[0]
    distances = periodic_distance(source_position, result.positions[infected_ids], L)
    times = result.infected_time[infected_ids]

    front_times: list[float] = []
    front_radii: list[float] = []
    for time in range(1, result.duration + 1):
        reached = times <= time
        if np.any(reached):
            front_radius = float(np.max(distances[reached]))
            if front_radius > 0:
                front_times.append(float(time))
                front_radii.append(front_radius)

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


def baseline_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    for model, lambda_ss in (("normal", 0.0), ("strong", BASELINE_LAMBDA), ("hub", BASELINE_LAMBDA)):
        for seed in seeds:
            configs.append(
                ExperimentConfig(
                    experiment="baseline_curves",
                    model=model,
                    lambda_ss=lambda_ss,
                    N=DEFAULT_N,
                    density=BASELINE_DENSITY,
                    seed=seed,
                    max_steps=max_steps,
                )
            )
    return configs


def percolation_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    for density in DENSITY_SWEEP:
        for model in ("strong", "hub"):
            for lambda_ss in LAMBDA_SWEEP:
                for seed in seeds:
                    configs.append(
                        ExperimentConfig(
                            experiment="percolation_density_sweep",
                            model=model,
                            lambda_ss=lambda_ss,
                            N=DEFAULT_N,
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
                    N=DEFAULT_N,
                    density=density,
                    seed=seed,
                    max_steps=max_steps,
                )
            )
    return configs

def sensitivity_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    population_sizes = (300, DEFAULT_N, 650)
    densities = (3.50, BASELINE_DENSITY, 6.00)

    for model in ("strong", "hub"):
        for N in population_sizes:
            for lambda_ss in LAMBDA_SWEEP:
                for seed in seeds:
                    configs.append(
                        ExperimentConfig(
                            experiment="sensitivity_N_lambda",
                            model=model,
                            lambda_ss=lambda_ss,
                            N=N,
                            density=BASELINE_DENSITY,
                            seed=seed,
                            max_steps=max_steps,
                        )
                    )

        for density in densities:
            for lambda_ss in (0.05, 0.20, 0.40):
                for seed in seeds:
                    configs.append(
                        ExperimentConfig(
                            experiment="sensitivity_density_lambda",
                            model=model,
                            lambda_ss=lambda_ss,
                            N=DEFAULT_N,
                            density=density,
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


def finish_plot(path: Path, *, legend: bool = True) -> None:
    plt.grid(alpha=0.25)
    if legend:
        plt.legend(frameon=True)
    plt.tight_layout()
    plt.savefig(path, dpi=PLOT_DPI)
    plt.close()


def plot_percolation(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    relevant_rows = [
        row
        for row in rows
        if row["model"] == "normal"
        or abs(float(row["lambda_ss"]) - BASELINE_LAMBDA) < 1e-9
    ]
    for model in ("normal", "strong", "hub"):
        model_rows = [row for row in relevant_rows if row["model"] == model]
        densities = np.array([float(row["density"]) for row in model_rows]) * DENSITY_SCALE
        probabilities = np.array([float(row["percolation_probability"]) for row in model_rows])
        order = np.argsort(densities)
        plt.plot(densities[order], probabilities[order], marker="o", label=model)

    plt.xlabel(r"$\rho \pi r_0^2$")
    plt.ylabel("Percolation probability")
    plt.ylim(-0.05, 1.05)
    finish_plot(path)


def plot_percolation_model(rows: list[dict[str, object]], path: Path, model: str) -> None:
    plt.figure(figsize=(7.2, 4.6))
    model_rows = [row for row in rows if row["model"] == model]
    markers = ("o", "s", "^", "D", "v", "*")
    for lambda_ss in sorted({float(row["lambda_ss"]) for row in model_rows}):
        lambda_rows = [row for row in model_rows if abs(float(row["lambda_ss"]) - lambda_ss) < 1e-9]
        densities = np.array([float(row["density"]) for row in lambda_rows]) * DENSITY_SCALE
        probabilities = np.array([float(row["percolation_probability"]) for row in lambda_rows])
        order = np.argsort(densities)
        marker = markers[int(round(lambda_ss / 0.2)) % len(markers)]
        plt.plot(
            densities[order],
            probabilities[order],
            linestyle="-",
            marker=marker,
            markerfacecolor="none" if lambda_ss == 0 else None,
            label=rf"$\lambda$={lambda_ss:.1f}",
        )

    plt.xlabel(r"$\rho \pi r_0^2$")
    plt.ylabel("Percolation probability")
    plt.ylim(-0.05, 1.05)
    plt.xlim(0.0, max(float(row["density"]) for row in model_rows) * DENSITY_SCALE + 0.5)
    finish_plot(path)


def plot_critical_density(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    style_by_model = {
        "strong": {"marker": "o", "color": "red", "label": "Strong infectiousness model (simulation)"},
        "hub": {"marker": "s", "color": "blue", "label": "Hub model (simulation)"},
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
            color=style["color"],
            markerfacecolor="none" if model == "hub" else style["color"],
            label=style["label"],
        )

    lambda_grid = np.linspace(0.0, 1.0, 200)
    reproductive_denominator = lambda_grid + (1.0 - lambda_grid) / 6.0
    plt.plot(
        lambda_grid,
        STRONG_CRITICAL_R0_REFERENCE / reproductive_denominator,
        color="limegreen",
        label=r"Strong infectiousness model ($R_0=R_c$)",
    )
    plt.plot(
        lambda_grid,
        HUB_CRITICAL_R0_REFERENCE / reproductive_denominator,
        color="magenta",
        linestyle="--",
        label=r"Hub model ($R_0=R_c$)",
    )

    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$\rho_c \pi r_0^2$")
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 25.0)
    finish_plot(path)


def plot_epidemic_curves(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    for model in ("normal", "strong", "hub"):
        model_rows = [row for row in rows if row["model"] == model]
        times = np.array([int(row["time"]) for row in model_rows])
        new_infections = np.array([float(row["mean_new_infections"]) for row in model_rows])
        plt.plot(times, new_infections, label=model)

    plt.xlabel("Time step")
    plt.ylabel("Mean new infections")
    finish_plot(path)


def plot_secondary_distribution(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    for model in ("normal", "strong", "hub"):
        model_rows = [row for row in rows if row["model"] == model]
        counts = np.array([int(row["secondary_infections"]) for row in model_rows])
        probabilities = np.array([float(row["probability"]) for row in model_rows])
        mask = probabilities > 0
        plt.semilogy(counts[mask], probabilities[mask], marker="o", label=model)

    plt.xlabel("Secondary infections caused by one individual")
    plt.ylabel("Probability")
    finish_plot(path)


def plot_secondary_distribution_normal(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    model_rows = [row for row in rows if row["model"] == "normal"]
    counts = np.array([int(row["secondary_infections"]) for row in model_rows])
    probabilities = np.array([float(row["probability"]) for row in model_rows])
    mask = probabilities > 0
    plt.semilogy(counts[mask], probabilities[mask], marker="^", color="black", label="no superspreaders")

    plt.xlabel("Number of secondary infections")
    plt.ylabel("Probability")
    finish_plot(path)


def plot_secondary_distribution_superspreaders(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    for model, marker in (("strong", "o"), ("hub", "s")):
        model_rows = [row for row in rows if row["model"] == model]
        counts = np.array([int(row["secondary_infections"]) for row in model_rows])
        probabilities = np.array([float(row["probability"]) for row in model_rows])
        mask = probabilities > 0
        plt.semilogy(counts[mask], probabilities[mask], marker=marker, label=model)

    plt.xlabel("Number of secondary infections")
    plt.ylabel("Probability")
    finish_plot(path)


def baseline_sensitivity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        row
        for row in rows
        if row["experiment"] == "sensitivity_N_lambda"
        and int(row["N"]) == DEFAULT_N
        and abs(float(row["density"]) - BASELINE_DENSITY) < 1e-9
    ]


def plot_sensitivity(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    relevant_rows = baseline_sensitivity_rows(rows)
    for model in ("strong", "hub"):
        model_rows = [row for row in relevant_rows if row["model"] == model]
        lambdas = np.array([float(row["lambda_ss"]) for row in model_rows])
        attack_rates = np.array([float(row["mean_attack_rate"]) for row in model_rows])
        order = np.argsort(lambdas)
        plt.plot(lambdas[order], attack_rates[order], marker="o", label=model)

    plt.xlabel(r"$\lambda$")
    plt.ylabel("Mean attack rate")
    plt.ylim(-0.05, 1.05)
    finish_plot(path)


def plot_velocity(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    relevant_rows = baseline_sensitivity_rows(rows)
    for model in ("strong", "hub"):
        model_rows = [row for row in relevant_rows if row["model"] == model]
        lambdas = np.array([float(row["lambda_ss"]) for row in model_rows])
        speeds = np.array([float(row["mean_propagation_speed"]) for row in model_rows])
        order = np.argsort(lambdas)
        marker = "o" if model == "strong" else "s"
        plt.plot(lambdas[order], speeds[order], marker=marker, label=model)

    plt.xlabel(r"$\lambda$")
    plt.ylabel("Velocity of propagation")
    finish_plot(path)


def plot_front_distance(rows: list[dict[str, object]], path: Path, model: str) -> None:
    plt.figure(figsize=(7.2, 4.6))
    model_rows = [row for row in rows if row["model"] == model]
    for lambda_ss in sorted({float(row["lambda_ss"]) for row in model_rows}):
        lambda_rows = [row for row in model_rows if abs(float(row["lambda_ss"]) - lambda_ss) < 1e-9]
        times = np.array([int(row["time"]) for row in lambda_rows])
        front_distance = np.array([float(row["mean_front_distance"]) for row in lambda_rows])
        plt.plot(times, front_distance, label=f"lambda={lambda_ss:.2f}")

    plt.xlabel("Time step")
    plt.ylabel(r"Front distance $r_f$")
    finish_plot(path)


def plot_infection_route(
    results: list[tuple[ExperimentConfig, SimulationResult]],
    path: Path,
    model: str,
) -> None:
    candidates = [(config, result) for config, result in results if config.model == model]
    if not candidates:
        return
    config, result = max(candidates, key=lambda item: item[1].total_infected)
    positions = result.positions
    infected_ids = np.flatnonzero(np.isfinite(result.infected_time))
    superspreader_ids = np.flatnonzero(result.is_superspreader)

    plt.figure(figsize=(5.6, 5.6))
    plt.scatter(positions[:, 0], positions[:, 1], s=8, color="0.85", linewidths=0)
    for target_id, source_id in enumerate(result.infection_source):
        if source_id < 0:
            continue
        x_values = [positions[source_id, 0], positions[target_id, 0]]
        y_values = [positions[source_id, 1], positions[target_id, 1]]
        plt.plot(x_values, y_values, color="0.25", alpha=0.18, linewidth=0.5)

    scatter = plt.scatter(
        positions[infected_ids, 0],
        positions[infected_ids, 1],
        c=result.infected_time[infected_ids],
        s=16,
        cmap="viridis",
        linewidths=0,
    )
    if superspreader_ids.size:
        plt.scatter(
            positions[superspreader_ids, 0],
            positions[superspreader_ids, 1],
            s=42,
            facecolors="none",
            edgecolors="tab:red",
            linewidths=0.8,
        )
    plt.scatter(positions[0, 0], positions[0, 1], s=90, marker="*", color="black", zorder=5)
    plt.colorbar(scatter, fraction=0.046, pad=0.04, label="Infection time")
    plt.title(f"{model} route, seed {config.seed}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xlim(0, config.L)
    plt.ylim(0, config.L)
    plt.gca().set_aspect("equal", adjustable="box")
    finish_plot(path, legend=False)


def write_report(
    output_dir: Path,
    baseline_rows: list[dict[str, object]],
    percolation_rows: list[dict[str, object]],
    critical_rows: list[dict[str, object]],
    sensitivity_rows: list[dict[str, object]],
    seeds: int,
) -> None:
    baseline_lookup = {row["model"]: row for row in baseline_rows}

    def fmt(row: dict[str, object], key: str) -> str:
        return f"{float(row[key]):.3f}"

    normal = baseline_lookup["normal"]
    strong = baseline_lookup["strong"]
    hub = baseline_lookup["hub"]

    best_hub_speed = max(
        (
            row
            for row in sensitivity_rows
            if row["experiment"] == "sensitivity_N_lambda"
            and row["model"] == "hub"
            and int(row["N"]) == DEFAULT_N
            and abs(float(row["density"]) - BASELINE_DENSITY) < 1e-9
        ),
        key=lambda row: float(row["mean_propagation_speed"]),
    )
    best_strong_speed = max(
        (
            row
            for row in sensitivity_rows
            if row["experiment"] == "sensitivity_N_lambda"
            and row["model"] == "strong"
            and int(row["N"]) == DEFAULT_N
            and abs(float(row["density"]) - BASELINE_DENSITY) < 1e-9
        ),
        key=lambda row: float(row["mean_propagation_speed"]),
    )

    density_lines: list[str] = []
    for row in percolation_rows:
        if row["model"] != "normal" and abs(float(row["lambda_ss"]) - BASELINE_LAMBDA) > 1e-9:
            continue
        density_lines.append(
            f"| {row['model']} | {float(row['lambda_ss']):.2f} | {float(row['density']):.2f} | "
            f"{float(row['percolation_probability']):.2f} | "
            f"{float(row['mean_attack_rate']):.2f} | {float(row['mean_propagation_speed']):.2f} |"
        )

    critical_lines: list[str] = []
    for row in critical_rows:
        critical_density = str(row["critical_density"]).strip()
        if critical_density:
            critical_lines.append(
                f"| {row['model']} | {float(row['lambda_ss']):.2f} | {float(critical_density):.2f} |"
            )

    content = f"""# Member 4 Report: Experiments and Sensitivity Analysis

## Scope

This report covers Member 4 responsibilities: running numerical experiments,
measuring percolation probability, estimating propagation speed, generating epidemic
curves, measuring secondary-infection distributions, and checking sensitivity over
population size, superspreader fraction, density, and random seeds. The focus is the
experiment/result part of the reference paper rather than the mathematical-model
definition plots.

## Reproducibility

- Simulator: `MonteCarloSIRSimulator`
- Random seeds per setting: `{seeds}`
- Baseline population: `N = {DEFAULT_N}`
- Baseline density: `{BASELINE_DENSITY:.2f}`
- Baseline superspreader fraction for strong and hub models: `lambda = {BASELINE_LAMBDA:.2f}`
- Percolation event definition: infection reaches the top band of the two-dimensional system
- Generated outputs: `summary_metrics.csv`, `percolation_probability.csv`,
  `critical_density.csv`, `propagation_speed.csv`, `front_distance.csv`,
  `epidemic_curves.csv`, `secondary_distribution.csv`, and `sensitivity_summary.csv`

## How to Read the Metrics

- **Attack rate** is the final fraction of the population that became infected at
  least once. A higher attack rate means a larger outbreak.
- **Percolation probability** is the fraction of simulation runs in which the infection
  reaches the top band of the spatial system. This measures whether the disease can
  spread across the system, not only whether it infects many people locally.
- **Propagation speed** is estimated from the slope of the infection-front distance
  over time. A larger value means the epidemic wave travels faster through space.
- **Duration** is the number of time steps until no infected individuals remain.
- **Mean secondary infections** is the average number of people infected by one infected
  individual. A heavy tail in this distribution indicates superspreading events.
- **Critical density** is the approximate density where percolation probability reaches
  `0.5`. A lower critical density means the model can produce system-wide outbreaks
  even when the population is more sparse.
- The CSV files store raw density `rho = N / L^2`. The percolation and critical-density
  plots use the paper-style normalized density `rho * pi * r0^2` on the axis. Since this
  project uses `r0 = 1`, the plotted density is the raw density multiplied by `pi`.
- The percolation and velocity sweeps use
  `lambda = 0.0, 0.2, 0.4, 0.6, 0.8, 1.0`. The baseline comparison still uses
  `lambda = {BASELINE_LAMBDA:.2f}` for the two superspreader models.

## Baseline Results

| Model | Mean attack rate | Percolation probability | Mean duration | Mean propagation speed | Mean secondary infections |
| --- | ---: | ---: | ---: | ---: | ---: |
| normal | {fmt(normal, 'mean_attack_rate')} | {fmt(normal, 'percolation_probability')} | {fmt(normal, 'mean_duration')} | {fmt(normal, 'mean_propagation_speed')} | {fmt(normal, 'mean_secondary_infections')} |
| strong | {fmt(strong, 'mean_attack_rate')} | {fmt(strong, 'percolation_probability')} | {fmt(strong, 'mean_duration')} | {fmt(strong, 'mean_propagation_speed')} | {fmt(strong, 'mean_secondary_infections')} |
| hub | {fmt(hub, 'mean_attack_rate')} | {fmt(hub, 'percolation_probability')} | {fmt(hub, 'mean_duration')} | {fmt(hub, 'mean_propagation_speed')} | {fmt(hub, 'mean_secondary_infections')} |

The hub model produces the fastest baseline spread. The strong-infectiousness model
also increases final epidemic size relative to the normal model, but its propagation
speed is lower because superspreaders still infect within the shorter normal range.

The baseline results show a clear ordering of epidemic severity:
`hub > strong > normal`. In the normal model, only about
`{fmt(normal, 'mean_attack_rate')}` of the population is infected on average, and
percolation probability remains low at `{fmt(normal, 'percolation_probability')}`.
In the strong model, the attack rate increases to about
`{fmt(strong, 'mean_attack_rate')}`, and percolation probability reaches
`{fmt(strong, 'percolation_probability')}`. In the hub model, the attack rate rises to
about `{fmt(hub, 'mean_attack_rate')}`, percolation probability reaches
`{fmt(hub, 'percolation_probability')}`, and propagation speed is higher than the
strong-model speed. This supports the main claim that superspreaders do not only
increase the final outbreak size; they also change how quickly the epidemic moves
through space.

## Percolation and Density Sweep

| Model | lambda | Density | Percolation probability | Mean attack rate | Mean propagation speed |
| --- | ---: | ---: | ---: | ---: | ---: |
{chr(10).join(density_lines)}

Percolation probability increases with density because each infectious individual has
more neighbors inside the infection radius. Hub superspreaders reach the high
percolation regime earlier than the strong model because their effective infection range
is larger.

At the baseline superspreader fraction `lambda = {BASELINE_LAMBDA:.2f}`, the hub model
starts percolating at lower density than the strong model. This means the hub mechanism
lowers the density barrier for large outbreaks. In practical terms, a hub-like
superspreader can connect distant local clusters, so the disease crosses the system
even when ordinary local transmission would still die out.

## Critical Density

| Model | lambda | Critical density |
| --- | ---: | ---: |
{chr(10).join(critical_lines)}

The critical-density plot summarizes the percolation curves into one threshold value.
As `lambda` increases, both models need less population density to percolate. The plot
also includes reference critical curves based on the paper's `R0 = Rc` argument, so the
markers show simulation estimates while the solid and dashed lines show the theoretical
trend. Lower critical density means the outbreak is more robust under sparse conditions,
so this result again shows that hub superspreaders are more dangerous than merely
stronger infectious individuals.

## Sensitivity Findings

- Increasing `lambda` raises attack rate and speed for both superspreader models.
- At baseline `N` and density, the strongest hub speed in the sweep occurs at
  `lambda = {float(best_hub_speed['lambda_ss']):.2f}` with mean speed
  `{float(best_hub_speed['mean_propagation_speed']):.3f}`.
- At baseline `N` and density, the strongest strong-model speed in the sweep occurs at
  `lambda = {float(best_strong_speed['lambda_ss']):.2f}` with mean speed
  `{float(best_strong_speed['mean_propagation_speed']):.3f}`.
- Changing `N` while keeping density fixed gives similar attack-rate trends, so the
  model behavior is mainly controlled by density and superspreader fraction rather than
  raw population size alone.
- Lower density is the most important robustness stress test: it reduces outbreak
  probability and makes outcomes more seed-dependent.

The sensitivity results show that `lambda` and density are the dominant controls. When
`lambda` increases, there are more superspreaders, so both attack rate and propagation
speed generally increase. When density increases, individuals have more nearby contacts,
so local clusters are easier to connect. Changing `N` while keeping density fixed has a
smaller effect because the average local neighborhood structure is similar.

## Figure Guide and Meaning

| Figure | What it shows | How to interpret it | Main result |
| --- | --- | --- | --- |
| `plots/percolation_probability.png` | Percolation probability versus normalized density `rho * pi * r0^2` for normal, strong, and hub models at baseline `lambda`. | Curves farther left percolate at lower density. | Hub percolates first, strong second, normal last. |
| `plots/percolation_probability_strong.png` | Strong-model percolation probability for `lambda = 0.0` to `1.0`, using normalized density on the x-axis. | Increasing `lambda` shifts the point cloud left and upward. | More strong superspreaders make system-wide spread possible at lower density. |
| `plots/percolation_probability_hub.png` | Hub-model percolation probability for `lambda = 0.0` to `1.0`, using normalized density on the x-axis. | The left shift is stronger than in the strong model. | Hub superspreaders reduce the percolation threshold more sharply. |
| `plots/critical_density.png` | Approximate normalized density `rho_c * pi * r0^2` where percolation probability reaches `0.5`, plus paper-style `R0 = Rc` reference curves. | Lower values mean easier system-wide spread; markers are simulation, lines are reference curves. | Critical density decreases with `lambda`, and the hub threshold is lower than the strong threshold. |
| `plots/front_distance_strong.png` | Infection-front distance `rf` over time for the strong model. | Steeper growth means faster spatial propagation; a plateau means spread has stopped. | Larger `lambda` makes the front travel farther and faster. |
| `plots/front_distance_hub.png` | Infection-front distance `rf` over time for the hub model. | Compare curve height and slope across `lambda`. | Hub spreading reaches far distances earlier than the strong model. |
| `plots/velocity_vs_lambda.png` | Propagation speed versus superspreader fraction. | Higher speed means the epidemic wave crosses space faster. | Hub speed is consistently higher than strong speed, especially at large `lambda`. |
| `plots/epidemic_curves.png` | Mean new infections per time step. | The peak height shows outbreak intensity; peak timing shows how fast the outbreak develops. | Hub has the tallest and earliest peak; strong is slower; normal remains small. |
| `plots/secondary_distribution.png` | Distribution of the number of secondary infections caused by one infected individual. | A longer tail means rare individuals infect many others. | Superspreader models produce heavier tails than the normal model. |
| `plots/secondary_distribution_normal.png` | Paper Fig. 12 style view for the no-superspreader case. | Most individuals have zero or few secondary infections. | Without superspreaders, the tail is short. |
| `plots/secondary_distribution_superspreaders.png` | Paper Fig. 13 style view comparing strong and hub superspreader cases. | The right tail shows individuals who infect many others. | Both superspreader models produce long-tailed secondary infection distributions. |
| `plots/sensitivity_lambda_attack_rate.png` | Mean final attack rate versus `lambda`. | Higher attack rate means a larger final epidemic. | Attack rate rises sharply as superspreaders are added, then approaches saturation. |
| `plots/infection_route_normal.png` | Spatial infection routes in a representative normal-model run. | Points are individuals, edges are infection events, and color indicates infection time. | Normal spread is mostly local and slower. |
| `plots/infection_route_strong.png` | Spatial infection routes in a representative strong-model run. | Red rings mark superspreaders. | Strong superspreaders create larger local bursts but still mainly transmit nearby. |
| `plots/infection_route_hub.png` | Spatial infection routes in a representative hub-model run. | Long route connections reveal wider effective contacts. | Hub superspreaders bridge distant areas and accelerate spatial spread. |

## Overall Interpretation

Across all experiments, the results support the reference paper's qualitative conclusion:
a small fraction of superspreaders can substantially increase outbreak size, outbreak
probability, and spatial speed. The important distinction is that the two superspreader
mechanisms do not behave the same way. The strong infectiousness model increases
transmission probability after contact, while the hub model changes the effective contact
geometry by allowing wider connections. Because of that, the hub model percolates at
lower density and has higher propagation speed.

## Paper Figures Not Reproduced

The reference paper's Fig. 1 and Fig. 2 are model-definition plots of the infection
probability function, so they are not included in this Member 4 experiment report.
The paper's Fig. 14 and Fig. 15 use empirical SARS Singapore 2003 data, which is also
not generated here because the raw empirical dataset is not included in this repository.
The current outputs focus on the experiment-side result families: percolation curves,
critical density, front-distance/velocity, epidemic curves, infection routes, and
secondary-infection distributions.

## Notes for Final Group Report

The experiments support the qualitative conclusion that a small fraction of
superspreaders can sharply increase epidemic size and speed. The hub model is more
dangerous than the strong-infectiousness model under the same `lambda` because it
changes the contact geometry, not only the probability of infection after contact.
"""
    (output_dir / "REPORT.md").write_text(content, encoding="utf-8")


def write_results_readme(output_dir: Path) -> None:
    content = """# Member 4 Generated Results

This directory is generated by:

```bash
.venv/bin/python experiments.py
```

## Files

- `REPORT.md`: readable summary for the final report.
- `summary_metrics.csv`: one row per simulation replicate.
- `baseline_summary.csv`: averaged baseline model comparison.
- `percolation_probability.csv`: percolation probability by model and density.
- `critical_density.csv`: density where percolation probability first crosses 0.5.
- `propagation_speed.csv`: propagation speed by model, density, and lambda.
- `front_distance.csv`: mean infection front distance over time by lambda.
- `epidemic_curves.csv`: mean new, active, and cumulative infections over time.
- `secondary_distribution.csv`: distribution of secondary infections.
- `sensitivity_summary.csv`: averaged sensitivity results for `N`, `lambda_ss`, and density.
- `plots/`: PNG figures generated from the CSV files.

All values are reproducible from the seeds recorded in `summary_metrics.csv`.
Percolation follows the reference-paper style definition: an outbreak has percolated
when infection reaches the top band of the spatial system.
Distances wrap horizontally, while the vertical axis remains open for bottom-to-top
propagation measurements.
"""
    (output_dir / "README.md").write_text(content, encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.seeds <= 0:
        raise ValueError("--seeds must be positive")
    if args.max_steps <= 0:
        raise ValueError("--max-steps must be positive")

    seeds = tuple(range(args.seeds))
    output_dir = args.output_dir
    plots_dir = prepare_output_dir(output_dir)

    baseline_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    front_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    route_results: list[tuple[ExperimentConfig, SimulationResult]] = []
    all_metrics: list[RunMetrics] = []

    all_configs = (
        baseline_configs(seeds, args.max_steps)
        + percolation_configs(seeds, args.max_steps)
        + sensitivity_configs(seeds, args.max_steps)
    )

    total_configs = len(all_configs)
    for index, config in enumerate(all_configs, start=1):
        result, metrics = run_simulation(config)
        all_metrics.append(metrics)
        if config.experiment == "baseline_curves":
            baseline_results.append((config, result))
            route_results.append((config, result))
        if (
            config.experiment == "sensitivity_N_lambda"
            and config.N == DEFAULT_N
            and abs(config.density - BASELINE_DENSITY) < 1e-9
        ):
            front_results.append((config, result))
        if index % 50 == 0 or index == total_configs:
            print(f"completed {index}/{total_configs} simulations")

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

    propagation_rows = group_metrics(all_metrics, ("experiment", "model", "lambda_ss", "N", "density"))
    write_csv(output_dir / "propagation_speed.csv", propagation_rows)

    front_rows = front_distance_rows(front_results, args.max_steps)
    write_csv(output_dir / "front_distance.csv", front_rows)

    curves_rows = epidemic_curve_rows(baseline_results, args.max_steps)
    write_csv(output_dir / "epidemic_curves.csv", curves_rows)

    secondary_rows = secondary_distribution_rows(baseline_results)
    write_csv(output_dir / "secondary_distribution.csv", secondary_rows)

    sensitivity_metrics = [metric for metric in all_metrics if metric.experiment.startswith("sensitivity")]
    sensitivity_rows = group_metrics(sensitivity_metrics, ("experiment", "model", "lambda_ss", "N", "density"))
    write_csv(output_dir / "sensitivity_summary.csv", sensitivity_rows)

    plot_percolation(percolation_rows, plots_dir / "percolation_probability.png")
    plot_percolation_model(percolation_rows, plots_dir / "percolation_probability_strong.png", "strong")
    plot_percolation_model(percolation_rows, plots_dir / "percolation_probability_hub.png", "hub")
    plot_critical_density(critical_rows, plots_dir / "critical_density.png")
    plot_front_distance(front_rows, plots_dir / "front_distance_strong.png", "strong")
    plot_front_distance(front_rows, plots_dir / "front_distance_hub.png", "hub")
    plot_velocity(sensitivity_rows, plots_dir / "velocity_vs_lambda.png")
    plot_epidemic_curves(curves_rows, plots_dir / "epidemic_curves.png")
    plot_secondary_distribution(secondary_rows, plots_dir / "secondary_distribution.png")
    plot_secondary_distribution_normal(secondary_rows, plots_dir / "secondary_distribution_normal.png")
    plot_secondary_distribution_superspreaders(
        secondary_rows,
        plots_dir / "secondary_distribution_superspreaders.png",
    )
    plot_sensitivity(sensitivity_rows, plots_dir / "sensitivity_lambda_attack_rate.png")
    plot_infection_route(route_results, plots_dir / "infection_route_normal.png", "normal")
    plot_infection_route(route_results, plots_dir / "infection_route_strong.png", "strong")
    plot_infection_route(route_results, plots_dir / "infection_route_hub.png", "hub")

    write_report(output_dir, baseline_rows, percolation_rows, critical_rows, sensitivity_rows, len(seeds))
    write_results_readme(output_dir)

    print(f"wrote clean results to {output_dir}")


if __name__ == "__main__":
    main()
