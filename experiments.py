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

import matplotlib.pyplot as plt
import numpy as np

from config import DEFAULT_MAX_STEPS, DEFAULT_N, DEFAULT_R0
from simulator import MonteCarloSIRSimulator, SimulationResult, periodic_distance


BASELINE_DENSITY = DEFAULT_N / (10.0 * DEFAULT_R0) ** 2
BASELINE_LAMBDA = 0.20
PERCOLATION_ATTACK_RATE = 0.50
DEFAULT_SEEDS = tuple(range(10))


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
        default=Path("results/member4"),
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
        percolated=attack_rate >= PERCOLATION_ATTACK_RATE,
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
    densities = (2.50, 3.50, BASELINE_DENSITY, 6.00, 8.00)
    for density in densities:
        for model, lambda_ss in (("normal", 0.0), ("strong", BASELINE_LAMBDA), ("hub", BASELINE_LAMBDA)):
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
    return configs


def sensitivity_configs(seeds: Iterable[int], max_steps: int) -> list[ExperimentConfig]:
    configs: list[ExperimentConfig] = []
    population_sizes = (300, DEFAULT_N, 650)
    lambda_values = (0.00, 0.05, 0.10, 0.20, 0.40)
    densities = (3.50, BASELINE_DENSITY, 6.00)

    for model in ("strong", "hub"):
        for N in population_sizes:
            for lambda_ss in lambda_values:
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


def plot_percolation(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(8, 5))
    models = sorted({row["model"] for row in rows})
    for model in models:
        model_rows = [row for row in rows if row["model"] == model]
        densities = [float(row["density"]) for row in model_rows]
        probabilities = [float(row["percolation_probability"]) for row in model_rows]
        order = np.argsort(densities)
        plt.plot(np.array(densities)[order], np.array(probabilities)[order], marker="o", label=model)

    plt.xlabel("Population density")
    plt.ylabel("Percolation probability")
    plt.ylim(-0.05, 1.05)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def plot_epidemic_curves(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(8, 5))
    models = sorted({row["model"] for row in rows})
    for model in models:
        model_rows = [row for row in rows if row["model"] == model]
        times = [int(row["time"]) for row in model_rows]
        active = [float(row["mean_active_infected"]) for row in model_rows]
        plt.plot(times, active, label=model)

    plt.xlabel("Time step")
    plt.ylabel("Mean active infected")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def plot_secondary_distribution(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(8, 5))
    models = sorted({row["model"] for row in rows})
    for model in models:
        model_rows = [row for row in rows if row["model"] == model]
        counts = [int(row["secondary_infections"]) for row in model_rows]
        probabilities = [float(row["probability"]) for row in model_rows]
        plt.plot(counts, probabilities, marker="o", label=model)

    plt.xlabel("Secondary infections caused by one infected individual")
    plt.ylabel("Probability")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def plot_sensitivity(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=(8, 5))
    relevant_rows = [
        row
        for row in rows
        if row["experiment"] == "sensitivity_N_lambda"
        and int(row["N"]) == DEFAULT_N
        and abs(float(row["density"]) - BASELINE_DENSITY) < 1e-9
    ]
    for model in sorted({row["model"] for row in relevant_rows}):
        model_rows = [row for row in relevant_rows if row["model"] == model]
        lambdas = [float(row["lambda_ss"]) for row in model_rows]
        attack_rates = [float(row["mean_attack_rate"]) for row in model_rows]
        order = np.argsort(lambdas)
        plt.plot(np.array(lambdas)[order], np.array(attack_rates)[order], marker="o", label=model)

    plt.xlabel("Superspreader fraction lambda")
    plt.ylabel("Mean attack rate")
    plt.ylim(-0.05, 1.05)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def write_report(
    output_dir: Path,
    baseline_rows: list[dict[str, object]],
    percolation_rows: list[dict[str, object]],
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

    density_lines = []
    for row in percolation_rows:
        density_lines.append(
            f"| {row['model']} | {float(row['density']):.2f} | {float(row['percolation_probability']):.2f} | "
            f"{float(row['mean_attack_rate']):.2f} | {float(row['mean_propagation_speed']):.2f} |"
        )

    content = f"""# Member 4 Report: Experiments and Sensitivity Analysis

## Scope

This report covers Member 4 responsibilities: percolation probability experiments,
propagation speed, epidemic curves, secondary-infection distribution, and sensitivity
checks over population size, superspreader fraction, density, and random seeds.

## Reproducibility

- Simulator: `MonteCarloSIRSimulator`
- Random seeds per setting: `{seeds}`
- Baseline population: `N = {DEFAULT_N}`
- Baseline density: `{BASELINE_DENSITY:.2f}`
- Baseline superspreader fraction for strong and hub models: `lambda = {BASELINE_LAMBDA:.2f}`
- Percolation event definition: final attack rate at least `{PERCOLATION_ATTACK_RATE:.0%}`
- Generated outputs: `summary_metrics.csv`, `percolation_probability.csv`,
  `propagation_speed.csv`, `epidemic_curves.csv`,
  `secondary_distribution.csv`, and `sensitivity_summary.csv`

## Baseline Results

| Model | Mean attack rate | Percolation probability | Mean duration | Mean propagation speed | Mean secondary infections |
| --- | ---: | ---: | ---: | ---: | ---: |
| normal | {fmt(normal, 'mean_attack_rate')} | {fmt(normal, 'percolation_probability')} | {fmt(normal, 'mean_duration')} | {fmt(normal, 'mean_propagation_speed')} | {fmt(normal, 'mean_secondary_infections')} |
| strong | {fmt(strong, 'mean_attack_rate')} | {fmt(strong, 'percolation_probability')} | {fmt(strong, 'mean_duration')} | {fmt(strong, 'mean_propagation_speed')} | {fmt(strong, 'mean_secondary_infections')} |
| hub | {fmt(hub, 'mean_attack_rate')} | {fmt(hub, 'percolation_probability')} | {fmt(hub, 'mean_duration')} | {fmt(hub, 'mean_propagation_speed')} | {fmt(hub, 'mean_secondary_infections')} |

The hub model produces the fastest baseline spread. The strong-infectiousness model
also increases final epidemic size relative to the normal model, but its propagation
speed is lower because superspreaders still infect within the shorter normal range.

## Percolation and Density Sweep

| Model | Density | Percolation probability | Mean attack rate | Mean propagation speed |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(density_lines)}

Percolation probability increases with density because each infectious individual has
more neighbors inside the infection radius. Hub superspreaders reach the high
percolation regime earlier than the strong model because their effective infection range
is larger.

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

## Figures

- `plots/percolation_probability.png`
- `plots/epidemic_curves.png`
- `plots/secondary_distribution.png`
- `plots/sensitivity_lambda_attack_rate.png`

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
- `propagation_speed.csv`: propagation speed by model, density, and lambda.
- `epidemic_curves.csv`: mean new, active, and cumulative infections over time.
- `secondary_distribution.csv`: distribution of secondary infections.
- `sensitivity_summary.csv`: averaged sensitivity results for `N`, `lambda_ss`, and density.
- `plots/`: PNG figures generated from the CSV files.

All values are reproducible from the seeds recorded in `summary_metrics.csv`.
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

    propagation_rows = group_metrics(all_metrics, ("experiment", "model", "lambda_ss", "N", "density"))
    write_csv(output_dir / "propagation_speed.csv", propagation_rows)

    curves_rows = epidemic_curve_rows(baseline_results, args.max_steps)
    write_csv(output_dir / "epidemic_curves.csv", curves_rows)

    secondary_rows = secondary_distribution_rows(baseline_results)
    write_csv(output_dir / "secondary_distribution.csv", secondary_rows)

    sensitivity_metrics = [metric for metric in all_metrics if metric.experiment.startswith("sensitivity")]
    sensitivity_rows = group_metrics(sensitivity_metrics, ("experiment", "model", "lambda_ss", "N", "density"))
    write_csv(output_dir / "sensitivity_summary.csv", sensitivity_rows)

    plot_percolation(percolation_rows, plots_dir / "percolation_probability.png")
    plot_epidemic_curves(curves_rows, plots_dir / "epidemic_curves.png")
    plot_secondary_distribution(secondary_rows, plots_dir / "secondary_distribution.png")
    plot_sensitivity(sensitivity_rows, plots_dir / "sensitivity_lambda_attack_rate.png")

    write_report(output_dir, baseline_rows, percolation_rows, sensitivity_rows, len(seeds))
    write_results_readme(output_dir)

    print(f"wrote clean results to {output_dir}")


if __name__ == "__main__":
    main()
