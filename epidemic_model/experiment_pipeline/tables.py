"""CSV row builders for simulation summaries and reference data."""

from datetime import date, timedelta
from statistics import mean
from typing import Iterable

import numpy as np

from epidemic_model.config import DEFAULT_R0, DEFAULT_W0
from epidemic_model.simulator import SimulationResult, infection_probability

from .metrics import active_case_curve, cumulative_case_curve, infection_front_curve, new_case_curve
from .records import SimulationConfig, SimulationMetrics
from .settings import (
    HUB_CRITICAL_REPRODUCTION_NUMBER,
    NORMALIZED_DENSITY_FACTOR,
    SARS_PROBABLE_CASE_COUNT,
    SARS_SECONDARY_CASE_FREQUENCIES,
    SARS_SIX_DAY_CASE_COUNTS,
    SARS_TIMESTEP_DAYS,
    STRONG_CRITICAL_REPRODUCTION_NUMBER,
)


def simulation_summary_rows(metrics: Iterable[SimulationMetrics]) -> list[dict[str, object]]:
    return [
        {
            "scenario_name": metric.scenario_name,
            "transmission_model": metric.transmission_model,
            "superspreader_fraction": f"{metric.superspreader_fraction:.4f}",
            "population_size": metric.population_size,
            "population_density": f"{metric.population_density:.4f}",
            "domain_size": f"{metric.domain_size:.4f}",
            "random_seed": metric.random_seed,
            "total_cases": metric.total_cases,
            "case_fraction": f"{metric.case_fraction:.6f}",
            "outbreak_duration": metric.outbreak_duration,
            "peak_new_cases": metric.peak_new_cases,
            "peak_active_cases": metric.peak_active_cases,
            "front_speed": f"{metric.front_speed:.6f}",
            "mean_secondary_cases": f"{metric.mean_secondary_cases:.6f}",
            "max_secondary_cases": metric.max_secondary_cases,
            "reached_top_boundary": int(metric.reached_top_boundary),
        }
        for metric in metrics
    ]


def infection_kernel_rows(points: int = 401) -> list[dict[str, object]]:
    max_radius = np.sqrt(6.0) * DEFAULT_R0
    distances = np.linspace(0.0, max_radius, points)
    rows: list[dict[str, object]] = []
    for transmission_model in ("strong", "hub"):
        for source_label, source_is_superspreader in (("normal", False), ("superspreader", True)):
            probabilities = infection_probability(
                distances,
                source_is_superspreader=source_is_superspreader,
                model=transmission_model,
                r0=DEFAULT_R0,
                w0=DEFAULT_W0,
            )
            for distance, probability in zip(distances, probabilities):
                rows.append(
                    {
                        "transmission_model": transmission_model,
                        "source_type": source_label,
                        "distance_over_r0": f"{distance / DEFAULT_R0:.6f}",
                        "probability_over_w0": f"{probability / DEFAULT_W0:.6f}",
                    }
                )
    return rows


def sars_secondary_case_rows() -> list[dict[str, object]]:
    total_frequency = sum(frequency for _count, frequency, _source in SARS_SECONDARY_CASE_FREQUENCIES)
    if total_frequency != SARS_PROBABLE_CASE_COUNT:
        raise ValueError(
            "SARS secondary-patient frequencies must account for all "
            f"{SARS_PROBABLE_CASE_COUNT} probable cases; got {total_frequency}"
        )

    rows: list[dict[str, object]] = []
    for secondary_count, frequency, source in SARS_SECONDARY_CASE_FREQUENCIES:
        rows.append(
            {
                "secondary_cases": secondary_count,
                "frequency": frequency,
                "probability": f"{frequency / SARS_PROBABLE_CASE_COUNT:.12f}",
                "denominator": SARS_PROBABLE_CASE_COUNT,
                "source": source,
            }
        )
    return rows


def sars_case_time_series_rows() -> list[dict[str, object]]:
    start = date(2003, 2, 13)
    total_cases = sum(SARS_SIX_DAY_CASE_COUNTS)
    rows: list[dict[str, object]] = []
    for time_step, cases in enumerate(SARS_SIX_DAY_CASE_COUNTS):
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


def grouped_metric_rows(metrics: Iterable[SimulationMetrics], keys: tuple[str, ...]) -> list[dict[str, object]]:
    buckets: dict[tuple[object, ...], list[SimulationMetrics]] = {}
    for metric in metrics:
        key = tuple(getattr(metric, name) for name in keys)
        buckets.setdefault(key, []).append(metric)

    rows: list[dict[str, object]] = []
    for key, values in sorted(buckets.items()):
        row = {name: value for name, value in zip(keys, key)}
        row.update(
            {
                "runs": len(values),
                "mean_case_fraction": f"{mean(item.case_fraction for item in values):.6f}",
                "std_case_fraction": f"{np.std([item.case_fraction for item in values], ddof=0):.6f}",
                "percolation_probability": f"{mean(float(item.reached_top_boundary) for item in values):.6f}",
                "mean_outbreak_duration": f"{mean(item.outbreak_duration for item in values):.6f}",
                "mean_peak_active_cases": f"{mean(item.peak_active_cases for item in values):.6f}",
                "mean_front_speed": f"{mean(item.front_speed for item in values):.6f}",
                "mean_secondary_cases": f"{mean(item.mean_secondary_cases for item in values):.6f}",
                "mean_max_secondary_cases": f"{mean(item.max_secondary_cases for item in values):.6f}",
            }
        )
        rows.append(row)
    return rows


def infection_curve_rows(
    results: list[tuple[SimulationConfig, SimulationResult]],
    max_time_steps: int,
) -> list[dict[str, object]]:
    buckets: dict[tuple[str, float], list[SimulationResult]] = {}
    for config, result in results:
        buckets.setdefault((config.transmission_model, config.superspreader_fraction), []).append(result)

    rows: list[dict[str, object]] = []
    for (transmission_model, superspreader_fraction), result_group in sorted(buckets.items()):
        new_curves = np.array([new_case_curve(result, max_time_steps) for result in result_group])
        active_curves = np.array([active_case_curve(result, max_time_steps) for result in result_group])
        cumulative_curves = np.array([cumulative_case_curve(result, max_time_steps) for result in result_group])

        for time in range(max_time_steps + 1):
            rows.append(
                {
                    "transmission_model": transmission_model,
                    "superspreader_fraction": f"{superspreader_fraction:.4f}",
                    "time": time,
                    "mean_new_cases": f"{np.mean(new_curves[:, time]):.6f}",
                    "mean_active_cases": f"{np.mean(active_curves[:, time]):.6f}",
                    "mean_cumulative_cases": f"{np.mean(cumulative_curves[:, time]):.6f}",
                }
            )
    return rows


def secondary_case_distribution_rows(
    results: list[tuple[SimulationConfig, SimulationResult]],
) -> list[dict[str, object]]:
    buckets: dict[tuple[str, float], list[int]] = {}
    for config, result in results:
        infected_counts = result.secondary_counts[np.isfinite(result.infected_time)]
        buckets.setdefault((config.transmission_model, config.superspreader_fraction), []).extend(int(value) for value in infected_counts)

    rows: list[dict[str, object]] = []
    for (transmission_model, superspreader_fraction), counts in sorted(buckets.items()):
        total = len(counts)
        for secondary_count in range(max(counts) + 1 if counts else 1):
            frequency = counts.count(secondary_count)
            rows.append(
                {
                    "transmission_model": transmission_model,
                    "superspreader_fraction": f"{superspreader_fraction:.4f}",
                    "secondary_cases": secondary_count,
                    "frequency": frequency,
                    "probability": f"{frequency / total:.6f}" if total else "0.000000",
                }
            )
    return rows


def infection_front_rows(
    results: list[tuple[SimulationConfig, SimulationResult]],
    max_time_steps: int,
) -> list[dict[str, object]]:
    buckets: dict[tuple[str, float], list[list[float]]] = {}
    for config, result in results:
        curve = infection_front_curve(result, config.domain_size, max_time_steps)
        buckets.setdefault((config.transmission_model, config.superspreader_fraction), []).append(curve)

    rows: list[dict[str, object]] = []
    for (transmission_model, superspreader_fraction), curves in sorted(buckets.items()):
        curve_array = np.array(curves, dtype=float)
        for time in range(max_time_steps + 1):
            rows.append(
                {
                    "transmission_model": transmission_model,
                    "superspreader_fraction": f"{superspreader_fraction:.4f}",
                    "time": time,
                    "mean_front_distance": f"{np.mean(curve_array[:, time]):.6f}",
                }
            )
    return rows


def critical_density_estimate_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    critical_rows: list[dict[str, object]] = []
    for transmission_model in ("strong", "hub"):
        model_rows = [row for row in rows if row["transmission_model"] == transmission_model]
        for superspreader_fraction in sorted({float(row["superspreader_fraction"]) for row in model_rows}):
            fraction_rows = [row for row in model_rows if abs(float(row["superspreader_fraction"]) - superspreader_fraction) < 1e-9]
            ordered = sorted(fraction_rows, key=lambda row: float(row["population_density"]))
            critical_density: float | None = None

            previous_density = float(ordered[0]["population_density"])
            previous_probability = float(ordered[0]["percolation_probability"])
            for row in ordered:
                population_density = float(row["population_density"])
                probability = float(row["percolation_probability"])
                if probability >= 0.5:
                    if probability == previous_probability:
                        critical_density = population_density
                    else:
                        fraction = (0.5 - previous_probability) / (probability - previous_probability)
                        critical_density = previous_density + fraction * (population_density - previous_density)
                    break
                previous_density = population_density
                previous_probability = probability

            critical_rows.append(
                {
                    "transmission_model": transmission_model,
                    "superspreader_fraction": f"{superspreader_fraction:.4f}",
                    "critical_density": (
                        f"{critical_density:.6f}" if critical_density is not None else ""
                    ),
                }
            )
    return critical_rows


def critical_density_theory_rows(points: int = 201) -> list[dict[str, object]]:
    """Return the analytical Eq. 3-5 reference curves used in Fig. 5."""
    if points < 2:
        raise ValueError("points must be at least 2")

    rows: list[dict[str, object]] = []
    for superspreader_fraction in np.linspace(0.0, 1.0, points):
        reproduction_factor = superspreader_fraction + (1.0 - superspreader_fraction) / 6.0
        for transmission_model, critical_r0 in (
            ("strong", STRONG_CRITICAL_REPRODUCTION_NUMBER),
            ("hub", HUB_CRITICAL_REPRODUCTION_NUMBER),
        ):
            rows.append(
                {
                    "transmission_model": transmission_model,
                    "superspreader_fraction": f"{superspreader_fraction:.6f}",
                    "reproduction_factor": f"{reproduction_factor:.12f}",
                    "critical_reproduction_number": f"{critical_r0:.6f}",
                    "normalized_critical_density": f"{critical_r0 / reproduction_factor:.6f}",
                }
            )
    return rows
