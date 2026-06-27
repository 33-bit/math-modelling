"""Filesystem helpers for generated scenario_name artifacts."""

import csv
from pathlib import Path


COLUMN_ALIASES = {
    "scenario_name": ("experiment",),
    "transmission_model": ("model",),
    "superspreader_fraction": ("lambda_ss",),
    "population_size": ("N",),
    "population_density": ("density",),
    "domain_size": ("L",),
    "random_seed": ("seed",),
    "max_time_steps": ("max_steps",),
    "case_fraction": ("attack_rate",),
    "outbreak_duration": ("duration",),
    "peak_new_cases": ("peak_new_infections",),
    "peak_active_cases": ("peak_active_infected",),
    "front_speed": ("propagation_speed",),
    "mean_case_fraction": ("mean_attack_rate",),
    "std_case_fraction": ("std_attack_rate",),
    "mean_outbreak_duration": ("mean_duration",),
    "mean_peak_active_cases": ("mean_peak_active_infected",),
    "mean_front_speed": ("mean_propagation_speed",),
    "mean_new_cases": ("mean_new_infections",),
    "mean_active_cases": ("mean_active_infected",),
    "mean_cumulative_cases": ("mean_cumulative_infected",),
    "mean_secondary_cases": ("mean_secondary_infections",),
    "mean_max_secondary_cases": ("mean_max_secondary_infections",),
    "max_secondary_cases": ("max_secondary_infections",),
    "secondary_cases": ("secondary_infections",),
    "critical_reproduction_number": ("critical_R0",),
    "reached_top_boundary": ("percolated",),
}


def normalize_csv_row(row: dict[str, object]) -> dict[str, object]:
    normalized = dict(row)
    for canonical_name, legacy_names in COLUMN_ALIASES.items():
        if canonical_name in normalized:
            continue
        for legacy_name in legacy_names:
            if legacy_name in normalized:
                normalized[canonical_name] = normalized[legacy_name]
                break
    return normalized


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
        return [normalize_csv_row(row) for row in csv.DictReader(file)]


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
