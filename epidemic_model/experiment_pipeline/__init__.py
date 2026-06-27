"""Purpose-named experiment pipeline API."""

from .cli import main
from .records import SimulationConfig, SimulationMetrics
from .settings import DENSITY_SWEEP_POPULATIONS
from .tables import (
    critical_density_estimate_rows,
    critical_density_theory_rows,
    sars_case_time_series_rows,
    sars_secondary_case_rows,
    simulation_summary_rows,
)
from .metrics import top_band_reached
from .plotting import route_configs_from_saved_rows

__all__ = [
    "DENSITY_SWEEP_POPULATIONS",
    "SimulationConfig",
    "SimulationMetrics",
    "critical_density_estimate_rows",
    "critical_density_theory_rows",
    "main",
    "route_configs_from_saved_rows",
    "sars_case_time_series_rows",
    "sars_secondary_case_rows",
    "simulation_summary_rows",
    "top_band_reached",
]
