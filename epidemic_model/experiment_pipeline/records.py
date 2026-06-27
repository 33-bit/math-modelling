"""Structured records used by scenario generation and reporting."""

from dataclasses import dataclass

import numpy as np

from epidemic_model.config import DEFAULT_MAX_STEPS


@dataclass(frozen=True)
class SimulationConfig:
    """Parameters that define one simulation replicate."""

    scenario_name: str
    transmission_model: str
    superspreader_fraction: float
    population_size: int
    population_density: float
    random_seed: int
    max_time_steps: int = DEFAULT_MAX_STEPS

    @property
    def domain_size(self) -> float:
        return float(np.sqrt(self.population_size / self.population_density))


@dataclass(frozen=True)
class SimulationMetrics:
    """Compact metrics collected from one simulation replicate."""

    scenario_name: str
    transmission_model: str
    superspreader_fraction: float
    population_size: int
    population_density: float
    domain_size: float
    random_seed: int
    total_cases: int
    case_fraction: float
    outbreak_duration: int
    peak_new_cases: int
    peak_active_cases: int
    front_speed: float
    mean_secondary_cases: float
    max_secondary_cases: int
    reached_top_boundary: bool
