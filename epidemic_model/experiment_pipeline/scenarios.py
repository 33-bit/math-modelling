"""Purpose-named simulation scenario builders."""

from typing import Iterable

from .records import SimulationConfig
from .settings import (
    DEFAULT_SUPERSPREADER_FRACTION,
    DENSITY_SWEEP_POPULATIONS,
    DOMAIN_SIDE_LENGTH,
    ROUTE_ANALYSIS_DENSITY,
    ROUTE_ANALYSIS_POPULATION,
    SARS_ANALYSIS_DENSITY,
    SARS_ANALYSIS_POPULATION,
    SARS_SUPERSPREADER_FRACTION,
    SPREAD_ANALYSIS_DENSITY,
    SPREAD_ANALYSIS_POPULATION,
    SUPERSPREADER_FRACTION_SWEEP,
)
from epidemic_model.config import DEFAULT_R0


def fixed_density_baseline_configs(seeds: Iterable[int], max_time_steps: int) -> list[SimulationConfig]:
    configs: list[SimulationConfig] = []
    for transmission_model, superspreader_fraction in (("normal", 0.0), ("strong", DEFAULT_SUPERSPREADER_FRACTION), ("hub", DEFAULT_SUPERSPREADER_FRACTION)):
        for random_seed in seeds:
            configs.append(
                SimulationConfig(
                    scenario_name="fixed_density_baseline",
                    transmission_model=transmission_model,
                    superspreader_fraction=superspreader_fraction,
                    population_size=SPREAD_ANALYSIS_POPULATION,
                    population_density=SPREAD_ANALYSIS_DENSITY,
                    random_seed=random_seed,
                    max_time_steps=max_time_steps,
                )
            )
    return configs


def density_sweep_configs(seeds: Iterable[int], max_time_steps: int) -> list[SimulationConfig]:
    configs: list[SimulationConfig] = []
    for population_size in DENSITY_SWEEP_POPULATIONS:
        population_density = population_size / (10.0 * DEFAULT_R0) ** 2
        for transmission_model in ("strong", "hub"):
            for superspreader_fraction in SUPERSPREADER_FRACTION_SWEEP:
                for random_seed in seeds:
                    configs.append(
                        SimulationConfig(
                            scenario_name="percolation_density_sweep",
                            transmission_model=transmission_model,
                            superspreader_fraction=superspreader_fraction,
                            population_size=population_size,
                            population_density=population_density,
                            random_seed=random_seed,
                            max_time_steps=max_time_steps,
                        )
                    )

        for random_seed in seeds:
            configs.append(
                SimulationConfig(
                    scenario_name="percolation_density_sweep",
                    transmission_model="normal",
                    superspreader_fraction=0.0,
                    population_size=population_size,
                    population_density=population_density,
                    random_seed=random_seed,
                    max_time_steps=max_time_steps,
                )
            )
    return configs


def epidemic_curve_configs(seeds: Iterable[int], max_time_steps: int) -> list[SimulationConfig]:
    configs: list[SimulationConfig] = []
    for transmission_model, superspreader_fraction in (("normal", 0.0), ("strong", DEFAULT_SUPERSPREADER_FRACTION), ("hub", DEFAULT_SUPERSPREADER_FRACTION)):
        for random_seed in seeds:
            configs.append(
                SimulationConfig(
                    scenario_name="epidemic_curve_comparison",
                    transmission_model=transmission_model,
                    superspreader_fraction=superspreader_fraction,
                    population_size=SPREAD_ANALYSIS_POPULATION,
                    population_density=SPREAD_ANALYSIS_DENSITY,
                    random_seed=random_seed,
                    max_time_steps=max_time_steps,
                )
            )
    return configs


def route_and_secondary_configs(seeds: Iterable[int], max_time_steps: int) -> list[SimulationConfig]:
    configs: list[SimulationConfig] = []
    for transmission_model, superspreader_fraction in (("normal", 0.0), ("strong", DEFAULT_SUPERSPREADER_FRACTION), ("hub", DEFAULT_SUPERSPREADER_FRACTION)):
        for random_seed in seeds:
            configs.append(
                SimulationConfig(
                    scenario_name="route_and_secondary_distribution",
                    transmission_model=transmission_model,
                    superspreader_fraction=superspreader_fraction,
                    population_size=ROUTE_ANALYSIS_POPULATION,
                    population_density=ROUTE_ANALYSIS_DENSITY,
                    random_seed=random_seed,
                    max_time_steps=max_time_steps,
                )
            )
    return configs


def sars_curve_comparison_configs(seeds: Iterable[int], max_time_steps: int) -> list[SimulationConfig]:
    configs: list[SimulationConfig] = []
    for transmission_model, superspreader_fraction in (("normal", 0.0), ("strong", SARS_SUPERSPREADER_FRACTION), ("hub", SARS_SUPERSPREADER_FRACTION)):
        for random_seed in seeds:
            configs.append(
                SimulationConfig(
                    scenario_name="sars_epidemic_comparison",
                    transmission_model=transmission_model,
                    superspreader_fraction=superspreader_fraction,
                    population_size=SARS_ANALYSIS_POPULATION,
                    population_density=SARS_ANALYSIS_DENSITY,
                    random_seed=random_seed,
                    max_time_steps=max_time_steps,
                )
            )
    return configs


def superspreader_fraction_sweep_configs(seeds: Iterable[int], max_time_steps: int) -> list[SimulationConfig]:
    configs: list[SimulationConfig] = []
    for transmission_model in ("strong", "hub"):
        for superspreader_fraction in SUPERSPREADER_FRACTION_SWEEP:
            for random_seed in seeds:
                configs.append(
                    SimulationConfig(
                        scenario_name="superspreader_fraction_sweep",
                        transmission_model=transmission_model,
                        superspreader_fraction=superspreader_fraction,
                        population_size=SPREAD_ANALYSIS_POPULATION,
                        population_density=SPREAD_ANALYSIS_DENSITY,
                        random_seed=random_seed,
                        max_time_steps=max_time_steps,
                    )
                )
    return configs
