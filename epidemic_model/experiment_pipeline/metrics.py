"""Simulation execution and derived outbreak metrics."""

import numpy as np

from epidemic_model.simulator import MonteCarloSIRSimulator, SimulationResult, periodic_distance

from .records import SimulationConfig, SimulationMetrics
from .settings import TOP_REACH_MARGIN


def active_case_curve(result: SimulationResult, max_time_steps: int) -> list[int]:
    """Return active infected counts for t = 0..max_time_steps."""
    curve: list[int] = []
    recovered_time = result.recovered_time
    infected_time = result.infected_time

    for time in range(max_time_steps + 1):
        infected = np.isfinite(infected_time) & (infected_time <= time)
        not_recovered = np.isnan(recovered_time) | (recovered_time > time)
        curve.append(int(np.count_nonzero(infected & not_recovered)))
    return curve


def cumulative_case_curve(result: SimulationResult, max_time_steps: int) -> list[int]:
    """Return cumulative infected counts for t = 0..max_time_steps."""
    curve: list[int] = []
    infected_time = result.infected_time
    for time in range(max_time_steps + 1):
        curve.append(int(np.count_nonzero(np.isfinite(infected_time) & (infected_time <= time))))
    return curve


def new_case_curve(result: SimulationResult, max_time_steps: int) -> list[int]:
    """Return new infections for t = 0..max_time_steps."""
    curve = [1]
    curve.extend(int(value) for value in result.new_infections_per_step)
    if len(curve) < max_time_steps + 1:
        curve.extend([0] * (max_time_steps + 1 - len(curve)))
    return curve[: max_time_steps + 1]


def infection_front_curve(result: SimulationResult, domain_size: float, max_time_steps: int) -> list[float]:
    """Return the furthest infected distance from the initial case over time."""
    curve: list[float] = []
    source_position = result.positions[0]
    infected_time = result.infected_time
    for time in range(max_time_steps + 1):
        infected_ids = np.flatnonzero(np.isfinite(infected_time) & (infected_time <= time))
        if infected_ids.size == 0:
            curve.append(0.0)
            continue
        distances = periodic_distance(source_position, result.positions[infected_ids], domain_size)
        curve.append(float(np.max(distances)))
    return curve


def top_band_reached(result: SimulationResult, domain_size: float) -> bool:
    """Return whether infection reaches the top boundary band."""
    infected_ids = np.flatnonzero(np.isfinite(result.infected_time))
    if infected_ids.size == 0:
        return False
    top_threshold = max(0.0, domain_size - TOP_REACH_MARGIN)
    return bool(np.any(result.positions[infected_ids, 1] >= top_threshold))


def estimate_front_speed(result: SimulationResult, domain_size: float) -> float:
    """Estimate front speed from advancing front points before the plateau."""
    infected_ids = np.flatnonzero(np.isfinite(result.infected_time))
    if infected_ids.size < 3 or result.duration < 2:
        return 0.0

    source_position = result.positions[0]
    distances = periodic_distance(source_position, result.positions[infected_ids], domain_size)
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


def run_simulation(config: SimulationConfig) -> tuple[SimulationResult, SimulationMetrics]:
    simulator = MonteCarloSIRSimulator(
        N=config.population_size,
        model=config.transmission_model,
        lambda_ss=config.superspreader_fraction,
        L=config.domain_size,
        max_steps=config.max_time_steps,
        seed=config.random_seed,
    )
    result = simulator.run()
    active_curve = active_case_curve(result, config.max_time_steps)
    secondary_counts = result.secondary_counts[np.isfinite(result.infected_time)]
    case_fraction = result.total_infected / config.population_size

    metrics = SimulationMetrics(
        scenario_name=config.scenario_name,
        transmission_model=config.transmission_model,
        superspreader_fraction=config.superspreader_fraction,
        population_size=config.population_size,
        population_density=config.population_density,
        domain_size=config.domain_size,
        random_seed=config.random_seed,
        total_cases=result.total_infected,
        case_fraction=case_fraction,
        outbreak_duration=result.duration,
        peak_new_cases=max(result.new_infections_per_step, default=0),
        peak_active_cases=max(active_curve, default=0),
        front_speed=estimate_front_speed(result, config.domain_size),
        mean_secondary_cases=float(np.mean(secondary_counts)) if secondary_counts.size else 0.0,
        max_secondary_cases=int(np.max(secondary_counts)) if secondary_counts.size else 0,
        reached_top_boundary=top_band_reached(result, config.domain_size),
    )
    return result, metrics


def run_configured_simulation(config: SimulationConfig) -> tuple[SimulationConfig, SimulationResult, SimulationMetrics]:
    result, metrics = run_simulation(config)
    return config, result, metrics
