"""Monte Carlo SIR simulator with Fujie-Odagaki superspreader variants."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from epidemic_model.config import (
    DEFAULT_GAMMA,
    DEFAULT_MAX_STEPS,
    DEFAULT_N,
    DEFAULT_R0,
    DEFAULT_W0,
    INFECTED,
    RECOVERED,
    SUSCEPTIBLE,
    VALID_MODELS,
    resolve_L,
)


def minimum_image_delta(
    source: np.ndarray,
    targets: np.ndarray,
    L: float,
    *,
    periodic_y: bool = False,
) -> np.ndarray:
    """Return signed displacement with cylindrical or fully periodic boundaries.

    The experiments use the cylindrical case because the paper's bottom-to-top
    percolation and front-speed definitions require an open vertical direction.
    Set ``periodic_y=True`` only for analyses that do not use those metrics.
    """
    source_array = np.asarray(source, dtype=float)
    targets_array = np.asarray(targets, dtype=float)
    delta = targets_array - source_array
    half_L = L / 2.0
    wrapped = delta.copy()
    wrapped[..., 0] = np.where(
        wrapped[..., 0] > half_L,
        wrapped[..., 0] - L,
        np.where(wrapped[..., 0] < -half_L, wrapped[..., 0] + L, wrapped[..., 0]),
    )
    if periodic_y:
        wrapped[..., 1] = np.where(
            wrapped[..., 1] > half_L,
            wrapped[..., 1] - L,
            np.where(wrapped[..., 1] < -half_L, wrapped[..., 1] + L, wrapped[..., 1]),
        )
    return wrapped


def periodic_distance(
    source: np.ndarray,
    targets: np.ndarray,
    L: float,
    *,
    periodic_y: bool = False,
) -> np.ndarray:
    """Return Euclidean minimum-image distance for the selected boundaries."""
    delta = minimum_image_delta(source, targets, L, periodic_y=periodic_y)
    return np.sqrt(np.sum(delta * delta, axis=-1))


def infection_probability(
    distances: np.ndarray,
    *,
    source_is_superspreader: bool,
    model: str,
    r0: float,
    w0: float,
) -> np.ndarray:
    """Return infection probabilities for one infectious source and many targets."""
    if model not in VALID_MODELS:
        raise ValueError(f"model must be one of {sorted(VALID_MODELS)}")
    if r0 <= 0:
        raise ValueError("r0 must be positive")
    if not 0 <= w0 <= 1:
        raise ValueError("w0 must be in [0, 1]")

    distances_array = np.asarray(distances, dtype=float)

    if model == "strong" and source_is_superspreader:
        probabilities = np.where(distances_array <= r0, w0, 0.0)
    elif model == "hub" and source_is_superspreader:
        rn = np.sqrt(6.0) * r0
        probabilities = np.where(distances_array <= rn, w0 * (1.0 - distances_array / rn) ** 2, 0.0)
    else:
        probabilities = np.where(distances_array <= r0, w0 * (1.0 - distances_array / r0) ** 2, 0.0)

    return np.clip(probabilities, 0.0, 1.0)


@dataclass
class SimulationResult:
    """Final state and infection history from one simulation run."""

    positions: np.ndarray
    states: np.ndarray
    is_superspreader: np.ndarray
    infected_time: np.ndarray
    recovered_time: np.ndarray
    infection_source: np.ndarray
    secondary_counts: np.ndarray
    new_infections_per_step: list[int]
    total_infected: int
    duration: int


class MonteCarloSIRSimulator:
    """Monte Carlo SIR simulator for normal, strong, and hub superspreader models."""

    def __init__(
        self,
        *,
        N: int = DEFAULT_N,
        model: str = "normal",
        lambda_ss: float = 0.0,
        r0: float = DEFAULT_R0,
        L: Optional[float] = None,
        w0: float = DEFAULT_W0,
        gamma: float = DEFAULT_GAMMA,
        max_steps: int = DEFAULT_MAX_STEPS,
        seed: Optional[int] = None,
        positions: Optional[np.ndarray] = None,
        is_superspreader: Optional[np.ndarray] = None,
    ) -> None:
        """Initialize population, superspreader flags, and SIR state."""
        if positions is None and N <= 0:
            raise ValueError("N must be positive")

        self.model = model
        self.lambda_ss = float(lambda_ss)
        self.r0 = float(r0)
        self.L = resolve_L(self.r0, L)
        self.w0 = float(w0)
        self.gamma = float(gamma)
        self.max_steps = int(max_steps)
        self.rng = np.random.default_rng(seed)

        if positions is not None:
            positions_array = np.asarray(positions, dtype=float)
            if positions_array.ndim != 2 or positions_array.shape[1] != 2:
                raise ValueError("positions must have shape (N, 2)")
            self.N = int(positions_array.shape[0])
            self.positions = positions_array.copy()
        else:
            self.N = int(N)
            self.positions = self.rng.uniform(0.0, self.L, size=(self.N, 2))

        self._validate_parameters()
        self._validate_positions()

        self.positions[0] = np.array([0.5 * self.L, 0.0])
        self.is_superspreader = self._initialize_superspreaders(is_superspreader)

        self.states = np.full(self.N, SUSCEPTIBLE, dtype=int)
        self.states[0] = INFECTED

        self.infected_time = np.full(self.N, np.nan, dtype=float)
        self.infected_time[0] = 0.0
        self.recovered_time = np.full(self.N, np.nan, dtype=float)
        self.infection_source = np.full(self.N, -1, dtype=int)
        self.secondary_counts = np.zeros(self.N, dtype=int)
        self.time = 0
        self.new_infections_per_step: list[int] = []

    def _validate_parameters(self) -> None:
        if self.N <= 0:
            raise ValueError("N must be positive")
        if self.model not in VALID_MODELS:
            raise ValueError(f"model must be one of {sorted(VALID_MODELS)}")
        if not 0 <= self.lambda_ss <= 1:
            raise ValueError("lambda_ss must be in [0, 1]")
        if self.r0 <= 0:
            raise ValueError("r0 must be positive")
        if self.L <= 0:
            raise ValueError("L must be positive")
        if not 0 <= self.w0 <= 1:
            raise ValueError("w0 must be in [0, 1]")
        if not 0 <= self.gamma <= 1:
            raise ValueError("gamma must be in [0, 1]")
        if self.max_steps < 0:
            raise ValueError("max_steps must be non-negative")

    def _validate_positions(self) -> None:
        if self.positions.shape != (self.N, 2):
            raise ValueError("positions must have shape (N, 2)")
        if np.any(self.positions < 0) or np.any(self.positions >= self.L):
            raise ValueError("positions must lie in [0, L)")

    def _initialize_superspreaders(self, is_superspreader: Optional[np.ndarray]) -> np.ndarray:
        if is_superspreader is not None:
            flags = np.asarray(is_superspreader, dtype=bool)
            if flags.shape != (self.N,):
                raise ValueError("is_superspreader must have shape (N,)")
            return flags.copy()

        flags = np.zeros(self.N, dtype=bool)
        if self.model == "normal" or self.lambda_ss == 0:
            return flags

        count = int(round(self.lambda_ss * self.N))
        if count > 0:
            selected = self.rng.choice(self.N, size=count, replace=False)
            flags[selected] = True
        return flags

    def step(self) -> int:
        """Run one Monte Carlo timestep and return the number of new infections."""
        infected_at_start = np.flatnonzero(self.states == INFECTED)
        next_time = self.time + 1
        new_infected_count = 0

        for source_id in self.rng.permutation(infected_at_start):
            susceptible_ids = np.flatnonzero(self.states == SUSCEPTIBLE)
            if susceptible_ids.size > 0:
                distances = periodic_distance(
                    self.positions[source_id],
                    self.positions[susceptible_ids],
                    self.L,
                )
                probabilities = infection_probability(
                    distances,
                    source_is_superspreader=bool(self.is_superspreader[source_id]),
                    model=self.model,
                    r0=self.r0,
                    w0=self.w0,
                )
                infected_mask = self.rng.random(susceptible_ids.size) < probabilities
                infected_target_ids = susceptible_ids[infected_mask]

                for target_id in infected_target_ids:
                    self.states[target_id] = INFECTED
                    self.infected_time[target_id] = float(next_time)
                    self.infection_source[target_id] = int(source_id)
                    self.secondary_counts[source_id] += 1
                    new_infected_count += 1

            if self.rng.random() < self.gamma:
                self.states[source_id] = RECOVERED
                self.recovered_time[source_id] = float(next_time)

        self.time = next_time
        self.new_infections_per_step.append(new_infected_count)
        return new_infected_count

    def run(self) -> SimulationResult:
        """Run until no infected individuals remain or max_steps is reached."""
        while self.time < self.max_steps and np.any(self.states == INFECTED):
            self.step()

        return SimulationResult(
            positions=self.positions.copy(),
            states=self.states.copy(),
            is_superspreader=self.is_superspreader.copy(),
            infected_time=self.infected_time.copy(),
            recovered_time=self.recovered_time.copy(),
            infection_source=self.infection_source.copy(),
            secondary_counts=self.secondary_counts.copy(),
            new_infections_per_step=list(self.new_infections_per_step),
            total_infected=int(np.count_nonzero(np.isfinite(self.infected_time))),
            duration=int(self.time),
        )
