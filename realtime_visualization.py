"""Realtime visualization for the spatial Monte Carlo SIR simulator."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
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
)
from epidemic_model.simulator import MonteCarloSIRSimulator


ROUTE_COLORS = {
    "edge": "#94A3B8",
    "susceptible_normal": "#CBD5E1",
    "susceptible_superspreader": "#F59E0B",
    "infected_normal": "#2563EB",
    "infected_superspreader": "#E11D48",
    "recovered_normal": "#93C5FD",
    "recovered_superspreader": "#FCA5A5",
    "initial": "#111827",
}


@dataclass(frozen=True)
class SimulationSettings:
    """Runtime parameters controlled by the visualization."""

    N: int = DEFAULT_N
    model: str = "strong"
    lambda_ss: float = 0.2
    r0: float = DEFAULT_R0
    L: float = 10.0 * DEFAULT_R0
    w0: float = DEFAULT_W0
    gamma: float = DEFAULT_GAMMA
    max_steps: int = DEFAULT_MAX_STEPS
    seed: Optional[int] = 42
    steps_per_tick: int = 1
    autoplay: bool = True


def wrapped_route_segments(
    positions: np.ndarray,
    infection_source: np.ndarray,
    L: float,
) -> list[list[tuple[float, float]]]:
    """Return infection tree line segments, split at the horizontal wrap boundary."""
    segments: list[list[tuple[float, float]]] = []
    for target_id, source_id in enumerate(infection_source):
        if source_id < 0:
            continue

        source_x, source_y = positions[source_id]
        target_x, target_y = positions[target_id]
        delta_x = target_x - source_x
        if abs(delta_x) <= L / 2.0:
            segments.append([(float(source_x), float(source_y)), (float(target_x), float(target_y))])
            continue

        if delta_x > 0:
            target_x_unwrapped = target_x - L
            boundary_x = 0.0
            wrapped_boundary_x = L
        else:
            target_x_unwrapped = target_x + L
            boundary_x = L
            wrapped_boundary_x = 0.0

        fraction = (boundary_x - source_x) / (target_x_unwrapped - source_x)
        boundary_y = source_y + fraction * (target_y - source_y)
        segments.append([(float(source_x), float(source_y)), (float(boundary_x), float(boundary_y))])
        segments.append([(float(wrapped_boundary_x), float(boundary_y)), (float(target_x), float(target_y))])
    return segments


def validate_settings(settings: SimulationSettings) -> None:
    """Validate settings before starting an interactive run."""
    if settings.N <= 0:
        raise ValueError("N must be positive")
    if settings.model not in VALID_MODELS:
        raise ValueError(f"model must be one of {sorted(VALID_MODELS)}")
    if not 0.0 <= settings.lambda_ss <= 1.0:
        raise ValueError("lambda must be in [0, 1]")
    if settings.r0 <= 0:
        raise ValueError("r0 must be positive")
    if settings.L <= 0:
        raise ValueError("L must be positive")
    if not 0.0 <= settings.w0 <= 1.0:
        raise ValueError("w0 must be in [0, 1]")
    if not 0.0 <= settings.gamma <= 1.0:
        raise ValueError("gamma must be in [0, 1]")
    if settings.max_steps < 0:
        raise ValueError("max_steps must be non-negative")
    if settings.steps_per_tick <= 0:
        raise ValueError("steps_per_tick must be positive")


class LiveSIRVisualization:
    """Matplotlib widget UI that animates one simulation run."""

    def __init__(self, settings: SimulationSettings) -> None:
        validate_settings(settings)

        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
        from matplotlib.collections import LineCollection
        from matplotlib.widgets import Button, RadioButtons, Slider

        self.plt = plt
        self.FuncAnimation = FuncAnimation
        self.LineCollection = LineCollection
        self.Button = Button
        self.RadioButtons = RadioButtons
        self.Slider = Slider

        self.settings = settings
        self.simulator = self._new_simulator(settings)
        self.running = settings.autoplay
        self.last_new_infections = 1
        self.times = [0]
        self.new_counts = [1]
        self.active_counts = [1]
        self.cumulative_counts = [1]

        self._setup_figure()
        self._draw()
        self.animation = self.FuncAnimation(
            self.figure,
            self._tick,
            interval=180,
            blit=False,
            cache_frame_data=False,
        )

    def _new_simulator(self, settings: SimulationSettings) -> MonteCarloSIRSimulator:
        return MonteCarloSIRSimulator(
            N=settings.N,
            model=settings.model,
            lambda_ss=settings.lambda_ss,
            r0=settings.r0,
            L=settings.L,
            w0=settings.w0,
            gamma=settings.gamma,
            max_steps=settings.max_steps,
            seed=settings.seed,
        )

    def _setup_figure(self) -> None:
        self.figure = self.plt.figure(figsize=(12.8, 7.4))
        self.route_ax = self.figure.add_axes([0.06, 0.11, 0.58, 0.82])
        self.curve_ax = self.figure.add_axes([0.71, 0.57, 0.25, 0.32])
        self.stats_ax = self.figure.add_axes([0.71, 0.42, 0.25, 0.11])
        self.stats_ax.axis("off")

        self.edge_collection = self.LineCollection(
            [],
            colors=ROUTE_COLORS["edge"],
            linewidths=0.65,
            alpha=0.52,
            zorder=1,
        )
        self.route_ax.add_collection(self.edge_collection)

        self.susceptible_normal = self.route_ax.scatter(
            [],
            [],
            s=18,
            facecolors="white",
            edgecolors=ROUTE_COLORS["susceptible_normal"],
            linewidths=0.55,
            label="S (normal)",
            zorder=2,
        )
        self.susceptible_superspreader = self.route_ax.scatter(
            [],
            [],
            s=42,
            facecolors="white",
            edgecolors=ROUTE_COLORS["susceptible_superspreader"],
            linewidths=0.9,
            label="S (superspreader)",
            zorder=3,
        )
        self.infected_normal = self.route_ax.scatter(
            [],
            [],
            s=28,
            color=ROUTE_COLORS["infected_normal"],
            linewidths=0,
            label="I (normal)",
            zorder=4,
        )
        self.infected_superspreader = self.route_ax.scatter(
            [],
            [],
            s=58,
            color=ROUTE_COLORS["infected_superspreader"],
            linewidths=0,
            label="I (superspreader)",
            zorder=5,
        )
        self.recovered_normal = self.route_ax.scatter(
            [],
            [],
            s=20,
            color=ROUTE_COLORS["recovered_normal"],
            alpha=0.62,
            linewidths=0,
            label="R (normal)",
            zorder=3,
        )
        self.recovered_superspreader = self.route_ax.scatter(
            [],
            [],
            s=48,
            color=ROUTE_COLORS["recovered_superspreader"],
            alpha=0.72,
            linewidths=0,
            label="R (superspreader)",
            zorder=4,
        )
        self.initial_source = self.route_ax.scatter(
            [],
            [],
            s=80,
            marker="*",
            color=ROUTE_COLORS["initial"],
            linewidths=0,
            label="initial",
            zorder=6,
        )

        self.route_ax.set_title("route of infection")
        self.route_ax.set_aspect("equal", adjustable="box")
        self.route_ax.grid(True, color="#E5E7EB", linewidth=0.7)
        self.route_ax.legend(
            frameon=True,
            framealpha=0.9,
            facecolor="white",
            edgecolor="0.86",
            fontsize=7.2,
            handlelength=1.0,
            loc="lower left",
            borderpad=0.35,
            labelspacing=0.25,
        )

        (self.new_line,) = self.curve_ax.plot([], [], color="#DC2626", linewidth=1.4, label="new")
        (self.active_line,) = self.curve_ax.plot([], [], color="#2563EB", linewidth=1.4, label="active")
        (self.cumulative_line,) = self.curve_ax.plot([], [], color="#111827", linewidth=1.3, label="total")
        self.curve_ax.set_title("epidemic curve", fontsize=10)
        self.curve_ax.set_xlabel("time step")
        self.curve_ax.set_ylabel("people")
        self.curve_ax.grid(True, color="#E5E7EB", linewidth=0.7)
        self.curve_ax.legend(frameon=False, fontsize=7.4, loc="upper left")

        model_axis = self.figure.add_axes([0.71, 0.25, 0.12, 0.12])
        model_labels = ("normal", "strong", "hub")
        model_index = model_labels.index(self.settings.model)
        self.model_radio = self.RadioButtons(model_axis, model_labels, active=model_index)

        self.N_slider = self._slider(
            [0.76, 0.220, 0.20, 0.020],
            "N",
            10,
            max(1500, self.settings.N),
            self.settings.N,
            1,
        )
        self.lambda_slider = self._slider(
            [0.76, 0.190, 0.20, 0.020],
            "lambda",
            0.0,
            1.0,
            self.settings.lambda_ss,
            0.01,
        )
        self.L_slider = self._slider(
            [0.76, 0.160, 0.20, 0.020],
            "L",
            1.0,
            max(30.0, self.settings.L),
            self.settings.L,
            0.1,
        )
        self.r0_slider = self._slider(
            [0.76, 0.130, 0.20, 0.020],
            "r0",
            0.1,
            max(5.0, self.settings.r0),
            self.settings.r0,
            0.05,
        )
        self.w0_slider = self._slider(
            [0.76, 0.100, 0.20, 0.020],
            "w0",
            0.0,
            1.0,
            self.settings.w0,
            0.01,
        )
        self.gamma_slider = self._slider(
            [0.76, 0.070, 0.20, 0.020],
            "gamma",
            0.0,
            1.0,
            self.settings.gamma,
            0.01,
        )
        self.speed_slider = self._slider(
            [0.76, 0.040, 0.20, 0.020],
            "speed",
            1,
            max(10, self.settings.steps_per_tick),
            self.settings.steps_per_tick,
            1,
        )

        self.run_button = self._button([0.68, 0.006, 0.06, 0.030], "Pause" if self.running else "Run", self._on_run)
        self.step_button = self._button([0.52, 0.006, 0.06, 0.030], "Step", self._on_step)
        self.apply_button = self._button([0.60, 0.006, 0.06, 0.030], "Apply", self._on_apply)
        self.seed_button = self._button([0.44, 0.006, 0.06, 0.030], "Seed", self._on_new_seed)

    def _slider(
        self,
        rect: list[float],
        label: str,
        valmin: float,
        valmax: float,
        valinit: float,
        valstep: float,
    ):
        axis = self.figure.add_axes(rect)
        return self.Slider(axis, label, valmin, valmax, valinit=valinit, valstep=valstep)

    def _button(self, rect: list[float], label: str, callback):
        axis = self.figure.add_axes(rect)
        button = self.Button(axis, label)
        button.on_clicked(callback)
        return button

    def _settings_from_controls(self) -> SimulationSettings:
        return replace(
            self.settings,
            N=int(round(self.N_slider.val)),
            model=str(self.model_radio.value_selected),
            lambda_ss=float(self.lambda_slider.val),
            L=float(self.L_slider.val),
            r0=float(self.r0_slider.val),
            w0=float(self.w0_slider.val),
            gamma=float(self.gamma_slider.val),
            steps_per_tick=int(round(self.speed_slider.val)),
            autoplay=self.running,
        )

    def _reset(self, settings: SimulationSettings) -> None:
        validate_settings(settings)
        self.settings = settings
        self.simulator = self._new_simulator(settings)
        self.last_new_infections = 1
        self.times = [0]
        self.new_counts = [1]
        self.active_counts = [1]
        self.cumulative_counts = [1]
        self._draw()

    def _is_finished(self) -> bool:
        return self.simulator.time >= self.simulator.max_steps or not np.any(self.simulator.states == INFECTED)

    def _on_run(self, _event) -> None:
        self.running = not self.running
        self.run_button.label.set_text("Pause" if self.running else "Run")

    def _on_step(self, _event) -> None:
        if self._step_once():
            self._draw()

    def _on_apply(self, _event) -> None:
        self._reset(self._settings_from_controls())

    def _on_new_seed(self, _event) -> None:
        next_seed = 1 if self.settings.seed is None else self.settings.seed + 1
        self._reset(replace(self._settings_from_controls(), seed=next_seed))

    def _step_once(self) -> bool:
        if self._is_finished():
            self.running = False
            self.run_button.label.set_text("Run")
            return False

        self.last_new_infections = self.simulator.step()
        self.times.append(self.simulator.time)
        self.new_counts.append(self.last_new_infections)
        self.active_counts.append(int(np.count_nonzero(self.simulator.states == INFECTED)))
        self.cumulative_counts.append(int(np.count_nonzero(np.isfinite(self.simulator.infected_time))))
        return True

    def _tick(self, _frame):
        if not self.running:
            return []

        changed = False
        for _ in range(max(1, int(round(self.speed_slider.val)))):
            changed = self._step_once() or changed
            if self._is_finished():
                break
        if changed:
            self._draw()
        return []

    def _draw(self) -> None:
        positions = self.simulator.positions
        states = self.simulator.states
        superspreaders = self.simulator.is_superspreader
        normal = ~superspreaders

        self.edge_collection.set_segments(
            wrapped_route_segments(positions, self.simulator.infection_source, self.simulator.L)
        )
        self._set_offsets(self.susceptible_normal, positions, (states == SUSCEPTIBLE) & normal)
        self._set_offsets(self.susceptible_superspreader, positions, (states == SUSCEPTIBLE) & superspreaders)
        self._set_offsets(self.infected_normal, positions, (states == INFECTED) & normal)
        self._set_offsets(self.infected_superspreader, positions, (states == INFECTED) & superspreaders)
        self._set_offsets(self.recovered_normal, positions, (states == RECOVERED) & normal)
        self._set_offsets(self.recovered_superspreader, positions, (states == RECOVERED) & superspreaders)
        self.initial_source.set_offsets(positions[[0]])

        self.route_ax.set_xlim(0.0, self.simulator.L)
        self.route_ax.set_ylim(0.0, self.simulator.L)
        self.route_ax.set_title(
            f"route of infection: {self.settings.model}, lambda={self.settings.lambda_ss:.2f}, "
            f"N={self.settings.N}, L={self.settings.L:.1f}"
        )

        self.new_line.set_data(self.times, self.new_counts)
        self.active_line.set_data(self.times, self.active_counts)
        self.cumulative_line.set_data(self.times, self.cumulative_counts)
        x_max = max(5, self.settings.max_steps, self.times[-1])
        y_max = max(1, max(self.new_counts), max(self.active_counts), max(self.cumulative_counts))
        self.curve_ax.set_xlim(0, x_max)
        self.curve_ax.set_ylim(0, y_max * 1.08)

        susceptible_count = int(np.count_nonzero(states == SUSCEPTIBLE))
        active_count = int(np.count_nonzero(states == INFECTED))
        recovered_count = int(np.count_nonzero(states == RECOVERED))
        total_infected = int(np.count_nonzero(np.isfinite(self.simulator.infected_time)))
        max_secondary = int(np.max(self.simulator.secondary_counts)) if self.simulator.secondary_counts.size else 0
        self.stats_ax.clear()
        self.stats_ax.axis("off")
        self.stats_ax.text(
            0.0,
            1.0,
            "\n".join(
                (
                    f"time: {self.simulator.time}/{self.settings.max_steps}",
                    f"S: {susceptible_count}   I: {active_count}   R: {recovered_count}",
                    f"new: {self.last_new_infections}   total: {total_infected}",
                    f"attack rate: {total_infected / self.settings.N:.3f}",
                    f"max secondary: {max_secondary}",
                    f"seed: {self.settings.seed}",
                )
            ),
            va="top",
            ha="left",
            fontsize=9,
            family="monospace",
        )
        self.figure.canvas.draw_idle()

    def _set_offsets(self, artist, positions: np.ndarray, mask: np.ndarray) -> None:
        offsets = positions[mask]
        artist.set_offsets(offsets if offsets.size else np.empty((0, 2)))

    def show(self) -> None:
        self.plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=DEFAULT_N, help="Number of people in the simulation.")
    parser.add_argument(
        "--model",
        choices=sorted(VALID_MODELS),
        default="strong",
        help="Infection model to visualize.",
    )
    parser.add_argument(
        "--lambda",
        dest="lambda_ss",
        type=float,
        default=0.2,
        help="Fraction of superspreaders for strong and hub models.",
    )
    parser.add_argument("--L", type=float, default=10.0 * DEFAULT_R0, help="Side length of the square simulation box.")
    parser.add_argument("--r0", type=float, default=DEFAULT_R0, help="Normal infection radius.")
    parser.add_argument("--w0", type=float, default=DEFAULT_W0, help="Maximum infection probability.")
    parser.add_argument("--gamma", type=float, default=DEFAULT_GAMMA, help="Recovery probability per timestep.")
    parser.add_argument("--max-steps", type=int, default=DEFAULT_MAX_STEPS, help="Maximum timesteps.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--speed", type=int, default=1, help="Simulation steps per animation tick.")
    parser.add_argument("--paused", action="store_true", help="Open the window without autoplay.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    settings = SimulationSettings(
        N=args.N,
        model=args.model,
        lambda_ss=args.lambda_ss,
        L=args.L,
        r0=args.r0,
        w0=args.w0,
        gamma=args.gamma,
        max_steps=args.max_steps,
        seed=args.seed,
        steps_per_tick=args.speed,
        autoplay=not args.paused,
    )
    app = LiveSIRVisualization(settings)
    app.show()


if __name__ == "__main__":
    main()
