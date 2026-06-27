"""Plot generation and saved route reconstruction."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from epidemic_model.simulator import SimulationResult

from .files import read_csv
from .metrics import run_simulation, top_band_reached
from .records import SimulationConfig
from .settings import (
    DEFAULT_SUPERSPREADER_FRACTION,
    INFECTION_SOURCE_COLORS,
    SUPERSPREADER_FRACTION_COLORS,
    SUPERSPREADER_FRACTION_MARKERS,
    TRANSMISSION_MODEL_COLORS,
    NORMALIZED_DENSITY_FACTOR,
    OUTPUT_PLOT_FILES,
    PLOT_LINE_WIDTH,
    PLOT_MARKER_SIZE,
    PLOT_RESOLUTION_DPI,
    ROUTE_COLORS,
    ROUTE_PLOT_SIZE,
    SARS_BAR_COLOR,
    SARS_SUPERSPREADER_FRACTION,
    SPREAD_ANALYSIS_DENSITY,
    SPREAD_ANALYSIS_POPULATION,
    STANDARD_PLOT_SIZE,
)


def apply_plot_axes() -> None:
    axis = plt.gca()
    axis.grid(False)
    axis.tick_params(direction="in", top=False, right=False, length=4, width=0.8)
    for spine in axis.spines.values():
        spine.set_linewidth(0.8)


def finish_plot(
    path: Path,
    *,
    legend: bool = True,
    legend_loc: str = "best",
    legend_ncol: int = 1,
    legend_fontsize: float = 7.6,
    legend_bbox: tuple[float, float] | None = None,
) -> None:
    apply_plot_axes()
    if legend:
        plt.legend(
            loc=legend_loc,
            bbox_to_anchor=legend_bbox,
            ncol=legend_ncol,
            frameon=True,
            framealpha=0.9,
            facecolor="white",
            edgecolor="0.86",
            fontsize=legend_fontsize,
            handlelength=2.0,
            borderpad=0.45,
            labelspacing=0.35,
            columnspacing=0.9,
        )
    plt.tight_layout()
    plt.savefig(path, dpi=PLOT_RESOLUTION_DPI, bbox_inches="tight", pad_inches=0.03)
    plt.close()


def plot_infection_kernel(rows: list[dict[str, object]], path: Path, transmission_model: str) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    model_rows = [row for row in rows if row["transmission_model"] == transmission_model]
    for source_type, style in (
        (
            "normal",
            {
                "linestyle": "--",
                "color": INFECTION_SOURCE_COLORS["normal"],
                "label": "Normal",
                "linewidth": PLOT_LINE_WIDTH + 0.2,
            },
        ),
        (
            "superspreader",
            {
                "linestyle": "-",
                "color": INFECTION_SOURCE_COLORS["superspreader"],
                "label": "Superspreader",
                "linewidth": PLOT_LINE_WIDTH + 0.2,
            },
        ),
    ):
        source_rows = [row for row in model_rows if row["source_type"] == source_type]
        distances = np.array([float(row["distance_over_r0"]) for row in source_rows])
        probabilities = np.array([float(row["probability_over_w0"]) for row in source_rows])
        plt.plot(distances, probabilities, **style)

    plt.xlabel(r"$r / r_0$")
    plt.ylabel(r"$w(r) / w_0$")
    plt.xlim(0.0, 1.0 if transmission_model == "strong" else np.sqrt(6.0))
    plt.ylim(-0.05, 1.05)
    finish_plot(path, legend_loc="upper right")


def plot_percolation(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    relevant_rows = [
        row
        for row in rows
        if row["transmission_model"] == "normal"
        or abs(float(row["superspreader_fraction"]) - DEFAULT_SUPERSPREADER_FRACTION) < 1e-9
    ]
    style_by_model = {
        "normal": {"marker": "^", "linestyle": ":", "label": "No superspreaders"},
        "strong": {"marker": "o", "linestyle": "-", "label": "Strong"},
        "hub": {"marker": "s", "linestyle": "--", "label": "Hub"},
    }
    for transmission_model in ("normal", "strong", "hub"):
        model_rows = [row for row in relevant_rows if row["transmission_model"] == transmission_model]
        densities = np.array([float(row["population_density"]) for row in model_rows]) * NORMALIZED_DENSITY_FACTOR
        probabilities = np.array([float(row["percolation_probability"]) for row in model_rows])
        order = np.argsort(densities)
        style = style_by_model[transmission_model]
        plt.plot(
            densities[order],
            probabilities[order],
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            marker=style["marker"],
            linestyle=style["linestyle"],
            markerfacecolor="white",
            markeredgecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            markersize=PLOT_MARKER_SIZE,
            linewidth=PLOT_LINE_WIDTH,
            label=style["label"],
        )

    plt.xlabel(r"$\rho \pi r_0^2$")
    plt.ylabel("Percolation probability")
    plt.xlim(0.0, 25.0)
    plt.ylim(-0.05, 1.05)
    finish_plot(path, legend_loc="lower right")


def plot_single_model_density_sweep(rows: list[dict[str, object]], path: Path, transmission_model: str) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    model_rows = [row for row in rows if row["transmission_model"] == transmission_model]
    for superspreader_fraction in sorted({float(row["superspreader_fraction"]) for row in model_rows}):
        fraction_rows = [row for row in model_rows if abs(float(row["superspreader_fraction"]) - superspreader_fraction) < 1e-9]
        densities = np.array([float(row["population_density"]) for row in fraction_rows]) * NORMALIZED_DENSITY_FACTOR
        probabilities = np.array([float(row["percolation_probability"]) for row in fraction_rows])
        order = np.argsort(densities)
        marker = SUPERSPREADER_FRACTION_MARKERS.get(round(superspreader_fraction, 1), "o")
        color = SUPERSPREADER_FRACTION_COLORS.get(round(superspreader_fraction, 1), "#111827")
        plt.plot(
            densities[order],
            probabilities[order],
            linestyle="-",
            marker=marker,
            color=color,
            linewidth=PLOT_LINE_WIDTH,
            markersize=PLOT_MARKER_SIZE,
            markerfacecolor="white",
            markeredgecolor=color,
            label=f"fraction={superspreader_fraction:.1f}",
        )

    plt.xlabel(r"$\rho \pi r_0^2$")
    plt.ylabel("Percolation probability")
    plt.ylim(-0.05, 1.05)
    plt.xlim(0.0, 25.0)
    finish_plot(
        path,
        legend_loc="upper center",
        legend_ncol=3,
        legend_fontsize=7.0,
        legend_bbox=(0.5, -0.16),
    )


def plot_critical_density(
    rows: list[dict[str, object]],
    reference_rows: list[dict[str, object]],
    path: Path,
) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    style_by_model = {
        "strong": {"marker": "o", "label": "Strong sim."},
        "hub": {"marker": "s", "label": "Hub sim."},
    }
    for transmission_model in ("strong", "hub"):
        model_rows = [row for row in rows if row["transmission_model"] == transmission_model and str(row["critical_density"]).strip()]
        fractions = np.array([float(row["superspreader_fraction"]) for row in model_rows])
        densities = np.array([float(row["critical_density"]) for row in model_rows]) * NORMALIZED_DENSITY_FACTOR
        order = np.argsort(fractions)
        style = style_by_model[transmission_model]
        plt.plot(
            fractions[order],
            densities[order],
            linestyle="None",
            marker=style["marker"],
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            markerfacecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            markeredgecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            markersize=PLOT_MARKER_SIZE + 0.5,
            label=style["label"],
        )

    for transmission_model, linestyle in (("strong", "-"), ("hub", "--")):
        model_reference_rows = [row for row in reference_rows if row["transmission_model"] == transmission_model]
        fractions = np.array([float(row["superspreader_fraction"]) for row in model_reference_rows])
        densities = np.array(
            [float(row["normalized_critical_density"]) for row in model_reference_rows]
        )
        order = np.argsort(fractions)
        plt.plot(
            fractions[order],
            densities[order],
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            linestyle=linestyle,
            linewidth=PLOT_LINE_WIDTH,
            alpha=0.85,
            label=rf"{transmission_model.title()} $R_0=R_c$",
        )

    plt.xlabel("Superspreader fraction")
    plt.ylabel(r"$\rho_c \pi r_0^2$")
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 25.0)
    finish_plot(path, legend_loc="upper right", legend_fontsize=7.0)


def plot_epidemic_curves(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    style_by_model = {
        "strong": {"marker": "o", "label": f"Strong (fraction={DEFAULT_SUPERSPREADER_FRACTION:.1f})"},
        "hub": {"marker": "s", "label": f"Hub (fraction={DEFAULT_SUPERSPREADER_FRACTION:.1f})"},
        "normal": {"marker": "^", "label": "No superspreaders"},
    }
    for transmission_model in ("strong", "hub", "normal"):
        model_rows = [row for row in rows if row["transmission_model"] == transmission_model]
        times = np.array([int(row["time"]) for row in model_rows])
        new_cases = np.array([float(row["mean_new_cases"]) for row in model_rows])
        style = style_by_model[transmission_model]
        plt.plot(
            times,
            new_cases,
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            linewidth=PLOT_LINE_WIDTH,
            marker=style["marker"],
            markersize=PLOT_MARKER_SIZE,
            markerfacecolor="white",
            markeredgecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            markevery=2,
            label=style["label"],
        )

    plt.xlabel("time step")
    plt.ylabel("the number of infected")
    plt.xlim(0.0, 40.0)
    finish_plot(path, legend_loc="upper right", legend_fontsize=7.4)


def plot_secondary_distribution(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    style_by_model = {
        "normal": {"marker": "^", "linestyle": ":", "label": "No superspreaders"},
        "strong": {"marker": "o", "linestyle": "-", "label": "Strong"},
        "hub": {"marker": "s", "linestyle": "--", "label": "Hub"},
    }
    for transmission_model in ("normal", "strong", "hub"):
        model_rows = [row for row in rows if row["transmission_model"] == transmission_model]
        counts = np.array([int(row["secondary_cases"]) for row in model_rows])
        probabilities = np.array([float(row["probability"]) for row in model_rows])
        mask = probabilities > 0
        style = style_by_model[transmission_model]
        plt.semilogy(
            counts[mask],
            probabilities[mask],
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            linestyle=style["linestyle"],
            marker=style["marker"],
            markersize=PLOT_MARKER_SIZE,
            markerfacecolor="white",
            markeredgecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            linewidth=PLOT_LINE_WIDTH,
            label=style["label"],
        )

    plt.xlabel("Secondary infections caused by one individual")
    plt.ylabel("Probability")
    finish_plot(path, legend_loc="upper right")


def plot_secondary_distribution_normal(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    model_rows = [row for row in rows if row["transmission_model"] == "normal"]
    counts = np.array([int(row["secondary_cases"]) for row in model_rows])
    probabilities = np.array([float(row["probability"]) for row in model_rows])
    mask = probabilities > 0
    plt.plot(
        counts[mask],
        probabilities[mask],
        marker="^",
        color=TRANSMISSION_MODEL_COLORS["normal"],
        markerfacecolor="white",
        markeredgecolor=TRANSMISSION_MODEL_COLORS["normal"],
        markersize=PLOT_MARKER_SIZE,
        linewidth=PLOT_LINE_WIDTH,
        label="No superspreaders",
    )

    plt.xlabel("the number of links")
    plt.ylabel("Probability")
    plt.xlim(0.0, 20.0)
    plt.ylim(0.0, 0.85)
    finish_plot(path, legend_loc="upper right")


def plot_secondary_distribution_superspreaders(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    for transmission_model, marker, linestyle in (("strong", "o", "-"), ("hub", "s", "--")):
        model_rows = [row for row in rows if row["transmission_model"] == transmission_model]
        counts = np.array([int(row["secondary_cases"]) for row in model_rows])
        probabilities = np.array([float(row["probability"]) for row in model_rows])
        mask = probabilities > 0
        plt.plot(
            counts[mask],
            probabilities[mask],
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            linestyle=linestyle,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            markersize=PLOT_MARKER_SIZE,
            linewidth=PLOT_LINE_WIDTH,
            label="Strong" if transmission_model == "strong" else "Hub",
        )

    plt.xlabel("the number of links")
    plt.ylabel("Probability")
    plt.xlim(0.0, 20.0)
    plt.ylim(0.0, 0.85)
    finish_plot(path, legend_loc="upper right")


def plot_sars_secondary_cases(
    empirical_rows: list[dict[str, object]],
    simulation_rows: list[dict[str, object]],
    path: Path,
) -> None:
    _ = simulation_rows
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    empirical_counts = np.array([int(row["secondary_cases"]) for row in empirical_rows])
    empirical_frequencies = np.array([int(row["frequency"]) for row in empirical_rows])
    plt.bar(
        empirical_counts,
        empirical_frequencies,
        width=0.8,
        color=SARS_BAR_COLOR,
        edgecolor="#92400E",
        hatch="//",
        linewidth=0.8,
    )

    plt.xlabel("number of direct secondary cases")
    plt.ylabel("number")
    plt.xlim(0, 40)
    plt.ylim(0, 180)
    finish_plot(path, legend=False)


def plot_sars_case_curve_comparison(
    empirical_rows: list[dict[str, object]],
    model_rows: list[dict[str, object]],
    path: Path,
) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    time_steps = np.array([int(row["time_step"]) for row in empirical_rows])
    cases = np.array([int(row["new_cases"]) for row in empirical_rows])
    plt.bar(
        time_steps,
        cases,
        width=0.85,
        color=SARS_BAR_COLOR,
        edgecolor="#92400E",
        hatch="//",
        linewidth=0.8,
        label="SARS data (digitized, approximate)",
    )

    for transmission_model, marker, linestyle in (("normal", "^", ":"), ("strong", "o", "-"), ("hub", "s", "--")):
        rows = [row for row in model_rows if row["transmission_model"] == transmission_model]
        times = np.array([int(row["time"]) for row in rows])
        new_cases = np.array([float(row["mean_new_cases"]) for row in rows])
        if transmission_model == "strong":
            label = f"Strong (fraction={SARS_SUPERSPREADER_FRACTION:.1f})"
        elif transmission_model == "hub":
            label = f"Hub (fraction={SARS_SUPERSPREADER_FRACTION:.1f})"
        else:
            label = "No superspreaders"
        plt.plot(
            times,
            new_cases,
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            linestyle=linestyle,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            markersize=PLOT_MARKER_SIZE,
            markevery=3,
            linewidth=PLOT_LINE_WIDTH,
            label=label,
        )

    plt.xlabel("time step")
    plt.ylabel("number of patients")
    plt.xlim(0, 25)
    plt.ylim(0, 80)
    finish_plot(path, legend_loc="upper right", legend_fontsize=7.2)


def spreading_sensitivity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        row
        for row in rows
        if row["scenario_name"] == "superspreader_fraction_sweep"
        and int(row["population_size"]) == SPREAD_ANALYSIS_POPULATION
        and abs(float(row["population_density"]) - SPREAD_ANALYSIS_DENSITY) < 1e-9
    ]


def plot_sensitivity(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    relevant_rows = spreading_sensitivity_rows(rows)
    for transmission_model in ("strong", "hub"):
        model_rows = [row for row in relevant_rows if row["transmission_model"] == transmission_model]
        fractions = np.array([float(row["superspreader_fraction"]) for row in model_rows])
        case_fractions = np.array([float(row["mean_case_fraction"]) for row in model_rows])
        order = np.argsort(fractions)
        marker = "o" if transmission_model == "strong" else "s"
        linestyle = "-" if transmission_model == "strong" else "--"
        plt.plot(
            fractions[order],
            case_fractions[order],
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            linestyle=linestyle,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            markersize=PLOT_MARKER_SIZE,
            linewidth=PLOT_LINE_WIDTH,
            label="Strong" if transmission_model == "strong" else "Hub",
        )

    plt.xlabel("Superspreader fraction")
    plt.ylabel("Mean case fraction")
    plt.ylim(-0.05, 1.05)
    finish_plot(path, legend_loc="lower right")


def plot_velocity(rows: list[dict[str, object]], path: Path) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    relevant_rows = spreading_sensitivity_rows(rows)
    for transmission_model in ("strong", "hub"):
        model_rows = [row for row in relevant_rows if row["transmission_model"] == transmission_model]
        fractions = np.array([float(row["superspreader_fraction"]) for row in model_rows])
        speeds = np.array([float(row["mean_front_speed"]) for row in model_rows])
        order = np.argsort(fractions)
        marker = "o" if transmission_model == "strong" else "s"
        linestyle = "-" if transmission_model == "strong" else "--"
        plt.plot(
            fractions[order],
            speeds[order],
            color=TRANSMISSION_MODEL_COLORS[transmission_model],
            linestyle=linestyle,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=TRANSMISSION_MODEL_COLORS[transmission_model],
            markersize=PLOT_MARKER_SIZE,
            linewidth=PLOT_LINE_WIDTH,
            label="Strong" if transmission_model == "strong" else "Hub",
        )

    plt.xlabel("Superspreader fraction")
    plt.ylabel(r"front speed (/ $r_0 \cdot s$)")
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.6)
    finish_plot(path, legend_loc="upper left")


def plot_front_distance(rows: list[dict[str, object]], path: Path, transmission_model: str) -> None:
    plt.figure(figsize=STANDARD_PLOT_SIZE)
    model_rows = [row for row in rows if row["transmission_model"] == transmission_model]
    for superspreader_fraction in sorted({float(row["superspreader_fraction"]) for row in model_rows}):
        fraction_rows = [row for row in model_rows if abs(float(row["superspreader_fraction"]) - superspreader_fraction) < 1e-9]
        times = np.array([int(row["time"]) for row in fraction_rows])
        front_distance = np.array([float(row["mean_front_distance"]) for row in fraction_rows])
        marker = SUPERSPREADER_FRACTION_MARKERS.get(round(superspreader_fraction, 1), "o")
        color = SUPERSPREADER_FRACTION_COLORS.get(round(superspreader_fraction, 1), "#111827")
        plt.plot(
            times,
            front_distance,
            color=color,
            marker=marker,
            markerfacecolor="white",
            markeredgecolor=color,
            markersize=PLOT_MARKER_SIZE,
            markevery=4,
            linewidth=PLOT_LINE_WIDTH,
            label=f"fraction={superspreader_fraction:.1f}",
        )

    plt.xlabel("time step")
    plt.ylabel(r"$r_f / r_0$")
    plt.xlim(0.0, 40.0)
    plt.ylim(0.0, 12.0)
    finish_plot(
        path,
        legend_loc="upper center",
        legend_ncol=3,
        legend_fontsize=7.0,
        legend_bbox=(0.5, -0.16),
    )


def plot_infection_route(
    results: list[tuple[SimulationConfig, SimulationResult]],
    path: Path,
    transmission_model: str,
) -> None:
    selection = select_representative_route(results, transmission_model)
    if selection is None:
        return
    config, result = selection

    positions = result.positions
    infected_ids = np.flatnonzero(np.isfinite(result.infected_time))
    superspreader_ids = np.flatnonzero(result.is_superspreader)
    infected_mask = np.isfinite(result.infected_time)
    normal_mask = ~result.is_superspreader
    susceptible_normal = normal_mask & ~infected_mask
    susceptible_superspreader = result.is_superspreader & ~infected_mask
    infected_normal = normal_mask & infected_mask
    infected_superspreader = result.is_superspreader & infected_mask

    plt.figure(figsize=ROUTE_PLOT_SIZE)
    for target_id, source_id in enumerate(result.infection_source):
        if source_id < 0:
            continue
        source_x, source_y = positions[source_id]
        target_x, target_y = positions[target_id]
        delta_x = target_x - source_x
        if abs(delta_x) <= config.domain_size / 2.0:
            plt.plot(
                [source_x, target_x],
                [source_y, target_y],
                color=ROUTE_COLORS["edge"],
                alpha=0.55,
                linewidth=0.45,
            )
        else:
            if delta_x > 0:
                target_x_unwrapped = target_x - config.domain_size
                boundary_x = 0.0
                wrapped_boundary_x = config.domain_size
            else:
                target_x_unwrapped = target_x + config.domain_size
                boundary_x = config.domain_size
                wrapped_boundary_x = 0.0
            fraction = (boundary_x - source_x) / (target_x_unwrapped - source_x)
            boundary_y = source_y + fraction * (target_y - source_y)
            plt.plot(
                [source_x, boundary_x],
                [source_y, boundary_y],
                color=ROUTE_COLORS["edge"],
                alpha=0.55,
                linewidth=0.45,
            )
            plt.plot(
                [wrapped_boundary_x, target_x],
                [boundary_y, target_y],
                color=ROUTE_COLORS["edge"],
                alpha=0.55,
                linewidth=0.45,
            )

    plt.scatter(
        positions[susceptible_normal, 0],
        positions[susceptible_normal, 1],
        s=7,
        facecolors="white",
        edgecolors=ROUTE_COLORS["susceptible_normal"],
        linewidths=0.45,
        label="S (normal)",
        zorder=2,
    )
    if superspreader_ids.size:
        plt.scatter(
            positions[susceptible_superspreader, 0],
            positions[susceptible_superspreader, 1],
            s=24,
            facecolors="white",
            edgecolors=ROUTE_COLORS["susceptible_superspreader"],
            linewidths=0.8,
            label="S (superspreader)",
            zorder=3,
        )
    plt.scatter(
        positions[infected_normal, 0],
        positions[infected_normal, 1],
        s=11,
        color=ROUTE_COLORS["infected_normal"],
        linewidths=0,
        label="I (normal)",
        zorder=3,
    )
    if superspreader_ids.size:
        plt.scatter(
            positions[infected_superspreader, 0],
            positions[infected_superspreader, 1],
            s=28,
            facecolors=ROUTE_COLORS["infected_superspreader"],
            edgecolors=ROUTE_COLORS["infected_superspreader"],
            linewidths=0.8,
            label="I (superspreader)",
            zorder=4,
        )
    plt.scatter(positions[0, 0], positions[0, 1], s=44, marker="*", color=ROUTE_COLORS["initial"], zorder=5)
    plt.title("route of infection", fontsize=10)
    plt.xlabel("")
    plt.ylabel("")
    plt.xlim(0, config.domain_size)
    plt.ylim(0, config.domain_size)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.legend(
        frameon=True,
        framealpha=0.88,
        facecolor="white",
        edgecolor="0.86",
        fontsize=6.8,
        handlelength=1.0,
        loc="lower left",
        borderpad=0.35,
        labelspacing=0.25,
    )
    finish_plot(path, legend=False)


def select_representative_route(
    results: list[tuple[SimulationConfig, SimulationResult]],
    transmission_model: str,
) -> tuple[SimulationConfig, SimulationResult] | None:
    """Select the deterministic representative used for one infection-route figure."""
    candidates = [(config, result) for config, result in results if config.transmission_model == transmission_model]
    if not candidates:
        return None
    if transmission_model == "normal":
        non_percolating = [
            (config, result)
            for config, result in candidates
            if not top_band_reached(result, config.domain_size)
        ]
        return max(non_percolating or candidates, key=lambda item: item[1].total_infected)
    return max(candidates, key=lambda item: item[1].total_infected)


def select_representative_routes(
    results: list[tuple[SimulationConfig, SimulationResult]],
) -> list[tuple[SimulationConfig, SimulationResult]]:
    selected = [select_representative_route(results, transmission_model) for transmission_model in ("normal", "strong", "hub")]
    return [item for item in selected if item is not None]


def route_selection_output_rows(
    results: list[tuple[SimulationConfig, SimulationResult]],
) -> list[dict[str, object]]:
    return [
        {
            "scenario_name": config.scenario_name,
            "transmission_model": config.transmission_model,
            "superspreader_fraction": f"{config.superspreader_fraction:.4f}",
            "population_size": config.population_size,
            "population_density": f"{config.population_density:.6f}",
            "domain_size": f"{config.domain_size:.6f}",
            "random_seed": config.random_seed,
            "max_time_steps": config.max_time_steps,
            "total_cases": result.total_infected,
            "reached_top_boundary": int(top_band_reached(result, config.domain_size)),
        }
        for config, result in results
    ]


def route_configs_from_saved_rows(rows: list[dict[str, object]]) -> list[SimulationConfig]:
    configs = [
        SimulationConfig(
            scenario_name=str(row["scenario_name"]),
            transmission_model=str(row["transmission_model"]),
            superspreader_fraction=float(row["superspreader_fraction"]),
            population_size=int(row["population_size"]),
            population_density=float(row["population_density"]),
            random_seed=int(row["random_seed"]),
            max_time_steps=int(row["max_time_steps"]),
        )
        for row in rows
    ]
    if {config.transmission_model for config in configs} != {"normal", "strong", "hub"}:
        raise ValueError("route_plot_selections.csv must contain normal, strong, and hub rows")
    return configs


def require_rows(name: str, rows: list[dict[str, object]] | None) -> list[dict[str, object]]:
    if rows is None:
        raise ValueError(f"{name} rows are required for the selected plot group")
    return rows


def write_selected_plots(
    plots_dir: Path,
    *,
    selected_experiments: set[str],
    infection_rows: list[dict[str, object]] | None = None,
    percolation_rows: list[dict[str, object]] | None = None,
    critical_rows: list[dict[str, object]] | None = None,
    critical_reference_rows: list[dict[str, object]] | None = None,
    front_rows: list[dict[str, object]] | None = None,
    sensitivity_rows: list[dict[str, object]] | None = None,
    curves_rows: list[dict[str, object]] | None = None,
    secondary_rows: list[dict[str, object]] | None = None,
    sars_secondary_rows: list[dict[str, object]] | None = None,
    sars_epicurve_rows: list[dict[str, object]] | None = None,
    sars_model_curve_rows: list[dict[str, object]] | None = None,
    route_results: list[tuple[SimulationConfig, SimulationResult]] | None = None,
) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    for path in plots_dir.glob("*.png"):
        path.unlink()

    if "percolation" in selected_experiments:
        infection_rows = require_rows("infection probability", infection_rows)
        percolation_rows = require_rows("percolation", percolation_rows)
        critical_rows = require_rows("critical density", critical_rows)
        critical_reference_rows = require_rows("critical density reference", critical_reference_rows)
        plot_infection_kernel(infection_rows, plots_dir / OUTPUT_PLOT_FILES["infection_probability_strong"], "strong")
        plot_infection_kernel(infection_rows, plots_dir / OUTPUT_PLOT_FILES["infection_probability_hub"], "hub")
        plot_percolation(percolation_rows, plots_dir / OUTPUT_PLOT_FILES["percolation_probability"])
        plot_single_model_density_sweep(
            percolation_rows,
            plots_dir / OUTPUT_PLOT_FILES["percolation_probability_strong"],
            "strong",
        )
        plot_single_model_density_sweep(
            percolation_rows,
            plots_dir / OUTPUT_PLOT_FILES["percolation_probability_hub"],
            "hub",
        )
        plot_critical_density(
            critical_rows,
            critical_reference_rows,
            plots_dir / OUTPUT_PLOT_FILES["critical_density"],
        )

    if "sensitivity" in selected_experiments:
        front_rows = require_rows("front distance", front_rows)
        sensitivity_rows = require_rows("sensitivity", sensitivity_rows)
        plot_front_distance(front_rows, plots_dir / OUTPUT_PLOT_FILES["front_distance_strong"], "strong")
        plot_front_distance(front_rows, plots_dir / OUTPUT_PLOT_FILES["front_distance_hub"], "hub")
        plot_velocity(sensitivity_rows, plots_dir / OUTPUT_PLOT_FILES["velocity_vs_superspreader_fraction"])
        plot_sensitivity(sensitivity_rows, plots_dir / OUTPUT_PLOT_FILES["sensitivity_superspreader_fraction_case_fraction"])

    if "epidemic-curves" in selected_experiments:
        curves_rows = require_rows("epidemic curves", curves_rows)
        plot_epidemic_curves(curves_rows, plots_dir / OUTPUT_PLOT_FILES["epidemic_curves"])

    if "routes" in selected_experiments:
        secondary_rows = require_rows("secondary distribution", secondary_rows)
        route_results = route_results or []
        plot_secondary_distribution(secondary_rows, plots_dir / OUTPUT_PLOT_FILES["secondary_distribution"])
        plot_secondary_distribution_normal(secondary_rows, plots_dir / OUTPUT_PLOT_FILES["secondary_distribution_normal"])
        plot_secondary_distribution_superspreaders(
            secondary_rows,
            plots_dir / OUTPUT_PLOT_FILES["secondary_distribution_superspreaders"],
        )
        plot_infection_route(route_results, plots_dir / OUTPUT_PLOT_FILES["infection_route_normal"], "normal")
        plot_infection_route(route_results, plots_dir / OUTPUT_PLOT_FILES["infection_route_strong"], "strong")
        plot_infection_route(route_results, plots_dir / OUTPUT_PLOT_FILES["infection_route_hub"], "hub")

    if "sars" in selected_experiments:
        sars_secondary_rows = require_rows("SARS secondary cases", sars_secondary_rows)
        sars_epicurve_rows = require_rows("SARS epidemic curve", sars_epicurve_rows)
        sars_model_curve_rows = require_rows("SARS model curves", sars_model_curve_rows)
        plot_sars_secondary_cases(
            sars_secondary_rows,
            secondary_rows or [],
            plots_dir / OUTPUT_PLOT_FILES["sars_secondary_patients"],
        )
        plot_sars_case_curve_comparison(
            sars_epicurve_rows,
            sars_model_curve_rows,
            plots_dir / OUTPUT_PLOT_FILES["sars_epidemic_curve_comparison"],
        )


def write_all_plots(
    plots_dir: Path,
    *,
    infection_rows: list[dict[str, object]],
    percolation_rows: list[dict[str, object]],
    critical_rows: list[dict[str, object]],
    critical_reference_rows: list[dict[str, object]],
    front_rows: list[dict[str, object]],
    sensitivity_rows: list[dict[str, object]],
    curves_rows: list[dict[str, object]],
    secondary_rows: list[dict[str, object]],
    sars_secondary_rows: list[dict[str, object]],
    sars_epicurve_rows: list[dict[str, object]],
    sars_model_curve_rows: list[dict[str, object]],
    route_results: list[tuple[SimulationConfig, SimulationResult]],
) -> None:
    write_selected_plots(
        plots_dir,
        selected_experiments={"percolation", "sensitivity", "epidemic-curves", "routes", "sars"},
        infection_rows=infection_rows,
        percolation_rows=percolation_rows,
        critical_rows=critical_rows,
        critical_reference_rows=critical_reference_rows,
        front_rows=front_rows,
        sensitivity_rows=sensitivity_rows,
        curves_rows=curves_rows,
        secondary_rows=secondary_rows,
        sars_secondary_rows=sars_secondary_rows,
        sars_epicurve_rows=sars_epicurve_rows,
        sars_model_curve_rows=sars_model_curve_rows,
        route_results=route_results,
    )


def regenerate_plots_from_outputs(output_dir: Path, selected_experiments: set[str] | None = None) -> None:
    selected_experiments = selected_experiments or {
        "percolation",
        "sensitivity",
        "epidemic-curves",
        "routes",
        "sars",
    }
    route_results: list[tuple[SimulationConfig, SimulationResult]] = []
    if "routes" in selected_experiments:
        selection_rows = read_csv(output_dir / "route_plot_selections.csv")
        for config in route_configs_from_saved_rows(selection_rows):
            result, _metrics = run_simulation(config)
            route_results.append((config, result))

    write_selected_plots(
        output_dir / "plots",
        selected_experiments=selected_experiments,
        infection_rows=(
            read_csv(output_dir / "infection_probability_functions.csv")
            if "percolation" in selected_experiments
            else None
        ),
        percolation_rows=(
            read_csv(output_dir / "percolation_probability.csv")
            if "percolation" in selected_experiments
            else None
        ),
        critical_rows=(
            read_csv(output_dir / "critical_density.csv")
            if "percolation" in selected_experiments
            else None
        ),
        critical_reference_rows=(
            read_csv(output_dir / "critical_density_reference_curves.csv")
            if "percolation" in selected_experiments
            else None
        ),
        front_rows=read_csv(output_dir / "front_distance.csv") if "sensitivity" in selected_experiments else None,
        sensitivity_rows=(
            read_csv(output_dir / "sensitivity_summary.csv")
            if "sensitivity" in selected_experiments
            else None
        ),
        curves_rows=read_csv(output_dir / "epidemic_curves.csv") if "epidemic-curves" in selected_experiments else None,
        secondary_rows=read_csv(output_dir / "secondary_distribution.csv") if "routes" in selected_experiments else None,
        sars_secondary_rows=(
            read_csv(output_dir / "sars_singapore_secondary_patients.csv")
            if "sars" in selected_experiments
            else None
        ),
        sars_epicurve_rows=(
            read_csv(output_dir / "sars_singapore_epidemic_curve.csv")
            if "sars" in selected_experiments
            else None
        ),
        sars_model_curve_rows=(
            read_csv(output_dir / "sars_epidemic_model_curves.csv")
            if "sars" in selected_experiments
            else None
        ),
        route_results=route_results,
    )
