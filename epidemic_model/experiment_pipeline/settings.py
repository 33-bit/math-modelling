"""Shared constants for the experiment pipeline."""

import numpy as np

from epidemic_model.config import DEFAULT_N, DEFAULT_R0

DEFAULT_SEED_INDICES = tuple(range(1000))
TOP_REACH_MARGIN = DEFAULT_R0
SUPERSPREADER_FRACTION_SWEEP = (0.00, 0.20, 0.40, 0.60, 0.80, 1.00)
DENSITY_SWEEP_POPULATIONS = (50, 75, 100, 125, *range(150, 901, 50))
PLOT_RESOLUTION_DPI = 300
STANDARD_PLOT_SIZE = (4.8, 3.35)
ROUTE_PLOT_SIZE = (3.8, 3.8)
PLOT_LINE_WIDTH = 1.0
PLOT_MARKER_SIZE = 4.2
TRANSMISSION_MODEL_COLORS = {
    "normal": "#6B7280",
    "strong": "#2563EB",
    "hub": "#E11D48",
}
INFECTION_SOURCE_COLORS = {
    "normal": "#6B7280",
    "superspreader": "#E11D48",
}
SUPERSPREADER_FRACTION_COLORS = {
    0.0: "#6B7280",
    0.2: "#2563EB",
    0.4: "#059669",
    0.6: "#D97706",
    0.8: "#7C3AED",
    1.0: "#DC2626",
}
ROUTE_COLORS = {
    "edge": "#94A3B8",
    "susceptible_normal": "#CBD5E1",
    "susceptible_superspreader": "#F59E0B",
    "infected_normal": "#2563EB",
    "infected_superspreader": "#E11D48",
    "initial": "#111827",
}
SARS_BAR_COLOR = "#FDE68A"
SUPERSPREADER_FRACTION_MARKERS = {
    0.0: "o",
    0.2: "s",
    0.4: "^",
    0.6: "D",
    0.8: "v",
    1.0: "*",
}
DOMAIN_SIDE_LENGTH = 10.0 * DEFAULT_R0
NORMALIZED_DENSITY_FACTOR = np.pi * DEFAULT_R0**2
SPREAD_ANALYSIS_NORMALIZED_DENSITY = 20.0
SPREAD_ANALYSIS_POPULATION = int(round(SPREAD_ANALYSIS_NORMALIZED_DENSITY * DOMAIN_SIDE_LENGTH**2 / NORMALIZED_DENSITY_FACTOR))
SPREAD_ANALYSIS_DENSITY = SPREAD_ANALYSIS_POPULATION / DOMAIN_SIDE_LENGTH**2
ROUTE_ANALYSIS_NORMALIZED_DENSITY = 15.0
ROUTE_ANALYSIS_POPULATION = int(round(ROUTE_ANALYSIS_NORMALIZED_DENSITY * DOMAIN_SIDE_LENGTH**2 / NORMALIZED_DENSITY_FACTOR))
ROUTE_ANALYSIS_DENSITY = ROUTE_ANALYSIS_POPULATION / DOMAIN_SIDE_LENGTH**2
SARS_ANALYSIS_NORMALIZED_DENSITY = 15.0
SARS_ANALYSIS_POPULATION = DEFAULT_N
SARS_ANALYSIS_DENSITY = SARS_ANALYSIS_POPULATION / DOMAIN_SIDE_LENGTH**2
DEFAULT_SUPERSPREADER_FRACTION = 0.20
SARS_SUPERSPREADER_FRACTION = 0.40
SARS_TIMESTEP_DAYS = 6
STRONG_CRITICAL_REPRODUCTION_NUMBER = 4.5
HUB_CRITICAL_REPRODUCTION_NUMBER = 3.2
SARS_PROBABLE_CASE_COUNT = 201
OUTPUT_PLOT_FILES = {
    "infection_probability_strong": "fig01_infection_probability_strong.png",
    "infection_probability_hub": "fig02_infection_probability_hub.png",
    "percolation_probability_strong": "fig03_percolation_probability_strong.png",
    "percolation_probability_hub": "fig04_percolation_probability_hub.png",
    "critical_density": "fig05_critical_density.png",
    "front_distance_strong": "fig06_front_distance_strong.png",
    "velocity_vs_superspreader_fraction": "fig07_velocity_vs_superspreader_fraction.png",
    "epidemic_curves": "fig08_epidemic_curves.png",
    "infection_route_strong": "fig09_infection_route_strong.png",
    "infection_route_hub": "fig10_infection_route_hub.png",
    "infection_route_normal": "fig11_infection_route_normal.png",
    "secondary_distribution_normal": "fig12_secondary_distribution_normal.png",
    "secondary_distribution_superspreaders": "fig13_secondary_distribution_superspreaders.png",
    "sars_secondary_patients": "fig14_sars_secondary_patients.png",
    "sars_epidemic_curve_comparison": "fig15_sars_epidemic_curve_comparison.png",
    "percolation_probability": "fig_extra_percolation_probability_comparison.png",
    "front_distance_hub": "fig_extra_front_distance_hub.png",
    "secondary_distribution": "fig_extra_secondary_distribution_comparison.png",
    "sensitivity_superspreader_fraction_case_fraction": "fig_extra_sensitivity_superspreader_fraction_case_fraction.png",
}
SARS_SECONDARY_CASE_FREQUENCIES = (
    (0, 162, "CDC MMWR: 162 probable SARS cases had no evidence of transmission"),
    (1, 19, "Digitized from CDC MMWR Figure 3"),
    (2, 8, "Digitized from CDC MMWR Figure 3"),
    (3, 6, "Digitized from CDC MMWR Figure 3"),
    (7, 1, "Digitized from CDC MMWR Figure 3"),
    (12, 1, "CDC MMWR Case 5 direct probable SARS links"),
    (21, 1, "CDC MMWR Case 1 direct probable SARS links"),
    (23, 2, "CDC MMWR Cases 2 and 3 direct probable SARS links"),
    (40, 1, "CDC MMWR Case 4 direct probable SARS links"),
)
SARS_SIX_DAY_CASE_COUNTS = (
    0,
    0,
    1,
    2,
    17,
    42,
    36,
    20,
    44,
    29,
    20,
    15,
    8,
    4,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
)


