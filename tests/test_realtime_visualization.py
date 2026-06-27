import unittest

import numpy as np

from realtime_visualization import SimulationSettings, validate_settings, wrapped_route_segments


class RealtimeVisualizationTests(unittest.TestCase):
    def test_wrapped_route_segments_keep_direct_edges_inside_half_box(self) -> None:
        positions = np.array([[1.0, 1.0], [3.0, 2.0]])
        infection_source = np.array([-1, 0])

        segments = wrapped_route_segments(positions, infection_source, L=10.0)

        self.assertEqual(segments, [[(1.0, 1.0), (3.0, 2.0)]])

    def test_wrapped_route_segments_split_edges_across_horizontal_boundary(self) -> None:
        positions = np.array([[9.0, 2.0], [1.0, 4.0]])
        infection_source = np.array([-1, 0])

        segments = wrapped_route_segments(positions, infection_source, L=10.0)

        self.assertEqual(
            segments,
            [[(9.0, 2.0), (10.0, 3.0)], [(0.0, 3.0), (1.0, 4.0)]],
        )

    def test_validate_settings_rejects_bad_lambda(self) -> None:
        with self.assertRaisesRegex(ValueError, "lambda"):
            validate_settings(SimulationSettings(lambda_ss=1.1))


if __name__ == "__main__":
    unittest.main()
