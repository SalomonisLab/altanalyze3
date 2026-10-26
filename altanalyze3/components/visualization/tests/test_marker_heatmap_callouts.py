"""Regression for sparse event callouts on a dense heatmap with GO terms."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from altanalyze3.components.visualization import marker_heatmap_h5ad as heatmap


class MarkerHeatmapCallouts(unittest.TestCase):
    def test_go_terms_keep_sparse_callouts_legible_and_spaced(self):
        features = [f"ENSG{i:011d}:E1.1-E2.1|ENSG{i:011d}:E1.1-E3.1"
                    for i in range(1600)]
        frame = pd.DataFrame(np.tile([-1., -0.5, 0.5, 1.], (1600, 1)),
                             index=features, columns=["a", "b", "c", "d"])
        original = frame.copy(deep=True)
        selected_rows = [10, 11, 12]
        observed = []
        native_callouts = heatmap._draw_row_callouts

        def capture(ax, labels, rows, n_rows, **kwargs):
            native_callouts(ax, labels, rows, n_rows, **kwargs)
            ax.figure.canvas.draw()
            renderer = ax.figure.canvas.get_renderer()
            texts = [t for t in ax.texts if t.get_text() in labels]
            observed[:] = [(t.get_fontsize(), t.get_window_extent(renderer)) for t in texts]

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "events.pdf"
            with patch.object(heatmap, "_draw_row_callouts", side_effect=capture) as spy:
                heatmap._plot_heatmap(
                    frame, str(output), [("aged", 2), ("MDS-pre", 2)],
                    ["aged", "MDS-pre"], ["aged", "aged", "MDS-pre", "MDS-pre"],
                    ["aged"] * 800 + ["MDS-pre"] * 800,
                    go_terms={"aged": [("RNA processing", 0.01)]},
                    label_events=[features[i] for i in selected_rows],
                    label_shared_exons=False, dpi=30, write_svg=False,
                )
            self.assertTrue(spy.called, "GO terms bypassed native leader-line callouts")
            self.assertEqual(len(observed), 3)
            self.assertTrue(all(size >= 3.5 for size, _ in observed))
            for i, (_, first) in enumerate(observed):
                for _, second in observed[i + 1:]:
                    self.assertFalse(first.overlaps(second), "Sparse event labels overlap")
            labels = pd.read_csv(output.with_name("events_label_events.tsv"), sep="\t")
            self.assertEqual(set(labels.heatmap_row), {11, 12, 13})
            self.assertTrue(output.is_file())
        pd.testing.assert_frame_equal(frame, original)


if __name__ == "__main__":
    unittest.main()
