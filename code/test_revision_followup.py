"""Check the intended weighting, rather than another copy of its implementation."""
import unittest
import numpy as np
from scipy.spatial.distance import cdist
from revision_followup import block_weights, weighted_jaccard


class BlockWeightTests(unittest.TestCase):
    def test_uniform_weights_reduce_to_jaccard(self):
        target = np.array([1, 0, 1, 0], dtype=bool)
        anchors = np.array([[0, 1, 1, 0], [1, 0, 0, 0]], dtype=bool)
        np.testing.assert_allclose(weighted_jaccard(target, anchors, np.ones(4)),
                                   cdist([target], anchors, metric="jaccard")[0])

    def test_duplicating_entire_block_preserves_distance(self):
        target = np.array([1, 0, 1])
        anchors = np.array([[0, 1, 1], [1, 0, 0]])
        blocks = np.array(["a", "a", "b"])
        expected = weighted_jaccard(target, anchors, block_weights(blocks))
        # Duplicate every column in block a, not just one chosen diagnostic.
        indices = [0, 1, 0, 1, 2]
        actual = weighted_jaccard(target[indices], anchors[:, indices], block_weights(blocks[indices]))
        np.testing.assert_allclose(actual, expected)
        for label in ["a", "b"]:
            self.assertAlmostEqual(block_weights(blocks[indices])[blocks[indices] == label].sum(), 1)

    def test_weighted_union_is_not_mean_of_block_distances(self):
        # Block a: mismatch 1/2, union 1/2; block b: mismatch 0, union 1.
        # Ratio of weighted sums is 1/3, whereas mean block distance is 1/2.
        weights = block_weights(["a", "a", "b"])
        self.assertAlmostEqual(weighted_jaccard([1, 0, 1], [[0, 0, 1]], weights)[0], 1/3)
        self.assertEqual(weighted_jaccard([0, 0, 0], [[0, 0, 0]], weights)[0], 0)


if __name__ == "__main__":
    unittest.main()
