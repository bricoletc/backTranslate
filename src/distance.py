import numpy as np
from typing import List, Dict, Tuple

from sampler import DNASample


class Distance(object):
    @staticmethod
    def factory(type):
        if type == "hamming":
            return Hamming()
        assert 0, f"{type} is not an implemented distance metric"


class Hamming(Distance):
    def __call__(self, seq1: np.array, seq2: np.array):
        assert len(seq1) == len(seq2), "ERROR: sequences must be of same length"
        return (seq1 != seq2).sum()


class DistMatrix:
    def __init__(
        self, samples: List[DNASample], distance_metric: str, min_distance: int = 0
    ):
        self.distance = Distance.factory(distance_metric)
        self.samples = samples
        self._vectorised_samples = self._vectorise(self.samples)
        self.min_dist = min_distance
        self.dimension = len(self.samples)
        self.matrix = np.full((self.dimension, self.dimension), np.inf)
        self._compute_dists()

    def _vectorise(self, samples: List[DNASample]):
        arrayed_strings = []
        for sample in self.samples:
            arrayed_strings.append(np.array([char for char in sample.sequence]))
        return np.array(arrayed_strings)

    def _compute_dists(self):
        for i in range(self.dimension):
            for j in range(i + 1):
                self.matrix[i][j] = self.distance(
                    self._vectorised_samples[i], self._vectorised_samples[j]
                )
        if self.dimension > 0:
            self.matrix /= len(self.samples[0])

    def to_Graph(self):
        graph = nx.Graph()
        graph.add_nodes_from(self.samples)
        for i in range(self.dimension):
            for j in range(i):
                if self.matrix[i][j] < self.min_dist:
                    graph.add_edge(self.samples[i], self.samples[j])
        return graph
