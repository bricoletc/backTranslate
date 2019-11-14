from support_objects import (
    DNASample,
    BackTranslater,
    AltCodons,
    Sampler,
    DistMatrix,
    Hamming,
    GraphPurifier,
)
import numpy as np
import networkx as nx
import unittest
from unittest.mock import patch, Mock
import random
import itertools

translation_table = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
}
translation_table = {key: sorted(value) for key, value in translation_table.items()}


class TestDNASample(unittest.TestCase):
    def test_EmptyConstructor_noSequence(self):
        sample = DNASample()
        self.assertEqual(sample.sequence, None)


class TestBackTranslater(unittest.TestCase):
    def tearDown(self):
        BackTranslater._reset_table()

    def test_callStaticMethod_populatedTranslationTable(self):
        BackTranslater._init_table()
        self.assertEqual(BackTranslater.translation_table, translation_table)

    def test_initialiseObject_populatedTranslationTable(self):
        bt = BackTranslater()
        self.assertEqual(bt.translation_table, translation_table)


class TestAltCodons(unittest.TestCase):
    def setUp(self):
        BackTranslater._init_table()

    def tearDown(self):
        BackTranslater._reset_table()

    def test_initAltCodonsWithoutBackTranslator_Fails(self):
        BackTranslater._reset_table()
        with self.assertRaises(AssertionError):
            alt_codons = AltCodons("A")

    def test_AltCodonsWithLowerCaseAA_Works(self):
        alt_codons = AltCodons("a")
        self.assertEqual("A", alt_codons.amino_acid)

    def test_AltCodonsUnknownAA_Fails(self):
        with self.assertRaises(AssertionError):
            alt_codons = AltCodons("Z")

    def test_sampleAltCodons_correctCodons(self):
        for aa in translation_table:
            alt_codons = AltCodons(aa)
            self.assertTrue(alt_codons.sample() in translation_table[aa])


class TestSampler(unittest.TestCase):
    def tearDown(self):
        BackTranslater._reset_table()

    def test_initSampler_populatedTranslationTable(self):
        self.assertEqual(len(BackTranslater.translation_table), 0)
        sampler = Sampler("AAAA")
        self.assertEqual(BackTranslater.translation_table, translation_table)

    def test_OneSampleDNASequence_CorrectObjectType(self):
        sampler = Sampler("ACAC")
        sampled = sampler._one_sample()
        self.assertTrue(type(sampled) is DNASample)

    def test_OneSampleDNASequence_CorrectLength(self):
        AA_seq = "ACAC"
        sampler = Sampler(AA_seq)
        sampled = sampler._one_sample()
        self.assertTrue(len(sampled) == 3 * len(AA_seq))

    def test_SampleDNASequences_CorrectUpperBound(self):
        make_seq = lambda keys, seq_size: "".join(
            [random.choice(keys) for time_repeated in range(seq_size)]
        )
        AAs = list(translation_table.keys())
        seq_list = [make_seq(AAs, 10) for seq_num in range(5)]

        for seq in seq_list:
            sampler = Sampler(seq)
            self.assertTrue(1 <= len(sampler.sample(10)) <= 10)

    def test_SampleForbiddenDNASequences_CorrectlyRejected(self):
        """
        Generate all DNA strings that can produce AA_seq,
        and pass those as forbidden.
        Then we expect all random samples to be rejected.
        """
        AA_seq = "ACAC"
        sampler = Sampler(AA_seq)
        all_choices = [alts.choices for alts in sampler.choices]  # list of lists
        product = itertools.product(*all_choices)  # '*': list unpacked
        product = ["".join(combi) for combi in product]
        forbidden = set(product)
        sampler = Sampler(AA_seq, forbidden=forbidden)

        samples = sampler.sample(50)
        self.assertEqual(len(samples), 0)


class TestDistMetrics(unittest.TestCase):
    def test_Hamming_correctVal(self):
        seq1, seq2 = "AGGACC", "AGGAAA"
        metric = Hamming()
        self.assertEqual(metric(seq1, seq2), 2)


class TestDistMatrix_Distances(unittest.TestCase):
    def setUp(self):
        self.DNA_sequences = ["AGGACC", "AGGCCC", "AGGTCC", "TTTTTT"]
        self.DNA_sequences = [DNASample(seq) for seq in self.DNA_sequences]

    def test_unknownDistMetric_Fails(self):
        with self.assertRaises(AssertionError):
            dm = DistMatrix(self.DNA_sequences, "weirdMetric")

    def test_HammingDistMetric_correctlyTyped(self):
        with patch.object(DistMatrix, "_compute_dists", return_value=None):
            dm = DistMatrix(self.DNA_sequences, "hamming")
            self.assertTrue(type(dm.distance) is Hamming)

    def test_HammingDistMetric_correctValues(self):
        dm = DistMatrix(self.DNA_sequences, "hamming")
        self.assertEqual(dm.dimension, 4)
        expected = np.array(
            [
                [0, np.inf, np.inf, np.inf],
                [1 / 6, 0, np.inf, np.inf],
                [1 / 6, 1 / 6, 0, np.inf],
                [1, 1, 5 / 6, 0],
            ]
        )
        self.assertTrue(np.array_equal(dm.matrix, expected))


class TestDistMatrix_vectorisedStrings(unittest.TestCase):
    def test_SequenceSet_correctNpArray(self):
        sequences = [DNASample("ACAG"), DNASample("AGGG")]
        dm = DistMatrix(sequences, "hamming")
        expected = np.array([["A", "C", "A", "G"], ["A", "G", "G", "G"]])
        self.assertTrue((dm._vectorised_samples == expected).all())


class TestDistMatrix_toGraph(unittest.TestCase):
    def setUp(self):
        distMat = np.array([[0, np.inf, np.inf], [1 / 3, 0, np.inf], [1 / 3, 1, 0]])
        dimension = distMat.shape[0]
        samples = [DNASample() for i in range(dimension)]
        self.dm = DistMatrix([], "hamming")
        self.dm.dimension, self.dm.matrix, self.dm.samples = dimension, distMat, samples

    def test_toGraphNoMinDistance_graphHasNoEdges(self):
        """
        The default maximum distance for an edge between nodes is 0;
        here we therefore expect no edges
        """
        induced_graph = self.dm.to_Graph()
        self.assertEqual(len(induced_graph.nodes), self.dm.dimension)
        self.assertEqual(len(induced_graph.edges), 0)

    def test_toGraphwithMinDistance_graphHasEdges(self):
        self.dm.min_dist = 2 / 3
        induced_graph = self.dm.to_Graph()
        self.assertEqual(len(induced_graph.edges), 2)


class TestGraphPurifier(unittest.TestCase):
    def test_ConnectedComponents_CorrectNumberAndType(self):
        graph = nx.Graph()
        graph.add_nodes_from([1, 2, 3, 4])
        graph.add_edges_from([(1, 2), (3, 4)])
        subgraphs = GraphPurifier.get_connected_components(graph)
        self.assertEqual(len(subgraphs), 2)
        self.assertTrue(all([type(s) == nx.classes.graph.Graph for s in subgraphs]))

    def test_TreeProximityGraph_twoRemainingSamples(self):
        samples = [
            DNASample("ACG"),
            DNASample("ATG"),
            DNASample("GCG"),
        ]
        dm = DistMatrix(samples, "hamming", min_distance=2 / 3)
        graph = dm.to_Graph()

        gp = GraphPurifier(graph)
        expected = sorted(samples[1:])
        self.assertTrue(1 <= len(gp.disjoint_samples) <= 2)

    def test_CliqueProximityGraph_singleRemainingSample(self):
        samples = [
            DNASample("ACA"),
            DNASample("ATA"),
            DNASample("AGA"),
        ]
        dm = DistMatrix(samples, "hamming", min_distance=2 / 3)
        graph = dm.to_Graph()
        gp = GraphPurifier(graph)
        self.assertEqual(len(gp.disjoint_samples), 1)

    def test_completeGraph_singleRemainingNode(self):
        graph = nx.complete_graph(10)
        gp = GraphPurifier(graph)
        self.assertEqual(len(gp.disjoint_samples), 1)


if __name__ == "__main__":
    unittest.main()
