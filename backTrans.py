from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
from typing import List, Dict, Tuple
import networkx as nx
from networkx.algorithms.flow import shortest_augmenting_path
import numpy as np
import random, logging, argparse, os
import time, csv


def CLI_parser():
    parser = argparse.ArgumentParser(
        description="Backtranslate amino acid sequences with pairwise distance constraints"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-f",
        "--file",
        dest="input_file",
        help="Path to fasta file containing amino acid sequences",
        type=str,
    )
    group.add_argument(
        "-s",
        "--sequence",
        dest="sequence",
        help="amino acid sequence passed on command-line",
        type=str,
    )

    parser.add_argument(
        "-n",
        "--num_samples",
        dest="num_samples",
        help="Maximum number of backtranslated DNA samples to produce",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "-d",
        "--min_distance",
        dest="min_dist",
        help="Minimum distance (Hamming, as fraction of distinct nucleotides) between all returned sequences",
        type=float,
        default=0.3,
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output",
        help="An existing directory where the output will go. This is a fasta of the sampled sequence(s) and a stats file",
        type=str,
    )

    parser.add_argument(
        "--output_prefix",
        dest="out_prefix",
        help="prefix for output files",
        type=str,
        default="results",
    )

    parser.add_argument(
        "--stats_header",
        dest="stats",
        action="store_true",
        help="Prints the header for the stats file",
    )

    return parser


def validate_args(args, parser):

    if args.input_file is not None and not os.path.exists(args.input_file):
        logging.error(f"file {args.input_file} does not exist")
        parser.print_help()
        exit(1)

    if args.output is not None and not os.path.exists(args.output):
        logging.error(f"directory {args.output} does not exist")
        parser.print_help()
        exit(1)

    if not 0 <= args.min_dist <= 1:
        logging.error(f"min distance {args.min_dist} is not between 0 and 1")
        parser.print_help()
        exit(1)


class DNASample:
    def __init__(self, sequence=None):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __repr__(self):
        return self.sequence

    def __lt__(self, other):
        return self.sequence < other.sequence


class Distance(object):
    @staticmethod
    def factory(type):
        if type == "hamming":
            return Hamming()
        assert 0, f"{type} is not an implemented distance metric"


class Hamming(Distance):
    def __call__(self, seq1: DNASample, seq2: DNASample):
        # return 0
        assert len(seq1) == len(seq2), "ERROR: sequences must be of same length"
        # num_diffs = 0
        # for i in range(len(seq1)):
        #     if seq1[i] != seq2[i]
        return sum([seq1[i] != seq2[i] for i in range(len(seq1))])


class BackTranslater:
    translation_table = defaultdict(list)

    def __init__(self):
        BackTranslater._init_table()

    @staticmethod
    def _init_table():
        if len(BackTranslater.translation_table) != 0:
            return
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
        for codon, aa in standard_table.items():
            BackTranslater.translation_table[aa].append(codon)
        BackTranslater.translation_table = {
            key: sorted(value)
            for key, value in BackTranslater.translation_table.items()
        }

    @staticmethod
    def _reset_table():
        BackTranslater.translation_table = defaultdict(list)


class AltCodons:
    def __init__(self, aa):
        assert (
            len(BackTranslater.translation_table) != 0
        ), "ERROR: back translation table uninitialised"

        self.amino_acid = aa.upper()
        assert (
            self.amino_acid in BackTranslater.translation_table
        ), f"ERROR: amino acid {self.amino_acid} does not exist"
        self.choices = BackTranslater.translation_table[self.amino_acid]

    def sample(self):
        return random.choice(self.choices)


class Sampler:
    def __init__(self, AA_seq: str):
        BackTranslater._init_table()
        self.sequence = AA_seq
        self.choices = [AltCodons(aa) for aa in self.sequence]

    def sample(self, num_samples):
        samples = set()
        while num_samples > 0:
            num_samples -= 1
            samples.add(self._one_sample())
        return list(samples)

    def _one_sample(self):
        dna_seq = ""
        for alt_codons in self.choices:
            dna_seq += alt_codons.sample()
        return DNASample(dna_seq)


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
            arrayed_strings.append([char for char in sample.sequence])
        return np.array(arrayed_strings)

    def _compute_dists(self):
        for i in range(self.dimension):
            for j in range(i + 1):
                self.matrix[i][j] = (
                    self._vectorised_samples[i] != self._vectorised_samples[j]
                ).sum()

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


class GraphPurifier:
    """
    'Purify' here means resolving the graph down to an edgeless graph
    with the greatest possible number of nodes (DNASamples)

    The _purify functions are ordered by how well they (seem to) perform, which
    is also how long they take to run.
    """

    def __init__(self, graph: nx.Graph):
        self.graph = graph
        initial_ccs = self.get_connected_components(self.graph)
        self._disjoint_samples = None
        # self._purify_degreeCentrality()
        self._purify_maxIndependentSet()
        # self._purify_minCut(initial_ccs)
        assert (
            len(self.graph.edges) == 0
        ), "ERROR: the graph should be edgeless after call to ._purify()"

    @property
    def disjoint_samples(self):
        if self._disjoint_samples is None:
            self._disjoint_samples = list(self.graph.nodes)
        return self._disjoint_samples

    @staticmethod
    def get_connected_components(graph: nx.Graph):
        """
        .subgraph() returns an immutable *view* on input graph
        """
        connected_components = list(nx.connected_components(graph))
        return [graph.subgraph(cc) for cc in connected_components]

    def _purify_degreeCentrality(self):
        degrees = nx.degree_centrality(self.graph)
        while self.graph.number_of_edges() > 0:
            max_node = max(degrees, key=degrees.get)
            degrees.pop(max_node)
            self.graph.remove_node(max_node)

    def _purify_maxIndependentSet(self):
        nodes = nx.maximal_independent_set(self.graph)
        self.graph = self.graph.subgraph(nodes).copy()

    def _purify_minCut(self, connected_components):
        for cc in connected_components:
            if len(cc.nodes) == 1:
                continue
            # spanning_tree = nx.minimum_spanning_tree(cc)
            # min_cut = list(nx.minimum_node_cut(spanning_tree))
            min_cut = list(nx.minimum_node_cut(cc, flow_func=shortest_augmenting_path))
            self.graph.remove_nodes_from(min_cut)
            new_ccs = self.get_connected_components(cc)
            self._purify_minCut(new_ccs)


class BackTranslate:
    current_stats_header = [
        "target_num_samples",
        "min_dist",
        "seqID",
        "DistMatrix_t",
        "GraphPuri_t",
        "num_sampled",
        "sequence_filepath",
    ]

    def __init__(
        self,
        aa_sequence: str,
        target_num_samples: int,
        min_distance: float,
        seq_ID: str,
        output_path=None,
    ):
        self.aa_seq = aa_sequence
        self.target_num = target_num_samples
        self.min_dist = min_distance
        self.output_path = output_path
        self.samples = []
        self.seq_ID = seq_ID

        self.stats = {
            "target_num_samples": self.target_num,
            "min_dist": self.min_dist,
            "seqID": seq_ID,
        }
        self._run()
        if self.output_path is not None:
            self._write_results()
        else:
            logging.info(self.stats)

    def _run(self):
        sampler = Sampler(self.aa_seq)
        logging.info(f"Sampling DNA sequences compatible with {self.aa_seq}")
        samples = sampler.sample(self.target_num)

        logging.info(f"Computing pairwise Hamming distances")
        t1 = time.clock()
        dm = DistMatrix(samples, "hamming", self.min_dist)
        t2 = time.clock()
        self.stats["DistMatrix_t"] = int(t2 - t1)

        logging.info(
            f"Pruning proximity graph to impose min pairwise distance of {self.min_dist} everywhere"
        )
        t1 = time.clock()
        graph = GraphPurifier(dm.to_Graph())
        t2 = time.clock()
        self.stats["GraphPuri_t"] = int(t2 - t1)

        self.samples = graph.disjoint_samples
        self.stats["num_sampled"] = len(self.samples)

    def _write_results(self):
        records = []
        for index, item in enumerate(self.samples):
            rec = SeqRecord(
                Seq(item.sequence),
                id=self.seq_ID,
                description=f"backtranslation_{index}",
            )

            records.append(rec)
        seqpath = self.output_path + ".fasta"
        SeqIO.write(records, seqpath, "fasta")
        self.stats["sequence_filepath"] = os.path.abspath(seqpath)

        with open(self.output_path + ".tsv", "w") as stats_out:
            fieldnames = self.stats.keys()
            writer = csv.DictWriter(stats_out, fieldnames=fieldnames, delimiter="\t")
            writer.writerow(self.stats)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO
    )

    parser = CLI_parser()
    args = parser.parse_args()
    if args.stats:
        header = "\t".join(BackTranslate.current_stats_header)
        logging.info(f"{header}")
        exit(0)
    validate_args(args, parser)
    output_path = (
        os.path.join(args.output, args.out_prefix) if args.output is not None else None
    )

    if args.sequence is not None:
        BackTranslate(
            args.sequence,
            args.num_samples,
            args.min_dist,
            seq_ID="No_ID",
            output_path=output_path,
        )

    else:
        num_records = 0
        this_output_path = output_path
        for seq_record in SeqIO.parse(args.input_file, "fasta"):
            if args.output is not None:
                this_output_path = f"{output_path}_{seq_record.id}"
            num_records += 1
            BackTranslate(
                seq_record.seq,
                args.num_samples,
                args.min_dist,
                seq_ID=seq_record.id,
                output_path=this_output_path,
            )

        if num_records == 0:
            logging.warning(f"no fasta records were parsed from {args.input_file}")
