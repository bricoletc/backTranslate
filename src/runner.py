from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from typing import List, Dict, Tuple
from collections import defaultdict
from networkx.algorithms.flow import shortest_augmenting_path
import os, logging
import time, csv

from src.sampler import Sampler, DNASample
from src.distance import DistMatrix
from src.graph import GraphPurifier


class DiffSeqFinder:
    def __init__(
        self,
        min_distance: float,
        output_path=None,
        forbidden=[],
        distance_measure="hamming",
    ):
        self.min_dist = min_distance
        self.output_path = output_path
        self.forbidden = forbidden
        self.samples = []

        self.stats = {
            "min_dist": self.min_dist,
        }
        self.distance_measure = distance_measure

    def _run_DiffFinder(self, samples: List[DNASample]):
        logging.info(f"Computing pairwise Hamming distances")
        t1 = time.clock()
        dm = DistMatrix(samples, self.distance_measure, self.min_dist)
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

    def _write_results(self, mode: str):
        records = []
        for index, item in enumerate(self.samples):
            description = f"{self.description}_{index}"
            rec = SeqRecord(Seq(item.sequence), id=item.id, description=description)
            records.append(rec)

        seqpath = self.output_path + ".fasta"
        SeqIO.write(records, seqpath, "fasta")
        self.stats["sequence_filepath"] = os.path.abspath(seqpath)

        with open(self.output_path + ".tsv", "w") as stats_out:
            fieldnames = self.current_stats_header
            for fn in fieldnames:
                if fn not in self.stats:
                    self.stats.__setitem__(fn, "NA")
            writer = csv.DictWriter(stats_out, fieldnames=fieldnames, delimiter="\t")
            if mode == "dna":
                writer.writeheader()
            writer.writerow(self.stats)


class fromProtein(DiffSeqFinder):
    current_stats_header = [
        "target_num_samples",
        "min_dist",
        "seq_ID",
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
        forbidden=[],
    ):
        super().__init__(min_distance, output_path, forbidden)
        self.seq_ID = seq_ID
        self.aa_seq = aa_sequence
        self.target_num = target_num_samples
        self.stats["target_num_samples"] = self.target_num
        self.stats["seq_ID"] = self.seq_ID
        self.description = "backtranslation"

        self._get_DNA_samples()
        super()._run_DiffFinder(self.samples)

        if self.output_path is not None:
            super()._write_results(mode="notDNA")
        else:
            logging.info(self.stats)

    def _get_DNA_samples(self):
        sampler = Sampler(self.aa_seq, self.seq_ID, forbidden=self.forbidden)
        logging.info(f"Sampling DNA sequences compatible with {self.aa_seq}")
        self.samples = sampler.sample(self.target_num)


class fromDNA(DiffSeqFinder):
    current_stats_header = [
        "min_dist",
        "DistMatrix_t",
        "GraphPuri_t",
        "num_sampled",
        "num_unique_sampled",
        "initial_perID",
        "sampled_perID",
        "sequence_filepath",
    ]

    def __init__(
        self, input_file: str, min_distance: float, output_path=None, forbidden=[],
    ):
        self.input_file = input_file
        self.description = "retained_seq"
        super().__init__(min_distance, output_path, forbidden)

        self._get_DNA_samples()
        self.initial_samples = self._count_uniques()
        super()._run_DiffFinder(self.samples)
        self.sampled_samples = self._count_uniques()

        self._compute_stats(self.initial_samples, self.sampled_samples)

        if self.output_path is not None:
            super()._write_results(mode="dna")
        else:
            logging.info(self.stats)

    def _get_DNA_samples(self):
        num_records = 0
        for seq_record in SeqIO.parse(self.input_file, "fasta"):
            num_records += 1
            if seq_record.seq.alphabet != IUPAC.unambiguous_dna:
                pass
                # print(seq_record.seq.alphabet)
                # logging.error(f"Sequence {str(seq_record.seq)} is not unambiguous DNA")
                # exit(1)
            self.samples.append(DNASample(str(seq_record.seq), ID=seq_record.id))

        if num_records == 0:
            logging.warning(f"no fasta records were parsed from {args.input_file}")

    def _count_uniques(self):
        unique_counts = defaultdict(int)
        for sample in self.samples:
            unique_counts[sample.id] += 1
        return unique_counts

    def _compute_stats(self, pre_sampling, post_sampling):
        total_pre, total_post = sum(post_sampling.values()), sum(pre_sampling.values())
        unique_pre, unique_post = len(post_sampling.keys()), len(pre_sampling.keys())
        breakdown_pre, breakdown_post = "", ""
        for key in pre_sampling.keys():
            breakdown_pre += str(pre_sampling[key]) + "/"
            if key in post_sampling:
                breakdown_post += str(post_sampling[key]) + "/"
            else:
                breakdown_post += "0/"
        self.stats.update(
            {
                "num_sampled": f"{total_pre}/{total_post}",
                "num_unique_sampled": f"{unique_pre}/{unique_post}",
                "initial_perID": breakdown_pre,
                "sampled_perID": breakdown_post,
            }
        )
