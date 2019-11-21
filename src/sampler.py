from Bio.Data import CodonTable
from collections import defaultdict
from typing import List, Dict, Tuple
import random


class AABackTranslationTable:
    translation_table = defaultdict(list)

    def __init__(self):
        AABackTranslationTable._init_table()

    @staticmethod
    def _init_table():
        if len(AABackTranslationTable.translation_table) != 0:
            return
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
        for codon, aa in standard_table.items():
            AABackTranslationTable.translation_table[aa].append(codon)
        AABackTranslationTable.translation_table = {
            key: sorted(value)
            for key, value in AABackTranslationTable.translation_table.items()
        }

    @staticmethod
    def _reset_table():
        AABackTranslationTable.translation_table = defaultdict(list)


class DNASample:
    def __init__(self, sequence=None, ID=None):
        self.sequence = sequence
        self.id = ID

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __repr__(self):
        return self.sequence

    def __lt__(self, other):
        return self.sequence < other.sequence


class AltCodons:
    def __init__(self, aa):
        assert (
            len(AABackTranslationTable.translation_table) != 0
        ), "ERROR: back translation table uninitialised"

        self.amino_acid = aa.upper()
        assert (
            self.amino_acid in AABackTranslationTable.translation_table
        ), f"ERROR: amino acid {self.amino_acid} does not exist"
        self.choices = AABackTranslationTable.translation_table[self.amino_acid]

    def sample(self):
        return random.choice(self.choices)


class Sampler:
    def __init__(self, AA_seq: str, forbidden: List[str] = []):
        AABackTranslationTable._init_table()
        self.sequence = AA_seq
        self.choices = [AltCodons(aa) for aa in self.sequence]
        self.forbidden = forbidden

    def sample(self, num_samples):
        samples = set()
        while num_samples > 0:
            invalid = False
            num_samples -= 1
            sample = self._one_sample()
            for _no_go in self.forbidden:
                if _no_go in sample.sequence:
                    invalid = True
                    break
            if invalid is False:
                samples.add(sample)

        return list(samples)

    def _one_sample(self):
        dna_seq = ""
        for alt_codons in self.choices:
            dna_seq += alt_codons.sample()
        return DNASample(dna_seq)
