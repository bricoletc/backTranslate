import sys
from Bio import SeqIO
from Bio.Data import CodonTable
from support_objects import DNASample, DistMatrix

try:
    seqfile = sys.argv[1]
    min_dist = float(sys.argv[2])

except Exception:
    print("Provide fasta file of sequences and then min distance as positional args")
    exit(1)

# standard_table = CodonTable.unambiguous_dna_by_name["Standard"].forward_table

translations = set()
records = []
for record in SeqIO.parse(sys.argv[1], "fasta"):
    translations.add(record.seq.translate().__str__())
    records.append(DNASample(record.seq))

assert len(translations) == 1, "More than one protein....."
print("......Checking all dna sequences produce the same protein.......OK")

dm = DistMatrix(records, "hamming", min_dist)
for i in range(dm.dimension):
    for j in range(i):
        assert (
            dm.matrix[i][j] > min_dist
        ), f"Entry {dm.matrix[i][j]} is below minimum distance {min_dist}"
print(
    f".............Checking pairwise distances all above {min_dist}.................OK"
)

print(f"\n \nProduced protein: {translations}")
