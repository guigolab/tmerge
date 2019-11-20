from importers import gtf
from mergers.contigs import build

data = gtf.parse("sample_data/chr1_short.gff")

contigs = build(data)

print(contigs[2])