from importers import gtf
from mergers.contigs import build as build_contigs
from mergers.transcripts import build as build_transcripts
from output.gtf import write

data = gtf.parse("sample_data/chr1_short.gff")

transcripts = build_transcripts(data)
print(f"Number of transcripts {len(transcripts)}")
contigs = build_contigs(transcripts)

print(f"Number of contigs {len(contigs)}")

write(contigs, "sample_data/out.gtf")