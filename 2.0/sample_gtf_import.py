from importers import gtf
from mergers.contigs import build  as build_contigs
from mergers.transcripts import build as build_transcripts
from output.gtf import write

data = gtf.parse("sample_data/chr1_short.gff")

transcripts = build_transcripts(data)
contigs = build_contigs(transcripts)

transcript_ids = [t.id for t in transcripts]
print(set([t for t in transcript_ids if transcript_ids.count(t) > 1]))

print(f"Number of transcripts {len(transcripts)}")
print(f"Number of contigs {len(contigs)}")
start_stops = [(t.start, t.end) for t in transcripts]
print(f"Transcript start and stops: {start_stops[0:5]}...")

write(contigs, "sample_data/out.gtf")