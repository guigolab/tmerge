from importers import gtf
from builders.contigs import build as build_contigs
from builders.transcripts import build as build_transcripts
from output.gtf import write

data = gtf.parse("sample_data/chr1_short.gff")

transcripts = build_transcripts(data)
print(f"Number of transcripts {len(transcripts)}")
contigs = build_contigs(transcripts)
print(f"Number of contigs {len(contigs)}")


def same_intron_set():
    # Not sure about this one
    pass

# If any of these rules return false, merging of the two transcripts fails
def exonic_overlap(transcript1, transcript2):
    # TODO: Ensure order of transcript1 and 2
    return transcript1.end >= transcript2.start - 1

def no_exon_intron_overlap(transcript1, transcript2):
    pass



# write(contigs, "sample_data/out.gtf")