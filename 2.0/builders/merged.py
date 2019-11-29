from models.contig import Contig
from models.transcript import Transcript
from functools import reduce

def build(contigs):
    return list(map(merge, contigs))

def merge(contig):
    transcripts = reduce(merge_monoexonic, contig.transcripts, [])
    return Contig(transcripts)

def merge_monoexonic(transcripts, to_merge):
    if len(transcripts) == 0:
        return [to_merge]

    last_transcript = transcripts[len(transcripts) -1]
    if last_transcript.end > to_merge.start:
        # Overlapping so merge
        last_transcript.add_exon(to_merge.exons[0])
        return transcripts

    transcripts.append(to_merge)
    return transcripts