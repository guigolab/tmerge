from models.transcript import Transcript

def build(exons):
    """
    Build a list of transcripts.
    """
    transcripts = [
        Transcript([exons[0]])
    ]

    cur_transcript_i = 0
    for exon in exons:
        try:
            transcripts[cur_transcript_i].add_exon(exon)
        except TypeError:
            transcripts.append(Transcript([exon]))
            cur_transcript_i += 1

    return transcripts