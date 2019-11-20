from models.contig import Contig

def build(transcripts):
    """
    Build a list of contigs from a list of transcripts.

    Contigs are overlapping transcripts. 

    Parameters
    ----------
    transcripts: list

    Returns
    -------
    list
        list of contigs
    """
    contigs = [
        Contig([transcripts[0]])
    ]

    cur_contig_i = 0
    for transcript in transcripts:
        try:
            contigs[cur_contig_i].add_transcript(transcript)
        except IndexError:
            contigs.append(Contig([transcript]))
            cur_contig_i += 1

    return contigs