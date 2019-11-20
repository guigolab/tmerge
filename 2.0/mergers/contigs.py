from models.contig import Contig

def build(alignments):
    """
    Build a list of contigs from a list of alignments.

    Contigs are overlapping alignments. 

    Parameters
    ----------
    alignments: list

    Returns
    -------
    list
        list of contigs
    """
    contigs = [
        Contig([alignments[0]], alignments[0].start, alignments[0].end)
    ]

    cur_contig_i = 0
    for alignment in alignments:
        if contigs[cur_contig_i].overlaps(alignment):
            contigs[cur_contig_i].add_alignment(alignment)
        else:
            contigs.append(Contig([alignment], alignment.start, alignment.end))
            cur_contig_i += 1

    return contigs