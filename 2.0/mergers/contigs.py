from models.contig import Contig

def build(exons):
    """
    Build a list of contigs from a list of exons.

    Contigs are overlapping exons. 

    Parameters
    ----------
    exons: list

    Returns
    -------
    list
        list of contigs
    """
    contigs = [
        Contig([exons[0]])
    ]

    cur_contig_i = 0
    for exon in exons:
        try:
            contigs[cur_contig_i].add_exon(exon)
        except IndexError:
            contigs.append(Contig([exon]))
            cur_contig_i += 1

    return contigs