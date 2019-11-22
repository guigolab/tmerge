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
    transcripts_copy = transcripts[0:]
    contigs = [ ]

    contig_i = 0
    while transcripts_copy:
        new_contig = build_contig(transcripts_copy)
        for t in new_contig.transcripts:
            transcripts_copy.remove(t)

        contigs.append(new_contig)
        contig_i += 1

    return contigs

def build_contig(transcripts, contig = None, i=0):
    """
    Recursively builds a contig from a given list of transcripts.

    Goes through the transcripts list sequentially and attempts to add it to the contig. 
    Recursion halts when contig does not except the transcript due to IndexError (i.e. no overlap and on same strand)
    """
    if not contig:
        contig = Contig( [transcripts[i] ])
        i+=1
    try:
        contig.add_transcript(transcripts[i])
        return build_contig(transcripts, contig, i+1)
    except TypeError:
        # If transcript[i] is not on same strand then the next transcript may be on same strand and overlap so try
        return build_contig(transcripts, contig, i+1)
    except IndexError:
        # If transcript[i] does not overlap then no overlap and on same strand so return
        return contig