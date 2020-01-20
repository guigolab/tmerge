from models.contig import Contig

def build(transcripts):
    """
    Build a list of contigs from a list of transcripts.

    Contigs are overlapping transcripts.
    
    Warning! Input transcript data is mutated

    Parameters
    ----------
    transcripts: list

    Returns
    -------
    Generator that yields one contig at a time
    """
    while transcripts:
        new_contig = build_contig(transcripts)
        for t in new_contig.transcripts:
            transcripts.remove(t)

        yield new_contig

def build_contig(transcripts):
    """
    Builds a contig from a given list of transcripts.

    Goes through the transcripts list sequentially and attempts to add it to the contig. 
    Halts when contig does not accept the transcript due to IndexError (i.e. no overlap and on same strand)
    """
    contig = Contig([transcripts[0]])
    for transcript in transcripts[1:]:
        try:
            contig.add_transcript(transcript)
        except TypeError:
            # If transcript[i] is not on same strand then the next transcript may be on same strand and overlap so try
            continue
        except IndexError:
            # If transcript[i] does not overlap then no overlap and on same strand so return
            break
    return contig