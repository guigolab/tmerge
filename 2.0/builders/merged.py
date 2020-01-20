from models.contig import Contig
from models.transcript import Transcript
from functools import reduce
from merge.rules import transcript_overlap, same_introns, no_exon_intron_overlap, monoexonic_overlap, ordered_subset
from merge.merge import merge as m

def build(contig):
    if len(contig.transcripts) == 1:
        return contig
    transcripts = contig.transcripts

    i = 0
    i_compare = 1
    while i < len(transcripts) -1:
        if i_compare >= len(transcripts) - 1:
            i += 1
            i_compare = 0

        if i == i_compare:
            i_compare += 1

        if ruleset(transcripts[i], transcripts[i_compare]):
            transcripts[i] = m(transcripts[i], transcripts[i_compare])
            transcripts.pop(i_compare)
        else:
            i_compare += 1

    return contig
            

def ruleset(t1, t2):
    if not transcript_overlap(t1, t2):
        return False
    if same_introns(t1, t2):
        return True
    if no_exon_intron_overlap(t1, t2) and monoexonic_overlap(t1, t2):
        return True
    if no_exon_intron_overlap(t1, t2) and ordered_subset(t1, t2):
        return True
    return False