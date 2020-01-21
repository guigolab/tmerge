from models.contig import Contig
from models.transcript import Transcript
from models.exon import Exon
from functools import reduce
from rules.rules import ruleset
from itertools import product

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
            transcripts[i] = merge(transcripts[i], transcripts[i_compare])
            transcripts.pop(i_compare)
        else:
            i_compare += 1

    return contig

def merge(t1, t2):
    shorter, longer = (t1, t2) if len(t1.exons) < len(t2.exons) else (t2, t1)

    # Create a new transcript from the coordinates of the input transcripts
    # Takes the longest overlapping exon from each transcript
    # I.e.
    # ======-----=====
    #    ===-----=========
    # becomes
    # ======-----=========
    exons = []
    for e1 in longer.exons:
        start = e1.start
        end = e1.end
        for e2 in shorter.exons:
            if e1.end > e2.start and e2.end > e1.start:
                start = start if start < e2.start else e2.start
                end = end if end > e2.end else e2.end

        exons.append(Exon(source_type="MERGED", chromosome=t1.chromosome, start=start, end=end, transcript_id=t1.id, gene_id=t1.id, strand=t1.strand))

    return Transcript(exons)