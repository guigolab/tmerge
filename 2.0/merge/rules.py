from itertools import product
from functools import reduce, partial, lru_cache
from utils import ranges

def transcript_overlap(t1, t2):
    return ranges.overlaps((t1.TSS, t1.TES), (t2.TSS, t2.TES))

def ordered_subset(t1, t2):
    return ranges.ordered_subset(t1.junctions, t2.junctions) or ranges.ordered_subset(t2.junctions, t1.junctions)

def no_TSS_TES_overlap(t1, t2):
    return (
        not ranges.within_any(t1.TSS, t2.junctions)
        and not ranges.within_any(t1.TES, t2.junctions)
    ) and (
        not ranges.within_any(t2.TSS, t1.junctions)
        and not ranges.within_any(t2.TES, t1.junctions)
    )

def no_first_last_exon_intron_overlap(t1, t2):
    if not t1.junctions or not t2.junctions:
        return True
    return (
        not ranges.overlaps_any((t1.junctions[len(t1.junctions) - 1][1], t1.TES), t2.junctions)
        and not ranges.overlaps_any((t2.junctions[len(t2.junctions) - 1][1], t2.TES), t1.junctions)
        and not ranges.overlaps_any((t1.TSS, t1.junctions[0][0]), t2.junctions)
        and not ranges.overlaps_any((t2.TSS, t2.junctions[0][0]), t1.junctions)
    )

@lru_cache
def ruleset(t1, t2):
    return (
        t1.chromosome == t2.chromosome
        and t1.strand == t2.strand
        and transcript_overlap(t1, t2)
        and ordered_subset(t1, t2)
        and no_TSS_TES_overlap(t1, t2)
        and no_first_last_exon_intron_overlap(t1,t2)
    )