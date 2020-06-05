from itertools import product
from functools import reduce, partial

from ..utils import ranges

def transcript_overlap(t1, t2):
    """
    Check if t1 and t2 overlap
    """
    return ranges.overlaps((t1.TSS, t1.TES), (t2.TSS, t2.TES))

def ordered_subset(t1, t2, tolerance = 0):
    """
    Check if the junctions of t1 are an ordered subset of the junctions of t2 or vice versa.
    """
    return ranges.ordered_subset(t1.sorted_junctions, t2.sorted_junctions, tolerance) or ranges.ordered_subset(t2.sorted_junctions, t1.sorted_junctions, tolerance)

def TSS_TES_overlap(t1, t2, tolerance = 0):
    """
    Check if the TSS (5') of t1 or the TES (3') of t1 is not within any of the junctions of t2 and vice versa.
    """
    return (
        ranges.within_any(t1.TSS + tolerance, t2.junctions, exclusive=False)
        or ranges.within_any(t1.TES - tolerance, t2.junctions, exclusive=False)
    ) or (
        ranges.within_any(t2.TSS + tolerance, t1.junctions, exclusive=False)
        or ranges.within_any(t2.TES - tolerance, t1.junctions, exclusive=False)
    )

def first_last_exon_intron_overlap(t1, t2, tolerance = 0):
    """
    Test whether the first/last exon overlaps with the first/last intron.
    Since junctions consists of splice junction coordinates, must construct the coordinates of the first/last exon from TSS/TES.
    """
    if not t1.junctions or not t2.junctions:
        return False
    
    # Invert tolerance since considering exons this time
    tolerance = -tolerance
    return (
        # First exon
        ranges.overlaps_any((t1.TSS, t1.sorted_junctions[0][0]), t2.sorted_junctions, tolerance)
        or ranges.overlaps_any((t2.TSS, t2.sorted_junctions[0][0]), t1.sorted_junctions, tolerance)
        # Last exon
        or ranges.overlaps_any((t1.sorted_junctions[len(t1.junctions) - 1][1], t1.TES), t2.sorted_junctions, tolerance)
        or ranges.overlaps_any((t2.sorted_junctions[len(t2.junctions) - 1][1], t2.TES), t1.sorted_junctions, tolerance)
    )

def ruleset(t1, t2, tolerance = 0):
    """
    Combines all the neccesary rules to determine if t1 and t2 can be merged.

    Parameters
    ----------
    t1: TranscriptModel
    t2: TranscriptModel
    tolerance: int
        Maximum number of nucleotides of terminal exon overhang of t1 allowed within an t2 (and vice versa).
        Setting this to a positive integer can correct mismapped splice junctions that sometimes occur when aligning very short, error-rich terminal read exons.
    """
    if tolerance:
        tolerance = tolerance + 1 # utils.ranges funcs are all exclusive. Want tolerance to be inclusive
    return (
        t1.chromosome == t2.chromosome
        and t1.strand == t2.strand
        and transcript_overlap(t1, t2)
        and not TSS_TES_overlap(t1, t2, tolerance)
        and ordered_subset(t1, t2)
        and not first_last_exon_intron_overlap(t1,t2, tolerance)
    )