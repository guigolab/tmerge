from functools import reduce

class ReadSupport():
    """
    Computes read support for each transcript and filters transcripts that do not meet min_isoform_fraction or min_read_support.

    isoform_fraction: read support as a fraction of the maximum read support in a contig.
    Read support: number of reads in the input (as defined by the exon/intron structure) that match a gien transcript.

    Parameters
    ----------
    hooks: dict
    end_fuzz: int
        Tolerated fuzziness of 5' and 3' ends for two reads to be considered equivalent when computing read support
    min_isoform_fraction: float
        Minimum isoform fraction for a transcript. Any transcripts with isoform_fraction below this threshold will be removed.
    min_read_support: int
        Minimum read support for a transcript
    """
    def __init__(self, hooks, end_fuzz = 0, min_isoform_fraction = 0, min_read_support = 1, no_merge = False, **kwargs):
        if min_isoform_fraction > 1:
            raise TypeError("min_isoform_fraction must be < 1.")
        if min_read_support < 1:
            raise TypeError("min_read_support must be >= 1")
        if min_isoform_fraction > 0 and min_read_support > 1:
            raise TypeError("You cannot use min_isoform_fraction and min_read_support at the same time.")

        self.end_fuzz = end_fuzz
        self.min_isoform_fraction = min_isoform_fraction
        self.min_read_support = min_read_support
        self.no_merge = no_merge

        hooks["transcript_added"].tap(self.add_meta)
        hooks["contig_built"].tap(self.calc_support)

    def add_meta(self, transcript):
        transcript.meta["read_support"] = 1
        transcript.meta["isoform_fraction"] = 0

    def supports_target(self, target, other):
        if (
            other.junctions == target.junctions and 
            other.TSS >= target.TSS - self.end_fuzz and
            other.TSS <= target.TSS + self.end_fuzz and
            other.TES >= target.TES - self.end_fuzz and
            other.TES <= target.TES + self.end_fuzz
        ):
            return True
        
        return False 

    def calc_support(self, transcripts):
        max_read_support = 1
        for i, target in enumerate(transcripts):
            if  target.removed:
                continue
            
            for other in transcripts[i+1:]:
                if other.TSS - self.end_fuzz > target.TSS:
                    # Since input is ordered, can break early if other comes after target
                    break
                
                if other.removed:
                    continue
                
                if self.supports_target(target, other):
                    target.meta["read_support"] += 1
                    if not self.no_merge:
                        # Remove if supports target as is the same TM and would be merged anyway
                        # Improves speed as reduces search space
                        # Don't run in no_merge mode as we want to keep all transcripts that would be merged.
                        other.remove()

            if target.meta["read_support"] > max_read_support:
                max_read_support = target.meta["read_support"]

            if target.meta["read_support"] < self.min_read_support:
                target.remove()

        for t in transcripts:
            t.meta["isoform_fraction"] = t.meta["read_support"] / max_read_support
            if t.meta["isoform_fraction"] < self.min_isoform_fraction:
                t.remove()
