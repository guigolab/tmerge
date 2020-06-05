from functools import reduce

class ReadSupport():
    """
    Compute read support for each transcript and filter by min_read_support.

    Read support is the number of reads in the input (as defined by the exon/intron structure) that match a gien transcript.

    Parameters
    ----------
    hooks: dict
    end_fuzz: int
        Tolerated fuzziness of 5' and 3' ends for two reads to be considered equivalent when computing read support
    min_read_support: int
        Minimum number of reads that must support a given transcript. Any transcripts with read support below this threshold will be removed.
    speed: bool
        When True, favour speed over sensitivity & precision
    """
    def __init__(self, hooks, end_fuzz = 0, min_read_support = 1, speed = False, **kwargs):
        self.end_fuzz = end_fuzz
        self.min_read_support = min_read_support
        self.speed = speed

        hooks["transcript_added"].tap(self.transcript_added)
        hooks["contig_built"].tap(self.calc_support)

    def add_full_length_count(self, transcript):
        if "full_length_count" not in transcript.meta:
            transcript.meta["full_length_count"] = 1

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

    def transcript_added(self, transcript):
        self.add_full_length_count(transcript)          

    def calc_support(self, transcripts):
        if self.speed and self.min_read_support == 1:
            # Don't bother computing if min read support is 1
            return
        supports = set()
        
        for i, target in enumerate(transcripts):
            if self.speed and target.id in supports:
                # Ignore transcripts that support another transcript. This forgoes precision slightly but with a great speed increase
                if target.meta["full_length_count"] < self.min_read_support:
                    target.remove()
                continue
            
            for other in transcripts[i+1:]:
                if other.TSS - self.end_fuzz > target.TSS:
                    # Since input is ordered, can break early if other comes after target
                    break
                
                if self.supports_target(target, other):
                    target.meta["full_length_count"] += 1
                    supports.add(other.id)

                if self.speed and target.meta["full_length_count"] >= self.min_read_support:
                    # If full length count already exceeds min read support, can break early
                    break
            
            if target.meta["full_length_count"] < self.min_read_support:
                target.remove()
