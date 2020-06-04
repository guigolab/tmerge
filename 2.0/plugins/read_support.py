from utils import ranges
from functools import reduce
from itertools import permutations, islice

class ReadSupport():
    def __init__(self, merger, end_fuzz = 0, min_read_support = 1, speed = False):
        self.end_fuzz = end_fuzz
        self.min_read_support = min_read_support
        self.speed = speed

        merger.hooks["transcript_added"].tap(self.transcript_added)
        merger.hooks["contig_built"].tap(self.calc_support)

    def add_full_length_count(self, transcript):
        if "full_length_count" not in transcript.meta:
            transcript.meta["full_length_count"] = 1

    def get_cum_intron_length(self, tm):
        return reduce(lambda x,y: x+y, map(lambda x: x[1] - x[0], tm), 0)

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
        
        # return
        # # Get the junctions from target that are before the TSS of other
        # TSS_junctions = [j for j in target.junctions if j[1] <= other.TSS]
        # TSS_intron_length = self.get_cum_intron_length(TSS_junctions)
        # # Calculate the spliced distance between target TSS and other TSS
        # dist_TSS = other.TSS - target.TSS - TSS_intron_length
        
        # # Get the junctions from target that are after the TES of other
        # TES_junctions = [j for j in target.junctions if j[0] >= other.TES]
        # TES_intron_length = self.get_cum_intron_length(TES_junctions)
        # # Calculate the spliced distance between other TES and target TES
        # dist_TES = target.TES - other.TES - TES_intron_length
        
        # if dist_TSS < self.end_fuzz and dist_TES < self.end_fuzz:
        #      target.meta["full_length_count"] += 1

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
