from utils import ranges
from functools import reduce

class ReadSupport():
    def __init__(self, merger, end_fuzz = 50, min_read_support = 1):
        self.end_fuzz = end_fuzz
        self.min_read_support = min_read_support
        merger.hooks["contig_built"].tap(self.add_full_length_count)
        merger.hooks["tm_merged"].tap(self.calc_full_length_count)
        merger.hooks["contig_merged"].tap(self.remove_unsupported)

    def add_full_length_count(self, transcripts):
        for t in transcripts:
            if "full_length_count" not in t.meta:
                t.meta["full_length_count"] = 1

    def get_cum_intron_length(self, tm):
        return reduce(lambda x,y: x+y, map(lambda x: x[1] - x[0], tm), 0)

    def calc_full_length_count(self, merge_result, merged):
        # Note: merge_result is always the result of the merge so the TM with lower TSS and higher TES

        # Get the junctions from merge_result that are before the TSS of merged
        TSS_junctions = [j for j in merge_result.junctions if j[1] <= merged.TSS]
        TSS_intron_length = self.get_cum_intron_length(TSS_junctions)
        # Calculate the spliced distance between merged_result TSS and merged TSS
        dist_TSS = merged.TSS - merge_result.TSS - TSS_intron_length
        
        # Get the junctions from merge_result that are after the TES of merged
        TES_junctions = [j for j in merge_result.junctions if j[0] >= merged.TES]
        TES_intron_length = self.get_cum_intron_length(TES_junctions)
        # Calculate the spliced distance between merged TES and merged_result TES
        dist_TES = merge_result.TES - merged.TES - TES_intron_length
        
        # merge_result.meta["full_length_count"] += merged.meta.get("full_length_count", 0)
        if dist_TSS <= self.end_fuzz and dist_TES <= self.end_fuzz:
            merge_result.meta["full_length_count"] += 1

    def remove_unsupported(self, transcripts):
        if self.min_read_support == 1:
            return

        to_remove = []
        for t in transcripts:
            if t.meta["full_length_count"] < self.min_read_support:
                to_remove.append(t)

        for t in to_remove:
            transcripts.remove(t)