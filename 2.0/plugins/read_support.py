from utils import ranges

class ReadSupport():
    def __init__(self, merger, end_fuzz = 0, min_read_support = 1):
        self.end_fuzz = end_fuzz
        self.min_read_support = min_read_support
        merger.hooks["contig_built"].tap(self.remove_unsupported)

    """
    Get the number of transcripts that support the given target along its full length
    """
    def get_full_length_support(self, target, transcripts):
        return len(
            [
                t for t in transcripts
                if (
                    ranges.within(t.TSS, (target.TSS - self.end_fuzz, target.TSS), exclusive=False)
                    and ranges.within(t.TES, (target.TES, target.TES + self.end_fuzz), exclusive=False)
                    and t.junctions == target.junctions
                )
            ]
        )

    def add_full_length_count(self, target, support):
        target.full_length_count = support
        return target

    def remove_unsupported(self, contig):
        for t in contig.transcripts:
            self.add_full_length_count(t, self.get_full_length_support(t, contig.transcripts))
        
        without_support = [
            t for t in contig.transcripts
            if t.full_length_count < self.min_read_support
        ]
        
        for t in without_support:
            contig.remove_transcript(t)        