from functools import reduce

class SplicedLengths():
    def __init__(self, hooks, **kwargs):
        hooks["contig_merged"].tap(self.add_spliced_lengths)

    def add_spliced_lengths(self, transcripts):
        for transcript in transcripts:
            # mRNA length. Add one to the total length and each of the junctions since junctions are coordinates of start and end of exons, not start and end of introns
            transcript.meta["mature_RNA_length"] = reduce(lambda l, x: l - (x[1] - x[0]) + 1, transcript.junctions, transcript.length + 1)
