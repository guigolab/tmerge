from src.utils import ranges
import statistics

class Depth:
    """
    Define depth for each transcript and trim TMs when read depth drops.
    """
    def __init__(self, hooks, trim = False, **kwargs):
        hooks["transcript_added"].tap(self.add_attributes)
        hooks["contig_merged"].tap(self.compute_depth)
        if trim:
            hooks["contig_merged"].tap(self.trim)

    def add_attributes(self, transcript):
        transcript.meta["depth"] = []

    def compute_depth(self, transcripts):
        for t in transcripts:
            for pos in range(t.TSS, t.TES + 1):
                t.meta["depth"].append(
                    len([
                        x for x in transcripts if ranges.within(pos, (x.TSS, x.TES))
                    ])
                )

    def trim(self, transcripts):
        for t in transcripts:
            self.trim_transcript(t)

    def trim_transcript(self, transcript):
        depths = transcript.meta["depth"]
        cutoff_depth = statistics.median_low(depths)

        for pos, depth in enumerate(depths):
            if depth >= cutoff_depth:
                continue

            if transcript.strand == "+":
                # Trim the 3' end
                to_remove = []
                for j in transcript.junctions:
                    if j[0] >= pos:
                        to_remove.append(j)
                for j in to_remove:
                    transcript.remove_junction(*j)
                transcript.TES = pos
                break

            elif transcript.strand == "-":
                # trim the 5' end
                to_remove = []
                for j in transcript.junctions:
                    if j[1] <= pos:
                        to_remove.append(j)
                for j in to_remove:
                    transcript.remove_junction(*j)
                
                transcript.TSS = pos
                break
