class MinLength:
    def __init__(self, hooks, min_length=0, **kwargs):
        """
        Remove TMs that do not satisfy min_length.
        """
        self.min_length = min_length
        if self.min_length > 0:
            hooks["contig_merged"].tap(self.remove_short_transcripts)

    def remove_short_transcripts(self, transcripts):
        for t in transcripts:
            if t.length < self.min_length:
                t.remove()