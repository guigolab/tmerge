class Stats:
    """
    Print basic stats to stdout
    """
    def __init__(self, hooks, stats, **kwargs):
        if not stats:
            return

        self.transcript_count = 0
        self.chromosome_count = 0
        self.contig_count = 0
        self.merged_count = 0

        hooks["chromosome_parsed"].tap(self.add_chromosome)
        hooks["contig_built"].tap(self.add_contig)
        hooks["contig_merged"].tap(self.add_merged_transcripts)
        hooks["complete"].tap(self.print)

    def add_chromosome(self, transcripts):
        self.chromosome_count += 1
        self.transcript_count += len(transcripts)

    def add_contig(self, *args):
        self.contig_count += 1

    def add_merged_transcripts(self, transcripts):
        self.merged_count += len(transcripts)

    def print(self):
        print(f"Found {self.transcript_count} in {self.chromosome_count} chromosomes.")
        print(f"Built {self.contig_count} contigs.")
        print(f"Final number of transcript models: {self.merged_count}")