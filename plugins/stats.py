class Stats:
    def __init__(self, hooks, stats, **kwargs):
        if not stats:
            return

        self.transcript_count = 0
        self.contig_count = 0
        self.merged_count = 0

        hooks["input_parsed"].tap(self.add_transcripts)
        hooks["contig_built"].tap(self.add_contig)
        hooks["contig_merged"].tap(self.add_merged_transcripts)
        hooks["complete"].tap(self.print)

    def add_transcripts(self, transcripts):
        self.transcript_count += len(transcripts)

    def add_contig(self, *args):
        self.contig_count += 1

    def add_merged_transcripts(self, contig):
        self.merged_count += len(contig.transcripts)

    def print(self):
        print(f"Number of contigs: {self.contig_count}")
        print(f"Number of transcripts: {self.transcript_count}")
        print(f"Number of merged transcripts: {self.merged_count}")