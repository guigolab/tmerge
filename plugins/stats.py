import logging

class Stats:
    """
    Print basic stats to stdout
    """
    def __init__(self, hooks, stats, stats_output,**kwargs):
        if not stats:
            return

        if stats_output:
            open(stats_output, "w")

        logging.basicConfig(
            filename=None if stats_output == "-" else stats_output,
            level=logging.INFO,
            format='%(levelname)s: %(message)s'
        )

        self.transcript_count = 0
        self.chromosome_count = 0
        self.contig_count = 0
        self.merged_count = 0

        hooks["chromosome_parsed"].tap(self.add_chromosome)
        hooks["contig_built"].tap(self.add_contig)
        hooks["transcripts_merged"].tap(self.log_merging_event)
        hooks["contig_complete"].tap(self.add_merged_transcripts)
        hooks["complete"].tap(self.print)

    def add_chromosome(self, transcripts):
        self.chromosome_count += 1
        self.transcript_count += len(transcripts)

    def add_contig(self, *args):
        self.contig_count += 1

    def add_merged_transcripts(self, transcripts):
        self.merged_count += len(transcripts)

    def log_merging_event(self, left, right):
        logging.info(f"{right.id} merged into {left.id}.")

    def print(self):
        logging.info(f"Found {self.transcript_count} transcripts in {self.chromosome_count} chromosomes.")
        logging.info(f"Built {self.contig_count} contigs.")
        logging.info(f"Final number of transcript models: {self.merged_count}")