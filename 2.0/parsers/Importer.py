from models.transcript_model import TranscriptModel
from collections import OrderedDict

class Importer():
    def __init__(self, parser):
        self.parser = parser

    def parse(self, path):
        transcript_models = OrderedDict() # Order is very important so that contigs are built properly
        with open(path, "r") as f:
            prev_chrom = None
            prev_start = None
            for num, line in enumerate(f):       
                try:
                    if line[0] == "#":
                        # Skip commented lines
                        continue

                    data = line.split("\t")
                    self.parser.set_data(data)
                    id = self.parser.get_transcript_id()
                    start = self.parser.get_start()
                    stop = self.parser.get_end()
                    chromosome = self.parser.get_chromosome()
                    strand = self.parser.get_strand()

                    if (prev_start and prev_chrom) and (prev_chrom == chromosome and start < prev_start):
                        raise TypeError("Input must be sorted by chromosome then start position.")

                    prev_start = start
                    prev_chrom = chromosome

                    if self.parser.get_transcript_id() not in transcript_models:
                        transcript_models[id] = TranscriptModel(
                            id,
                            chromosome,
                            strand,
                            start,
                            stop
                        )

                    else:
                        # Add the splice junction to the transcript
                        # Assumes ordered input
                        tm = transcript_models[id]

                        if not stop >= tm.TES:
                            raise IndexError("Input not sorted.")
                        if chromosome != tm.chromosome:
                            raise TypeError("Transcripts spread across multiple chromosomes!")
                        if strand != tm.strand:
                            raise TypeError("Transcripts spread across two strands!")

                        old_TES = tm.TES
                        tm.TES = stop
                        tm.add_junction(old_TES, start)

                except Exception as e:
                    # TODO: Handle this exception
                    raise IOError(f"Error in file parsing at line {num}: {e}")
        
        return transcript_models