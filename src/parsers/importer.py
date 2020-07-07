from collections import OrderedDict

from ..models.transcript_model import TranscriptModel

class Importer():
    def __init__(self, parser):
        self.parser = parser

    def parse(self, path):
        transcript_models = OrderedDict() # Maintains order as well as uniqueness. Like an ordered set
        with open(path, "r") as f:
            prev_chrom = None
            prev_start = None
            for num, line in enumerate(f):       
                try:
                    if line[0] == "#":
                        # Skip commented lines
                        continue

                    self.parser.set_data(line)
                    
                    if self.parser.get_type() != "exon":
                        continue

                    id = self.parser.get_transcript_id()
                    start = self.parser.get_start()
                    stop = self.parser.get_end()
                    chromosome = self.parser.get_chromosome()
                    strand = self.parser.get_strand()

                    if (prev_start and prev_chrom) and (prev_chrom == chromosome and start < prev_start):
                        raise TypeError("Input must be sorted by chromosome then start position.")

                    if not prev_chrom:
                        prev_chrom = chromosome

                    if prev_chrom != chromosome:
                        # Can yield the chromosome and start afresh
                        yield transcript_models
                        prev_chrom = chromosome
                        transcript_models = OrderedDict()


                    prev_start = start

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
                    raise IOError(f"Error in file parsing at line {num}: {e} \n Line: {line}")

        yield transcript_models