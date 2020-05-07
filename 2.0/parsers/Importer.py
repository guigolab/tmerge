from models.transcript_model import TranscriptModel
from collections import OrderedDict

class Importer():
    def __init__(self, parser):
        self.parser = parser

    def parse(self, path):
        transcript_models = OrderedDict() # Order is very important so that contigs are built properly
        with open(path, "r") as f:
            for line in f:            
                # try:
                data = line.split("\t")
                self.parser.set_data(data)
                id = self.parser.get_transcript_id()
                start = self.parser.get_start()
                stop = self.parser.get_end()

                if self.parser.get_transcript_id() not in transcript_models:
                    transcript_models[id] = TranscriptModel(
                        id,
                        self.parser.get_chromosome(),
                        self.parser.get_strand(),
                        start,
                        stop
                    )

                else:
                    # Add the splice junction to the transcript
                    # Assumes ordered input
                    tm = transcript_models[id]

                    if not stop >= tm.TES:
                        raise IndexError("Input not sorted.")
                    if self.parser.get_chromosome() != tm.chromosome:
                        raise TypeError("Transcripts spread across multiple chromosomes!")
                    if self.parser.get_strand() != tm.strand:
                        raise TypeError("Transcripts spread across two strands!")

                    old_TES = tm.TES
                    tm.TES = stop
                    tm.add_junction(old_TES, start)

                # except Exception as e:
                #     # TODO: Handle this exception
                #     print(e)
                #     pass
        
        return transcript_models