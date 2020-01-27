from models.exon import Exon
from models.transcript import Transcript

class Importer():
    def __init__(self, parser):
        self.parser = parser

    def parse(self, path):
        transcripts = {}
        prev_chr = ""
        with open(path, "r") as f:
            for line in f:            
                try:
                    data = line.split("\t")
                    exon = Exon(
                        self.parser.get_source(),
                        self.parser.get_chromosome(data),
                        self.parser.get_start(data),
                        self.parser.get_end(data),
                        self.parser.get_strand(data),
                        self.parser.get_transcript_id(data),
                        self.parser.get_gene_id(data)
                    )
                    prev_transcript = transcripts.get(exon.transcript_id, None)
                    
                    if prev_transcript:
                        prev_transcript.add_exon(exon)
                    else:
                        transcripts[exon.transcript_id] = Transcript([exon])

                    # Yield values one chromosome at a time
                    if prev_chr and exon.chromosome != prev_chr:
                        yield list(transcripts.values())
                        transcripts = {}

                    prev_chr = exon.chromosome

                except Exception as e:
                    # TODO: Handle this exception
                    print(e)
                    pass
        
        yield list(transcripts.values())