from models.transcript import Transcript
from models.exon import Exon
from models.contig import Contig
import random

class Faker():
    def __init__(self, source_type, chromosome, chromosome_length = 10000):
        self.source_type = source_type
        self.chromosome = chromosome
        self.chromosome_length = chromosome_length

    def exon(self, strand="+", start=None, length=50, transcript_suffix="", gene_suffix=""):
        if not start:
            start = random.randrange(0, self.chromosome_length)

        transcript_id = f"transcript_start{start}_length{length}_{transcript_suffix}"
        gene_id= f"gene_start{start}_length{length}_{gene_suffix}"
        return Exon(self.source_type, self.chromosome, start, start+length, strand, transcript_id, gene_id)

    def transcript(self, num_exons=5, exons_seed=None, **kwargs):
        if not exons_seed:
            exons_seed = [self.exon(**kwargs)]
        transcript = Transcript(exons_seed)
        for i in range(0, num_exons - len(exons_seed)):
            transcript.add_exon(self.exon(**kwargs))

        return transcript

    def contig(self, num_transcripts=10, transcripts_seed=None, **kwargs):
        if not transcripts_seed:
            transcripts_seed = [self.transcript(**kwargs)]
        contig = Contig(transcripts_seed)
        for i in range(0, num_transcripts - len(transcripts_seed)):
            contig.add_transcript(self.transcript(**kwargs))

        return contig
