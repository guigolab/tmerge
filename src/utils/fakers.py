import uuid
import random

from ..models.transcript_model import TranscriptModel

class Faker():
    """
    Used in unit tests. Creates "fake" transcript models
    """
    def __init__(self, source_type, chromosome, chromosome_length = 10000):
        self.source_type = source_type
        self.chromosome = chromosome
        self.chromosome_length = chromosome_length

    def tm(self, num_introns, intron_length = 50, start = 0, length = 9*9999, strand="+"):
        transcript_id = f"transcript_{uuid.uuid1()}"
        tm = TranscriptModel(transcript_id, self.chromosome, strand, start, start + length)

        for i in range(0, num_introns):
            start = start + random.randrange(2, 10) * 50  # Ensures introns aren't right next to eachother
            tm.add_junction(start, start+intron_length)

        return tm

    def contig(self, num_transcripts=10, transcript_length = 100, introns_per_transcript=0):
        tm = self.tm(introns_per_transcript, length=transcript_length)
        contig = [tm]
        tm_start = 0
        for i in range(0, num_transcripts - 1):
            contig.append(self.tm(introns_per_transcript, start=tm_start))
            tm_start = tm_start + transcript_length

        return contig

    def contig_overlapping_monoexonic(self, num_transcripts=10, **kwargs):
        transcripts_seed = []
        start = 0
        length = kwargs.get("length") or 50
        for i in range(0, num_transcripts):
            transcripts_seed.append(self.transcript(1, None, start=start, **kwargs))
            start = start + length - 1
        
        return self.contig(num_transcripts, transcripts_seed)