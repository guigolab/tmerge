from dataclasses import dataclass
from utils import ranges

@dataclass
class Contig:
    __slots__ = ["transcripts", "start", "end", "strand"]

    transcripts: dict
    start: int
    end: int
    strand: str
    
    def __init__(self, transcripts):
        if not all([t.strand for t in transcripts.values()]):
            raise TypeError("Exons must be on the same strand.")

        self.start = min((a.TSS for a in transcripts.values()))
        self.end = max((a.TES for a in transcripts.values()))
        self.transcripts = transcripts
        self.strand = next(iter(transcripts.values())).strand

    def add_transcript(self, transcript):
        if not self.are_same_strand(transcript):
            raise TypeError("Transcripts not on same strand.")
        if not self.overlaps(transcript):
            raise IndexError("Transcript not within contig.")
        
        self.transcripts[transcript.id] = transcript
        
        if transcript.TES > self.end:
            self.end = transcript.TES

    def overlaps(self, transcript):
        return ranges.overlaps((self.start, self.end), (transcript.TSS, transcript.TES))

    def is_plus(self):
        return self.strand == "+"

    def is_minus(self):
        return self.strand == "-"
        
    def are_same_strand(self, transcript):
        return self.strand == transcript.strand
        
    def has_transcript(self, transcript):
        return transcript in self.transcripts