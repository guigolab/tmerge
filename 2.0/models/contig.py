from dataclasses import dataclass
from utils import ranges

@dataclass
class Contig:
    __slots__ = ["_transcripts", "start", "end", "strand"]

    _transcripts: dict
    start: int
    end: int
    strand: str
    
    def __init__(self, seed_transcript):
        self.start = seed_transcript.TSS
        self.end = seed_transcript.TES
        self.strand = seed_transcript.strand
        self._transcripts = {}

        self.add_transcript(seed_transcript)

    @property
    def transcripts(self):
        return list(self._transcripts.values())

    @transcripts.setter
    def transcripts(self, transcript):
        raise Error("Use add_transcript instead")

    def add_transcript(self, transcript):
        if not self.are_same_strand(transcript):
            raise TypeError("Transcripts not on same strand.")
        if not self.overlaps(transcript):
            raise IndexError("Transcript not within contig.")
        
        self._transcripts[transcript.id] = transcript
        
        if transcript.TES > self.end:
            self.end = transcript.TES

    def remove_transcript(self, transcript):
        del self._transcripts[transcript.id]

    def overlaps(self, transcript):
        return ranges.overlaps((self.start, self.end), (transcript.TSS, transcript.TES))

    def is_plus(self):
        return self.strand == "+"

    def is_minus(self):
        return self.strand == "-"
        
    def are_same_strand(self, transcript):
        return self.strand == transcript.strand
        
    def has_transcript(self, transcript):
        return transcript in self._transcripts