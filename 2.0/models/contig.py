from dataclasses import dataclass

@dataclass
class Contig:
    __slots__ = ["transcripts", "start", "end", "strand"]

    transcripts: list
    start: int
    end: int
    strand: str
    
    def __init__(self, transcripts):
        if not all([t.strand for t in transcripts]):
            raise TypeError("Exons must be on the same strand.")

        self.start = min((a.start for a in transcripts))
        self.end = max((a.end for a in transcripts))
        self.transcripts = transcripts
        self.strand = transcripts[0].strand

    def add_transcript(self, transcript):
        if not self.are_same_strand(transcript):
            raise TypeError("Transcripts not on same strand.")
        if not self.overlaps(transcript):
            raise IndexError("Transcript not within contig.")
        
        self.transcripts.append(transcript)
        
        if transcript.end > self.end:
            self.end = transcript.end
        if transcript.start < self.start:
            self.start = transcript.start

    def overlaps(self, transcript):
        return self.end >= transcript.start and transcript.end >= self.start

    def is_plus(self):
        return self.strand == "+"

    def is_minus(self):
        return self.strand == "-"
        
    def are_same_strand(self, transcript):
        return self.strand == transcript.strand
        
    def has_transcript(self, transcript):
        return transcript in self.transcripts