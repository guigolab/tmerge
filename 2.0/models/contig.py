from dataclasses import dataclass

@dataclass
class Contig:
    __slots__ = ["transcripts", "start", "end"]

    transcripts: list
    start: int
    end: int
    
    def __init__(self, transcripts):
        self.start = min((a.start for a in transcripts))
        self.end = max((a.end for a in transcripts))
        self.transcripts = transcripts

    def add_transcript(self, transcript):
        if not self.overlaps(transcript):
            raise IndexError("Transcript not within contig.")
        
        self.transcripts.append(transcript)

        if transcript.start < self.start:
            self.start = transcript.start

        if transcript.end > self.end:
            self.end = transcript.end

    def overlaps(self, transcript):
        return (self.start < transcript.end and transcript.start < self.end) or (transcript.start < self.end and self.start < transcript.end)