from dataclasses import dataclass

@dataclass
class Transcript:
    __slots__ = ["exons", "start", "end", "id"]

    exons: list
    start: int
    end: int
    id: str

    
    def __init__(self, exons):
        if len(exons) > 1 and not all(a.transcript_id for a in exons):
            raise TypeError("Exons are not from the same transcript.")

        self.start = min((a.start for a in exons))
        self.end = max((a.end for a in exons))
        self.exons = exons
        self.id = exons[0].transcript_id

    def add_exon(self, exon):
        if not self.are_same_transcript(exon):
            raise TypeError(f"Not same transcript. {self.id} and {exon.transcript_id}")

        self.exons.append(exon)

        if exon.start < self.start:
            self.start = exon.start

        if exon.end > self.end:
            self.end = exon.end

    def are_same_transcript(self, exon):
        return self.id == exon.transcript_id