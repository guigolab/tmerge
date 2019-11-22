from dataclasses import dataclass

@dataclass
class Transcript:
    __slots__ = ["exons", "start", "end", "id", "strand"]

    exons: list
    start: int
    end: int
    id: str
    strand: str

    
    def __init__(self, exons):
        if not all([a.transcript_id for a in exons]):
            raise TypeError("Exons are not from the same transcript.")
        if not all(e.strand for e in exons):
            raise TypeError("Exons must be on the same strand.")

        self.start = min((a.start for a in exons))
        self.end = max((a.end for a in exons))
        self.exons = exons
        self.id = exons[0].transcript_id
        self.strand = exons[0].strand

    def add_exon(self, exon):
        if not self.are_same_transcript(exon):
            raise TypeError(f"Not same transcript. {self.id} and {exon.transcript_id}")
        if not self.are_same_transcript(exon):
            raise TypeError("Not same strand.")

        self.exons.append(exon)

        if exon.start < self.start:
            self.start = exon.start

        if exon.end > self.end:
            self.end = exon.end

    def are_same_transcript(self, exon):
        return self.id == exon.transcript_id

    def are_same_strand(self, exon):
        return self.strand == exon.strand