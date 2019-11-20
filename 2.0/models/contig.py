from dataclasses import dataclass

@dataclass
class Contig:
    __slots__ = ["exons", "start", "end"]

    exons: list
    start: int
    end: int
    
    def __init__(self, exons):
        self.start = min((a.start for a in exons))
        self.end = max((a.end for a in exons))
        self.exons = exons

    def add_exon(self, exon):
        if not self.overlaps(exon):
            raise IndexError("Exon not within contig.")
        
        self.exons.append(exon)

        if exon.start < self.start:
            self.start = exon.start

        if exon.end > self.end:
            self.end = exon.end

    def overlaps(self, exon):
        return (self.start <= exon.end and exon.start <= self.end) or (exon.start <= self.end and self.start <= exon.end)