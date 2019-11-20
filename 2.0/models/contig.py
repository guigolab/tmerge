from dataclasses import dataclass

@dataclass
class Contig:
    __slots__ = ["alignments", "start", "end"]

    alignments: list
    start: int
    end: int

    def add_alignment(self, alignment):
        if not self.overlaps(alignment):
            raise IndexError("Alignment not within contig.")
        
        self.alignments.append(alignment)

        if alignment.start < self.start:
            self.start = alignment.start

        if alignment.end > self.end:
            self.end = alignment.end

    def overlaps(self, alignment):
        return (self.start <= alignment.end and alignment.start <= self.end) or (alignment.start <= self.end and self.start <= alignment.end)