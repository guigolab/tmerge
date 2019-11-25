from dataclasses import dataclass
from typing import List
from exon import Exon
from dataclass_type_validator import TypeValidationError


@dataclass
class Transcript:
    __slots__ = ["exons", "start", "end", "id", "strand"]

    exons: List[Exon]
    start: int
    end: int
    id: str
    strand: str

    
    def __init__(self, exons):

        self.start = min((a.start for a in exons))
        self.end = max((a.end for a in exons))
        self.exons = exons
        self.id = exons[0].transcript_id
        self.strand = exons[0].strand

        if not all([self.are_same_transcript(e) for e in exons]):
            raise TypeError("Exons are not from the same transcript.")
        if not all([self.are_same_strand(e) for e in exons]):
            raise TypeError("Exons must be on the same strand.")


    def __post_init__(self):
        dataclass_type_validator(self)
        
        if(self.start >= self.end):
            raise IndexError("start must be less than end.")

    def add_exon(self, exon):
        if not self.are_same_transcript(exon):
            raise TypeError(f"Not same transcript. {self.id} and {exon.transcript_id}")
        if not self.are_same_strand(exon):
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