from dataclasses import dataclass
from typing import List
from models.exon import Exon
from dataclass_type_validator import TypeValidationError


@dataclass
class Transcript:
    __slots__ = ["exons", "start", "end", "id", "strand", "chromosome"]

    exons: List[Exon]
    start: int
    end: int
    id: str
    strand: str
    chromosome: str

    
    def __init__(self, exons):

        self.start = min((a.start for a in exons))
        self.end = max((a.end for a in exons))
        self.exons = exons
        self.id = exons[0].transcript_id
        self.strand = exons[0].strand
        self.chromosome = exons[0].chromosome

        def check_all(func):
            return not all([func(e) for e in exons])

        if check_all(self.are_same_strand):
            raise TypeError("Exons must be on the same strand.")
        if check_all(self.are_same_chromosome):
            raise TypeError("Exons must be on the same chromosome.")
        

    def __post_init__(self):
        dataclass_type_validator(self)
        
        if(self.start >= self.end):
            raise IndexError("start must be less than end.")

    def add_exon(self, exon):
        if not self.are_same_strand(exon):
            raise TypeError("Not same strand.")
        if not self.are_same_chromosome(exon):
            raise TypeError("Not same chromosome.")

        self.exons.append(exon)

        if exon.start < self.start:
            self.start = exon.start

        if exon.end > self.end:
            self.end = exon.end

    def are_same_strand(self, exon):
        return self.strand == exon.strand

    def are_same_chromosome(self, exon):
        return self.chromosome == exon.chromosome