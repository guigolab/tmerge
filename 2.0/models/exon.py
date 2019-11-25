from dataclasses import dataclass
from dataclass_type_validator import dataclass_type_validator

@dataclass
class Exon:
    # Use slots to save memory
    __slots__ = ["source_type", "chromosome", "start", "end", "strand", "transcript_id", "gene_id"]

    source_type: str
    chromosome: str
    start: int
    end: int
    strand: str
    transcript_id: str
    gene_id: str

    def __post_init__(self):
        dataclass_type_validator(self)
        
        if(self.start >= self.end):
            raise IndexError("start must be less than end.")