from dataclasses import dataclass
from dataclass_type_validator import dataclass_type_validator

@dataclass
class Exon:
    # Use slots to save memory
    __slots__ = ["source_type", "seq_name", "start", "end", "strand", "transcript_id", "gene_id"]

    source_type: str
    seq_name: str
    start: int
    end: int
    strand: str
    transcript_id: str
    gene_id: str

    def __post_init__(self):
        dataclass_type_validator(self)