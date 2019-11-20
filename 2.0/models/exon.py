from dataclasses import dataclass

@dataclass
class Exon:
    # Use slots to save memory
    __slots__ = ["source_type", "seq_name", "start", "end", "strand", "attributes"]

    source_type: str
    seq_name: str
    start: int
    end: int
    strand: str
    attributes: str