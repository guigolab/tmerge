from dataclasses import dataclass

@dataclass
class Alignment:
    # Use slots to save memory
    __slots__ = ["source_type", "seq_name", "start", "end", "strand"]

    source_type: str
    seq_name: str
    start: int
    end: int
    strand: str