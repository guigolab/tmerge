from dataclasses import dataclass

@dataclass
class Transcript:
    __slots__ = ["exons", "start", "end"]

    exons: list
    start: int
    end: int