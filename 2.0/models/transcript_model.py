from dataclasses import dataclass
from typing import List, MutableSet, Tuple
from dataclass_type_validator import TypeValidationError
from utils import ranges
from functools import partial


class TranscriptModel:
    __slots__ = ["id", "chromosome", "strand", "_TSS", "_TES", "_junctions", "transcript_count", "contains", "meta"]

    id: str
    strand: str
    chromosome: str
    _TSS: int
    _TES: int
    _junctions: MutableSet[Tuple[int, int]]
    transcript_count: int
    contains: [str]
    meta: dict
    
    def __init__(self, id, chromosome, strand, TSS, TES):
        if TSS >= TES:
            raise IndexError("TSS must be less than TES.")
        
        self.id = id
        self.chromosome = chromosome
        self.strand = strand
        self._TSS = TSS
        self._TES = TES
        self._junctions = set()
        self.transcript_count = 1
        self.contains = [self.id]
        self.meta = {}

    @property
    def length(self):
        return self.TES - self.TSS

    @property
    def monoexonic(self):
        return len(self._junctions) == 0

    @property
    def TSS(self):
        return self._TSS

    @TSS.setter
    def TSS(self, TSS):
        if not ranges.lte_all(TSS, self.junctions):
            raise IndexError("TSS is not less than or equal to the models junctions.")

        self._TSS = TSS

    @property
    def TES(self):
        return self._TES

    @TES.setter
    def TES(self, TES):
        if not ranges.gte_all(TES, self.junctions):
            raise IndexError("TES is not greater than or equal to the models junctions.")

        self._TES = TES

    @property
    def junctions(self):
        return sorted(list(self._junctions), key=lambda x: x[0]) # Sort by start

    @junctions.setter
    def junctions(self):
        raise Exception("Junctions cannot be set. Please use junction methods.")

    def add_junction(self, start, stop):
        if not ranges.overlaps((start, stop), (self.TSS, self.TES)):
            raise IndexError("Junction out of range of TSS and TES.")
        self._junctions.update({(start, stop)})

    def remove_junction(self, start, stop):
        self._junctions.remove((start, stop))


