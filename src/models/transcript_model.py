from typing import Tuple, Set

from ..utils import ranges

class TranscriptModel:
    __slots__ = ["id", "chromosome", "strand", "_TSS", "_TES", "_junctions", "transcript_count", "contains", "meta", "_removed"]

    id: str
    strand: str
    chromosome: str
    _TSS: int
    _TES: int
    _junctions: Set[Tuple[int, int]] # TODO: used sortedcontainers here
    transcript_count: int
    contains: [str]
    meta: dict
    _removed: bool
    
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
        self._removed = False

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
        return self._junctions

    @junctions.setter
    def junctions(self):
        raise Exception("Junctions cannot be set. Please use junction methods.")

    def add_junction(self, start, stop):
        if not ranges.overlaps((start, stop), (self.TSS, self.TES)):
            raise IndexError("Junction out of range of TSS and TES.")
        if (start, stop) not in self._junctions:
            self._junctions.add((start, stop))

    def remove_junction(self, start, stop):
        if (start, stop) in self._junctions:
            self._junctions.remove((start, stop))

    @property
    def sorted_junctions(self):
        return sorted(self._junctions, key=lambda x: x[0])

    @property
    def removed(self):
        return self._removed

    # Used in plugins. Plugins can only alter transcripts in place so this is used a removal flag
    def remove(self):
        self._removed = True



