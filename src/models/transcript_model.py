from typing import Tuple, Set

from ..utils import ranges

class TranscriptModel:
    """
    Represents a transcript as it flows through tmerge.

    At first this is just the transcript from the input reads but after merging it "contains" other transcripts, making it a "transcript model". 
    When two transcripts are merged, the junctions of the merged in transcript is added to the host and the TSS/TES are updated to be the minimum/maximum of the two.
    """
    __slots__ = ["id", "chromosome", "strand", "_TSS", "_TES", "_junctions", "transcript_count", "meta", "_removed"]

    id: str
    strand: str
    chromosome: str
    # 5' end
    _TSS: int
    # 3' end
    _TES: int
    _junctions: Set[Tuple[int, int]] # TODO: used sortedcontainers here
    # Any key/value pairs added to meta will be added to the attributes column of the output GTF
    meta: dict
    _removed: bool
    
    def __init__(self, id, chromosome, strand, TSS, TES):
        if TSS > TES:
            raise IndexError("TSS must be less than TES.")
        
        self.id = id
        self.chromosome = chromosome
        self.strand = strand
        self._TSS = TSS
        self._TES = TES
        self._junctions = set()
        self.transcript_count = 1
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
        """
        Add a junction to the TranscriptModel

        Parameters
        ----------
        start: int
            5' coordinate
        end: init
            3' coordinate
        """
        if not ranges.overlaps((start, stop), (self.TSS, self.TES)):
            raise IndexError("Junction out of range of TSS and TES.")
        if (start, stop) not in self._junctions:
            self._junctions.add((start, stop))

    def remove_junction(self, start, stop):
        """
        Remove a junction from the TranscriptModel

        Parameters
        ----------
        start: int
            5' coordinate
        end: init
            3' coordinate
        """
        if (start, stop) in self._junctions:
            self._junctions.remove((start, stop))

    @property
    def sorted_junctions(self):
        # Sets can't be sorted so if a sorted list is needed, convert it to a list and sort it by the 5' end of the junction
        return sorted(self._junctions, key=lambda x: x[0])

    @property
    def removed(self):
        return self._removed

    def remove(self):
        """
        Flag the TranscriptModel for removal.

        This is how plugins remove transcripts. In Merge.merge, transcripts are filtered out that have been flagged as removed.
        """
        self._removed = True



