from functools import partial
from itertools import combinations
import os
from multiprocessing import Pool, Manager, Process, cpu_count
from threading import Thread
from queue import Queue

from ..parsers import Importer, Gtf
from ..models import TranscriptModel
from ..output import gtf as write
from ..utils import ranges
from .rules import ruleset
from .hook import Hook

# Any hooks added here should also be updated in the docs
HOOKS = ["chromosome_parsed", "transcript_added", "contig_built", "transcripts_merged", "contig_merged", "contig_complete", "merging_complete", "pre_sort", "post_sort", "complete"]

# For now only support gtf
gtf_importer = Importer(Gtf())

class Merge:
    """
    Controls merging. This is the main class of tmerge.
    """
    def __init__(self, inputPath, outputPath, tolerance = 0, processes=None, no_merge=False):
        self._add_hooks()
        
        self.inputPath = inputPath
        self.outputPath = outputPath
        self.tolerance = tolerance
        self.processes = processes
        self.no_merge = no_merge

        # Overwrite file contents first
        open(self.outputPath, 'w').close()

    def _add_hooks(self):
        self.hooks = {}
        for hook in HOOKS:
            self.hooks[hook] = Hook()

    def build_contigs(self, transcripts):
        """
        Build sets of overlapping transcripts. 

        Reduces the search space for merging.

        Parameters
        ----------
        transcripts: list
            List of TranscriptModel

        Returns
        -------
        generator
            A generator that yield lists of TranscriptModel (contigs)
        """
        while transcripts:
            t1 = next(iter(transcripts.values()))
            cur_TSS = t1.TSS
            cur_TES = t1.TES
            cur_strand = t1.strand
            cur_contig = [t1]
            self.hooks["transcript_added"].exec(t1)
            
            for t2 in transcripts.values():
                if cur_strand != t2.strand:
                    # If t2 is not on same strand then the next transcript may be on same strand and overlap so try
                    continue
                if not ranges.overlaps((cur_TSS, cur_TES), (t2.TSS, t2.TES)):
                    # If t2 does not overlap then no overlap and on same strand so stop iteration
                    break
                
                cur_contig.append(t2)
                self.hooks["transcript_added"].exec(t2)

                if t2.TES > cur_TES:
                    cur_TES = t2.TES

            for t in cur_contig[1:]:
                del transcripts[t.id]

            yield cur_contig
            
    def merge_transcripts(self, left, right):
        """
        Attempts to merge right into left in-place.
        """
        if ruleset(left, right, self.tolerance):
            left.TSS = min([left.TSS, right.TSS])
            left.TES = max([left.TES, right.TES])
            for j in right.junctions:
                left.add_junction(*j)

            self.hooks["transcripts_merged"].exec(left, right)

            return True

        return False

    def merge_contig(self, transcripts):
        """
        Merge a given contig sequentially from left to right until no more merging can be done.

        Parameters
        ----------
        transcripts: list
            List of TranscriptModel

        Returns
        -------
        generator
            Generator that yields a merged contig (list of TranscriptModel)
        """
        # This might be another method that may be quicker and cleaner but needs investigating as doesnt work quite right. Skips some transcripts
        # NOTE: transcripts should be a dequeue not a list with this method
        # if len(transcripts) == 1:
        #     # No need to merge if <= one in contig
        #     yield transcripts[0]
        #     return

        # merged = set()
        # last = None
        # for (left, right) in combinations(transcripts, 2):
        #     if left in merged or right in merged:
        #         continue
        #     if not last:
        #         last = left
        #     if left is not last:
        #         yield last
        #         last = left
        #     if self.merge_transcripts(left, right):
        #         merged.add(right)

        # if last:
        #     yield last

        i = 0
        i_compare = 1
        while i < len(transcripts) - 1:
            if i == i_compare:
                i_compare += 1
                continue
            
            t1 = transcripts[i]

            if i_compare > len(transcripts) - 1:
                i += 1
                i_compare = 0
                continue

            t2 = transcripts[i_compare]

            if self.merge_transcripts(t1, t2):
                transcripts.pop(i_compare)
            else:
                i_compare += 1

        return transcripts

    def merge(self):
        """
        The main function that begins tmerge.
        """
        q = Queue()

        # Use multithreading and a FIFO queue for slight speed increase
        Thread(target=write, args=(q, self.outputPath), daemon=True).start()

        for chromosome in gtf_importer.parse(self.inputPath):
            self.hooks["chromosome_parsed"].exec(chromosome)

            for contig in self.build_contigs(chromosome):
                self.hooks["contig_built"].exec(contig)

                unremoved = [x for x in contig if not x.removed]
                if self.no_merge:
                    merged = unremoved
                else:
                    merged = self.merge_contig(unremoved)

                self.hooks["contig_merged"].exec(merged)

                unremoved_and_merged = [x for x in merged if not x.removed]
                for transcript in unremoved_and_merged:
                    q.put(transcript)

                self.hooks["contig_complete"].exec(unremoved_and_merged)
        
        q.join()
        self.hooks["merging_complete"].exec()

        # Method for using multiprocessing
        # NOTE: It doesn't work very well at the moment but has potential to create vast speed improvements
        # mg = Manager()
        # q = mg.Queue()

        # writer = Process(target=write, args=(q, self.outputPath))

        # num_processes = self.processes if self.processes else cpu_count()
        # print(f"Running on {num_processes} threads.")
        
        # contigs = list(self.build_contigs(transcripts))
        # print(f"Built {len(contigs)} contigs")

        # pool = []
        # for contig in contigs:
        #     p = Process(target=self.merge_contig, args = (contig,q))
        #     p.start()
        #     pool.append(p)

        # for p in pool:
        #     p.join()

        # q.put("DONE")
        # writer.join()

        # # with Pool(processes=num_processes) as p:
        # #     p.map(partial(self.merge_contig, q=q), contigs)
        # #     q.put("DONE")
            
        # # write(q, self.outputPath)

        self.hooks["pre_sort"].exec()
        self._sort()
        self.hooks["post_sort"].exec()
        self.hooks["complete"].exec()

    def _sort(self):
        os.system(f"sort -k 1,1 -k4,4n -o \"{self.outputPath}\" \"{self.outputPath}\"")
