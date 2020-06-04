from parsers import Importer
from parsers import gtf
from models.transcript_model import TranscriptModel
from output.gtf import write
from functools import reduce, partial
from itertools import combinations, permutations, filterfalse, chain, groupby
from models.contig import Contig
from merge.hook import Hook
import os
from utils import ranges, iterators
from collections import OrderedDict
from merge.rules import ruleset
from multiprocessing import Pool, Manager, Process, cpu_count
from threading import Thread
from queue import Queue
from collections import deque

HOOKS = ["chromosome_parsed", "contig_built", "contig_merged" "merging_complete", "pre_sort", "post_sort", "complete"]
gtf_importer = Importer.Importer(gtf.Gtf())

class Merge:
    def __init__(self, inputPath, outputPath, tolerance = 0, processes=None, speed=False):
        self._add_hooks()
        
        self.inputPath = inputPath
        self.outputPath = outputPath
        self.tolerance = tolerance
        self.processes = processes

        # Overwrite file contents first
        open(self.outputPath, 'w').close()

    def _add_hooks(self):
        self.hooks = {}
        for hook in HOOKS:
            self.hooks[hook] = Hook()

    def build_contigs(self, transcripts):
        while transcripts:
            t1 = next(iter(transcripts.values()))
            cur_TSS = t1.TSS
            cur_TES = t1.TES
            cur_strand = t1.strand
            cur_contig = [t1]
            
            for t2 in transcripts.values():
                if cur_strand != t2.strand:
                    # If t2 is not on same strand then the next transcript may be on same strand and overlap so try
                    continue
                if not ranges.overlaps((cur_TSS, cur_TES), (t2.TSS, t2.TES)):
                    # If t2 does not overlap then no overlap and on same strand so stop iteration
                    break
                
                cur_contig.append(t2)

                if t2.TES > cur_TES:
                    cur_TES = t2.TES

            for t in cur_contig[1:]:
                del transcripts[t.id]

            yield cur_contig
        
            
    """
    Attempts to merge right into left in-place.
    """
    def merge_transcripts(self, left, right):
        if ruleset(left, right, self.tolerance):
            left.TSS = min([left.TSS, right.TSS])
            left.TES = max([left.TES, right.TES])
            for j in right.junctions:
                left.add_junction(*j)

            left.transcript_count = left.transcript_count + right.transcript_count
            left.contains.extend(right.contains)
            left.contains.append(right.id)

            return True

        return False

    def merge_contig(self, transcripts):
        # This might be another method that may be quicker and cleaner but needs investigating as doesnt work quite right. Skips some transcripts
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
        q = Queue()

        # Use multithreading and a FIFO queue for slight speed increase
        Thread(target=write, args=(q, self.outputPath), daemon=True).start()

        for chromosome in gtf_importer.parse(self.inputPath):
            self.hooks["chromosome_parsed"].exec(chromosome)

            for contig in self.build_contigs(chromosome):
                self.hooks["contig_built"].exec(contig)

                unremoved = [x for x in contig if not x.removed]
                merged = self.merge_contig(unremoved)

                self.hooks["contig_merged"].exec(merged)
                
                for transcript in self.merge_contig(unremoved):
                    q.put(transcript)
        
        q.join()
        self.hooks["merging_complete"].exec()

        # Method for using multiprocessing
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
