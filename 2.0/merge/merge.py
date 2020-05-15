from parsers import Importer
from parsers import gtf
from models.transcript_model import TranscriptModel
from output.gtf import write
from functools import reduce, partial
from itertools import product, combinations
from models.contig import Contig
from merge.hook import Hook
import os
from utils import ranges, iterators
from collections import OrderedDict
from merge.rules import ruleset
from multiprocessing import Pool, Manager, Process, cpu_count
from collections import deque

HOOKS = ["input_parsed", "contig_built", "tm_pre_merge", "tm_merged", "tm_not_merged", "contig_merged", "pre_sort", "post_sort", "complete"]
gtf_importer = Importer.Importer(gtf.Gtf())

class Merge:
    def __init__(self, inputPath, outputPath, tolerance = 0, processes=None):
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
        added_ids = []
        items = transcripts.items()

        for first_id, first_transcript in items:
            if first_id in added_ids:
                continue

            cur_contig = [first_transcript]
            cur_TSS = first_transcript.TSS
            cur_TES = first_transcript.TES
            
            added_ids.append(first_id)

            for t_id, transcript in items:
                if t_id in added_ids:
                    continue
                if first_transcript.strand != transcript.strand:
                    # If transcript[i] is not on same strand then the next transcript may be on same strand and overlap so try
                    continue
                if not ranges.overlaps((cur_TSS, cur_TES), (transcript.TSS, transcript.TES)):
                     # If transcript[i] does not overlap then no overlap and on same strand so return
                    break

                cur_contig.append(transcript)
                added_ids.append(t_id)
                if transcript.TES > cur_TES:
                    cur_TES = transcript.TES
                
            self.hooks["contig_built"].exec(cur_contig)
            yield cur_contig
        
        return
            
    """
    Attempts to merge right into left in-place.
    Returns True if merge succeeded, otherwise returns False.
    """
    def merge_transcripts(self, left, right):
        self.hooks["tm_pre_merge"].exec(left, right)
        
        if ruleset(left, right, self.tolerance):
            left.TSS = min([left.TSS, right.TSS])
            left.TES = max([left.TES, right.TES])
            for j in right.junctions:
                left.add_junction(*j)

            left.transcript_count = left.transcript_count + right.transcript_count
            left.contains.extend(right.contains)
            left.contains.append(right.id)

            self.hooks["tm_merged"].exec(left, right)

            return (left, None)

        self.hooks["tm_not_merged"].exec(left, right)

        return (left, right)


    def merge_contig(self, transcripts):
        # first_id = None
        # first_fail_merged_id = None
        # merged = []
        # for left in transcripts:
        #     for right in transcripts:
        #         if left not in merged and right not in merged:
        #             self.merge_transcripts(left, right)
        #     merged.append(left)

        # # while True:
        # #     try:
        # #         left = transcripts.popleft()
        # #         right = transcripts.popleft()
        # #         if not first_id:
        # #             first_id = left.id

        # #         if left.id == first_id and right.id == first_fail_merged_id:
        # #             break

        # #         (merge_result, fail) = self.merge_transcripts(left, right)
        # #         if not fail:
        # #             transcripts.appendleft(merge_result)
        # #         else:
        # #             if not first_fail_merged_id:
        # #                 first_fail_merged_id = fail.id
        # #             transcripts.appendleft(merge_result)
        # #             transcripts.append(fail)
        # #     except IndexError:
        # #         break
        
        # for t in transcripts:
        #     if t.meta["full_length_count"] > 1:
        #         print(t.meta)

        # self.hooks["contig_merged"].exec(transcripts)
        # return transcripts

        # TODO: refactor to use itertools
        i = 0
        i_compare = 1
        while i < len(transcripts) - 1:
            if i == i_compare:
                i_compare += 1
                continue
            
            if i_compare > len(transcripts) - 1:
                i += 1
                i_compare = 0
                continue

            t1 = transcripts[i]
            t2 = transcripts[i_compare]

            (merged, fail) = self.merge_transcripts(t1, t2)

            if not fail:
                transcripts.pop(i_compare)
            else:
                i_compare += 1

        self.hooks["contig_merged"].exec(transcripts)
        return transcripts


    def merge(self):
        transcripts = gtf_importer.parse(self.inputPath)
        self.hooks["input_parsed"].exec(transcripts)
        mg = Manager()
        q = mg.Queue()
        # Put the writer in its own thread
        writer = Process(target=write, args=(q, self.outputPath))
        writer.start()

        with Pool(processes=self.processes) as p:
            contigs = p.imap_unordered(self.merge_contig, self.build_contigs(transcripts))
            
            for contig in contigs:
                # Put the contig to be written in the queue to avoid collisions
                q.put(contig)
            
            p.close()
            p.join()

        q.put("KILL")
        writer.join()
        
        self.hooks["pre_sort"].exec()
        self._sort()
        self.hooks["post_sort"].exec()
        self.hooks["complete"].exec()

    def _sort(self):
        os.system(f"sort -n -k4 -o \"{self.outputPath}\" \"{self.outputPath}\"")
