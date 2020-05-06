from parsers import Importer
from parsers import gtf
from models.transcript_model import TranscriptModel
from output.gtf import write
from functools import reduce, partial
from itertools import product
from models.contig import Contig
from merge.hook import Hook
import os
from utils import ranges, iterators
from collections import OrderedDict
from merge.rules import ruleset
from multiprocessing import Pool, Manager, Process, cpu_count

HOOKS = ["input_parsed", "contig_built", "contig_merged", "pre_sort", "post_sort", "complete"]
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
        # Note: transripts parameter here is a dict but contig.transcripts is a list
        while transcripts:
            id, transcript = next(iter(transcripts.items()))
            new_contig = Contig(transcript)

            for i, transcript in enumerate(transcripts.values(), start=1):
                try:
                    new_contig.add_transcript(transcript)
                except TypeError:
                    # If transcript[i] is not on same strand then the next transcript may be on same strand and overlap so try
                    continue
                except IndexError:
                    # If transcript[i] does not overlap then no overlap and on same strand so return
                    break

            # TODO: Should not be accessing _transcripts
            for k, v in new_contig._transcripts.items():
                del transcripts[k]

            self.hooks["contig_built"].exec(new_contig)
            yield new_contig
        
        return
            

    def merge_contig(self, contig):
        # TODO: refactor to use itertools
        i = 0
        i_compare = 1
        merged_transcripts = contig.transcripts
        while i < len(merged_transcripts) - 1:
            if i == i_compare:
                i_compare += 1
                continue
            
            if i_compare > len(merged_transcripts) - 1:
                i += 1
                i_compare = 0
                continue
            
            t1 = merged_transcripts[i]
            t2 = merged_transcripts[i_compare]

            if ruleset(t1, t2, self.tolerance):
                t1.TSS = min([t1.TSS, t2.TSS])
                t1.TES = max([t1.TES, t2.TES])
                for j in t2.junctions:
                    t1.add_junction(*j)

                t1.transcript_count = t1.transcript_count + t2.transcript_count

                # Will replace
                contig.add_transcript(t1)

                contig.remove_transcript(t2)
                merged_transcripts.pop(i_compare)
            else:
                i_compare += 1

        self.hooks["contig_merged"].exec(contig)

        return contig

    def merge(self):
        transcripts = gtf_importer.parse(self.inputPath) # dict
        # TODO: refactor here so don't have to convert to list. Best to use an iterator for transcripts everywhere
        self.hooks["input_parsed"].exec(list(transcripts.values()))
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
