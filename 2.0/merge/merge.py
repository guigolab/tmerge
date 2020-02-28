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

HOOKS = ["input_parsed", "contig_built", "contig_merged", "contig_written", "pre_sort", "post_sort", "complete"]
gtf_importer = Importer.Importer(gtf.Gtf())

class Merge:
    def __init__(self, inputPath, outputPath, tolerance = 0, end_fuzz = 0, min_read_support = 1):
        self._add_hooks()
        
        self.inputPath = inputPath
        self.outputPath = outputPath
        self.tolerance = tolerance
        self.end_fuzz = end_fuzz
        self.min_read_support = min_read_support

        # Overwrite file contents first
        open(self.outputPath, 'w').close()

    def _add_hooks(self):
        self.hooks = {}
        for hook in HOOKS:
            self.hooks[hook] = Hook()

    def build_contigs(self, transcripts):
        while transcripts:
            id, transcript = next(iter(transcripts.items()))
            contig_transcripts = OrderedDict()
            contig_transcripts[id] = transcript
            new_contig = Contig(contig_transcripts)
            
            for i, transcript in enumerate(transcripts.values()):
                try:
                    new_contig.add_transcript(transcript)
                except TypeError:
                    # If transcript[i] is not on same strand then the next transcript may be on same strand and overlap so try
                    continue
                except IndexError:
                    # If transcript[i] does not overlap then no overlap and on same strand so return
                    break
            
            for k, v in new_contig.transcripts.items():
                del transcripts[k]

            yield new_contig

    def merge(self):
        transcripts = gtf_importer.parse(self.inputPath)
        self.hooks["input_parsed"].exec(transcripts)

        for contig in self.build_contigs(transcripts):
            self.hooks["contig_built"].exec(contig)

            # TODO; refactor to use itertools
            i = 0
            i_compare = 1
            transcripts = list(contig.transcripts.values())
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

                if ruleset(t1, t2, self.tolerance):
                    t1.full_length_count = max([t1, t2], key=lambda x: x.length).full_length_count
                    t1.TSS = min([t1.TSS, t2.TSS])
                    t1.TES = max([t1.TES, t2.TES])
                    for j in t2.junctions:
                        t1.add_junction(*j)

                    t1.transcript_count = t1.transcript_count + t2.transcript_count

                    if ranges.within(t2.TSS, (t1.TSS - self.end_fuzz, t1.TSS), exclusive=False) and ranges.within(t2.TES, (t1.TES, t1.TES + self.end_fuzz), exclusive=False):
                        t1.full_length_count += 1

                    transcripts.pop(i_compare)
                else:
                    i_compare += 1

            transcripts = [t for t in transcripts if t.full_length_count >= self.min_read_support]
            contig.transcripts = transcripts

            self.hooks["contig_merged"].exec(contig)

            write(transcripts, self.outputPath)

            self.hooks["contig_written"].exec()

        
        self.hooks["pre_sort"].exec()
        self._sort()
        self.hooks["post_sort"].exec()
        self.hooks["complete"].exec()

    def _sort(self):
        os.system(f"sort -n -k4 -o \"{self.outputPath}\" \"{self.outputPath}\"")
