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

HOOKS = ["chromosome_parsed", "contig_built", "contig_merged", "contig_written", "pre_sort", "post_sort", "complete"]
gtf_importer = Importer.Importer(gtf.Gtf())

class Merge:
    def __init__(self, inputPath, outputPath):
        self._add_hooks()
        
        self.inputPath = inputPath
        self.outputPath = outputPath

        # Overwrite file contents first
        open(self.outputPath, 'w').close()

    def _add_hooks(self):
        self.hooks = {}
        for hook in HOOKS:
            self.hooks[hook] = Hook()

    @staticmethod
    def ruleset(transcript1, transcript2):
        return (
            transcript1.chromosome == transcript2.chromosome
            and transcript1.strand == transcript2.strand
            and ranges.overlaps((transcript1.TSS, transcript1.TES), (transcript2.TSS, transcript2.TES)) # The transcripts overlap
            and (
                ranges.ordered_subset(transcript1.junctions, transcript2.junctions) # Ordered subset?
                or ranges.ordered_subset(transcript2.junctions, transcript1.junctions)
            )
            and (
                not ranges.within_any(transcript1.TSS, transcript2.junctions)
                and not ranges.within_any(transcript1.TES, transcript2.junctions)
            )
            and (
                not ranges.within_any(transcript2.TSS, transcript1.junctions)
                and not ranges.within_any(transcript2.TES, transcript1.junctions)
            )
        )

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
        for contig in self.build_contigs(transcripts):
            # TODO; refactor to use itertools
            merged = {}
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
                
                if self.ruleset(transcripts[i], transcripts[i_compare]):
                    transcripts[i].TSS = min([transcripts[i].TSS, transcripts[i_compare].TSS])
                    transcripts[i].TES = max([transcripts[i].TES, transcripts[i_compare].TES])
                    for j in transcripts[i_compare].junctions:
                        transcripts[i].add_junction(*j)                    
                    transcripts[i].transcript_count = transcripts[i].transcript_count + transcripts[i_compare].transcript_count
                    transcripts.pop(i_compare)
                else:
                    i_compare += 1

            write(transcripts, self.outputPath)
        
        self.hooks["pre_sort"].exec()
        self._sort()
        self.hooks["post_sort"].exec()
        self.hooks["complete"].exec()

    def _sort(self):
        os.system(f"sort -n -k4 -o \"{self.outputPath}\" \"{self.outputPath}\"")
