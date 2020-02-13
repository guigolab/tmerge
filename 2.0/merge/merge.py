from parsers import Importer
from parsers import gtf
from builders.contigs import build as build_contigs
from builders.merged import build as build_merged
from models.transcript_model import TranscriptModel
from output.gtf import write
from functools import reduce
from models.contig import Contig
from merge.hook import Hook
import os
from utils import ranges

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
            and ranges.ordered_subset(transcript1.junctions, transcript2.junctions) # Ordered subset?
            and (not ranges.within_set(transcript1.TSS, transcript2.junctions) or transcript2.monoexonic) # Monoexonics have no junctions so must specify
            and (not ranges.within_set(transcript1.TES, transcript2.junctions) or transcript2.monoexonic)
        )

    def merge(self):
        parsed = gtf_importer.parse(self.inputPath)
        merged = {}

        for tm in parsed:
            if tm.id in merged:
                continue
            # Check if transcript can be merged
            mergeable = [
                x for x in parsed
                    if self.ruleset(tm, x)
                ]

            if not mergeable:
                merged[tm.id] = tm
                continue


            lowest_TSS = min([x.TSS for x in mergeable])
            highest_TES = max([x.TES for x in mergeable])
            longest = max([x.length for x in mergeable])
            new_tm_id = [x for x in mergeable if x.length == longest][0].id # TODO; This may cause collisions
            
            # Build the merged model
            new_tm = TranscriptModel(new_tm_id, tm.chromosome, tm.strand, lowest_TSS, highest_TES)
            for m in mergeable:
                for j in m.junctions:
                    new_tm.add_junction(*j)

            new_tm.transcript_count = len(mergeable)
            if new_tm_id not in merged:
                merged[new_tm_id] = new_tm

        write(list(merged.values()), self.outputPath)
        
        self.hooks["pre_sort"].exec()
        self._sort()
        self.hooks["post_sort"].exec()
        self.hooks["complete"].exec()

    def _sort(self):
        os.system(f"sort -n -k4 -o \"{self.outputPath}\" \"{self.outputPath}\"")
