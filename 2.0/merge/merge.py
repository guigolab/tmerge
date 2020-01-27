from parsers import Importer
from parsers import gtf
from builders.contigs import build as build_contigs
from builders.merged import build as build_merged
from output.gtf import write
from functools import reduce
from models.contig import Contig
from merge.hook import Hook
import os

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

    def merge(self):
        for transcripts in gtf_importer.parse(self.inputPath):
            self.hooks["chromosome_parsed"].exec(transcripts)

            for i, contig in enumerate(build_contigs(transcripts)):
                self.hooks["contig_built"].exec(contig)
                merged = build_merged(contig)
                self.hooks["contig_merged"].exec(merged)
                write(merged, f"contig_{i}", self.outputPath)
                self.hooks["contig_written"].exec(contig)
        
        self.hooks["pre_sort"].exec()
        self._sort()
        self.hooks["post_sort"].exec()
        self.hooks["complete"].exec()

    def _sort(self):
        os.system(f"sort -n -k4 -o \"{self.outputPath}\" \"{self.outputPath}\"")
