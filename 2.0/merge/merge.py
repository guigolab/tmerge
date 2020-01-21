from importers import gtf
from builders.contigs import build as build_contigs
from builders.merged import build as build_merged
from output.gtf import write
from functools import reduce
from models.contig import Contig

class Merge:
    def __init__(self, inputPath, outputPath):
        self.inputPath = inputPath
        self.outputPath = outputPath

        # Overwrite file contents first
        open(self.outputPath, 'w').close()

    def merge(self):
        for transcripts in gtf.parse(self.inputPath):
            for i, contig in enumerate(build_contigs(transcripts)):
                merged = build_merged(contig)
                write(merged, f"contig_{i}", self.outputPath)
