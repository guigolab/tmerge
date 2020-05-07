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
        added_ids = []
        for i, transcript in enumerate(transcripts):
            if transcript.id in added_ids:
                continue

            new_contig = Contig(transcript)
            added_ids.append(transcript.id)

            for i, transcript in enumerate([t for t in transcripts if t.id not in added_ids], start=1):
                try:
                    new_contig.add_transcript(transcript)
                    added_ids.append(transcript.id)
                except TypeError:
                    # If transcript[i] is not on same strand then the next transcript may be on same strand and overlap so try
                    continue
                except IndexError:
                    # If transcript[i] does not overlap then no overlap and on same strand so return
                    break

            self.hooks["contig_built"].exec(new_contig)
            yield new_contig
        
        return
            
    """
    Attempts to merge right into left in-place.
    Returns True if merge succeeded, otherwise returns False.
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


    def merge_contig(self, contig):
        # From left to right, merge all the transcripts
        # Go again until contig.transcripts is exhausted
        # Hint: this works because contig.transcripts returns a new iterator each time

        merged = set()
        for left in contig.transcripts:
            for right in contig.transcripts:
                if left.id is not right.id and right.id not in merged and left.id not in merged:
                    if self.merge_transcripts(left, right):
                        merged.add(right.id)

        for t_id in merged:
            contig.remove_transcript_by_id(t_id)
            
        self.hooks["contig_merged"].exec(contig)

        return contig

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
                print("contig")
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
