from importers import gtf
from builders.contigs import build as build_contigs
from builders.merged import build as build_merged
from output.gtf import write
from utils.fakers import Faker
from functools import reduce
from models.contig import Contig
import getopt, argparse

unix_options = "ho:i:"
gnu_options = ["help", "output=", "input="]

description = "Tmerge 2.0 Beta"
parser = argparse.ArgumentParser(description=description)

parser.add_argument("-i", "--input", help="Input GTF file")
parser.add_argument("-o", "--output", help="Output GTF file")
args = parser.parse_args()

transcripts = gtf.parse(args.input)

print(f"Number of transcripts {len(transcripts)}")
contigs = build_contigs(transcripts)
print(f"Number of contigs {len(contigs)}")
merged = build_merged(contigs)
print(f"Number of merged contigs: {len(merged)}")
merged_transcripts = [c.transcripts for c in merged]
num_merged = sum([len(t) for t in merged_transcripts])
print(f"Number of transcripts in total: {num_merged} ")

write(merged, args.output)