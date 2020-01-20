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

# Overwrite file contents first
open(args.output, 'w').close()

for i, contig in enumerate(build_contigs(transcripts)):
    merged = build_merged(contig)
    write(merged, f"contig_{i}", args.output)