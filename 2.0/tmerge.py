#!/usr/bin/env python3

from merge.merge import Merge
from plugins.stats import Stats
import getopt, argparse

unix_options = "ho:i:s"
gnu_options = ["help", "output=", "input=", "stats"]

description = "tmerge 2.0 Beta"
parser = argparse.ArgumentParser(description=description)

parser.add_argument("-i", "--input", help="Input GTF file")
parser.add_argument("-o", "--output", help="Output GTF file")
parser.add_argument("-s", "--stats", action="store_true", help="Provide statistics for merged transcripts.")

args = parser.parse_args()

merger = Merge(args.input, args.output)

if args.stats:
    Stats(merger)

merger.merge()