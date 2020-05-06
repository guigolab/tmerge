#!/usr/bin/env python3

from merge.merge import Merge
from plugins.stats import Stats
from plugins.read_support import ReadSupport
import getopt, argparse

unix_options = "ho:i:s"
gnu_options = ["help", "output=", "input=", "stats"]

description = "tmerge 2.0 Beta"
parser = argparse.ArgumentParser(description=description)

parser.add_argument("-i", "--input", help="Input GTF file")
parser.add_argument("-o", "--output", help="Output GTF file")
parser.add_argument("-s", "--stats", action="store_true", help="Provide statistics for merged transcripts.")
parser.add_argument("-t", "--tolerance", help="Exon overhang tolerance", type=int, default=0)
parser.add_argument("-f", "--support", help="Minimum read support", type=int, default=0)
parser.add_argument("-e", "--fuzz", help="end fuzz", type=int, default=0)
parser.add_argument("-p", "--processes", help="Processes", type=int, default=None)


args = parser.parse_args()

merger = Merge(args.input, args.output, args.tolerance, args.processes)

ReadSupport(merger, args.fuzz, args.support)

if args.stats:
    Stats(merger)

merger.merge()