#!/usr/bin/env python3

from merge.merge import Merge
from plugins.stats import Stats
from plugins.read_support import ReadSupport
from plugins.spliced_lengths import SplicedLengths
from plugins.merged_info import MergedInfo
from plugins.splice_site_scoring import SpliceSiteScoring
import getopt, argparse

unix_options = "ho:i:s"
gnu_options = ["help", "output=", "input=", "stats"]

description = "tmerge 2.0 Beta"
parser = argparse.ArgumentParser(description=description)

parser.add_argument("-i", "--input", help="Input GTF file")
parser.add_argument("-o", "--output", help="Output GTF file")
parser.add_argument("-s", "--stats", action="store_true", help="Provide statistics for merged transcripts.")
parser.add_argument("--tolerance", help="Exon overhang tolerance", type=int, default=0)
parser.add_argument("--support", help="Minimum read support", type=int, default=0)
parser.add_argument("--fuzz", help="end fuzz", type=int, default=0)
parser.add_argument("--speed", help="Speed mode. Enables options that forgoe sensitivity and precision for faster merge time.", action="store_true")
parser.add_argument("--processes", help="The number of processes (threads) allowed to run. Must be greater than 2. If left unspecified, then will use the maximum number available.", type=int, default=None)
parser.add_argument("--spliceScoring", action="store_true")
parser.add_argument("--acceptor")
parser.add_argument("--donor")
parser.add_argument("--fasta")
parser.add_argument("--validAcceptor", type=int, default=4)
parser.add_argument("--validDonor", type=int, default=4)


args = parser.parse_args()

merger = Merge(args.input, args.output, args.tolerance, args.processes, args.speed)

"""
Load plugins
============
Order matters here for speed. 
Functions that tap into hooks are ran in the order they are tapped. I.e. the first one to tap, will be the first to execute.
If plugins are removing transcripts than they should run early to reduce the search space.
"""
if args.spliceScoring:
    SpliceSiteScoring(merger, args.donor, args.acceptor, args.fasta, args.validDonor, args.validAcceptor)

ReadSupport(merger, args.fuzz, args.support, args.speed)
SplicedLengths(merger, args.speed)
MergedInfo(merger, args.speed)

if args.stats:
    Stats(merger)

"""
Begin merging
=============
"""
merger.merge()