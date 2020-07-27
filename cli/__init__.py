import getopt, argparse
from plugins import Stats, ReadSupport, SplicedLengths, MergedInfo, SpliceSiteScoring, MinLength, Depth
from src import merge

"""
tmerge CLI
==========
Runs tmerge as a CLI with stats, read support, spliced lengths, merged info and splice site scoring plugins enabled.
"""
def main():
    description = "tmerge. Merge transcriptome read-to-genome alignments into non-redundant transcript models."
    parser = argparse.ArgumentParser(description=description)

    # Basic options passed to merger and not for plugins
    parser.add_argument("-i", "--input", help="Input GTF file")
    parser.add_argument("-o", "--output", help="Output GTF file")
    parser.add_argument("--processes", help="The number of processes (threads) allowed to run. Must be greater than 2. If left unspecified, then will use the maximum number available.", type=int, default=None)
    parser.add_argument("--tolerance", help="Maximum number of nucleotides of terminal exon overhang allowed within an intron of another transcript during the merging of input reads", type=int, default=0)
    parser.add_argument("--no_merge", action="store_true", default=False, help="Disable merging. Only pre-merge filtering is done. Useful for debugging.")

    # Stats
    parser.add_argument("-s", "--stats", action="store_true", help="Provide statistics for merged transcripts.")
    parser.add_argument("--stats_output", help="Path to the output stats. Set to '-' for stderr.", default="tmerge.log")

    # Options for MinLength
    parser.add_argument("--min_length", help="The minimum length required for output transcripts.", default=0, type=int)

    # Options for Depth
    # Depth disabled
    # parser.add_argument("--trim", action="store_true", help="(experimental) Trim output transcript ends? Ends are trimmed when there is a sudden drop in depth.", default=False)

    # Options for ReadSupport
    # TODO: Change --min_isoform_fraction/--min_read_support to one option
    # If it is < 1 it is min_isoform_fraction if it is > 1 it's min_read_support
    parser.add_argument("--min_isoform_fraction", help="Minimum number of times a read alignment (by exon/intron structure) needs to be present in the input expressed as a fraction of the maximum value at a gene locus.", type=float, default=0)
    parser.add_argument("--min_read_support", help="Minimum number of times a read alignment (by exon/intron structure) needs to be present in the input.", type=int, default=1)
    parser.add_argument("--end_fuzz", help="Tolerated fuzziness of 5' and 3' ends for two reads to be considered equivalent when calculating read support", type=int, default=0)
    
    # Options for SpliceSiteScoring
    parser.add_argument("--splice_scoring", action="store_true", help="Enable splice site scoring.", default=False)
    parser.add_argument("--acceptor_path", help="Only if splice_scoring enabled. Path to the acceptor file.")
    parser.add_argument("--donor_path", help="Only if splice_scoring enabled. Path to the donor file.")
    parser.add_argument("--fasta_path", help="Only if splice_scoring enabled. Path to the FASTA genome file.")
    parser.add_argument("--valid_acceptor", type=int, default=-24, help="Only if splice_scoring enabled. Threshold at which transcripts are removed.")
    parser.add_argument("--valid_donor", type=int, default=-19, help="Only if splice_scoring enabled. Threshold at which transcripts are removed.")

    # Flags for ONT/PacBio
    parser.add_argument("--ont", action="store_true", default=False, help="Use pre-defined optimum options for ONT datasets.")
    parser.add_argument("--pacbio", action="store_true", default=False, help="Use pre-defined optimum options for PacBio datasets.")

    args = parser.parse_args()

    if args.ont:
        args.min_length = 200
        args.min_isoform_fraction = 0.08
        args.end_fuzz = 2
        args.tolerance = 4
    elif args.pacbio:
        args.min_length = 200
        args.end_fuzz = 4
        args.tolerance = 2
        args.min_isoform_fraction = 0.02

    if (args.ont or args.pacbio) and args.splice_scoring:
        args.min_isoform_fraction = 0
        args.min_read_support = 1


    if args.splice_scoring and (args.min_isoform_fraction > 0 or args.min_read_support > 1):
        raise TypeError("You cannot use splice scoring and min_isoform_fraction/min_read_support at the same time")

    if args.ont and args.pacbio:
        raise TypeError("ONT and PacBio flags cannot be used together")

    """
    Load plugins
    ============
    Order matters here for speed. 
    Functions that tap into hooks are ran in the order they are tapped. I.e. the first one to tap, will be the first to execute.
    If plugins are removing transcripts than they should run early to reduce the search space.
    """
    plugins = [
        SpliceSiteScoring,
        ReadSupport,
        SplicedLengths,
        MergedInfo,
        MinLength,
        # Disabled depth as is experimental
        # Depth,
        Stats
    ]

    """
    Begin merging
    =============
    """
    merge(
        input_path=args.input,
        output_path=args.output,
        plugins=plugins,
        **vars(args)
    )

if __name__ == "main":
    main()