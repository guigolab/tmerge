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

    # Stats
    parser.add_argument("-s", "--stats", action="store_true", help="Provide statistics for merged transcripts.")

    # Options for MinLength
    parser.add_argument("--min_length", help="The minimum length required for output transcripts.", default=200, type=int)

    # Options for Depth
    parser.add_argument("--trim", help="Trim output transcript ends? Ends are trimmed when there is a sudden drop in depth.", default=True)

    # Options for ReadSupport
    parser.add_argument("--min_abundance", help="Minimum abundance for transcripts", type=float, default=0.15)
    parser.add_argument("--min_read_support", help="Minimum number of times a read alignment (by exon/intron structure) needs to be present in the input.", type=int, default=1)
    parser.add_argument("--end_fuzz", help="Tolerated fuzziness of 5' and 3' ends for two reads to be considered equivalent when calculating read support", type=int, default=0)
    
    # Options for SpliceSiteScoring
    parser.add_argument("--splice_scoring", action="store_true", help="Enable splice site scoring.", default=False)
    parser.add_argument("--acceptor_path", help="Only if splice_scoring enabled. Path to the acceptor file.")
    parser.add_argument("--donor_path", help="Only if splice_scoring enabled. Path to the donor file.")
    parser.add_argument("--fasta_path", help="Only if splice_scoring enabled. Path to the FASTA genome file.")
    parser.add_argument("--valid_acceptor", type=int, default=4, help="Only if splice_scoring enabled. Threshold at which transcripts are removed.")
    parser.add_argument("--valid_donor", type=int, default=4, help="Only if splice_scoring enabled. Threshold at which transcripts are removed.")


    args = parser.parse_args()

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