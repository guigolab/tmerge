# NAME

tmerge

# SYNOPSIS

A utility to merge transcriptome read-to-genome alignments into non-redundant transcript models.

`tmerge` compares transcript structures (or read-to-genome alignments) present in the input and attempts to reduce transcript redundancy, _i.e._, merge compatible input transcripts into non-redundant transcript models.

See DESCRIPTION below for more details.

**Usage example**:

`tmerge -cpu <number of CPUs> --tmPrefix <string to prepend to output transcript_ids> <input GTF file> > <output file>`

## INPUT

GTF file of read-to-genome alignments, sorted by chromosome and start position.
Only `exon` records are considered.
Read alignments need to be uniquely identified with the `transcript_id` GTF attribute. `transcript_id` is the only mandatory GTF attribute in input records.

## OPTIONS

- **cpu** (integer) = Number of CPUs to use

    Default: 1

- **tmPrefix** (string) = Prefix string for transcript\_id identifiers in the output

    Default: '' (empty string)

## OUTPUT

`tmerge` outputs GTF. The output contains a set of transcript models
chr7    tmerge  exon    98207754        98209634        .       +       .       gene\_id "blahTM\_000000000007"; transcript\_id "blahTM\_000000000007"; contains "23aa004d-2eab-42ea-b731-3ebe13820c2e/436f4b3f50a7179009afed1bea7c0a631b23f514,3fa67222-f04e-4932-a69c-7585a142ddb8/436f4b3f50a7179009afed1bea7c0a631b23f514,900811a4-ec46-4bf7-b8ea-5d5fbead752e/436f4b3f50a7179009afed1bea7c0a631b23f514,d086ce93-926e-4ad5-90c8-a8a1844b1832/436f4b3f50a7179009afed1bea7c0a631b23f514,dbb204c8-5b77-4b9c-a37d-6334e4db58f6/436f4b3f50a7179009afed1bea7c0a631b23f514,f7b7c67e-c65e-4014-a4c7-8f784541fd1b/436f4b3f50a7179009afed1bea7c0a631b23f514";

# DESCRIPTION

# DEPENDENCIES

CPAN:

- UUID::Generator::PurePerl
- Parallel::ForkManager

# AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com
