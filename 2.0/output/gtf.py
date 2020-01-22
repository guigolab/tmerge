def write(contig, contig_name, output_path):
    with open(output_path, "a") as f:
        for transcript in contig.transcripts:
            for exon in transcript.exons:
                data = [
                    exon.chromosome,
                    contig_name,
                    "exon",
                    str(exon.start),
                    str(exon.end),
                    "0",
                    exon.strand,
                    ".",
                    f"\gene_id \"{exon.gene_id}\"; transcript_id \"{exon.transcript_id}\";"
                ]

                f.write("\t".join(data))
                f.write("\n")