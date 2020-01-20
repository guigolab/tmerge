def write(contigs, output_path):
    with open(output_path, "w") as f:
        for i, contig in enumerate(contigs):
            for transcript in contig.transcripts:
                    for exon in transcript.exons:
                        data = [
                            exon.chromosome,
                            f"contig{i}",
                            "exon",
                            str(exon.start),
                            str(exon.end),
                            ".",
                            exon.strand,
                            ".",
                            f"\"gene_id\" {exon.gene_id}; \"transcript_id\" {exon.transcript_id}"
                        ]

                        f.write("\t".join(data))
                        f.write("\n")