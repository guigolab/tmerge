def write(contigs, output_path):
    for i, contig in enumerate(contigs):
        with open(output_path + f"_contig{i}.gtf", "w") as f:
            for transcript in contig.transcripts:
                    for exon in transcript.exons:
                        data = [
                            "chr1",
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