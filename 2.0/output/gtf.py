"""
Reads from a queue of contigs
"""

def write(queue, output_path):
    with open(output_path, "w") as f:
        while True:
            msg = queue.get()

            if msg == "KILL":
                return

            transcripts = msg.transcripts

            for transcript in transcripts:
                if transcript.monoexonic:
                    start = transcript.TSS
                    end = transcript.TES

                    data = [
                        transcript.chromosome,
                        "TMERGE",
                        "exon",
                        str(start),
                        str(end),
                        "0",
                        transcript.strand,
                        ".",
                        f"gene_id \"{transcript.id}\"; transcript_id \"{transcript.id}\"; MERGED_COUNT \"{transcript.transcript_count}\"; FULL_LENGTH_SUPPORT_COUNT \"{transcript.full_length_count}\";"
                    ]

                    f.write("\t".join(data))
                    f.write("\n")
                    f.flush()

                else:
                    for idx in range(0, len(transcript.junctions) + 1):
                        if idx == 0:
                            start = transcript.TSS
                            end = transcript.junctions[0][0]
                        elif idx >= len(transcript.junctions):
                            start = transcript.junctions[idx - 1][1]
                            end = transcript.TES
                        else:
                            start = transcript.junctions[idx - 1][1]
                            end = transcript.junctions[idx][0]

                        data = [
                            transcript.chromosome,
                            "TMERGE",
                            "exon",
                            str(start),
                            str(end),
                            "0",
                            transcript.strand,
                            ".",
                            f"gene_id \"{transcript.id}\"; transcript_id \"{transcript.id}\"; MERGED_COUNT \"{transcript.transcript_count}\"; FULL_LENGTH_SUPPORT_COUNT \"{transcript.full_length_count}\";"
                        ]

                        f.write("\t".join(data))
                        f.write("\n")
                        f.flush()