from functools import reduce

def write_exon(f, tm_id, transcript, start, end):
    meta_info = reduce(lambda x,y: f"{x} {y}", map(lambda i: f"{i[0]} \"{i[1]}\";", transcript.meta.items())) if transcript.meta else ""

    data = [
        transcript.chromosome,
        "TMERGE",
        "exon",
        str(start),
        str(end),
        "0",
        transcript.strand,
        ".",
        f"gene_id \"{tm_id}\"; transcript_id \"{tm_id}\"; {meta_info}"
    ]

    f.write("\t".join(data))
    f.write("\n")
    f.flush()

def write(q, output_path):
    """
    Read from a queue of transcripts and write the transcript to the output GTF.
    """
    with open(output_path, "w") as f:
        id_count = 0
        while True:
            transcript = q.get()

            if transcript.monoexonic:
                start = transcript.TSS
                end = transcript.TES
                write_exon(f, f"TM_{id_count}", transcript, start, end)

            else:
                t_junctions = transcript.sorted_junctions
                for idx in range(0, len(t_junctions) + 1):
                    if idx == 0:
                        start = transcript.TSS
                        end = t_junctions[0][0]
                    elif idx >= len(transcript.junctions):
                        start = t_junctions[idx - 1][1]
                        end = transcript.TES
                    else:
                        start = t_junctions[idx - 1][1]
                        end = t_junctions[idx][0]

                    write_exon(f, f"TM_{id_count}", transcript, start, end)
            
            id_count += 1

            q.task_done()