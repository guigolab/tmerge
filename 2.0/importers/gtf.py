def parse(path, lines=False):
    sequences = []
    with open(path, "r") as f:
        line_count = 0
        while not lines or line_count < lines:
            try:
                line = f.readline().split("\t")
                sequence = {
                    "seqname": line[0],
                    "source": line[1],
                    "feature": line[2],
                    "start": line[3],
                    "end": line[4],
                    "score": line[5],
                    "strand": line[6],
                    "frame": line[7],
                    "attributes": line[8]
                }
                sequences.append(sequence)
            except:
                # TODO: Handle this exception
                pass
            finally:
                line_count += 1
    return sequences