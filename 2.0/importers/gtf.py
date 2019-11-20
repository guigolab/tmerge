from Models.Alignment import Alignment

def parse(path, lines=False):
    sequences = []
    with open(path, "r") as f:
        line_count = 0
        while not lines or line_count < lines:
            try:
                line = f.readline().split("\t")
                sequence = Alignment("gtf", line[0], line[3], line[4])
                sequences.append(sequence)
            except Exception as e:
                # TODO: Handle this exception
                print(e)
                pass
            finally:
                line_count += 1
    return sequences