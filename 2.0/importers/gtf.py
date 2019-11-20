from Models.Alignment import Alignment

def parse(path):
    sequences = []
    with open(path, "r") as f:
        for line in f:            
            try:
                data = line.split("\t")
                sequence = Alignment("gtf", data[0], data[3], data[4])
                sequences.append(sequence)
            except Exception as e:
                # TODO: Handle this exception
                print(e)
                pass
    return sequences