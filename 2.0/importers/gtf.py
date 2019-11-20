from models.exon import Exon

def parse(path):
    sequences = []
    with open(path, "r") as f:
        for line in f:            
            try:
                data = line.split("\t")
                transcript_id = data[8].split(";")[1].replace("transcript_id", "").replace("\"", "").replace(" ", "")
                gene_id = data[8].split(";")[0].replace("gene_id", "").replace("\"", "").replace(" ", "")
                sequence = Exon("gtf", data[0], data[3], data[4], data[6], transcript_id, gene_id)
                sequences.append(sequence)
            except Exception as e:
                # TODO: Handle this exception
                print(e)
                pass
    return sequences