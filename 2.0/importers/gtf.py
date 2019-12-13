from models.exon import Exon
from models.transcript import Transcript

def parse(path):
    transcripts = {}
    with open(path, "r") as f:
        for line in f:            
            try:
                data = line.split("\t")
                transcript_id = data[8].split(";")[1].replace("transcript_id", "").replace("\"", "").replace(" ", "")
                gene_id = data[8].split(";")[0].replace("gene_id", "").replace("\"", "").replace(" ", "")
                exon = Exon("gtf", data[0], int(data[3]), int(data[4]), data[6], transcript_id, gene_id)
                prev_transcript = transcripts.get(exon.transcript_id, None)
                
                if prev_transcript:
                    prev_transcript.add_exon(exon)
                else:
                    transcripts[exon.transcript_id] = Transcript([exon])

            except Exception as e:
                # TODO: Handle this exception
                print(e)
                pass
    return list(transcripts.values())