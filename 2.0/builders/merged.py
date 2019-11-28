from models.contig import Contig
from models.transcript import Transcript

def build(contigs):
    return list(map(merge, contigs))

def merge(contig):
    contig = Contig(merge_monoexonic(contig.transcripts))
    return contig

def merge_monoexonic(transcripts, i=0):
    if i >= len(transcripts) -1:
        return transcripts
    
    if transcripts[i].exons[0].end > transcripts[i+1].exons[0].start:
        # Overlapping so merge
        transcripts[i].add_exon(transcripts[i+1].exons[0])

        # Remove i+1 from list
        transcripts.pop(i+1)
    else: 
        i=i+1

    return merge_monoexonic(transcripts, i)    