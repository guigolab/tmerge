from models.transcript import Transcript

def merge(t1, t2):
    first, last = (t1, t2) if t1.start <= t2.start else (t2, t1)
    shorter, longer = (first, last) if first.end >= last.end else (first, last)
    
    if first is longer:
        return first

    # Get the exons that are overhanging.
    # Also gets an extra exon in case the shorter transcript is cut mid-exon
    exon_overhang = longer.exons[len(longer.exons) - len(shorter.exons) - 2:]
    
    # Remove last exon
    first.exons.pop()
    first.exons.extend(exon_overhang)
    return Transcript(first.exons)