def overlap(t1, t2):
    return t1.end > t2.start

def _get_intron_set(t):
    return list(map(lambda exons: (exons[0].end, exons[1].start), zip(t.exons, t.exons[1:])))

def same_introns(t1, t2):
    return _get_intron_set(t1) == _get_intron_set(t2)