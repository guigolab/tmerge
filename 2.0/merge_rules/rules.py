def overlap(t1, t2):
    return t1.end > t2.start

def _get_intron_set(t):
    introns = []
    for idx, item in enumerate(t.exons):
        if idx >= len(t.exons) - 1:
            break
        introns.append((item.end, t.exons[idx+1].start))
    return introns

def same_introns(t1, t2):
    return _get_intron_set(t1) == _get_intron_set(t2)