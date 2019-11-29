from itertools import product
from functools import reduce, partial

def _coordinates_to_tuple(t):
    return (t.start, t.end)

_objs_to_coordinate_tuples = partial(map, _coordinates_to_tuple)

def _overlap(t1, t2):
    return t1[1] > t2[0] and t2[1] > t1[0]

def _get_intron_set(t):
    return map(lambda exons: (exons[0].end, exons[1].start), zip(t.exons, t.exons[1:]))

def transcript_overlap(t1, t2):
    return _overlap(_coordinates_to_tuple(t1), _coordinates_to_tuple(t2))

def same_introns(t1, t2):
    return list(_get_intron_set(t1)) == list(_get_intron_set(t2))

def no_exon_intron_overlap(t1, t2):
    longer = t1 if len(t1.exons) > len(t2.exons) else t2
    shorter = t1 if len(t1.exons) <= len(t2.exons) else t2

    return all(
        map(
            lambda x: not _overlap(*x),
            product(
                _get_intron_set(shorter),
                _objs_to_coordinate_tuples(longer.exons)
            )
        )
    )