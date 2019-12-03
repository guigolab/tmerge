from itertools import product
from functools import reduce, partial

def _make_coordinates_tuple(t):
    return (t.start, t.end)

_objs_to_coordinate_tuples = partial(map, _make_coordinates_tuple)

def _overlap(t1, t2):
    return t1[1] > t2[0] and t2[1] > t1[0]

def _get_intron_set(t):
    return map(lambda exons: (exons[0].end, exons[1].start), zip(t.exons, t.exons[1:]))

def transcript_overlap(t1, t2):
    return _overlap(_make_coordinates_tuple(t1), _make_coordinates_tuple(t2))

def same_introns(t1, t2):
    return list(_get_intron_set(t1)) == list(_get_intron_set(t2))

def no_exon_intron_overlap(t1, t2):
    shorter, longer = (t1, t2) if len(t1.exons) < len(t2.exons) else (t2, t1)

    return all(
        map(
            lambda x: not _overlap(*x),
            product(
                _get_intron_set(shorter),
                _objs_to_coordinate_tuples(longer.exons)
            )
        )
    )

def monoexonic_overlap(t1, t2):
    # TODO: Ask julien. What if both are monoexonic?
    if (len(t1.exons) != 1 and len(t2.exons) > 1) or (len(t2.exons) != 1 and len(t1.exons) > 1):
        raise TypeError("At least one of the transcripts must be monoexonic")

    if len(t1.exons) == 1 and len(t2.exons) == 1:
        raise TypeError("At least one of the transcripts must be polyexonic")
    
    monoexon, polyexon = (t1, t2) if len(t1.exons) == 1 else (t2, t1)
    overlaps_monoexon = partial(_overlap, _make_coordinates_tuple(monoexon.exons[0]))

    return any(
        map(
            overlaps_monoexon,
            _objs_to_coordinate_tuples(polyexon.exons)
        )
    )

    