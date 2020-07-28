from functools import partial

# Exclusive
def overlaps(range1, range2, tolerance = 0):
    return (range2[1] + tolerance) > (range1[0] - tolerance) and (range1[1] + tolerance) > (range2[0] - tolerance)

# Exclusive
def within(x, range1, tolerance = 0, exclusive = True):
    if len(range1) == 0:
        return False
    if exclusive:
        return x > range1[0] - tolerance and x < range1[1] + tolerance
    else:
        return x >= range1[0] - tolerance and x <= range1[1] + tolerance


def lte(x, range1, tolerance = 0):
    return x <= range1[0] + tolerance

def gte(x, range1, tolerance = 0):
    return x >= range1[1] - tolerance

def find_overlapping(range1, range_set, tolerance = 0):
    return [x for x in range_set if overlaps(range1, x, tolerance)]

# Operations on sets of ranges
def ordered_subset(base, compare, tolerance = 0):
    if len(compare) == 0:
        return True
    if len(compare) > len(base):
        return False
    if base == compare:
        return True

    for i in range(len(base)):
        if base[i:i+len(compare)] == compare:
            return True

        if not tolerance:
            continue
            
        with_tolerance = True
        for compare_i in range(len(compare)):
            tolerated_range_start = (compare[compare_i][0] - tolerance, compare[compare_i][0] + tolerance)
            tolerated_range_end = (compare[compare_i][1] - tolerance, compare[compare_i][1] + tolerance)
            base_i = i + compare_i
            if base_i > len(base) - 1:
                break
            if not (within(base[base_i][0], tolerated_range_start) and within(base[base_i][1], tolerated_range_end)):
                with_tolerance = False
        if with_tolerance:
            return True
        
    return False

def _to_set_all(func):
    return lambda x, range_set, tolerance = 0, *args, **kwargs: all(map(partial(func, x, tolerance=tolerance, *args, **kwargs), range_set))

def _to_set_any(func):
    return lambda x, range_set, tolerance = 0, *args, **kwargs: any(map(partial(func, x, tolerance=tolerance, *args, **kwargs), range_set))

overlaps_all = _to_set_all(overlaps)
within_all = _to_set_all(within)
lte_all = _to_set_all(lte)
gte_all = _to_set_all(gte)

overlaps_any = _to_set_any(overlaps)
within_any = _to_set_any(within)
lte_any = _to_set_any(lte)
gte_any = _to_set_any(gte)