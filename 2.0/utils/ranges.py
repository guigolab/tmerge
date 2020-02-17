from functools import partial

# Operations on two ranges
def overlaps(range1, range2):
    return range2[1] > range1[0] and range1[1] > range2[0]

def within(x, range1):
    if len(range1) == 0:
        return False
    return x > range1[0] and x < range1[1]

def lte(x, range1):
    return x <= range1[0]

def gte(x, range1):
    return x >= range1[1]

# Operations on sets of ranges
def ordered_subset(base, compare):
    if len(compare) == 0:
        return True
    if len(compare) > len(base):
        return False
    if base == compare:
        return True

    for i in range(len(base)):
        if base[i:i+len(compare)] == compare:
            return True
    return False

def _to_set_all(func):
    return lambda x, range_set: all(map(partial(func, x), range_set))

def _to_set_any(func):
    return lambda x, range_set: any(map(partial(func, x), range_set))

overlaps_all = _to_set_all(overlaps)
within_all = _to_set_all(within)
lte_all = _to_set_all(lte)
gte_all = _to_set_all(gte)

overlaps_any = _to_set_any(overlaps)
within_any = _to_set_any(within)
lte_any = _to_set_any(lte)
gte_any = _to_set_any(gte)