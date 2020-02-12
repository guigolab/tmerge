from functools import partial

# Operations on two ranges
def overlaps(range1, range2):
    return range2[1] > range1[0] and range1[1] > range2[0]

def within(x, range1):
    if len(range1) == 0:
        return False
    return x >= range1[0] and x <= range1[1]

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

def _to_set(func):
    return lambda x, range_set: all(map(partial(func, x), range_set))

overlaps_set = _to_set(overlaps)
within_set = _to_set(within)
lte_set = _to_set(lte)
gte_set = _to_set(gte)