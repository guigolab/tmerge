import unittest

from . import ranges

class TestRanges(unittest.TestCase):
    def test_overlaps(self):
        self.assertFalse(ranges.overlaps((0, 5), (8, 10)))
        self.assertFalse(ranges.overlaps((0,5), (5, 19)))
        self.assertFalse(ranges.overlaps((-5,0),(-10,-5)))
        self.assertTrue(ranges.overlaps((5,10),(9,15)))
        self.assertTrue(ranges.overlaps((5,10), (5,10)))
        self.assertTrue(ranges.overlaps((3,10),(0,4)))
        self.assertFalse(ranges.overlaps((0,5), (-5, 0)))
        self.assertTrue(ranges.overlaps((5,10),(11,15), 1))
        self.assertTrue(ranges.overlaps((5,10),(15,20), 5))
        self.assertFalse(ranges.overlaps((5,10),(11,15), 0))

    def test_overlaps_any(self):
        range1 = [
            (5,10),
            (15,20),
            (25,30),
            (30,35)
        ]
        self.assertTrue(ranges.overlaps_any((10,30), range1))
        self.assertFalse(ranges.overlaps_any((0,2), range1))
        self.assertTrue(ranges.overlaps_any((0,4), range1, 2))
        self.assertFalse(ranges.overlaps_any((0,3), range1, 1))
        self.assertTrue(ranges.overlaps_any((0, 399), [(400, 450), (550, 600), (650, 700)], 2))

    def test_overlaps_all(self):
        range1 = [
            (5,10),
            (15,20),
            (25,30),
            (30,35)
        ]
        self.assertTrue(ranges.overlaps_all((0,50), range1))
        self.assertFalse(ranges.overlaps_all((10,30), range1))
        self.assertTrue(ranges.overlaps_all((0,3), range1, 10000))

    def test_within(self):
        r = (10,20)
        self.assertTrue(ranges.within(15, r))
        self.assertFalse(ranges.within(10, r))
        self.assertFalse(ranges.within(20,r))
        self.assertFalse(ranges.within(9,r))
        self.assertFalse(ranges.within(21,r))
        self.assertFalse(ranges.within(5,r,5))
        self.assertFalse(ranges.within(25,r,5))
        self.assertTrue(ranges.within(5,r,6))
        self.assertTrue(ranges.within(25,r,6))

        # Test inclusive
        self.assertTrue(ranges.within(10, r, exclusive=False))
        self.assertTrue(ranges.within(25,r,5, exclusive=False))

    def test_within_any(self):
        r = [
            (5,10),
            (15,20),
            (25,30),
            (30,35)
        ]
        self.assertTrue(ranges.within_any(6, r))
        self.assertFalse(ranges.within_any(40, r))
        self.assertTrue(ranges.within_any(37,r,3))

    def test_ordered_subset(self):
        r1 = [
            (5,10),
            (15,20),
            (25,30),
            (35,40),
            (45,50)
        ]

        r2 = [
            (55, 60),
            (65, 70),
            (75, 80)
        ]

        self.assertTrue(ranges.ordered_subset(r1,r1[2:]))
        # Test is unidirectional. I.e. r1 is subset of r2 but r2 not subset of r1 
        self.assertFalse(ranges.ordered_subset(r1[2:], r1))

        self.assertTrue(ranges.ordered_subset(r1, r1))
        self.assertFalse(ranges.ordered_subset(r1, r2))
        self.assertFalse(ranges.ordered_subset(r2,r1))

        # Test tolerance
        # Tolerance should be tolerated from either end of the range and plus or minus from either end 
        self.assertTrue(ranges.ordered_subset(
            r1,
            list(map(lambda x: (x[0] - 1, x[1]), r1)),
            2
        ))
        self.assertTrue(ranges.ordered_subset(
            r1,
            list(map(lambda x: (x[0] + 1, x[1]), r1)),
            2
        ))
        self.assertTrue(ranges.ordered_subset(
            r1,
            list(map(lambda x: (x[0], x[1] + 1), r1)),
            2
        ))
        self.assertTrue(ranges.ordered_subset(
            r1,
            list(map(lambda x: (x[0], x[1] - 1), r1)),
            2
        ))
        self.assertTrue(ranges.ordered_subset(
            r1,
            list(map(lambda x: (x[0] + 1, x[1] + 1), r1)),
            2
        ))
        self.assertTrue(ranges.ordered_subset(
            r1,
            list(map(lambda x: (x[0] - 1, x[1] - 1), r1)),
            2
        ))
        # Test exclusivity
        self.assertTrue(ranges.ordered_subset(
            r1,
            list(map(lambda x: (x[0] - 0, x[1]), r1)),
            1
        ))
        # Test when only one of ranges is altered
        r3 = r1[::]
        r3[0] = (r3[0][0] - 5, r3[0][1] + 5)
        self.assertTrue(ranges.ordered_subset(r1, r3, 6))