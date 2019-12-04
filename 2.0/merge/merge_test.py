import unittest
from merge.merge import merge
from utils.fakers import Faker
import copy

class TestMerge(unittest.TestCase):
    """
    Test the merge function.

    Note that only test real world scenarios. Merge is dumb and given any two transcripts, it will try and merge.
    There is no checking if the two transcripts should be merged. That is done by the rules given in merge.rules
    """
    def setUp(self):
        self.faker = Faker("gtf", "chr1")

    def test_monoexonic_merge(self):
        t1 = self.faker.transcript(1, start=0)
        t2 = self.faker.transcript(1, start=t1.end - 10)

        merged = merge(t1, t2)
        self.assertEquals(merged.end, 90)

        merged2 = merge(t2 ,t1)
        self.assertEquals(merged2.end, 90)

    def test_simple_polyexonic_merge(self):
        t1 = self.faker.transcript(5, start=0)
        t2 = copy.deepcopy(t1)

        t2.pop_exon()
        merged = merge(t1, t2)
        self.assertEquals(merged.end, t1.end)
        self.assertEquals(merged.start, 0)

    def test_exonic_overhang_merge(self):
        t1 = self.faker.transcript(5, start=0)
        t2 = copy.deepcopy(t1)

        t2.pop_exon()
        t2.exons[len(t2.exons) - 2].end -= 10

        merged = merge(t1, t2)
        self.assertEquals(merged.end, t1.end)
        self.assertEquals(merged.start, 0)

    def test_missing_start_exon_merge(self):
        t1 = self.faker.transcript(5, start = 0)
        t2 = copy.deepcopy(t1)

        t2.exons.pop(0)
        t2.start = t2.exons[0].start

        merged = merge(t1, t2)
        self.assertEquals(merged.end, t1.end)
        self.assertEquals(merged.start, t1.start)

    def test_exonic_underhang_merge(self):
        t1 = self.faker.transcript(5, start=0)
        t2 = copy.deepcopy(t1)

        # Cut off start of first exon
        t2.exons[0].start += 10
        t2.start += 10

        merged = merge(t1, t2)
        self.assertEquals(merged.start, 0)
        self.assertEquals(merged.end, t1.end)