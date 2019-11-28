import unittest
from merge_rules.rules import overlap, same_introns
from utils.fakers import Faker
import copy

class TestRules(unittest.TestCase):
    def setUp(self):
        self.faker = Faker("gtf", "chr1")

    def test_overlap(self):
        t1 = self.faker.transcript(1, None, start=0, length=50)
        t2 = self.faker.transcript(1, None, start=49, length=50)

        self.assertTrue(overlap(t1, t1))

        t3 = self.faker.transcript(1, None, start=200)
        self.assertFalse(overlap(t1, t3))

    def test_same_introns(self):
        t1 = self.faker.transcript(10)
        t2 = copy.deepcopy(t1)

        self.assertTrue(same_introns(t1, t2))

        t2.add_exon(self.faker.exon())
        self.assertFalse(same_introns(t1, t2))