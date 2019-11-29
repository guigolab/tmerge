import unittest
from merge_rules.rules import transcript_overlap, same_introns, no_exon_intron_overlap
from utils.fakers import Faker
import copy

class TestRules(unittest.TestCase):
    def setUp(self):
        self.faker = Faker("gtf", "chr1")

    def test_transcript_overlap(self):
        t1 = self.faker.transcript(1, None, start=0, length=50)
        t2 = self.faker.transcript(1, None, start=49, length=50)

        self.assertTrue(transcript_overlap(t1, t1))

        t3 = self.faker.transcript(1, None, start=200)
        self.assertFalse(transcript_overlap(t1, t3))

    def test_same_introns(self):
        t1 = self.faker.transcript(10)
        t2 = copy.deepcopy(t1)

        self.assertTrue(same_introns(t1, t2))

        t2.add_exon(self.faker.exon())
        self.assertFalse(same_introns(t1, t2))

    def test_no_exon_intron_overlap(self):
        t1 = self.faker.transcript(2)
        t2 = copy.deepcopy(t1)

        self.assertTrue(no_exon_intron_overlap(t1, t2))

        t3 = copy.deepcopy(t1)
        t3.exons[0].end += 50
        t3.exons[1].start += 50
        t3.exons[1].end += 50

        self.assertFalse(no_exon_intron_overlap(t1, t3))

        self.assertFalse(no_exon_intron_overlap(t3, t1))