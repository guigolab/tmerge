import unittest
from rules.rules import transcript_overlap, same_introns, no_exon_intron_overlap, monoexonic_overlap, ordered_subset
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

        # test with a monoexon
        t4 = copy.deepcopy(t3)
        t4.exons = [t4.exons[0]]
        t4.exons[0].end += 50

        self.assertFalse(no_exon_intron_overlap(t3, t4))
        self.assertFalse(no_exon_intron_overlap(t4, t3))

        # Test with two monoexons
        t5 = copy.deepcopy(t4)
        t5.exons[0].end += 50

        self.assertFalse(no_exon_intron_overlap(t4, t5))
        self.assertFalse(no_exon_intron_overlap(t5,t4))

    def test_monoexonic_overlap(self):
        polyexon = self.faker.transcript(5)

        # Test overlaps with last exon
        monoexon = self.faker.transcript(1, start=polyexon.end - 1)
        self.assertTrue(monoexonic_overlap(polyexon, monoexon))
        self.assertTrue(monoexonic_overlap(monoexon, polyexon))

        # Test overlaps with middle exon
        monoexon = self.faker.transcript(1, start=polyexon.exons[3].end - 1)
        self.assertTrue(monoexonic_overlap(polyexon, monoexon))
        self.assertTrue(monoexonic_overlap(monoexon, polyexon))

        # Returns true on two monoexons that overlap
        monoexon2 = self.faker.transcript(1, start=monoexon.end - 3)
        self.assertTrue(monoexonic_overlap(monoexon, monoexon2))
        self.assertTrue(monoexonic_overlap(monoexon2, monoexon))
        
        # Test fails on no monoexon
        self.assertFalse(monoexonic_overlap(polyexon, copy.deepcopy(polyexon)))

    def test_ordered_subset(self):
        t1 = self.faker.transcript(5)
        t2 = copy.deepcopy(t1)

        # Knock off an exon from t1
        t1.exons = t1.exons[:-1]
        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

        # Knock off an exon from start of t1
        t1.exons = t1.exons[1:]
        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

        # Add exon back to end of t1
        t1 = copy.deepcopy(t2)
        t1.exons = t1.exons[1:]
        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

        # Add two exons to t1. Hence, adding two introns making them not ordered subsets
        t1.add_exon(self.faker.exon(start=t1.end + 100))
        t1.add_exon(self.faker.exon(start=t1.end + 100))

        self.assertFalse(ordered_subset(t1, t2))
        self.assertFalse(ordered_subset(t2, t1))

        # Compare two transcripts with ordered subset of introns but with very large difference in exon length
        t1 = self.faker.transcript(15)
        t2 = copy.deepcopy(t1)
        t2.exons = t2.exons[5:]

        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))
        

