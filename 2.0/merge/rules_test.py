import unittest
from merge.rules import transcript_overlap, no_TSS_TES_overlap, ordered_subset, ruleset
from utils.fakers import Faker
import copy

class TestRules(unittest.TestCase):
    def setUp(self):
        self.faker = Faker("gtf", "chr1")

    # Test individual rules
    # =====================
    def test_transcript_overlap(self):
        t1 = self.faker.tm(0, 0, 0, 1000)
        t2 = self.faker.tm(0, 0, 20, 1000)
        t3 = self.faker.tm(0, 0, 2000, 8000)

        self.assertTrue(transcript_overlap(t1, t2))
        self.assertTrue(transcript_overlap(t2,t1))
        self.assertFalse(transcript_overlap(t1,t3))


    def test_ordered_subset(self):
        t1 = self.faker.tm(5)
        t2 = copy.deepcopy(t1)

        # Knock off an intron from end of t1
        t1.remove_junction(*t1.junctions[len(t1.junctions) - 1])
        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

        # Knock off an intron from start of t1
        t1.remove_junction(*t1.junctions[0])
        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

        # Add exon back to end of t1
        t1 = copy.deepcopy(t2)
        t1.remove_junction(*t1.junctions[0])
        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

        # Add two exons to t1. Hence, adding two introns making them not ordered subsets
        t1.add_junction(t1.TES - 100, t1.TES - 50)
        t1.add_junction(t1.TES - 200, t1.TES - 150)

        self.assertFalse(ordered_subset(t1, t2))
        self.assertFalse(ordered_subset(t2, t1))

        # Compare two transcripts with ordered subset of introns but with very large difference in num introns
        t1 = self.faker.tm(15)
        t2 = copy.deepcopy(t1)
        for i in range(0,5):
            t2.remove_junction(*t1.junctions[i])

        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

    def test_no_TSS_TES_overlap(self):
        t1 = self.faker.tm(2)
        t2 = copy.deepcopy(t1)

        self.assertTrue(no_TSS_TES_overlap(t1, t2))

        t3 = copy.deepcopy(t1)
        t3.remove_junction(*t3.junctions[0])
        t3.TSS = t1.junctions[0][0] + 5

        self.assertFalse(no_TSS_TES_overlap(t1, t3))
        self.assertFalse(no_TSS_TES_overlap(t3, t1))

        t4 = copy.deepcopy(t1)
        t4.remove_junction(*t4.junctions[len(t4.junctions) - 1])
        t4.TES = t1.junctions[len(t1.junctions) - 1][1] - 5

        self.assertFalse(no_TSS_TES_overlap(t1, t4))
        self.assertFalse(no_TSS_TES_overlap(t4, t1))

        # test with a monoexon
        t5 = self.faker.tm(0, 0, t1.junctions[0][0] + 5, 50)

        self.assertFalse(no_TSS_TES_overlap(t1, t5))
        self.assertFalse(no_TSS_TES_overlap(t1, t5))

        # Test with two monoexons
        t6 = copy.deepcopy(t5)
        t6.TES = 200
        t6.TSS = 150

        self.assertTrue(no_TSS_TES_overlap(t6, t5))
        self.assertTrue(no_TSS_TES_overlap(t5,t6))

    # Test ruleset
    # ============
    def test_same_introns(self):
        t1 = self.faker.tm(10)
        t2 = copy.deepcopy(t1)

        self.assertTrue(ruleset(t1, t2))

        t2.add_junction(t1.junctions[2][1] + 50, t1.junctions[2][1] + 100)
        self.assertFalse(ruleset(t1, t2))

    def test_monoexonic_overlap(self):
        polyexon = self.faker.tm(5)

        # Test overlaps with last exon
        monoexon = self.faker.tm(0, start=polyexon.TES - 1)
        self.assertTrue(ruleset(polyexon, monoexon))
        self.assertTrue(ruleset(monoexon, polyexon))

        # Test overlaps with middle exon
        monoexon = self.faker.tm(0, start=polyexon.junctions[2][1] + 1, length=5)
        self.assertTrue(ruleset(polyexon, monoexon))
        self.assertTrue(ruleset(monoexon, polyexon))

        # Returns true on two monoexons that overlap
        monoexon2 = self.faker.tm(0, start=monoexon.TES - 3)
        self.assertTrue(ruleset(monoexon, monoexon2))
        self.assertTrue(ruleset(monoexon2, monoexon))
