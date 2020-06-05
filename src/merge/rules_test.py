import unittest
import copy

from .rules import transcript_overlap, TSS_TES_overlap, ordered_subset, ruleset, first_last_exon_intron_overlap
from ..utils.fakers import Faker

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
        t1.remove_junction(*t1.sorted_junctions[len(t1.sorted_junctions) - 1])
        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

        # Knock off an intron from start of t1
        t1.remove_junction(*t1.sorted_junctions[0])
        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

        # Add exon back to end of t1
        t1 = copy.deepcopy(t2)
        t1.remove_junction(*t1.sorted_junctions[0])
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
            t2.remove_junction(*t1.sorted_junctions[i])

        self.assertTrue(ordered_subset(t1, t2))
        self.assertTrue(ordered_subset(t2, t1))

    def test_no_TSS_TES_overlap(self):
        t1 = self.faker.tm(2)
        t2 = copy.deepcopy(t1)

        self.assertFalse(TSS_TES_overlap(t1, t2))

        t3 = copy.deepcopy(t1)
        t3.remove_junction(*t3.sorted_junctions[0])
        t3.TSS = t1.sorted_junctions[0][0] + 5

        self.assertTrue(TSS_TES_overlap(t1, t3))
        self.assertTrue(TSS_TES_overlap(t3, t1))

        t4 = copy.deepcopy(t1)
        t4.remove_junction(*t4.sorted_junctions[len(t4.sorted_junctions) - 1])
        t4.TES = t1.sorted_junctions[len(t1.sorted_junctions) - 1][1] - 5

        self.assertTrue(TSS_TES_overlap(t1, t4))
        self.assertTrue(TSS_TES_overlap(t4, t1))

        # test with a monoexon
        t5 = self.faker.tm(0, 0, t1.sorted_junctions[0][0] + 5, 50)

        self.assertTrue(TSS_TES_overlap(t1, t5))
        self.assertTrue(TSS_TES_overlap(t1, t5))

        # Test with two monoexons
        t6 = copy.deepcopy(t5)
        t6.TES = 200
        t6.TSS = 150

        self.assertFalse(TSS_TES_overlap(t6, t5))
        self.assertFalse(TSS_TES_overlap(t5,t6))

        # Test boundary case where TES has same coordinate as junction
        t7 = self.faker.tm(0, 0, 0, t1.sorted_junctions[1][1])
        self.assertTrue(TSS_TES_overlap(t1, t7))

        # Test boundary case where TSS has same coordinate as junction
        t8 = self.faker.tm(0, 0, t1.sorted_junctions[1][0], 50)
        self.assertTrue(TSS_TES_overlap(t1, t8))

    # Test ruleset
    # ============
    def test_same_introns(self):
        t1 = self.faker.tm(4)
        t2 = copy.deepcopy(t1)

        self.assertTrue(ruleset(t1, t2))

        t2.add_junction(t1.sorted_junctions[2][1] + 1, t1.sorted_junctions[2][1] + 1)
        self.assertFalse(ruleset(t1, t2))
        self.assertFalse(ruleset(t2, t1))

    def test_monoexonic_overlap(self):
        polyexon = self.faker.tm(5)

        # Test overlaps with last exon
        monoexon = self.faker.tm(0, start=polyexon.TES - 1)
        self.assertTrue(ruleset(polyexon, monoexon))
        self.assertTrue(ruleset(monoexon, polyexon))

        # Test overlaps with middle exon
        monoexon = self.faker.tm(0, start=polyexon.sorted_junctions[2][1] + 1, length=5)
        self.assertTrue(ruleset(polyexon, monoexon))
        self.assertTrue(ruleset(monoexon, polyexon))

        # Returns true on two monoexons that overlap
        monoexon2 = self.faker.tm(0, start=monoexon.TES - 3)
        self.assertTrue(ruleset(monoexon, monoexon2))
        self.assertTrue(ruleset(monoexon2, monoexon))

    def test_last_exon_overlap(self):
        # Since only considering introns (sorted_junctions), it is possible for one transcript's junction chain to be an ordered subset of anothers
        # but the last exon to overlap many introns of the others and still have the same TES
        # so:
        # ====---====----====----====---====
        # ====---====----===================
        t1 = self.faker.tm(10)
        t2 = copy.deepcopy(t1)

        for i in range(3):
            t2.remove_junction(*t2.sorted_junctions[len(t2.sorted_junctions) - 1])

        self.assertFalse(ruleset(t1, t2))

    def test_first_exon_overlap(self):
        # same as test_last_exon_overlap but for the first exon
        # so:
        # ====---====----====----====---====
        # ===================----====---====
        t1 = self.faker.tm(10)
        t2 = copy.deepcopy(t1)

        for i in range(3):
            t2.remove_junction(*t2.sorted_junctions[0])

        self.assertFalse(ruleset(t1, t2))

    def test_tolerance(self):
        t1 = self.faker.tm(3)
        t2 = copy.deepcopy(t1)

        for junction in t2.sorted_junctions:
            t2.add_junction(junction[0] - 1, junction[1] + 1)
            t2.remove_junction(*junction)
        
        self.assertFalse(ruleset(t1, t2, 0))
        self.assertFalse(ruleset(t1, t2, 5))

        # Should merge when one of the terminal exons overhangs
        # ======----========----====
        #            =======---====
        t3 = self.faker.tm(0, 50, t1.sorted_junctions[0][1] + 5)
        t3.add_junction(t1.sorted_junctions[1][0], t1.sorted_junctions[1][1])
        t3.add_junction(*t1.sorted_junctions[2])
        self.assertTrue(ruleset(t1,t3,0))
        self.assertTrue(ruleset(t1,t3,5))

        # Should merge when one of the terminal exons underhangs
        # ======----========----====
        #          =========----====
        t4 = self.faker.tm(0, 50, t1.sorted_junctions[1][1] - 5)
        t4.add_junction(t1.sorted_junctions[1][0], t1.sorted_junctions[1][1])
        t4.add_junction(*t1.sorted_junctions[2])
        self.assertFalse(ruleset(t1,t4,0))
        self.assertTrue(ruleset(t1,t4,5))
        

        # Test monoexon merges
        # E.g. should merge when tolerance set to > zero value:
        # =====-----=====-----=====-----====
        #          ======
        t4 = self.faker.tm(0, 0, t1.sorted_junctions[1][1] - 5, 20)
        self.assertTrue(ruleset(t1, t4, 5))
        self.assertFalse(ruleset(t1, t4, 0))
        # =====-----=====-----=====-----====
        #             ====
        t5 = self.faker.tm(0, 0, t1.sorted_junctions[1][0] - 15, 20)
        self.assertTrue(ruleset(t1, t5, 5))
        self.assertFalse(ruleset(t1, t5, 0))

        # Test merge with TSS "staircase"
        t6 = self.faker.tm(4)
        t7 = copy.deepcopy(t6)

        t7.TSS = t6.TSS + 20
        self.assertTrue(ruleset(t6,t7, 5))