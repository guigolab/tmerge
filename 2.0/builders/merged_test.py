from utils.fakers import Faker
from builders.merged import build
import unittest

class MergedTest(unittest.TestCase):
    def setUp(self):
        self.faker = Faker("gtf", "chr1")

    def test_overlapping_monoexonic(self):
        contig = self.faker.contig_overlapping_monoexonic()
        merged = build([contig])
        self.assertEqual(len(merged), 1)
        self.assertEqual(len(merged[0].transcripts), 1)
        self.assertEqual(contig.start, merged[0].start)
        self.assertEqual(contig.end, merged[0].end)
        self.assertEqual(len(contig.transcripts[0].exons), len(merged[0].transcripts[0].exons))