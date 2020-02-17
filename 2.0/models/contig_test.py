import unittest
from utils.fakers import Faker
class TestContigModel(unittest.TestCase):
    def setUp(self):
        self.faker = Faker("gtf", "chr1")
        self.contig = self.faker.contig(10)

    def test_creates_data(self):
        self.assertEqual(len(self.contig.transcripts), 10)

    def test_adds_transcript(self):
        tm = self.faker.tm(5, 50, 10)
        self.contig.add_transcript(tm)

        self.assertEqual(self.contig.start, 0)
        self.assertEqual(len(self.contig.transcripts), 11)

    def test_non_overlapping(self):
        with self.assertRaises(IndexError):
            self.contig.add_transcript(self.faker.tm(0, 0, self.contig.end))

    def test_not_same_strand(self):
        with self.assertRaises(TypeError):
            self.contig.add_transcript(self.faker.tm(strand="-"))