import unittest
from models.transcript import Transcript
from models.exon import Exon
from models.contig import Contig

class TestContigModel(unittest.TestCase):
    def setUp(self):
        self.sample_transcripts = [
            Transcript([Exon("gtf", "chr1", 0, 100, "+", "abc", "123")] )
        ]

        self.contig = Contig(self.sample_transcripts)

    def test_creates_data(self):
        self.assertEqual(self.contig.start, 0)
        self.assertEqual(self.contig.end, 100)
        self.assertEqual(len(self.contig.transcripts), 1)

    def test_adds_transcript(self):
        self.contig.add_transcript(Transcript([Exon("gtf", "chr1", 50, 300, "+", "def", "456")]))

        self.assertEqual(self.contig.start, 0)
        self.assertEqual(self.contig.end, 300)
        self.assertEqual(len(self.contig.transcripts), 2)

    def test_non_overlapping(self):
        with self.assertRaises(IndexError):
            self.contig.add_transcript(Transcript([Exon("gtf", "chr1", 400, 600, "+", "def", "456")]))

    def test_not_same_strand(self):
        with self.assertRaises(TypeError):
            self.contig.add_transcript(Transcript([Exon("gtf", "chr1", 150, 300, "-", "def", "456")]))

    def test_has_transcript(self):
        self.assertTrue(self.contig.has_transcript(self.sample_transcripts[0]))
        self.assertFalse(self.contig.has_transcript(Transcript([Exon("gtf", "chr1", 800, 900, "+", "abc", "123")])))