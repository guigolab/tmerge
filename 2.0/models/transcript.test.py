import unittest
from exon import Exon
from transcript import Transcript

class TestTranscriptModel(unittest.TestCase):
    def setUp(self):
        self.sample_exons = [
            Exon("gtf", "chr1", 0, 42, "+", "abc", "def")
        ]
        self.transcript = Transcript(self.sample_exons)

    def test_creates_data(self):
        self.assertEqual(self.transcript.exons, self.sample_exons)
        self.assertEqual(self.transcript.start, 0)
        self.assertEqual(self.transcript.end, 42)
        self.assertEqual(self.transcript.strand, "+")
        self.assertEqual(self.transcript.id, "abc")

    def test_adds_exon(self):
        exon = Exon("gtf", "chr1", 43, 48, "+", "abc", "gfa")
        self.transcript.add_exon(exon)

        self.assertEqual(len(self.transcript.exons), 2)

    def test_not_same_strand(self):
        exon = Exon("gtf", "chr1", 43, 48, "-", "abc", "gfa")
        with self.assertRaises(TypeError):
            self.transcript.add_exon(exon)

    def test_not_same_transcript(self):
        exon = Exon("gtf", "chr1", 43, 48, "+", "def", "gfa")
        with self.assertRaises(TypeError):
            self.transcript.add_exon(exon)

    def test_multiple_exons(self):
        self.sample_exons.extend([
            Exon("gtf", "chr1", 43, 49, "+", "abc", "def"),
            Exon("gtf", "chr1", 58, 981, "+", "abc", "def")
        ])

        self.transcript = Transcript(self.sample_exons)
        self.assertEqual(self.transcript.end, 981)

    def test_multiple_exons_failure(self):
        self.sample_exons.extend([
            Exon("gtf", "chr1", 43, 49, "+", "abc", "def"),
            Exon("gtf", "chr1", 58, 981, "-", "abc", "def")
        ])

        with self.assertRaises(TypeError):
            Transcript(self.sample_exons)

        self.sample_exons = [
            Exon("gtf", "chr1", 43, 49, "+", "abc", "def"),
            Exon("gtf", "chr1", 84, 983, "+", "defg", "def")
        ]

        with self.assertRaises(TypeError):
            Transcript(self.sample_exons)

if __name__ == '__main__':
    unittest.main()