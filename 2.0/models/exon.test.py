import unittest
from exon import Exon

sample_data = {
    "source_type": "liff",
    "seq_name": "zaphod",
    "start": 0,
    "end": 42,
    "strand": "Betelgeuse",
    "transcript_id": "vogons",
    "gene_id": "constructor"
}

class TestExon(unittest.TestCase):

    def test_creates_data(self):
        exon = Exon(**sample_data)
        for k, v in sample_data.items():
            self.assertEqual(getattr(exon, k), v)

if __name__ == '__main__':
    unittest.main()