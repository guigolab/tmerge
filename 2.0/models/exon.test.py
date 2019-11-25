import unittest
from exon import Exon
from dataclass_type_validator import TypeValidationError

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


    def test_start_stop_ints(self):
        sample_data["start"] = "0"
        sample_data["end"] = "42"

        with self.assertRaises(TypeValidationError):
            exon = Exon(**sample_data)

    def test_start_lt_stop(self):
        sample_data["start"] = 42
        sample_data["end"] = 0

        with self.assertRaises(IndexError):
            Exon(**sample_data)

if __name__ == '__main__':
    unittest.main()