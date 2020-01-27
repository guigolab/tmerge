from models.exon import Exon
from models.transcript import Transcript
from parsers import AbstractParser

class Gtf(AbstractParser.AbstractParser):
    def get_source(self):
        return "gtf"

    def get_chromosome(self, data):
            return data[0]

    def get_start(self, data):
        return int(data[3])

    def get_end(self, data):
        return int(data[4])

    def get_strand(self, data):
        return data[6]

    def get_transcript_id(self, data):
        return data[8].split(";")[1].replace("transcript_id", "").replace("\"", "").replace(" ", "")

    def get_gene_id(self, data):
        return data[8].split(";")[0].replace("gene_id", "").replace("\"", "").replace(" ", "")

    def get_meta(self, data):
        return ""
