from .abstract_parser import AbstractParser

class Gtf(AbstractParser):
    def set_data(self, line):
        self.data = line.split("\t")

    def get_source(self):
        return "gtf"

    def get_chromosome(self):
        return self.data[0]

    def get_type(self):
        return self.data[2]

    def get_start(self):
        return int(self.data[3])

    def get_end(self):
        return int(self.data[4])

    def get_strand(self):
        return self.data[6]

    def get_transcript_id(self):
        return self.data[8].split(";")[1].replace("transcript_id", "").replace("\"", "").replace(" ", "")

    def get_gene_id(self):
        return self.data[8].split(";")[0].replace("gene_id", "").replace("\"", "").replace(" ", "")

    def get_meta(self):
        return ""
