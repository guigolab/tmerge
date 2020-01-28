class AbstractParser():
    def get_source(self):
        raise NotImplementedError()

    def get_chromosome(self, line):
        raise NotImplementedError()

    def get_start(self, line):
        raise NotImplementedError()

    def get_end(self, line):
        raise NotImplementedError()

    def get_strand(self, line):
        raise NotImplementedError()

    def get_transcript_id(self, line):
        raise NotImplementedError()

    def get_gene_id(self, line):
        raise NotImplementedError()

    def get_meta(self, line):
        raise NotImplementedError()