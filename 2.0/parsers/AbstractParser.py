class AbstractParser():
    def set_data(self, data):
        self.data = data

    def get_source(self):
        raise NotImplementedError()

    def get_chromosome(self):
        raise NotImplementedError()

    def get_start(self):
        raise NotImplementedError()

    def get_end(self):
        raise NotImplementedError()

    def get_strand(self):
        raise NotImplementedError()

    def get_transcript_id(self):
        raise NotImplementedError()

    def get_gene_id(self):
        raise NotImplementedError()

    def get_meta(self):
        raise NotImplementedError()