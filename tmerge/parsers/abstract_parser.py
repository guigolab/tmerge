class AbstractParser():
    def set_data(self, line):
        # Should seperate the line here so that it's accesible via other methods.
        # I.e. seperate by comma or tab for TSVs and CSVs
        raise NotImplementedError()

    def get_source(self):
        raise NotImplementedError()

    def get_type(self):
        # I.e. column 3 of GFF
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