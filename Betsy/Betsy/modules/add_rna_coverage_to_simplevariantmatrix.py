from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import filelib
        import add_coverage_to_simplevariantmatrix

        simple_node, coverage_node = antecedents
        filelib.assert_exists_nz(simple_node.identifier)
        filelib.assert_exists_nz(coverage_node.identifier)

        # Figure out if I'm adding coverage data from DNA or RNA.
        #in_attrs = simple_node.data.attributes
        #out_attrs = out_attributes
        #name = "with_rna_coverage"
        #assert name in in_attrs and name in out_attrs
        #is_rna_cov = False
        #if in_attrs[name] == "no" and out_attrs[name] == "yes":
        #    is_rna_cov = True
        add_coverage_to_simplevariantmatrix.add_coverage_to_svm(
            simple_node.identifier, coverage_node.identifier, out_filename,
            True)


    def name_outfile(self, antecedents, user_options):
        return "matrix.txt"
