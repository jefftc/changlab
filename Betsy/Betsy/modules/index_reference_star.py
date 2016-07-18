from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib

        ref_node, gene_node = antecedents
        ref = alignlib.standardize_reference_genome(
            ref_node.identifier, out_path, use_symlinks=True)
        filelib.safe_mkdir(out_path)

        x = alignlib.make_STAR_index_command(
            ref.fasta_file_full, out_path, gtf_file=gene_node.identifier,
            num_cores=num_cores)
        x = "%s >& out.txt" % x
        parallel.sshell(x, path=out_path)

        # Check to make sure index was created successfully.
        alignlib.assert_is_STAR_reference(out_path)


    def name_outfile(self, antecedents, user_options):
        return "reference.star"
