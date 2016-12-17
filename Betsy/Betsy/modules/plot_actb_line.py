from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils as mlib
        
        in_data = antecedents
        metadata = {}

        lineplot = mlib.get_config("lineplot", which_assert_file=True)

        gene_names = [
            "ACTB",
            "TUBB",
            ]

        infile = in_data.identifier
        
        sq = parallel.quote
        cmd = [
            sq(lineplot),
            "--gene_names", ",".join(gene_names),
            sq(infile),
            sq(outfile),
            ]
        cmd = " ".join(cmd)
        parallel.sshell(cmd)
        metadata["commands"] = [cmd]

        filelib.assert_exists_nz(outfile)
        return metadata
            

    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'actb_plot' + original_file + '.png'
        #return filename
        return "actb_plot.png"


