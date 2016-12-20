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
            "ACTB", 60,       # Human beta actin.
            "TUBB", 203068,   # Human beta tubulin.
            "Actb", 22461,    # Mouse beta actin.
            "Tubb4a", 22153,  # Mouse beta tubulin.
            ]

        infile = in_data.identifier
        
        sq = parallel.quote
        cmd = [
            sq(lineplot),
            "--gene_names", ",".join(map(str, gene_names)),
            "--mar_bottom", 1.50,
            sq(infile),
            sq(outfile),
            ]
        cmd = " ".join(map(str, cmd))
        parallel.sshell(cmd)
        metadata["commands"] = [cmd]

        filelib.assert_exists_nz(outfile)
        return metadata
            

    def name_outfile(self, antecedents, user_options):
        return "actb_plot.png"


