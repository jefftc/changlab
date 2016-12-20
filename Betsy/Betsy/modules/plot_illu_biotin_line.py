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

        #module_utils.plot_line_keywd(in_data.identifier, 'biotin', outfile)
        lineplot = mlib.get_config("lineplot", which_assert_file=True)

        sq = parallel.quote
        cmd = [
            sq(lineplot),
            "--gene_names", "biotin",
            "--mar_bottom", 1.50,
            "--yaxis_starts_at_0",
            sq(in_data.identifier),
            sq(outfile),
            ]
        cmd = " ".join(map(str, cmd))
        parallel.sshell(cmd)
        metadata["commands"] = [cmd]
        filelib.assert_exists_nz(outfile)
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "biotin.png"


