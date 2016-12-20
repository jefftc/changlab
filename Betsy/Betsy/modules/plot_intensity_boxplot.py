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


        #M = arrayio.read(in_data.identifier)
        #data = jmath.transpose(M._X)
        #tickname = M._col_names['_SAMPLE_NAME']
        #fig = mplgraph.boxplot(
        #    data,
        #    xlabel='Sample Name',
        #    ylabel='Signal',
        #    title='Signal Intensity',
        #    box_label=tickname)
        #fig.savefig(outfile)

        boxplot = mlib.get_config("boxplot", which_assert_file=True)

        sq = parallel.quote
        cmd = [
            sq(boxplot),
            sq(in_data.identifier),
            sq(outfile),
            ]
        cmd = " ".join(map(str, cmd))
        parallel.sshell(cmd)

        metadata["commands"] = [cmd]
        
        filelib.assert_exists_nz(outfile)
        
        return metadata
        

    def name_outfile(self, antecedents, user_options):
        return "signal_distribution.png"


