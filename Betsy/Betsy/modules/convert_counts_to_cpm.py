from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import config

        signal_node = in_data
        signal_file = signal_node.identifier
        assert os.path.exists(signal_file)
        
        slice_matrix = filelib.which_assert(config.slice_matrix)

        sq = parallel.quote
        cmd = [
            sq(slice_matrix),
            "--cpm",
            signal_file,
            ]
        cmd = " ".join(cmd)
        cmd = "%s >& %s" % (cmd, outfile)

        parallel.sshell(cmd)
        filelib.assert_exists_nz(outfile)


    def name_outfile(self, antecedents, user_options):
        return "signal.cpm"


    
