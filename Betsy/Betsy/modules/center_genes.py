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
        from Betsy import module_utils as mlib

        metadata = {}
        
        center_alg = {
            'mean': 'a',
            'median': 'm',
            }
        assert "gene_center" in out_attributes
        center = out_attributes['gene_center']
        assert center in center_alg, "Invalid center option: %s" % center
        center_parameter = center_alg[center]

        cluster = mlib.get_config("cluster", which_assert_file=True)
        sq = parallel.quote
        cmd = [
            sq(cluster),
            "-f", sq(in_data.identifier),
            "-cg", center_parameter,
            "-u", outfile,
            ]
        cmd = " ".join(map(str, cmd))
        parallel.sshell(cmd)
        metadata["commands"] = [cmd]
        
        outputfile = outfile + '.nrm'
        filelib.assert_exists_nz(outputfile)
        os.rename(outputfile, outfile)

        return metadata


    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'signal_center_' + original_file + '.tdf'
        #return filename
        return "signal.tdf"

    
    #def set_out_attributes(self, antecedents, out_attributes):
    #    new_parameters = out_attributes.copy()
    #    new_parameters['format'] = 'tdf'
    #    return new_parameters



