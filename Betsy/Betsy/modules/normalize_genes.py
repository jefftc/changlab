from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import os
        import arrayio
        from genomicode import jmath
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils as mlib

        metadata = {}

        norm_para = ["variance", "sum_of_squares"]
        assert "gene_normalize" in out_attributes
        normalize = out_attributes["gene_normalize"]
        assert normalize in norm_para, \
               "Invalid normalize option: %s" % normalize

        if normalize == "variance":
            f = file(outfile, 'w')
            M = arrayio.read(in_data.identifier, format=arrayio.pcl_format)
            M_n = jmath.safe_norm_mv(M.slice())
            M._X = M_n
            M_c = arrayio.convert(M, to_format=arrayio.pcl_format)
            arrayio.pcl_format.write(M_c, f)
            f.close()
        elif normalize == "sum_of_squares":
            cluster = mlib.get_config("cluster", which_assert_file=True)
            sq = parallel.quote
            cmd = [
                sq(cluster),
                "-f", sq(in_data.identifier),
                "-ng", 
                "-u", outfile,
                ]
            parallel.sshell(cmd)
            metadata["command"] = cmd
            outputfile = outfile + '.nrm'
            filelib.assert_exists_nz(outputfile)
            os.rename(outputfile, outfile)
    
        filelib.assert_exists_nz(outfile)
        return metadata


    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'signal_normalize_' + original_file + '.tdf'
        #return filename
        return "signal.tdf"
    

    #def set_out_attributes(self, antecedents, out_attributes):
    #    new_parameters = out_attributes.copy()
    #    new_parameters['format'] = 'tdf'
    #    return new_parameters



