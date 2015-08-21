from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """mean or median"""
        import os
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        CLUSTER_BIN = config.cluster
        cluster = module_utils.which(CLUSTER_BIN)
        assert cluster, 'cannot find the %s' % CLUSTER_BIN
        center_alg = {'mean': 'a', 'median': 'm'}
        try:
            center_parameter = center_alg[out_attributes['gene_center']]
        except:
            raise ValueError("Centering parameter is not recognized")
    
        
        process = subprocess.Popen([CLUSTER_BIN, '-f', in_data.identifier, '-cg',
                                    center_parameter, '-u', outfile],
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        outputfile = outfile + '.nrm'
        os.rename(outputfile, outfile)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for centering fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_center_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        new_parameters = out_attributes.copy()
        new_parameters['format'] = 'tdf'
        return new_parameters



