from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """clustering the input file"""
        import os
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        CLUSTER_BIN = config.cluster
        cluster = module_utils.which(CLUSTER_BIN)
        assert cluster, 'cannot find the %s' % CLUSTER_BIN
        distance_para = {'correlation': '1', 'euclidean': '7'}
        dist = distance_para[out_attributes['distance']]
        com_parameter = ["-g", dist, '-s', '-e', '1']
        command = [CLUSTER_BIN, '-f', in_data.identifier, '-u', outfile]
        for i in com_parameter:
            command.append(i)
    
        
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        result_files = os.listdir(".")
        result_format = 'txt'
        for result_file in result_files:
            if result_file.endswith(result_format):
                os.rename(result_file, outfile)
    
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for cluster_genes_by_som fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'cluster_file_' + original_file + '.cdt'
        return filename



