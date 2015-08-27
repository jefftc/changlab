from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        data_node, rename_node = antecedents
        rename_path = config.slice_matrix
        rename_BIN = module_utils.which(rename_path)
        assert rename_BIN, 'cannot find the %s' % rename_path
        command = ['python', rename_BIN, data_node.identifier, '--relabel_col_ids',
                   rename_node.identifier + ',NewName']
        f = file(outfile, 'w')
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=f,
                                   stderr=subprocess.PIPE)
        f.close()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)

        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for relabel_samples does not exist' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, rename_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'signal_rename_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        data_node, rename_node = antecedents
        new_parameters = data_node.data.attributes.copy()
        new_parameters['rename_sample'] = 'yes'
        return new_parameters


    
