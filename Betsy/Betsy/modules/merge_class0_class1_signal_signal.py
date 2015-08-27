from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """merge three signal file to generate a joined signal file"""
        import os
        from Betsy import module_utils
        
        merge_node1, merge_node2 = antecedents
        assert os.path.exists(merge_node1.identifier), \
            'File not found: %s' % merge_node1.identifier
        assert os.path.exists(merge_node2.identifier), \
            'File not found: %s' % merge_node2.identifier
        
        file1, file2 = module_utils.convert_to_same_platform(
            merge_node1.identifier, merge_node2.identifier)
        f = file(outfile, 'w')
        module_utils.merge_two_files(file1, file2, f)
        f.close()


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node1, data_node2 = antecedents
        original_file = module_utils.get_inputid(data_node1.identifier)
        filename = 'signal_merge' + original_file + '.tdf'
        return filename


    
