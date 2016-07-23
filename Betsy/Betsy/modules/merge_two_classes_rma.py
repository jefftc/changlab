from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """merge three signal file to generate a joined signal file"""
        import os
        from genomicode import filelib
        from Betsy import module_utils
        merge_node1, merge_node2 = antecedents
        assert os.path.exists(merge_node1.identifier), (
            'the merge_file1 %s in merge_data does not exist' % merge_node1.identifier
        )
        assert os.path.exists(merge_node2.identifier), (
            'the merge_file2 %s in merge_data does not exist' % merge_node2.identifier
        )
        file1, file2 = module_utils.convert_to_same_platform(
            merge_node1.identifier, merge_node2.identifier)
        f = file(outfile, 'w')
        module_utils.merge_two_files(file1, file2, f)
        f.close()
        assert filelib.exists_nz(outfile), (
            'the output file %s for merge_data fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node1, data_node2 = antecedents
        original_file = module_utils.get_inputid(data_node1.identifier)
        filename = 'signal_merge' + original_file + '.tdf'
        return filename


    
