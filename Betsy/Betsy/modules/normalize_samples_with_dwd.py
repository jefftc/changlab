from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import dwdnorm
        import arrayio
        from Betsy import read_label_file
        from Betsy import module_utils
        data_node, cls_node = antecedents
        # BUG: What happens if no antecedents?
        if data_node and cls_node:
            M = arrayio.read(data_node.identifier)
            result, label_line, second_line = read_label_file.read(
                cls_node.identifier)
            assert len(result) == 2, 'for dwd,there should be only 2 classes'
            assert [i in ['0', '1']
                    for i in label_line] == [True] * len(label_line), (
                        'the label of class shoul be 0 and 1'
                    )
            y = [i.replace('0', '-1') for i in label_line]
            M_y = dwdnorm.normalize(M, y)
            f = file(outfile, 'w')
            arrayio.tab_delimited_format.write(M_y, f)
            f.close()
            assert module_utils.exists_nz(outfile), (
                'the output file %s for dwd fails' % outfile
            )
            #new_parameters = out_attributes.copy()
        #return False


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'signal_dwd_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        data_node, cls_node = antecedents
        new_parameters = data_node.data.attributes.copy()
        new_parameters['dwd_norm'] = 'yes'
        return new_parameters

    
