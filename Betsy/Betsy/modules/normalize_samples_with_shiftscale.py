from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import shiftscalenorm
        import arrayio
        from Betsy import read_label_file
        from genomicode import filelib
        data_node, cls_node = antecedents
        if data_node and cls_node:
            result, label_line, second_line = read_label_file.read(
                cls_node.identifier)
            assert len(result) == 2, 'for shiftscale,there should be only 2 classes'
            M = arrayio.read(data_node.identifier)
            index1 = result[0][0]
            index2 = result[1][0]
            M_1 = M.matrix(None, index1)
            M_2 = M.matrix(None, index2)
            M_y = shiftscalenorm.normalize(M_1, M_2)
            for i in range(M_y.dim()[0]):
                for j in range(M_y.dim()[1]):
                    if str(M_y._X[i][j]) == 'nan':
                        M_y._X[i][j] = M_2._X[i][0]
            for j in range(M.nrow()):
                for i in range(len(index1)):
                    M._X[j][index1[i]] = M_y._X[j][i]

            f = file(outfile, 'w')
            arrayio.tab_delimited_format.write(M, f)
            f.close()
            assert filelib.exists_nz(outfile), (
                'the output file %s for shiftscale fails' % outfile
            )
        
        
        return False


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'signal_shiftscale_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        data_node, cls_node = antecedents
        new_parameters = data_node.data.attributes.copy()
        new_parameters['shiftscale_norm'] = 'yes'
        return new_parameters


    
