from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """compare signature predictions"""
        import os
        import arrayio
        from genomicode import filelib
        from Betsy import module_utils
        data_node1, data_node2, data_node3 = antecedents
        assert 'probabilities.pcl' in os.listdir(data_node1.identifier)
        assert 'probabilities.pcl' in os.listdir(data_node2.identifier)
        assert 'probabilities.pcl' in os.listdir(data_node3.identifier)
        f1 = os.path.join(data_node1.identifier, 'probabilities.pcl')
        f2 = os.path.join(data_node2.identifier, 'probabilities.pcl')
        f3 = os.path.join(data_node3.identifier, 'probabilities.pcl')
        M1 = arrayio.read(f1)
        M2 = arrayio.read(f2)
        M3 = arrayio.read(f3)
        total_sample = list(set(M1.col_names('_SAMPLE_NAME') + M2.col_names(
            '_SAMPLE_NAME') + M3.col_names('_SAMPLE_NAME')))
        d1 = get_new_matrix(total_sample, f1)
        d2 = get_new_matrix(total_sample, f2)
        d3 = get_new_matrix(total_sample, f3)
        f = file(outfile, 'w')
        header = ['SigID', 'NAME'] + total_sample
        f.write('\t'.join(header) + '\n')
        assert len(d1) == len(d2)
        assert len(d1) == len(d3)
        for i in range(len(d1)):
            d1[i][1] = 'Affy_' + d1[i][1]
            d2[i][1] = 'Agilent_' + d2[i][1]
            d3[i][1] = 'RSEM_' + d3[i][1]
            f.write('\t'.join(d1[i]) + '\n')
            f.write('\t'.join(d2[i]) + '\n')
            f.write('\t'.join(d3[i]) + '\n')
        
        f.close()
        assert filelib.exists_nz(outfile), (
            'the output file %s for compare_signature_predictions fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node1, data_node2, data_node3 = antecedents
        original_file = module_utils.get_inputid(data_node1.identifier)
        filename = 'compare_three_signature_predictions_' + original_file + '.txt'
        return filename


def get_new_matrix(total_sample, filename):
    import arrayio
    M = arrayio.read(filename)
    samples = M.col_names('_SAMPLE_NAME')
    data = []
    indexes = []
    for i in total_sample:
        index = samples.index(i) if i in samples else None
        indexes.append(index)
    
    for i in range(M.nrow()):
        line = [M.row_names('SigID')[i], M.row_names('NAME')[i]]
        for j in indexes:
            if j is None:
                line.append('NA')
            else:
                line.append(str(M._X[i][j]))
        data.append(line)
    
    return data
