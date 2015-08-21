from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from Betsy import read_label_file
        from Betsy import module_utils
        data_node, cls_node_test = antecedents
        f = file(data_node.identifier, 'r')
        text = f.readlines()
        f.close()
        a, test_label, second_line = read_label_file.read(cls_node_test.identifier)
        actual_label = [second_line[int(i)] for i in test_label]
        f = file(outfile, 'w')
        header = text[0].replace('\n', '').split('\t')
        header.extend(['Actual_class', 'Correct?'])
        f.write('\t'.join(header) + '\n')
        for index in range(1, len(text)):
            line = text[index].replace('\n', '')
            line = line.split('\t')
            correct = 'no'
            if line[1] == actual_label[index - 1]:
                correct = 'yes'
            line.extend([actual_label[index - 1], correct])
            f.write('\t'.join(line) + '\n')
        
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for evaluate_prediction' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node_test = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'prediction_' + original_file + '.txt'
        return filename


    def find_antecedents(
        self, network, module_id, out_attributes, user_attributes, pool):
        from Betsy import module_utils
        filter1 = module_utils.AntecedentFilter(datatype_name='ClassifyFile')
        filter2 = module_utils.AntecedentFilter(
            datatype_name='ClassLabelFile', contents="test")
        x = module_utils.find_antecedents(
            network, module_id, user_attributes, pool, filter1, filter2)
        return x
