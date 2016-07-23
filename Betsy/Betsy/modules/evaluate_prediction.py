from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from Betsy import read_label_file
        from Betsy import module_utils
        
        cls_node_test, data_node = antecedents
        text = open(data_node.identifier).readlines()
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
        #assert filelib.exists_nz(outfile), (
        #    'the output file %s for evaluate_prediction' % outfile
        #)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node_test = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'prediction_' + original_file + '.txt'
        return filename
