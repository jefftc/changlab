from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from Betsy import rule_engine_bie3
        import os
        import arrayio
        from Betsy import rulebase
        from Betsy import read_label_file
        from Betsy import module_utils
        data_node, cls_node = antecedents
        M = arrayio.read(data_node.identifier)
        a, training_label, second_line = read_label_file.read(cls_node.identifier)
        full_index = range(M.ncol())
        predict_model = __import__(
            'Betsy.modules.' + 'classify_with_random_forest', globals(), locals(),
            ['classify_with_random_forest'], -2)
        evaluate_model = __import__('Betsy.modules.' + 'evaluate_prediction',
                                    globals(), locals(), ['evaluate_prediction'],
                                    -2)
        f = file(outfile, 'w')
        f.write('\t'.join(['sample_name', 'Predicted_class', 'Confidence',
                           'Actual_class', 'Correct?']))
        f.write('\n')
        for i in range(M.ncol()):
            test_index = i
            train_index = full_index[:]
            train_index.remove(test_index)
            merge_index = train_index + [test_index]
            y_training = [training_label[x] for x in train_index]
            y_test = [training_label[test_index]]
            M_merge = M.matrix(None, merge_index)
            merge_file = 'merge' + '_' + str(i)
            f_out = file(merge_file, 'w')
            arrayio.gct_format.write(M_merge, f_out)
            f_out.close()
            train_label = 'train_label' + '_' + str(i)
            test_label = 'test_label' + '_' + str(i)
            read_label_file.write(train_label, second_line, y_training)
            read_label_file.write(test_label, second_line, y_test[0])
            merge_node = rulebase.SignalFile.output(format='gct',
                                                    contents='class0,class1,test')
            merge_data = rule_engine_bie3.IdentifiedDataNode(merge_node,
                                                     identifier=merge_file)
            train_label_node = rulebase.ClassLabelFile.output(
                contents='class0,class1')
            train_label_data = rule_engine_bie3.IdentifiedDataNode(train_label_node,
                                                           identifier=train_label)
            test_label_node = rulebase.ClassLabelFile.output(contents='test')
            test_label_data = rule_engine_bie3.IdentifiedDataNode(test_label_node,
                                                          identifier=test_label)
            new_parameters = out_attributes.copy()
            x = merge_data, train_label_data
            out_node = predict_model.run(x, out_attributes, user_options, network,
                                         num_cores)
            out_node_label = evaluate_model.run(
                (out_node,
                 test_label_data), new_parameters, user_options, network, num_cores)
            f1 = open(out_node_label.identifier, 'r')
            lines = f1.readlines()
            f1.close()
            f.write(lines[1])
            os.remove(merge_file)
            os.remove(train_label)
            os.remove(test_label)
            os.remove(out_node.identifier)
            os.remove(out_node_label.identifier)
        
        
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for loocv fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'predication_loocv_random_forest' + original_file + '.txt'
        return filename


    
