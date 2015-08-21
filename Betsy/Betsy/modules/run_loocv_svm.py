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
        predict_model = __import__('Betsy.modules.' + 'classify_with_svm',
                                   globals(), locals(), ['classify_with_svm'], -2)
        train_model = __import__('Betsy.modules.' + 'train_svm_model', globals(),
                                 locals(), ['train_svm_model'], -2)
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
            del new_parameters['loocv']
            del new_parameters['actual_label']
            del new_parameters['wv_feature_stat']
            x1 = merge_data, train_label_data
            svm_model = train_model.run(x1, new_parameters, user_options, network)
            x = svm_model, merge_data, train_label_data
            out_node = predict_model.run(x, out_attributes, user_options, network)
            out_node_label = evaluate_model.run(
                (out_node, test_label_data), out_attributes, user_options, network)
            f1 = open(out_node_label.identifier, 'r')
            lines = f1.readlines()
            f1.close()
            f.write(lines[1])
            os.remove(merge_file)
            os.remove(train_label)
            os.remove(test_label)
            os.remove(out_node.identifier)
            os.remove(out_node_label.identifier)
            os.remove(svm_model.identifier)
        
        
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for loocv fails' % outfile
        )


    def find_antecedents(
        self, network, module_id, out_attributes, user_attributes, pool):
        from Betsy import module_utils
        filter1 = module_utils.AntecedentFilter(
            datatype_name='SignalFile', contents="class0,class1")
        filter2 = module_utils.AntecedentFilter(
            datatype_name='ClassLabelFile', contents="class0,class1")
        x = module_utils.find_antecedents(
            network, module_id, user_attributes, pool, filter1, filter2)
        return x


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'predication_loocv_svm' + original_file + '.txt'
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        identifier = data_node.identifier
        return module_utils.hash_input(identifier, pipeline, out_attributes,
                                             user_options)
