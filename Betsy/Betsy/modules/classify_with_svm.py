from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import svmutil
        import arrayio
        from Betsy import read_label_file
        from Betsy import module_utils
        svm_model, data_node_test, cls_node_train = antecedents
        a, train_label, second_line = read_label_file.read(
            cls_node_train.identifier)
        M = arrayio.read(data_node_test.identifier)
        # convert to the format libsvm accept
        test = M.matrix(None, range(len(train_label), M.dim()[1]))
        x_test = module_utils.format_convert(test)
        model = svmutil.svm_load_model(svm_model.identifier)
        a, train_label, second_line = read_label_file.read(
            cls_node_train.identifier)
        y_test = [0] * len(x_test)
        p_label, p_acc, p_val = svmutil.svm_predict(y_test, x_test, model)
        prediction_index = [int(i) for i in p_label]
        prediction = [second_line[i] for i in prediction_index]
        name = test._col_names.keys()[0]
        sample_name = test._col_names[name]
        result = [['Sample_name', 'Predicted_class', 'Confidence']]
        for i in range(len(sample_name)):
            result.append([str(sample_name[i]), prediction[i], str(p_val[i][0])])
        
        f = file(outfile, 'w')
        for i in result:
            f.write('\t'.join(i))
            f.write('\n')
        
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for classify_with_svm fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        svm_model, data_node_test, cls_node_train = antecedents
        original_file = module_utils.get_inputid(svm_model.identifier)
        filename = 'svm_result' + original_file + '.txt'
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        svm_model, data_node_test, cls_node_train = antecedents
        identifier = svm_model.identifier
        return module_utils.hash_input(identifier, pipeline, out_attributes,
                                             user_options)


    def find_antecedents(
        self, network, module_id, out_attributes, user_attributes, pool):
        from Betsy import module_utils
        filter1 = module_utils.AntecedentFilter(datatype_name='SvmModel')
        filter2 = module_utils.AntecedentFilter(
            datatype_name='SignalFile', contents="class0,class1,test")
        filter3 = module_utils.AntecedentFilter(
            datatype_name='ClassLabelFile', contents="class0,class1")
        x = module_utils.find_antecedents(
            network, module_id, user_attributes, pool, filter1, filter2, filter3)
        return x
