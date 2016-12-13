from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from Betsy import read_label_file
        from genomicode import jmath
        
        cls_node_train, data_node = antecedents
        result, label_line, second_line = read_label_file.read(
            cls_node_train.identifier)
        y = [second_line[int(i)] for i in label_line]
        R = jmath.start_R()
        M = arrayio.read(data_node.identifier)
        M_train = M.matrix(None, range(0, len(label_line)))
        M_test = M.matrix(None, range(len(label_line), M.dim()[1]))
        M1 = M_train.slice()
        M_train = jmath.transpose(M1)
        jmath.R_equals_matrix(M_train, 'data')
        M2 = M_test.slice()
        M2 = jmath.transpose(M2)
        jmath.R_equals_matrix(M2, 'test')
        jmath.R_equals(y, 'y')
        R('y<-as.factor(y)')
        R('require(randomForest, quietly=TRUE)')
        R('library(randomForest)')
        R('model <- randomForest(data,y=y,importance=TRUE)')
        R('predict_result <- predict(model, test)')
        predict_result = R['predict_result']
        levels = predict_result.levels
        predict_labels = predict_result[:]
        predict_labels = [levels[i - 1] for i in predict_labels]
        name = M_test._col_names.keys()[0]
        sample_name = M_test._col_names[name]
        result = [['Sample_name', 'Predicted_class', 'Confidence']]
        for i in range(len(sample_name)):
            result.append([str(sample_name[i]), predict_labels[i], ''])
        
        f = file(outfile, 'w')
        for i in result:
            f.write('\t'.join(i))
            f.write('\n')
        f.close()


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node_train = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'random_forest_' + original_file + '.txt'
        return filename


    
