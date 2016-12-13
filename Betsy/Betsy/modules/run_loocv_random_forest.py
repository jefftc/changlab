from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import arrayio
        from genomicode import filelib
        from Betsy import bie3
        from Betsy import rulebase
        from Betsy import read_label_file
        
        cls_node, data_node = antecedents
        M = arrayio.read(data_node.identifier)
        x = read_label_file.read(cls_node.identifier)
        a, training_label, second_line = x
        
        predict_model = __import__(
            'Betsy.modules.' + 'classify_with_random_forest',
            globals(), locals(),
            ['classify_with_random_forest'], -2)
        evaluate_model = __import__(
            'Betsy.modules.' + 'evaluate_prediction',
            globals(), locals(), ['evaluate_prediction'], -2)
        
        full_index = range(M.ncol())
        
        f = file(outfile, 'w')
        f.write('\t'.join(['sample_name', 'Predicted_class', 'Confidence',
                           'Actual_class', 'Correct?']))
        f.write('\n')
        for i in range(M.ncol()):
            # Make filenames
            # gene expression for N samples.
            merge_file = 'merge' + '_' + str(i)
            # class label file for the training samples (samples 1-(N-1)).
            train_label = 'train_label' + '_' + str(i)
            # class label file for the test sample (sample N).
            test_label = 'test_label' + '_' + str(i)
            # Save the output of the prediction and evaluation.
            predict_file = "predict.txt"
            evaluate_file = "evaluate.txt"

            test_index = i
            train_index = full_index[:]
            train_index.remove(test_index)
            merge_index = train_index + [test_index]
            y_training = [training_label[x] for x in train_index]
            y_test = [training_label[test_index]]

            # Write the files for this iteration.
            M_merge = M.matrix(None, merge_index)
            arrayio.gct_format.write(M_merge, open(merge_file, 'w'))
            read_label_file.write(train_label, second_line, y_training)
            read_label_file.write(test_label, second_line, y_test[0])

            # Make objects to be used in this analysis.
            x = rulebase.SignalFile.output(
                format='gct', contents='class0,class1,test')
            merge_data = bie3.IdentifiedDataNode(x, identifier=merge_file)
            x = rulebase.ClassLabelFile.output(contents='class0,class1')
            train_label_data = bie3.IdentifiedDataNode(
                x, identifier=train_label)
            x = rulebase.ClassLabelFile.output(contents='test')
            test_label_data = bie3.IdentifiedDataNode(x, identifier=test_label)

            # Make a fake object to pass to evaluate_model.run.
            out_node = filelib.GenericObject()
            out_node.identifier = predict_file

            # Run the predictions.
            x = train_label_data, merge_data
            predict_model.Module().run(
                network, x, out_attributes, user_options, num_cores,
                predict_file)

            # Run the evaluation.
            new_parameters = out_attributes.copy()
            x = test_label_data, out_node
            evaluate_model.Module().run(
                network, x, new_parameters, user_options, num_cores,
                evaluate_file)

            # Is this the right line?
            lines = open(evaluate_file).readlines()
            
            f.write(lines[1])
            os.remove(merge_file)
            os.remove(train_label)
            os.remove(test_label)
            os.remove(predict_file)
            os.remove(evaluate_file)
        
        f.close()
        

    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'predication_loocv_random_forest' + original_file + '.txt'
        return filename


    
