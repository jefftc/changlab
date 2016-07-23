from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import svmutil
        import arrayio
        from genomicode import filelib
        from Betsy import read_label_file
        from Betsy import module_utils
        
        data_node_train, cls_node_train = antecedents
        a, training_label, second_line = read_label_file.read(
            cls_node_train.identifier)
        training = arrayio.read(data_node_train.identifier)
        x_training = module_utils.format_convert(
            training.matrix(None, range(0, len(training_label)))
        )
          #convert to the format libsvm accept
        y_training = [int(x) for x in training_label]
        svm_kernel = ['linear', 'polynomial', 'RBF', 'sigmoid',
                      'precomputed_kernel']
        #if 'svm_kernel' in parameters.keys():
        kernel_type = svm_kernel.index(out_attributes['svm_kernel'])
        command = '-t ' + str(kernel_type)
        # else:
        #    command = '-t 0'
        param = svmutil.svm_parameter(command)
        prob = svmutil.svm_problem(y_training, x_training)
        model = svmutil.svm_train(prob, param)
        svmutil.svm_save_model(outfile, model)
        filelib.assert_exists_nz(outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node_train, cls_node_train = antecedents
        original_file = module_utils.get_inputid(data_node_train.identifier)
        filename = 'svm_model_' + original_file + '.txt'
        return filename


    
