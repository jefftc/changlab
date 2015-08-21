from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        import arrayio
        from Betsy import read_label_file
        from Betsy import module_utils
        from genomicode import config
        data_node_train, data_node_test, cls_node_train = antecedents
        module_name = 'WeightedVoting'
        gp_parameters = dict()
        file1, file2 = module_utils.convert_to_same_platform(
            data_node_train.identifier, data_node_test.identifier)
        result, label_line, class_name = read_label_file.read(
            cls_node_train.identifier)
        M = arrayio.read(data_node_test.identifier)
        label_line = ['0'] * M.dim()[1]
        read_label_file.write('temp_test.cls', class_name, label_line)
        gp_parameters['train.filename'] = file1
        gp_parameters['train.class.filename'] = cls_node_train.identifier
        gp_parameters['test.filename'] = file2
        gp_parameters['test.class.filename'] = 'temp_test.cls'
        if 'wv_num_features' in user_options:
            gp_parameters['num.features'] = str(user_options['wv_num_features'])
        
        if 'wv_minstd' in user_options:
            assert module_utils.is_number(
                user_options['wv_minstd']), 'the sv_minstd should be number'
            gp_parameters['min.std'] = str(user_options['wv_minstd'])

        
        wv_feature_stat = ['wv_snr', 'wv_ttest', 'wv_snr_median',
                           'wv_ttest_median', 'wv_snr_minstd', 'wv_ttest_minstd',
                           'wv_snr_median_minstd', 'wv_ttest_median_minstd']

        assert out_attributes['wv_feature_stat'] in wv_feature_stat, (
            'the wv_feature_stat is invalid'
        )
        gp_parameters['feature.selection.statistic'] = str(
            wv_feature_stat.index(out_attributes['wv_feature_stat']))
        gp_path = config.genepattern
        gp_module = module_utils.which(gp_path)
        assert gp_module, 'cannot find the %s' % gp_path
        download_directory = os.path.join(".", 'wv_result')
        command = [gp_module, module_name, '-o', download_directory]
        for key in gp_parameters.keys():
            a = ['--parameters', key + ':' + gp_parameters[key]]
            command.extend(a)
        
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        
        assert os.path.exists(download_directory), (
            'there is no output directory for weightedVoting'
        )
        result_files = os.listdir(download_directory)
        assert 'stderr.txt' not in result_files, 'gene_pattern get error'
        gp_files = os.listdir(download_directory)
        for gp_file in gp_files:
            if gp_file.endswith('pred.odf'):
                gp_file = os.path.join(download_directory, gp_file)
                f = file(gp_file, 'r')
                text = f.readlines()
                f.close()
                os.rename(os.path.join(download_directory, gp_file),
                          os.path.join(download_directory, 'prediction.odf'))
                assert text[1][0:12] == 'HeaderLines='
                start = int(text[1][12:-1])
                newresult = [['Sample_name', 'Predicted_class', 'Confidence']]
                for i in text[start + 2:]:
                    line = i.split()
                    n = len(line)
                    newline = [' '.join(line[0:n - 4]), line[n - 3], line[n - 2]]
                    newresult.append(newline)
                f = file(outfile, 'w')
                for i in newresult:
                    f.write('\t'.join(i))
                    f.write('\n')
                f.close()
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for classify_with_weighted_voting fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node_train, data_node_test, cls_node_train = antecedents
        original_file = module_utils.get_inputid(data_node_train.identifier)
        filename = 'weighted_voting_' + original_file + '.cdt'
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        data_node_train, data_node_test, cls_node_train = antecedents
        identifier = data_node_train.identifier
        return module_utils.hash_input(identifier, pipeline, out_attributes,
                                             user_options)


    def find_antecedents(
        self, network, module_id, out_attributes, user_attributes, pool):
        from Betsy import module_utils
        filter1 = module_utils.AntecedentFilter(
            datatype_name='SignalFile', contents="class0,class1")
        filter2 = module_utils.AntecedentFilter(
            datatype_name='SignalFile', contents="test")
        filter3 = module_utils.AntecedentFilter(
            datatype_name='ClassLabelFile', contents="class0,class1")
        x = module_utils.find_antecedents(
            network, module_id, user_attributes, pool, filter1, filter2, filter3)
        return x
