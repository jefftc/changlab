from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        #from Betsy import read_label_file
        from genomicode import filelib
        from genomicode import config
        #from genomicode import arraysetlib
        from Betsy import module_utils
        
        cls_node_train, data_node_train = antecedents
        #result, label_line, class_name = read_label_file.read(
        #    cls_node_train.identifier)
        #x = arraysetlib.read_cls_file(cls_node_train.identifier)
        #class_names, classes = x

        metadata = {}
        
        module_name = 'WeightedVotingXValidation'
        module_id_version = '00028:2'
        
        gp_params = dict()
        gp_params['data.filename'] = data_node_train.identifier
        gp_params['class.filename'] = cls_node_train.identifier
        if 'wv_num_features' in user_options:
            gp_params['num.features'] = str(user_options['wv_num_features'])
    ##    if 'wv_minstd' in user_input:	
    ##    	assert module_utils.is_number(
    ##            user_input['wv_minstd']), 'the sv_minstd should be number'
    ##        gp_parameters['min.std'] = str(user_input['wv_minstd'])
    ##        
    ##    wv_feature_stat = ['wv_snr', 'wv_ttest', 'wv_snr_median',
    ##                       'wv_ttest_median', 'wv_snr_minstd',
    ##                       'wv_ttest_minstd', 'wv_snr_median_minstd',
    ##                       'wv_ttest_median_minstd']
    ##    
    ##    assert parameters['wv_feature_stat'] in wv_feature_stat, (
    ##            'the wv_feature_stat is invalid')
    ##    gp_parameters['feature.selection.statistic'] = str(
    ##            wv_feature_stat.index(parameters[
    ##                'wv_feature_stat']))

        
        gp_path = config.genepattern
        gp_module = filelib.which(gp_path)
        assert gp_module, 'cannot find the %s' % gp_path
        
        download_directory = os.path.join(".", 'wv_result')
        command = [
            gp_module,
            module_name,
            '--id_and_version', module_id_version,
            '-o', download_directory,
            ]
        for key in gp_params:
            x = ['--parameters', "%s:%s" % (key, gp_params[key])]
            command.extend(x)
        x = " ".join(map(str, x))
        metadata["commands"] = [x]

        # DEBUG: If this is already run, don't run it again.
        if not os.path.exists(download_directory):
            process = subprocess.Popen(
                command,
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()[1]

            # Ignore warning:
            # /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/
            #   python2.7/site-packages/rpy2/rinterface/__init__.py:185:
            #   RRuntimeWarning: Loading required package: rJava
            # 
            # warnings.warn(x, RRuntimeWarning)
            x = error_message.strip()
            if x.endswith("warnings.warn(x, RRuntimeWarning)"):
                pass
            else:
                assert not x, error_message
        filelib.assert_exists(download_directory)

        # Find the prediction file.
        x = os.listdir(download_directory)
        assert 'stderr.txt' not in x, 'gene_pattern get error'
        x = [x for x in x if x.endswith("pred.odf")]
        assert x, "Missing predictions file"
        assert len(x) == 1, "Too many prediction files: %s" % repr(x)
        gp_file = x[0]

        gp_file = os.path.join(download_directory, gp_file)
        text = open(gp_file).readlines()
        #os.rename(os.path.join(download_directory, gp_file),
        #          os.path.join(download_directory, 'prediction.odf'))
        assert text[1][0:12] == 'HeaderLines='
        start = int(text[1][12:-1])
        newresult = [['Sample_name', 'Predicted_class', 'Confidence',
                      'Actual_class', 'Correct?']]
        for i in text[start + 2:]:
            line = i.split()
            n = len(line)
            newline = [' '.join(line[0:n - 4]), line[n - 3], line[n - 2],
                       line[n - 4], line[n - 1]]
            newresult.append(newline)

        f = file(outfile, 'w')
        for i in newresult:
            f.write('\t'.join(i))
            f.write('\n')
        f.close()

        return metadata


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'predication_loocv_wv' + original_file + '.txt'
        return filename
