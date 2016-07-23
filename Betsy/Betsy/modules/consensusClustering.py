from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        from genomicode import filelib
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        module_name = 'ConsensusClustering'
        gp_parameters = dict()
        file1 = in_data.identifier

        gp_parameters['input.filename'] = file1

        if 'cc_kmax' in user_options:
            assert module_utils.is_number(
                user_options['cc_kmax']), 'the cc_kmax should be number'
            gp_parameters['cc_kmax'] = str(user_options['cc_kmax'])
    
        
        if 'cc_resampling_iter' in user_options:
            assert module_utils.is_number(user_options['cc_resampling_iter'
                          ]), 'the cc_resampling_iter should be number'
            gp_parameters['resampling.iterations'] = str(
                user_options['cc_resampling_iter'])
    
        
        if 'cc_seed_value' in user_options:
            assert module_utils.is_number(
                user_options['cc_seed_value']), 'the cc_seed_value should be number'
            gp_parameters['seed.value'] = str(user_options['cc_seed_value'])
    
        
        if 'cc_decent_iter' in user_options:
            assert module_utils.is_number(
                user_options['cc_decent_iter']), 'the cc_decent_iter should be number'
            gp_parameters['descent.iterations'] = str(user_options['cc_decent_iter'])
    
        
        if 'cc_norm_iter' in user_options:
            assert module_utils.is_number(
                user_options['cc_norm_iter']), 'the cc_norm_iter should be number'
            gp_parameters['normalization.iterations'] = str(
                user_options['cc_norm_iter'])
    
        
        if 'cc_heatmap_size' in user_options:
            assert module_utils.is_number(
                user_options['cc_heatmap_size']), 'the cc_heatmap_size should be number'
            gp_parameters['heat.map.size'] = str(user_options['cc_heatmap_size'])

    
        
        gp_parameters['clustering.algorithm'
                     ] = out_attributes['Consensus_algorithm'].upper()
        gp_parameters['distance.measure'] = out_attributes['cc_distance'].upper()
        gp_parameters['merge.type'] = out_attributes['merge_type']
        gp_parameters['create.heat.map'] = out_attributes['create_heatmap']
        gp_parameters['normalize.type'] = out_attributes['normalize_type']
        gp_parameters['cluster.by'] = out_attributes['clusterby']
        gp_parameters['resample'] = out_attributes['cc_resample']

        gp_path = config.genepattern
        gp_module = module_utils.which(gp_path)
        assert gp_module, 'cannot find the %s' % gp_path
        download_directory = outfile
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
            'there is no output directory for consensusClustering'
        )
        result_files = os.listdir(download_directory)
        assert 'stderr.txt' not in result_files, 'gene_pattern get error'
        assert filelib.exists_nz(outfile), (
            'the output file %s for consensusClustering fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'ConsensusClusteringFolder' + original_file
        return filename



