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
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        bfrm_path = config.bfrmnorm
        bfrm_BIN = module_utils.which(bfrm_path)
        assert bfrm_BIN, 'cannot find the %s' % bfrm_path
        num_factor = 1
        #num_factor = 10
        if 'num_factors' in user_options.keys():
            num_factor = int(user_options['num_factors'])
            assert num_factor >= 1, 'the num_factor should be >=1'
            # What is single_object?
            #M = arrayio.read(single_object.identifier)
            M = arrayio.read(in_data.identifier)
            col_num = M.ncol()
            assert num_factor <= col_num, (
                'the num_factor should be less than %d' % col_num
            )
    
        
        tmp = 'tmp_dir'
        command = ['python', bfrm_BIN, in_data.identifier, '-f', str(num_factor),
                   '-o', tmp]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        assert module_utils.exists_nz(tmp), (
            'the output dir %s for bfrm_normalize fails' % tmp
        )
        assert module_utils.exists_nz(os.path.join(tmp, 'normalized.gct')), (
            'the output gct file for bfrm_normalize fails'
        )
        out = os.path.join(tmp, 'normalized.gct')
        M = arrayio.read(out)
        M_new = arrayio.convert(M, to_format=arrayio.pcl_format)
        f = file(outfile, 'w')
        arrayio.tab_delimited_format.write(M_new, f)
        f.close()


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_bfrm_' + original_file + '.tdf'
        return filename


    def set_out_attributes(self, antecedents, out_attributes):
        new_parameters = antecedents.data.attributes.copy()
        new_parameters['bfrm_norm'] = 'yes'
        return new_parameters



