from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        data_node, cls_node = antecedents
        diffexp_bin = config.find_diffexp_genes
        assert os.path.exists(diffexp_bin)
        cmd = ['python', diffexp_bin, data_node.identifier, '--cls_file',
               cls_node.identifier, '--algorithm', 'ebayes']
        if 'diffexp_foldchange_value' in user_options:
            foldchange = float(user_options['diffexp_foldchange_value'])
            cmd = cmd + ['--fold_change', str(foldchange)]
        
        handle = open(outfile, 'w')
        try:
            process = subprocess.Popen(cmd,
                                       shell=False,
                                       stdout=handle,
                                       stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)
        finally:
            handle.close()
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for calc_diffexp_with_ebayes fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'ebayes_' + original_file + '.txt'
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        identifier = data_node.identifier
        return module_utils.hash_input(
            identifier, pipeline, out_attributes, user_options)
