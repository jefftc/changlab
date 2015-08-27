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
        cmd = [
            'python',
            diffexp_bin,
            data_node.identifier,
            '--cls_file', cls_node.identifier,
            '--algorithm', 'ttest',
            ]
        ## x = user_options.get("fold_change_cutoff")
        ## if x:
        ##     cmd = cmd + ['--fold_change', x]
        ## x = user_options.get("p_cutoff")
        ## if x:
        ##     cmd = cmd + ['--p_cutoff', x]
        ## x = user_options.get("bonf_cutoff")
        ## if x:
        ##     cmd = cmd + ['--bonf_cutoff', x]
        ## x = user_options.get("fdr_cutoff")
        ## if x:
        ##     cmd = cmd + ['--fdr_change', x]
        
        handle = open(outfile, 'w')
        try:
            process = subprocess.Popen(
                cmd, shell=False, stdout=handle, stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)
        finally:
            handle.close()


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 't_test_' + original_file + '.txt'
        return filename


    
