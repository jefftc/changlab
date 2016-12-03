from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
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
        ## fc = user_options.get("diffexp_foldchange_value")
        ## if fc:
        ##     cmd = cmd + ['--fold_change', fc]
        
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
        filename = 'ebayes_' + original_file + '.txt'
        return filename


    
