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
        from Betsy import module_utils

        data_node, cel_node = antecedents
        #out_attributes = set_out_attributes(data_node, out_attributes)
        phenotype_BIN = config.analyze_phenotype
        assert os.path.exists(phenotype_BIN)
        assert "geneset_value" in user_options, 'no geneset are provided'
        if not os.path.exists(outfile):
            os.mkdir(outfile)

        command = [
            'python', phenotype_BIN, '--phenotype', 'EMT',
            '--ignore_samples', 'shCDH1,1', '--gene',
            user_options['geneset_value'], '-o', outfile + '/EMT',
            data_node.identifier, cel_node.identifier]
        process = subprocess.Popen(
            command, shell=False, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        x = process.communicate()
        error_message = x[1]
        assert not error_message, error_message
        assert module_utils.exists_nz(outfile), (
            'the output file %s for analyze_phenotype fails' % outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cel_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'EMT_' + original_file
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        data_node, cel_node = antecedents
        identifier = data_node.identifier
        return module_utils.hash_input(
            identifier, pipeline, out_attributes, user_options)
