from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        rma_node, mas5_node = antecedents
        scoresig_path = config.scoresig
        scoresig_BIN = module_utils.which(scoresig_path)
        assert scoresig_BIN, 'cannot find the %s' % scoresig_path
        file1, file2 = module_utils.convert_to_same_platform(rma_node.identifier,
                                                             mas5_node.identifier)
        command = ['python', scoresig_BIN, '-r', file1, '-m', file2, '-j', '20',
                   '-o', outfile]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for score_path_with_scoresig_affymetrix does not exists'
            % outfile)


    def find_antecedents(
        self, network, module_id, out_attributes, user_attributes, pool):
        from Betsy import module_utils
        filter1 = module_utils.AntecedentFilter(
            datatype_name='SignalFile', preprocess="rma")
        filter2 = module_utils.AntecedentFilter(
            datatype_name='SignalFile', preprocess="mas5")
        x = module_utils.find_antecedents(
            network, module_id, user_attributes, pool, filter1, filter2)
        return x


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        rma_node, mas5_node = antecedents
        original_file = module_utils.get_inputid(rma_node.identifier)
        filename = 'signature_score' + original_file
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        rma_node, mas5_node = antecedents
        identifier = rma_node.identifier
        return module_utils.hash_input(identifier, pipeline, out_attributes,
                                             user_options)
