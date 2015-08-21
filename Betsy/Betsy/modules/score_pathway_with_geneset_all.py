from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """analyze geneset"""
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        data_node, geneset_node = antecedents
        score_geneset_path = config.score_geneset
        score_geneset_BIN = module_utils.which(score_geneset_path)
        assert score_geneset_BIN, 'cannot find the %s' % score_geneset_path
        automatch = out_attributes['automatch']
        command = ['python', score_geneset_BIN, '-o', outfile, '--geneset_file',
                   geneset_node.identifier, data_node.identifier, '--all']
        if automatch == 'yes':
            command.append('--automatch')
        
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for score_pathway_with_geneset fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'score_geneset_' + original_file + '.txt'
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        identifier = data_node.identifier
        return module_utils.hash_input(identifier, pipeline, out_attributes,
                                             user_options)
