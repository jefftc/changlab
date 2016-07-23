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
        from genomicode import filelib
        in_data = antecedents
        scoresig_path = config.scoresig
        scoresig_BIN = module_utils.which(scoresig_path)
        assert scoresig_BIN, 'cannot find the %s' % scoresig_path
        command = ['python', scoresig_BIN, '-r', in_data.identifier, '-m',
                   in_data.identifier, '-j', '20', '-o', outfile]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for run_scoresig does not exists' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signature_score' + original_file
        return filename



