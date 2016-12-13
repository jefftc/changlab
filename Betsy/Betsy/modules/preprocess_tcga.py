from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import subprocess
        from genomicode import config
        from genomicode import filelib
        in_data = antecedents
        #out_attributes = set_out_attributes(in_data, out_attributes)
        TCGA_BIN = config.download_tcga
        command = ['python', TCGA_BIN, '--input', in_data.identifier, '--data',
                   in_data.data.attributes['preprocess'], '--process_only',
                   outfile]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        x = process.communicate()
        error_message = x[1]
        if error_message:
            raise ValueError(error_message)
    
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for preprocess_tcga fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signalfile_tcga_' + original_file + '.tdf'
        return filename



