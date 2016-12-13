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
        from genomicode import filelib
        #out_attributes = set_out_attributes(in_data, out_attributes)
        TCGA_BIN = config.download_tcga
        assert 'disease' in user_options
        if 'date' in user_options:
            x = ['--date', user_options['date']]
        else:
            x = []
    
        
        command = ['python', TCGA_BIN, '--disease', user_options['disease'],
                   '--data', out_attributes['preprocess'], '--download_and_extract'] + x
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        result_files = os.listdir(".")
        result_format = '.txt'
        for result_file in result_files:
            if result_file.endswith(result_format) and result_file != 'MANIFEST.txt':
                os.rename(result_file, outfile)

    
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for download_tcga_agilent fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'tcga_file' + original_file + '.txt'
        return filename
