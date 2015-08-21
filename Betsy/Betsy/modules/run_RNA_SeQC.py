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
        in_data = antecedents
        directory = module_utils.unzip_if_zip(in_data.identifier)
        filenames = os.listdir(directory)
        assert filenames, 'The input folder or zip file is empty.'
        if not os.path.exists(outfile):
            os.mkdir(outfile)
    
        
        RNA_SeQC_path = config.RNASeQC
        REF = user_options['RNA_ref']
        GTF = user_options['RNA_gtf']
        assert os.path.exists(RNA_SeQC_path), ('cannot find the %s' % RNA_SeQC_path
                                             )
        assert os.path.exists(REF), ('cannot find the %s' % REF)
        assert os.path.exists(GTF), ('cannot find the %s' % GTF)
        for filename in filenames:
            infile = os.path.join(directory, filename)
            outname = os.path.splitext(filename)[0]

            command = ['java', '-jar', RNA_SeQC_path, '-o', outfile, '-r', REF,
                       '-s', outname + '|' + infile + '|NA', '-t', GTF]
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()
            # guess if there is error message in the output
            if 'error' in error_message[1]:
                raise ValueError(error_message)
    
    
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'bamFolder_' + original_file
        return filename
