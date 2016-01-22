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
    
        
        mark_duplicates_path = config.Mark_duplicates
        assert os.path.exists(
            mark_duplicates_path), 'cannot find the %s' % mark_duplicates_path
        for filename in filenames:
            infile = os.path.join(directory, filename)
            outname = os.path.splitext(filename)[-2] + '.bam'
            outname = os.path.join(outfile, outname)
            command = ['java', '-Xmx5g', '-jar', mark_duplicates_path,
                       'I=' + infile, 'O=' + outname, 'METRICS_FILE=metricsFile',
                       'VALIDATION_STRINGENCY=LENIENT', 'REMOVE_DUPLICATES=true']
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()
            if 'error' in error_message[1]:
                raise ValueError(error_message)
            assert module_utils.exists_nz(outname), (
                'the output file %s for flag_dups_in_bam_folder does not exist' %
                outname)
    



    
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'bamFiles_dups' + original_file
        return filename

