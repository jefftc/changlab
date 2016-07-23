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
        from genomicode import filelib
        in_data = antecedents
        sortsam_BIN = config.sortsam
        assert os.path.exists(sortsam_BIN), 'cannot find the %s' % sortsam_BIN
        command = ['java', '-Xmx5g', '-jar', sortsam_BIN,
                   'I=' + in_data.identifier, 'O=' + outfile, 'SO=coordinate',
                   'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true']
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()[1]
        if 'error' in error_message:
            raise ValueError(error_message)
    
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for sort_sam_file does not exist' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        return 'sorted_' + original_file + '.bam'

