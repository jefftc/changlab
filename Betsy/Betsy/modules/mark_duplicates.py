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
        mark_duplicates_path = config.Mark_duplicates
        assert os.path.exists(
            mark_duplicates_path), 'cannot find the %s' % mark_duplicates_path
        command = ['java', '-Xmx5g', '-jar', mark_duplicates_path,
                   'I=' + in_data.identifier, 'O=' + outfile,
                   'METRICS_FILE=metricsFile', 'VALIDATION_STRINGENCY=LENIENT',
                   'REMOVE_DUPLICATES=true']
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        #error_message = process.communicate()[1]
        #if 'error' in error_message:
        #    raise ValueError(error_message)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for mark_duplcates does not exist' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'marked_duplicates_' + original_file + '.bam'
        return filename
