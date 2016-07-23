from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        import shutil
        from genomicode import config
        from genomicode import filelib
        from Betsy import module_utils
        
        data_node, clinical_node = antecedents
        analyze_path = config.analyze_clinical
        analyze_BIN = module_utils.which(analyze_path)
        assert analyze_BIN, 'cannot find the %s' % analyze_path
        x = []
        if 'rank_cutoff' in user_options:
            x.extend(['--rank_cutoff', user_options['rank_cutoff']])
        if 'zscore_cutoff' in user_options:
            x.extend(['--zscore_cutoff', user_options['zscore_cutoff']])
        command = [
            'python', analyze_BIN, data_node.identifier,
            clinical_node.identifier, '--outcome', user_options['outcome'] +
            ',' + user_options['dead'], '--gene', user_options['genename'],
            '-o', 'clin'] + x
        process = subprocess.Popen(
            command, shell=False, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        outputfiles = os.listdir(".")
        os.mkdir(outfile)
        for i in outputfiles:
            shutil.copyfile(i, os.path.join(outfile, i))
        assert filelib.exists_nz(outfile), (
            'the output file %s for analyze_clinical_outcome fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, clinical_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'survial_analysis_' + original_file
        return filename
