from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import shutil
        from Betsy import module_utils
        in_data = antecedents
        result_files = os.listdir(in_data.identifier)
        for result_file in result_files:
            if '-controls' not in result_file:
                goal_file = os.path.join(in_data.identifier, result_file)
                shutil.copyfile(goal_file, outfile)
    
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for illu_signal fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_illumina_' + original_file + '.gct'
        return filename



