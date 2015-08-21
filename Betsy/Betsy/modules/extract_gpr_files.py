from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outpath):
        """extract the files that are gpr format"""
        import os
        import shutil
        from Betsy import gpr_module
        from Betsy import module_utils

        directory = module_utils.unzip_if_zip(antecedents.identifier)
        x = os.listdir(directory)
        x = [x for x in x if x != ".DS_Store"]
        files = x
        assert files, 'The input folder or zip file is empty.'

        files = [
            x for x in files
            if gpr_module.check_gpr(os.path.join(directory, x))]
        assert files, 'There are no gpr files in the input.'
        
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        for file in files:
            x1 = os.path.join(directory, file)
            x2 = os.path.join(outpath, file)
            shutil.copyfile(x1, x2)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'gpr_files_' + original_file
        return filename
