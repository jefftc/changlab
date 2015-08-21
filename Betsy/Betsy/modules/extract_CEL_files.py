from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """extract the cel files with cc or v3_4"""
        import os
        import shutil
        from Betsy import module_utils
        from genomicode import affyio
        in_data = antecedents
        directory = module_utils.unzip_if_zip(in_data.identifier)
        filenames = os.listdir(directory)
        assert filenames, 'The input folder or zip file is empty.'
        ver_list = []
        if not os.path.exists(outfile):
            os.mkdir(outfile)
    
        
        for filename in filenames:
            if filename == '.DS_Store':
                pass
            else:
                fileloc = os.path.join(directory, filename)
                cel_v = affyio.guess_cel_version(fileloc)
                if cel_v in ['cc1', 'v3', 'v4']:
                    shutil.copyfile(fileloc, os.path.join(outfile, filename))
                    ver_list.append(True)
                else:
                    ver_list.append(False)

    
        
        if True in ver_list:
            assert module_utils.exists_nz(outfile), (
                'the output file %s for extract_CEL_files fails' % outfile
            )
        else:
            assert ValueError('There is no cel file in the input.')



    
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'cel_files_' + original_file
        return filename



