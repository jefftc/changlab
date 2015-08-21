from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """convert the cel file with ccl or v3_4 to v3_4"""
        import os
        import shutil
        from genomicode import affyio
        from Betsy import module_utils
        in_data = antecedents
        filenames = os.listdir(in_data.identifier)
        if not os.path.exists(outfile):
            os.mkdir(outfile)
    
        #ver_list = []
        
        for filename in filenames:
            if filename == '.DS_Store':
                pass
            else:
                fileloc = os.path.join(in_data.identifier, filename)
                cel_v = affyio.guess_cel_version(fileloc)
                if fileloc.endswith('.gz'):
                    newcelfname = os.path.splitext(filename)[0]
                    cel_file = module_utils.gunzip(fileloc)
                else:
                    cel_file = fileloc
                    newcelfname = filename
                if cel_v == 'cc1':
                    f = file(os.path.join(outfile, newcelfname), 'w')
                    affyio.convert_cel_cc1_to_3(cel_file, f)
                    f.close()
                elif cel_v in ('v3', 'v4'):
                    shutil.copyfile(cel_file, os.path.join(outfile, newcelfname))
                if fileloc.endswith('.gz'):
                    os.remove(cel_file)
    
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for convert_CEL_to_v3 fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        return 'cel_files_' + original_file


