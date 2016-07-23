from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """convert the cel file with ccl or v3_4 to v3_4"""
        import shutil
        from genomicode import filelib
        from Betsy import module_utils
        in_data = antecedents
        #new_parameters = set_out_attributes(in_data, out_attributes)
        shutil.copytree(in_data.identifier, outfile)
        assert filelib.exists_nz(outfile), (
            'the output file %s for detect_CEL_version' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'cel_files_' + original_file
        return filename


    def set_out_attributes(self, antecedents, out_attributes):
        import os
        from genomicode import affyio
        filenames = os.listdir(antecedents.identifier)
        ver_list = []
        for filename in filenames:
            if filename == '.DS_Store':
                pass
            else:
                fileloc = os.path.join(antecedents.identifier, filename)
                cel_v = affyio.guess_cel_version(fileloc)
                ver_list.append(cel_v)
        
        new_parameters = out_attributes.copy()
        if 'cc1' in ver_list:
            new_parameters['version'] = 'cc'
        elif 'v3' in ver_list:
            new_parameters['version'] = 'v3_v4'
        elif 'v4' in ver_list:
            new_parameters['version'] = 'v3_v4'
        else:
            raise ValueError('the cel file can only be cc,v3,v4')
        
        return new_parameters



