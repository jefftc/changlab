from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import subprocess
        import shutil
        import arrayio
        from genomicode import config
        from genomicode import arrayplatformlib
        from genomicode import filelib
        from Betsy import module_utils
        
        in_data = antecedents
        DATA = arrayio.read(in_data.identifier)
        chipname = arrayplatformlib.identify_platform_of_matrix(DATA)
        platform = user_options['platform_name']
        assert arrayplatformlib.get_bm_attribute(platform), (
            'the desire platform %s is not recognized by Betsy' % platform
        )
        if chipname == platform:
            shutil.copyfile(in_data.identifier, outfile)
        else:
            Annot_path = config.annotate_matrix
            Annot_BIN = module_utils.which(Annot_path)
            assert Annot_BIN, 'cannot find the %s' % Annot_path
            command = ['python', Annot_BIN, in_data.identifier, "--platform",
                       platform]
            f = file(outfile, 'w')
            try:
                process = subprocess.Popen(command,
                                           shell=False,
                                           stdout=f,
                                           stderr=subprocess.PIPE)
            finally:
                f.close()
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)

        filelib.assert_exists_nz(outfile)



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_' + original_file + '.tdf'
        return filename



