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
        from genomicode import filelib
        in_data = antecedents
        directory = module_utils.unzip_if_zip(in_data.identifier)
        agilent_files = []
        filenames = os.listdir(directory)
        assert filenames, 'The input folder or zip file is empty.'
        for filename in filenames:
            if filename in ['.DS_Store', '._.DS_Store', '.Rapp.history']:
                continue
            if os.path.isdir(os.path.join(directory, filename)):
                continue
            postag = []
            fline = []
            f = open(os.path.join(directory, filename), 'r')
            for i in range(10):
                line = f.readline()
                words = line.split()
                if len(words) > 0:
                    postag.append(words[0])
                    if words[0] == 'FEATURES':
                        fline = set(words)
            f.close()
            signal_tag = set(['gProcessedSignal', 'rProcessedSignal'])
            if signal_tag.issubset(fline):
                if postag == ['TYPE', 'FEPARAMS', 'DATA', '*', 'TYPE', 'STATS',
                              'DATA', '*', 'TYPE', 'FEATURES']:
                    agilent_files.append(filename)

    
        
        if agilent_files:
            if not os.path.exists(outfile):
                os.mkdir(outfile)
            for filename in agilent_files:
                old_file = os.path.join(directory, filename)
                new_file = os.path.join(outfile, filename)
                shutil.copyfile(old_file, new_file)
            assert filelib.exists_nz(outfile), (
                'the output file %s for extract_agilent_files fails' % outfile
            )
        else:
            raise ValueError('There is no agilent file in the input.')




    
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'agilent_files' + original_file
        return filename


