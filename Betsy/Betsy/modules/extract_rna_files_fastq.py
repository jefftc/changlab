from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """extract the fastq rna seq files"""
        import os
        import shutil
        from Betsy import module_utils
        in_data = antecedents
        directory = module_utils.unzip_if_zip(in_data.identifier)
        filenames = os.listdir(directory)
        assert filenames, 'The input folder or zip file is empty.'
        if not os.path.exists(outfile):
            os.mkdir(outfile)
    
        
        format_types = ['fa', 'fastq']
        for format_type in format_types:
            for filename in filenames:
                if filename == '.DS_Store':
                    continue
                fileloc = os.path.join(in_data.identifier, filename)
                if fileloc.endswith(format_type + '.gz'):
                    newfname = os.path.splitext(filename)[0]
                    new_file = module_utils.gunzip(fileloc)
                elif fileloc.endswith(format_type):
                    new_file = fileloc
                    newfname = filename
                    shutil.copyfile(new_file, os.path.join(outfile, newfname))
                if fileloc.endswith('.gz'):
                    os.remove(new_file)
    
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_rna_files_fastq fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'fastq_files_' + original_file
        return filename
