from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """select tumor sample only """
        import os
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        infile = in_data.identifier
        if in_data.identifier.endswith('tar.gz'):
            infile = extract_files(in_data.identifier)
    
        
        slice_matrix_BIN = config.slice_matrix
        slice_matrix = module_utils.which(slice_matrix_BIN)
        assert slice_matrix, 'cannot find the %s' % slice_matrix_BIN
        tempfile = 'temp.txt'

        process = subprocess.Popen([slice_matrix_BIN,
                                    '--reorder_col_alphabetical',
                                    infile, ],
                                   shell=False,
                                   stdout=file(tempfile, 'w'),
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        assert module_utils.exists_nz(tempfile), (
            'the temp file %s for select_tumor_only fails' % tempfile
        )
        tempfile1 = 'temp1.txt'
        process = subprocess.Popen([slice_matrix_BIN,
                                    '--tcga_solid_tumor_only',
                                    tempfile, ],
                                   shell=False,
                                   stdout=file(tempfile1, 'w'),
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        assert module_utils.exists_nz(tempfile1), (
            'the temp file %s for select_tumor_only fails' % tempfile1
        )
        tempfile2 = 'temp2.txt'
        process = subprocess.Popen([slice_matrix_BIN,
                                    '--tcga_relabel_patient_barcodes',
                                    tempfile1, ],
                                   shell=False,
                                   stdout=file(tempfile2, 'w'),
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        assert module_utils.exists_nz(tempfile2), (
            'the output file %s for select_tumor_only fails' % tempfile2
        )
        process = subprocess.Popen([slice_matrix_BIN,
                                    '--remove_duplicate_cols',
                                    tempfile2, ],
                                   shell=False,
                                   stdout=file(outfile, 'w'),
                                   stderr=subprocess.PIPE)

        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)

    
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for select_tumor_only fails' % outfile
        )

        os.remove(tempfile)
        os.remove(tempfile1)
        os.remove(tempfile2)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        original_file = original_file.replace('.rar', '')
        filename = 'TCGAFile_' + original_file + '.txt'
        return filename


def extract_files(gzfile):
    import os
    import tarfile
    
    assert gzfile.lower().endswith('tar.gz')
    tfile = tarfile.open(gzfile, 'r:gz')
    gzname = os.path.split(gzfile)[-1]
    newdir = os.path.join(".", gzname[:-7])
    tfile.extractall(newdir)
    folder = os.listdir(newdir)
    directory = os.path.join(newdir, folder[0])
    assert os.path.exists(directory)
    files = os.listdir(directory)
    for filename in files:
        if filename.lower().endswith('.txt') and filename != 'MANIFEST.txt':
            return os.path.join(directory, filename)
    return None
