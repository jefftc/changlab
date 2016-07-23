from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """extract the matrix file from the expression files"""
        import os
        from Betsy import module_utils
        from genomicode import filelib
        #from genomicode import affyio
        in_data = antecedents
        directory = module_utils.unzip_if_zip(in_data.identifier)
        filenames = os.listdir(directory)
        assert filenames, 'The input folder or zip file is empty.'
        for filename in filenames:
            if 'series_matrix.txt' in filename:
                fileloc = os.path.join(directory, filename)
                outname = os.path.join(outfile, filename)
                extract_expression_file(fileloc, outname)
    
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for extract_matrix_file fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'UnprocessedSignalFile_' + original_file
        return filename



def extract_expression_file(fileloc, outname):
    """extract the matrix from the series_matrix file"""
    text_list = open(fileloc, 'r').readlines()
    start = text_list.index('!series_matrix_table_begin')
    end = text_list.index('!series_matrix_table_end')
    f = file(outname, 'w')
    for i in range(start + 1, end):
        line = text_list[i].replace('"', '')
        f.write(line)
    
    f.close()
