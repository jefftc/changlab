from Module import AbstractModule

# TODO: Have attribute that indicates what kind of array data to pull
# out of this folder.  If multiple types of files are given, then this
# needs to be set.

# Multiple platforms detected (cel, gpr).  Please specify desired
# platform with the --mattr option.


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
        out_path = outfile
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        
        in_path = module_utils.unzip_if_zip(in_data.identifier)
        assert in_path == in_data.identifier
        filenames = os.listdir(in_path)
        assert filenames, "The input folder or zip file is empty."
    
        x = guess_datatype(in_path)
        datatype, filenames = x
        for in_filename in filenames:
            in_path, in_file = os.path.split(in_filename)
            out_filename = os.path.join(out_path, in_file)
            shutil.copyfile(in_filename, out_filename)

    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'expression_' + original_file
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        in_data = antecedents
        attrs = out_attributes.copy()
        datatype, filenames = guess_datatype(in_data.identifier)
        attrs['filetype'] = datatype
        return attrs


def guess_datatype(path):
    import os
    
    # TODO: What is folder contains files with multiple data types?
    assert os.path.isdir(path)

    # Make a list of all the files in the directory.
    filenames = []
    for x in os.walk(path):
        dirpath, dirnames, files = x
        x = [os.path.join(dirpath, x) for x in files]
        filenames.extend(x)

    typed_files = []  # list of (filename, type)
    for filename in filenames:
        stem, ext = os.path.splitext(filename)
        uext = ext.upper()

        filetype = None
        if uext == ".CEL":
            filetype = "cel"
        elif uext == ".IDAT":
            filetype = "idat"
        elif uext == ".GPR":
            filetype = "gpr"
        elif is_agilent_file(filename):
            filetype = "agilent"
        if not filetype:
            continue
        x = filename, filetype
        typed_files.append(x)

    assert typed_files, "No known microarray file types: %s" % path

    type2files = {}
    for filename, ftype in typed_files:
        if ftype not in type2files:
            type2files[ftype] = []
        type2files[ftype].append(filename)

    if 'cel' in type2files:
        return 'cel', type2files['cel']
    elif 'idat' in type2files:
        return 'idat', type2files['idat']
    elif 'gpr' in type2files:
        return 'gpr', type2files['gpr']
    elif 'agilent' in type2files:
        return 'agilent', type2files['agilent']
    # Should not get here.
    raise AssertionError



    ## def guess_datatype(folder, matrix_folder):
    ##     # TODO: What is folder contains files with multiple data types?
    ##     directory = module_utils.unzip_if_zip(folder)
    ##     filenames = os.listdir(directory)
    ##     assert filenames, 'The input folder or zip file is empty.'
    ##     result_files = dict()
    ##     for filename in filenames:
    ##         if '.CEL' or '.cel' in filename:
    ##             if 'cel' not in result_files:
    ##                 result_files['cel'] = []
    ##             result_files['cel'].append(filename)
    ##         elif '.idat' or '.IDAT' in filename:
    ##             if 'idat' not in result_files:
    ##                 result_files['idat'] = []
    ##             result_files['idat'].append(filename)
    ##         elif is_agilent_file(os.path.join(folder, filename)):
    ##             if 'agilent' or 'AGILENT' not in result_files:
    ##                 result_files['agilent'] = []
    ##             result_files['agilent'].append(filename)
    ##         elif '.gpr' or '.gpr' in filename:
    ##             if 'gpr' not in result_files:
    ##                 result_files['gpr'] = []
    ##             result_files['gpr'].append(filename)
    ##     if not result_files:
    ##         matrix_files = os.listdir(matrix_folder)
    ##         for filename in matrix_files:
    ##             if 'series_matrix.txt' in filename:
    ##                 result_files['matrix'] = [filename]
    ##     if not result_files:
    ##         raise ValueError(
    ##             'we cannot guess the datatype in the folder %s' % folder)
    ##     if 'cel' in result_files:
    ##         return 'cel', result_files['cel']
    ##     elif 'idat' in result_files:
    ##         return 'idat', result_files['idat']
    ##     elif 'gpr' in result_files:
    ##         return 'gpr', result_files['gpr']
    ##     elif 'agilent' in result_files:
    ##         return 'agilent', result_files['agilent']
    ##     elif 'matrix' in result_files:
    ##         return 'matrix', result_files['matrix']
    ##     else:
    ##         return None


def is_agilent_file(filename):
    # Return a boolean indicating whether this is an agilent file.

    # TODO: Figure out this code is supposed to do and test it.
    postag = []
    fline = []
    f = open(filename, 'r')
    for i in range(10):
        # Bug: What if nothing found in first 10 lines?
        line = f.readline()
        words = line.split()
        if len(words) > 0:
            postag.append(words[0])
            if words[0] == 'FEATURES':
                fline = set(words)
    
    f.close()
    signal_tag = set(['gProcessedSignal', 'rProcessedSignal'])
    if signal_tag.issubset(fline):
        if postag == ['TYPE', 'FEPARAMS', 'DATA', '*', 'TYPE', 'STATS', 'DATA',
                      '*', 'TYPE', 'FEATURES']:
            return True
    
    return False
