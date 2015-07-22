import os
import shutil
from Betsy import module_utils, bie3, rulebase

# TODO: rename all check_ modules.

def run(network, antecedents, out_attributes, user_options, num_cores):
    """check the data type from the expression file"""
    #from genomicode import affyio
    data_node, matrix_node = antecedents
    outfile = name_outfile(antecedents, user_options)
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    directory = module_utils.unzip_if_zip(data_node.identifier)
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    x = guess_datatype(data_node.identifier, matrix_node.identifier)
    datatype, filenames = x
    for filename in filenames:
        fileloc = os.path.join(directory, filename)
        shutil.copyfile(fileloc, os.path.join(outfile, filename))
    assert module_utils.exists_nz(outfile), (
        'the output file %s for check_geo_file_type fails' % outfile
    )
    out_node = bie3.Data(rulebase.ExpressionFiles, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    data_node, matrix_node = antecedents
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'expression_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    data_node, matrix_node = antecedents
    new_parameters = out_attributes.copy()
    datatype, filenames = guess_datatype(data_node.identifier,
                                         matrix_node.identifier)
    new_parameters['filetype'] = datatype
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    data_node, matrix_node = antecedents
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    filter1 = module_utils.AntecedentFilter(datatype_name='ExpressionFiles')
    filter2 = module_utils.AntecedentFilter(
        datatype_name='GEOSeriesMatrixFile')
    x = module_utils.find_antecedents(
        network, module_id, user_attributes, pool, filter1, filter2)
    return x


def guess_datatype(folder, matrix_folder):
    #should be consider a folder contain multiple datatype case in the future
    "guess the datatype of the files in the input folder"
    directory = module_utils.unzip_if_zip(folder)
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    result_files = dict()
    for filename in filenames:
        if '.CEL' or '.cel' in filename:
            if 'cel' not in result_files:
                result_files['cel'] = []
            result_files['cel'].append(filename)
        elif '.idat' or '.IDAT' in filename:
            if 'idat' not in result_files:
                result_files['idat'] = []
            result_files['idat'].append(filename)
        elif is_agilent_file(os.path.join(folder, filename)):
            if 'agilent' or 'AGILENT' not in result_files:
                result_files['agilent'] = []
            result_files['agilent'].append(filename)
        elif '.gpr' or '.gpr' in filename:
            if 'gpr' not in result_files:
                result_files['gpr'] = []
            result_files['gpr'].append(filename)
    if not result_files:
        matrix_files = os.listdir(matrix_folder)
        for filename in matrix_files:
            if 'series_matrix.txt' in filename:
                result_files['matrix'] = [filename]
    if not result_files:
        raise ValueError(
            'we cannot guess the datatype in the folder %s' % folder)
    if 'cel' in result_files:
        return 'cel', result_files['cel']
    elif 'idat' in result_files:
        return 'idat', result_files['idat']
    elif 'gpr' in result_files:
        return 'gpr', result_files['gpr']
    elif 'agilent' in result_files:
        return 'agilent', result_files['agilent']
    elif 'matrix' in result_files:
        return 'matrix', result_files['matrix']
    else:
        return None


def is_agilent_file(filename):
    "test a file is a agilent or not"
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
