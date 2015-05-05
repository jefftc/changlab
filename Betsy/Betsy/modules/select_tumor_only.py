#select_tumor_only.py
import os
import subprocess
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils
from genomicode import config


def extract_files(gzfile):
    assert gzfile.lower().endswith('tar.gz')
    import tarfile
    tfile = tarfile.open(gzfile, 'r:gz')
    gzname = os.path.split(gzfile)[-1]
    newdir = os.path.join(os.getcwd(), gzname[:-7])
    tfile.extractall(newdir)
    folder = os.listdir(newdir)
    directory = os.path.join(newdir, folder[0])
    assert os.path.exists(directory)
    files = os.listdir(directory)
    for filename in files:
        if filename.lower().endswith('.txt') and filename != 'MANIFEST.txt':
            return os.path.join(directory, filename)
    return None


def run(network, antecedents, out_attributes, user_options, num_cores):
    """select tumor sample only """
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
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

    out_node = bie3.Data(rulebase.TCGAFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    original_file = original_file.replace('.rar', '')
    filename = 'TCGAFile_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node
