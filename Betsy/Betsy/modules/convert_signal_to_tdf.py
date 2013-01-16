#convert_signal_to_tdf.py
import os
from Betsy import module_utils
import shutil
import xlrd 
import openpyxl
import arrayio
import subprocess
from genomicode import config, arrayplatformlib, arrayannot

def run(parameters, objects, pipeline):
    """check an input file is xls or xlsx format"""
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    if single_object.identifier.endswith('.gz'):
        newfile = module_utils.gunzip(single_object.identifier)
    else:
        newfile = single_object.identifier
    M = None
    tmp_file = None
    tmp1_file = newfile
    try:
        xlrd.open_workbook(newfile)
        tmp_file = 'tmp.xls'
    except Exception, XLRDError:
        try:
            book = openpyxl.load_workbook(newfile)
            tmp_file = 'tmp.xlsx'
        except Exception, InvalidFileException:
            tmp_file = None
        except (SystemError, MemoryError, KeyError), x:
            raise
    if tmp_file:
        shutil.copyfile(newfile, tmp_file)
        xls2txt_path = config.xls2txt
        xls2txt_BIN = module_utils.which(xls2txt_path)
        assert xls2txt_BIN, 'cannot find the %s' % xls2txt_path
        f = file('tmp1.txt', 'w')
        command = ['python', xls2txt_BIN, tmp_file]
        process = subprocess.Popen(command, shell=False,
                                   stdout=f,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        os.remove(tmp_file)
        f.close()
        tmp1_file = 'tmp1.txt'
    try:
        M = arrayio.choose_format(tmp1_file)
    except Exception, x:
        raise
    except (SystemError, MemoryError, KeyError), x:
        raise
    if M:
        M = check_gct_header(tmp1_file)
        M_c = arrayio.convert(M, to_format=arrayio.tab_delimited_format)
        f = file(outfile, 'w')
        arrayio.tab_delimited_format.write(M_c, f)
        f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_signal_to_tdf does not exists'
        % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(parameters, single_object,
                                             pipeline, outfile)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = module_utils.find_object(
        parameters, objects, 'input_signal_file', 'contents')
    if single_object:
        parameters = module_utils.renew_parameters(parameters, ['status'])
        attributes = parameters.values()
        new_object = module_utils.DataObject('signal_file',
                                             attributes, outfile)
        new_objects = objects[:]
        new_objects.append(new_object)
        return new_objects
    else:
        single_object = get_identifier(parameters, objects)
        return module_utils.get_newobjects(outfile, 'signal_file',
                                           parameters, objects, single_object)


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'input_signal_file', 'contents')
    if not single_object:
        single_object = module_utils.find_object(parameters, objects,
                                                 'signal_file',
                                                 'contents,preprocess')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for convert_signal_to_tdf does not exist'
        % single_object.identifier)
    return single_object


def check_gct_header(filename):
    M_name = arrayio.choose_format(filename)
    M = arrayio.read(filename)
    ids = M._row_order
    if M_name.__name__ == 'arrayio.gct_format':
        all_platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
        if all_platforms:
             probe_ids = M._row_names[all_platforms[0][0]]
             for platform in all_platforms:
                if 'Entrez' not in platform[1]:
                    old_header = platform[0]
                    new_header = 'Probe ID'
                    M, ids = module_utils.replace_matrix_header(
                        M, old_header, new_header)
                elif 'Entrez' in platform[1]:
                    old_header = platform[0]
                    new_header = 'Entrez ID'
                    M, ids = module_utils.replace_matrix_header(
                        M, old_header, new_header)
             annotate_header = arrayplatformlib.annotate_header
             dictionary = arrayannot.annotate_probes_multiple(
                 probe_ids, annotate_header)      
             old_ids = ids
             for id in old_ids:
                flag = True
                for key in dictionary.keys():
                    flag = True
                    value_list = dictionary[key]
                    gene_list = M._row_names[id]
                    for gene in gene_list:
                        if gene not in value_list:
                            flag = False
                            break
                    if flag:
                        M, ids = module_utils.replace_matrix_header(
                            M, id, key)
    return M
