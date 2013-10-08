#convert_signal_to_tdf.py
import os
#from Betsy import module_utils
import module_utils
import shutil
import xlrd 
import openpyxl
import arrayio
import subprocess
from genomicode import config, arrayplatformlib
import bie
import rulebase
        
def run(data_node,parameters):
    """check an input file is xls or xlsx format"""
    starttime = strftime(module_utils.FMT, localtime())
    outfile = name_outfile(data_node)
    if data_node.attributes['filename'].endswith('.gz'):
        unzip_file = module_utils.gunzip(data_node.attributes['filename'])
    else:
        unzip_file = data_node.attributes['filename']
    M = None
    xls_file = None
    txt_file = unzip_file
    try:
        xlrd.open_workbook(unzip_file)
        xls_file = 'tmp.xls'
    except Exception, XLRDError:
        try:
            book = openpyxl.load_workbook(unzip_file)
            xls_file = 'tmp.xlsx'
        except Exception, InvalidFileException:
            xls_file = None
        except (SystemError, MemoryError, KeyError), x:
            raise
    if xls_file:
        shutil.copyfile(unzip_file, xls_file)
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
        txt_file = 'tmp1.txt'
    M = arrayio.choose_format(txt_file)
    M = guess_and_change_gct_header(txt_file)
    M_c = arrayio.convert(M, to_format=arrayio.tab_delimited_format)
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_signal_to_tdf does not exists'
        % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node




def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'signal_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def get_out_attributes(parameters,data_node):
    return parameters


def guess_and_change_gct_header(filename):
    M_name = arrayio.choose_format(filename)
    M = arrayio.read(filename)
    ids = M._row_order
    if not M_name.__name__ == 'arrayio.gct_format':
        return M
    all_platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
    if not all_platforms:
        return M
    old_header = all_platforms[0][0]
    platform = all_platforms[0][1]
    new_header = platform2header[platform]
    M = module_utils.replace_matrix_header(M, old_header, new_header)
    return M

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)

    return data_node

platform2header = {'Agilent_Human1A':'Probe ID',
                   'HG_U133A_2':'Probe ID',
                   'HG_U133A':'Probe ID',
                   'HG_U133B':'Probe ID',
                   'HG_U133_Plus_2':'Probe ID',
                   'HG_U95A':'Probe ID',
                   'HG_U95Av2':'Probe ID',
                   'Hu35KsubA':'Probe ID',
                   'Hu35KsubB':'Probe ID',
                   'Hu35KsubC':'Probe ID',
                   'Hu35KsubD':'Probe ID',
                   'Hu6800':'Probe ID',
                   'HumanHT_12_control':'Probe ID',
                   'HumanHT_12':'Probe ID',
                   'MG_U74Av2':'Probe ID',
                   'MG_U74Bv2':'Probe ID',
                   'MG_U74Cv2':'Probe ID',
                   'Mouse430_2':'Probe ID',
                   'Mouse430A_2':'Probe ID',
                   'MouseRef_8_control':'Probe ID',
                   'MouseRef_8':'Probe ID',
                   'Mu11KsubA':'Probe ID',
                   'Mu11KsubB':'Probe ID',
                   'RAE230A':'Probe ID',
                   'RG_U34A':'Probe ID',
                   'Entrez_ID_human':'Gene ID',
                   'Entrez_ID_mouse':'Gene ID',
                   'Entrez_symbol_human':'Gene symbol',
                   'Entrez_symbol_mouse':'Gene symbol'}
