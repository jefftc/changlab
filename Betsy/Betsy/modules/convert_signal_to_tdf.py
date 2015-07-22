#convert_signal_to_tdf.py
import os
from Betsy import module_utils, userfile
import shutil
import xlrd
import openpyxl
import arrayio
import subprocess
from genomicode import config, arrayplatformlib
from Betsy import bie3
from Betsy import rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """check an input file is xls or xlsx format"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    real_name = in_data.identifier
    try:
        x = userfile._unhash_storefile(in_data.identifier)
        real_name = x[1]
    except:
        pass
    
    if (in_data.identifier.endswith('.gz') or real_name.endswith('.gz')):
        unzip_file = module_utils.gunzip(in_data.identifier)
    else:
        unzip_file = in_data.identifier
    
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
        command = ['python', xls2txt_BIN, xls_file]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=f,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        os.remove(xls_file)
        f.close()
        txt_file = 'tmp1.txt'
    
    M = guess_and_change_gct_header(txt_file)
    M_c = arrayio.convert(M, to_format=arrayio.tab_delimited_format)
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_signal_to_tdf does not exists' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)

    return data_node


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


platform2header = {
    'Agilent_Human1A': 'Probe ID',
    'HG_U133A_2': 'Probe ID',
    'HG_U133A': 'Probe ID',
    'HG_U133B': 'Probe ID',
    'HG_U133_Plus_2': 'Probe ID',
    'HG_U95A': 'Probe ID',
    'HG_U95Av2': 'Probe ID',
    'Hu35KsubA': 'Probe ID',
    'Hu35KsubB': 'Probe ID',
    'Hu35KsubC': 'Probe ID',
    'Hu35KsubD': 'Probe ID',
    'Hu6800': 'Probe ID',
    'HumanHT_12_control': 'Probe ID',
    'HumanHT_12': 'Probe ID',
    'MG_U74Av2': 'Probe ID',
    'MG_U74Bv2': 'Probe ID',
    'MG_U74Cv2': 'Probe ID',
    'Mouse430_2': 'Probe ID',
    'Mouse430A_2': 'Probe ID',
    'MouseRef_8_control': 'Probe ID',
    'MouseRef_8': 'Probe ID',
    'Mu11KsubA': 'Probe ID',
    'Mu11KsubB': 'Probe ID',
    'RAE230A': 'Probe ID',
    'RG_U34A': 'Probe ID',
    'Entrez_ID_human': 'Gene ID',
    'Entrez_ID_mouse': 'Gene ID',
    'Entrez_symbol_human': 'Gene symbol',
    'Entrez_symbol_mouse': 'Gene symbol'
}
