#download_geo_supplement.py
import os
import shutil
from genomicode import affyio
from ftplib import FTP
from Betsy import module_utils,bie3,rulebase
import string
import gzip


def run(data_node, parameters, user_input, network, num_cores):
    """given a database ID and GPLID, get the cel files"""
    outfile = name_outfile(data_node,user_input)
    GSEID = user_input['GSEID']
    GPLID = None
    if 'GPLID' in user_input:
        GPLID = user_input['GPLID']
    assert GSEID.startswith('GSE'), 'GSEID %s is not correct' % GSEID
    if not GPLID:
        download_geo_with_GSEID(GSEID, outfile)
    else:
        assert GPLID.startswith('GPL'), 'GPLID %s is not correct' % GPLID
        download_geo_with_GPLID(GSEID, GPLID, outfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for download_geo_supplement fails' % outfile)
    out_node = bie3.Data(rulebase.ExpressionFiles,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(user_input['GSEID'])
    filename = 'expression_files_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(data_node,pipeline,parameters,user_input):
    GPLID = ''
    if  'GPLID' in user_input:
        GPLID = user_input['GPLID']
    identifier = user_input['GSEID'] + GPLID
    return identifier

def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node

def clean_cel_filename(cel_file):
    """clean the cel_file name"""
    if cel_file.upper().endswith('CEL'):
        cel_file = cel_file[0: -4]
        punc = string.punctuation
        indicate = [x in cel_file for x in punc]
        if True in indicate:
            punct_index = []
            start = 0
            while start < len(indicate) - 1:
                try:
                    start = indicate.index(True, start + 1)
                    punct_index.append(start)
                except (SystemError, MemoryError, KeyError), x:
                    raise
                except Exception:
                    break
            break_point = [cel_file.index(punc[x]) for x in punct_index]
            return cel_file[0: min(break_point)] + '.CEL'
        else:
            return cel_file + '.CEL'
    else:
        return cel_file

def download_geo_with_GSEID(GSEID, outfile):
    file_folder = os.path.join(os.getcwd(), GSEID)
    module_utils.download_dataset(GSEID)
    #get chip name
    cel_files = os.listdir(file_folder)
    unknown_folder = os.path.join(os.getcwd(), 'unknown_folder')
    chip_name_list = []
    for cel_file in cel_files:
            fileloc = os.path.join(file_folder, cel_file)
            if fileloc.endswith('.gz'):
                newcelfname = clean_cel_filename(os.path.splitext(cel_file)[0])
                #unzip the gz data
                unzipfile = os.path.splitext(fileloc)[0]
                fileObj = gzip.GzipFile(fileloc, 'rb')
                fileObjOut = file(unzipfile, 'wb')
                while 1:
                    line = fileObj.readline()
                    if line == '':
                        break
                    fileObjOut.write(line)
                fileObj.close()
                fileObjOut.close()
                assert os.path.exists(unzipfile), (
                    'the unzip %s fails' % unzipfile)
            else:
                unzipfile = fileloc
                newcelfname = clean_cel_filename(cel_file)
            #get chip_name and copy into different folder
            chip_name = None
            try:
                chip_name = affyio.extract_chip_name(unzipfile)
            except (SystemError, MemoryError, KeyError), x:
                raise
            except Exception, x:
                if not os.path.exists(unknown_folder):
                    os.mkdir(unknown_folder)
                shutil.copyfile(unzipfile,
                                os.path.join(unknown_folder, newcelfname))
            if chip_name is not None:
                if chip_name not in chip_name_list:
                    chip_name_list.append(chip_name)
                    os.mkdir(os.path.join(os.getcwd(), chip_name))
                chip_folder = os.path.join(os.getcwd(), chip_name)
                shutil.copyfile(unzipfile,
                                os.path.join(chip_folder, newcelfname))
            if fileloc.endswith('.gz'):
                os.remove(unzipfile)
    #determine the one to preprocess
    if len(chip_name_list) == 1:
        out_filename = os.path.join(os.getcwd(), chip_name_list[0])
    elif len(chip_name_list) > 1:
        size_list = [os.path.getsize(os.path.join(os.getcwd(), x))
                     for x in chip_name_list]
        #makesure there is no multiple folder have the same maximum size
        maxsize = max(size_list)
        new_size_list = size_list[:]
        new_size_list.remove(maxsize)
        #only one folder is maximum size
        if  maxsize > max(new_size_list):
            out_chip_name = chip_name_list[size_list.index(maxsize)]
            out_filenanme = os.path.join(os.getcwd(), out_chip_name)
        #multiple has same maximum size
        elif maxsize == max(new_size_list):
            start = -1
            folder_index = []
            while start < len(size_list) - 1:
                try:
                    start = size_list.index(maxsize, start + 1)
                    folder_index.append(start)
                except (SystemError, MemoryError, KeyError), x:
                    raise
                except Exception:
                    break
            folder_names = [chip_name_list[x] for x in folder_index]
            Is_HG = [x.startswith('HG') for x in folder_names]
            a = []
            for i in Is_HG:
                if i:
                    a.append(1)
                else:
                    a.append(0)
            #choose the human platform
            if sum(a) == 1:
                out_chip_name = folder_names[a.index(1)]
                out_filename = os.path.join(os.getcwd(), out_chip_name)
            #multipld human paltforms
            elif sum(a) > 1:
                if 'HG-U133_Plus_2' in folder_names:
                    out_filename = os.path.join(os.getcwd(), 'HG-U133_Plus_2')
                elif 'HG-U133A' in folder_names:
                    out_filename = os.path.join(os.getcwd(), 'HG-U133A')
                elif 'HG-U95A' in folder_names:
                    out_filename = os.path.join(os.getcwd(), 'HG-U95A')
                else:
                    raise ValueError('does not recognazie the platform')
    os.rename(out_filename, outfile)
    
def download_geo_with_GPLID(GSEID, GPLID, outfile):
    GSEID_path = module_utils.download_dataset(GSEID)
    platform_txtfiles = get_seriesmatrix_file(GSEID, GPLID)
    #get the cel file name for the GPL platform
    if not os.path.exists(outfile):
            os.mkdir(outfile)
    if len(platform_txtfiles) > 0:
        for platform_txtfile in platform_txtfiles:
            cel_list = open(platform_txtfile, 'r').readlines()
            cel_line = None
            for linecontent in cel_list:
                if linecontent.startswith('!Sample_geo_accession'):
                    cel_line = linecontent
                    break
            assert cel_line, (
                'the file %s does not contain "!Sample_geo_accession"'
                % platform_txtfile)
            filecontent = os.listdir(GSEID_path)
            cel_names = []
            for x in linecontent.split()[1:]:
                x = x.strip()
                assert x.startswith('\"') and x.endswith('\"')
                x = x[1: -1]
                cel_names.append(x)
            #check if the GSM Id cannot found in the data set
            file_name_string = ' '.join(filecontent)
            for cel_name in cel_names:
                if cel_name not in file_name_string:
                    raise ValueError('The GSM ID %s cannot find in data set'
                                     % cel_name)
                else:
                    for cel_file in filecontent:
                        if cel_file.upper().startswith(cel_name.upper()):
                            if cel_file.lower().endswith('gz'):
                                cel_file = clean_cel_filename(
                                    os.path.splitext(cel_file)[0]) + '.gz'
                            outfilename = os.path.join(outfile, cel_file)
                            shutil.copyfile(
                                os.path.join(GSEID_path, cel_file),
                                outfilename)
    else:
        os.rename(GSEID_path, outfile)

def get_seriesmatrix_file(GSEID, GPLID):
    'download series matrix and unzip'
    try:
        ftp = FTP('ftp.ncbi.nih.gov')
        ftp.login()
    except Exception, e:
        raise ValueError(e)
    try:
        ftp.cwd('pub/geo/DATA/SeriesMatrix/' + GSEID)
    except FTP.error_perm, x:
        if str(x).find('No such file') >= 0:
            raise AssertionError('cannot find the %s' % path)
    entry = []
    ftp.retrlines('NLST', entry.append)
    platform_txtfiles = []
    for platform_filename in entry:
        if GPLID in platform_filename:
            f = open(platform_filename, 'wb')
            ftp.retrbinary('RETR ' + platform_filename, f.write)
            f.close()
            platform_txtfile = platform_filename[: -3]
            assert not os.path.exists(platform_txtfile), (
                'the seriesmatrix file %s already exists' % platform_txtfile)
            #unzip the gz data
            import gzip
            fileObj = gzip.GzipFile(platform_filename, 'rb')
            fileObjOut = file(platform_txtfile, 'wb')
            while 1:
                line = fileObj.readline()
                if line == '':
                    break
                fileObjOut.write(line)
            fileObj.close()
            fileObjOut.close()
            os.remove(platform_filename)
            assert os.path.exists(platform_txtfile), (
                'the unzip %s in download_geo_dataset_GPL fails'
                % platform_txtfile)
            platform_txtfiles.append(platform_txtfile)
    ftp.close()
    return platform_txtfiles
