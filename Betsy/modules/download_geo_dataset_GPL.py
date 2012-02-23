#download_geo_dataset_GPL.py
import os
import tarfile
import shutil
import rule_engine
from ftplib import FTP
import module_utils
import string
def download_dataset(GSEID):
    #download the tar folder from geo
    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()
    GSE_directory = 'pub/geo/DATA/supplementary/series/'+GSEID
    ftp.cwd(GSE_directory)    
    download_rarfile = GSEID+'_RAW.tar'
    f = open(download_rarfile,'wb')
    ftp.retrbinary('RETR '+GSEID+'_RAW.tar',f.write)
    f.close()
    ftp.close()
    #untar the data folder
    GSEID_path = GSEID
    if not tarfile.is_tarfile(download_rarfile):
        raise ValueError('download file is not tar file')
    else:
        tar = tarfile.open(download_rarfile)
        tar.extractall(path=GSEID_path)
        tar.close()
    os.remove(download_rarfile)
    assert os.path.exists(GSEID_path),'the download file does not exist'
    assert len(os.listdir(GSEID_path))>0, 'the untar fails'
    return GSEID_path

def get_seriesmatrix_file(GSEID,GPLID):
    'download series matrix and unzip'
    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()
    GPL_directory = 'pub/geo/DATA/SeriesMatrix/'+GSEID
    ftp.cwd(GPL_directory)
    entry = []
    ftp.retrlines('NLST',entry.append)
    platform_txtfiles = []
    for platform_filename in entry:
        if GPLID in platform_filename:
            f = open(platform_filename,'wb')
            ftp.retrbinary('RETR '+platform_filename,f.write)
            f.close()
            platform_txtfile = platform_filename[:-3]
            assert not os.path.exists(platform_txtfile),'the seriesmatrix file already exists'
            #unzip the gz data
            import gzip
            fileObj = gzip.GzipFile(platform_filename, 'rb');
            fileObjOut = file(platform_txtfile, 'wb');
            while 1:
                line = fileObj.readline()
                if line == '':
                    break
                fileObjOut.write(line)
            fileObj.close()
            fileObjOut.close()
            os.remove(platform_filename)
            assert os.path.exists(platform_txtfile),'the unzip fails'
            platform_txtfiles.append(platform_txtfile)
    ftp.close()
    return platform_txtfiles

def run(parameters,objects,pipeline):
    """given a database ID and GPLID, get the cel files"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    assert len(identifier.split(','))==2,'the identifier should contain GSEID and GPLID'
    GSEID = identifier.split(',')[0]
    GPLID = identifier.split(',')[1]
    assert GSEID.startswith('GSE'),'GSEID is not correct'
    assert GPLID.startswith('GPL'),'GPLID is not correct'
    GSEID_path = download_dataset(GSEID)
    platform_txtfiles = get_seriesmatrix_file(GSEID,GPLID)
    #get the cel file name for the GPL platform
    if not os.path.exists(outfile):
            os.mkdir(outfile)
    for platform_txtfile in platform_txtfiles:
        cel_list = open(platform_txtfile,'r').readlines()
        cel_line = None
        for linecontent in cel_list:
            if linecontent.startswith('!Sample_geo_accession'):
                cel_line = linecontent
                break
        assert cel_line,'the file does not contain "!Sample_geo_accession"'
        filecontent = os.listdir(GSEID_path)
        cel_names = []
        for x in linecontent.split()[1:]:
            x = x.strip()
            assert x.startswith('\"') and x.endswith('\"')
            x = x[1:-1]
            cel_names.append(x)
        #check if the GSM Id cannot found in the data set
        file_name_string = ' '.join(filecontent)
        for cel_name in cel_names:
            if cel_name not in file_name_string:
                raise ValueError('The GSM ID %s cannot find in data set'%cel_name)
            else:
                for cel_file in filecontent:
                    if cel_file.upper().startswith(cel_name.upper()):
                        if cel_file.lower().endswith('gz'):
                            cel_file = clean_cel_filename(os.path.splitext(cel_file)[0])+'.gz'
                        outfilename=os.path.join(outfile,cel_file)
                        shutil.copyfile(os.path.join(GSEID_path,cel_file),outfilename)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
         
def make_unique_hash(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    return identifier.split(',')[0]+'.'+identifier.split(',')[1]
    
def get_outfile(parameters,objects):
    hash_string = make_unique_hash(parameters,objects)
    identifier,single_object = get_identifier(parameters,objects)
    filename = hash_string
    outfile = os.path.join(os.getcwd(),filename)
    new_object = rule_engine.DataObject('geo_dataset',[parameters['DatasetId'],
                                        'unknown',parameters['Contents']],outfile)
    new_objects = objects[:]
    new_objects.remove(single_object)
    new_objects.append(new_object)
    return outfile,new_objects

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,
                                    'gse_dataset_and_platform','Contents,DatasetId')

def clean_cel_filename(cel_file):
    """clean the cel_file name"""
    if cel_file.upper().endswith('CEL'):
        cel_file=cel_file[0:-4]
        punc = string.punctuation
        indicate = [x in cel_file for x in punc]
        if True in indicate:
            punct_index = []
            start=0
            while start<len(indicate)-1:
                try:
                    start = indicate.index(True,start+1)
                    punct_index.append(start)
                except (SystemError,MemoryError,KeyError),x:
                      raise
                except Exception:
                    break
            break_point=[cel_file.index(punc[x]) for x in punct_index]
            return cel_file[0:min(break_point)]+'.CEL'
        else:
            return cel_file+'.CEL'
    else:
        return cel_file
    
