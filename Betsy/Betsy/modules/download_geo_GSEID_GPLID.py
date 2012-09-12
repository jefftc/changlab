#download_geo_GSEID_GPLID.py
import os
import shutil
from ftplib import FTP
import module_utils
import string
import gzip

def get_seriesmatrix_file(GSEID,GPLID):
    'download series matrix and unzip'
    try:
        ftp = FTP('ftp.ncbi.nih.gov')
        ftp.login()
    except Exception,e:
        raise ValueError(e)
    try:
        ftp.cwd('pub/geo/DATA/SeriesMatrix/'+GSEID)
    except FTP.error_perm,x:
        if str(x).find('No such file')>=0:
            raise AssertionError,'cannot find the %s' %path
    entry = []
    ftp.retrlines('NLST',entry.append)
    platform_txtfiles = []
    for platform_filename in entry:
        if GPLID in platform_filename:
            f = open(platform_filename,'wb')
            ftp.retrbinary('RETR '+platform_filename,f.write)
            f.close()
            platform_txtfile = platform_filename[:-3]
            assert not os.path.exists(platform_txtfile),(
                'the seriesmatrix file %s already exists'%platform_txtfile)
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
            assert os.path.exists(platform_txtfile),(
                'the unzip %s in download_geo_dataset_GPL fails'
                %platform_txtfile)
            platform_txtfiles.append(platform_txtfile)
    ftp.close()
    return platform_txtfiles

def run(parameters,objects,pipeline):
    """given a database ID and GPLID, get the cel files"""
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    GSEID = single_object.identifier.split(',')[0]
    GPLID = single_object.identifier.split(',')[1]
    assert GSEID.startswith('GSE'),'GSEID %s is not correct'%GSEID
    assert GPLID.startswith('GPL'),'GPLID %s is not correct'%GPLID
    GSEID_path = module_utils.download_dataset(GSEID)
    platform_txtfiles = get_seriesmatrix_file(GSEID,GPLID)
    #get the cel file name for the GPL platform
    if not os.path.exists(outfile):
            os.mkdir(outfile)
    if len(platform_txtfiles)>0:
        for platform_txtfile in platform_txtfiles:
            cel_list = open(platform_txtfile,'r').readlines()
            cel_line = None
            for linecontent in cel_list:
                if linecontent.startswith('!Sample_geo_accession'):
                    cel_line = linecontent
                    break
            assert cel_line,(
                'the file %s does not contain "!Sample_geo_accession"'
                %platform_txtfile)
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
    else:
        os.rename(GSEID_path,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for download_geo_dataset_GPL fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects
         
def make_unique_hash(identifier,pipeline,parameters):
    return identifier.split(',')[0]+'.'+identifier.split(',')[1]
    
def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    hash_string = make_unique_hash(single_object.identifier,pipeline,parameters)
    filename = hash_string
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    newobjecttype = parameters['filetype']
    new_object = module_utils.DataObject(newobjecttype,[
                        'unknown_version',parameters['contents']],outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects


def get_identifier(parameters,objects):
    single_object = module_utils.find_object(parameters,objects,
                        'gse_id_and_platform','contents')
    assert len(single_object.identifier.split(','))==2,(
        'the identifier %s should contain GSEID and GPLID'
        %single_object.identifier)
    return single_object

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
    
