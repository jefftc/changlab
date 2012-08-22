#download_geo_GSEID.py
import os
import rule_engine
from genomicode import affyio
import shutil
import gzip
import string
import module_utils
def run(parameters,objects,pipeline):
    """given a database ID and database,download and untar the data folder"""
    #download the tar folder from geo
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    file_folder = os.path.join(os.getcwd(),single_object.identifier)
    module_utils.download_dataset(single_object.identifier)
    #get chip name
    cel_files = os.listdir(file_folder)
    unknown_folder = os.path.join(os.getcwd(),'unknown_folder')
    chip_name_list = []
    for cel_file in cel_files:
            fileloc = os.path.join(file_folder,cel_file)
            if fileloc.endswith('.gz'):
                newcelfname = clean_cel_filename(os.path.splitext(cel_file)[0])
                #unzip the gz data
                unzipfile = os.path.splitext(fileloc)[0]
                fileObj = gzip.GzipFile(fileloc, 'rb');
                fileObjOut = file(unzipfile, 'wb');
                while 1:
                    line = fileObj.readline()
                    if line == '':
                        break
                    fileObjOut.write(line)
                fileObj.close()
                fileObjOut.close()
                assert os.path.exists(unzipfile),'the unzip %s fails' %unzipfile
            else:
                unzipfile = fileloc
                newcelfname = clean_cel_filename(cel_file)
            #get chip_name and copy into different folder
            chip_name = None
            try:
                chip_name = affyio.extract_chip_name(unzipfile)
            except (SystemError,MemoryError,KeyError),x:
                    raise 
            except Exception,x:
                if not os.path.exists(unknown_folder):
                    os.mkdir(unknown_folder)
                shutil.copyfile(unzipfile,
                                os.path.join(unknown_folder,newcelfname))
            if chip_name is not None:
                if chip_name not in chip_name_list:
                    chip_name_list.append(chip_name)
                    os.mkdir(os.path.join(os.getcwd(),chip_name))
                chip_folder = os.path.join(os.getcwd(),chip_name)
                shutil.copyfile(unzipfile,
                                os.path.join(chip_folder,newcelfname))
            if fileloc.endswith('.gz'):
                os.remove(unzipfile)
    #determine the one to preprocess
    if len(chip_name_list) == 1:
        out_filename = os.path.join(os.getcwd(),chip_name_list[0])
        
    elif len(chip_name_list) > 1:
        size_list = [os.path.getsize(os.path.join(os.getcwd(),x)) for x in chip_name_list]
        #makesure there is no multiple folder have the same maximum size  
        maxsize = max(size_list)
        new_size_list = size_list[:]
        new_size_list.remove(maxsize)
        #only one folder is maximum size
        if  maxsize>max(new_size_list):
            out_chip_name = chip_name_list[size_list.index(maxsize)]
            out_filenanme = os.path.join(os.getcwd(),out_chip_name)
            
        #multiple has same maximum size
        elif maxsize == max(new_size_list):
            start = -1
            folder_index = []
            while start<len(size_list)-1:
                try:
                    start = size_list.index(maxsize,start+1)
                    folder_index.append(start)
                except (SystemError,MemoryError,KeyError),x:
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
            if sum(a)== 1:
                out_chip_name = folder_names[a.index(1)]
                out_filename = os.path.join(os.getcwd(),out_chip_name)
            #multipld human paltforms
            elif sum(a)>1:
                if 'HG-U133_Plus_2' in folder_names:
                    out_filename = os.path.join(os.getcwd(),'HG-U133_Plus_2')
                elif 'HG-U133A' in folder_names:
                    out_filename = os.path.join(os.getcwd(),'HG-U133A')
                elif 'HG-U95A' in folder_names:
                    out_filename = os.path.join(os.getcwd(),'HG-U95A')
                else:
                    raise ValueError('does not recognazie the platform')
    os.rename(out_filename,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for download_geo_dataset fails' %outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    hash_string = identifier
    return hash_string

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    filename = 'after_select_' + single_object.identifier
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    newobjecttype = parameters['filetype']
    new_object = rule_engine.DataObject(newobjecttype,[
                        'unknown_version',parameters['contents']],outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'gse_id','contents')
    assert single_object.identifier.startswith('GSE'),(
        'the GSEID %s should start with GSE'%single_object.identifier)
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
    
    