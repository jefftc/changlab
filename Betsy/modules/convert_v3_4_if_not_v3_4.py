#convert_v3_4_if_not_V3_4.py
import os
import module_utils
import rule_engine
import shutil
import gzip
import module_utils
import hash_method
def run(parameters,objects,pipeline):
    """check if all the cel file are v3_4"""
    from genomicode import affyio
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    filenames = os.listdir(identifier)
    os.mkdir(outfile)
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(identifier,filename)
            cel_v = affyio.guess_cel_version(fileloc)
            if fileloc.endswith('.gz'):
                newcelfname =  os.path.splitext(filename)[0]
                #unzip the gz data
                cel_file = fileloc[:-3]
                fileObj = gzip.GzipFile(fileloc, 'rb');
                fileObjOut = file(cel_file, 'wb');
                while 1:
                    line = fileObj.readline()
                    if line == '':
                        break
                    fileObjOut.write(line)
                fileObj.close()
                fileObjOut.close()
                assert os.path.exists(cel_file),'unzip the cel_file %s fails'%cel_file
            else:
                    cel_file = fileloc
                    newcelfname = filename
            if cel_v == 'cc1':
                f = file(os.path.join(outfile,newcelfname),'w')
                affyio.convert_cel_cc1_to_3(cel_file,f)
                f.close()
            elif cel_v in ('v3','v4'):
                shutil.copyfile(cel_file,os.path.join(outfile,newcelfname))
            if fileloc.endswith('.gz'):
                os.remove(cel_file)
    assert module_utils.exists_nz(outfile),'the output \
                          file %s for convert_v3_4_if_not_v3_4 fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects
   

def make_unique_hash(identifier,pipeline,parameters):
    inputid = module_utils.get_inputid(identifier)
    hash_profile={'version': 'v3_4',
                   'number of files':str(len(os.listdir(identifier)))}
    hash_result=hash_method.hash_parameters(
                    inputid,pipeline,**hash_profile)
    return hash_result

def get_outfile(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(identifier)
    hash_string = make_unique_hash(identifier,pipeline,parameters)
    filename = original_file + '_BETSYHASH1_' + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'cel_files',parameters,objects,single_object)
    return new_objects

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'cel_files','contents')
    assert os.path.exists(identifier),'folder %s \
           for convert_v3_4_if_not_v3_4 does not exit.' % identifier
    assert os.path.isdir(identifier),"input is not a folder"
    return identifier,single_object