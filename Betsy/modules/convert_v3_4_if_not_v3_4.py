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
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    filenames = os.listdir(single_object.identifier)
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    ver_list = []
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(single_object.identifier,filename)
            cel_v = affyio.guess_cel_version(fileloc)
            if cel_v in ['cc1','v3','v4']:
                ver_list.append(True)
            else:
                ver_list.append(False)
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
    if True in ver_list:
        assert module_utils.exists_nz(outfile),'the output \
                              file %s for convert_v3_4_if_not_v3_4 fails'%outfile
        new_objects = get_newobjects(parameters,objects,pipeline)
        module_utils.write_Betsy_parameters_file(
            parameters,single_object,pipeline,outfile)
        return new_objects
    else:
        return None

def make_unique_hash(identifier,pipeline,parameters):
    inputid = module_utils.get_inputid(identifier)
    hash_profile={'version': 'v3_4',
                   'number of files':str(len(os.listdir(identifier))),
                  'filenames':str(os.listdir(identifier))}
    hash_result=hash_method.hash_parameters(
                    inputid,pipeline,**hash_profile)
    return hash_result

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    hash_result = make_unique_hash(single_object.identifier,pipeline,parameters)
    filename = 'cel_files_' + hash_result + '*'+ original_file
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'cel_files',parameters,objects,single_object)
    return new_objects

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'cel_files','contents')
    assert os.path.exists(single_object.identifier),'folder %s \
           for convert_v3_4_if_not_v3_4 does not exit.' % single_object.identifier
    assert os.path.isdir(single_object.identifier),"input is not a folder"
    return single_object
