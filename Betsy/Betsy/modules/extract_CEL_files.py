#extract_CEL_files.py
import os
from Betsy import module_utils
import shutil
import gzip

def run(parameters,objects,pipeline):
    """extract tge cel files with cc or v3_4"""
    from genomicode import affyio
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    filenames = os.listdir(single_object.identifier)
    ver_list = []
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    for filename in filenames:
        if filename == '.DS_Store':
            pass
        else:
            fileloc = os.path.join(single_object.identifier,filename)
            cel_v = affyio.guess_cel_version(fileloc)
            if cel_v in ['cc1','v3','v4']:
                shutil.copyfile(fileloc,os.path.join(outfile,filename))
                ver_list.append(True)
            else:
                ver_list.append(False)
    if True in ver_list:
        assert module_utils.exists_nz(outfile),(
            'the output file %s for extract_CEL_files fails'%outfile)
        new_objects = get_newobjects(parameters,objects,pipeline)
        module_utils.write_Betsy_parameters_file(
            parameters,single_object,pipeline,outfile)
        return new_objects
    else:
        return None
               

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

##    inputid = module_utils.get_inputid(identifier)
##    hash_profile={'version': 'cc_or_v3_4',
##                   'number of files':str(len(os.listdir(identifier))),
##                  'filenames':str(os.listdir(identifier))}
##    hash_result=hash_method.hash_parameters(
##                    inputid,pipeline,**hash_profile)
##    return hash_result

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'cel_files_' + original_file
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
    assert os.path.exists(single_object.identifier),(
        'folder %s for extract_CEl_files does not exit.'
        % single_object.identifier)
    assert os.path.isdir(single_object.identifier),"input is not a folder"
    return single_object
