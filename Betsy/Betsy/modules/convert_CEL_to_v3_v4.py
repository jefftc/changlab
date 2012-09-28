#convert_CEL_to_v3_v4.py
import os
from Betsy import module_utils
import shutil
import gzip


def run(parameters,objects,pipeline):
    """convert the cel file with ccl or v3_4 to v3_4"""
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
            if fileloc.endswith('.gz'):
                newcelfname =  os.path.splitext(filename)[0]
                cel_file = module_utils.gunzip(fileloc)
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
    assert module_utils.exists_nz(outfile),(
        'the output file %s for convert_CEL_to_v3_v4 fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects
    

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    hash_result = make_unique_hash(single_object.identifier,pipeline,parameters)
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
        'folder %s for convert_CEL_to_v3_v4 does not exit.'
        % single_object.identifier)
    assert os.path.isdir(single_object.identifier),"input is not a folder"
    return single_object
