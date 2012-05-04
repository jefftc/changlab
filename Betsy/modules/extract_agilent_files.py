#extract_agilent_files.py

import module_utils
import shutil
import os
import hash_method
import rule_engine
import zipfile

def extract_all(zipName):
    z = zipfile.ZipFile(zipName)
    for f in z.namelist():
        if f.endswith('/'):
            os.makedirs(f)
        else:
            z.extract(f)

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    if zipfile.is_zipfile(identifier):
        directory = os.path.split(identifier)[-1]
        directory = os.path.splitext(directory)[0]
        directory = os.path.join(os.getcwd(),directory)
        extract_all(identifier)
    else:
        directory = identifier
    agilent_files = []
    for filename in os.listdir(directory):
        if filename in ['.DS_Store','._.DS_Store','.Rapp.history']:
            continue
        postag = []
        fline =[]
        f = open(os.path.join(directory,filename),'r')
        for i in range(10):
            line = f.readline()
            words =line.split()
            postag.append(words[0])
            if words[0] == 'FEATURES':
                fline = set(words)
        f.close()
        signal_tag =set(['gMeanSignal','rMeanSignal','gMedianSignal','rMedianSignal'])
        if (postag == ['TYPE', 'FEPARAMS', 'DATA', '*',
                      'TYPE', 'STATS', 'DATA', '*', 'TYPE', 'FEATURES'] and
                    signal_tag.issubset(fline)):
           agilent_files.append(filename)
    if agilent_files:
        os.mkdir(outfile)
        for filename in agilent_files: 
            old_file = os.path.join(directory,filename)
            new_file = os.path.join(outfile,filename)
            shutil.copyfile(old_file,new_file)
        assert module_utils.exists_nz(outfile),'the output \
                        file %s for extract_agilent_files fails'%outfile
        new_objects = get_newobjects(parameters,objects,pipeline)
        module_utils.write_Betsy_parameters_file(parameters,
                                                 single_object,pipeline)
        return new_objects
    else:
        raise ValueError('no agilent file in the input %s'&identifier)
def make_unique_hash(identifier,pipeline,parameters):
    original_file = module_utils.get_inputid(identifier)
    hash_profile={'version': 'agilent',
                   'number of files':str(len(os.listdir(identifier)))}
    hash_result=hash_method.hash_parameters(
        original_file,pipeline,**hash_profile)
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
        outfile,'agilent_files',parameters,objects,single_object)
    return new_objects

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'agilent_files','contents')
    assert os.path.exists(identifier),'the input \
                file %s for extract_agilent_files does not exist'%identifier
    return identifier,single_object
