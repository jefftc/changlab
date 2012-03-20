#convert_v3_4.py

import os
import module_utils
import rule_engine
import shutil
import gzip

def run(parameters,objects,pipeline):
    """convert from cc to v3_4"""
    from genomicode import affyio
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    filenames = os.listdir(identifier)
    os.mkdir(outfile)
    for filename in filenames:
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
            assert os.path.exists(cel_file),'the unzip fails'
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
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'geo_dataset','Contents,DatasetId',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object =  module_utils.find_object(
        parameters,objects,'geo_dataset','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'geo_dataset',parameters,objects,single_object)
    return new_objects
