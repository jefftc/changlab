#agilent.py
import module_utils
import shutil
import os
from genomicode import jmath

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    
    cwd = os.getcwd()
    R = jmath.start_R()
    R('library(marray)')
    os.chdir(identifier)
    try:
        R('dir<-getwd()')
        R('files<-list.files(dir)')
        R('x.read<-read.Agilent(files)')
    finally:
        os.chdir(cwd)
    R('xnorm.loc <- maNorm(x.read, norm = "loess")')
    R('x.norm <- maNormScale(xnorm.loc, norm = "p")')
    tmpfile = 'tmp.txt'
    jmath.R_equals('\"' + tmpfile + '\"','tmpfile')
    R('write.marray(x.norm,tmpfile)')
    f=open(tmpfile,'r')
    text=f.readlines()
    firstline=text[0].split()
    f.close()
    firstindex = firstline.index('"ProbeName"')
    secondindex = firstline.index('"Sequence"')
    sample = range(secondindex+1,len(firstline))
    f=open(outfile,'w')
    for i in text:
        line = i.split()
        f.write(line[firstindex]+'\t')
        for j in sample:
            f.write(line[j]+'\t')
        f.write('\n')
    f.close()
    os.remove(tmpfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
    parameters,single_object,pipeline)
    return new_objects
    

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'geo_dataset','Contents,DatasetId',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(parameters,objects,'geo_dataset','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
