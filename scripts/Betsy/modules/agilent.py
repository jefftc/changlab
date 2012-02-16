#agilent.py
import module_utils
import shutil
import os
import Betsy_config
from genomicode import jmath

def run(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
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
    module_utils.write_Betsy_parameters_file(
    parameters,single_object)
    return new_objects
    

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'geo_dataset','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'geo_dataset','Contents,DatasetId','signal_file')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'geo_dataset','Contents,DatasetId')
    
