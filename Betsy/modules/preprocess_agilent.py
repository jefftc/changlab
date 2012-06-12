#agilent.py
import module_utils
import shutil
import os
from genomicode import jmath
import rule_engine
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    
    cwd = os.getcwd()
    R = jmath.start_R()
    R('library(marray)')
    os.chdir(single_object.identifier)
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
    assert module_utils.exists_nz(outfile),'the output\
                 file %s for preprocess_agilent fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
    parameters,single_object,pipeline)
    return new_objects
    

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_agilent' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
            parameters,objects,'agilent_files','contents')
    assert os.path.exists(single_object.identifier),'the input \
        file %s for preprocess_agilent does not exist'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    parameters = module_utils.renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
