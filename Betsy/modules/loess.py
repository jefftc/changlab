#loess.py

import os
import module_utils
import shutil
import gzip
from genomicode import smarray
import gpr_module
import rule_engine
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    filenames=os.listdir(single_object.identifier)
    keep=[]
    red_sig_matrix=[]
    green_sig_matrix=[]
    red_back_matrix=[]
    green_back_matrix=[]
    sample=[]
    for filename in filenames:
        fileloc=os.path.join(single_object.identifier,filename)
        if not filename.endswith('gpr.gz') and not filename.endswith('gpr'):
            continue
        if filename.endswith('gpr.gz'):   
            samplename=filename[:-7] 
        elif filename.endswith('gpr'):
            samplename=filename.s[:-4]
        sample.append(samplename)
        red_sig,green_sig,red_back,green_back,keep=gpr_module.extract_multiple(fileloc,keep)
        red_sig_matrix.append(red_sig)
        green_sig_matrix.append(green_sig)
        red_back_matrix.append(red_back)
        green_back_matrix.append(green_back)
    red_signal=[[0]*len(red_sig_matrix) for i in range(len(red_sig_matrix[0]))]
    green_signal=[[0]*len(red_sig_matrix) for i in range(len(red_sig_matrix[0]))]
    red_back=[[0]*len(red_sig_matrix) for i in range(len(red_sig_matrix[0]))]
    green_back=[[0]*len(red_sig_matrix) for i in range(len(red_sig_matrix[0]))]
    for i in range(len(red_sig_matrix[0])):
        for j in range(len(red_sig_matrix)):
            red_signal[i][j]=red_sig_matrix[j][i]
            green_signal[i][j]=green_sig_matrix[j][i]
            red_back[i][j]=red_back_matrix[j][i]
            green_back[i][j]=green_back_matrix[j][i]
    x = smarray.correct_background(
       red_signal, green_signal, red_back, green_back,
       method="normexp",offset=50)
    R, G = x
    SIGNAL = smarray.normalize_within_arrays(R, G, method="loess")
    keep[0][1]=keep[0][1].upper() #convert the 'Name' to 'NAME'
    f=open(outfile,'w')
    f.write('\t'.join(keep[0][0:2]))
    f.write('t'.join(sample))
    for i in range(len(keep)-1):
        f.write('\t'.join(keep[i+1][0:2]))
        for j in range(len(SIGNAL[0])):
            f.write('\t')
            f.write(str(SIGNAL[i][j]))
        f.write('\n')
    f.close()
    assert module_utils.exists_nz(outfile),'the output file %s\
                                              for loess fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects


def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_loess_'+original_file+'.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'geo_dataset','contents',['gpr'])
    assert os.path.exists(single_object.identifier),'the input file %s\
                   for loess does not exist'%single_object.identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    parameters = module_utils.renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
