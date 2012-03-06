#preprocess_gpr.py
import os
import module_utils
import hash_method
import rule_engine
import gzip
import gpr_module

def run(parameters,objects,pipeline):
    """preprocess the input gpr files, generate a signal file"""
    #preprocess the cel file to text signal file
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
    filenames=os.listdir(identifier)
    keep=[]
    logmatrix=[]
    for filename in filenames:
        fileloc=os.path.join(identifier,filename)
        if not filename.endswith('gpr.gz') and not filename.endswith('gpr'):
            continue
        logratio,keep=gpr_module.extract_gpr(fileloc,keep)
        logmatrix.append(logratio)
    keep[0][1]=keep[0][1].upper() #convert the 'Name' to 'NAME'
    f=open(outfile,'w')
    for i in range(len(keep)):
        f.write('\t'.join(keep[i][0:2]))
        for j in range(len(logmatrix)):
            f.write('\t')
            f.write(logmatrix[j][i])
        f.write('\n')
    f.close()
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'geo_dataset','Contents,DatasetId',pipeline)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'geo_dataset','Contents,DatasetId','signal_file',pipeline)

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'geo_dataset','Contents,DatasetId')
   
