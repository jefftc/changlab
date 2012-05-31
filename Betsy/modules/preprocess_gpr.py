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
    outfile = get_outfile(parameters,objects,pipeline) 
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
    assert module_utils.exists_nz(outfile),'the output\
                        file %s for preprocess_gpr fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
                            parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'gpr_files','contents',pipeline)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
                     parameters,objects,'gpr_files','contents')
    assert os.path.exists(identifier),'the input \
                file %s for preprocess_gpr does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    parameters = renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('signal_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
