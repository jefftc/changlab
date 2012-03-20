#filter_by_fold_change.py
import os
import module_utils
import math
def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    import arrayio
    min_fold_change = 1.5
    f_out = file(outfile,'w')
    M = arrayio.read(identifier)
    I_good = []
    for i in range(M.dim()[0]):
        med = median(M.slice(i,None)[0])
        new = [x/med for x in M.slice(i,None)[0]]
        flag = [j > min_fold_change for j in new]
        if True in flag:
            I_good.append(i)
    M_c = M.matrix(I_good,None)
    arrayio.pcl_format.write(M_c,f_out)
    f_out.close()
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects

def median(numberList):
    numberList.sort()
    middle = len(numberList)/2
    if len(numberList)%2 == 1:
        return numberList[middle]
    else:
        right = int(math.ceil(middle))
        left = right - 1
        return (numberList[right]+numberList[left])/2.0

