#median_fill_if_missing.py
import os
import shutil
import module_utils
import arrayio
def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    if not is_missing(identifier):
            shutil.copyfile(identifier,outfile)
    else:
        M = arrayio.read(identifier)
        f_out = file(outfile,'w')
        X = M.slice()
        for i in range(M.dim()[0]):
           med = X[i]
           for j in range(M.dim()[1]):
               if M._X[i][j] == None:
                        M._X[i][j] = str(med)
        arrayio.pcl_format.write(M,f_out)    
        f_out.close()
    assert module_utils.exists_nz(outfile),'the output\
        file %s for median_fill_if_missing does not exist'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents,preprocess',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(identifier),'the input \
        file %s for median_fill_if_missing does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
    
def is_missing(identifier):
    import arrayio
    M = arrayio.read(identifier)
    has_missing = False
    for i in range(M.dim()[0]):
       for j in range(M.dim()[1]):
           if M._X[i][j] == None:
               has_missing = True
               break
       if has_missing:
            break
    return has_missing

def get_median(input_list):
    new_list = sorted(input_list)
    if len(new_list) % 2 == 1:
        return new_list[(len(new_list)+1)/2-1]
    else:
        lower = new_list[len(new_list)/2-1]
        upper = new_list[len(new_list)/2]
        return (float(lower + upper)) / 2  


