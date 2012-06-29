#filter_by_fold_change.py
import os
import module_utils
import math
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    folder_number = parameters['filter_fc']
    try:
        fold_change = float(folder_number)
    except ValueError:
        raise 'The filter_fc should be a number'
    assert fold_change > 0, 'the filter_fc should be positive'
    import arrayio
    min_fold_change = math.log(fold_change,2)
    f_out = file(outfile,'w')
    M = arrayio.read(single_object.identifier)
    I_good = []
    X = M.slice()
    for i in range(M.nrow()):
        gene = X[i]
        fold_change = max(gene)-min(gene)
        if fold_change > min_fold_change:
            I_good.append(i)
    M_c = M.matrix(I_good,None)
    arrayio.pcl_format.write(M_c,f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile),'the output\
                              file for filter_by_fold_change fails'
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_foldchange_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents',pipeline)
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),'the input\
                file for filter_by_fold_change does not exist'
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects

