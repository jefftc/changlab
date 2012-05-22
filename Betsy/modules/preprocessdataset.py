#preprocessdataset.py
import os
import module_utils

def run(parameters,objects,pipeline):
    """run preprocessdataset """
    import arrayio
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    threshold = 20
    ceiling = 16000
    min_fold_change = 5
    min_delta = 100.0
    M = arrayio.read(identifier)
    X = M.slice()
    I_good = []
    for i in range(M.nrow()):
        for j in range(len(X[i])):
            if X[i][j]<threshold:
                M._X[i][j] = threshold
            if X[i][j]>ceiling:
                M._X[i][j]=ceiling
        gene = M._X[i]
        fold_change = max(gene)/float(min(gene))
        delta = max(gene)-min(gene)
        if fold_change >= min_fold_change and delta>=min_delta:
            I_good.append(i)
    f = file(outfile,'w')
    M_c = M.matrix(I_good,None)
    arrayio.pcl_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),'the output\
                        file %s for preprocessdataset fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents',
        pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier),'the input \
                file %s for preprocessdataset does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects