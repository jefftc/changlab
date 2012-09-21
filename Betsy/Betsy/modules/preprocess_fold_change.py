#preprocess_fold_change.py
import os
import module_utils

def run(parameters,objects,pipeline):
    """run preprocessdataset """
    import arrayio
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    threshold = 20
    ceiling = 16000
    min_fold_change = 5
    min_delta = 100.0
    M = arrayio.read(single_object.identifier)
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
    arrayio.tab_delimited_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for preprocessdataset fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_preprocessdataset_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for preprocessdataset does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
