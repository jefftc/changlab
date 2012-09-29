# select_num_genes.py
import os
from Betsy import module_utils

def run(parameters,objects,pipeline):
    """select a num of genes according to num_features"""
    import arrayio
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    f = file(outfile,'w')
    num_features = int(parameters['num_features'])
    assert num_features>0, 'the num_features should be >0'
    M = arrayio.read(single_object.identifier)
    M_c = M.matrix(range(0,num_features),None)
    arrayio.tab_delimited_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for select_num_genes fails'%outfile)
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
    filename = 'signal_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for select_num_genes does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
