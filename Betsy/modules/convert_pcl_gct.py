#convert_pcl_gct.py
import os
import module_utils

def run(parameters,objects,pipeline):
    """convert pcl signal file to gct format"""
    import arrayio
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    f = file(outfile,'w')
    M = arrayio.read(identifier)
    M_c = arrayio.convert(M,to_format=arrayio.gct_format)
    arrayio.gct_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),'the output file %s\
                             for convert_pcl_gct fails'%outfile
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
    assert os.path.exists(identifier),'the input file %s\
                   for convert_pcl_gct does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
