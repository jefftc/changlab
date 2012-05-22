#convert_to_pcl_if_not_pcl.py
import os
import module_utils
import shutil

def run(parameters,objects,pipeline):
    """convert not xls signal file to pcl format"""
    import arrayio
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    f = file(outfile,'w')
    M = arrayio.choose_format(identifier)
    if M.__name__[8:-7] == 'pcl':
        shutil.copyfile(identifier,outfile)
        f.close()
    else:
        M = arrayio.read(identifier)
        M_c = arrayio.convert(M,to_format = arrayio.pcl_format)
        arrayio.pcl_format.write(M_c,f)
        f.close()
    assert module_utils.exists_nz(outfile),'the output\
    file %s for convert_to_pcl_if_not_pcl fails'%outfile
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
    assert os.path.exists(identifier),'the input file \
        %s for convert_to_pcl_if_not_pcl does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects