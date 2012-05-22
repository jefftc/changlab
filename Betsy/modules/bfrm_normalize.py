#bfrm_normalize.py
import os
import hash_method
import subprocess
import module_utils
import rule_engine
import shutil
def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile =  get_outfile(parameters,objects,pipeline)
    import Betsy_config
    bfrm_path = Betsy_config.BFRMNORM
    bfrm_BIN = module_utils.which(bfrm_path)
    assert bfrm_BIN,'cannot find the %s' %bfrm_path
    num_factor = 10
    if 'num_factors' in parameters.keys():
        num_factor = int(parameters['num_factors'])
        assert num_factor >= 1, 'the num_factor should be >=1'
        import arrayio
        M = arrayio.read(identifier)
        col_num = M.ncol()
        assert num_factor <= col_num,'the num_factor should be less than %d'%col_num
    tmp = 'tmp_dir'
    command = ['python', bfrm_BIN,identifier,'-f',str(num_factor), '-o',tmp]
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(tmp),'the output dir %s\
             for bfrm_normalize fails' %tmp
    assert module_utils.exists_nz(os.path.join(tmp,'normalized.gct')),'the output gct file \
             for bfrm_normalize fails'
    shutil.copyfile(os.path.join(tmp,'normalized.gct'),outfile)
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
        %s for bfrm_normalize does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
    

