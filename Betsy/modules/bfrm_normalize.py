#bfrm_normalize.py
import os
import hash_method
import subprocess
import module_utils
import rule_engine
import shutil
import arrayio
import tempfile
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile =  get_outfile(parameters,objects,pipeline)
    import Betsy_config
    bfrm_path = Betsy_config.BFRMNORM
    bfrm_BIN = module_utils.which(bfrm_path)
    assert bfrm_BIN,'cannot find the %s' %bfrm_path
    num_factor = 10
    if 'num_factors' in parameters.keys():
        num_factor = int(parameters['num_factors'])
        assert num_factor >= 1, 'the num_factor should be >=1'
        M = arrayio.read(single_object.identifier)
        col_num = M.ncol()
        assert num_factor <= col_num,(
            'the num_factor should be less than %d'%col_num)
    tmp = 'tmp_dir'
    command = ['python', bfrm_BIN,single_object.identifier,'-f',str(num_factor), '-o',tmp]
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert module_utils.exists_nz(tmp),(
        'the output dir %s for bfrm_normalize fails' %tmp)
    assert module_utils.exists_nz(os.path.join(tmp,'normalized.gct')),(
        'the output gct file for bfrm_normalize fails')
    out = os.path.join(tmp,'normalized.gct')
    M = arrayio.read(out)
    M_new = arrayio.convert(M,to_format = arrayio.pcl_format)
    f = file(outfile,'w')
    arrayio.pcl_format.write(M_new,f)
    f.close()
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
    filename = 'signal_bfrm_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for bfrm_normalize does not exist'%single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
    

