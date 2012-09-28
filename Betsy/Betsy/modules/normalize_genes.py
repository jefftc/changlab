#normalize_genes.py
import os
from Betsy import module_utils
import subprocess
def run(parameters,objects,pipeline):
    """variance or sum_of_square"""
    norm_para = ["variance","sum_of_squares"]
    if parameters['gene_normalize'] not in norm_para:
        raise ValueError("Cannot recognizd the normalize parameter")
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    if parameters['gene_normalize'] == "variance":
        import arrayio
        from genomicode import jmath
        f = file(outfile,'w')
        M = arrayio.read(single_object.identifier,format = arrayio.pcl_format)
        M_n = jmath.safe_norm_mv(M.slice())
        M._X = M_n
        M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
        arrayio.pcl_format.write(M_c,f)
        f.close()
    elif parameters['gene_normalize'] == "sum_of_squares":
        CLUSTER_BIN = 'cluster'
        process = subprocess.Popen([CLUSTER_BIN,'-f',single_object.identifier,'-ng',
                                    '-u',outfile],shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        
        outputfile = outfile+'.nrm'
        os.rename(outputfile,outfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for normalize fails'%outfile)
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
    filename = 'signal_normalize_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for normalize does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
