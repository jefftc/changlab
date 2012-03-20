#normalize.py
import os
import module_utils
import subprocess
def run(parameters,objects,pipeline):
    """variance or sum_of_square"""
    norm_para = ["variance","sum_of_squares"]
    if parameters['gene_normalize'] not in norm_para:
        raise ValueError("Cannot recognizd the normalize parameter")
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    if parameters['gene_normalize'] == "variance":
        import arrayio
        from genomicode import jmath
        f = file(outfile,'w')
        M = arrayio.read(identifier,format = arrayio.pcl_format)
        M_n = jmath.safe_norm_mv(M.slice())
        M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
        arrayio.pcl_format.write(M_c,f)
        f.close()
    elif parameters['gene_normalize'] == "sum_of_squares":
        CLUSTER_BIN = 'cluster'
        process = subprocess.Popen([CLUSTER_BIN,'-f',identifier,'-ng',
                                    '-u',outfile],shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        
        outputfile = outfile+'.nrm'
        os.rename(outputfile,outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId',
        pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
