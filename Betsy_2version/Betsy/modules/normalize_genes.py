#normalize_genes.py
import os
#from Betsy
import module_utils
import subprocess
import bie
import rulebase

def run(data_node,parameters):
    """variance or sum_of_square"""
    norm_para = ["variance","sum_of_squares"]
    if parameters['gene_normalize'] not in norm_para:
        raise ValueError("Cannot recognizd the normalize parameter")
    outfile = name_outfile(data_node)
    if parameters['gene_normalize'] == "variance":
        import arrayio
        from genomicode import jmath
        f = file(outfile,'w')
        M = arrayio.read(data_node.attributes['filename'],
                         format = arrayio.pcl_format)
        M_n = jmath.safe_norm_mv(M.slice())
        M._X = M_n
        M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
        arrayio.pcl_format.write(M_c,f)
        f.close()
    elif parameters['gene_normalize'] == "sum_of_squares":
        CLUSTER_BIN = 'cluster'
        process = subprocess.Popen([CLUSTER_BIN,'-f',data_node.attributes['filename'],'-ng',
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
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile1,**new_parameters)
    return out_node


def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'signal_normalize_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    new_parameters = parameters.copy()
    new_parameters['format']='tdf'
    return new_parameters
    

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
