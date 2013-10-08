#convert_signal_to_gct.py
import os
#from Betsy
import module_utils
import bie
import rulebase

def run(data_node,parameters):
    """convert signal file to gct format"""
    import arrayio
    outfile = name_outfile(data_node)
    f = file(outfile, 'w')
    M = arrayio.read(data_node.attributes['filename'])
    M_c = arrayio.convert(M, to_format=arrayio.gct_format)
    arrayio.gct_format.write(M_c, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_pcl_gct fails' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile2,**new_parameters)
    return out_node


def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_' + original_file + '.gct'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):

    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

