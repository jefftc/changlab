#fill_missing_with_median.py
import os
import shutil
from Betsy import module_utils,bie3,rulebase
import arrayio
import numpy
 

def run(data_node,parameters, user_input,network):
    outfile = name_outfile(data_node,user_input)
    assert module_utils.is_missing(
        data_node.identifier),'no missing values'
    M = arrayio.read(data_node.identifier)
    f_out = file(outfile, 'w')
    X = M.slice()
    for i in range(M.dim()[0]):
        med = numpy.median([j for j in X[i] if j])
        for j in range(M.dim()[1]):
            if M._X[i][j] is None:
                M._X[i][j] = med
    arrayio.tab_delimited_format.write(M, f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for median_fill_if_missing does not exist'
        % outfile)
    out_node = bie3.Data(rulebase.SignalFile_Impute,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_median_fill_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)


