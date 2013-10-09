#fill_missing_with_median.py
import os
import shutil
#from Betsy
import module_utils
import arrayio
import numpy
import bie
import rulebase

def run(data_node,parameters):
    outfile = name_outfile(data_node)
    assert module_utils.is_missing(data_node.attributes['filename']),'no missing values'
    M = arrayio.read(data_node.attributes['filename'])
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
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_median_fill_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters


def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


