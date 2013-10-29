#filter_genes_by_missing_values.py
import os
from Betsy import module_utils
from Betsy import bie, rulebase

def run(data_node,parameters, network):
    outfile = name_outfile(data_node)
    import arrayio
    f_out = file(outfile,'w')
    M = arrayio.read(data_node.attributes['filename'])
    I_good = []
    #get the percentage of gene filter
    percent=float(parameters[str('filter')])/100
    for i in range(M.dim()[0]):
       missing_count = 0
       for j in range(M.dim()[1]):
           if M._X[i][j] in [None,'NA']:
                missing_count = missing_count + 1
       if float(missing_count)/M.dim()[1]<percent:
            I_good.append(i)
    M_c = M.matrix(I_good,None)
    arrayio.tab_delimited_format.write(M_c,f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for gene_filter fails'%outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_filter_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

