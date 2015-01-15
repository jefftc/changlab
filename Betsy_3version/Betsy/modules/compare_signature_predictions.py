#compare_signature_predictions.py
import os
import arrayio
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(in_nodes,parameters, user_input, network,num_cores):
    """compare signature predictions"""
    outfile =  name_outfile(in_nodes,user_input)
    data_node1, data_node2, data_node3 = in_nodes
    assert 'probabilities.pcl' in os.listdir(data_node1.identifier)
    assert 'probabilities.pcl' in os.listdir(data_node2.identifier)
    assert 'probabilities.pcl' in os.listdir(data_node3.identifier)
    f1 = os.path.join(data_node1.identifier,'probabilities.pcl')
    f2 = os.path.join(data_node2.identifier,'probabilities.pcl')
    f3 = os.path.join(data_node3.identifier,'probabilities.pcl')
    M1=arrayio.read(f1)
    M2=arrayio.read(f2)
    M3=arrayio.read(f3)
    total_sample = list(set(M1.col_names('_SAMPLE_NAME')+
             M2.col_names('_SAMPLE_NAME')+
             M3.col_names('_SAMPLE_NAME')))
    d1=get_new_matrix(total_sample,f1)
    d2=get_new_matrix(total_sample,f2)
    d3=get_new_matrix(total_sample,f3)
    f = file(outfile,'w')
    header = ['SigID','NAME']+total_sample
    f.write('\t'.join(header)+'\n')
    for line in d1:
        f.write('\t'.join(line)+'\n')
    for line in d2:
        f.write('\t'.join(line)+'\n')
    for line in d3:
        f.write('\t'.join(line)+'\n')
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for compare_signature_predictions fails' %outfile)
    out_node = bie3.Data(rulebase.ScoreCompareReportFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def name_outfile(in_nodes,user_input):
    data_node1, data_node2, data_node3 = in_nodes
    original_file = module_utils.get_inputid(data_node1.identifier)
    filename = 'compare_three_signature_predictions_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node1, data_node2, data_node3 = in_nodes
    identifier = data_node1.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node1 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            datatype='SignatureScore',
                                            optional_key='preprocess',optional_value='affymetrix')
    data_node2 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            datatype='SignatureScore',
                                            optional_key='preprocess',optional_value='agilent')
    data_node3 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            datatype='SignatureScore',
                                            optional_key='preprocess',optional_value='RSEM_genes')
    
    return data_node1, data_node2, data_node3
    
def get_new_matrix(total_sample,filename):
    M=arrayio.read(filename)
    samples = M.col_names('_SAMPLE_NAME')
    data = []
    indexes = []
    for i in total_sample:
        index = samples.index(i) if i in samples else None
        indexes.append(index)
    for i in range(M.nrow()):
        line = [M.row_names('SigID')[i],M.row_names('NAME')[i]]
        for j in indexes:
           if j is None:
               line.append('NA')
           else:
               line.append(str(M._X[i][j]))
        data.append(line)
    return data
