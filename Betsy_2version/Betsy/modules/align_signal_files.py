#select_common_genes.py
import os
import sys
import arrayio
#from Betsy
import module_utils
from genomicode import matrixlib
import bie
import rulebase

def run(in_nodes,parameters):
    train_node,test_node = in_nodes
    outfile = name_outfile(in_nodes)
    file1, file2 = module_utils.convert_to_same_platform(
        train_node.attributes['filename'],
        test_node.attributes['filename'])
    training = arrayio.read(file1)
    test = arrayio.read(file2)
    [M_A, M_B] = matrixlib.align_rows(training, test)
    assert M_A.nrow() > 0, (
        'there is no common genes betwee %s and %s'
        % (train_node.attributes['filename'], test_node.attributes['filename']))
##    f=file(outfile,'w')
##    arrayio.tab_delimited_format.write(M_A,f)
##    f.close()
##    test_filename = os.path.splitext(outfile)[0]+'_test.txt'
##    f = file(test_filename,'w')
##    arrayio.tab_delimited_format.write(M_B, f)
##    f.close()
    f = file(outfile, 'w')
    if parameters['contents'] == 'class0,class1':
        arrayio.tab_delimited_format.write(M_A, f)
    elif parameters['contents'] == 'test':
        arrayio.tab_delimited_format.write(M_B, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for align_signal_file fails' % outfile)
##    assert module_utils.exists_nz(test_filename), (
##        'the output file %s for align_signal_file fails' % test_filename)
    new_parameters = parameters.copy()
##    new_parameters['contents'] = 'class0,class1'
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile2,**new_parameters)
##    new_parameters['contents'] = 'test'
##    new_parameters['filename'] = os.path.split(test_filename)[-1]
##    out_node_test = bie.Data(SignalFile2_rule.SignalFile2,**new_parameters)
    return out_node

def name_outfile(in_nodes):
    train_data_node,test_data_node = in_nodes
    original_file = module_utils.get_inputid(train_data_node.attributes['filename'])
    filename = 'signal_common_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    train_data_node,test_data_node = in_nodes
    identifier = train_data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes):
    train_data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,contents='class0,class1',
                                            datatype='SignalFile2')
    test_data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,contents='test',
                                            datatype='SignalFile2')
    return train_data_node,test_data_node



