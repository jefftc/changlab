#normalize_samples_with_dwd.py

import os
import shutil
from genomicode import dwdnorm
import arrayio
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils, read_label_file

def run(in_nodes, parameters, user_input,network):
    data_node,cls_node = in_nodes
    if data_node and cls_node:
        outfile = name_outfile(in_nodes,user_input)
        M = arrayio.read(data_node.identifier)
        result,label_line,second_line=read_label_file.read(
            cls_node.identifier)
        assert len(result) == 2, 'for dwd,there should be only 2 classes'
        assert [i in ['0','1'] for i in label_line] == [True]*len(label_line),(
            'the label of class shoul be 0 and 1')
        y = [i.replace('0','-1') for i in label_line]
        M_y = dwdnorm.normalize(M,y)
        f = file(outfile,'w')
        arrayio.tab_delimited_format.write(M_y,f)
        f.close()
        assert module_utils.exists_nz(outfile),(
            'the output file %s for dwd fails'%outfile)
        new_parameters = parameters.copy()
        out_node = bie3.Data(rulebase.SignalFile_Merge,**parameters)
        out_object = module_utils.DataObject(out_node,outfile)
        return out_object
    return False


def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            datatype='SignalFile_Merge')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           user_attributes,datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_dwd_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    data_node, cls_node = in_nodes
    new_parameters = data_node.data.attributes.copy()
    new_parameters['dwd_norm']='yes'
    return new_parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)
