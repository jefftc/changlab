#normalize_samples_with_dwd.py

import os
import module_utils, read_label_file
import shutil
from genomicode import dwdnorm
import arrayio
import bie
import rulebase

def run(in_nodes, parameters):
    data_node,cls_node = in_nodes
    if data_node and cls_node:
        outfile = name_outfile(in_nodes)
        M = arrayio.read(data_node.attributes['filename'])
        result,label_line,second_line=read_label_file.read(
            cls_node.attributes['filename'])
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
        new_parameters['filename'] = os.path.split(outfile)[-1]
        out_node = bie.Data(rulebase.SignalFile,**new_parameters)
        return out_node
    return False


def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_dwd_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
