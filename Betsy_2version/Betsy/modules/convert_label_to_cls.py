#convert_label_to_cls.py

import os
from Betsy import module_utils, read_label_file
import shutil
import arrayio
from Betsy import bie, rulebase

def run(in_nodes, parameters, network):
    data_node,cls_node = in_nodes
    if data_node and cls_node:
        outfile = name_outfile(in_nodes)
        f = file(cls_node.attributes['filename'],'rU')
        text = f.readlines()
        f.close()
        text = [i.rstrip() for i in text]
        label_dict = {}
        for line in text:
            words = line.split('\t')
            if words[1] in label_dict:
                label_dict[words[1]].append(words[0])
            else:
                label_dict[words[1]] = [words[0]]
        class_names = label_dict.keys()
        M = arrayio.read(data_node.attributes['filename'])
        column_names = M.col_names('_SAMPLE_NAME')
        label_line = [0]*len(column_names)
        for i in range(len(class_names)):
            sample_names = label_dict[class_names[i]]
            for sample_name in sample_names:
                index = column_names.index(sample_name)
                label_line[index] = str(i)
        read_label_file.write(outfile,class_names,label_line)
        assert module_utils.exists_nz(outfile),(
            'the output file %s for convert_label_to_cls fails'%outfile)
        new_parameters = parameters.copy()
        new_parameters['filename'] = os.path.split(outfile)[-1]
        out_node = bie.Data(rulebase.ClassLabelFile,**new_parameters)
        return out_node
    return False


def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'class_label_' + original_file + '.cls'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)
