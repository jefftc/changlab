#normalize_sampels_with_shiftscale.py

import os
#from Betsy
import module_utils
import shutil
#from Betsy
import read_label_file
from genomicode import shiftscalenorm
import arrayio
import bie
import rulebase

def run(in_nodes,parameters):
    data_node, cls_node = in_nodes
    if data_node and cls_node:
        outfile = name_outfile(in_nodes)
        result,label_line,second_line=read_label_file.read(
            cls_node.attributes['filename'])
        assert len(result) == 2, 'for shiftscale,there should be only 2 classes'
        M = arrayio.read(data_node.attributes['filename'])
        index1=result[0][0]
        index2=result[1][0]
        M_1=M.matrix(None,index1)
        M_2=M.matrix(None,index2)
        M_y = shiftscalenorm.normalize(M_1,M_2)
        for i in range(M_y.dim()[0]):
            for j in range(M_y.dim()[1]):
                if str(M_y._X[i][j]) == 'nan':         
                    M_y._X[i][j] = M_2._X[i][0]
        for j in range(M.nrow()):
            for i in range(len(index1)):
                M._X[j][index1[i]]=M_y._X[j][i]
            
        f = file(outfile,'w')
        arrayio.tab_delimited_format.write(M,f)
        f.close()
        assert module_utils.exists_nz(outfile),(
            'the output file %s for shiftscale fails'%outfile)
        new_parameters = parameters.copy()
        new_parameters['filename'] = os.path.split(outfile)[-1]
        out_node = bie.Data(SignalFile_rule.SignalFile,**new_parameters)
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
    filename = 'signal_shiftscale_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

    


