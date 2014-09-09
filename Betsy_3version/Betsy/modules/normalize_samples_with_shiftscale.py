#normalize_sampels_with_shiftscale.py

import os
import shutil
from genomicode import shiftscalenorm
import arrayio
from Betsy import bie3, rulebase
from Betsy import module_utils
from Betsy import read_label_file

def run(in_nodes,parameters, user_input, network):
    data_node, cls_node = in_nodes
    if data_node and cls_node:
        outfile = name_outfile(in_nodes,user_input)
        result,label_line,second_line=read_label_file.read(
            cls_node.identifier)
        assert len(result) == 2, 'for shiftscale,there should be only 2 classes'
        M = arrayio.read(data_node.identifier)
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
        out_node = bie3.Data(rulebase._SignalFile_Merge,**parameters)
        out_object = module_utils.DataObject(out_node,outfile)
        return out_object
    return False

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            datatype='_SignalFile_Merge')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           user_attributes,datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_shiftscale_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    data_node, cls_node = in_nodes
    new_parameters = data_node.data.attributes.copy()
    new_parameters['shiftscale_norm']='yes'
    return new_parameters
    

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

    


