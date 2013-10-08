#normalize_samples_with_combat.py
import os
import subprocess
#from Betsy
import module_utils, read_label_file
from genomicode import config
import bie
import rulebase

def run(in_nodes,parameters):
    data_node, cls_node = in_nodes
    if data_node and cls_node:
        outfile = name_outfile(in_nodes)
        result,label_line,second_line=read_label_file.read(
            cls_node.attributes['filename'])
        assert len(result) >= 2, 'for combat,there should be equal or larger than 2 classes'
        combat_path = config.combatnorm
        combat_BIN = module_utils.which(combat_path)
        assert combat_BIN,'cannot find the %s' %combat_path
        command = ['python', combat_BIN,'-f',data_node.attributes['filename'],'-o',
                   outfile,'-label',cls_node.attributes['filename']]
        process = subprocess.Popen(command,shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        assert module_utils.exists_nz(outfile),(
            'the output file %s for combat fails' %outfile)
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
    filename = 'signal_combat_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

    


