#normalize_with_loess.py
import os
from Betsy import module_utils
import shutil
import gzip
from genomicode import smarray
from Betsy import gpr_module, bie3,rulebase

def run(data_node, parameters, user_input,network):
    parameters = get_out_attributes(parameters,data_node)
    filenames=os.listdir(data_node.identifier)
    outfile = name_outfile(data_node,user_input)
    keep=[]
    red_sig_matrix=[]
    green_sig_matrix=[]
    red_back_matrix=[]
    green_back_matrix=[]
    sample=[]
    for filename in filenames:
        fileloc=os.path.join(data_node.identifier,filename)
        if not filename.endswith('gpr.gz') and not filename.endswith('gpr'):
            continue
        if filename.endswith('gpr.gz'):   
            samplename=filename[:-7] 
        elif filename.endswith('gpr'):
            samplename=filename.s[:-4]
        sample.append(samplename)
        red_sig,green_sig,red_back,green_back,keep=gpr_module.extract_multiple(fileloc,keep)
        red_sig_matrix.append(red_sig)
        green_sig_matrix.append(green_sig)
        red_back_matrix.append(red_back)
        green_back_matrix.append(green_back)
    red_signal=[[0]*len(red_sig_matrix) for i in range(len(red_sig_matrix[0]))]
    green_signal=[[0]*len(red_sig_matrix) for i in range(len(red_sig_matrix[0]))]
    red_back=[[0]*len(red_sig_matrix) for i in range(len(red_sig_matrix[0]))]
    green_back=[[0]*len(red_sig_matrix) for i in range(len(red_sig_matrix[0]))]
    for i in range(len(red_sig_matrix[0])):
        for j in range(len(red_sig_matrix)):
            red_signal[i][j]=red_sig_matrix[j][i]
            green_signal[i][j]=green_sig_matrix[j][i]
            red_back[i][j]=red_back_matrix[j][i]
            green_back[i][j]=green_back_matrix[j][i]
    x = smarray.correct_background(
       red_signal, green_signal, red_back, green_back,
       method="normexp",offset=50)
    R, G = x
    SIGNAL = smarray.normalize_within_arrays(R, G, method="loess")
    keep[0][1]=keep[0][1].upper() #convert the 'Name' to 'NAME'
    f=open(outfile,'w')
    f.write('\t'.join(keep[0][0:2]))
    f.write('\t')
    f.write('\t'.join(sample))
    f.write('\n')
    for i in range(len(keep)-1):
        f.write('\t'.join(keep[i+1][0:2]))
        for j in range(len(SIGNAL[0])):
            f.write('\t')
            f.write(str(SIGNAL[i][j]))
        f.write('\n')
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for loess fails'%outfile)
    out_node = bie3.Data(rulebase.SignalFile_Postprocess,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object

def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_loess_'+original_file+'.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,
                                         parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node
