#agilent.py

import shutil
import os
from genomicode import jmath
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(data_node,parameters, user_input, network,num_cores):
    outfile = name_outfile(data_node,user_input)
    cwd = os.getcwd()
    R = jmath.start_R()
    R('require(limma,quietly=TRUE)')
    R('library(marray)')
    os.chdir(data_node.identifier)
    try:
        R('dir<-getwd()')
        R('files<-list.files(dir)')
        R('x.read<-read.Agilent(files)')
    finally:
        os.chdir(cwd)
    R('xnorm.loc <- maNorm(x.read, norm = "loess")')
    R('x.norm <- maNormScale(xnorm.loc, norm = "p")')
    tmpfile = 'tmp.txt'
    jmath.R_equals(tmpfile,'tmpfile')
    R('write.marray(x.norm,tmpfile)')
    f=open(tmpfile,'r')
    text=f.readlines()
    firstline=text[0].split()
    f.close()
    firstindex = firstline.index('"ProbeName"')
    if '"Sequence"' in firstline:
        secondindex = firstline.index('"Sequence"')
    else:
        secondindex = firstline.index('"ControlType"')
    sample = range(secondindex+1,len(firstline))
    f=open(outfile,'w')
    for i in text:
        line = i.split()
        f.write(line[firstindex]+'\t')
        for j in sample:
            f.write(line[j]+'\t')
        f.write('\n')
    f.close()
    os.remove(tmpfile)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for preprocess_agilent fails' %outfile)
    out_node = bie3.Data(rulebase._SignalFile_Postprocess,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'signal_agilent' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node
