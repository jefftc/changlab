#agilent.py

import shutil
import os
from genomicode import jmath
from Betsy import bie
from Betsy import rulebase
from Betsy import module_utils

def run(data_node,parameters, network):
    outfile = name_outfile(data_node)
    cwd = os.getcwd()
    R = jmath.start_R()
    R('require(limma,quietly=TRUE)')
    R('library(marray)')
    os.chdir(data_node.attributes['filename'])
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
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.SignalFile,**new_parameters)
    return out_node
    


def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_agilent' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
