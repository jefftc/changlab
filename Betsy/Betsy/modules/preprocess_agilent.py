#agilent.py

import os
from genomicode import jmath
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    cwd = os.getcwd()
    R = jmath.start_R()
    R('require(limma,quietly=TRUE)')
    R('library(marray)')
    os.chdir(in_data.identifier)
    try:
        R('dir<-getwd()')
        R('files<-list.files(dir)')
        R('x.read<-read.Agilent(files)')
    finally:
        os.chdir(cwd)
    
    R('xnorm.loc <- maNorm(x.read, norm = "loess")')
    R('x.norm <- maNormScale(xnorm.loc, norm = "p")')
    tmpfile = 'tmp.txt'
    jmath.R_equals(tmpfile, 'tmpfile')
    R('write.marray(x.norm,tmpfile)')
    f = open(tmpfile, 'r')
    text = f.readlines()
    firstline = text[0].split()
    f.close()
    firstindex = firstline.index('"ProbeName"')
    if '"Sequence"' in firstline:
        secondindex = firstline.index('"Sequence"')
    else:
        secondindex = firstline.index('"ControlType"')
    
    sample = range(secondindex + 1, len(firstline))
    f = open(outfile, 'w')
    for i in text:
        line = i.split()
        f.write(line[firstindex] + '\t')
        for j in sample:
            f.write(line[j] + '\t')
        f.write('\n')
    
    f.close()
    os.remove(tmpfile)
    assert module_utils.exists_nz(outfile), (
        'the output file %s for preprocess_agilent fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Postprocess, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_agilent' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
