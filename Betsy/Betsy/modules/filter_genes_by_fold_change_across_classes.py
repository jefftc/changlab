#filter_genes_by_fold_change_across_classes.py
import os
from Betsy import module_utils, read_label_file
from genomicode import jmath
import math
from time import strftime,localtime

def run(parameters,objects,pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    import arrayio
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_object = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    # obtain the class label
    label, label_line, second_line = read_label_file.read(
        label_object.identifier)
    class_num = len(label)
    assert class_num == 2, 'the number of class is not 2'
    fc = int(parameters['group_fc'])
    M = arrayio.read(single_object.identifier)
    first = M.slice(None, label[0][0])
    second = M.slice(None, label[1][0])
    X = M.slice()
    I_good = []
    for i in range(M.nrow()):
        fold_change = abs(jmath.mean(first[i])-jmath.mean(second[i]))
        if fold_change >= math.log(fc,2):
            I_good.append(i)
    assert I_good, 'there is no gene is significant in fold change with 2'
    f = file(outfile,'w')
    M_c = M.matrix(I_good,None)
    arrayio.tab_delimited_format.write(M_c,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for filter_genes_by_fold_change_across_classes fails'
         % outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,[single_object,label_object],pipeline,outfile,
        starttime,user,jobname)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_group_fc' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for filter_genes_by_fold_change_across_classes does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
