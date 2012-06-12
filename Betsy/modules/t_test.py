#t_test.py
import hash_method
import gene_ranking
import module_utils
import shutil
import os
import rule_engine
import arrayio
import read_label_file
from genomicode import jmath
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    assert os.path.exists(label_file.identifier),'cannot find label_file %s'%label_file.identifier
    label,label_line,second_line = read_label_file.read(label_file.identifier)
    M = arrayio.read(single_object.identifier)
    assert len(label) == 2, ' the length of label in %s should be 2'%label_file.identifier
    assert len(label[0]) == 2
    assert len(label[1]) == 2
    first = M.slice(None,label[0][0])
    second = M.slice(None,label[1][0])
    t,p = gene_ranking.t_test(first,second)
    higher_group = get_higherexpression(M,label,second_line)
    bonf = jmath.cmh_bonferroni(p)
    fdr = jmath.cmh_fdr_bh(p)
    f=file(outfile,'w')
    header = M._row_order[:]
    header.extend(['p_value','cmh_bonferroni','cmh_fdr','higher_expression'])
    f.write('\t'.join(header))
    f.write('\n')
    for i in range(len(p)):
        for key in M._row_order:
            f.write(M._row_names[key][i])
            f.write('\t')
        f.write(str(p[i])+'\t')
        f.write(str(bonf[i])+'\t')
        f.write(str(fdr[i])+'\t')
        f.write(str(higher_group[i])+'\n')
    f.close()
    assert module_utils.exists_nz(outfile),'the output file %s\
                             for t_test fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 't_test' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),'the input file %s \
                for t_test does not exist'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(parameters,['status'])
    new_object = rule_engine.DataObject(
        'differential_expressed_genes',[parameters['contents']],outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects

def get_higherexpression(M,label,second_line):
    higher_group = []
    assert len(label) == 2, ' the length of label should be 2'
    assert len(label[0]) == 2
    assert len(label[1]) == 2
    first = M.slice(None,label[0][0])
    second = M.slice(None,label[1][0])
    for i in range(M.nrow()):
        group1 = sum(first[i])/float(len(first[i]))
        group2 = sum(second[i])/float(len(second[i]))
        if group1 >= group2:
            higher_group.append(second_line[int(label[0][1])])
        else:
            higher_group.append(second_line[int(label[1][1])])
    return higher_group
        
