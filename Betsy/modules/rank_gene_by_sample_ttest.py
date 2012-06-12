#rank_gene_by_sample_ttest.py
import hash_method
import gene_ranking
import module_utils
import shutil
import os
import rule_engine
from genomicode import jmath
import arrayio
import read_label_file
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
    gene_list=[]
    key = M._row_order[0]
    threshold = 0.05
    if 'gene_select_threshold' in parameters.keys():
        threshold = float(parameters['gene_select_threshold'])
    if parameters['gene_order'] == 't_test_p':
        for i in range(len(p)):
            if float(p[i]) < threshold:
                gene_list.append(M._row_names[key][i])
    elif parameters['gene_order'] == 't_test_fdr':
        fdr = jmath.cmh_fdr_bh(p)
        for i in range(len(fdr)):
            if float(fdr[i]) < threshold:
                gene_list.append(M._row_names[key][i])
    f = open(outfile,'w')
    f.write('\t'.join(gene_list))
    f.close()
    assert len(gene_list)>0,'there is no significant genes can be found in ttest'
    assert module_utils.exists_nz(outfile),'the output\
                file %s for rank_gene_by_sample_ttest fails'%outfile
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
    filename = 'gene_list_file_rank_ttest' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

    

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),'the input\
               file %s for rank_gene_by_sample_ttest does not exist'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(parameters,['status'])
    new_object = rule_engine.DataObject(
        'gene_list_file',[parameters['contents']],outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
   
