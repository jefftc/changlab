#rank_gene_by_sample_ttest.py
import hash_method
import gene_ranking
import module_utils
import shutil
import os
import rule_engine

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    label_file,onj=module_utils.find_object(parameters,objects,'class_label_file','Contents,DatasetId')
    assert os.path.exists(label_file),'cannot find label_file'
    gene_list = gene_ranking.t_test_for_file(identifier,label_file)
    f = open(outfile,'w')
    f.write('\t'.join(gene_list))
    f.close()
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')


def get_outfile(parameters,objects):
    outfile = os.path.join(os.getcwd(),'gene_list_file.txt')
    new_object=rule_engine.DataObject('gene_list_file',[parameters['DatasetId']],outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return outfile,new_objects

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
