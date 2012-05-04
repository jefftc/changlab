#rank_gene_by_sample_ttest.py
import hash_method
import gene_ranking
import module_utils
import shutil
import os
import rule_engine

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file,onj=module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    assert os.path.exists(label_file),'cannot find label_file %s'%label_file
    gene_list = gene_ranking.t_test_for_file(identifier,label_file)
    f = open(outfile,'w')
    f.write('\t'.join(gene_list))
    f.close()
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
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents',pipeline)

    

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier),'the input\
               file %s for rank_gene_by_sample_ttest does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(parameters,['status'])
    new_object = rule_engine.DataObject(
        'gene_list_file',[parameters['contents']],outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
   
