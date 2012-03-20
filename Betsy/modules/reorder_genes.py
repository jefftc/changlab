#reorder_genes.py
import hash_method
import gene_ranking
import module_utils
import os
def run(parameters,objects,pipeline):
    #also if input is other kind of file
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    #read the gene order list
    gene_list_file,obj=module_utils.find_object(parameters,
                                objects,'gene_list_file','DatasetId')
    assert os.path.exists(gene_list_file),'cannot find gene_list_file'    
    gene_list = open(gene_list_file,'r').read().split()
    #read the pcl signal file
    f_signal= open(identifier,'r')
    content = f_signal.readlines()
    f_signal.close()
    #get the original gene list
    original_list = []
    for i in range(1,len(content)):
        original_list.append(content[i].split()[0])
    #get the order index and write to the outout file
    indexlist = gene_ranking.find_sorted_index(original_list,gene_list)
    f = open(outfile,'w')
    f.write(content[0]+'\n')
    for i in range(len(indexlist)):
        f.write(content[indexlist[i]+1]+'\n')
    f.close()
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return  module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
    assert os.path.exists(identifier),'the input file does not exist'
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
