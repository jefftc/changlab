#reorder_genes.py
import hash_method
import gene_ranking
import module_utils
import os
def run(parameters,objects,pipeline):
    #also if input is other kind of file
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
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
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)

def get_outfile(parameters,objects,pipeline):
    return  module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId','signal_file',pipeline)

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
