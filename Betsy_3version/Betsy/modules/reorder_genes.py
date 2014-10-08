#reorder_genes.py

import os
import arrayio
from genomicode import arrayplatformlib, config
from Betsy import bie3
from Betsy import rulebase
from Betsy import gene_ranking
from Betsy import module_utils

def run(in_nodes, parameters, user_input,network):
    data_node,gene_node = in_nodes
    outfile = name_outfile(in_nodes,user_input)
    #read the gene order list
    gene_list = open(gene_node.identifier, 'r').read().split()
    M = arrayio.read(data_node.identifier)
    x = arrayplatformlib.identify_all_platforms_of_matrix(M)
    if x:
        id = x[0][0]
        platform = x[0][1]
        chip = arrayplatformlib.identify_platform_of_annotations(gene_list)
        if not chip:
            chip = []
        signal_file = data_node.identifier
        if platform == chip:
            tmpfile = data_node.identifier
        else:
            platform_name = 'unknown_platform'
            if 'platform_name' in user_input:
                platform_name = user_input['platform_name']
            if platform_name in chip:#, 'unknown_platform':
                import subprocess
                Annot_path = config.annotate_matrix
                Annot_BIN = module_utils.which(Annot_path)
                assert Annot_BIN, 'cannot find the %s' % Annot_path
                signal_file = 'tmp'    
                command = ['python', Annot_BIN, '-f', single_object.identifier,
                           '-o', signal_file, "--platform", chip]
                process = subprocess.Popen(command, shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
                error_message = process.communicate()[1]
                if error_message:
                    raise ValueError(error_message)
                assert module_utils.exists_nz(signal_file), 'the platform conversion fails'
                id = parameters['platform']
                M = arrayio.read(signal_file)
            elif platform_name == platform:
                infile = gene_list_file.identifier
                f=file(infile, 'rU')
                genes=f.readlines()
                f.close()
                gene_list = module_utils.convert_gene_list_platform(
                    genes, platform)
    else:
        id = M._row_order[0]
    #id = M._row_order[3]
    #print id
    original_list = M._row_names[id]
    #get the order index and write to the outout file
    indexlist = gene_ranking.find_sorted_index(original_list, gene_list)
    M_new = M.matrix(indexlist, None)
    f = open(outfile, 'w')
    arrayio.tab_delimited_format.write(M_new, f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for reorder_genes fails' % outfile)
    out_node = bie3.Data(rulebase._SignalFile_Order,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                            datatype='_SignalFile_Order')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           user_attributes,datatype='GeneListFile')
    return data_node, cls_node

def name_outfile(in_nodes,user_input):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'gene_reorder' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node,cls_node = in_nodes
    identifier = data_node.identifier
    newparameters = parameters.copy()
    newparameters['genelistfile']=cls_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,newparameters,user_input)


   




