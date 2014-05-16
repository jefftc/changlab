#annotate_probes.py
from genomicode import arrayannot,arrayplatformlib
import arrayio
import os
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils

def run(data_node,parameters,user_input, network):
    outfile = name_outfile(data_node,user_input)
    M = arrayio.read(data_node.identifier)
    all_platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
    if not all_platforms:
        raise ValueError('we cannot guess the platform and annotate the file')
    ids = M._row_order
    probe_header = all_platforms[0][0]
    probe_id = M._row_names[probe_header]
    new_ids = ids[:]
    new_ids.remove(probe_header)
    #annotate_type = parameters['annotate_type']
    #if annotate_type == 'all':
    annotate_header = arrayplatformlib.annotate_header
    #elif annotate_type == 'gene_id':
    #    annotate_header = ['Gene ID']
    dictionary = arrayannot.annotate_probes_multiple(probe_id, annotate_header)
    column=[]
    for id in new_ids:
        flag = True
        for key in dictionary.keys():
            flag = True
            value_list = dictionary[key]
            gene_list = M._row_names[id]
            for gene in gene_list:
                if gene not in value_list:
                    flag = False
                    break
            if flag:
                column.append((key,id))
    header = [i[0] for i in column]
    miss_header = list(set(annotate_header).difference(set(header)))
    original_ids = ids[:]
    for col in miss_header:
        col_2 = col
        if col in original_ids:
            col_1 = col + '_1'
            col_2 = col + '_2'
            M = module_utils.replace_matrix_header(M,col,col_1)
            ids = M._row_order
        ids.append(col_2)
        M._row_order = ids
        M._row_names[col_2] = dictionary[col]
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for annot_probes fails' % outfile)
    out_node = bie3.Data(rulebase.SignalFile_Annotate,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'signal_annot_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters
    

def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
