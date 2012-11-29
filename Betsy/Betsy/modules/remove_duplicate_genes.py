#remove_duplicate_genes.py
import os
from Betsy import module_utils
from genomicode import jmath,arrayplatformlib,arrayannot
import arrayio
import re


def run(parameters,objects,pipeline):
    """remove duplicate genes"""
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    if parameters['unique_genes'] == 'average_genes':
        M_new = average_genes(M)
    elif parameters['unique_genes'] == 'high_var':
        M_new = get_high_variance(M)
    elif parameters['unique_genes'] == 'first_gene':
        M_new = pick_first_one(M)
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_new,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for remove_duplicate_genes fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects


def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)


def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile


def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for remove_duplicate_genes does not exist'
        %single_object.identifier)
    return single_object


def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects


def get_duplicated_genes(M):
    ProbeId,GeneID = guess_gene_header(M)
    gene_names = M._row_names[GeneID]
    name2indexes=dict()
    for i in range(M.nrow()):
        b = re.match("^[-]*",gene_names[i])
        if len(gene_names[i])==0:
            continue
        elif len(b.group())>0:
            continue
        else:
            if gene_names[i] not in name2indexes:
                name2indexes[gene_names[i]]=[]
            name2indexes[gene_names[i]].append(i)
    return name2indexes


def pick_first_one(M):
    name2indexes = get_duplicated_genes(M)
    ProbeId,GeneID = guess_gene_header(M)
    probe_names = M._row_names[ProbeId]
    for key in name2indexes.keys():
        if len(name2indexes[key])>1:
            a = [(probe_names[i],i) for i in name2indexes[key]]
            a.sort()
            index = a[0][1]
            name2indexes[key]=[index]
    all_index = name2indexes.values()
    all_index.sort()
    all_index = [i[0] for i in all_index if len(i)==1]
    M_new = M.matrix(all_index,None)
    return M_new


def get_high_variance(M):
    name2indexes = get_duplicated_genes(M)
    for key in name2indexes.keys():
        if len(name2indexes[key])>1:
            a=[(jmath.var(M._X[i]),i) for i in name2indexes[key]]
            a.sort()
            index = a[-1][1]
            name2indexes[key]=[index]
    all_index = name2indexes.values()
    all_index.sort()
    all_index = [i[0] for i in all_index if len(i)==1]
    M_new = M.matrix(all_index,None)
    return M_new


def average_genes(M):
    name2indexes = get_duplicated_genes(M)
    for key in name2indexes.keys():
        if len(name2indexes[key])>1:
            newmatrix=[M._X[i] for i in name2indexes[key]]
            new = jmath.mean_matrix(newmatrix,byrow=0)
            M._X[name2indexes[key][0]]=new
            name2indexes[key]=[name2indexes[key][0]]
    all_index = name2indexes.values()
    all_index.sort()
    all_index = [i[0] for i in all_index if len(i)==1]
    M_new = M.matrix(all_index,None)
    return M_new

def guess_gene_header(M):
    all_platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
    ids = M._row_order
    probe_header = all_platforms[0][0]
    probe_id = M._row_names[probe_header]
    annotate_header = 'Gene_ID'
    value_list = arrayannot.annotate_probes(probe_id, annotate_header)
    value_list = [i.lower() for i in value_list]
    new_ids = ids[:]
    new_ids.remove(probe_header)
    column=[]
    for id in new_ids:
        flag = True
        gene_list = M._row_names[id]
        for gene in gene_list:
            if gene.lower() not in value_list:
                flag = False
                break
        if flag:
            ProbeID = probe_header
            GeneID = id
            return ProbeID, GeneID
    assert flag,'we cannot guess the header of this file'