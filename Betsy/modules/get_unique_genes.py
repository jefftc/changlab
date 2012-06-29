#get_unique_genes.py
import os
import module_utils
from genomicode import jmath
import arrayio
import re

def run(parameters,objects,pipeline):
    """get unique genes"""
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
    arrayio.pcl_format.write(M_new,f)
    f.close()
    assert module_utils.exists_nz(outfile),'the output file %s\
                             for convert_pcl_gct fails'%outfile
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
    filename = 'signal_' + original_file + '.gct'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),'the input file %s\
                   for convert_pcl_gct does not exist'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects

def average_genes(M):
    a = M._row_order
    gene_header = a[1]
    gene_names = M._row_names[gene_header]
    index = range(M.nrow())
    unique_genes=[]
    dup_genes =[]
    unique_index=[]
    dup_dict=dict()
    for i in range(M.nrow()):
        b=re.match("^[-]*",gene_names[i])
        if len(gene_names[i])==0:
            continue
        elif len(b.group())>0:
            continue
        elif gene_names[i] not in unique_genes:
            unique_genes.append(gene_names[i])
            unique_index.append(i)
        else:
            dup_genes.append(i)
    for i in dup_genes:
            index=gene_names.index(gene_names[i])
            if index not in dup_dict.keys():
                dup_dict[index]=[index]
            dup_dict[index].append(i)
    for key in dup_dict.keys():
        newmatrix=[M._X[i] for i in dup_dict[key]]
        new = jmath.mean_matrix(newmatrix,byrow=0)
        M._X[key]=new
    M_new = M.matrix(unique_index,None)
    return M_new

def get_high_variance(M):
    a = M._row_order
    gene_header = a[1]
    gene_names = M._row_names[gene_header]
    index = range(M.nrow())
    unique_genes=[]
    dup_genes =[]
    unique_index=[]
    dup_dict=dict()
    for i in range(M.nrow()):
        b=re.match("^[-]*",gene_names[i])
        if len(gene_names[i])==0:
            continue
        elif len(b.group())>0:
            continue
        elif gene_names[i] not in unique_genes:
            unique_genes.append(gene_names[i])
            unique_index.append(i)
        else:
            dup_genes.append(i)
    for i in dup_genes:
            index=gene_names.index(gene_names[i])
            if index not in dup_dict.keys():
                dup_dict[index]=[index]
            dup_dict[index].append(i)
    for key in dup_dict.keys():
        a=[(jmath.var(M._X[i]),i) for i in dup_dict[key]]
        a.sort()
        dup_dict[key]=a[-1][1]
    newindex=[]
    for i in unique_index:
        if i in dup_dict.keys():
            newindex.append(dup_dict[i])
        else:
            newindex.append(i)
    M_new = M.matrix(newindex,None)
    return M_new

def pick_first_one(M):
    a = M._row_order
    gene_header = a[1]
    probe_header = a[0]
    gene_names = M._row_names[gene_header]
    index = range(M.nrow())
    unique_genes=[]
    dup_genes =[]
    unique_index=[]
    dup_dict=dict()
    for i in range(M.nrow()):
        b=re.match("^[-]*",gene_names[i])
        if len(gene_names[i])==0:
            continue
        elif len(b.group())>0:
            continue
        elif gene_names[i] not in unique_genes:
            unique_genes.append(gene_names[i])
            unique_index.append(i)
        else:
            dup_genes.append(i)
    for i in dup_genes:
            index=gene_names.index(gene_names[i])
            if index not in dup_dict.keys():
                dup_dict[index]=[index]
            dup_dict[index].append(i)
    probe_names = M._row_names[probe_header]
    for key in dup_dict.keys():
        a = [(probe_names[i],i) for i in dup_dict[key]]
        a.sort()
        dup_dict[key]=a[0][1]
    newindex = []
    for i in unique_index:
        if i in dup_dict.keys():
            newindex.append(dup_dict[i])
        else:
            newindex.append(i)
    M_new = M.matrix(newindex,None)
    return M_new
