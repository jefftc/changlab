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
                             for get)unqiue_genes fails'%outfile
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
    filename = 'signal_' + original_file + '.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),'the input file %s\
                   for get_unique_genes does not exist'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects

def average_genes(M):
    unique_index,dup_gene_list = get_duplicated_genes(M)
    for gene_list in dup_gene_list:
        newmatrix=[M._X[i] for i in gene_list]
        new = jmath.mean_matrix(newmatrix,byrow=0)
        M._X[gene_list[0]]=new
    M_new = M.matrix(unique_index,None)
    return M_new

def get_high_variance(M):
    unique_index,dup_gene_list = get_duplicated_genes(M)
    for gene_list in dup_gene_list:
        a=[(jmath.var(M._X[i]),i) for i in gene_list]
        a.sort()
        index = a[-1][1]
        first_index = unique_index.index(gene_list[0])
        unique_index[first_index] = index
    M_new = M.matrix(unique_index,None)
    return M_new

def pick_first_one(M):
    unique_index,dup_gene_list = get_duplicated_genes(M)
    ProbeId,GeneID = guess_gene_header(M)
    probe_names = M._row_names[ProbeId]
    for gene_list in dup_gene_list:
        a = [(probe_names[i],i) for i in gene_list]
        a.sort()
        index = a[0][1]
        first_index = unique_index.index(gene_list[0])
        unique_index[first_index] = index
    M_new = M.matrix(unique_index,None)
    return M_new

def get_duplicated_genes(M):
    ProbeId,GeneID = guess_gene_header(M)
    gene_names = M._row_names[GeneID]
    unique_genes=[]
    unique_index=[]
    dup_genes =[]
    dup_dict=dict()
    for i in range(M.nrow()):
        b = re.match("^[-]*",gene_names[i])
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
        index = gene_names.index(gene_names[i])
        if index not in dup_dict.keys():
            dup_dict[index]=[index]
        dup_dict[index].append(i)
    return unique_index,dup_dict.values()

def guess_gene_header(M):
    headers = M._row_order
    if headers == ['Probe.Set.ID','NAME']:# the signal file header when running rma or mas5 preprocess
        ProbeID = 'Probe.Set.ID'
        GeneID = 'NAME'
        return ProbeID,GeneID
    elif headers == ['NAME','Description']: #the control file header when runnning illumina 
        ProbeID = 'NAME'
        GeneID = 'Description'
        return ProbeID,GeneID
    elif headers == ['GeneID','NAME']: #the signal file header when running illumina
        ProbeID = 'GeneID'
        GeneID ='NAME'
        return ProbeID,GeneID
    #try to guess the header
    for header in headers: 
        if header in ['Gene ID','Gene.ID','GeneID',
                          'GeneSymbol','Gene Symbol','Gene.Symbol']:
            GeneID = header
            ProbeID = headers[0]
            return ProbeID,GeneID
    if len(headers) == 2:
        newlist1 = list(set(M._row_names[headers[0]]))
        newlist2 = list(set(M._row_names[headers[1]]))
        #the first one is unique and second one is not unique
        if len(newlist1) == M.nrow() and len(newlist2) < M.nrow(): 
            ProbeID = headers[0]
            GeneID = headers[1]
            return ProbeID,GeneID
        #the first second is unique and first one is not unique
        elif len(newlist2) == M.nrow() and len(newlist1) < M.nrow():
            ProbeID = headers[1]
            GeneID = headers[0]
            return ProbeID,GeneID
        else:
            raise AssertionError, 'we cannot guess the gene header'
    else:   
        raise AssertionError, 'we cannot guess the gene header'
