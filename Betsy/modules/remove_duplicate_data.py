#remove_duplicate_data.py
import os
import module_utils
import arrayio

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    M_new = pick_first_one(M)
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_new,f)
    f.close()
    assert module_utils.exists_nz(outfile),'the output file %s\
                                for remove_duplicate_data fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_remove_duplicate_data_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),'the input file %s\
                          for remove_duplicate_data'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects


def get_duplicated_probe(M):
    probeID = M._row_order[1]
    probe_names = M._row_names[probeID]
    unique_probe=[]
    unique_index=[]
    dup_probe =[]
    dup_dict=dict()
    for i in range(M.nrow()):
        if probe_names[i] not in unique_probe:
            unique_probe.append(probe_names[i])
            unique_index.append(i)
        else:
            dup_probe.append(i)
    for i in dup_probe:
        index = probe_names.index(probe_names[i])
        if index not in dup_dict.keys():
            dup_dict[index]=[index]
        dup_dict[index].append(i)
    return unique_index,dup_dict.values()

def pick_first_one(M):
    unique_index,dup_gene_list = get_duplicated_probe(M)
    probeId = M._row_order[0]
    probe_names = M._row_names[probeId]
    for gene_list in dup_gene_list:
        a = [(probe_names[i],i) for i in gene_list]
        a.sort()
        index = a[0][1]
        first_index = unique_index.index(gene_list[0])
        unique_index[first_index] = index
    M_new = M.matrix(unique_index,None)
    return M_new


