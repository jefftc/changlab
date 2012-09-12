#remove_duplicate_data.py
import os
import module_utils
import arrayio
from genomicode import arrayannot

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    M_new = pick_first_one(M)
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_new,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for remove_duplicate_data fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_remove_duplicate_data_'+original_file+'.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for remove_duplicate_data'%single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects


def get_duplicated_probe(M):
    platform_list = arrayannot.identify_all_platforms_of_matrix(M)
    assert len(platform_list)>=2,'the number of header identified should be at least 2'
    if platform_list[0][0] == M._row_order[-1]:
        probeID = platform_list[1][0]
    else:
        probeID = platform_list[0][0]
    probe_names = M._row_names[probeID]
    name2indexes = dict()
    for i in range(M.nrow()):
        if probe_names[i] not in name2indexes:
            name2indexes[probe_names[i]] = []
        name2indexes[probe_names[i]].append(i)
    return name2indexes

def pick_first_one(M):
    name2indexes = get_duplicated_probe(M)
    probeId = M._row_order[-1]
    probe_names = M._row_names[probeId]
    for key in name2indexes.keys():
        if len(name2indexes[key]) > 1:
            a = [(probe_names[i],i) for i in name2indexes[key]]
            a.sort()
            index = a[0][1]
            name2indexes[key] = [index]
    all_index = name2indexes.values()
    all_index.sort()
    all_index = [i[0] for i in all_index if len(i)==1]
    M_new = M.matrix(all_index,None)
    return M_new

##def get_duplicated_probe(M):
##    assert M._row_names>=2,'the number of header should be at least 2'
##    probeID = M._row_order[1]
##    probe_names = M._row_names[probeID]
##    name2indexes = dict()
##    for i in range(M.nrow()):
##        if probe_names[i] not in name2indexes:
##            name2indexes[probe_names[i]] = []
##        name2indexes[probe_names[i]].append(i)
##    return name2indexes
##
##def pick_first_one(M):
##    name2indexes = get_duplicated_probe(M)
##    probeId = M._row_order[0]
##    probe_names = M._row_names[probeId]
##    for key in name2indexes.keys():
##        if len(name2indexes[key]) > 1:
##            a = [(probe_names[i],i) for i in name2indexes[key]]
##            a.sort()
##            index = a[0][1]
##            name2indexes[key] = [index]
##    all_index = name2indexes.values()
##    all_index.sort()
##    all_index = [i[0] for i in all_index if len(i)==1]
##    M_new = M.matrix(all_index,None)
##    return M_new

