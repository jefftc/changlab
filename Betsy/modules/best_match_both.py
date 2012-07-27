#best_match_both.py
import os
import module_utils
import arrayio
import Betsy_config

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    mapfile = Betsy_config.MAPPING
    assert os.path.exists(mapfile),'mapping file %s does not exist'%mapfile
    f = file(mapfile,'r')
    lines = f.readlines()
    f.close()
    #lines = text.split('\r')
    lines = [line for line in lines if len(line)>0]
    result=[]
    for line in lines[1:]:
        words = line.split('\t')
        if int(words[5])<=1000 and 'Best for Both' in words[6]:
            result.append((words[1],words[3]))
    M = arrayio.read(single_object.identifier)
    headers = M._row_order
    probe_id = M._row_names[headers[0]]
    illu_id = M._row_names[headers[1]]
    index = []
    for i in range(M.nrow()):
        if (probe_id[i],illu_id[i]) in result:
            index.append(i)
    M_new = M.matrix(index,None)
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_new,f)
    f.close()
    assert module_utils.exists_nz(outfile),'the output file %s\
                                for best_match_both fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_best_match_both_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),'the input file %s\
                          for best_match_both'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
