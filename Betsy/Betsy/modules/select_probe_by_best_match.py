#select_probe_by_best_match.py
import os
from Betsy import module_utils
import arrayio
from Betsy import config
from genomicode import filelib,arrayannot

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    mapfile = config.MAPPING
    assert os.path.exists(mapfile),'mapping file %s does not exist'%mapfile
    result=[]
    for d in filelib.read_row(mapfile,header=True):
        if int(d.Distance)<=1000 and d.Match=='Best for Both':
            result.append((d.Affymetrix_Probe_Set_ID,d.Illumina_Probe_ID))
    M = arrayio.read(single_object.identifier)
    platform_list = arrayannot.identify_all_platforms_of_matrix(M)
    illu_id = None
    probe_id = None
    for platform in platform_list:
        if 'HumanHT_12' in platform:
            illu_id = M._row_names[platform[0]]
        if 'HG_U133_Plus_2' in platform:
            probe_id = M._row_names[platform[0]]
    if not illu_id or not probe_id:
        return None
    index = []
    for i in range(M.nrow()):
        if (probe_id[i],illu_id[i]) in result:
            index.append(i)
    if len(index)>0:
        M_new = M.matrix(index,None)
        f = file(outfile,'w')
        arrayio.tab_delimited_format.write(M_new,f)
        f.close()
        assert module_utils.exists_nz(outfile),(
            'the output file %s for best_match_both fails'%outfile)
        new_objects = get_newobjects(parameters,objects,pipeline)
        module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
        return new_objects
    else:
        return None

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_best_match_both_'+original_file+'.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for best_match_both'%single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
