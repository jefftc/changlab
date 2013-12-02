#select_probe_by_best_match.py
import os
import arrayio
from genomicode import filelib,arrayplatformlib, config
from Betsy import bie
from Betsy import rulebase
from Betsy import module_utils

def run(data_node, parameters, network):
    outfile = name_outfile(data_node)
    mapfile = config.HumanHT_12_to_HG_u133_Plus_2
    assert os.path.exists(mapfile),'mapping file %s does not exist'%mapfile
    result=[]
    for d in filelib.read_row(mapfile,header=True):
        if int(d.Distance)<=1000 and d.Match=='Best for Both':
            result.append((d.Affymetrix_Probe_Set_ID,d.Illumina_Probe_ID))
    M = arrayio.read(data_node.attributes['filename'])
    #platform_list = arrayplatformlib.identify_all_platforms_of_matrix(M)
    platform_list = arrayplatformlib.score_all_platforms_of_matrix(M)
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
        new_parameters = parameters.copy()
        new_parameters['filename'] = os.path.split(outfile)[-1]
        out_node = bie.Data(rulebase.SignalFile2,**new_parameters)
        return out_node
    else:
        return None

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'signal_best_match_both_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,data_node):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)


