#select_probe_by_best_match.py
import os
import arrayio
from genomicode import filelib, arrayplatformlib, config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    mapfile = config.HumanHT_12_to_HG_u133_Plus_2
    assert os.path.exists(mapfile), 'mapping file %s does not exist' % mapfile
    result = []
    for d in filelib.read_row(mapfile, header=True):
        if int(d.Distance) <= 1000 and d.Match == 'Best for Both':
            result.append((d.Affymetrix_Probe_Set_ID, d.Illumina_Probe_ID))
    
    M = arrayio.read(in_data.identifier)
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
        if (probe_id[i], illu_id[i]) in result:
            index.append(i)
    
    if len(index) > 0:
        M_new = M.matrix(index, None)
        f = file(outfile, 'w')
        arrayio.tab_delimited_format.write(M_new, f)
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for best_match_both fails' % outfile
        )
        out_node = bie3.Data(rulebase._SignalFile_Filter, **out_attributes)
        out_object = module_utils.DataObject(out_node, outfile)
        return out_object
    else:
        return None



def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_best_match_both_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
