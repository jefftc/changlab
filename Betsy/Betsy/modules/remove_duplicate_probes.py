#remove_duplicate_probes.py
import os
import arrayio
from genomicode import jmath, arrayplatformlib
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    M = arrayio.read(in_data.identifier)
    M_new = remove_duplicate_probes_var(M)
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M_new, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for remove_duplicate_probes fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Filter, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes, pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_select_probe_var_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def get_high_variance(M, name2indexes, empty_item):
    for key in name2indexes.keys():
        if len(name2indexes[key]) > 1:
            a = [(jmath.var(M._X[i]), i) for i in name2indexes[key]]
            a.sort()
            index = a[-1][1]
            name2indexes[key] = [index]
    all_index = name2indexes.values()
    all_index.sort()
    all_index = [i[0] for i in all_index if len(i) == 1]
    all_index.extend(empty_item)
    all_index.sort()
    M_new = M.matrix(all_index, None)
    return M_new


def remove_duplicate_probes_var(M):
    probe_headers = arrayplatformlib.identify_all_platforms_of_matrix(M)
    ids = [probe_header[0] for probe_header in probe_headers]
    for id in ids:
        probe_names = M._row_names[id]
        name2indexes = dict()
        empty_item = []
        for i in range(M.nrow()):
            if len(probe_names[i]) == 0 or probe_names[i] == 'NA':
                empty_item.append(i)
            elif probe_names[i] in name2indexes:
                name2indexes[probe_names[i]].append(i)
            elif probe_names[i] not in name2indexes:
                name2indexes[probe_names[i]] = [i]
        M = get_high_variance(M, name2indexes, empty_item)
    return M
