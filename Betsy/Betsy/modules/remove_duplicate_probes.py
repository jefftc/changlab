from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from genomicode import filelib
        from Betsy import module_utils
        in_data = antecedents
        M = arrayio.read(in_data.identifier)
        M_new = remove_duplicate_probes_var(M)
        f = file(outfile, 'w')
        arrayio.tab_delimited_format.write(M_new, f)
        f.close()
        assert filelib.exists_nz(outfile), (
            'the output file %s for remove_duplicate_probes fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_select_probe_var_' + original_file + '.tdf'
        return filename



def get_high_variance(M, name2indexes, empty_item):
    from genomicode import jmath
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
    from genomicode import arrayplatformlib
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
