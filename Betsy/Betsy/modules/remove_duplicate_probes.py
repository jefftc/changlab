from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import shutil
        import arrayio
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils as mlib

        filename = in_data.identifier
        filelib.assert_exists_nz(filename)

        # De-duplicate by every single header.  Not sure if this is
        # right.
        MATRIX = arrayio.read(filename)
        # Figure out which columns has duplicates.
        has_dup = []
        for name in MATRIX.row_names():
            annots = MATRIX.row_names(name)
            assert name not in has_dup
            seen = {}
            for annot in annots:
                if annot in seen:
                    has_dup.append(name)
                    break
                seen[annot] = 1
        if not has_dup:
            shutil.copy2(filename, outfile)
            return
        
        sq = parallel.quote
        slice_matrix = mlib.get_config("slice_matrix", which_assert_file=True)
        for i, name in enumerate(has_dup):
            f = "outfile.%d.txt" % i
            x = [
                sq(slice_matrix),
                "--dedup_row_by_var", sq(name),
                sq(filename),
                ">&",
                sq(f),
                ]
            x = " ".join(map(str, x))
            parallel.sshell(x)
        shutil.copy2(f, outfile)
            


        # slice_matrix.py --dedup_row_by_var "gene_short_name" exp31.txt > exp32.txt

        #M = arrayio.read(in_data.identifier)
        #M_new = remove_duplicate_probes_var(M)
        #arrayio.tab_delimited_format.write(M_new, open(outfile, 'w'))
        
        #assert filelib.exists_nz(outfile), (
        #    'the output file %s for remove_duplicate_probes fails' % outfile
        #)


    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'signal_select_probe_var_' + original_file + '.tdf'
        #return filename
        return "signal.tdf"



## def get_high_variance(M, name2indexes, empty_item):
##     from genomicode import jmath
##     for key in name2indexes.keys():
##         if len(name2indexes[key]) > 1:
##             a = [(jmath.var(M._X[i]), i) for i in name2indexes[key]]
##             a.sort()
##             index = a[-1][1]
##             name2indexes[key] = [index]
    
##     all_index = name2indexes.values()
##     all_index.sort()
##     all_index = [i[0] for i in all_index if len(i) == 1]
##     all_index.extend(empty_item)
##     all_index.sort()
##     M_new = M.matrix(all_index, None)
##     return M_new


## def remove_duplicate_probes_var(M):
##     from genomicode import arrayplatformlib
##     probe_headers = arrayplatformlib.identify_all_platforms_of_matrix(M)
##     ids = [probe_header[0] for probe_header in probe_headers]
##     for id in ids:
##         probe_names = M._row_names[id]
##         name2indexes = dict()
##         empty_item = []
##         for i in range(M.nrow()):
##             if len(probe_names[i]) == 0 or probe_names[i] == 'NA':
##                 empty_item.append(i)
##             elif probe_names[i] in name2indexes:
##                 name2indexes[probe_names[i]].append(i)
##             elif probe_names[i] not in name2indexes:
##                 name2indexes[probe_names[i]] = [i]
##         M = get_high_variance(M, name2indexes, empty_item)
##     return M
