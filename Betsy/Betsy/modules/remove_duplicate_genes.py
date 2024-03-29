from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import arrayplatformlib as apl
        from Betsy import module_utils as mlib

        in_data = antecedents
        metadata = {}

        M = arrayio.read(in_data.identifier)
        cat2header = apl.categorize_headers(M)
        header = cat2header.get(apl.GENE_SYMBOL)
        if header is None:
            header = cat2header.get(apl.GENE_ID)
        assert header is not None, "I could not find gene IDs or symbols: %s" \
               % in_data.identifier
        metadata["dedup_header"] = header

        slice_matrix = mlib.get_config("slice_matrix", which_assert_file=True)

        sq = parallel.quote
        algorithm = out_attributes['unique_genes']
        if algorithm == "average_genes":
            raise NotImplementedError
        elif algorithm == "high_var":
            dedup_cmd = ["--dedup_row_by_var", sq(header)]
            pass
        elif algorithm == "first_gene":
            raise NotImplementedError
        else:
            raise AssertionError, "Unknown algorithm: %s" % algorithm

        cmd = [
            sq(slice_matrix),
            ]
        cmd += dedup_cmd
        cmd += [
            sq(in_data.identifier)
            ]
        cmd = " ".join(cmd)
        cmd = "%s >& %s" % (cmd, outfile)
        parallel.sshell(cmd)

        filelib.assert_exists_nz(outfile)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "signal.tdf"



## def get_duplicated_genes(M):
##     import re
##     ProbeId, GeneID = guess_gene_header(M)
##     gene_names = M._row_names[GeneID]
##     name2indexes = dict()
##     for i in range(M.nrow()):
##         b = re.match("^[-]*", gene_names[i])
##         if len(gene_names[i]) == 0:
##             continue
##         elif len(b.group()) > 0:
##             continue
##         else:
##             if gene_names[i] not in name2indexes:
##                 name2indexes[gene_names[i]] = []
##             name2indexes[gene_names[i]].append(i)
    
##     return name2indexes


## def pick_first_one(M):
##     name2indexes = get_duplicated_genes(M)
##     ProbeId, GeneID = guess_gene_header(M)
##     probe_names = M._row_names[ProbeId]
##     for key in name2indexes.keys():
##         if len(name2indexes[key]) > 1:
##             a = [(probe_names[i], i) for i in name2indexes[key]]
##             a.sort()
##             index = a[0][1]
##             name2indexes[key] = [index]
##     all_index = name2indexes.values()
##     all_index.sort()
##     all_index = [i[0] for i in all_index if len(i) == 1]
##     M_new = M.matrix(all_index, None)
##     return M_new


## def get_high_variance(M):
##     from genomicode import jmath
##     name2indexes = get_duplicated_genes(M)
##     for key in name2indexes.keys():
##         if len(name2indexes[key]) > 1:
##             a = [(jmath.var(M._X[i]), i) for i in name2indexes[key]]
##             a.sort()
##             index = a[-1][1]
##             name2indexes[key] = [index]
    
##     all_index = name2indexes.values()
##     all_index.sort()
##     all_index = [i[0] for i in all_index if len(i) == 1]
##     M_new = M.matrix(all_index, None)
##     return M_new


## def average_genes(M):
##     from genomicode import jmath
##     name2indexes = get_duplicated_genes(M)
##     for key in name2indexes.keys():
##         if len(name2indexes[key]) > 1:
##             newmatrix = [M._X[i] for i in name2indexes[key]]
##             new = jmath.mean_matrix(newmatrix, byrow=0)
##             M._X[name2indexes[key][0]] = new
##             name2indexes[key] = [name2indexes[key][0]]
    
##     all_index = name2indexes.values()
##     all_index.sort()
##     all_index = [i[0] for i in all_index if len(i) == 1]
##     M_new = M.matrix(all_index, None)
##     return M_new


## def guess_gene_header(M):
##     from genomicode import arrayannot
##     from genomicode import arrayplatformlib
    
##     all_platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
##     ids = M._row_order
##     probe_header = all_platforms[0][0]
##     probe_id = M._row_names[probe_header]
##     annotate_header = 'Gene ID'
##     value_list = arrayannot.annotate_probes(probe_id, annotate_header)
##     value_list = [i.lower() for i in value_list]
##     new_ids = ids[:]
##     new_ids.remove(probe_header)
##     #column = []
##     for id in new_ids:
##         flag = True
##         gene_list = M._row_names[id]
##         for gene in gene_list:
##             if gene.lower() not in value_list:
##                 flag = False
##                 break
##         if flag:
##             ProbeID = probe_header
##             GeneID = id
##             return ProbeID, GeneID
    
##     assert flag, 'we cannot guess the header of this file'
