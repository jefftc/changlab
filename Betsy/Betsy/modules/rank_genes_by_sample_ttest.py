from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import numpy
        import arrayio
        from genomicode import jmath
        from Betsy import gene_ranking
        from Betsy import read_label_file
        from Betsy import module_utils as mlib
        
        data_node, cls_node = antecedents
        M = arrayio.read(data_node.identifier)
        x = read_label_file.read(cls_node.identifier)
        label, label_line, second_line = x
        assert len(label) == 2, \
               'the length of label in %s should be 2' % cls_node.identifier
        assert len(label[0]) == 2
        assert len(label[1]) == 2

        threshold = mlib.get_user_option(
            user_options, "gene_select_threshold", not_empty=True, type=float)

        first = M.slice(None, label[0][0])
        second = M.slice(None, label[1][0])
        t, p = gene_ranking.t_test(first, second)
        for i in range(len(p)):
            if p[i] is None:
                p[i] = 1
        
        gene_list = []
        key = M._row_order[0]

        gene_order = out_attributes["gene_order"]
        if gene_order == 'ttest_p':
            sort_p = [(p[index], index) for index in range(len(p))]
            sort_p.sort()
            for i in range(len(sort_p)):
                if sort_p[i][0] < threshold:
                    gene_list.append(M._row_names[key][sort_p[i][1]])
        elif gene_order == 'ttest_fdr':
            #for i in range(len(p)):
            #    if p[i] == 10:
            #        p[i] = ''
            fdr = jmath.cmh_fdr_bh(p)
            #for i in range(len(fdr)):
            #    if numpy.isnan(fdr[i]):
            #        fdr[i] = 10
            sort_fdr = [(fdr[index], index) for index in range(len(fdr))]
            sort_fdr.sort()
            for i in range(len(fdr)):
                if sort_fdr[i][0] < threshold:
                    gene_list.append(M._row_names[key][sort_fdr[i][1]])
        else:
            raise AssertionError, "Unknown gene_order: %s" % gene_order

        assert gene_list, 'there is no significant genes can be found in ttest'

        f = open(outfile, 'w')
        f.write('\t'.join(gene_list))
        f.close()


    def name_outfile(self, antecedents, user_options):
        return "gene_list.txt"
