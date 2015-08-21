from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from Betsy import module_utils
        from genomicode import pcalib
        
        in_data = antecedents
        M = arrayio.read(in_data.identifier)
        X = M._X
        N = 500
        if 'pca_gene_num' in user_options:
            N = int(user_options['pca_gene_num'])
    
        
        N = min(N, M.nrow())
        index = pcalib.select_genes_var(X, N)
        M_new = M.matrix(index, None)
        f = file(outfile, 'w')
        arrayio.tab_delimited_format.write(M_new, f)
        f.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for analyze_samples_pca fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'pca_matrix_' + original_file + '.tdf'
        return filename



