from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """run preprocessdataset """
        from genomicode import filelib
        from Betsy import module_utils
        import arrayio
        in_data = antecedents
        threshold = 20
        ceiling = 16000
        min_fold_change = 5
        min_delta = 100.0
        M = arrayio.read(in_data.identifier)
        X = M.slice()
        I_good = []
        for i in range(M.nrow()):
            for j in range(len(X[i])):
                if X[i][j] < threshold:
                    M._X[i][j] = threshold
                if X[i][j] > ceiling:
                    M._X[i][j] = ceiling
            gene = M._X[i]
            fold_change = max(gene) / float(min(gene))
            delta = max(gene) - min(gene)
            if fold_change >= min_fold_change and delta >= min_delta:
                I_good.append(i)
    
        
        f = file(outfile, 'w')
        M_c = M.matrix(I_good, None)
        arrayio.tab_delimited_format.write(M_c, f)
        f.close()
        assert filelib.exists_nz(outfile), (
            'the output file %s for filter_and_threshold_genes fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_preprocessdataset_' + original_file + '.tdf'
        return filename



