from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from Betsy import module_utils
        in_data = antecedents
        import arrayio
        f_out = file(outfile, 'w')
        M = arrayio.read(in_data.identifier)
        I_good = []
        #get the percentage of gene filter
        percent = float(user_options['filter_value']) / 100
        for i in range(M.dim()[0]):
            missing_count = 0
            for j in range(M.dim()[1]):
                if M._X[i][j] in [None, 'NA']:
                    missing_count = missing_count + 1
            if float(missing_count) / M.dim()[1] < percent:
                I_good.append(i)
    
        
        M_c = M.matrix(I_good, None)
        arrayio.tab_delimited_format.write(M_c, f_out)
        f_out.close()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for gene_filter fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_filter_' + original_file + '.tdf'
        return filename



