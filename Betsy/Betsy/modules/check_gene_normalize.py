from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """check gene normalize"""
        import shutil
        from Betsy import module_utils
        in_data = antecedents
        #out_attributes = set_out_attributes(in_data, out_attributes)
        shutil.copyfile(in_data.identifier, outfile)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for check_gene_normalize fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_check_normalize_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        import arrayio
        new_parameters = out_attributes.copy()
        M = arrayio.read(antecedents.identifier)
        if is_gene_normalize_varaince(M):
            new_parameters['gene_normalize'] = 'variance'
        elif is_gene_normalize_ss(M):
            new_parameters['gene_normalize'] = 'sum_of_squares'
        else:
            new_parameters['gene_normalize'] = 'no'
        
        return new_parameters



def is_gene_normalize_varaince(M):
    import numpy
    for line in M.slice():
        if abs(numpy.var(line) - 1) > 0.000001:
            return False
    
    return True


def is_gene_normalize_ss(M):
    import numpy
    for line in M.slice():
        if abs(numpy.sum([(x - numpy.mean(line)) ** 2
                          for x in line]) - 1) > 0.000001:
            return False
    
    return True
