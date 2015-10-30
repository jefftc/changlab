from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import shutil
        shutil.copyfile(in_data.identifier, outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_check_normalize_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        import arrayio
        new_parameters = out_attributes.copy()
        M = arrayio.read(antecedents.identifier)
        if is_gene_normalize_variance(M):
            new_parameters['gene_normalize'] = 'variance'
        elif is_gene_normalize_ss(M):
            new_parameters['gene_normalize'] = 'sum_of_squares'
        else:
            new_parameters['gene_normalize'] = 'no'
        
        return new_parameters



def is_gene_normalize_variance(M):
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
