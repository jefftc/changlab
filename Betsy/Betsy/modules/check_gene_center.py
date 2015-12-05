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
        filename = 'signal_check_center_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        import arrayio
        new_parameters = out_attributes.copy()
        M = arrayio.read(antecedents.identifier)
        if is_gene_center_mean(M):
            new_parameters['gene_center'] = 'mean'
        elif is_gene_center_median(M):
            new_parameters['gene_center'] = 'median'
        else:
            new_parameters['gene_center'] = 'no'
        
        return new_parameters



def is_gene_center_mean(M):
    import numpy
    for line in M.slice():
        if numpy.mean(line) > 0.0000001:
            return False
    
    return 'mean'


def is_gene_center_median(M):
    import numpy
    for line in M.slice():
        if numpy.median(line) > 0.0000001:
            return False
    
    return 'median'
