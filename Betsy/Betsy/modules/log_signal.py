from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        """log the input file"""
        import math
        import arrayio
        from genomicode import filelib
        from genomicode import binreg
        from Betsy import module_utils

        signal_file = in_data.identifier
        filelib.assert_exists_nz(signal_file)
        
        M = arrayio.read(signal_file)
        assert not binreg.is_logged_array_data(M), 'the file is logged'
        # Change the matrix in place.
        X = M._X
        for i in range(len(X)):
            for j in range(len(X[i])):
                x = X[i][j]
                if x is None:
                    continue
                x = float(x)
                if x < 1:
                    x = 1
                x = math.log(x, 2)
                X[i][j] = x
    
        M_c = arrayio.convert(M, to_format=arrayio.tab_delimited_format)
        
        handle = open(outfile, 'w')
        arrayio.tab_delimited_format.write(M_c, handle)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_log_' + original_file + '.tdf'
        return filename




