from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from genomicode import binreg

        M = arrayio.read(antecedents.identifier)
        X = M._X
        for i in range(len(X)):
            for j in range(len(X[i])):
                if X[i][j] is None:
                    continue
                X[i][j] = 2 ** float(X[i][j])

        # If all the values should be integers, convert them back to
        # ints.
        EPS = 1E-4
        # 628494.000012.  1E-5 is too strict.
        is_int = True
        for i in range(len(X)):
            for j in range(len(X[i])):
                if X[i][j] is None:
                    continue
                if abs(X[i][j] - int(round(X[i][j]))) > EPS:
                    is_int = False
                    break
            if not is_int:
                break
        if is_int:
            for i in range(len(X)):
                for j in range(len(X[i])):
                    if X[i][j] is None:
                        continue
                    X[i][j] = int(round(X[i][j]))
                    
        #assert binreg.is_logged_array_data(M), (
        #    'the input file %s should be logged' % antecedents.identitifer)
    
        f = file(outfile, 'w')
        arrayio.tab_delimited_format.write(M, f)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        return 'signal_unlog' + original_file + '.tdf'
