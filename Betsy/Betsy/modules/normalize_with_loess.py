from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import smarray
        from genomicode import filelib
        from Betsy import gpr_module
        from Betsy import module_utils
        
        in_data = antecedents
        #out_attributes = set_out_attributes(in_data, out_attributes)
        filenames = os.listdir(in_data.identifier)
        keep = []
        red_sig_matrix = []
        green_sig_matrix = []
        red_back_matrix = []
        green_back_matrix = []
        sample = []
        for filename in filenames:
            fileloc = os.path.join(in_data.identifier, filename)
            if not filename.endswith('gpr.gz') and not filename.endswith('gpr'):
                continue
            if filename.endswith('gpr.gz'):
                samplename = filename[:-7]
            elif filename.endswith('gpr'):
                samplename = filename.s[:-4]
            sample.append(samplename)
            x = gpr_module.extract_multiple(fileloc, keep)
            red_sig, green_sig, red_back, green_back, keep = x
            red_sig_matrix.append(red_sig)
            green_sig_matrix.append(green_sig)
            red_back_matrix.append(red_back)
            green_back_matrix.append(green_back)
    
        
        red_signal = [[0] * len(red_sig_matrix)
                      for i in range(len(red_sig_matrix[0]))]
        green_signal = [[0] * len(red_sig_matrix)
                        for i in range(len(red_sig_matrix[0]))]
        red_back = [[0] * len(red_sig_matrix)
                    for i in range(len(red_sig_matrix[0]))]
        green_back = [[0] * len(red_sig_matrix)
                      for i in range(len(red_sig_matrix[0]))]
        for i in range(len(red_sig_matrix[0])):
            for j in range(len(red_sig_matrix)):
                red_signal[i][j] = red_sig_matrix[j][i]
                green_signal[i][j] = green_sig_matrix[j][i]
                red_back[i][j] = red_back_matrix[j][i]
                green_back[i][j] = green_back_matrix[j][i]
    
        
        x = smarray.correct_background(
            red_signal, green_signal, red_back, green_back,
            method="normexp", offset=50)
        R, G = x
        SIGNAL = smarray.normalize_within_arrays(R, G, method="loess")
        keep[0][1] = keep[0][1].upper()
          #convert the 'Name' to 'NAME'
        f = open(outfile, 'w')
        f.write('\t'.join(keep[0][0:2]))
        f.write('\t')
        f.write('\t'.join(sample))
        f.write('\n')
        for i in range(len(keep) - 1):
            f.write('\t'.join(keep[i + 1][0:2]))
            for j in range(len(SIGNAL[0])):
                f.write('\t')
                f.write(str(SIGNAL[i][j]))
            f.write('\n')
        f.close()
        assert filelib.exists_nz(outfile), (
            'the output file %s for loess fails' % outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        return 'signal_loess_' + original_file + '.tdf'



