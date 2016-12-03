from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import jmath
        from genomicode import filelib
        in_data = antecedents
        cwd = os.getcwd()
        R = jmath.start_R()
        R('require(limma,quietly=TRUE)')
        R('library(marray)')
        os.chdir(in_data.identifier)
        try:
            R('dir<-getwd()')
            R('files<-list.files(dir)')
            R('x.read<-read.Agilent(files)')
        finally:
            os.chdir(cwd)
    
        
        R('xnorm.loc <- maNorm(x.read, norm = "loess")')
        R('x.norm <- maNormScale(xnorm.loc, norm = "p")')
        tmpfile = 'tmp.txt'
        jmath.R_equals(tmpfile, 'tmpfile')
        R('write.marray(x.norm,tmpfile)')
        f = open(tmpfile, 'r')
        text = f.readlines()
        firstline = text[0].split()
        f.close()
        firstindex = firstline.index('"ProbeName"')
        if '"Sequence"' in firstline:
            secondindex = firstline.index('"Sequence"')
        else:
            secondindex = firstline.index('"ControlType"')
    
        
        sample = range(secondindex + 1, len(firstline))
        f = open(outfile, 'w')
        for i in text:
            line = i.split()
            f.write(line[firstindex] + '\t')
            for j in sample:
                f.write(line[j] + '\t')
            f.write('\n')
    
        
        f.close()
        os.remove(tmpfile)
        assert filelib.exists_nz(outfile), (
            'the output file %s for preprocess_agilent fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_agilent' + original_file + '.tdf'
        return filename


