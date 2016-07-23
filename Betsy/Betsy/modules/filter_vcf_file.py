from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import subprocess
        from genomicode import filelib
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        bcftools_BIN = config.bcftools
        bcftools_module = module_utils.which(bcftools_BIN)
        assert bcftools_module, 'cannot find the %s' % bcftools_BIN
        vcfutils_BIN = config.vcfutils
        #vcfutils_module = module_utils.which(vcfutils_BIN)
        #assert bcftools_module, 'cannot find the %s' % bcftools_BIN
        command = [bcftools_BIN, 'view', in_data.identifier, '|', vcfutils_BIN,
                   'varFilter', '-D500']
        #command = ['vcfutils.pl','varFilter','-D100',single_object.identifier]
        f = file(outfile, 'w')
        try:
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=f,
                                       stderr=subprocess.PIPE)
        finally:
            f.close()
    
        
        error_message = process.communicate()[1]
        if 'error' in error_message:
            raise ValueError(error_message)
    
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for filter_vcf_file does not exist' % outfile
        )

    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'filter_' + original_file + '.vcf'
        return filename

