from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        from genomicode import config
        from Betsy import module_utils
        in_data = antecedents
        species = out_attributes['ref']
        if species == 'hg18':
            ref_file = config.hg18_ref
        elif species == 'hg19':
            ref_file = config.hg19_ref
        elif species == 'dm3':
            ref_file = config.dm3_ref
        elif species == 'mm9':
            ref_file = config.mm9_ref
    
        
        assert os.path.exists(ref_file), 'the ref file %s does not exsits' % ref_file
        #command = ['samtools','mpileup','-uf',ref,single_object.identifier,'|',
        #           'bcftools','view', '-bvcg','-','>',outfile]
        samtools_BIN = config.samtools
        samtools_module = module_utils.which(samtools_BIN)
        assert os.path.exists(samtools_module), \
               'cannot find the %s' % samtools_BIN
        command = [
            samtools_BIN, 'mpileup', '-uf', ref_file, in_data.identifier]
        f = file(outfile, 'w')
        try:
            process = subprocess.Popen(
                command, shell=False, stdout=f, stderr=subprocess.PIPE)
        finally:
            f.close()
    
        
        error_message = process.communicate()[1]
        if 'error' in error_message:
            raise ValueError(error_message)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'mpileup_' + original_file + '.bcf'
        return filename
