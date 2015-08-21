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
        GATK_path = config.Gatk
        GATK_BIN = module_utils.which(GATK_path)
        assert os.path.exists(GATK_path), ('cannot find the %s' % GATK_path)
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
        command = ['java', '-jar', GATK_path, '-T', 'UnifiedGenotyper', '-R',
                   ref_file, '-I', in_data.identifier, '-o', outfile, '-rf',
                   'BadCigar', '-stand_call_conf', '50.0', '-stand_emit_conf',
                   '10.0']
        process = subprocess.Popen(
            command, shell=False, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        #process.wait()
        error_message = process.communicate()[1]
        #print error_message
        if 'error' in error_message:
            raise ValueError(error_message)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'GATK_' + original_file + '.vcf'
        return filename

