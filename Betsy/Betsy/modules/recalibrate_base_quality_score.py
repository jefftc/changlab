from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        GATK_path = config.Gatk
        assert GATK_path, 'cannot find the %s' % GATK_path
        species = out_attributes['ref']
        if species == 'hg18':
            ref_file = config.hg18_ref
            dbsnp_file = config.hg18_dbsnp
            indels_file = config.hg18_indels
        elif species == 'hg19':
            ref_file = config.hg19_ref
            dbsnp_file = config.hg19_dbsnp
            indels_file = config.hg19_indels
        else:
            raise ValueError('we cannot process %s' % species)
    
        
        assert os.path.exists(ref_file), 'the ref file %s does not exsits' % ref_file
        assert os.path.exists(
            dbsnp_file), 'the dbsnp file %s does not exsits' % dbsnp_file
        assert os.path.exists(
            indels_file), 'the indels file %s does not exsits' % indels_file
        command = ['java', '-jar', GATK_path, '-T', 'BaseRecalibrator', '-R',
                   ref_file, '-I', in_data.identifier, '-knownSites', dbsnp_file,
                   '-knownSites', indels_file, '-o', 'recal.bam']
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        #error_message = process.communicate()[1]
        #if 'error' in error_message:
        #    raise ValueError(error_message)
        process.wait()
        command = ['java', '-jar', GATK_path, '-T', 'PrintReads', '-R', ref_file,
                   '-BQSR', 'recal.bam', '-I', in_data.identifier, '-o', outfile]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        assert module_utils.exists_nz(outfile), (
            'the output file %s for base_quality_score_recalibration does not exist'
            % outfile)




    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'base_quality_score_recalibration_' + original_file + '.bam'
        return filename

