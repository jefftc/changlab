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
        from genomicode import filelib
        from genomicode import filelib
        from Betsy import module_utils
        
        in_data = antecedents
        #out_attributes = set_out_attributes(in_data, out_attributes)
        species = out_attributes['ref']
        if species == 'hg18':
            ref_file = config.hg18_ref
        elif species == 'hg19':
            ref_file = config.hg19_ref
        elif species == 'dm3':
            ref_file = config.dm3_ref
        elif species == 'mm9':
            ref_file = config.mm9_ref
        else:
            raise ValueError('cannot handle %s' % species)
    
        assert os.path.exists(ref_file), "File not found: %s" % ref_file
        bwa_BIN = config.bwa
        bwa_module = module_utils.which(bwa_BIN)
        assert bwa_module, 'cannot find the %s' % bwa_BIN
        command = [bwa_BIN, 'aln', ref_file, in_data.identifier]
        f = file(outfile, 'w')
        try:
            process = subprocess.Popen(
                command, shell=False, stdout=f, stderr=subprocess.PIPE)
        finally:
            f.close()
    
        error_message = process.communicate()[1]
        if 'error' in error_message:
            raise ValueError(error_message)
    
        assert filelib.exists_nz(outfile), (
            'the output file %s for align_sequence does not exist' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        return 'align_sequence' + original_file
