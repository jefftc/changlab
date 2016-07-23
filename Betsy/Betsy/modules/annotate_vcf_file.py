from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import subprocess
        from genomicode import config
        from genomicode import filelib
        from Betsy import module_utils
        in_data = antecedents
        species = out_attributes['ref']
        annotate_BIN = config.annotate_vcf
        command = ['python', annotate_BIN, in_data.identifier, '-o', outfile,
                   '-species', species]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()[1]
        if 'error' in error_message:
            raise ValueError(error_message)
    
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for annot_vcf_file fails' % outfile
        )




    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'vcf_annot_' + original_file + '.txt'
        return filename
