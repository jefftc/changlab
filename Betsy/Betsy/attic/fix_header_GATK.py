from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        from genomicode import filelib
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        AddOrReplaceReadGroups_path = config.AddOrReplaceReadGroups
        assert os.path.exists(AddOrReplaceReadGroups_path), (
            'cannot find the %s' % AddOrReplaceReadGroups_path
        )
        command = ['java', '-Xmx5g', '-jar', AddOrReplaceReadGroups_path,
                   'I=' + in_data.identifier, 'O=' + outfile, 'PL=illumina',
                   'ID=Group1', 'LB=Al_chrom3', 'PU=Al_chrom3', 'SM=Al_chrom3',
                   'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=LENIENT']
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        error_message = process.communicate()[1]
        #if 'error' in error_message:
        #    raise ValueError(error_message)
        assert filelib.exists_nz(outfile), (
            'the output file %s for fix_header_GATK does not exist' % outfile
        )




    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'fix_header_' + original_file + '.bam'
        return filename

