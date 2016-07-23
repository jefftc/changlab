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
        from Betsy import module_utils
        in_data = antecedents
        directory = module_utils.unzip_if_zip(in_data.identifier)
        filenames = os.listdir(directory)
        assert filenames, 'The input folder or zip file is empty.'
        if not os.path.exists(outfile):
            os.mkdir(outfile)
    
        
        samtools_BIN = config.samtools
        assert os.path.exists(samtools_BIN), 'cannot find the %s' % samtools_BIN
        for filename in filenames:
            infile = os.path.join(directory, filename)
            outname = os.path.splitext(filename)[-2] + '.bam'
            outname = os.path.join(outfile, outname)
            command = [samtools_BIN, 'view', '-S', '-b', '-o', outname, infile]
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()
            if 'error' in error_message[1]:
                raise ValueError(error_message)
            assert filelib.exists_nz(outname), (
                'the output file %s for convert_sam_to_bam does not exist' % outname
            )
    

    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        return 'bamFiles_' + original_file
