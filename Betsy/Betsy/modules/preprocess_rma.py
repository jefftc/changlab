from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """preprocess the inputfile with RMA 
           using preprocess.py will generate a output file"""
        import os
        import subprocess
        from Betsy import module_utils
        from genomicode import config
        from genomicode import filelib
        in_data = antecedents
        #preprocess the cel file to text signal file
        PREPROCESS_path = config.preprocess
        PREPROCESS_BIN = module_utils.which(PREPROCESS_path)
        assert PREPROCESS_BIN, 'cannot find the %s' % PREPROCESS_path
        command = ['python', PREPROCESS_BIN, 'RMA', in_data.identifier]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            if not "Loading required package: Biobase" in error_message:
                raise ValueError(error_message)
    
        
        outputfiles = os.listdir(".")
        outputfile = None
        for i in outputfiles:
            if i.endswith('.rma'):
                outputfile = i
        
        assert outputfile, "No output file created."

        os.rename(outputfile, outfile)
        assert filelib.exists_nz(outfile), (
            'the output file %s for preprocess_rma fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_rma_' + original_file + '.jeffs'
        return filename


