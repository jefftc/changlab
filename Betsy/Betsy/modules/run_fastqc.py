from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils

        filelib.safe_mkdir(out_path)
        filenames = module_utils.find_fastq_files(in_data.identifier)
        assert filenames, "FASTQ files not found: %s" % in_data.identifier
        metadata = {}

        commands = [
            "fastqc --outdir=%s --extract %s" % (out_path, x)
            for x in filenames]
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        #commands = ["ls > %s" % x for x in filenames]
        parallel.pshell(commands, max_procs=num_cores)

        # Fastqc generates files:
        # <file>_fastqc/
        # <file>_fastqc.zip
        # The contents of the .zip file are identical to the directories.
        # If this happens, then delete the .zip files because they are
        # redundant.
        files = os.listdir(out_path)
        filenames = [os.path.join(out_path, x) for x in files]
        for filename in filenames:
            zip_filename = "%s.zip" % filename
            if os.path.exists(zip_filename):
                os.unlink(zip_filename)

        
    def name_outfile(self, antecedents, user_options):
        return "fastqc_results"


