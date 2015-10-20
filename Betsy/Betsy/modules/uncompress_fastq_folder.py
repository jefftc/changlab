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
        assert filenames, "I could not find any FASTQ files."

        REMOVE = [".gz", ".bz2", ".xz"]

        # Uncompress the files to the new directory in parallel.
        commands = []
        for in_filename in filenames:
            in_path, in_file = os.path.split(in_filename)
            x = in_file
            for r in REMOVE:
                if x.lower().endswith(r):
                    x = x[:-len(r)]
            out_file = x
            out_filename = os.path.join(out_path, out_file)
            
            args = in_filename, out_filename
            keywds = {}
            x = uncompress_file, args, keywds
            commands.append(x)
        parallel.run(commands, num_cores)

        
    def name_outfile(self, antecedents, user_options):
        return "fastq_uncompressed"


def uncompress_file(in_filename, out_filename):
    from genomicode import filelib
    CHUNK_SIZE = 16*1024*1024
    
    in_handle = filelib.openfh(in_filename)
    out_handle = open(out_filename, 'w')
    while True:
        x = in_handle.read(CHUNK_SIZE)
        if not x:
            break
        out_handle.write(x)
