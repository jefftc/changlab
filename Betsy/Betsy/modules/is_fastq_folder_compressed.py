from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        #import shutil
        #from genomicode import filelib
        #from Betsy import module_utils

        os.symlink(in_data.identifier, out_path)

        #filelib.safe_mkdir(out_path)

        #filenames = module_utils.find_fastq_files(antecedents.identifier)
        #assert filenames, "I could not find any FASTQ files."

        # Symlink the files to the new directory.
        #for in_filename in filenames:
        #    in_path, in_file = os.path.split(in_filename)
        #    out_filename = os.path.join(out_path, in_file)
        #    assert not os.path.exists(out_filename)
        #    os.symlink(in_filename, out_filename)
        #    #shutil.copyfile(in_filename, out_filename)

        
    def set_out_attributes(self, antecedents, out_attributes):
        import os
        from Betsy import module_utils
        
        filenames = module_utils.find_fastq_files(antecedents.identifier)
        assert filenames, "I could not find any FASTQ files."

        ext2count = {}  # ".gz" -> num files
        is_compressed = False
        for filename in filenames:
            x, ext = os.path.splitext(filename)
            ext = ext.lower()
            if ext in [".fq", ".fastq"]:
                pass
            elif ext in [".gz", ".bz2", ".xz"]:
                if ext not in ext2count:
                    ext2count[ext] = 0
                ext2count[ext] += 1
            else:
                raise AssertionError, \
                      "Unrecognized extention for fastq file: %s" % ext

        # Figure out which compression was used the most.
        schwartz = [(-ext2count[x], x) for x in ext2count]
        schwartz.sort()

        compressed = "no"
        if schwartz:
            x, ext = schwartz[0]
            assert ext[0] == "."
            compressed = ext[1:]
        out_attributes["compressed"] = compressed
        return out_attributes


    def name_outfile(self, antecedents, user_options):
        return "fastq"
    
