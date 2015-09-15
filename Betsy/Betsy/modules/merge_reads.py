from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from Betsy import module_utils

        CHUNK_SIZE = 1024*1024

        module_utils.safe_mkdir(out_path)

        fastq_node, group_node = antecedents
        fastq_path = fastq_node.identifier
        sample_group_file = group_node.identifier

        module_utils.assert_sample_group_file(sample_group_file, fastq_path)
        x = module_utils.read_sample_group_file(group_node.identifier)
        x = module_utils.fix_sample_group_filenames(x, fastq_path)
        sample_groups = x

        # The new files should be named:
        # <Sample>.fastq          # if single end
        # <Sample>_<Pair>.fastq   # if paired end
        for x in sample_groups:
            file_, sample, pair = x
            in_filename = os.path.join(fastq_path, file_)
            assert os.path.exists(in_filename)
            
            out_file = "%s.fastq" % sample
            if pair:
                out_file = "%s_%s.fastq" % (sample, pair)
            out_filename = os.path.join(out_path, out_file)

            #if not os.path.exists(out_filename):
            #    # Create an empty outfile that I can append to.
            #    open(out_filename, 'w')
            #else:
            if os.path.exists(out_filename):
                raise NotImplementedError, "Haven't debugged merging."

            
            in_handle = filelib.openfh(in_filename)
            out_handle = open(out_filename, 'a')
            while True:
                x = in_handle.read(CHUNK_SIZE)
                if not x:
                    break
                out_handle.write(x)
        

    def name_outfile(self, antecedents, user_options):
        return "fastq_merged"

