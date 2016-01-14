from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils

        fastq_node, group_node = antecedents
        fastq_path = fastq_node.identifier
        sample_group_file = group_node.identifier
        filelib.safe_mkdir(out_path)

        module_utils.assert_sample_group_file(sample_group_file, fastq_path)
        x = module_utils.read_sample_group_file(group_node.identifier)
        x = module_utils.fix_sample_group_filenames(x, fastq_path)
        sample_groups = x

        # For merging, the order of the files in the sample_group_file
        # must be maintainted.  Otherwise, will be merged out of order.
        
        # The new files should be named:
        # <Sample>.fastq          # if single end
        # <Sample>_<Pair>.fastq   # if paired end
        jobs = []  # list of (in_filename, out_filename)
        for x in sample_groups:
            file_, sample, pair = x
            in_filename = os.path.join(fastq_path, file_)
            assert os.path.exists(in_filename)
            
            out_file = "%s.fastq" % sample
            if pair:
                out_file = "%s_%s.fastq" % (sample, pair)
            out_filename = os.path.join(out_path, out_file)
            x = in_filename, out_filename
            jobs.append(x)

        out2ins = {}  # out_filename -> list of in_filenames
        for (in_filename, out_filename) in jobs:
            if out_filename not in out2ins:
                out2ins[out_filename] = []
            out2ins[out_filename].append(in_filename)

        commands = []
        for out_filename, in_filenames in out2ins.iteritems():
            args = in_filenames, out_filename
            keywds = {}
            x = merge_or_symlink_files, args, keywds
            commands.append(x)

        parallel.run(commands, num_cores)

    def name_outfile(self, antecedents, user_options):
        return "fastq_merged"


def merge_or_symlink_files(in_filenames, out_filename):
    # If only 1 file, then just symlink it rather than copy.
    # out_filename must not exist.
    import os
    from genomicode import filelib
    
    CHUNK_SIZE = 1024*1024
    assert not os.path.exists(out_filename)
    
    # If only one file, then just symlink it.
    if len(in_filenames) == 1:
        in_filename = in_filenames[0]
        os.symlink(in_filename, out_filename)
        return

    # Create an empty outfile that I can append to.
    open(out_filename, 'w')

    # Append the files in order.
    for in_filename in in_filenames:
        in_handle = filelib.openfh(in_filename)
        out_handle = open(out_filename, 'a')
        while True:
            x = in_handle.read(CHUNK_SIZE)
            if not x:
                break
            out_handle.write(x)
