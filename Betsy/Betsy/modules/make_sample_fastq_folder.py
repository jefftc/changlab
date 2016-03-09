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
        from Betsy import module_utils as mlib

        # This this is I/O heavy, don't use so many cores.
        MAX_CORES = 2

        filenames = mlib.find_fastq_files(in_data.identifier)
        assert filenames, "I could not find any FASTQ files."
        filelib.safe_mkdir(out_path)
        metadata = {}
        
        num_samples = mlib.get_user_option(
            user_options, "num_samples", not_empty=True, type=int)
        metadata["num_samples"] = num_samples

        jobs = []
        for in_filename in filenames:
            p, f = os.path.split(in_filename)
            out_filename = os.path.join(out_path, f)
            x = in_filename, out_filename
            jobs.append(x)

        cmds = []
        for x in jobs:
            in_filename, out_filename = x
            x = copy_fastq_file, (in_filename, out_filename, num_samples), {}
            cmds.append(x)

        nc = min(MAX_CORES, num_cores)
        metadata["num cores"] = nc
        parallel.pyfun(cmds, num_procs=nc)

        return metadata

        
    def name_outfile(self, antecedents, user_options):
        return "sample.fastq"


def copy_fastq_file(in_filename, out_filename, num_samples):
    from genomicode import genomelib

    outhandle = open(out_filename, 'w')
    for i, x in enumerate(genomelib.read_fastq(in_filename)):
        if i >= num_samples:
            break
        genomelib.write_fastq(*x, handle=outhandle)
