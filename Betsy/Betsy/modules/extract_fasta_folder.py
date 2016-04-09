from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import parallel
        from genomicode import filelib
        from Betsy import module_utils as mlib

        bam_filenames = mlib.find_bam_files(in_data.identifier)
        filelib.safe_mkdir(out_path)

        metadata = {}
        metadata["tool"] = "bam2fastx (unknown version)"

        jobs = []
        for bam_filename in bam_filenames:
            p, f, e = mlib.splitpath(bam_filename)
            fa_filename = os.path.join(out_path, "%s.fa" % f)
            jobs.append((bam_filename, fa_filename))

        bam2fastx = mlib.findbin("bam2fastx")

        commands = []
        for x in jobs:
            bam_filename, fa_filename = x
            # bam2fastx -A --fasta -o rqc14.fa rqc11.bam
            x = [
                mlib.sq(bam2fastx),
                "-A",
                "--fasta",
                "-o", mlib.sq(fa_filename),
                mlib.sq(bam_filename),
                ]
            x = " ".join(x)
            commands.append(x)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        parallel.pshell(commands, max_procs=num_cores)

        x = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "fasta"
