from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib

        # If align_with_star is run with two_pass=yes, this will leave
        # two BAM files for every sample.
        # p1.<sample>.Aligned.out.bam    pass 1
        # <sample>.Aligned.out.bam       pass 2
        # Make sure to ignore the pass1 files.
        x = filelib.list_files_in_path(
            in_data.identifier, endswith=".Aligned.out.bam",
            file_not_startswith="p1.")
        bam_filenames = x
        if not bam_filenames:
            x = filelib.list_files_in_path(
                in_data.identifier, endswith=".Aligned.out.sam")
            sam_filenames = x
            if sam_filenames:
                assert bam_filenames, \
                       "No .Aligned.out.bam files.  Looks like .sam generated."
            assert bam_filenames, "No .Aligned.out.bam files."
        filelib.safe_mkdir(out_path)

        jobs = []  # list of (in_filename, out_filename)
        for in_filename in bam_filenames:
            # in_filename has format:
            # <path>/<sample>.Aligned.out.sam
            path, f = os.path.split(in_filename)
            sample, x = f.split(".", 1)
            assert x == "Aligned.out.bam", f
            out_filename = os.path.join(out_path, "%s.bam" % sample)
            assert in_filename != out_filename
            jobs.append((in_filename, out_filename))

        # Make sure outfiles are unique.
        x = [x[-1] for x in jobs]
        x = {}.fromkeys(x)
        assert len(jobs) == len(x), "Duplicate sample names."

        for x in jobs:
            in_filename, out_filename = x
            os.symlink(in_filename, out_filename)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)


    def name_outfile(self, antecedents, user_options):
        return "star.bam"

