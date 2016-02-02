from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib

        x = filelib.list_files_in_path(
            in_data.identifier, endswith=".Aligned.out.sam")
        sam_filenames = x
        assert sam_filenames, "No .Aligned.out.sam files."
        filelib.safe_mkdir(out_path)

        jobs = []  # list of (in_filename, out_filename)
        for in_filename in sam_filenames:
            # in_filename has format:
            # <path>/<sample>.Aligned.out.sam
            path, x = os.path.split(in_filename)
            sample, x = x.split(".", 1)
            assert x == "Aligned.out.sam"
            out_filename = os.path.join(out_path, "%s.sam" % sample)
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
        return "star.sam"

