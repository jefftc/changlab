from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import shell
        from genomicode import alignlib
        from Betsy import module_utils

        align_node = in_data
        x = module_utils.find_bam_files(align_node.identifier)
        x = [x for x in x if x.endswith("accepted_hits.bam")]
        bam_filenames = x
        assert bam_filenames, "No accepted_hits.bam files."
        filelib.safe_mkdir(out_path)

        jobs = []  # list of (in_filename, out_filename)
        for in_filename in bam_filenames:
            # Names must in the format:
            # <path>/<sample>.tophat/accepted_hits.bam
            # full_path   <path>/<sample>.tophat
            # path        <path>
            # tophat_dir  <sample>.tophat
            # file_       accepted_hits.bam
            # sample      <sample>
            
            full_path, file_ = os.path.split(in_filename)
            path, tophat_dir = os.path.split(full_path)

            assert file_ == "accepted_hits.bam"
            assert tophat_dir.endswith(".tophat")
            sample = tophat_dir[:-7]
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
        return "tophat.bam"

