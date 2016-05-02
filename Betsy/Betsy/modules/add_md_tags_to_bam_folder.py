from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils

        ## Importing pysam is hard!
        #import sys
        #sys_path_old = sys.path[:]
        #sys.path = [x for x in sys.path if x.find("RSeQC") < 0]
        #import pysam
        #sys.path = sys_path_old

        bam_node, ref_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        # list of (in_filename, err_filename, out_filename)
        jobs = []
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            s, ext = os.path.splitext(f)
            log_filename = os.path.join(out_path, "%s.log" % s)
            out_filename = os.path.join(out_path, f)
            assert in_filename != out_filename
            x = in_filename, log_filename, out_filename
            jobs.append(x)

        # Don't do this.  Need MD, NM, NH in
        # summarize_alignment_cigar.  To be sure, just redo it.
        ## If the files already have MD tags, then just symlink the
        ## files.  Don't add again.
        #i = 0
        #while i < len(jobs):
        #    in_filename, out_filename = jobs[i]
        #
        #    handle = pysam.AlignmentFile(in_filename, "rb")
        #    align = handle.next()
        #    tag_dict = dict(align.tags)
        #    if "MD" not in tag_dict:
        #        i += 1
        #        continue
        #    # Has MD tags.  Just symlink and continue.
        #    os.symlink(in_filename, out_filename)
        #    del jobs[i]

        # Make a list of samtools commands.
        # Takes ~200 Mb per process, so should not be a big issue.
        samtools = filelib.which_assert(config.samtools)
        sq = parallel.quote
        commands = []
        for x in jobs:
            in_filename, log_filename, out_filename = x

            # samtools calmd -b <in.bam> <ref.fasta> > <out.bam>

            # May generate error:
            # [bam_fillmd1] different NM for read
            #   'ST-J00106:118:H75L3BBXX:3:2128:21846:47014': 0 -> 19
            # Pipe stderr to different file.
            x = [
                samtools,
                "calmd", "-b",
                sq(in_filename),
                sq(ref.fasta_file_full),
                ]
            x = " ".join(x)
            x = "%s 2> %s 1> %s" % (x, sq(log_filename), sq(out_filename))
            commands.append(x)
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)

    
    def name_outfile(self, antecedents, user_options):
        return "md.bam"
