from Module import AbstractModule


class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import itertools
        from genomicode import config
        from genomicode import parallel
        from genomicode import filelib

        signal_node, annotation_node = antecedents
        signal_filename = signal_node.identifier
        annotation_filename = annotation_node.identifier
        filelib.assert_exists_nz(signal_filename)
        filelib.assert_exists_nz(annotation_filename)
        metadata = {}

        align_matrices = filelib.which_assert(config.align_matrices)

        # Make sure the signal_filename has an ID_REF header.
        header = filelib.read_cols(signal_filename).next()
        assert header[0] == "ID_REF", "Missing ID_REF header: %s" % \
               signal_filename

        signal_align_file = "signal.aligned.txt"
        annot_align_file = "annot.aligned.txt"

        # First, align the two files.
        sq = parallel.quote
        cmd = [
            sq(align_matrices),
            "--annot_file", signal_filename,
            "--header", "ID_REF",
            "--annot_file", annotation_filename,
            "--left_join",
            signal_align_file, annot_align_file,
            ]
        cmd = " ".join(cmd)
        parallel.sshell(cmd)
        metadata["command"] = cmd

        # Now merge them.  Take the first column of the expression
        # file (should be ID_REF), the whole annotation file, then the
        # remainder of the expression file.
        signal_handle = filelib.read_cols(signal_align_file)
        annot_handle = filelib.read_cols(annot_align_file)
        outhandle = open(outfile, 'w')
        for x1, x2 in itertools.izip(signal_handle, annot_handle):
            x = [x1[0]] + x2 + x1[1:]
            print >>outhandle, "\t".join(x)
        outhandle.close()

        #cmd = "paste %s %s > %s" % (
        #    annot_align_file, signal_align_file, outfile)
        #shell.single(cmd)

        filelib.assert_exists_nz(outfile)


    def name_outfile(self, antecedents, user_options):
        return "signal.txt"
