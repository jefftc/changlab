from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import parallel
        from genomicode import filelib
        from genomicode import alignlib
        from Betsy import module_utils

        bam_folder, sample_node = antecedents
        bam_path = bam_folder.identifier
        filelib.safe_mkdir(out_path)

        assert os.path.exists(bam_path)
        assert os.path.isdir(bam_path)

        attr2order = {
            "name" : "name",
            "coordinate" : "pos",
            }
        x = bam_folder.data.attributes["sorted"]
        sort_order = attr2order.get(x)
        assert sort_order, "Cannot handle sorted: %s" % x

        attr2stranded = {
            "single" : "no", 
            "paired" : "no",
            "paired_ff" : None,
            "paired_fr" : "yes",
            "paired_rf" : "reverse",
            }
        x = sample_node.data.attributes["orientation"]
        stranded = attr2stranded.get(x)
        assert stranded, "Cannot handle orientation: %s" % x

        gtf_file = module_utils.get_user_option(
            user_options, "gtf_file", not_empty=True)
        assert os.path.exists(gtf_file), "File not found: %s" % gtf_file

        mode = module_utils.get_user_option(
            user_options, "htseq_count_mode", allowed_values=
            ["union", "intersection-strict", "intersection-nonempty"])

        # Make a list of the jobs to run.
        jobs = []
        for bam_filename in filelib.list_files_in_path(
            bam_path, endswith=".bam", case_insensitive=True):
            x = os.path.split(bam_filename)[1]
            x = os.path.splitext(x)[0]
            x = "%s.count" % x
            out_file = x
            x = bam_filename, out_file
            jobs.append(x)
        
        # Generate commands for each of the files.
        commands = []
        for x in jobs:
            bam_filename, out_file = x
            x = alignlib.make_htseq_count_command(
                bam_filename, gtf_file, sort_order, stranded, mode=mode)
            x = "%s >& %s" % (x, out_file)
            commands.append(x)

        parallel.pshell(commands, max_procs=num_cores, path=out_path)

        # Make sure the analysis completed successfully.
        x = [x[1] for x in jobs]
        x = [os.path.join(out_path, x) for x in x]
        output_filenames = x
        filelib.assert_exists_nz_many(output_filenames)


    def name_outfile(self, antecedents, user_options):
        return "signal.counts"
