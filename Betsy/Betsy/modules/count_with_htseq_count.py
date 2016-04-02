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
        from Betsy import module_utils as mlib

        bam_folder, sample_node, gene_node, strand_node = antecedents
        bam_path = bam_folder.identifier
        assert mlib.dir_exists(bam_path)
        gtf_file = gene_node.identifier
        filelib.assert_exists_nz(gtf_file)
        stranded = mlib.read_stranded(strand_node.identifier)
        filelib.safe_mkdir(out_path)
        
        metadata = {}

        attr2order = {
            "name" : "name",
            "coordinate" : "pos",
            }
        x = bam_folder.data.attributes["sorted"]
        sort_order = attr2order.get(x)
        assert sort_order, "Cannot handle sorted: %s" % x

        #attr2stranded = {
        #    "single" : "no", 
        #    "paired" : "no",
        #    "paired_ff" : None,
        #    "paired_fr" : "yes",
        #    "paired_rf" : "reverse",
        #    }
        #x = sample_node.data.attributes["orientation"]
        #stranded = attr2stranded.get(x)
        #assert stranded, "Cannot handle orientation: %s" % x

        ht_stranded = None
        if stranded.stranded == "unstranded":
            ht_stranded = "no"
        elif stranded.stranded == "firststrand":
            ht_stranded = "reverse"
        elif stranded.stranded == "secondstrand":
            ht_stranded = "yes"
        assert ht_stranded is not None


        #gtf_file = mlib.get_user_option(
        #    user_options, "gtf_file", not_empty=True)
        #assert os.path.exists(gtf_file), "File not found: %s" % gtf_file

        mode = mlib.get_user_option(
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
                bam_filename, gtf_file, sort_order, ht_stranded, mode=mode)
            x = "%s >& %s" % (x, out_file)
            commands.append(x)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        parallel.pshell(commands, max_procs=num_cores, path=out_path)

        # Make sure the analysis completed successfully.
        x = [x[1] for x in jobs]
        x = [os.path.join(out_path, x) for x in x]
        output_filenames = x
        filelib.assert_exists_nz_many(output_filenames)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "signal.counts"
