from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import parallel
        from genomicode import alignlib
        from genomicode import filelib
        from Betsy import module_utils
        
        bam_node, ref_node, pos_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        positions_filename = pos_node.identifier
        filelib.assert_exists_nz(positions_filename)
        filelib.safe_mkdir(out_path)

        # list of (in_filename, err_filename, out_filename)
        jobs = []
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            sample, ext = os.path.splitext(f)
            err_filename = os.path.join(out_path, "%s.log" % sample)
            out_filename = os.path.join(out_path, "%s.vcf" % sample)
            x = in_filename, err_filename, out_filename
            jobs.append(x)

        ## Get possible positions file.
        #positions_filename = module_utils.get_user_option(
        #    user_options, "positions_file", check_file=True)
        
        # Figure out whether the purpose is to get coverage.  Change
        # the parameters if it is.
        assert "vartype" in out_attributes
        vartype = out_attributes["vartype"]
        assert vartype in ["all", "snp", "indel", "consensus"]
        #if cov == "yes":
        #    assert positions_filename, "Missing: positions_file"

        # samtools mpileup -l freq04.txt -R -B -q 0 -Q 0 -d10000000 \
        #   -f genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
        #   $i > $j"
        samtools = filelib.which_assert(config.samtools)

        if vartype == "consensus":
            args = [
                "-R",        # Ignore read group tags.
                "-B",        # Disable BAQ (base quality) computation.
                "-q", 0,     # Skip bases with mapQ smaller than this.
                "-Q", 0,     # Skip bases with BAQ smaller than this.
                "-d10000000",  # Allow deep reads.
                ]
        else:
            raise NotImplementedError

        sq = parallel.quote
        commands = []
        for x in jobs:
            in_filename, err_filename, out_filename = x
            
            x = [
                sq(samtools),
                "mpileup",
                "-f", sq(ref.fasta_file_full),
                ]
            if positions_filename:
                x.extend(["-l", positions_filename])
            x.extend(args)
            x.append(sq(in_filename))
            x = " ".join(map(str, x))
            x = "%s 2> %s 1> %s" % (x, err_filename, out_filename)
            commands.append(x)

        #for x in commands:
        #    print x
            
        parallel.pshell(commands, max_procs=num_cores)

        x = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)
        

    def name_outfile(self, antecedents, user_options):
        return "consensus.mpileup"
