from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        import stat
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from genomicode import config

        pileup_filenames = filelib.list_files_in_path(
            in_data.identifier, endswith=".pileup")
        assert pileup_filenames, "No .pileup files."
        #ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        # TODO: get version of varscan

        # Figure out whether the purpose is to get coverage.  Change
        # the parameters if it is.
        assert "vartype" in out_attributes
        vartype = out_attributes["vartype"]
        assert vartype in ["all", "snp", "indel", "consensus"]

        #if vartype == "consensus":
        #    # Figure out the consensus-specific arguments.
        #    pass
        #else:
        #    raise NotImplementedError

        # list of (sample, in_filename, tmp1_filename, tmp2_filename,
        #          out_filename)
        jobs = []
        for in_filename in pileup_filenames:
            p, f = os.path.split(in_filename)
            sample, ext = os.path.splitext(f)
            tmp1_filename = os.path.join(out_path, "%s.tmp1" % sample)
            tmp2_filename = os.path.join(out_path, "%s.tmp2" % sample)
            out_filename = os.path.join(out_path, "%s.vcf" % sample)
            x = sample, in_filename, tmp1_filename, tmp2_filename, out_filename
            jobs.append(x)

        # VarScan will generate a "Parsing Exception" if there are 0
        # reads in a location.  Filter those out.
        sq = parallel.quote
        commands = []
        for x in jobs:
            sample, in_filename, tmp1_filename, tmp2_filename, out_filename = x
            x = "awk -F'\t' '$4 != 0 {print}' %s > %s" % (
                in_filename, tmp1_filename)
            commands.append(x)
        parallel.pshell(commands, max_procs=num_cores)
        # Some files may be empty if there are no reads.
        x = [x[2] for x in jobs]
        filelib.assert_exists_many(x)

        # java -jar /usr/local/bin/VarScan.jar mpileup2cns $i \
        #   --min-coverage 0 --min-reads2 0 --min-avg-qual 0 --min-var-freq 0 \
        #   --p-value 1.0 --strand-filter 0 --output-vcf 1 > $j
        varscan = filelib.which_assert(config.varscan_jar)

        # If there are no reads in the pileup file, then just generate
        # an empty file and don't call varscan.
        i = 0
        while i < len(jobs):
            x = jobs[i]
            sample, in_filename, tmp1_filename, tmp2_filename, out_filename = x
            if os.stat(in_filename)[stat.ST_SIZE] == 0:
                open(out_filename, 'w')
                del jobs[i]
            else:
                i += 1
        
        # Make a list of commands.
        commands = []
        for x in jobs:
            sample, in_filename, tmp1_filename, tmp2_filename, out_filename = x
            x = [
                "java", "-jar", sq(varscan),
                "mpileup2cns",
                tmp1_filename,
                "--min-coverage", 0,
                "--min-reads2", 0,
                "--min-avg-qual", 0,
                "--min-var-freq", 0,
                "--p-value", 1.0,
                "--strand-filter", 0,
                "--output-vcf", 1,
                ]
            x = " ".join(map(str, x))
            x = "%s >& %s" % (x, tmp2_filename)
            commands.append(x)

        parallel.pshell(commands, max_procs=num_cores)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores

        x = [x[3] for x in jobs]
        filelib.assert_exists_nz_many(x)

        # Clean up the VCF files.  VarScan leaves extraneous lines
        # there.
        for x in jobs:
            sample, in_filename, tmp1_filename, tmp2_filename, out_filename = x
            alignlib.clean_varscan_vcf(sample, tmp2_filename, out_filename)

        # Files in out_path can get very big.  Clean them up.
        for x in jobs:
            sample, in_filename, tmp1_filename, tmp2_filename, out_filename = x
            if os.path.exists(tmp1_filename):
                os.unlink(tmp1_filename)
            if os.path.exists(tmp2_filename):
                os.unlink(tmp2_filename)

        x = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "varscan.vcf"
