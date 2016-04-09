from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from genomicode import config
        from Betsy import module_utils as mlib

        mpileup_node = in_data
        mpileup_filenames = filelib.list_files_in_path(
            mpileup_node.identifier, endswith=".pileup")
        assert mpileup_filenames, "No .pileup files."
        nc_match = mlib.read_normal_cancer_file(nc_node.identifier)
        #ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        # Figure out whether the purpose is to get coverage.  Change
        # the parameters if it is.
        assert "vartype" in out_attributes
        vartype = out_attributes["vartype"]
        assert vartype in ["snp", "indel"]
        tool = "mpileup2snp"
        if vartype == "indel":
            tool = "mpileup2indel"


        # list of (sample, in_filename, tmp1_filename, tmp2_filename,
        #          out_filename)
        jobs = []
        for in_filename in mpileup_filenames:
            p, sample, ext = mlib.splitpath(in_filename)
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
        x = [x[2] for x in jobs]
        filelib.assert_exists_nz_many(x)
        

        # java -jar /usr/local/bin/VarScan.jar <tool> $i --output_vcf 1 > $j
        varscan = filelib.which_assert(config.varscan_jar)
        
        # Make a list of commands.
        commands = []
        for x in jobs:
            sample, in_filename, tmp1_filename, tmp2_filename, out_filename = x
            x = [
                "java", "-jar", sq(varscan),
                tool,
                tmp1_filename,
                "--output-vcf", 1,
                ]
            x = " ".join(map(str, x))
            x = "%s >& %s" % (x, tmp2_filename)
            commands.append(x)

        #for x in commands:
        #    print x
        #import sys; sys.exit(0)

        parallel.pshell(commands, max_procs=num_cores)
        x = [x[3] for x in jobs]
        filelib.assert_exists_nz_many(x)

        # Clean up the VCF files.  VarScan leaves extraneous lines
        # there.
        for x in jobs:
            sample, in_filename, tmp1_filename, tmp2_filename, out_filename = x
            alignlib.clean_varscan_vcf(sample, tmp2_filename, out_filename)
        x = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)


    def name_outfile(self, antecedents, user_options):
        return "varscan.vcf"


