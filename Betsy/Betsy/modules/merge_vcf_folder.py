from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        import os
        import shutil
        from genomicode import filelib
        from genomicode import shell
        from Betsy import module_utils

        vcf_node = in_data
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf")
        assert vcf_filenames, "No .vcf files."

        bgzip = module_utils.get_config("bgzip")
        tabix = module_utils.get_config("tabix")
        bcftools = module_utils.get_config("bcftools")
        sq = shell.quote
        
        tmp_path = "indexed.vcf"
        tmp_path = os.path.realpath(tmp_path)
        filelib.safe_mkdir(tmp_path)
        

        # list of (in_filename, tmp_filename)
        jobs = []
        for in_filename in vcf_filenames:
            p, f = os.path.split(in_filename)
            #sample, ext = os.path.splitext(f)
            tmp_filename = os.path.join(tmp_path, f)
            x = in_filename, tmp_filename
            jobs.append(x)

        # Copy all the VCF files to a temporary directory.
        for x in jobs:
            in_filename, tmp_filename = x
            shutil.copy2(in_filename, tmp_filename)

        # Compress the VCF files.
        # bgzip file.vcf
        commands = []
        for x in jobs:
            in_filename, tmp_filename = x
            x = "%s %s" % (sq(bgzip), sq(tmp_filename))
            commands.append(x)
        shell.parallel(commands, max_procs=num_cores, path=tmp_path)
        x = ["%s.gz" % x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)

        # Index the VCF files.
        # tabix -p vcf file.vcf
        commands = []
        for x in jobs:
            in_filename, tmp_filename = x
            x = "%s -p vcf %s.gz" % (sq(tabix), sq(tmp_filename))
            commands.append(x)
        shell.parallel(commands, max_procs=num_cores, path=tmp_path)
        x = ["%s.gz.tbi" % x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)
        
        # Run bcftools
        cmd = [
            sq(bcftools),
            "merge", 
            "-o %s" % sq(out_filename),
            "-O v",
            ]
        for x in jobs:
            in_filename, tmp_filename = x
            cmd.append("%s.gz" % tmp_filename)
        x = " ".join(cmd)
        shell.single(x)
        filelib.assert_exists_nz(out_filename)


    def name_outfile(self, antecedents, user_options):
        return "merged.vcf"
