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
        from genomicode import parallel
        from genomicode import config

        in_filename = in_data.identifier
        filelib.assert_exists_nz(in_filename)

        vcftools = filelib.which_assert(config.vcftools)
        
        # vcftools --vcf test31.txt --remove-indels --recode --recode-INFO-all
        #   --out test32
        # Writes stuff to console.  Should capture in log file.
        # Saves file test32.recode.vcf

        p, f = os.path.split(in_filename)
        s, ext = os.path.splitext(in_filename)
        sample = s

        out_stem = "%s.filtered" % sample
        log_filename = "%s.log" % sample
        # Should create file <out_stem>.recode.vcf
        outfile = "%s.recode.vcf" % out_stem

        sq = parallel.quote
        cmd = [
            sq(vcftools),
            "--vcf", sq(in_filename),
            "--remove-indels",
            "--recode",
            "--recode-INFO-all",
            "--out", out_stem,
            ]
        cmd = " ".join(cmd)
        cmd = "%s >& %s" % (cmd, log_filename)
        parallel.sshell(cmd)

        filelib.assert_exists_nz(outfile)
        shutil.copy2(outfile, out_filename)


    def name_outfile(self, antecedents, user_options):
        return "filtered.vcf"
