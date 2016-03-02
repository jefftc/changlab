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
        from genomicode import alignlib
        from Betsy import module_utils

        mvcf_node = in_data
        in_filename = mvcf_node.identifier
        filelib.assert_exists_nz(in_filename)

        buildver = module_utils.get_user_option(
            user_options, "buildver", allowed_values=["hg19"], not_empty=True)

        # Annovar takes a filestem, without the ".vcf".
        p, f = os.path.split(in_filename)
        f, exp = os.path.splitext(f)
        log_filename = "%s.log" % f
        
        p, f = os.path.split(out_filename)
        f, exp = os.path.splitext(f)
        out_filestem = f

        cmd = alignlib.make_annovar_command(
            in_filename, log_filename, out_filestem, buildver)
        parallel.sshell(cmd)

        # Make sure the analysis completed successfully.
        x = "%s.%s_multianno.vcf" % (out_filestem, buildver)
        filelib.assert_exists_nz(x)
        if os.path.realpath(x) != os.path.realpath(out_filename):
            shutil.copy2(x, out_filename)


    def name_outfile(self, antecedents, user_options):
        return "annotated.vcf"


