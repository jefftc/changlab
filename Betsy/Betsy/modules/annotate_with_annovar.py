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
        from Betsy import module_utils

        vcf_node = in_data
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf")
        assert vcf_filenames, "No .vcf files."
        filelib.safe_mkdir(out_path)

        buildver = module_utils.get_user_option(
            user_options, "buildver", allowed_values=["hg19"], not_empty=True)

        jobs = []  # list of (in_filename, log_filename, out_filestem)
        for in_filename in vcf_filenames:
            # Annovar takes a filestem, without the ".vcf".
            p, f = os.path.split(in_filename)
            f, exp = os.path.splitext(f)
            log_filename = os.path.join(out_path, "%s.log" % f)
            out_filestem = os.path.join(out_path, f)
            x = in_filename, log_filename, out_filestem
            jobs.append(x)
            
        # Make a list of commands.
        commands = []
        for x in jobs:
            in_filename, log_filename, out_filestem = x

            x = alignlib.make_annovar_command(
                in_filename, log_filename, out_filestem, buildver)
            commands.append(x)
            
        #for x in commands:
        #    print x
        #import sys; sys.exit(0)

        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x = [x[-1] for x in jobs]  # out_filestems
        x = ["%s.%s_multianno.vcf" % (x, buildver) for x in x]
        filelib.assert_exists_nz_many(x)


    def name_outfile(self, antecedents, user_options):
        return "annovar.vcf"


