from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import shell
        from genomicode import alignlib
        from genomicode import config
        from Betsy import module_utils

        vcf_node = in_data
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf")
        assert vcf_filenames, "No .vcf files."
        filelib.safe_mkdir(out_path)

        jobs = []  # list of (in_filename, log_filename, out_filestem)
        for in_filename in vcf_filenames:
            # Annovar takes a filestem, without the ".vcf".
            p, f = os.path.split(in_filename)
            f, exp = os.path.splitext(f)
            log_filename = os.path.join(out_path, "%s.log" % f)
            out_filestem = os.path.join(out_path, f)
            x = in_filename, log_filename, out_filestem
            jobs.append(x)
            
        buildver = module_utils.get_user_option(
            user_options, "buildver", allowed_values=["hg19"], not_empty=True)

        # list of (name, operation).
        # These are just for buildver hg19.
        protocols = [
            ("refGene", "g"), ("cytoBand", "r"), ("genomicSuperDups", "r"),
            ("esp6500siv2_all", "f"), ("snp138", "f"),
            ("ljb26_all", "f"),
            ("1000g2015aug_all", "f"), ("1000g2015aug_afr", "f"),
            ("1000g2015aug_eas", "f"), ("1000g2015aug_eur", "f"),
            ]
            
        # P1=refGene,cytoBand,genomicSuperDups
        # P2=esp6500siv2_all,snp138,ljb26_all
        #P3=1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur
        # table_annovar.pl -buildver hg19 -remove -vcfinput
        #   -protocol $P1,$P2,$P3 -operation g,r,r,f,f,f,f,f,f,f
        #   -out $j $i humandb/

        table_annovar = filelib.which_assert(config.table_annovar)
        annodb = config.annovar_db
        assert os.path.exists(annodb)
        assert os.path.isdir(annodb)

        # Make a list of commands.
        sq = shell.quote
        commands = []
        for x in jobs:
            in_filename, log_filename, out_filestem = x

            x1 = [x[0] for x in protocols]
            x2 = [x[1] for x in protocols]

            x = [
                sq(table_annovar),
                "-buildver", buildver,
                "-remove",
                "-vcfinput",
                "-protocol", ",".join(x1),
                "-operation", ",".join(x2),
                "-out", sq(out_filestem),
                sq(in_filename),
                sq(annodb),
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, log_filename)
            commands.append(x)
            
        #for x in commands:
        #    print x
        #import sys; sys.exit(0)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x = [x[-1] for x in jobs]  # out_filestems
        x = ["%s.%s_multianno.vcf" % (x, buildver) for x in x]
        filelib.assert_exists_nz_many(x)


    def name_outfile(self, antecedents, user_options):
        return "annovar.vcf"

