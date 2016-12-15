from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils as mlib

        pindel_node, ref_node = antecedents
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        # Should come up with better defaults for these.
        REF_NAME = "reference"
        REF_DATE = "date"

        # pindel2vcf -p <pindel_output_file> -r <reference_file>
        #   -R <name_and_version_of_reference_genome>
        #   -d <date_of_reference_genome_version>
        #   [-v <vcf_output_file>]
        #
        # -P <out_prefix>   Instead of -p.

        # Find the <out_prefix> in the pindel path.  <out_prefix>
        # should be the sample names.  Identify them by looking for
        # the config files:
        #   <sample>.config.txt
        pindel_path = pindel_node.identifier
        x = os.listdir(pindel_path)
        x = [x for x in x if x.endswith(".config.txt")]
        assert x, "I could not find config files: %s" % pindel_path
        x = [x.replace(".config.txt", "") for x in x]
        samples = x

        jobs = []
        for sample in samples:
            out_prefix = os.path.join(pindel_path, sample)
            # Make sure out_prefix is right.
            filelib.assert_exists("%s_D" % out_prefix)
            filelib.assert_exists("%s_SI" % out_prefix)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            out_filename = os.path.join(out_path, "%s.vcf" % sample)
            x = filelib.GenericObject(
                sample=sample,
                out_prefix=out_prefix,
                log_filename=log_filename,
                out_filename=out_filename)
            jobs.append(x)

        # Make a list of commands.
        pindel2vcf = mlib.get_config("pindel2vcf", which_assert_file=True)
        sq = parallel.quote
        commands = []
        for j in jobs:
            x = [
                sq(pindel2vcf),
                "-P", sq(j.out_prefix),
                "-r", sq(ref.fasta_file_full),
                "-R", REF_NAME,
                "-d", REF_DATE,
                "-v", sq(j.out_filename),
                "--gatk_compatible",
                ]
            x = " ".join(map(str, x))
            x = "%s >& %s" % (x, j.log_filename)
            commands.append(x)
        parallel.pshell(commands, max_procs=num_cores)
        metadata["num_cores"] = num_cores
        metadata["commands"] = commands

        # Make sure the analysis completed successfully.
        x1 = [x.log_filename for x in jobs]
        x2 = [x.out_filename for x in jobs]
        filelib.assert_exists_nz_many(x1+x2)

        raise NotImplementedError

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "pindel.vcf"
    
