from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import shell
        from genomicode import alignlib
        from Betsy import module_utils

        bam_node, ref_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        jobs = []  # list of (in_filename, log_filename, out_filename)
        for in_filename in bam_filenames:
            p, f = os.path.split(in_filename)
            sample, ext = os.path.splitext(f)
            out_filename = os.path.join(out_path, "%s.vcf" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = in_filename, log_filename, out_filename
            jobs.append(x)
        
        # java -Xmx5g -jar /usr/local/bin/GATK/GenomeAnalysisTK.jar
        #   -T HaplotypeCaller -R ucsc.hg19.fasta
        #   -dontUseSoftClippedBases -stand_call_conf 20.0
        #   -stand_emit_conf 20.0 -I $i -o $j
               
        # Make a list of commands.
        commands = []
        for x in jobs:
            in_filename, log_filename, out_filename = x
            x = alignlib.make_GATK_command(
                T="HaplotypeCaller", R=ref.fasta_file_full,
                dontUseSoftClippedBases=None, stand_call_conf=20.0,
                stand_emit_conf=20.0, I=in_filename, o=out_filename)
            x = "%s >& %s" % (x, log_filename)
            commands.append(x)

        shell.parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)


    def name_outfile(self, antecedents, user_options):
        return "GATK.vcf"

