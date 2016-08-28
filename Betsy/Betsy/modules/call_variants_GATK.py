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
        from Betsy import module_utils

        bam_node, ref_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        # TODO: Figure out GATK version.

        ## Figure out whether the user wants SNPs or INDELs.
        #assert "vartype" in out_attributes
        #vartype = out_attributes["vartype"]
        #assert vartype in ["all", "snp", "indel"]

        jobs = []
        for bam_filename in bam_filenames:
            p, f = os.path.split(bam_filename)
            sample, ext = os.path.splitext(f)
            #raw_outfile = os.path.join(out_path, "%s.raw" % sample)
            vcf_outfile = os.path.join(out_path, "%s.vcf" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = filelib.GenericObject(
                bam_filename=bam_filename, vcf_outfile=vcf_outfile,
                log_filename=log_filename)
            jobs.append(x)
        
        # java -Xmx5g -jar /usr/local/bin/GATK/GenomeAnalysisTK.jar
        #   -T HaplotypeCaller -R ucsc.hg19.fasta
        #   -dontUseSoftClippedBases -stand_call_conf 20.0
        #   -stand_emit_conf 20.0 -I $i -o $j
               
        # Make a list of commands.
        commands = []
        for j in jobs:
            # For debugging.  If exists, don't do it again.
            #if filelib.exists_nz(j.raw_outfile):
            if filelib.exists_nz(j.vcf_outfile):
                continue
            x = alignlib.make_GATK_command(
                T="HaplotypeCaller", R=ref.fasta_file_full,
                dontUseSoftClippedBases=None, stand_call_conf=20.0,
                stand_emit_conf=20.0, I=j.bam_filename, o=j.vcf_outfile)
            x = "%s >& %s" % (x, j.log_filename)
            commands.append(x)

        parallel.pshell(commands, max_procs=num_cores)

        # Filter each of the VCF files.
        #for j in jobs:
        #    filter_by_vartype(vartype, j.raw_outfile, j.vcf_outfile)
        #metadata["filter"] = vartype

        # Make sure the analysis completed successfully.
        x = [j.vcf_outfile for j in jobs]
        filelib.assert_exists_nz_many(x)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "GATK.vcf"


