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

        bam_node, ref_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        metadata["tool"] = "samtools %s" % alignlib.get_samtools_version()
        # bcfutils is same version as samtools.
        # Cannot get vcfutils version.

        # Ignore max_read_depth for now.
        #max_read_depth = mlib.get_user_option(
        #    user_options, "max_read_depth", type=int)
        #assert max_read_depth > 0 and max_read_depth <= 1E6

        # samtools mpileup -uf ref.fa aln1.bam aln2.bam |
        #   bcftools view -bvcg - > var.raw.bcf
        # bcftools view var.raw.bcf | vcfutils.pl varFilter -D100
        #   > var.flt.vcf
        # Actually, just do in one step.  Don't use vcfutils.pl to
        # filter.

        # list of (sample, bam_infilename, bcf_filename, vcf_outfilename)
        jobs = []
        for in_filename in bam_filenames:
            path, sample, ext = mlib.splitpath(in_filename)
            # bcf_filename not stored.
            bcf_filename = os.path.join(out_path, "%s.bcf" % sample)
            vcf_outfilename = os.path.join(out_path, "%s.vcf" % sample)
            x = sample, in_filename, bcf_filename, vcf_outfilename
            jobs.append(x)

        # Generate the bcf files.
        sq = mlib.sq
        samtools = mlib.findbin("samtools")
        bcftools = mlib.findbin("bcftools")
        #vcfutils = mlib.findbin("vcfutils")
        commands = []
        for x in jobs:
            sample, bam_infilename, bcf_filename, vcf_outfilename = x
            x = [
                sq(samtools), "mpileup", "-g", "-f", sq(ref.fasta_file_full),
                sq(bam_infilename),
                "|",
                #sq(bcftools), "view", "-bvcg", "-",
                sq(bcftools), "call", "-c", "-v",
                ">",
                sq(vcf_outfilename),
                ]
            x = " ".join(x)
            commands.append(x)
        parallel.pshell(commands, max_procs=num_cores)
        x = [x[3] for x in jobs]
        filelib.assert_exists_nz_many(x)

        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "mpileup.vcf"
