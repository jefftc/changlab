from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import parallel
        from genomicode import filelib
        from genomicode import genomelib
        from genomicode import config
        from Betsy import module_utils as mlib

        fasta_node, bam_node, sample_node, orient_node = antecedents
        fasta_data = mlib.find_merged_fastq_files(
            sample_node.identifier, fasta_node.identifier, find_fasta=True)
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        orient = mlib.read_orientation(orient_node.identifier)
        filelib.safe_mkdir(out_path)
        
        # TODO: Try to figure out version.
        metadata = {}
        metadata["tool"] = "RSeQC (unknown version)"

        pyrseqc = mlib.findbin("pyrseqc")
        
        gene_model = mlib.get_user_option(
            user_options, "gene_model", not_empty=True,
            allowed_values=["hg19"])
        if gene_model == "hg19":
            gene_path = config.rseqc_hg19
        else:
            raise AssertionError, "Unhandled: %s" % gene_model
        
        filelib.dir_exists(gene_path)
        gene_model_bed = os.path.join(gene_path, "RefSeq.bed12")
        housekeeping_model_bed = os.path.join(
            gene_path, "HouseKeepingGenes.bed")

        sample2fastadata = {}
        for x in fasta_data:
            sample, f1, f2 = x
            sample2fastadata[sample] = x

        is_paired = orient.orientation.startswith("paired")

        # Guess the read length.  Read the first fasta.
        assert sample2fastadata
        x = sample2fastadata.keys()[0]
        filename = sample2fastadata[x][1]
        lengths = {}  # length -> count
        for i, x in enumerate(genomelib.read_fasta_many(filename)):
            if i >= 100:
                break
            title, sequence = x
            l = len(sequence)
            lengths[l] = lengths.get(l, 0) + 1
        # Use the most common length.
        c_length = c_count = None
        for (l, c) in lengths.iteritems():
            if c_count is None or c > c_count:
                c_length, c_count = l, c
        assert c_length
        read_length = c_length

        jobs = []  # sample, bam_filename, fasta_file1, fasta_file2, outdir
        for bam_filename in bam_filenames:
            # <path>/<sample>.bam
            p, sample, e = mlib.splitpath(bam_filename)
            assert sample in sample2fastadata
            x, f1, f2 = sample2fastadata[sample]
            outdir = os.path.join(out_path, sample)
            x = sample, bam_filename, f1, f2, outdir
            jobs.append(x)

        # Some of the modules of RSeQC uses a lot of memory.  Have
        # seen a Python process take 33 Gb, and an R process take 200
        # Gb.  However, most of the modules use much less memory.  So
        # run one pyrseqc at a time, and run each one of those
        # processes in parallel.  Is probably slower than running
        # multiple pyrseqc, but takes less memory.
        commands = []
        for x in jobs:
            sample, bam_filename, fasta_filename1, fasta_filename2, outdir = x

            # pyrseqc.py -j 20 --paired_end rqc11.bam rqc14.fa 76 \
            #   mod07.txt hg19.HouseKeepingGenes.bed rqc21 --dry_run
            x = [
                mlib.sq(pyrseqc),
                "-j", str(num_cores),
                ]
            if is_paired:
                x += ["--paired_end"]
            x += [
                mlib.sq(bam_filename), 
                mlib.sq(fasta_filename1),
                str(read_length),
                mlib.sq(gene_model_bed),
                mlib.sq(housekeeping_model_bed),
                mlib.sq(outdir),
                ]
            x = " ".join(x)
            commands.append(x)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        # pyrseqc takes up to ~40 Gb per process.
        # read_distribution.py takes 33 Gb.
        # read_quality.py spins off an R process that takes ~200 Gb.
        # Make sure we don't use up more memory than is available on
        # the machine.
        #nc = mlib.calc_max_procs_from_ram(60, upper_max=num_cores)
        #metadata["num cores"] = nc
        #x = parallel.pshell(commands, max_procs=nc)

        # Because of memory, just run one at a time, but each one, use
        # multiple cores.
        for cmd in commands:
            x = parallel.sshell(cmd)
            assert x.find("Traceback") < 0, x
            
        filelib.assert_exists_nz(out_path)
        
        return metadata
        
    def name_outfile(self, antecedents, user_options):
        return "rseqc"
