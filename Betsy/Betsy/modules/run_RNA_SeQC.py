from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import filelib
        from genomicode import shell
        from genomicode import alignlib
        from genomicode import hashlib
        from Betsy import module_utils

        bam_node, ref_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        # java -jar /usr/local/bin/RNA-SeQC_v1.1.8.jar \
        #   -o <sample> -r <reference_file> -s "<sample>|<in_filename>|NA"
        #   -t <gtf_file> >& <log_filename>"
        # <out_path>        Output directory.  Will be created if not exists.
        # <in_filename>     BAM file
        # <reference_file>  /data/biocore/genomes/UCSC/mm10.fa
        # <gtf_file>   /data/biocore/rsem/mouse_refseq_mm10/UCSC_knownGenes.gtf
        #
        # <reference_file> must be indexed and have a dict file.

        rna_seqc_jar = filelib.which_assert(config.rna_seqc_jar)

        GTF = module_utils.get_user_option(
            user_options, "rna_seqc_gtf_file", not_empty=True)
        assert os.path.exists(GTF), "File not found: %s" % GTF

        # list of infile, out_path, ref_file, gtf_file, sample, log_file
        jobs = []
        for in_filename in bam_filenames:
            p, file_ = os.path.split(in_filename)
            f, e = os.path.splitext(file_)
            sample = hashlib.hash_var(f)
            out_path_rna_seqc = os.path.join(out_path, sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)

            x = in_filename, out_path_rna_seqc, ref.fasta_file_full, GTF, \
                sample, log_filename
            jobs.append(x)

        sq = shell.quote
        commands = []
        for x in jobs:
            (in_filename, out_path_rna_seqc, ref_filename, gtf_filename, \
             sample, log_filename) = x

            x = [sample, in_filename, "NA"]
            x = "|".join(x)
            x = [
                'java',
                '-jar', rna_seqc_jar,
                '-o', sq(out_path_rna_seqc),
                '-r', sq(ref_filename),
                '-s', "'%s'" % x,
                '-t', gtf_filename,
                ]
            x = " ".join(x)
            cmd = "%s >& %s" % (x, log_filename)
            commands.append(cmd)

        # XXX
        for cmd in commands:
            print cmd
        import sys; sys.exit(0)
        
        x = shell.parallel(commands, max_procs=num_cores)
        run_log = os.path.join(out_path, "run.log")
        open(run_log, 'w').write(x)

        # Check for outfile.
        # Make sure the analysis completed successfully.
        for x in jobs:
            (in_filename, out_path_rna_seqc, ref_filename, gtf_filename, \
             sample, log_filename) = x
            filelib.assert_exists_nz(out_path_rna_seqc)
    
    
    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'bamFolder_' + original_file
        #return filename
        return "rna_seqc"
