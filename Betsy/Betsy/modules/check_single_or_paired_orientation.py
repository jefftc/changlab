from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import shutil
        
        group_node, fastq_node, reference_node = antecedents
        shutil.copy2(group_node.identifier, outfile)
        
        
    def set_out_attributes(self, antecedents, out_attributes):
        import os
        #from genomicode import config
        #from genomicode import filelib
        from genomicode import alignlib
        from Betsy import module_utils

        group_node, fastq_node, reference_node = antecedents
        sample_group_file = group_node.identifier
        fastq_path = fastq_node.identifier
        assert os.path.exists(fastq_path)
        assert os.path.isdir(fastq_path)
        ref = alignlib.create_reference_genome(reference_node.identifier)

        # FASTQ files are merged.  Find the merged files.
        x = module_utils.find_merged_fastq_files(sample_group_file, fastq_path)
        fastq_files = x
        assert fastq_files, "No fastq files."

        # Possibilities:
        # 1.  All single.
        # 2.  All paired.
        # 3.  Mixed.
        attrs = out_attributes.copy()
        all_pair2 = [x[-1] for x in fastq_files]
        uniq_pair2 = {}.fromkeys(all_pair2).keys()
        if uniq_pair2 == [None]:
            # All single.
            attrs["orientation"] = "single"
            return attrs
        if None in all_pair2:
            # Mixed.
            raise AssertionError, "Mixed single and paired-end."
        # All paired.

        # Optimization: check just the first group of FASTQ files and
        # assume they all have the same orientation.
        sample, pair1_filename, pair2_filename = fastq_files[0]
        x = get_paired_orientation(
            ref.fasta_file_full, pair1_filename, pair2_filename)
        orientation = "paired"
        if x:
            orientation = "paired_%s" % x
        attrs["orientation"] = orientation
        return attrs
    

    def name_outfile(self, antecedents, user_options):
        # Maintain the same extension so I know what file type this is.
        import os
        group_node, fastq_node, reference_node = antecedents
        f, ext = os.path.splitext(group_node.identifier)
        return "sample_group%s" % ext


def copy_fastq(in_filename, out_filename, MAX_READS=None):
    from genomicode import genomelib

    in_iter = genomelib.read_fastq(in_filename)
    out_handle = open(out_filename, 'w')
    i = 0
    while MAX_READS is None or i < MAX_READS:
        i += 1
        x = in_iter.next()
        if not x:  # no more reads
            assert i  # make sure at least 1 read
            break
        genomelib.write_fastq(*x, **{"handle" : out_handle})


def get_paired_orientation(
    reference_genome, filename1, filename2, outpath=None):
    # Return None, "ff", "fr", or "rf".  None means unstranded.
    import os
    import shutil
    import tempfile
    from genomicode import genomelib
    from genomicode import alignlib
    from genomicode import shell

    # Strategy: run bowtie2 in all orientations.  Return the one with
    # most reads aligned.  Do with a subset of the data, so this
    # doesn't take a long time.

    NUM_READS = 100

    # If outpath is None, then put everything into a temporary
    # directory.
    path = outpath   # where to write the results
    tempdir = None   # temporary directory to be deleted
    try:
        if path is None:
            tempdir = tempfile.mkdtemp(dir=".")
            path = tempdir   # write into a temporary directory

        #short_filename1 = os.path.join(path, "short_1.fq")
        #short_filename2 = os.path.join(path, "short_2.fq")
        #copy_fastq(filename1, short_filename1, NUM_READS)
        #copy_fastq(filename2, short_filename2, NUM_READS)

        sam_ff = os.path.join(path, "orient_ff.sam")
        sam_fr = os.path.join(path, "orient_fr.sam")
        sam_rf = os.path.join(path, "orient_rf.sam")
        log_ff = os.path.join(path, "orient_ff.log")
        log_fr = os.path.join(path, "orient_fr.log")
        log_rf = os.path.join(path, "orient_rf.log")

        x1 = alignlib.make_bowtie2_command(
            reference_genome, fastq_file1=filename1,
            fastq_file2=filename2, sam_file=sam_ff, orientation="ff",
            max_reads=NUM_READS)
        x2 = alignlib.make_bowtie2_command(
            reference_genome, fastq_file1=filename1,
            fastq_file2=filename2, sam_file=sam_fr, orientation="fr",
            max_reads=NUM_READS)
        x3 = alignlib.make_bowtie2_command(
            reference_genome, fastq_file1=filename1,
            fastq_file2=filename2, sam_file=sam_rf, orientation="rf",
            max_reads=NUM_READS)
        x1 += " >& %s" % log_ff
        x2 += " >& %s" % log_fr
        x3 += " >& %s" % log_rf
        commands = [x1, x2, x3]

        shell.parallel(commands)

        # Read the results.
        output_ff = alignlib.parse_bowtie2_output(log_ff)
        output_fr = alignlib.parse_bowtie2_output(log_fr)
        output_rf = alignlib.parse_bowtie2_output(log_rf)

        reads_ff = output_ff["concordant_reads"]
        reads_fr = output_fr["concordant_reads"]
        reads_rf = output_rf["concordant_reads"]

        orient = [
            (reads_ff, "ff"),
            (reads_fr, "fr"),
            (reads_rf, "rf"),
            ]
        orient.sort()

        # If all are within 10%, then it's probably un-stranded.
        cutoff = int(output_ff["reads_processed"] * 0.10)
        if orient[2][0] - orient[0][0] <= cutoff:
            return None
        return orient[2][-1]
    finally:
        if tempdir is not None and os.path.exists(tempdir):
            shutil.rmtree(tempdir)

