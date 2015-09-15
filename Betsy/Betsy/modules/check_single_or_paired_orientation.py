from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import shutil
        group_node, fastq_node, reference_node = antecedents

        shutil.copyfile(group_node.identifier, outfile)
        
        
    def set_out_attributes(self, antecedents, out_attributes):
        import os
        #from genomicode import config
        #from genomicode import filelib
        from Betsy import module_utils

        attrs = out_attributes.copy()
        group_node, fastq_node, reference_node = antecedents

        sample_group_file = group_node.identifier
        fastq_path = fastq_node.identifier
        reference_path = reference_node.identifier

        assert os.path.exists(fastq_path)
        assert os.path.exists(reference_path)
        assert os.path.isdir(fastq_path)
        assert os.path.isdir(reference_path)

        ## module_utils.assert_sample_group_file(sample_group_file, fastq_path)
        ## x = module_utils.read_sample_group_file(sample_group_file)
        ## x = module_utils.fix_sample_group_filenames(x, fastq_path)
        ## sample_groups = x

        ## x = [x[2] for x in sample_groups]
        ## x = {}.fromkeys(x)
        ## pairs = sorted(x)
        ## for x in pairs:
        ##     assert x in [None, 1, 2], "Invalid pair: %s" % repr(x)
        ## # If pair is all None, then orientation is single.
        ## if len(pairs) == 1:
        ##     pair = pairs[0]
        ##     if pair in [None, 1]:
        ##         attrs["orientation"] = "single"
        ##         return attrs
        ##     raise AssertionError, "Weird pairing."
        ## assert len(pairs) == 2, "Invalid pairs."
        ## if pairs == [None, 1]:
        ##     # A big weird, but acceptable.
        ##     attrs["orientation"] = "single"
        ##     return attrs
        ## if pairs == [None, 2]:
        ##     raise AssertionError, "Weird pairing."
        ## if pairs == [None, 1, 2]:
        ##     raise AssertionError, "Multiple pairings."
        ## assert pairs == [1, 2]
        
        # Paired ends.  Need to figure out which orientation.
        reference_genome = module_utils.find_bowtie2_reference(reference_path)

        # FASTQ files are merged.  Find the merged files.
        x = module_utils.find_merged_fastq_files(sample_group_file, fastq_path)
        fastq_files = x
        assert fastq_files, "No fastq files."

        # Possibilities:
        # 1.  All single.
        # 2.  All paired.
        # 3.  Mixed.
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
        x = _get_paired_orientation(
            reference_genome, pair1_filename, pair2_filename)
        orientation = "paired_%s" % x
        attrs["orientation"] = orientation
        return attrs
    

    def name_outfile(self, antecedents, user_options):
        # Maintain the same extension so I know what file type this is.
        import os
        group_node, fastq_node, reference_node = antecedents
        f, ext = os.path.splitext(group_node.identifier)
        return "sample_group%s" % ext




def _get_paired_orientation(reference_genome, filename1, filename2):
    # Return "ff", "fr", or "rf".
    import os
    import shutil
    import tempfile
    from genomicode import genomelib
    from Betsy import module_utils

    # Strategy: run bowtie2 in all orientations.  Return the one with
    # most reads aligned.  Do with a subset of the data, so this
    # doesn't take a long time.

    NUM_READS = 1000

    tempdir = None
    try:
        tempdir = tempfile.mkdtemp(dir=".")

        short_filename1 = os.path.join(tempdir, "short_1.fq")
        short_filename2 = os.path.join(tempdir, "short_2.fq")

        in_iter1 = genomelib.read_fastq(filename1)
        in_iter2 = genomelib.read_fastq(filename2)
        out_handle1 = open(short_filename1, 'w')
        out_handle2 = open(short_filename2, 'w')
        for i in range(NUM_READS):
            # How does this work?  StopIteration?
            x1 = in_iter1.next()
            x2 = in_iter2.next()
            if not x1 or not x2:  # doesn't have NUM_READS
                assert i  # make sure at least 1 read
                break
            genomelib.write_fastq(*x1, **{"handle" : out_handle1})
            genomelib.write_fastq(*x2, **{"handle" : out_handle2})
        out_handle1.close()
        out_handle2.close()

        sam_ff = os.path.join(tempdir, "orient_ff.sam")
        sam_fr = os.path.join(tempdir, "orient_fr.sam")
        sam_rf = os.path.join(tempdir, "orient_rf.sam")
        log_ff = os.path.join(tempdir, "orient_ff.log")
        log_fr = os.path.join(tempdir, "orient_fr.log")
        log_rf = os.path.join(tempdir, "orient_rf.log")

        x1 = module_utils.make_bowtie2_command(
            reference_genome, fastq_file1=short_filename1,
            fastq_file2=short_filename2, sam_file=sam_ff, orientation="ff")
        x2 = module_utils.make_bowtie2_command(
            reference_genome, fastq_file1=short_filename1,
            fastq_file2=short_filename2, sam_file=sam_fr, orientation="fr")
        x3 = module_utils.make_bowtie2_command(
            reference_genome, fastq_file1=short_filename1,
            fastq_file2=short_filename2, sam_file=sam_rf, orientation="rf")
        x1 += " >& %s" % log_ff
        x2 += " >& %s" % log_fr
        x3 += " >& %s" % log_rf
        commands = [x1, x2, x3]
        
        module_utils.run_parallel(commands)

        # Read the results.
        output_ff = module_utils.parse_bowtie2_output(log_ff)
        output_fr = module_utils.parse_bowtie2_output(log_fr)
        output_rf = module_utils.parse_bowtie2_output(log_rf)

        reads_ff = output_ff["concordant_reads"]
        reads_fr = output_fr["concordant_reads"]
        reads_rf = output_rf["concordant_reads"]

        reads, orientation = reads_ff, "ff"
        if reads_fr > reads:
            reads, orientation = reads_fr, "fr"
        if reads_rf > reads:
            reads, orientation = reads_rf, "rf"
        return orientation
    finally:
        if tempdir is not None and os.path.exists(tempdir):
            shutil.rmtree(tempdir)


