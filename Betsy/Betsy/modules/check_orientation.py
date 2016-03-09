from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)
        self.cache = {}

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import shutil
        
        fastq_node, group_node, reference_node = antecedents
        metadata = {}
        
        key = (group_node.identifier, fastq_node.identifier,
               reference_node.identifier)
        x = self.cache.get(key)
        if x:
            reads_ns, reads_fr, reads_rf, reads_ff = x
            metadata["reads_ns"] = reads_ns
            metadata["reads_fr"] = reads_fr
            metadata["reads_rf"] = reads_rf
            metadata["reads_ff"] = reads_ff

        os.symlink(fastq_node.identifier, out_path)
        return metadata
        
        
    def set_out_attributes(self, antecedents, out_attributes):
        #import os
        #from genomicode import config
        #from genomicode import filelib
        from genomicode import alignlib
        from Betsy import module_utils

        group_node, fastq_node, reference_node = antecedents
        fastq_files = module_utils.find_merged_fastq_files(
            group_node.identifier, fastq_node.identifier)
        assert fastq_files, "No fastq files."
        ref = alignlib.create_reference_genome(reference_node.identifier)

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
        orient, reads_ns, reads_fr, reads_rf, reads_ff = x
        #orientation = "paired"
        #if x:
        #    orientation = "paired_%s" % x
        #attrs["orientation"] = orientation
        attrs["orientation"] = "paired_%s" % orient
        key = (group_node.identifier, fastq_node.identifier,
               reference_node.identifier)
        self.cache[key] = reads_ns, reads_fr, reads_rf, reads_ff
        return attrs
    

    def name_outfile(self, antecedents, user_options):
        # Maintain the same extension so I know what file type this is.
        import os
        group_node, fastq_node, reference_node = antecedents
        f, ext = os.path.splitext(group_node.identifier)
        return "sample_group%s" % ext


def get_paired_orientation_rseqc(reference_bed, bam_filename):
    from genomicode import alignlib
    from genomicode import filelib
    from genomicode import parallel
    from Betsy import module_lib as mlib

    script = alignlib.find_rseqc_script("infer_experiment.py")
    filelib.assert_exists_nz(reference_bed)
    filelib.assert_exists_nz(bam_filename)

    # RSeQC scripts use #!/usr/bin/python, which may not be the right
    # one.  Use the python on the path.
    cmd = [
        "python",
        mlib.sq(script),
        "-r", mlib.sq(reference_bed)
        "-i", mlib.sq(bam_filename)
        ]
    cmd = " ".join(cmd)
    x = parallel.sshell(cmd)

    print x
    import sys; sys.exit(0)


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


def get_paired_orientation_bowtie2(
    reference_genome, filename1, filename2, outpath=None):
    # Return tuple of ("ff", "fr", or "rf"; reads_ns; reads_fr;
    # reads_rf; reads_ff).
    import os
    import shutil
    import tempfile
    import multiprocessing
    #from genomicode import genomelib
    from genomicode import alignlib
    from genomicode import parallel

    # Strategy: run bowtie2 in all orientations.  Return the one with
    # most reads aligned.  Do with a subset of the data, so this
    # doesn't take a long time.

    # 100 is too low.  Gave is wrong result (fr instead of rf) on the
    # Thunderbolts miSEQ data.  Minimum number that gives right answer
    # is 250.
    # NUM_READS  Time (s)  1 core
    #     50       2.4
    #    100       2.4
    #    250       2.5
    #    500       2.7
    #   1000       2.9
    #   2000       3.3
    #   5000       4.6
    #  10000       7.8
    # 100000      52.6
    NUM_READS = 1000

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
        sam_ns = os.path.join(path, "orient_ns.sam")
        log_ff = os.path.join(path, "orient_ff.log")
        log_fr = os.path.join(path, "orient_fr.log")
        log_rf = os.path.join(path, "orient_rf.log")
        log_ns = os.path.join(path, "orient_ns.log")


        nc = multiprocessing.cpu_count()
        nc = int(max(nc/4.0, 1))
        nc = 1
        x1 = alignlib.make_bowtie2_command(
            reference_genome, fastq_file1=filename1,
            fastq_file2=filename2, sam_file=sam_ff, orientation="ff",
            max_reads=NUM_READS, num_threads=nc)
        x2 = alignlib.make_bowtie2_command(
            reference_genome, fastq_file1=filename1,
            fastq_file2=filename2, sam_file=sam_fr, orientation="fr",
            max_reads=NUM_READS, num_threads=nc)
        x3 = alignlib.make_bowtie2_command(
            reference_genome, fastq_file1=filename1,
            fastq_file2=filename2, sam_file=sam_rf, orientation="rf",
            max_reads=NUM_READS, num_threads=nc)
        x4 = alignlib.make_bowtie2_command(
            reference_genome, fastq_file1=filename1,
            fastq_file2=filename2, sam_file=sam_ns, orientation=None,
            max_reads=NUM_READS, num_threads=nc)
        x1 += " >& %s" % log_ff
        x2 += " >& %s" % log_fr
        x3 += " >& %s" % log_rf
        x4 += " >& %s" % log_ns
        commands = [x1, x2, x3, x4]

        parallel.pshell(commands)

        # Read the results.
        output_ff = alignlib.parse_bowtie2_output(log_ff)
        output_fr = alignlib.parse_bowtie2_output(log_fr)
        output_rf = alignlib.parse_bowtie2_output(log_rf)
        output_ns = alignlib.parse_bowtie2_output(log_ns)
    finally:
        if tempdir is not None and os.path.exists(tempdir):
            shutil.rmtree(tempdir)

    reads_ff = output_ff["concordant_reads"]
    reads_fr = output_fr["concordant_reads"]
    reads_rf = output_rf["concordant_reads"]
    reads_ns = output_ns["concordant_reads"]
    assert type(reads_ff) is type(0)

    orient = [
        (reads_ff, "ff"),
        (reads_fr, "fr"),
        (reads_rf, "rf"),
        #(reads_ns, None),
        ]
    orient.sort(reverse=True)

    # Debug:
    if False:
        print orient
        raise AssertionError

    # If highest is within 10% of the un-stranded one, then it's
    #cutoff = reads_ns * 0.10
    #if reads_ns >= orient[3][0] - reads_ns*0.10:
    #    return None
    return orient[0][-1], reads_ns, reads_fr, reads_rf, reads_ff
