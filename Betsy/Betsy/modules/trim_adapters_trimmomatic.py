from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from Betsy import module_utils

        fastq_node, sample_node = antecedents

        fastq_path = fastq_node.identifier
        assert os.path.exists(fastq_path)
        assert os.path.isdir(fastq_path)

        adapters_filename = user_options.get("adapters_fasta")
        assert adapters_filename, "MISSING: adapters_fasta"
        assert filelib.exists_nz(adapters_filename), \
               "Not found: %s" % adapters_filename
               #"Not found: %s" % os.path.realpath(adapters_filename)
        
        module_utils.safe_mkdir(out_path)

        # Find the merged fastq files.
        x = module_utils.find_merged_fastq_files(
            sample_node.identifier, fastq_path)
        fastq_files = x
        assert fastq_files, "I could not find any FASTQ files."

        jobs = []
        for x in fastq_files:
            sample, pair1, pair2 = x
            p1, f1 = os.path.split(pair1)
            trimmed1 = os.path.join(out_path, f1)
            trimmed2 = None
            if pair2:
                p2, f2 = os.path.split(pair2)
                trimmed2 = os.path.join(out_path, f2)
            unpaired1 = os.path.join(out_path, "unpaired_1.fasta")
            unpaired2 = os.path.join(out_path, "unpaired_2.fasta")
            log_filename = os.path.join(out_path, "%s.log" % sample)
            x = sample, pair1, pair2, trimmed1, trimmed2, \
                unpaired1, unpaired2, log_filename
            jobs.append(x)

        sq = module_utils.shellquote
        commands = []
        for x in jobs:
            sample, pair1, pair2, trimmed1, trimmed2, unpaired1, unpaired2, \
                    log_filename = x
            nc = max(1, num_cores/len(jobs))
            x = _make_trimmomatic_cmd(
                pair1, pair2, trimmed1, trimmed2, unpaired1, unpaired2,
                adapters_filename, num_threads=nc)
            x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)

        # XXX test what if this is broken.
        module_utils.run_parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            sample, pair1, pair2, trimmed1, trimmed2, unpaired1, unpaired2, \
                    log_filename = x
            # Make sure outfile created.
            assert module_utils.exists_nz(trimmed1), \
                   "Missing: %s" % trimmed1
            if trimmed2:
                assert module_utils.exists_nz(trimmed2), \
                       "Missing: %s" % trimmed2
            x = open(log_filename).read()
            assert not x.startswith("Usage:"), "usage problem"

        
    def name_outfile(self, antecedents, user_options):
        return "fastq_trimmed"


def _make_trimmomatic_cmd(
    pair1, pair2, trimmed1, trimmed2, unpaired1, unpaired2, adapters_filename,
    num_threads=None):
    import os
    from genomicode import config
    from Betsy import module_utils

    assert os.path.exists(adapters_filename)
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 200

    #java -jar /usr/local/bin/trimmomatic.jar SE -phred33 fastq01.txt\
    #  fastq21.txt ILLUMINACLIP:NEBNext_SR.fa:2:30:10 LEADING:3 TRAILING:3 \
    #  SLIDINGWINDOW:4:15 MINLEN:15

    #TRIM=/usr/local/bin/trimmomatic.jar
    #CLIP=/home/jchang/src/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa
    #java -jar $TRIM PE -phred33 \
    #  ${b}_R1_001.fastq ${b}_R2_001.fastq \
    #  ${b}_R1_001.fastq.paired ${b}_R1_001.fastq.unpaired \
    #  ${b}_R2_001.fastq.paired ${b}_R2_001.fastq.unpaired \
    #  ILLUMINACLIP:$CLIP:2:30:10 \
    #  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    orientation = "SE"
    if pair2:
        orientation = "PE"
        assert pair2 and trimmed2 and unpaired1 and unpaired2

    trimmomatic = module_utils.which_assert(config.trimmomatic_jar)
    sq = module_utils.shellquote
    cmd = [
        "java",
        "-jar", sq(trimmomatic),
        orientation,
        ]
    if num_threads:
        cmd += ["-threads", str(num_threads)]
    cmd += ["-phred33"]
    if orientation == "SE":
        cmd += [pair1, trimmed1]
    else:
        cmd += [
            pair1, pair2,
            trimmed1, unpaired1,
            trimmed2, unpaired2,
            ]
    cmd += [
        # What if there are spaces in adapters_filename?
        "ILLUMINACLIP:%s:2:30:10" % adapters_filename,
        "LEADING:3",
        "TRAILING:3",
        "SLIDINGWINDOW:4:15",
        "MINLEN:15",
        ]
    return " ".join(cmd)
