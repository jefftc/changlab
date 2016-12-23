"""Next generation sequencing genome alignment.

Classes:
ReferenceGenome
HTSeqCountOutput

Methods:
create_reference_genome
standardize_reference_genome

find_bai_file

get_samtools_version
get_bcfutils_version
get_vcfutils_version
get_radia_version

call_samtools_idxstats

get_bowtie1_version
make_bowtie1_command
parse_bowtie1_output

get_bowtie2_version
make_bowtie2_command
parse_bowtie2_output

get_tophat_version
make_tophat_command
parse_tophat_align_summary

get_bwa_version
make_bwa_mem_command
make_bwa_aln_command

get_rsem_version
make_rsem_command
find_rsem_result_files
find_reference_stem

get_STAR_version
make_STAR_command
make_STAR_index_command
assert_is_STAR_reference

make_htseq_count_command
parse_htseq_count_output

find_picard_jar
make_GATK_command

make_platypus_command

make_annovar_command

get_muse_version

find_rseqc_script

gtf_to_bed

clean_mutect_vcf
clean_varscan_vcf
clean_strelka_vcf
clean_muse_vcf
clean_radia_vcf

# identify_caller     Figure out caller of a VCF file.

"""
# _create_reference_genome_path
# _is_subset

class ReferenceGenome:
    # Contains information about a ReferenceGenome.
    # Should only be instantiated by create_reference_genome
    #
    # Members:
    # name             Name of reference.
    # fasta_file_full  Full path of fasta file.  <name>.fa or <name>.fasta
    #
    # path             Path where the files exist.
    # fasta_file       Fasta file, relative to path.
    # dict_file        dict file.
    # samtools_index   File for samtools index.
    # bowtie1_indexes  Files for bowtie1 index, relative to path.
    #                  Only provided if indexes are complete.
    # bowtie2_indexes  Files for bowtie2 index, relative to path.
    #                  Only provided if indexes are complete.
    # bwa_indexes      Files for bwa index, relative to path.
    #                  Only provided if indexes are complete.
    # rsem_indexes     Files for rsem index, relative to path.
    #                  Only provided if indexes are complete.
    # other_files      Other unclassified files.

    def __init__(
        self, name, fasta_file_full, dict_file, samtools_index,
        bowtie1_indexes, bowtie2_indexes, bwa_indexes, rsem_indexes,
        other_files):
        import os

        p, f = os.path.split(fasta_file_full)

        self.name = name
        self.fasta_file_full = fasta_file_full

        self.path = p
        self.fasta_file = f
        self.dict_file = dict_file
        self.samtools_index = samtools_index
        self.bowtie1_indexes = bowtie1_indexes[:]
        self.bowtie2_indexes = bowtie2_indexes[:]
        self.bwa_indexes = bwa_indexes[:]
        self.rsem_indexes = rsem_indexes[:]
        self.other_files = other_files[:]


def _create_reference_genome_path(path, name=None):
    # path should be the path where all the files are contained.
    # If there are multiple reference genomes in the same path,
    # please provide name to specify which one to get.
    import os

    opj = os.path.join

    # Files:
    # <path>/<name>.[fa|fasta]       Fasta file.
    ## <path>/<name>.[fa|fasta].dict  dict file (GATK wants <path>/<name>.dict)
    # <path>/<name>.dict             dict file (GATK wants this)
    # <path>/<name>.[fa|fasta].fai   Samtools index.
    # <path>/<name>.[1234].ebwt      Bowtie1 index.
    # <path>/<name>.rev.[12].ebwt
    # <path>/<name>.[1234].bt2       Bowtie2 index.
    # <path>/<name>.rev.[12].bt2
    # <path>/<name>.[fa|fasta].amb   BWA
    # <path>/<name>.[fa|fasta].ann
    # <path>/<name>.[fa|fasta].bwt
    # <path>/<name>.[fa|fasta].pac
    # <path>/<name>.[fa|fasta].sa
    # <path>/<name>.chrlist          RSEM
    # <path>/<name>.grp
    # <path>/<name>.idx.fa
    # <path>/<name>.n2g.idx.fa
    # <path>/<name>.seq
    # <path>/<name>.ti
    # <path>/<name>.transcripts.fa
    # RSEM may also include Bowtie1 and Bowtie2 indexes.

    files = os.listdir(path)
    ## Make sure there are no directories.
    #for x in files:
    #    filename = opj(path, x)
    #    assert os.path.isfile(filename)  # allow symlinks?
    # Ignore the directories.
    files = [x for x in files if not os.path.isdir(opj(path, x))]

    # Find the fasta files.
    x = [x for x in files
         if x.lower().endswith(".fa") or x.lower().endswith(".fasta")]
    # Filter out the known RSEM index files with ".fa" extensions.
    x = [x for x in x if not x.lower().endswith(".idx.fa")]
    x = [x for x in x if not x.lower().endswith(".transcripts.fa")]
    fasta_files = x
    # Find the names of the fasta files.
    x = [os.path.splitext(x)[0] for x in fasta_files]
    names = x
    assert fasta_files, "No fasta files found: %s" % path
    if name is None:
        uniq_names = {}.fromkeys(names).keys()
        assert len(uniq_names) == 1, "Multiple fasta files found."
        name = names[0]
    assert name in names, "Reference genome not found: %s" % name

    # May have multiple fasta_files if same file ends with ".fa" and
    # ".fasta".  This might be necessary because different software
    # have different requirements.
    i = names.index(name)
    fasta_file = fasta_files[i]
    fasta_file_full = os.path.join(path, fasta_file)

    # Make a list of all known index files.
    all_bowtie1 = [
        "%s.1.ebwt" % name, "%s.2.ebwt" % name, "%s.3.ebwt" % name,
        "%s.4.ebwt" % name, "%s.rev.1.ebwt" % name, "%s.rev.2.ebwt" % name,
        ]
    all_bowtie2 = [
        "%s.1.bt2" % name, "%s.2.bt2" % name, "%s.3.bt2" % name,
        "%s.4.bt2" % name,  "%s.rev.1.bt2" % name,  "%s.rev.2.bt2" % name,
        ]
    all_rsem = [
        "%s.chrlist" % name, "%s.grp" % name, "%s.idx.fa" % name,
        "%s.n2g.idx.fa" % name, "%s.seq" % name, "%s.ti" % name,
        "%s.transcripts.fa" % name,
        ]
    all_bwa = [
        "%s.amb" % fasta_file, "%s.ann" % fasta_file, "%s.bwt" % fasta_file,
        "%s.pac" % fasta_file, "%s.sa" % fasta_file,
        ]

    # Identify each one of the files.
    samtools_index = None
    dict_file = None
    bowtie1_indexes = []
    bowtie2_indexes = []
    bwa_indexes = []
    rsem_indexes = []
    other_files = []
    for file_ in files:
        # Ignore files not from this reference genome.
        if not file_.startswith("%s." % name):
            continue
        #lfile = file_.lower()
        #if lfile.endswith(".fa") or lfile.endswith(".fasta"):
        #    # Already handled fasta files.
        #    continue
        if file_ == fasta_file:
            continue
        #elif file_ == "%s.dict" % fasta_file:
        elif file_ == "%s.dict" % name:
            assert dict_file is None
            dict_file = file_
        elif file_ == "%s.fai" % fasta_file:
            assert samtools_index is None
            samtools_index = file_
        elif file_ in all_bowtie1:
            bowtie1_indexes.append(file_)
        elif file_ in all_bowtie2:
            bowtie2_indexes.append(file_)
        elif file_ in all_bwa:
            bwa_indexes.append(file_)
        elif file_ in all_rsem:
            rsem_indexes.append(file_)
        else:
            other_files.append(file_)
    bowtie1_indexes.sort()
    bowtie2_indexes.sort()
    bwa_indexes.sort()
    rsem_indexes.sort()

    # Make sure indexes are complete.  If not, put incomplete
    # files in other_files.
    if bowtie1_indexes != sorted(all_bowtie1):
        other_files.extend(bowtie1_indexes)
        bowtie1_indexes = []
    if bowtie2_indexes != sorted(all_bowtie2):
        other_files.extend(bowtie2_indexes)
        bowtie2_indexes = []
    if bwa_indexes != sorted(all_bwa):
        other_files.extend(bwa_indexes)
        bwa_indexes = []
    if rsem_indexes != sorted(all_rsem):
        other_files.extend(rsem_indexes)
        rsem_indexes = []

    x = ReferenceGenome(
        name, fasta_file_full, dict_file, samtools_index, bowtie1_indexes,
        bowtie2_indexes, bwa_indexes, rsem_indexes, other_files)
    return x


def _is_subset(small_list, full_list):
    # Return whether small_list is a subset of full_list.
    return set(small_list).issubset(full_list)


def create_reference_genome(file_or_path, name=None):
    # Return a ReferenceGenome object.
    import os

    # If file_or_path is a file, then find the path of this file.
    path = file_or_path
    if not os.path.isdir(path):
        p, f = os.path.split(file_or_path)
        path = p
    return _create_reference_genome_path(path, name=name)


def standardize_reference_genome(
    in_file_or_path, out_path, use_symlinks=False):
    # Put a cleaner version of the reference genome in out_path.
    # 
    # in_file_or_path can be a reference FASTA file or path that
    # contains a reference FASTA file.
    # <file>         ->  <out_path>/<file>
    # <path>/<file>  ->  <out_path>/<file>
    #
    # Return a ReferenceGenome object.
    import os
    import filelib

    # If file_or_path is a file, then find the path of this file.
    path = in_file_or_path
    if os.path.isfile(path):
        p, f = os.path.split(in_file_or_path)
        path = p

    if use_symlinks:
        filelib.symlink_file_or_path_to_path(path, out_path)
    else:
        filelib.copy_file_or_path_to_path(path, out_path)
    x = create_reference_genome(out_path)
    return x


def find_bai_file(bam_filename):
    # Given a bam file, find the bai index file.  Return None if not
    # found.
    import os
    from genomicode import filelib
    
    filelib.assert_exists_nz(bam_filename)
    path, filename = os.path.split(bam_filename)
    stem, ext = os.path.splitext(filename)
    assert ext.lower() == ".bam"

    # <stem>.bam
    # <stem>.bam.bai
    # <stem>.bai
    file1 = os.path.join(path, "%s.bam.bai" % stem)
    file2 = os.path.join(path, "%s.bai" % stem)
    if filelib.exists(file1):
        return file1
    elif filelib.exists(file2):
        return file2
    return None
    

def get_radia_version():
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    opj = os.path.join
    sq = parallel.quote
    python = filelib.which_assert(config.python)
    radia_path = filelib.which_assert(config.radia_path)
    radia_py = opj(radia_path, "scripts", "radia.py")
    filelib.assert_exists_nz(radia_py)

    x = [
        sq(python),
        sq(radia_py),
        "--version",
        ]
    x = " ".join(x)
    x = parallel.sshell(x)
    x = x.strip()
    # v1.1.3
    return x


def get_samtools_version():
    # Returns full version string.
    # 1.2 (using htslib 1.2.1)
    # 0.1.18 (r982:295)
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    samtools = filelib.which_assert(config.samtools)
    x = parallel.sshell(samtools, ignore_nonzero_exit=True)
    x = x.strip()
    # Version: 1.2 (using htslib 1.2.1)
    # Version: 0.1.18 (r982:295)
    m = re.search(r"Version: ([\w\. \(\)-:]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def get_bcftools_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    bcftools = filelib.which_assert(config.bcftools)
    x = parallel.sshell(bcftools, ignore_nonzero_exit=True)
    x = x.strip()
    # Version: 1.2 (using htslib 1.2.1)
    m = re.search(r"Version: ([\w\. \(\)-]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def get_vcfutils_version():
    # Cannot get version.
    raise NotImplementedError


def call_samtools_idxstats(bam_filename):
    # Return list of (seqname, length, mapped, unmapped).
    import StringIO

    import parallel
    import filelib
    import config

    # BAM file should be indexed for maximum performance.
    filelib.assert_exists_nz(bam_filename)

    sq = parallel.quote
    samtools = filelib.which_assert(config.samtools)
    cmd = [
        sq(samtools),
        "idxstats",
        sq(bam_filename),
        ]
    # How to check if this is broken?
    x = parallel.sshell(cmd)
    handle = StringIO.StringIO(x)
    data = []
    for x in filelib.read_cols(handle):
        seqname, length, mapped, unmapped = x
        length = int(length)
        mapped = int(mapped)
        unmapped = int(unmapped)
        x = seqname, length, mapped, unmapped
        data.append(x)
    return data


def get_bowtie1_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    bowtie = filelib.which_assert(config.bowtie)
    x = parallel.sshell("%s --version" % bowtie, ignore_nonzero_exit=True)
    x = x.strip()
    # bowtie version 1.1.1
    m = re.search(r"version ([\w\.]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def make_bowtie1_command(
    reference_fa, sam_file, fastq_file1, fastq_file2=None,
    orientation=None, num_threads=None):
    # reference_fa is the full path to the fasta file.
    # Orientation must be None, "ff", "fr", "rf"
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    assert os.path.exists(fastq_file1)
    if fastq_file2:
        assert os.path.exists(fastq_file2)
    if orientation:
        assert orientation in ["ff", "fr", "rf"]
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100

    bowtie1 = filelib.which_assert(config.bowtie)

    # bowtie --sam -p <nthreads> <reference_base> <sample.fq>
    #   <sample.sam>
    # bowtie --sam --fr -p <nthreads> <reference_base> -1 <sample_1.fq> -2
    #   <sample_2.fq> <sample.sam>

    sq = parallel.quote
    cmd = [
        sq(bowtie1),
        "--sam",
        ]
    if orientation:
        cmd += ["--%s" % orientation]
    if num_threads:
        cmd += ["-p", str(num_threads)]
    stem = find_reference_stem(reference_fa)
    cmd += [sq(stem)]
    if not fastq_file2:
        cmd += [sq(fastq_file1)]
    else:
        cmd += [
            "-1", sq(fastq_file1),
            "-2", sq(fastq_file2),
            ]
    cmd += [sq(sam_file)]
    return " ".join(cmd)


def parse_bowtie1_output(filename):
    # Return a dictionary with keys:
    # reads_processed
    # aligned_reads
    #
    # Warning: Exhausted best-first chunk memory ...; skipping read
    # # reads processed: 62830141
    # # reads with at least one reported alignment: 28539063 (45.42%)
    # # reads that failed to align: 34291078 (54.58%)
    # Reported 28539063 paired-end alignments to 1 output stream(s)
    from genomicode import filelib

    results = {}
    for line in filelib.openfh(filename):
        if line.startswith("Warning:"):
            continue
        elif line.startswith("# reads processed"):
            x = line.split(":")
            assert len(x) == 2
            results["reads_processed"] = int(x[1])
        elif line.startswith("# reads with at least one reported alignment"):
            x = line.split(":")
            assert len(x) == 2
            x = x[1].strip()
            x = x.split(" ")
            assert len(x) == 2
            results["aligned_reads"] = int(x[0])
        elif line.startswith("# reads that failed to align"):
            pass
        elif line.startswith("Reported "):
            pass
        else:
            raise AssertionError, "Unknown line: %s" % line.strip()
    return results


def get_bowtie2_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    bowtie = filelib.which_assert(config.bowtie2)
    x = parallel.sshell("%s --version" % bowtie, ignore_nonzero_exit=True)
    x = x.strip()
    # /usr/local/bin/bowtie2-align-s version 2.2.4
    m = re.search(r"version ([\w\.]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def make_bowtie2_command(
    reference_fa, fastq_file1, fastq_file2=None, orientation=None,
    sam_file=None, bam_file=None, num_threads=None, skip_reads=None,
    max_reads=None):
    # Orientation must be None, "ff", "fr", "rf".
    # Will save to either sam file or bam file.  Either sam_file or
    # bam_file, or neither, must be given.
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    assert os.path.exists(fastq_file1)
    assert not (sam_file and bam_file)
    if fastq_file2:
        assert os.path.exists(fastq_file2)
    if orientation:
        assert orientation in ["ff", "fr", "rf"]
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100
    if skip_reads is not None:
        assert skip_reads >= 0
    if max_reads is not None:
        assert max_reads > 0

    bowtie2 = filelib.which_assert(config.bowtie2)

    # bowtie2 -p <nthreads> -x <reference_base> -U <sample.fq>
    #   -S <sample.sam>
    # bowtie2 -p <nthreads> -x <reference_base> -1 <sample_1.fq>
    #   -2 <sample_2.fq> --fr -S <sample.sam>

    sq = parallel.quote
    cmd = [
        sq(bowtie2),
        ]
    if num_threads:
        cmd += ["-p", str(num_threads)]
    stem = find_reference_stem(reference_fa)
    cmd += ["-x", sq(stem)]
    if skip_reads is not None:
        cmd += ["-s", str(skip_reads)]
    if max_reads is not None:
        cmd += ["--upto", str(max_reads)]
    if not fastq_file2:
        cmd += [
            "-U", sq(fastq_file1),
            ]
    else:
        cmd += [
            "-1", sq(fastq_file1),
            "-2", sq(fastq_file2),
            ]
    if orientation:
        cmd += ["--%s" % orientation]
    if sam_file:
        cmd += [
            "-S", sq(sam_file),
            ]
    cmd = " ".join(cmd)
    if bam_file:
        # samtools view -bS -o <bam_filename> <sam_filename>
        samtools = filelib.which_assert(config.samtools)
        x = [
            sq(samtools),
            "view",
            "-bS",
            "-o", sq(bam_file),
            "-",
            ]
        x = " ".join(x)
        cmd = "%s | %s" % (cmd, x)
    return cmd


def parse_bowtie2_output(filename):
    # Return a dictionary with keys:
    # reads_processed    BOTH
    # aligned_reads      BOTH    (aligned >= 1 time)
    # concordant_reads   paired  (concordant >= 1 time)
    #
    # 1000 reads_processed means there are 1000 possible pairs.

    # For single-end reads.
    # 20000 reads; of these:
    #   20000 (100.00%) were unpaired; of these:
    #      1247 (6.24%) aligned 0 times
    #      18739 (93.69%) aligned exactly 1 time
    #      14 (0.07%) aligned >1 times
    # 93.77% overall alignment rate

    # For paired-end reads.
    # 62830141 reads; of these:
    #   62830141 (100.00%) were paired; of these:
    #     6822241 (10.86%) aligned concordantly 0 times
    #     48713485 (77.53%) aligned concordantly exactly 1 time
    #     7294415 (11.61%) aligned concordantly >1 times
    #     ----
    #     6822241 pairs aligned concordantly 0 times; of these:
    #       3004782 (44.04%) aligned discordantly 1 time
    #     ----
    #     3817459 pairs aligned 0 times concordantly or discordantly; of these:
    #       7634918 mates make up the pairs; of these:
    #         3135692 (41.07%) aligned 0 times
    #         1870375 (24.50%) aligned exactly 1 time
    #         2628851 (34.43%) aligned >1 times
    # 97.50% overall alignment rate
    from genomicode import filelib

    concordant_1 = concordant_more = None
    results = {}
    for line in filelib.openfh(filename):
        if line.find("reads; of these:") >= 0:
            x = line.strip().split()
            results["reads_processed"] = int(x[0])
        elif line.find("overall alignment rate") >= 0:
            x = line.strip().split()
            # 93.77% overall alignment rate
            x = x[0]            # 93.77%
            assert x.endswith("%")
            x = x[:-1]          # 93.77
            x = float(x) / 100  # 0.9377
            assert "reads_processed" in results
            aligned = int(x * results["reads_processed"])
            results["aligned_reads"] = aligned
        elif line.find("aligned concordantly exactly 1 time") >= 0:
            x = line.strip().split()
            x = x[0]
            concordant_1 = int(x)
        elif line.find("aligned concordantly >1 times") >= 0:
            x = line.strip().split()
            x = x[0]
            concordant_more = int(x)
    x = open(filename).read()
    assert "reads_processed" in results, x
    assert "aligned_reads" in results, x

    if concordant_1 is not None:
        assert concordant_more is not None
    if concordant_more is not None:
        assert concordant_1 is not None
    if concordant_1 is not None:
        results["concordant_reads"] = concordant_1 + concordant_more
    return results


def get_tophat_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    tophat = filelib.which_assert(config.tophat)
    x = parallel.sshell("%s --version" % tophat, ignore_nonzero_exit=True)
    x = x.strip()
    # TopHat v2.0.13
    m = re.search(r"TopHat v([\w\.]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def make_tophat_command(
    reference_fa, out_path, fastq_file1, fastq_file2=None,
    gtf_file=None, transcriptome_fa=None, library_type=None,
    num_threads=None):
    # reference_fa is the full path to the fasta file.
    # transcriptome_fa should be tophat indexed fasta file of
    # transcriptome.
    # library_type must be None, "fr-unstranded", "fr-firststrand", or
    # "fr-secondstrand".
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    assert os.path.exists(fastq_file1)
    if fastq_file2:
        assert os.path.exists(fastq_file2)
    if gtf_file:
        assert os.path.exists(gtf_file)
    transcriptome_index = None
    if transcriptome_fa:
        # Make sure seems reasonable and bowtie indexed.
        p, f = os.path.split(transcriptome_fa)
        f, ext = os.path.splitext(f)
        name = f
        ref = create_reference_genome(transcriptome_fa, name=name)
        assert ref.bowtie2_indexes, \
               "transcriptome index must be indexed with bowtie2"
        x = ref.fasta_file_full
        if x.endswith(".fa"):
            x = x[:-3]
        elif x.endswith(".fasta"):
            x = x[:-6]
        else:
            raise AssertionError, "Unknown fasta extension: %s" % \
                  ref.fasta_file
        transcriptome_index = x
    assert gtf_file or transcriptome_index, \
           "Either gtf_file or transcriptome_index must be provided."
    assert library_type in [
        None, "fr-unstranded", "fr-firststrand", "fr-secondstrand"]
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100

    tophat = filelib.which_assert(config.tophat)

    # tophat [options]* <stem> <reads_1.fq> [<reads_2.fa>]
    # --bowtie1  Use bowtie1 instead of bowtie2.
    # -o <outdir>
    # -r <inner dict>
    # -p <num_threads>
    # --library-type fr-unstranded, fr-firststrand, fr-secondstrand

    sq = parallel.quote
    cmd = [
        sq(tophat),
        "-o", sq(out_path),
        ]
    if gtf_file:
        cmd += ["--GTF", sq(gtf_file)]
    if transcriptome_index:
        cmd += ["--transcriptome-index", sq(transcriptome_index)]
    if num_threads:
        cmd += ["-p", str(num_threads)]
    if library_type:
        cmd += ["--library-type", library_type]

    stem = find_reference_stem(reference_fa)
    cmd += [sq(stem)]
    cmd += [sq(fastq_file1)]
    if fastq_file2:
        cmd += [sq(fastq_file2)]
    return " ".join(cmd)


def parse_tophat_align_summary(filename):
    # Return a dictionary with keys:
    # reads_processed
    # aligned_reads
    from genomicode import filelib
    filelib.assert_exists_nz(filename)

    # This parses a paired end read file.  Have not implemented a
    # parser for single end results.
    # Format for paired end reads: (Spacing not preserved.)
    # Left reads:
    #   Input     :  78337847
    #   Mapped   :  78243095 (99.9% of input)
    #     of these:   6068262 ( 7.8%) have multiple alignments (83814 have >20)
    # Right reads:
    #   Input     :  78337847
    #   Mapped   :  77663143 (99.1% of input)
    #     of these:   5987549 ( 7.7%) have multiple alignments (83425 have >20)
    # 99.5% overall read mapping rate.
    #
    # Aligned pairs:  77575761
    #  of these:   5980080 ( 7.7%) have multiple alignments
    #              4017526 ( 5.2%) are discordant alignments
    # 93.9% concordant pair alignment rate.

    # Pull out just the aligned pairs, and concordant pair alignment
    # rate.
    aligned_pairs = None
    concordant_pair_alignment_rate = None
    for line in filelib.openfh(filename):
        line = line.strip()
        if not line:
            continue
        if line.startswith("Aligned pairs:"):
            x = line.split()
            assert len(x) == 3, repr(x)
            aligned_pairs = int(x[-1])
        elif line.find("concordant pair alignment rate") >= 0:
            x = line.split()
            assert len(x) == 5
            x = x[0]
            assert x.endswith("%")
            x = x[:-1]
            concordant_pair_alignment_rate = float(x)
    assert aligned_pairs, "Not found: Aligned pair"
    assert concordant_pair_alignment_rate, "Not found: alignment rate"

    results = {}
    results["reads_processed"] = aligned_pairs
    x = int(round(aligned_pairs*concordant_pair_alignment_rate/100.0))
    results["aligned_reads"] = x
    return results


def get_bwa_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    bwa = filelib.which_assert(config.bwa)
    x = parallel.sshell(bwa, ignore_nonzero_exit=True)
    x = x.strip()
    # Version: 0.7.12-r1039
    m = re.search(r"Version: ([\w\.-]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def make_bwa_mem_command(
    reference_fa, err_filename, fastq_file1,
    fastq_file2=None, sam_filename=None, bam_filename=None, num_threads=None):
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    assert os.path.exists(reference_fa)
    assert os.path.exists(fastq_file1)
    if fastq_file2:
        assert os.path.exists(fastq_file2)
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100
    assert sam_filename or bam_filename
    assert not (sam_filename and bam_filename)

    bwa = filelib.which_assert(config.bwa)

    # bwa mem -t <num_cores> ref.fa read1.fq read2.fq > aln-pe.sam
    # bwa mem -t <num_cores> ref.fa read1.fq read2.fq 2> <err_file>
    #   | samtools view -bS -o aln_pe.bam -
    sq = parallel.quote
    cmd = [sq(bwa), "mem"]
    if num_threads:
        cmd += ["-t", str(num_threads)]
    cmd += [sq(reference_fa)]
    cmd += [sq(fastq_file1)]
    if fastq_file2:
        cmd += [sq(fastq_file2)]
    cmd += ["2>", sq(err_filename),]
    if sam_filename:
        cmd += ["1>", sq(sam_filename),]
    elif bam_filename:
        samtools = filelib.which_assert(config.samtools)
        cmd += [
            "|", sq(samtools), "view", "-bS", "-o", sq(bam_filename), "-",
            ]
    else:
        raise AssertionError
    return " ".join(cmd)


def make_bwa_aln_command(
    reference_fa, fastq_filename, sai_filename, err_filename,
    num_threads=None):
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    assert os.path.exists(reference_fa)
    assert os.path.exists(fastq_filename)
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100
    bwa = filelib.which_assert(config.bwa)

    # bwa aln -t <num_cores> <reference.fa> <input.fastq> > <output.sai>
    sq = parallel.quote
    cmd = [sq(bwa), "aln"]
    if num_threads:
        cmd += ["-t", str(num_threads)]
    cmd += [sq(reference_fa)]
    cmd += [sq(fastq_filename)]
    cmd += [
        "1>", sq(sai_filename),
        "2>", sq(err_filename),
        ]
    return " ".join(cmd)


## def find_rsem_reference(search_path):
##     # Find the indexed reference genome at:
##     # <search_path> should be:
##     #   /ref/
##     # Where the actual index files are:
##     #   /ref/hg19.transcripts.fa (etc.)
##     #   <search_path>/<assembly>.transcripts.fa
##     # Return <search_path>/.../<assembly>
##     #
##     # This can be used as the <index_stem> parameter for
##     # rsem-calculate-expression.
##     from genomicode import filelib

##     x = filelib.list_files_in_path(search_path)
##     x = [x for x in x if x.lower().endswith(".idx.fa")]
##     x = [x for x in x if not x.lower().endswith(".n2g.idx.fa")]
##     assert x, "Cannot find rsem index."
##     assert len(x) == 1, "Found multiple potential rsem indexes."
##     x = x[0]
##     assert x.endswith(".idx.fa")
##     x = x[:-len(".idx.fa")]
##     index_stem = x

##     index_files = [
##         "%s.chrlist" % index_stem,
##         "%s.grp" % index_stem,
##         "%s.idx.fa" % index_stem,
##         "%s.n2g.idx.fa" % index_stem,
##         "%s.seq" % index_stem,
##         "%s.ti" % index_stem,
##         "%s.transcripts.fa" % index_stem,
##         ]
##     for filename in index_files:
##         assert exists_nz(filename), "Missing: %s" % filename
##     return index_stem


def get_rsem_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    rsem = filelib.which_assert(config.rsem_calculate)
    x = parallel.sshell("%s --version" % rsem, ignore_nonzero_exit=True)
    x = x.strip()
    # Current version is RSEM v1.2.19
    m = re.search(r"RSEM v([\w\.]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def make_rsem_command(
    reference_fa, sample_name, fastq_file1, fastq_file2=None,
    forward_prob=None, num_threads=None, output_genome_bam=True,
    bowtie_chunkmbs=None, align_with_bowtie2=False, align_with_star=False):
    # <sample_name> is prefix for all output files,
    #   e.g. <sample_name>.genes.results
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    # rsem-calculate-expression -p <num_cores> --output-genome-bam
    #   --paired-end <file1.fastq> <file2.fastq>
    #   <index_stem> <sample_name> >& <sample_name>.log"
    # For strand-specific:
    #   --forward-prob 1.0 (upstream from forward strand; secondstrand)
    #   --forward-prob 0.0 (upstream from reverse strand; firststrand)

    # If there are lots of memory warnings in the log files:
    #   Warning: Exhausted best-first chunk memory for read
    #   ST-J00106:110:H5NY5BBXX:6:1101:18203:44675 1:N:0:1/1 (patid
    #   2076693); skipping read
    # May need to increase --chunkmbs.
    
    assert os.path.exists(fastq_file1)
    if fastq_file2:
        assert os.path.exists(fastq_file2)
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100
    assert forward_prob in [None, 0.0, 0.5, 1.0]
    if bowtie_chunkmbs:
        assert bowtie_chunkmbs > 0 and bowtie_chunkmbs < 100*1024
    assert not (align_with_bowtie2 and align_with_star)

    rsem_calculate = filelib.which_assert(config.rsem_calculate)

    sq = parallel.quote
    cmd = [
        sq(rsem_calculate),
        ]
    if num_threads:
        cmd += ["-p", str(num_threads)]
    if bowtie_chunkmbs:
        cmd += ["--bowtie-chunkmbs", str(bowtie_chunkmbs)]
    if align_with_bowtie2:
        cmd += ["--bowtie2"]
        bowtie2 = filelib.which_assert(config.bowtie2)
        # RSEM wants the path that contains the bowtie2 executable.
        bowtie2 = os.path.split(bowtie2)[0]
        cmd += ["--bowtie2-path", sq(bowtie2)]
    if align_with_star:
        cmd += ["--star"]
        STAR = filelib.which_assert(config.STAR)
        # RSEM wants the path that contains the executable.
        STAR = os.path.split(STAR)[0]
        cmd += ["--star-path", sq(STAR)]
    if output_genome_bam:
        cmd += ["--output-genome-bam"]
    if forward_prob is not None:
        cmd += ["--forward-prob", str(forward_prob)]
    if not fastq_file2:
        cmd += [sq(fastq_file1)]
    else:
        cmd += [
            "--paired-end",
            sq(fastq_file1),
            sq(fastq_file2),
            ]
    stem = find_reference_stem(reference_fa)
    cmd += [sq(stem), sq(sample_name)]
    return " ".join(cmd)


def find_rsem_result_files(search_path):
    # Find the files with rsem gene expression estimates.  Return list
    # of (sample, gene_filename, isoform_filename).  gene_filename or
    # isoform_filename can be None.
    import os
    from genomicode import filelib

    # Look for files in the format:
    # <sample>.genes.results
    # <sample>.isoforms.results
    x0 = filelib.list_files_in_path(search_path)
    x1 = [x for x in x0 if x.endswith(".genes.results")]
    x2 = [x for x in x0 if x.endswith(".isoforms.results")]
    sample2files = {}  # sample -> gene_filename, isoform_filename
    for gene_filename in x1:
        p, f = os.path.split(gene_filename)
        sample = f.replace(".genes.results", "")
        assert sample not in sample2files
        sample2files[sample] = gene_filename, None
    for isoform_filename in x2:
        p, f = os.path.split(isoform_filename)
        sample = f.replace(".isoforms.results", "")
        gene_filename, x = sample2files.get(sample, (None, None))
        assert x is None
        sample2files[sample] = gene_filename, isoform_filename
    data = []  # list of (sample, gene_filename, isoform_filename)X
    for sample in sorted(sample2files):
        gene_filename, isoform_filename = sample2files[sample]
        x = sample, gene_filename, isoform_filename
        data.append(x)
    return data


def find_reference_stem(ref_fasta):
    # ref_fasta is the full path to the fasta file for the reference
    # genome.
    #
    # <path>/<name>.[fa|fasta]
    # <path>/<name>.[1234].ebwt
    # <path>/<name>.rev.[12].ebwt
    # [etc...]
    #
    # <reference_stem> is <path>/<name>.
    import os

    assert os.path.exists(ref_fasta)
    ref_fasta = os.path.realpath(ref_fasta)
    stem, ext = os.path.splitext(ref_fasta)
    assert ext in [".fa", ".fasta"]
    return stem


def get_STAR_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    star = filelib.which_assert(config.STAR)
    x = parallel.sshell(star, ignore_nonzero_exit=True)
    x = x.strip()
    # ### versions
    # versionSTAR             020201
    # int>0: STAR release numeric ID. Please do not change this value!
    m = re.search(r"versionSTAR\s+([\w\. \(\)-]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def make_STAR_command(
    reference_path, out_prefix, num_cores, is_stranded, pair1, pair2,
    log_filename):
    # out_prefix should be full path.
    from genomicode import filelib
    from genomicode import parallel
    from genomicode import config

    STAR = filelib.which_assert(config.STAR)

    sq = parallel.quote
    x = [
        sq(STAR),
        "--genomeDir", sq(reference_path),
        "--outFileNamePrefix", sq(out_prefix),
        "--runThreadN", num_cores,
        "--outSAMtype", "BAM Unsorted",
        ]
    if not is_stranded:
        x += ["--outSAMstrandField", "intronMotif"]
    x += ["--readFilesIn", sq(pair1)]
    if pair2:
        x += [sq(pair2)]
    x = " ".join(map(str, x))
    x = "%s >& %s" % (x, sq(log_filename))
    return x


def make_STAR_index_command(
    fasta_filename, out_path, gtf_file=None, sjdb_files=None, num_cores=None):
    # Either gtf_file or sjdb_file must be given, but not both.
    # sjdb_files should be a list of the files for multi-sample 2-pass
    # mapping.
    import os
    from genomicode import filelib
    from genomicode import config
    from genomicode import parallel

    filelib.assert_exists_nz(fasta_filename)
    #filelib.assert_exists(out_path)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    assert gtf_file or sjdb_files
    assert not (gtf_file and sjdb_files)
    if gtf_file:
        filelib.assert_exists_nz(gtf_file)
    if sjdb_files:
        assert type(sjdb_files) is type([]), \
               "sjdb_files should be a list of files."
        for x in sjdb_files:
            filelib.assert_exists_nz(x)

    STAR = filelib.which_assert(config.STAR)

    # STAR --runThreadN 40 --runMode genomeGenerate --genomeDir test05
    #  --genomeFastaFiles <file.fasta> \
    #  --sjdbGTFfile $GTF

    # --sjdbOverhang.  According to manual:
    # In most cases, the default value of 100 will work as well as the
    # ideal value.

    sq = parallel.quote

    x = [
        sq(STAR),
        ]
    if num_cores is not None:
        x += ["--runThreadN", num_cores]
    x += [
        "--runMode", "genomeGenerate",
        "--genomeDir", sq(out_path),
        "--genomeFastaFiles", sq(fasta_filename),
        ]
    if gtf_file:
        x += ["--sjdbGTFfile", sq(gtf_file)]
    else:
        assert sjdb_files
        x += ["--sjdbFileChrStartEnd"] + [sq(x) for x in sjdb_files]
    x = " ".join(map(str, x))
    #x = "%s >& out.txt" % " ".join(map(str, x))
    return x


def assert_is_STAR_reference(reference_path):
    import os
    from genomicode import filelib

    # Check to make sure index was created successfully.
    files = [
        "chrLength.txt",
        "chrNameLength.txt",
        "chrName.txt",
        "chrStart.txt",
        #"exonGeTrInfo.tab",    # gtf_file only.  Missing for sjdb_files.
        #"exonInfo.tab",        # gtf_file only.  Missing for sjdb_files.
        #"geneInfo.tab",        # gtf_file only.  Missing for sjdb_files.
        "Genome",
        "genomeParameters.txt",
        "SA",
        "SAindex",
        "sjdbInfo.txt",
        #"sjdbList.fromGTF.out.tab",  # gtf_file only.  Missing for sjdb_files.
        "sjdbList.out.tab",
        #"transcriptInfo.tab",        # gtf_file only.  Missing for sjdb_files.
        ]
    x = [os.path.join(reference_path, x) for x in files]
    filelib.assert_exists_nz_many(x)


def make_htseq_count_command(
    bam_file, gtf_file, sort_order, stranded, mode=None):
    # sort_order should be "name" or "pos".
    # stranded should be "yes", "no", "reverse".
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    assert sort_order in ["name", "pos"]
    assert stranded in ["yes", "no", "reverse"]
    assert mode is None or mode in \
           ["union", "intersection-strict", "intersection-nonempty"]

    assert os.path.exists(bam_file)
    assert os.path.exists(gtf_file)

    # htseq-count [options] <alignment_file> <gff_file>
    # -f <format>    bam
    # -r <order>     name (default), pos.  How sorted.
    # -s <stranded>  yes, no, reverse
    # -t <feature type>
    # -m <mode>
    htseq_count = filelib.which_assert(config.htseq_count)

    sq = parallel.quote
    cmd = [
        sq(htseq_count),
        "-f", "bam",
        "-r", sort_order,
        "-s", stranded,
        ]
    if mode:
        cmd += ["-m", mode]
    cmd += [
        bam_file,
        gtf_file,
        ]
    return " ".join(cmd)


class HTSeqCountOutput:
    def __init__(
        self, counts, no_feature, ambiguous, too_low_aQual,
        not_aligned, alignment_not_unique, warnings, errors):
        # counts is a dictionary of gene_id -> count.
        # warnings is a list of warning lines
        # errors is a list of error lines
        self.counts = counts.copy()
        self.no_feature = no_feature        # count not be assigned to feature
        self.ambiguous = ambiguous          # multiple features
        self.too_low_aQual = too_low_aQual  # low quality
        self.not_aligned = not_aligned      # no alignment in SAM file
        self.alignment_not_unique = alignment_not_unique  # multiple alignment
        self.warnings = warnings[:]
        self.errors = errors[:]


def parse_htseq_count_output(file_or_handle):
    # Return an HTSeqCountOutput object.
    import os

    #filename = None
    handle = file_or_handle
    if type(file_or_handle) is type(""):
        #filename = file_or_handle
        assert os.path.exists(file_or_handle)
        handle = open(file_or_handle)

    counts = {}
    meta = {}
    warnings = []
    errors = []

    # Default values, in case there are warnings or errors.
        # __no_feature    343859
        # __ambiguous     279435
        # __too_low_aQual 247071
        # __not_aligned   9744
        # __alignment_not_unique  0


    for line in handle:
        # 100000 GFF lines processed.
        if line.rstrip().endswith("processed."):
            continue
        # Warning: Mate records missing for 128 records; first such
        # record: <SAM_Alignment object: Paired-end read
        # 'SRR988443.108873869' aligned to MT:[11627,11677)/->.
        if line.startswith("Warning:"):
            warnings.append(line.strip())
            continue
        #err_msg = line.strip()
        #if filename:
        #    err_msg = "%s: %s" % line.strip()
        #assert not line.startswith("Error occured"), err_msg
        if line.startswith("Error"):
            errors.append(line.strip())
            continue
        x = line.rstrip("\r\n").split("\t")
        assert len(x) == 2, "Unrecognized line: %s" % repr(line.rstrip("\r\n"))
        gene_id, count = x
        gene_id, count = gene_id.strip(), int(count)
        # __no_feature    343859
        # __ambiguous     279435
        # __too_low_aQual 247071
        # __not_aligned   9744
        # __alignment_not_unique  0
        if gene_id.startswith("__"):
            name = gene_id[2:]
            meta[name] = count
            continue
        counts[gene_id] = count
    assert len(meta) == 5

    x = HTSeqCountOutput(counts, warnings=warnings, errors=errors, **meta)
    return x


def find_picard_jar(jar_name):
    # jar_name should not include ".jar",
    # e.g. "AddOrReplaceReadGroups".
    import os
    from genomicode import config

    picard_path = config.picard
    assert os.path.exists(picard_path)
    jar_filename = os.path.join(picard_path, "%s.jar" % jar_name)
    assert os.path.exists(jar_filename), "File not found: %s" % jar_filename
    return jar_filename


def _make_java_command(config_name, params, num_dashes, memory=None,
                       java_path=None):
    # params a a dictionary of <key> : <value> to pass to the java
    # command.
    # <value> of None means no argument.
    # A dash is prepended to each key.
    # Special key:
    # _UNHASHABLE  list of (key, value) tuples
    #
    # memory is the number of gb of RAM to allocate for java.  Should
    # be an integer.  Default is 5.
    from genomicode import config
    from genomicode import parallel
    from genomicode import filelib

    UNHASHABLE = "_UNHASHABLE"
    assert num_dashes >= 1 and num_dashes <= 2

    if memory is None:
        memory = 5
    assert type(memory) is type(0)
    assert memory >= 1 and memory <= 256

    jarfile = getattr(config, config_name)
    filelib.assert_exists_nz(jarfile)

    sq = parallel.quote
    java = "java"
    if java_path is not None:
        filelib.which_assert(java_path)
        java = sq(java_path)

    cmd = [
        java,
        "-Xmx%dg" % memory,
        "-Djava.io.tmpdir=.", 
        "-jar", sq(jarfile),
        ]
    x1 = params.get(UNHASHABLE, [])
    assert type(x1) in [type([]), type(())], \
           "%s should be list of key-value tuples" % UNHASHABLE
    x2 = list(params.iteritems())
    all_params = x1 + x2

    DASH = "-"*num_dashes
    for (key, value) in all_params:
        if key == UNHASHABLE:
            continue
        #assert key[0] in "A-Za-z_-", "Possibly bad name: %s" % key
        if value is None:
            cmd.append("%s%s" % (DASH, key))
        else:
            cmd.extend(["%s%s" % (DASH, key), sq(value)])
    return " ".join(map(str, cmd))


def make_GATK_command(**params):
    # SplitNCigarReads complains with only 5 Gb RAM:
    # There was a failure because you did not provide enough memory to
    # run this program.  See the -Xmx JVM argument to adjust the
    # maximum heap size provided to Java
    # GATK commands by default just use one dash.
    
    return _make_java_command("gatk_jar", params, 1, memory=20)


def make_MuTect_command(**params):
    from genomicode import config

    # Mutect 1.1.4 only works with Java 6.
    # Mutect 1.1.7 (from Github) only works with Java 7.
    name = "mutect_java"
    assert hasattr(config, name), "Missing from genomicode config: %s" % name
    java_path = getattr(config, name)

    return _make_java_command("mutect_jar", params, 2, java_path=java_path)


def make_platypus_command(
    bam_file, ref_file, log_file=None, out_file=None, buffer_size=None,
    max_reads=None):
    from genomicode import config
    from genomicode import parallel
    from genomicode import filelib

    if max_reads is not None:
        assert max_reads > 2 and max_reads < 1E10
    if buffer_size is not None:
        assert buffer_size > 1E4 and buffer_size < 1E8

    # nCPU command seems to be broken.  May cut off last line.
    num_cores = None
    if num_cores is not None:
        assert num_cores >= 1 and num_cores < 256

    # /usr/local/bin/Platypus/Platypus.py callVariants
    #   --bamFiles $i
    #   --refFile ../index/erdman.fa
    #   --output $j
    #   --logFileName <filename>
    platypus = filelib.which_assert(config.platypus)

    sq = parallel.quote
    cmd = [
        sq(platypus),
        "callVariants",
        "--bamFiles", sq(bam_file),
        "--refFile", sq(ref_file),
        ]
    if log_file:
        cmd += ["--logFileName", log_file]
    if out_file:
        cmd += ["--output", out_file]
    if num_cores:
        cmd += ["--nCPU", num_cores]
    if buffer_size:
        cmd += ["--bufferSize", buffer_size]
    if max_reads:
        cmd += ["--maxReads", max_reads]
    return " ".join(map(str, cmd))


def make_annovar_command(
    in_filename, log_filename, out_filestem, buildver, vcf_input=True):
    # in_filename is an input VCF file.  out_filestem is just the file
    # with no extension.  buildver, e.g. "hg19"
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    # list of (name, operation).
    # These are just for buildver hg19.
    protocols = [
        ("refGene", "g"), ("cytoBand", "r"), ("genomicSuperDups", "r"),
        ("esp6500siv2_all", "f"), ("snp138", "f"),
        ("ljb26_all", "f"),
        ("1000g2015aug_all", "f"), ("1000g2015aug_afr", "f"),
        ("1000g2015aug_eas", "f"), ("1000g2015aug_eur", "f"),
        ]
    if buildver != "hg19":
        raise NotImplementedError

    # P1=refGene,cytoBand,genomicSuperDups
    # P2=esp6500siv2_all,snp138,ljb26_all
    # P3=1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur
    # table_annovar.pl -buildver hg19 -remove -vcfinput
    #   -protocol $P1,$P2,$P3 -operation g,r,r,f,f,f,f,f,f,f
    #   -out $j $i humandb/

    table_annovar = filelib.which_assert(config.table_annovar)
    annodb = config.annovar_db
    assert os.path.exists(annodb)
    assert os.path.isdir(annodb)

    sq = parallel.quote
    x1 = [x[0] for x in protocols]
    x2 = [x[1] for x in protocols]
    x = [
        sq(table_annovar),
        "-buildver", buildver,
        "-remove",
        ]
    if vcf_input:
        x += ["-vcfinput"]
    x += [
        "-protocol", ",".join(x1),
        "-operation", ",".join(x2),
        "-out", sq(out_filestem),
        sq(in_filename),
        sq(annodb),
        ]
    x = " ".join(x)
    x = "%s >& %s" % (x, log_filename)
    return x


def find_rseqc_script(name):
    # example name "infer_experiment.py"
    import os
    from genomicode import config
    from genomicode import filelib

    rseqc_path = filelib.which_assert(config.rseqc_path)
    filename = os.path.join(rseqc_path, name)
    filelib.assert_exists_nz(filename)
    return filename


def gtf_to_bed(gtf_filename, bed_filename):
    # Convert a GTF file to a BED file.  Currently only exports exons.
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    filelib.assert_exists_nz(gtf_filename)
    gtfutils = filelib.which_assert(config.gtfutils)

    sq = parallel.quote
    cmd = [
        sq(gtfutils), "tobed",
        "-exons",
        sq(gtf_filename),
        ]
    cmd = " ".join(cmd)
    cmd = "%s > %s" % (cmd, sq(bed_filename))
    parallel.sshell(cmd)

    filelib.assert_exists_nz(bed_filename)



def _get_samtools_SM(filename):
    # @RG
    # ID:Q-Zhang-QZ-Cap475-PM78-NA-Cap475-5983-19_141209-SN746-0319-AC5K...
    # PL:Illumina_HiSeq_2000
    # PU:141209_SN746_0319_AC5KJFACXX
    # LB:Cap475-5983-19
    # SM:Cap475-5983-19
    # CN:IACS-GM

    from genomicode import config
    from genomicode import parallel
    from genomicode import filelib

    samtools = filelib.which_assert(config.samtools)

    # samtools view -H <filename>
    sq = parallel.quote
    cmd = [
        sq(samtools),
        "view",
        "-H",
        sq(filename),
        ]
    cmd = " ".join(cmd)
    x = parallel.sshell(cmd)
    lines = x.split("\n")
    x = [x for x in lines if x.startswith("@RG")]
    assert x, "I could not find the @RG header."
    assert len(x) == 1, "Multiple @RG lines."
    rg_line = x[0]

    x = rg_line.split("\t")
    info = {}
    for x in x:
        if x == "@RG":
            continue
        key, value = x.split(":", 1)
        info[key] = value
    assert "SM" in info, "Misssing SM: %s" % rg_line.strip()

    return info["SM"]


def clean_mutect_vcf(
    normal_bam, cancer_bam, normal_sample, cancer_sample,
    in_filename, out_filename):
    #from genomicode import hashlib

    # MuTect uses the SM field from the BAM file header.  Change this
    # to the name of the samples given by the user.
    # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <normal> <tumor>

    normal_sm = _get_samtools_SM(normal_bam)
    cancer_sm = _get_samtools_SM(cancer_bam)

    # Not sure why hashing this?  Other callers don't hash.
    #normal_sample_h = hashlib.hash_var(normal_sample)
    #cancer_sample_h = hashlib.hash_var(cancer_sample)

    outhandle = open(out_filename, 'w')
    for line in open(in_filename):
        if line.startswith("#CHROM") and line.find(normal_sm) >= 0:
            line = line.replace(normal_sm, normal_sample)
        if line.startswith("#CHROM") and line.find(cancer_sm) >= 0:
            line = line.replace(cancer_sm, cancer_sample)
        outhandle.write(line)
    outhandle.close()


def clean_varscan_vcf(sample, in_filename, out_filename):
    #from genomicode import hashlib

    # Remove lines that contain these values.
    BAD_LINES = [
        "Min coverage:",
        "Min reads2:",
        "Min var freq:",
        "Min avg qual:",
        "P-value thresh:",
        "Reading input from",
        "bases in pileup file",
        "variant positions",
        "were failed by the strand-filter",
        "variant positions reported",
        "Input file was not ready after 100 5-second cycles!",
        "Only SNPs will be reported",
        "Only indels will be reported",
        ("Warning: No p-value threshold provided, so p-values will " 
         "not be calculated"),
        ]

    # Gives error message (and no results) if the machine is out of
    # memory.  Detect this case and raise an exception.
    #
    #There is insufficient memory for the Java Runtime Environment to continue.
    # Cannot create GC thread. Out of system resources.
    # An error report file with more information is saved as:
    # /tmp/jvm-23978/hs_error.log

    # Generate an exception if these lines are found.
    ERROR_LINES = [
        ("There is insufficient memory for the Java Runtime Environment to "
         "continue."),
        "Cannot create GC thread. Out of system resources.",
        "An error report file with more information is saved as:",
        ]

    # Varscan calls each sample "Sample1".  Doesn't have the read
    # group info from the BAM file.  Change this to the proper sample
    # name.
    # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1

    # Not sure why hashing this?  Other callers don't hash.
    #sample_h = hashlib.hash_var(sample)
    
    outhandle = open(out_filename, 'w')
    for line in open(in_filename):
        found = False
        for x in BAD_LINES:
            if line.find(x) >= 0:
                found = True
                break
        if found:
            continue
        for x in ERROR_LINES:
            if line.find(x) >= 0:
                raise AssertionError(line.strip())
        if line.startswith("#CHROM") and line.find("Sample1") >= 0:
            #line = line.replace("Sample1", sample_h)
            line = line.replace("Sample1", sample)
        outhandle.write(line)
    outhandle.close()


def clean_strelka_vcf(normal_sample, cancer_sample, in_filename, out_filename):
    #from genomicode import hashlib

    # Strelka calls the samples "NORMAL" and "TUMOR".  Change this to
    # the proper sample name.
    # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NORMAL TUMOR
    #normal_sample_h = hashlib.hash_var(normal_sample)
    #cancer_sample_h = hashlib.hash_var(cancer_sample)
    # Why need hashing?
    normal_sample_h = normal_sample
    cancer_sample_h = cancer_sample

    outhandle = open(out_filename, 'w')
    for line in open(in_filename):
        if line.startswith("#CHROM") and line.find("NORMAL") >= 0:
            line = line.replace("NORMAL", normal_sample_h)
        if line.startswith("#CHROM") and line.find("TUMOR") >= 0:
            line = line.replace("TUMOR", cancer_sample_h)
        outhandle.write(line)
    outhandle.close()


def clean_muse_vcf(normal_sample, cancer_sample, in_filename, out_filename):
    #from genomicode import hashlib

    # MuSE calls the samples "NORMAL" and "TUMOR".  Change this to
    # the proper sample name.
    # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT TUMOR NORMAL
    #normal_sample_h = hashlib.hash_var(normal_sample)
    #cancer_sample_h = hashlib.hash_var(cancer_sample)
    # Why need hashing?
    normal_sample_h = normal_sample
    cancer_sample_h = cancer_sample

    outhandle = open(out_filename, 'w')
    for line in open(in_filename):
        if line.startswith("#CHROM") and line.find("NORMAL") >= 0:
            line = line.replace("NORMAL", normal_sample_h)
        if line.startswith("#CHROM") and line.find("TUMOR") >= 0:
            line = line.replace("TUMOR", cancer_sample_h)
        outhandle.write(line)
    outhandle.close()


def clean_radia_vcf(normal_sample, cancer_sample, in_filename, out_filename):
    # Rename the samples:
    # DNA_NORMAL   <normal_sample>
    # DNA_TUMOR    <cancer_sample>
    # RNA_TUMOR    <cancer_sample>_RNA

    normal_name = normal_sample
    dna_tumor_name = cancer_sample
    rna_tumor_name = "%s_RNA" % cancer_sample

    # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT DNA_NORMAL
    # DNA_TUMOR RNA_TUMOR
    outhandle = open(out_filename, 'w')
    for line in open(in_filename):
        if line.startswith("#CHROM") and line.find("DNA_NORMAL") >= 0:
            line = line.replace("DNA_NORMAL", normal_name)
        if line.startswith("#CHROM") and line.find("DNA_TUMOR") >= 0:
            line = line.replace("DNA_TUMOR", dna_tumor_name)
        if line.startswith("#CHROM") and line.find("RNA_TUMOR") >= 0:
            line = line.replace("RNA_TUMOR", rna_tumor_name)
        outhandle.write(line)
    outhandle.close()


def get_muse_version():
    import re
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    muse = filelib.which_assert(config.muse)
    x = parallel.sshell(muse, ignore_nonzero_exit=True)
    x = x.strip()
    # Version: v1.0rc
    m = re.search(r"Version: ([\w\.]+)", x)
    assert m, "Missing version string:\n%s" % x
    return m.group(1)
