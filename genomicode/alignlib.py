"""Next generation sequencing genome alignment.

Classes:
ReferenceGenome
HTSeqCountOutput

Methods:
create_reference_genome
standardize_reference_genome

make_bowtie1_command
parse_bowtie1_output

make_bowtie2_command
parse_bowtie2_output

make_tophat_command

make_bwa_mem_command
make_bwa_aln_command

make_rsem_command
find_rsem_result_files

make_htseq_count_command
parse_htseq_count_output

"""

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
    # <path>/<name>.[fa|fasta].dict  dict file
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
    # Filter out the known RSEM index files.
    x = [x for x in x if not x.lower().endswith(".idx.fa")]
    x = [x for x in x if not x.lower().endswith(".transcripts.fa")]
    fasta_files = x
    # Find the names of the fasta files.
    x = [os.path.splitext(x)[0] for x in fasta_files]
    names = x
    assert fasta_files, "No fasta files found."
    if name is None:
        uniq_names = {}.fromkeys(names).keys()
        assert len(uniq_names) == 1, "Multiple fasta files found."
        name = names[0]
    assert name in names, "Reference genome not found: %s" % name
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
    all_bwa = [
        "%s.amb" % fasta_file, "%s.ann" % fasta_file, "%s.bwt" % fasta_file,
        "%s.pac" % fasta_file, "%s.sa" % fasta_file,
        ]
    all_rsem = [
        "%s.chrlist" % name, "%s.grp" % name, "%s.idx.fa" % name,
        "%s.n2g.idx.fa" % name, "%s.seq" % name, "%s.ti" % name,
        "%s.transcripts.fa" % name,
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
        elif file_ == "%s.dict" % fasta_file:
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


def create_reference_genome(file_or_path, name=None):
    import os

    # If file_or_path is a file, then find the path of this file.
    path = file_or_path
    if not os.path.isdir(path):
        p, f = os.path.split(file_or_path)
        path = p
    return _create_reference_genome_path(path, name=name)


def standardize_reference_genome(
    in_file_or_path, out_path, use_symlinks=False):
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


def make_bowtie1_command(
    reference_fa, sam_file, fastq_file1, fastq_file2=None,
    orientation=None, num_threads=None):
    # reference_fa is the full path to the fasta file.
    # Orientation must be None, "ff", "fr", "rf"
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

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

    sq = shell.quote
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


def make_bowtie2_command(
    reference_fa, fastq_file1, fastq_file2=None, orientation=None,
    sam_file=None, num_threads=None, skip_reads=None, max_reads=None):
    # Orientation must be None, "ff", "fr", "rf"
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

    assert os.path.exists(fastq_file1)
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

    sq = shell.quote
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
    return " ".join(cmd)


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


def make_tophat_command(
    reference_fa, out_path, fastq_file1, fastq_file2=None,
    gtf_file=None, transcriptome_fa=None, orientation=None,
    num_threads=None):
    # reference_fa is the full path to the fasta file.
    # transcriptome_fa should be tophat indexed fasta file of
    # transcriptome.
    # Orientation must be None, "ff", "fr", "rf"
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell
    
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
    if orientation:
        assert orientation in ["ff", "fr", "rf"]
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100

    tophat = filelib.which_assert(config.tophat)

    assert orientation != "ff", "Not handled."
    orientation2ltype = {
        "fr" : "fr-firststrand",
        "rf" : "fr-secondstrand",
        }

    # tophat [options]* <stem> <reads_1.fq> [<reads_2.fa>]
    # --bowtie1  Use bowtie1 instead of bowtie2.
    # -o <outdir>
    # -r <inner dict>
    # -p <num_threads>
    # --library-type fr-unstranded, fr-firststrand, fr-secondstrand

    sq = shell.quote
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
    if orientation:
        ltype = orientation2ltype[orientation]
        cmd += ["--library-type", ltype]
        
    stem = find_reference_stem(reference_fa)
    cmd += [sq(stem)]
    cmd += [sq(fastq_file1)]
    if fastq_file2:
        cmd += [sq(fastq_file2)]
    return " ".join(cmd)
    

def make_bwa_mem_command(
    reference_fa, sam_filename, err_filename, fastq_file1,
    fastq_file2=None, num_threads=None):
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

    assert os.path.exists(reference_fa)
    assert os.path.exists(fastq_file1)
    if fastq_file2:
        assert os.path.exists(fastq_file2)
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100

    bwa = filelib.which_assert(config.bwa)

    # bwa mem -t <num_cores> ref.fa read1.fq read2.fq > aln-pe.sam
    sq = shell.quote
    cmd = [sq(bwa), "mem"]
    if num_threads:
        cmd += ["-t", str(num_threads)]
    cmd += [sq(reference_fa)]
    cmd += [sq(fastq_file1)]
    if fastq_file2:
        cmd += [sq(fastq_file2)]
    cmd += [
        "1>", sq(sam_filename),
        "2>", sq(err_filename),
        ]
    return " ".join(cmd)


def make_bwa_aln_command(
    reference_fa, fastq_filename, sai_filename, err_filename,
    num_threads=None):
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

    assert os.path.exists(reference_fa)
    assert os.path.exists(fastq_filename)
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100
    bwa = filelib.which_assert(config.bwa)

    # bwa aln -t <num_cores> <reference.fa> <input.fastq> > <output.sai>
    sq = shell.quote
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


def make_rsem_command(
    reference_fa, sample_name, fastq_file1, fastq_file2=None,
    num_threads=None):
    # <sample_name> is prefix for all output files,
    #   e.g. <sample_name>.genes.results
    # Orientation must be None, "ff", "fr", "rf"
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

    # rsem-calculate-expression -p <num_cores> --output-genome-bam
    #   --paired-end <file1.fastq> <file2.fastq>
    #   <index_stem> <sample_name> >& <sample_name>.log"
    # For strand-specific:
    #   --forward-prob 1.0
    #   --forward-prob 0.0

    assert os.path.exists(fastq_file1)
    if fastq_file2:
        assert os.path.exists(fastq_file2)
    if num_threads is not None:
        assert num_threads >= 1 and num_threads < 100

    rsem_calculate = filelib.which_assert(config.rsem_calculate)

    sq = shell.quote
    cmd = [
        sq(rsem_calculate),
        ]
    if num_threads:
        cmd += ["-p", str(num_threads)]
    cmd += ["--output-genome-bam"]
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
    x = filelib.list_files_in_path(search_path)
    x1 = [x for x in x if x.endswith(".genes.results")]
    x2 = [x for x in x if x.endswith(".isoforms.results")]
    sample2files = {}  # sample -> gene_filename, isoform_filename
    for gene_filename in x1:
        p, f = os.path.split(gene_filename)
        sample = f.replace(".genes.results", "")
        assert sample not in sample2files
        sample2files[sample] = gene_filename, None
    for isoform_filename in x2:
        p, f = os.path.split(gene_filename)
        sample = f.replace(".genes.results", "")
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


def make_htseq_count_command(
    bam_file, gtf_file, sort_order, stranded, mode=None):
    # sort_order should be "name" or "pos".
    # stranded should be "yes", "no", "reverse".
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

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

    sq = shell.quote
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
    # Return an HTSeqCountResults object.
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
