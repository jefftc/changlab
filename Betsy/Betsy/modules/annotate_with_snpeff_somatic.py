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
        from genomicode import vcflib
        from Betsy import module_utils as mlib

        vcf_node, nc_node = antecedents
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf")
        assert vcf_filenames, "No .vcf files."
        nc_match = mlib.read_normal_cancer_file(nc_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        # Filenames:
        # <caller>.vcf

        wgs_or_wes = mlib.get_user_option(
            user_options, "wgs_or_wes", not_empty=True,
            allowed_values=["wgs", "wes"])
        genome = mlib.get_user_option(
            user_options, "snpeff_genome", not_empty=True)
        databases = list_snpeff_databases()
        assert genome in databases, "Unknown genome database: %s" % genome

        # For each caller, do the SnpEFF calls.  Some callers include
        # the somatic information, others do not.  If germline samples
        # are present, then do with _cancer.  Otherwise, do not.

        # java -Xmx16g -jar $SNPEFF -v -cancer -cancerSamples vcf03.txt
        #   GRCh37.75 vcf02.txt 1> test03.txt 2> test03.log

        # Don't bother annotating positions that do not pass filter.
        # Filter them out first based on FILTER column.

        opj = os.path.join
        jobs = []
        for in_filename in vcf_filenames:
            path, stem, ext = mlib.splitpath(in_filename)
            samples_file = opj(out_path, "%s.cancerSamples.txt" % stem)
            filtered_filename = opj(out_path, "%s.filtered_input" % stem)
            out_filename = opj(out_path, "%s.vcf" % stem)
            log_filename = opj(out_path, "%s.log" % stem)
            x = filelib.GenericObject(
                in_filename=in_filename,
                samples_file=samples_file,
                filtered_filename=filtered_filename,
                out_filename=out_filename,
                log_filename=log_filename)
            jobs.append(x)

        # First, filter each of the VCF files.
        commands = []
        for j in jobs:
            args = j.in_filename, j.filtered_filename, wgs_or_wes
            x = vcflib.filter_vcf_file, args, {}
            commands.append(x)
        parallel.pyfun(commands, num_procs=num_cores)

        # Make the cancer_samples files.
        for j in jobs:
            # Will generate this if there are cancer samples.
            make_cancer_samples_file(
                j.in_filename, nc_match, j.samples_file)

        # Make a list of commands.
        commands = []
        for j in jobs:
            cancer = False
            if os.path.exists(j.samples_file):
                cancer = True
            x = make_snpeff_command(
                j.filtered_filename, genome, j.out_filename, j.log_filename,
                is_cancer=cancer, cancer_samples_file=j.samples_file)
            commands.append(x)

        nc = mlib.calc_max_procs_from_ram(16, upper_max=num_cores)
        parallel.pshell(commands, max_procs=nc)
        metadata["commands"] = commands
        metadata["num_cores"] = nc
            
        # Make sure the analysis completed successfully.
        x = [x.out_filename for x in jobs]
        filelib.assert_exists_nz_many(x)

        # Log files should be empty.
        for j in jobs:
            filelib.assert_exists(j.log_filename)
            assert not filelib.exists_nz(j.log_filename), \
                   "Error with %s.\n%s" % (j.stem, j.log_filename)
            filelib.safe_unlink(j.log_filename)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "snpeff_annotated.vcf"


def list_snpeff_databases():
    import os
    import StringIO
    from genomicode import parallel
    from genomicode import filelib
    from Betsy import module_utils as mlib
    
    path = mlib.get_config("snp_eff_path", which_assert_file=True)
    snpeff = os.path.join(path, "snpEff.jar")
    filelib.assert_exists_nz(snpeff)

    # Genome    Organism    Status    Bundle    Database download link
    # ------    --------    ------    ------    ----------------------
    sq = parallel.quote
    cmd = [
        "java", "-Xmx16g", "-jar", sq(snpeff),
        "databases",
        ]
    output = parallel.sshell(cmd)
    header = i_db = None
    databases = []
    for cols in filelib.read_cols(StringIO.StringIO(output)):
        cols = [x.strip() for x in cols]
        if header is None:
            header = cols
            assert "Genome" in header
            i_db = header.index("Genome")
            continue
        assert len(cols) == len(header)
        if cols[0].startswith("---"):
            continue
        db_name = cols[i_db]
        databases.append(db_name)
    return databases


def make_cancer_samples_file(vcf_file, nc_match, outfile):
    # Two column tab-delimited text.  No headers.
    # <germline>  <tumor>
    from genomicode import vcflib
    from genomicode import jmath

    # vcf samples (joined with bcftools).
    # PIM005_G   peak1   2:PIM001_G      peak2   3:PIM001_G   [...]

    germline_samples = [x[0] for x in nc_match]
    tumor_samples = [x[1] for x in nc_match]

    # Hopefully should be able to find the samples in the first 1000
    # rows.
    vcf = vcflib.read(vcf_file, nrows=1000)

    # Get the samples from the VCF file.
    samples = vcf.samples

    # HACK: Radia has calls from RNA.  Ignore them.
    # <tumor_sample>_RNA
    rna = {}.fromkeys(["%s_RNA" % x for x in tumor_samples])
    samples = [x for x in samples if x not in rna]
    
    # Clean up samples.
    clean = []   # list of tuples ("G" or "T", sample_name)
    for sample in samples:
        if sample in germline_samples:
            x = "G", sample
        elif sample in tumor_samples:
            x = "T", sample
        else:
            # <num>:<germline sample name>
            x = sample.split(":", 1)
            assert len(x) == 2, "Unknown sample name: %s" % sample
            assert jmath.is_int(x[0]), "Unknown sample name: %s" % sample
            s = x[1]
            assert s in germline_samples, "Unknown sample name: %s" % sample
            x = "G", s
        clean.append(x)
    samples = clean

    # If there are no germline samples, then don't make a file.
    x1 = [x for x in samples if x[0] == "G"]
    x2 = [x for x in samples if x[0] == "T"]
    if not x1:
        return None
    # Make sure there are the same number of germline samples.
    assert len(x1) == len(x2), "Germline/Tumor mismatch: %s" % vcf_file
    assert len(samples) % 2 == 0

    # Pairs should contain one "G" and one "T".
    for i in range(0, len(samples), 2):
        t1, s1 = samples[i]
        t2, s2 = samples[i+1]
        assert t1 != t2, "Bad Germline/Tumor ordering: %s" % vcf_file

    lines = []
    for i in range(0, len(samples), 2):
        t1, s1 = samples[i]
        t2, s2 = samples[i+1]
        # Want germline, then tumor.
        if t1 == "T" and t2 == "G":
            t1, s1, t2, s2 = t2, s2, t1, s1
        assert t1 == "G" and t2 == "T"
        x = "%s\t%s\n" % (s1, s2)
        lines.append(x)
    open(outfile, 'w').writelines(lines)


def make_snpeff_command(in_file, genome, out_file, log_file, is_cancer=False,
                        cancer_samples_file=None):
    import os
    from genomicode import filelib
    from genomicode import parallel
    from Betsy import module_utils as mlib
    
    if is_cancer:
        filelib.assert_exists_nz(cancer_samples_file)

    path = mlib.get_config("snp_eff_path", which_assert_file=True)
    snpeff = os.path.join(path, "snpEff.jar")
    filelib.assert_exists_nz(snpeff)

    sq = parallel.quote
    cmd = [
        "java", "-Xmx16g",
        "-jar", sq(snpeff),
        ]
    if is_cancer:
        cmd += [
            "-cancer",
            "-cancerSamples", sq(cancer_samples_file),
            ]
    cmd += [
        sq(genome),
        sq(in_file),
        ]
    cmd = " ".join(cmd)
    cmd = "%s 1> %s 2> %s" % (cmd, sq(out_file), sq(log_file))
    return cmd
