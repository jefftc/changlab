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

        # For debugging.
        RUN_VARIANT_CALLING = True
        FILTER_CALLS = True
        MERGE_CALLS = True
        FIX_VCF_FILES = True

        dna_bam_node, rna_bam_node, nc_node, ref_node = antecedents
        dna_bam_filenames = mlib.find_bam_files(dna_bam_node.identifier)
        assert dna_bam_filenames, "No DNA .bam files."
        rna_bam_filenames = mlib.find_bam_files(rna_bam_node.identifier)
        assert rna_bam_filenames, "No RNA .bam files."
        nc_match = mlib.read_normal_cancer_file(nc_node.identifier)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}
        metadata["tool"] = "Radia %s" % alignlib.get_radia_version()

        ## Make sure the BAM files do not contain spaces in the
        ## filenames.  Radia doesn't work well with spaces.
        #filenames = dna_bam_filenames + rna_bam_filenames
        #has_spaces = []
        #for filename in filenames:
        #    if filename.find(" ") >= 0:
        #        has_spaces.append(filename)
        #x = has_spaces
        #if len(x) > 5:
        #    x = x[:5] + ["..."]
        #x = ", ".join(x)
        #msg = "Radia breaks if there are spaces in filenames: %s" % x
        #assert not has_spaces, msg

        # sample -> bam filename
        dnasample2bamfile = mlib.root2filename(dna_bam_filenames)
        rnasample2bamfile = mlib.root2filename(rna_bam_filenames)
        # Make sure files exist for all the samples.  The DNA-Seq
        # should have both normal and cancer.  RNA is not needed for
        # normal sample.
        mlib.assert_normal_cancer_samples(nc_match, dnasample2bamfile)
        mlib.assert_normal_cancer_samples(
            nc_match, rnasample2bamfile, ignore_normal_sample=True)

        # Make sure Radia and snpEff are configured.
        radia_genome_assembly = mlib.get_user_option(
            user_options, "radia_genome_assembly", not_empty=True)
        assert radia_genome_assembly == "hg19", "Only hg19 handled."
        snp_eff_genome = mlib.get_user_option(
            user_options, "snp_eff_genome", not_empty=True)

        radia_path = mlib.get_config("radia_path", assert_exists=True)
        snp_eff_path = mlib.get_config("snp_eff_path", assert_exists=True)
        radia_files = get_radia_files(radia_path, radia_genome_assembly)

        # Make a list of the chromosomes to use.  Pick an arbitrarily
        # BAM file.  Look at only the chromosomes that are present in
        # all files.
        all_bamfiles = dnasample2bamfile.values() + rnasample2bamfile.values()
        chroms = list_common_chromosomes(all_bamfiles)
        assert chroms, "No chromosomes found in all files."
        # Only use the chromosomes that can be filtered by Radia.
        chroms = filter_radia_chromosomes(chroms, radia_files)

        # Make output directories.
        radia_outpath = "radia1.tmp"
        filter_outpath = "radia2.tmp"
        merge_outpath = "radia3.tmp"

        if not os.path.exists(radia_outpath):
            os.mkdir(radia_outpath)
        if not os.path.exists(filter_outpath):
            os.mkdir(filter_outpath)
        if not os.path.exists(merge_outpath):
            os.mkdir(merge_outpath)

        
        # Steps:
        # 1.  Call variants (radia.py)
        #     -o <file.vcf>
        # 2.  Filter variants (filterRadia.py)
        #     <outpath>
        #     Creates a file: <filter_outpath>/<patient_id>_chr<chrom>.vcf
        # 3.  Merge (mergeChroms.py)
        #     Takes as input: <filter_outpath>
        #     Produces: <merge_outpath>/<patient_id>.vcf

        # list of (normal_sample, cancer_sample, chrom, 
        #   normal_bamfile, dna_tumor_bamfile, rna_tumor_bamfile,
        #   radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile,
        #   final_vcf_outfile,
        #   radia_logfile, filter_logfile, merge_logfile)
        opj = os.path.join
        jobs = []
        for i, (normal_sample, cancer_sample) in enumerate(nc_match):
            normal_bamfile = dnasample2bamfile[normal_sample]
            dna_tumor_bamfile = dnasample2bamfile[cancer_sample]
            rna_tumor_bamfile = rnasample2bamfile[cancer_sample]

            merge_vcf_outfile = opj(merge_outpath, "%s.vcf" % cancer_sample)
            merge_logfile = opj(merge_outpath, "%s.log" % cancer_sample)
            final_vcf_outfile = opj(out_path, "%s.vcf" % cancer_sample)
            
            for chrom in chroms:
                radia_vcf_outfile = opj(
                    radia_outpath, "%s_chr%s.vcf" % (cancer_sample, chrom))
                filter_vcf_outfile = opj(
                    filter_outpath, "%s_chr%s.vcf" % (cancer_sample, chrom))
                radia_logfile = opj(
                    radia_outpath, "%s_chr%s.log" % (cancer_sample, chrom))
                filter_logfile = opj(
                    filter_outpath, "%s_chr%s.log" % (cancer_sample, chrom))
                x = normal_sample, cancer_sample, chrom, \
                    normal_bamfile, dna_tumor_bamfile, rna_tumor_bamfile, \
                    radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile, \
                    final_vcf_outfile, \
                    radia_logfile, filter_logfile, merge_logfile
                jobs.append(x)

        # Since Radia doesn't work well if there are spaces in the
        # filenames, symlink these files here to guarantee that there
        # are no spaces.
        normal_path = "normal.bam"
        dna_path = "dna.bam"
        rna_path = "rna.bam"
        if not os.path.exists(normal_path):
            os.mkdir(normal_path)
        if not os.path.exists(dna_path):
            os.mkdir(dna_path)
        if not os.path.exists(rna_path):
            os.mkdir(rna_path)
        for i, x in enumerate(jobs):
            normal_sample, cancer_sample, chrom, \
                normal_bamfile, dna_tumor_bamfile, rna_tumor_bamfile, \
                radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile, \
                final_vcf_outfile, \
                radia_logfile, filter_logfile, merge_logfile = x
            x1 = hash_and_symlink_bamfile(normal_bamfile, normal_path)
            x2 = hash_and_symlink_bamfile(dna_tumor_bamfile, dna_path)
            x3 = hash_and_symlink_bamfile(rna_tumor_bamfile, rna_path)
            clean_normal, clean_dna, clean_rna = x1, x2, x3
            x = normal_sample, cancer_sample, chrom, \
                clean_normal, clean_dna, clean_rna, \
                radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile, \
                final_vcf_outfile, \
                radia_logfile, filter_logfile, merge_logfile
            jobs[i] = x
            

        # Generate the commands for doing variant calling.
        python = mlib.get_config("python", which_assert_file=True)

        # filterRadia.py calls the "blat" command, and there's no way
        # to set the path.  Make sure "blat" is executable.
        if not filelib.which("blat"):
            # Find "blat" in the configuration and add it to the path.
            x = mlib.get_config("blat", which_assert_file=True)
            path, x = os.path.split(x)
            if os.environ["PATH"]:
                path = "%s:%s" % (os.environ["PATH"], path)
            os.environ["PATH"] = path
            # Make sure it's findable now.
            filelib.which_assert("blat")
        

        # STEP 1.  Call variants with radia.py.
        # python radia.py test31 5 \
        # -n bam04/PIM001_G.bam \
        # -t bam04/196B-MG.bam \
        # -r bam34/196B-MG.bam \
        # -f genomes/Broad.hg19/Homo_sapiens_assembly19.fa \
        # -o test32.vcf
        # --dnaTumorMitochon MT \
        # --rnaTumorMitochon MT \
        sq = mlib.sq
        commands = []
        for x in jobs:
            normal_sample, cancer_sample, chrom, \
                normal_bamfile, dna_tumor_bamfile, rna_tumor_bamfile, \
                radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile, \
                final_vcf_outfile, \
                radia_logfile, filter_logfile, merge_logfile = x

            x = [
                sq(python),
                sq(radia_files.radia_py),
                cancer_sample,
                chrom,
                "-n", sq(normal_bamfile),
                "-t", sq(dna_tumor_bamfile),
                "-r", sq(rna_tumor_bamfile),
                "-f", sq(ref.fasta_file_full),
                "-o", radia_vcf_outfile,
                ]
            if "MT" in chroms:
                x += [
                    "--dnaNormalMitochon MT",
                    "--dnaTumorMitochon MT",
                    "--rnaTumorMitochon MT",
                    ]
            x = " ".join(x)
            x = "%s >& %s" % (x, radia_logfile)
            commands.append(x)
        assert len(commands) == len(jobs)
        # Only uses ~200 Mb of ram.
        if RUN_VARIANT_CALLING:
            parallel.pshell(commands, max_procs=num_cores)
        metadata["num_cores"] = num_cores
        metadata["commands"] = commands

        # Make sure log files are empty.
        logfiles = [x[10] for x in jobs]
        filelib.assert_exists_z_many(logfiles)



        # STEP 2.  Filter variants with filterRadia.py.
        commands = []
        for x in jobs:
            normal_sample, cancer_sample, chrom, \
                normal_bamfile, dna_tumor_bamfile, rna_tumor_bamfile, \
                radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile, \
                final_vcf_outfile, \
                radia_logfile, filter_logfile, merge_logfile = x

            x = [
                sq(python),
                sq(radia_files.filterRadia_py),
                cancer_sample,
                chrom,
                sq(radia_vcf_outfile),
                sq(filter_outpath),
                sq(radia_files.scripts_dir),
                "-b", sq(radia_files.blacklist_dir),
                "-d", sq(radia_files.snp_dir),
                "-r", sq(radia_files.retro_dir),
                "-p", sq(radia_files.pseudo_dir),
                "-c", sq(radia_files.cosmic_dir),
                "-t", sq(radia_files.target_dir),
                "-s", sq(snp_eff_path),
                "-e", snp_eff_genome,
                "--rnaGeneBlckFile", sq(radia_files.rnageneblck_file),
                "--rnaGeneFamilyBlckFile", sq(
                    radia_files.rnagenefamilyblck_file),
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, filter_logfile)
            commands.append(x)
        assert len(commands) == len(jobs)

        # Sometimes samtools crashes in the middle of a run.  Detect
        # this case, and re-run the analysis if needed.
        assert len(commands) == len(jobs)
        py_commands = []
        for x, cmd in zip(jobs, commands):
            normal_sample, cancer_sample, chrom, \
                normal_bamfile, dna_tumor_bamfile, rna_tumor_bamfile, \
                radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile, \
                final_vcf_outfile, \
                radia_logfile, filter_logfile, merge_logfile = x
            args = cmd, cancer_sample, chrom, filter_logfile
            x = _run_filterRadia_with_restart, args, {}
            py_commands.append(x)
        # Takes ~10 Gb each.
        nc = mlib.calc_max_procs_from_ram(25, upper_max=num_cores)
        if FILTER_CALLS:
            parallel.pyfun(py_commands, num_procs=nc)
        metadata["commands"] += commands

        # Make sure log files are empty.
        logfiles = [x[11] for x in jobs]
        filelib.assert_exists_z_many(logfiles)

        # Make sure filter_vcf_outfile exists.
        outfiles = [x[7] for x in jobs]
        filelib.assert_exists_nz_many(outfiles)


        # STEP 3.  Merge the results.
        commands = []
        for x in jobs:
            normal_sample, cancer_sample, chrom, \
                normal_bamfile, dna_tumor_bamfile, rna_tumor_bamfile, \
                radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile, \
                final_vcf_outfile, \
                radia_logfile, filter_logfile, merge_logfile = x

            # python /usr/local/radia/scripts/mergeChroms.py 196B-MG \
            #   radia2.tmp/ radia3.tmp
            # The "/" after radia2.tmp is important.  If not given,
            # will generate some files with only newlines.

            fo = filter_outpath
            if not fo.endswith("/"):
                fo = "%s/" % fo
            x = [
                sq(python),
                sq(radia_files.mergeChroms_py),
                cancer_sample,
                fo,
                merge_outpath,
                ]
            x = " ".join(x)
            x = "%s >& %s" % (x, merge_logfile)
            commands.append(x)
        assert len(commands) == len(jobs)
        # Since the chromosomes were separated for the previous steps,
        # this will generate one merge for each chromosome.  This is
        # unnecessary, since we only need to merge once per sample.
        # Get rid of duplicates.
        commands = sorted({}.fromkeys(commands))
        if MERGE_CALLS:
            parallel.pshell(commands, max_procs=num_cores)
        metadata["commands"] += commands

        # Make sure log files are empty.
        logfiles = [x[12] for x in jobs]
        logfiles = sorted({}.fromkeys(logfiles))
        filelib.assert_exists_z_many(logfiles)

        # Fix the VCF files.
        commands = []
        for x in jobs:
            normal_sample, cancer_sample, chrom, \
                normal_bamfile, dna_tumor_bamfile, rna_tumor_bamfile, \
                radia_vcf_outfile, filter_vcf_outfile, merge_vcf_outfile, \
                final_vcf_outfile, \
                radia_logfile, filter_logfile, merge_logfile = x
            args = normal_sample, cancer_sample, \
                   merge_vcf_outfile, final_vcf_outfile
            x = alignlib.clean_radia_vcf, args, {}
            commands.append(x)
        if FIX_VCF_FILES:
            parallel.pyfun(commands, num_procs=num_cores)

        # Make sure output VCF files exist.
        x = [x[9] for x in jobs]
        filelib.assert_exists_nz_many(x)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "radia.vcf"


class RadiaFiles:
    def __init__(
        self, radia_py, filterRadia_py, mergeChroms_py,
        scripts_dir,
        blacklist_dir, snp_dir, retro_dir, pseudo_dir, cosmic_dir, target_dir,
        rnageneblck_file, rnagenefamilyblck_file):
        self.radia_py = radia_py
        self.filterRadia_py = filterRadia_py
        self.mergeChroms_py = mergeChroms_py

        self.scripts_dir = scripts_dir
        
        self.blacklist_dir = blacklist_dir
        self.snp_dir = snp_dir
        self.retro_dir = retro_dir
        self.pseudo_dir = pseudo_dir
        self.cosmic_dir = cosmic_dir
        self.target_dir = target_dir

        self.rnageneblck_file = rnageneblck_file
        self.rnagenefamilyblck_file = rnagenefamilyblck_file
    

def get_radia_files(radia_path, assembly):
    import os
    from genomicode import filelib

    opj = os.path.join

    radia_py = opj(radia_path, "scripts", "radia.py")
    filterRadia_py = opj(radia_path, "scripts", "filterRadia.py")
    mergeChroms_py = opj(radia_path, "scripts", "mergeChroms.py")

    # For hg19 only.
    scripts_dir = opj(radia_path, "scripts")
    blacklist_dir = opj(
        radia_path, "data/%s/blacklists/1000Genomes/phase1" % assembly)
    snp_dir = opj(radia_path, "data/%s/snp135" % assembly)
    retro_dir = opj(radia_path, "data/%s/retroGenes" % assembly)
    pseudo_dir = opj(radia_path, "data/%s/pseudoGenes" % assembly)
    cosmic_dir = opj(radia_path, "data/%s/cosmic" % assembly)
    target_dir = opj(radia_path, "data/%s/gaf/2_1" % assembly)

    rnageneblck_file = opj(radia_path, "data/rnaGeneBlacklist.tab")
    rnagenefamilyblck_file = opj(radia_path, "data/rnaGeneFamilyBlacklist.tab")

    files = [
        radia_py,
        filterRadia_py,
        mergeChroms_py,
        rnageneblck_file,
        rnagenefamilyblck_file,
        ]
    paths = [
        scripts_dir,
        blacklist_dir,
        snp_dir,
        retro_dir,
        pseudo_dir,
        cosmic_dir,
        target_dir,
        ]
    filelib.assert_exists_nz_many(files)
    filelib.assert_exists_many(paths)

    x = RadiaFiles(
        radia_py, filterRadia_py, mergeChroms_py,
        scripts_dir,
        blacklist_dir, snp_dir, retro_dir, pseudo_dir, cosmic_dir, target_dir,
        rnageneblck_file, rnagenefamilyblck_file)
    return x


def list_common_chromosomes(bam_filenames):
    from genomicode import alignlib
    
    common_chroms = None
    for filename in bam_filenames:
        x = alignlib.call_samtools_idxstats(filename)
        x = [x[0] for x in x]
        x = [x for x in x if x != "*"]
        chroms = sorted(x)

        if common_chroms is None:
            common_chroms = chroms
        common_chroms = [x for x in common_chroms if x in chroms]
    return common_chroms

def filter_radia_chromosomes(chroms, radia_files):
    # Only want the ones that exists in the black files.
    # Files in blacklist_dir.
    # - <blacklist_dir>/chr1.bed.gz
    import os

    x = os.listdir(radia_files.blacklist_dir)
    x = [x.replace(".gz", "") for x in x]
    x = [x.replace(".bed", "") for x in x]
    blacklist_chroms = x

    good_chroms = []
    for chrom in chroms:
        if chrom in blacklist_chroms:
            good_chroms.append(chrom)
            continue
        if not chrom.startswith("chr"):
            c = "chr%s" % chrom
            if c in blacklist_chroms:
                good_chroms.append(chrom)
            continue
    return good_chroms


def hash_and_symlink_bamfile(bam_filename, out_path):
    import os
    from genomicode import hashlib

    p, f = os.path.split(bam_filename)
    f = hashlib.hash_alnum(f)
    outfile = os.path.join(out_path, f)
    if not os.path.exists(outfile):
        os.symlink(bam_filename, outfile)

    # Also symlink the index.
    index_filename = "%s.bai" % bam_filename
    index_outfile = "%s.bai" % outfile
    if os.path.exists(index_filename) and not os.path.exists(index_outfile):
        os.symlink(index_filename, index_outfile)
    return outfile


def _run_filterRadia_with_restart(cmd, cancer_sample, chrom, logfile):
    # Sometimes samtools crashes in the middle of a run.  Detect this
    # case, and re-run the analysis if needed.
    from genomicode import parallel
    from genomicode import filelib

    num_tries = 0
    while num_tries <= 3:
        num_tries += 1
        parallel.sshell(cmd, ignore_nonzero_exit=True)
        filelib.assert_exists(logfile)
        log = open(logfile).read()
        # Empty logfile means cmd completed successfully.
        if not log.strip():
            break
        # Look for evidence that samtools died.  If this occurs, try again.
        # 06/29/2016 09:57:16 AM  ERROR   The return code of '1' from the
        #   following filter command indicates an error.
        # 06/29/2016 09:57:16 AM  ERROR   Error from /usr/bin/python
        #   /usr/local/radia/scripts/createBlatFile.pyc 196C-lung2
        #   radia2.tmp/196C-lung2_dnaFiltered_chr1.vcf
        #   radia2.tmp/196C-lung2_mpileup_rna_origin_chr1.vcf
        #   -o radia2.tmp/196C-lung2_blatInput_chr1.fa
        #   --allVCFCalls --blatRnaNormalReads --blatRnaTumorReads:
        # <Traceback>
        # [...]
        #   samtoolsCall.kill()
        # [...]
        # OSError: [Errno 3] No such process
        if log.find("samtoolsCall.kill") >= 0 \
               and log.find("No such process") >= 0:
            continue
        # Otherwise, the process failed for some other reason.  Raise
        # an exception.
        raise AssertionError, "Problem filtering: %s %s\n%s" % (
            cancer_sample, chrom, log)
