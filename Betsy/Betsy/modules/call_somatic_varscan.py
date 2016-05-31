from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        import shutil
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import parselib
        from Betsy import module_utils as mlib

        mpileup_node, nc_node = antecedents
        mpileup_filenames = filelib.list_files_in_path(
            mpileup_node.identifier, endswith=".pileup")
        assert mpileup_filenames, "No .pileup files."
        nc_match = mlib.read_normal_cancer_file(nc_node.identifier)
        #ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)

        # Figure out whether the purpose is to get coverage.  Change
        # the parameters if it is.
        assert "vartype" in out_attributes
        vartype = out_attributes["vartype"]
        assert vartype in ["snp", "indel"]

        sample2pufile = {}  # sample -> mpileup filename
        for filename in mpileup_filenames:
            path, sample, ext = mlib.splitpath(filename)
            sample2pufile[sample] = filename
        
        # Make sure files exist for all the samples.
        all_samples = []
        for (normal_sample, cancer_sample) in nc_match:
            if normal_sample not in all_samples:
                all_samples.append(normal_sample)
            if cancer_sample not in all_samples:
                all_samples.append(cancer_sample)
        missing = [x for x in all_samples if x not in sample2pufile]
        x = parselib.pretty_list(missing, max_items=5)
        assert not missing, "Missing BAM files for samples: %s" % x

        # list of (sample, normal_pileup, cancer_pileup,
        #          tmp1_normal, tmp1_cancer, tmp2_filename, out_filename)
        opj = os.path.join
        jobs = []
        for (normal_sample, cancer_sample) in nc_match:
            normal_pileup = sample2pufile[normal_sample]
            cancer_pileup = sample2pufile[cancer_sample]
            p, sample, ext = mlib.splitpath(cancer_pileup)
            tmp1_normal = opj(out_path, "%s.normal.tmp1" % sample)
            tmp1_cancer = opj(out_path, "%s.cancer.tmp1" % sample)
            tmp2_filename = opj(out_path, "%s.tmp2" % sample)
            out_filename = opj(out_path, "%s.vcf" % sample)
            x = sample, normal_sample, cancer_sample, \
                normal_pileup, cancer_pileup, \
                tmp1_normal, tmp1_cancer, tmp2_filename, out_filename
            jobs.append(x)

        # VarScan will generate a "Parsing Exception" if there are 0
        # reads in a location.  Will be either "0" or blank.  Filter
        # those lines out.
        sq = parallel.quote
        commands = []
        for x in jobs:
            sample, normal_sample, cancer_sample, \
                    normal_pileup, cancer_pileup, \
                    tmp1_normal, tmp1_cancer, tmp2_filename, out_filename = x
            x1 = "awk -F'\t' '$4 >= 1 {print}' %s > %s" % (
                normal_pileup, tmp1_normal)
            x2 = "awk -F'\t' '$4 >= 1 {print}' %s > %s" % (
                cancer_pileup, tmp1_cancer)
            commands.extend([x1, x2])
        parallel.pshell(commands, max_procs=num_cores)
        x = [x[3] for x in jobs] + [x[4] for x in jobs]
        filelib.assert_exists_nz_many(x)
        

        # java -jar VarScan.jar somatic [normal_pileup] [tumor_pileup]
        #   [output] OPTIONS
        varscan = mlib.findbin("varscan_jar")

        # Use parameters from:
        # Using VarScan 2 for Germline Variant Calling and Somatic
        # Mutation Detection
        
        # Make a list of commands.
        commands = []
        for x in jobs:
            sample, normal_sample, cancer_sample, \
                    normal_pileup, cancer_pileup, \
                    tmp1_normal, tmp1_cancer, tmp2_filename, out_filename = x
            x = [
                "java", "-jar", sq(varscan),
                "somatic",
                sq(tmp1_normal), sq(tmp1_cancer),
                sample,
                "--min-coverage", 10,
                "--min-avg-qual", 15,
                "--min-normal-coverage", 10,
                "--min-tumor-coverage", 10,
                "--min-var-freq", 0.05,
                "--somatic-p-value", 0.05,
                "--output-vcf", 1,
                ]
            x = " ".join(map(str, x))
            x = "%s >& %s" % (x, tmp2_filename)
            commands.append(x)

        parallel.pshell(commands, max_procs=num_cores)
        x = [x[5] for x in jobs]
        filelib.assert_exists_nz_many(x)

        # Copy the final file to the right place.
        for x in jobs:
            sample, normal_sample, cancer_sample, \
                    normal_pileup, cancer_pileup, \
                    tmp1_normal, tmp1_cancer, tmp2_filename, out_filename = x
            # Will be written in current directory.
            varscan_out = "%s.snp.vcf" % sample
            if vartype == "indel":
                varscan_out = "%s.indel.vcf" % sample
            filelib.assert_exists(varscan_out)
            shutil.copy2(varscan_out, out_filename)

        # VarScan names the samples "NORMAL" and "TUMOR".  Replace
        # them with the actual names.
        for x in jobs:
            sample, normal_sample, cancer_sample, \
                    normal_pileup, cancer_pileup, \
                    tmp1_normal, tmp1_cancer, tmp2_filename, out_filename = x
            _fix_normal_cancer_names(
                out_filename, normal_sample, cancer_sample)

            
    def name_outfile(self, antecedents, user_options):
        return "varscan.vcf"


def _fix_normal_cancer_names(filename, normal_sample, cancer_sample):
    # VarScan calls the samples NORMAL, TUMOR.  Fix them.
    from genomicode import hashlib
    
    lines, header_i, samples = _read_vcf(filename)

    normal_sample_h = hashlib.hash_var(normal_sample)
    cancer_sample_h = hashlib.hash_var(cancer_sample)

    header = lines[header_i]
    header1 = header[:-len(samples)]
    header2 = header[-len(samples):]

    assert sorted(header2) == ["NORMAL", "TUMOR"], header2
    for i in range(len(header2)):
        if header2[i] == "NORMAL":
            header2[i] = normal_sample_h
        if header2[i] == "TUMOR":
            header2[i] = cancer_sample_h
    lines[header_i] = header1 + header2

    handle = open(filename, 'w')
    for x in lines:
        print >>handle, "\t".join(x)
    handle.close()


def _read_vcf(filename):
    # Return a tuple of:
    # - a list of lines.  Each line is a list of columns.
    # - the index of the header row (or None)
    # - the sample names
    from genomicode import filelib
    
    lines = [x for x in filelib.read_cols(filename)]
    header_i = None
    for i, cols in enumerate(lines):
        if cols[0] == "#CHROM":
            header_i = i
            break
    assert header_i is not None, "Could not find #CHROM: %s" % filename

    header = lines[header_i]
    x = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
         "INFO", "FORMAT"]
    assert header[:len(x)] == x, "Unknown format: %s" % header
    samples = header[len(x):]
    return lines, header_i, samples

    
