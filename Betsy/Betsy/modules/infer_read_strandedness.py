from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)
        self.cache = {}

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib
        from genomicode import alignlib
        from Betsy import module_utils as mlib

        bam_node, gene_node = antecedents
        bam_filenames = mlib.find_bam_files(bam_node.identifier)
        gtf_file = gene_node.identifier
        filelib.assert_exists_nz(gtf_file)
        assert bam_filenames, "No bam files found."
        metadata = {}

        # Make output filenames.
        p, r, e = mlib.splitpath(gtf_file)
        bed_file = "%s.bed" % r

        # Make bed file.
        alignlib.gtf_to_bed(gtf_file, bed_file)
        #bed_file = "/data/jchang/biocore/gtf02.txt"

        # Figure out the orientation.
        x = get_paired_stranded_rseqc(bed_file, bam_filenames[0])
        single_or_paired, stranded, frac_failed, frac_first, frac_second = x

        x = mlib.Stranded(
            single_or_paired, stranded, frac_failed, frac_first, frac_second)
        mlib.write_stranded(x, outfile)
        return metadata
    
        
    def name_outfile(self, antecedents, user_options):
        return "stranded.json"


def parse_rseqc_infer_experiment(output):
    # Return tuple of (
    #   "single" or "paired",
    #   "unstranded", "firststrand", "secondstrand"
    #   <fraction failed>,
    #   <fraction firststrand>,
    #   <fraction secondstrand>
    #   )

    # Reading reference gene model gencode.v19.annotation.bed ... Done
    # Loading SAM/BAM file ...  Finished
    # Total 184900 usable reads were sampled
    #
    #
    # This is PairEnd Data
    # Fraction of reads failed to determine: 0.0889
    # Fraction of reads explained by "1++,1--,2+-,2-+": 0.0229
    # Fraction of reads explained by "1+-,1-+,2++,2--": 0.8882

    # This is SingleEnd Data
    # Fraction of reads failed to determine: 0.0170
    # Fraction of reads explained by "++,--": 0.9669
    # Fraction of reads explained by "+-,-+": 0.0161

    # MISSING FILE ERROR:
    # /data/jchang/biocore/gtf02.txt does NOT exists.
    
    x = output.split("\n")
    x = [x.strip() for x in x]
    x = [x for x in x if x]
    lines = x

    # Look for missing file error.
    for x in lines:
        if x.find("does NOT exists") >= 0:
            raise AssertionError, "RSeQC Error: %s" % x.strip()
    
    # Might be correct.  Do more checks.
    assert len(lines) >= 4

    # Look for the "This is" line.
    for i in range(len(lines)):
        if lines[i].startswith("This is"):
            break
    else:
        raise AssertionError, 'Cannot find "This is ..." line:\n%s' % output

    single_or_paired = None
    if lines[i].find("PairEnd") >= 0:
        single_or_paired = "paired"
    elif lines[i].find("SingleEnd") >= 0:
        single_or_paired = "single"
    else:
        raise AssertionError, "Unknown output: %s" % lines[i]

    for i in range(i+1, len(lines)):
        if lines[i].find("failed to determine") >= 0:
            x = lines[i].split()
            frac_failed = float(x[-1])
        elif lines[i].find("1++,1--,2+-,2-+") >= 0:
            assert single_or_paired == "paired"
            x = lines[i].split()
            frac_secondstrand = float(x[-1])
        elif lines[i].find("1+-,1-+,2++,2--") >= 0:
            assert single_or_paired == "paired"
            x = lines[i].split()
            frac_firststrand = float(x[-1])
        elif lines[i].find("++,--") >= 0:
            assert single_or_paired == "single"
            x = lines[i].split()
            frac_secondstrand = float(x[-1])
        elif lines[i].find("+-,-+") >= 0:
            assert single_or_paired == "single"
            x = lines[i].split()
            frac_firststrand = float(x[-1])
        else:
            raise AssertionError, "Unknown line: %s" % lines[i]
    assert frac_failed is not None
    assert frac_firststrand is not None
    assert frac_secondstrand is not None

    stranded = "unstranded"
    if frac_firststrand >= 0.7:
        stranded = "firststrand"
    elif frac_secondstrand >= 0.7:
        stranded = "secondstrand"
    x = single_or_paired, stranded, frac_failed, \
        frac_firststrand, frac_secondstrand
    return x

    
def get_paired_stranded_rseqc(reference_bed, bam_filename):
    from genomicode import alignlib
    from genomicode import filelib
    from genomicode import parallel
    from Betsy import module_utils as mlib

    script = alignlib.find_rseqc_script("infer_experiment.py")
    filelib.assert_exists_nz(reference_bed)
    filelib.assert_exists_nz(bam_filename)

    # RSeQC scripts use #!/usr/bin/python, which may not be the right
    # one.  Use the python on the path.
    cmd = [
        "python",
        mlib.sq(script),
        "-r", mlib.sq(reference_bed),
        "-i", mlib.sq(bam_filename),
        ]
    cmd = " ".join(cmd)
    print "HERE 3", cmd
    x = parallel.sshell(cmd)
    x = parse_rseqc_infer_experiment(x)
    #single_or_paired, stranded, frac_failed, frac_first, frac_second = x
    return x
