from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parselib
        from genomicode import parallel
        from Betsy import module_utils as mlib

        in_vcf_node, bf_vcf_node = antecedents
        in_vcf_filenames = filelib.list_files_in_path(
            in_vcf_node.identifier, endswith=".vcf", toplevel_only=True)
        bf_vcf_filenames = filelib.list_files_in_path(
            bf_vcf_node.identifier, endswith=".vcf", toplevel_only=True)
        filelib.safe_mkdir(out_path)
        metadata = {}

        common_only = mlib.get_user_option(
            user_options, "backfill_common_only", allowed_values=["no", "yes"],
            not_empty=True)

        in_vcf_samples = [mlib.splitpath(x)[1] for x in in_vcf_filenames]
        bf_vcf_samples = [mlib.splitpath(x)[1] for x in bf_vcf_filenames]

        # Make sure there are no duplicate sample names.
        x1 = {}.fromkeys(in_vcf_samples).keys()
        x2 = {}.fromkeys(bf_vcf_samples).keys()
        assert len(in_vcf_samples) == len(x1), "Duplicate samples"
        assert len(bf_vcf_samples) == len(x2), "Duplicate samples"

        # Find the samples.
        common = [x for x in in_vcf_samples if x in bf_vcf_samples]
        in_only = [x for x in in_vcf_samples if x not in common]
        bf_only = [x for x in bf_vcf_samples if x not in common]
        assert common, "No common samples."
        
        pretty_in = parselib.pretty_list(in_only, max_items=5)
        pretty_bf = parselib.pretty_list(bf_only, max_items=5)
        if common_only == "no":
            assert not (in_only and bf_only), \
                   "Extra samples in both sets:\n%s\n%s" % (
                pretty_in, pretty_bf)
            assert not in_only, "Target VCF file has extra samples: %s" % \
                   pretty_in
            assert not bf_only, "Source VCF file has extra samples: %s." % \
                   pretty_bf
        SAMPLES = common
        
        # list of sample, in_vcf_filename, bf_vcf_filename, out_filename
        jobs = []
        for sample in SAMPLES:
            assert sample in in_vcf_samples
            assert sample in bf_vcf_samples
            i = in_vcf_samples.index(sample)
            j = bf_vcf_samples.index(sample)
            in_filename = in_vcf_filenames[i]
            bf_filename = bf_vcf_filenames[j]
            out_filename = os.path.join(out_path, "%s.vcf" % sample)
            x = sample, in_filename, bf_filename, out_filename
            jobs.append(x)

        jobs2 = []
        for x in jobs:
            sample, in_filename, bf_filename, out_filename = x
            fn = backfill_vcf
            args = in_filename, bf_filename, out_filename
            keywds = {}
            jobs2.append((fn, args, keywds))
        #num_cores = 1
        parallel.pyfun(jobs2, num_procs=num_cores)
        metadata["num_cores"] = num_cores

        return metadata
        

    def name_outfile(self, antecedents, user_options):
        return "backfilled.vcf"



def backfill_vcf(in_file, bf_file, out_file):
    import copy
    from genomicode import vcflib
    
    #print in_mvcf_node.identifier
    #print back_mvcf_node.identifier
    in_vcf = vcflib.read(in_file)
    bf_vcf = vcflib.read(bf_file)

    # May have multiple samples, e.g. germline and tumor.
    #assert len(in_vcf.samples) == 1, "Too many samples: %s" % in_vcf.samples
    x = [x for x in in_vcf.samples if x in bf_vcf.samples]
    SAMPLES = x

    # Parse out the read counts from the backfill vcf.
    bf_variants = {}   # (sample, chrom, pos) -> ref, alt, Variant, Call
    for i in range(bf_vcf.num_variants()):
        var = vcflib.get_variant(bf_vcf, i)
        for sample in SAMPLES:
            call = vcflib.get_call(var, sample)

            if call.num_ref is None and call.num_alt is None and \
               call.total_reads is None and call.vaf is None:
                continue
            x = sample, var.chrom, var.pos
            assert x not in bf_variants, "Duplicate: %s %s %s" % x
            bf_variants[x] = var.ref, var.alt, var, call

    # Find the variants that can be backfilled.
    # List of (chrom, pos, in_var_num, sample, in_call, bf_var, bf_call)
    matches = []
    for i in range(in_vcf.num_variants()):
        in_var = vcflib.get_variant(in_vcf, i)
        for sample in SAMPLES:
            # Skip if there is no backfill information.
            key = sample, in_var.chrom, in_var.pos
            if key not in bf_variants:
                continue
            bf_ref, bf_alt, bf_var, bf_call = bf_variants[key]
            # Don't worry if the variants match.  Just want a
            # rough estimate of the coverage at this location.
            ## Make sure the variants match.
            ##if not is_same_variants(ref, alt, bf_ref, bf_alt):
            ##    continue
            in_call = vcflib.get_call(in_var, sample)
            x = in_var.chrom, in_var.pos, i, sample, in_call, bf_var, bf_call
            matches.append(x)

    # Update the read counts from annotated VCF file.
    out_vcf = copy.deepcopy(in_vcf)
    add_backfill_genotypes(out_vcf)

    seen = {}
    for x in matches:
        chrom, pos, var_num, sample, in_call, bf_var, bf_call = x
        seen[(sample, chrom, pos)] = 1

        var = vcflib.get_variant(out_vcf, var_num)
        GD = var.sample2genodict[sample]

        mapping = [
            ("BFILL_REF", "num_ref"),
            ("BFILL_ALT", "num_alt"),
            ("BFILL_COV", "total_reads"),
            ("BFILL_VAF", "vaf"),
            ]
        changed = False
        for gt_key, call_attr in mapping:
            x = getattr(bf_call, call_attr)
            if x is None:
                continue
            if type(x) is type([]):  # arbitrarily use max
                x = max(x)
            GD[gt_key] = vcflib._format_vcf_value(x)
            changed = True
        if changed:
            vcflib.set_variant(out_vcf, var_num, var)

    # Add the variants that are in bf_file, but not in in_file.
    for x in bf_variants:
        # sample, chrom, pos = x
        if x in seen:
            continue
        bf_ref, bf_alt, bf_var, bf_call = bf_variants[x]
        # VarScan sets the filter_ to "PASS" for everything.  Get rid
        # of this.
        bf_var.filter_ = "BACKFILL"
        vcflib.add_variant(out_vcf, bf_var)

    vcflib.write(out_file, out_vcf)
    


def add_backfill_genotypes(vcf):
    # Will add genotype columns for backfill in place.
    from genomicode import vcflib
    
    # FORMAT      GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
    # <genotype>  0/1:12:28:28:24:4:14.29%:5.5746E-2:37:36:14:10:3:1

    # Columns to add.
    COLUMNS = ["BFILL_REF", "BFILL_ALT", "BFILL_COV", "BFILL_VAF"]

    for i in range(vcf.num_variants()):
        var = vcflib.get_variant(vcf, i)

        changed = False
        for col in COLUMNS:
            if col not in var.genotype_names:
                var.genotype_names.append(col)
                changed = True
            for genodict in var.sample2genodict.itervalues():
                if col in genodict:
                    continue
                genodict[col] = "."
                changed = True
        if changed:
            vcflib.set_variant(vcf, i, var)


def is_same_variants(ref1, alt1, ref2, alt2):
    alphabet = "ACGT"
    # One of them is mutation, the other is deletion.
    # SEQ1  A G
    # SEQ2  A .
    if ref1 == ref2 and alt1 in alphabet and alt2 == ".":
        return False
    if ref1 == ref2 and alt2 in alphabet and alt1 == ".":
        return False
    # One of them has more alleles than the other.
    # SEQ1  A G
    # SEQ2  A G,T
    alts1 = sorted(alt1.split(","))
    alts2 = sorted(alt2.split(","))
    num_alts_1 = len(alts1)
    num_alts_2 = len(alts2)
    if ref1 == ref2 and num_alts_1 != num_alts_2:
        return False
    # One of them has different alt alleles than the other.
    # SEQ1  C G,A
    # SEQ2  C G,T
    if ref1 == ref2 and alts1 != alts2:
        return False
    # Something weird going on.
    # SEQ1  GC G
    # SEQ2  G  .
    if ref1 != ref2 and alt1 != alt2:
        return False
    return True
