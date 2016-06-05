from Module import AbstractModule

# Columns:
# Caller
# Sample
# Chrom
# Pos
# Ref
# Alt
# Num Ref
# Num Alt
# Total Reads
# VAF
# Filter
# Call
# GQ

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        #import os
        import shutil
        #from genomicode import filelib
        #from genomicode import parallel
        #from genomicode import hashlib
        from genomicode import vcflib
        from genomicode import hashlib
        #from Betsy import module_utils as mlib

        #CALLER2FILTERFN = {
        #    "mutect" : filter_mutect,
        #    "varscan" : filter_varscan,
        #    "strelka" : filter_strelka,
        #    "somaticsniper" : filter_somaticsniper,
        #    "jointsnvmix" : filter_jointsnvmix,
        #    }

        # # Make sure I can handle each of the callers.
        # x = [x[0] for x in vcf_files]
        # x = [x for x in x if x not in CALLER2FILTERFN]  # not handled
        # assert not x, "Unknown caller(s): %s" % ", ".join(x)

        vcf_folder = in_data
        vcf_files = find_vcf_files(vcf_folder.identifier)
        metadata = {}

        TEMPFILE = "temp.vcf"
        handle = open(TEMPFILE, 'w')
        header = (
            "Caller", "File", "Sample", "Chrom", "Pos",
            "Ref", "Alt", "Num Ref", "Num Alt", "Total Reads", "VAF",
            "Filter", "Call", "GQ")
        print >>handle, "\t".join(header)

        # Write out data from each of the VCF files.w
        for x in vcf_files:
            caller, filestem, filename = x

            vcf = vcflib.read(filename)
            
            for i in range(vcf.num_variants()):
                var = vcflib.get_variant(vcf, i)
                ref = var.ref
                if type(ref) is not type(""):
                    ref = ",".join(ref)
                alt = var.alt
                if type(alt) is not type(""):
                    alt = ",".join(alt)
                filt = ",".join(var.filter_)
                filt = clean_vcf_value(filt)

                for sample in var.samples:
                    # If sample begins with an integer, there may be a
                    # "X" pre-pended to it.  Try to detect this case
                    # and fix it.
                    clean_sample = sample
                    if sample == hashlib.hash_var(filestem):
                        clean_sample = filestem
                        
                    genodict = var.sample2genodict[sample]
                    GQ = genodict.get("GQ", "")
                    if GQ in [None, "."]:
                        GQ = ""
                    call = vcflib.get_call(var, sample)

                    num_ref = vcflib._format_vcf_value(call.num_ref)
                    num_alt = vcflib._format_vcf_value(call.num_alt)
                    total_reads = vcflib._format_vcf_value(call.total_reads)
                    vaf = vcflib._format_vcf_value(call.vaf)

                    num_ref = clean_vcf_value(num_ref)
                    num_alt = clean_vcf_value(num_alt)
                    total_reads = clean_vcf_value(total_reads)
                    vaf = clean_vcf_value(vaf)

                    call = call.call
                    if call is None:
                        call = ""

                    x = caller, filestem, clean_sample, var.chrom, var.pos, \
                        ref, alt, num_ref, num_alt, total_reads, vaf, \
                        filt, call, GQ
                    assert len(x) == len(header)
                    print >>handle, "\t".join(map(str, x))
                    
            
            #filter_fn = CALLER2FILTERFN[caller]
            #vcf_f = filter_fn(vcf)
        handle.close()

        

        ## # Make a list of each of the variant positions.
        ## variants = {}  # (chrom, pos, ref, alt) -> 1
        ## assert len(vcf_files) == len(vcf_filtered)
        ## for j in range(len(vcf_files)):
        ##     caller, sample, filename = vcf_files[j]
        ##     vcf = vcf_filtered[j]
        ##     for i in range(vcf.num_variants()):
        ##         var = vcflib.get_variant(vcf, i)

        ##         assert type(var.ref) is type(""), "Multiple ref: %s %s %s %d" \
        ##                % (caller, sample, var.chrom, var.pos)
        ##         #assert type(var.alt) is type(""), "Multiple alt: %s %s %s %d"\
        ##         #       % (caller, sample, var.chrom, var.pos)

        ##         # Alt may be a list of possible alternatives.
        ##         alt = var.alt
        ##         if type(var.alt) is type([]):
        ##             alt = tuple(alt)
        ##         x = var.chrom, var.pos, var.ref, alt
        ##         variants[x] = 1


        ## # Make sure there are no conflicting positions.
        ## variants = sorted(variants)
        shutil.move(TEMPFILE, out_filename)

        return metadata
            
    
    def name_outfile(self, antecedents, user_options):
        return "variants.txt"


def find_vcf_files(vcf_path):
    # Return list of (<caller>, <sample>, <filename>).
    
    import os
    from genomicode import filelib
    from genomicode import alignlib
    from genomicode import vcflib

    filenames = filelib.list_files_in_path(
        vcf_path, endswith=".vcf", case_insensitive=True)
    
    # Format:
    # <path>/<sample>.vcf
    vcf_files = []
    for filename in filenames:
        p, f = os.path.split(filename)
        sample = os.path.splitext(f)[0]
        caller = vcflib.identify_caller(filename)
        assert caller is not None, "Unknown caller: %s" % filename
        x = caller, sample, filename
        vcf_files.append(x)
    return vcf_files
    ## vcf_files = []
    ## for caller_path in os.listdir(vcf_path):
    ##     assert caller_path.endswith(".vcf")
    ##     caller = caller_path[:-4]
    ##     for sample_file in os.listdir(opj(vcf_path, caller_path)):
    ##         assert sample_file.endswith(".vcf")
    ##         sample = sample_file[:-4]
    ##         filename = opj(vcf_path, caller_path, sample_file)
    ##         x = caller, sample, filename
    ##         vcf_files.append(x)
    ## return vcf_files


def filter_mutect(vcf):
    # FILTER:
    # PASS
    # REJECT
    from genomicode import AnnotationMatrix
    from genomicode import vcflib

    # Keep FILTER==PASS
    I = []
    for i in range(vcf.num_variants()):
        var = vcflib.get_variant(vcf, i)
        assert len(var.filter_) == 1
        filt = var.filter_[0]
        assert filt in ["PASS", "REJECT"]
        if filt == "PASS":
            I.append(i)

    x = AnnotationMatrix.rowslice(vcf.matrix, I)
    return vcflib.VCFFile(x)


def filter_varscan(vcf):
    # FILTER:
    # PASS
    # indelError
    from genomicode import AnnotationMatrix
    from genomicode import vcflib

    # Keep FILTER==PASS
    I = []
    for i in range(vcf.num_variants()):
        var = vcflib.get_variant(vcf, i)
        assert len(var.filter_) == 1
        filt = var.filter_[0]
        assert filt in ["PASS", "indelError"]
        if filt == "PASS":
            I.append(i)

    x = AnnotationMatrix.rowslice(vcf.matrix, I)
    return vcflib.VCFFile(x)


def filter_strelka(vcf):
    # FILTER:
    # PASS
    # BCNoise
    # BCNoise;DP
    # BCNoise;QSS_ref
    # BCNoise;QSS_ref;DP
    # BCNoise;QSS_ref;SpanDel
    # DP
    # DP;SpanDel
    # QSS_ref
    # QSS_ref;DP
    # QSS_ref;DP;SpanDel
    # QSS_ref;SpanDel
    from genomicode import AnnotationMatrix
    from genomicode import vcflib

    # Keep FILTER==PASS
    I = []
    for i in range(vcf.num_variants()):
        var = vcflib.get_variant(vcf, i)
        if len(var.filter_) > 1:
            continue
        filt = var.filter_[0]
        if filt == "PASS":
            I.append(i)

    x = AnnotationMatrix.rowslice(vcf.matrix, I)
    return vcflib.VCFFile(x)


def filter_somaticsniper(vcf):
    # FILTER:
    # .

    # Assume filtered already.
    return vcf


def filter_jointsnvmix(vcf):
    # FILTER:
    # PASS

    # Assume filtered already.
    return vcf

def clean_vcf_value(value):
    # In VCF files, empty values are typically repesented as ".".
    # Convert these to empty strings.
    value_list = value.split(",")
    for i in range(len(value_list)):
        if value_list[i] == ".":
            value_list[i] = ""
    value = ",".join(value_list)
    return value

    
