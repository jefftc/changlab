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
        import shutil
        from genomicode import parallel

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
        handle.close()

        # Write out data from each of the VCF files.
        jobs = []
        for x in vcf_files:
            caller, filestem, filename = x

            args = filename, caller, filestem, header, TEMPFILE
            x = summarize_vcf_file, args, {}
            jobs.append(x)
        parallel.pyfun(jobs, num_procs=num_cores, lock_keyword="lock")
        metadata["num_cores"] = num_cores

        shutil.move(TEMPFILE, out_filename)

        return metadata
            
    
    def name_outfile(self, antecedents, user_options):
        return "variants.txt"


def summarize_vcf_file(filename, caller, filestem, header, outfilename, lock):
    from genomicode import hashlib
    from genomicode import vcflib

    CALLER2FILTERFN = {
        "MuTect" : get_filter_mutect,
        "VarScan2" : get_filter_varscan,
        "Strelka" : get_filter_strelka,
        "SomaticSniper" : get_filter_somaticsniper,
        "mutationSeq" : get_filter_jointsnvmix,
        }

    vcf = vcflib.read(filename)

    lines = []
    for i in range(vcf.num_variants()):
        var = vcflib.get_variant(vcf, i)
        ref = var.ref
        if type(ref) is not type(""):
            ref = ",".join(ref)
        alt = var.alt
        if type(alt) is not type(""):
            alt = ",".join(alt)

        # Some callers may need more processing.
        assert caller in CALLER2FILTERFN, "Unknown caller: %s" % caller
        filter_fn = CALLER2FILTERFN[caller]
        filt = filter_fn(var)

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

            num_ref = vcflib._format_vcf_value(
                call.num_ref, char_for_None="")
            num_alt = vcflib._format_vcf_value(
                call.num_alt, char_for_None="")
            total_reads = vcflib._format_vcf_value(
                call.total_reads, char_for_None="")
            vaf = vcflib._format_vcf_value(
                call.vaf, char_for_None="")

            call = call.call
            if call is None:
                call = ""

            x = caller, filestem, clean_sample, var.chrom, var.pos, \
                ref, alt, num_ref, num_alt, total_reads, vaf, \
                filt, call, GQ
            assert len(x) == len(header)
            x = "\t".join(map(str, x))
            lines.append(x)

            if len(lines) >= 100000:
                x = "\n".join(lines)+"\n"
                lock.acquire()
                handle = open(outfilename, 'a')
                handle.write(x)
                handle.close()
                lock.release()
                lines = []

    x = "\n".join(lines)+"\n"
    lock.acquire()
    handle = open(outfilename, 'a')
    handle.write(x)
    handle.close()
    lock.release()
    


def find_vcf_files(vcf_path):
    # Return list of (<caller>, <sample>, <filename>).
    
    import os
    from genomicode import filelib
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


def get_filter_mutect(var):
    # Mutect uses "PASS" to indicate a confident somatic variant.
    # INFO
    # SOMATIC  If present, indicate somatic mutation.
    # VT       SNP,INS,DEL
    return ",".join(var.filter_)


def get_filter_varscan(var):
    # Varscan uses "PASS" to mean that the variant could be called
    # correctly, and indicates the somatic variant in the SS info.
    # SS 0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown
    assert "SS" in var.infodict
    SS = var.infodict["SS"]
    assert SS in ["0", "1", "2", "3", "5"]
    if SS == "2":
        return "PASS"
    elif SS == "0":
        return "REFERENCE"
    elif SS == "1":
        return "GERMLINE"
    elif SS == "3":
        return "LOH"
    return "UNKNOWN"


def get_filter_strelka(var):
    # INFO
    # SOMATIC   Present in every line
    
    # FILTER seems to be pretty good.
    # ##FILTER=<ID=DP,Description="Greater than 3.0x chromosomal mean depth
    #   in Normal sample">
    # ##FILTER=<ID=BCNoise,Description="Fraction of basecalls filtered at
    #   this site in either sample is at or above 0.4">
    # ##FILTER=<ID=SpanDel,Description="Fraction of reads crossing site
    #   with spanning deletions in either sample exceeeds 0.75">
    # ##FILTER=<ID=QSS_ref,Description="Normal sample is not homozygous ref
    #   or ssnv Q-score < 15, ie calls with NT!=ref or QSS_NT < 15">
    return ",".join(var.filter_)


def get_filter_somaticsniper(var):
    # FILTER   .
    # Information embedded into FORMAT.
    # ##FORMAT=<ID=SS,Number=1,Type=Integer,
    #   Description="Variant status relative to non-adjacent Normal,
    #   0=wildtype,1=germline,2=somatic,3=LOH,4=unknown">

    # Order of samples should be:
    # <normal> <tumor>
    assert len(var.samples) == 2, "SomaticSniper should have two samples"
    s0, s1 = var.samples
    g0, g1 = var.sample2genodict[s0], var.sample2genodict[s1]
    assert "SS" in g0
    assert "SS" in g1
    ss0 = g0["SS"]
    ss1 = g1["SS"]
    # First sample should be germline.
    assert ss0 in ["0", "1", "2", "3", "4"]
    assert ss1 in ["0", "1", "2", "3", "4"]
    assert ss0 in ["0", "1"], "Invalid SS for germline: %s" % ss0

    if ss1 == "0":
        return "WILDTYPE"
    elif ss1 == "1":
        return "GERMLINE"
    elif ss1 == "2":
        return "PASS"
    elif ss1 == "3":
        return "LOH"
    return "UNKNOWN"


def get_filter_jointsnvmix(var):
    # Only PASS calls.  Only shows tumor data.  No information for
    # germline.
    return ",".join(var.filter_)


## def filter_mutect(vcf):
##     # FILTER:
##     # PASS
##     # REJECT
##     from genomicode import AnnotationMatrix
##     from genomicode import vcflib

##     # Keep FILTER==PASS
##     I = []
##     for i in range(vcf.num_variants()):
##         var = vcflib.get_variant(vcf, i)
##         assert len(var.filter_) == 1
##         filt = var.filter_[0]
##         assert filt in ["PASS", "REJECT"]
##         if filt == "PASS":
##             I.append(i)

##     x = AnnotationMatrix.rowslice(vcf.matrix, I)
##     return vcflib.VCFFile(x)


## def filter_varscan(vcf):
##     # FILTER:
##     # PASS
##     # indelError
##     from genomicode import AnnotationMatrix
##     from genomicode import vcflib

##     # Keep FILTER==PASS
##     I = []
##     for i in range(vcf.num_variants()):
##         var = vcflib.get_variant(vcf, i)
##         assert len(var.filter_) == 1
##         filt = var.filter_[0]
##         assert filt in ["PASS", "indelError"]
##         if filt == "PASS":
##             I.append(i)

##     x = AnnotationMatrix.rowslice(vcf.matrix, I)
##     return vcflib.VCFFile(x)


## def filter_strelka(vcf):
##     # FILTER:
##     # PASS
##     # BCNoise
##     # BCNoise;DP
##     # BCNoise;QSS_ref
##     # BCNoise;QSS_ref;DP
##     # BCNoise;QSS_ref;SpanDel
##     # DP
##     # DP;SpanDel
##     # QSS_ref
##     # QSS_ref;DP
##     # QSS_ref;DP;SpanDel
##     # QSS_ref;SpanDel
##     from genomicode import AnnotationMatrix
##     from genomicode import vcflib

##     # Keep FILTER==PASS
##     I = []
##     for i in range(vcf.num_variants()):
##         var = vcflib.get_variant(vcf, i)
##         if len(var.filter_) > 1:
##             continue
##         filt = var.filter_[0]
##         if filt == "PASS":
##             I.append(i)

##     x = AnnotationMatrix.rowslice(vcf.matrix, I)
##     return vcflib.VCFFile(x)


## def filter_somaticsniper(vcf):
##     # FILTER:
##     # .

##     # Assume filtered already.
##     return vcf


## def filter_jointsnvmix(vcf):
##     # FILTER:
##     # PASS

##     # Assume filtered already.
##     return vcf


## def clean_vcf_value(value):
##     # In VCF files, empty values are typically repesented as ".".
##     # Convert these to empty strings.
##     value_list = value.split(",")
##     for i in range(len(value_list)):
##         if value_list[i] == ".":
##             value_list[i] = ""
##     value = ",".join(value_list)
##     return value

    
