"""
Classes:
VCFFile
Variant
Call

Caller
NextGENeCaller
samtoolsCaller
GATKCaller
PlatypusCaller
MuTectCaller
MuTect2Caller
VarScan2Caller
SomaticSniperCaller
StrelkaCaller
JointSNVMixCaller
MuSECaller
RadiaCaller
PindelCaller
CALLERS              Global variable with all Caller classes.

Functions:
read
write

get_caller        Return the caller based on a vcf_filename

get_variant      Get Variant from VCFFile.
set_variant      Set a Variant to a VCFFile.
add_variant      Add a Variant to a VCFFile.
get_call         Get a Call from a Variant.
set_call         Set a Call to a Variant.

is_pass_filter   Whether the variant passes a filter.
select_variants  Select a subset of the variants.

make_coverage_matrix
make_vaf_matrix

simplify_call

"""
# _parse_info
# _parse_genotype
# _parse_vcf_value
# _format_info
# _format_genotype
# _format_vcf_value
#
# has_info
# has_all_info
# has_format
# has_all_format
#
# _safe_int
# _safe_float
# _safe_add
# _percent_to_decimal


# VCF Format:
# INFO        ADP=28;WT=0;HET=1;HOM=0;NC=0
# FORMAT      GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
# <genotype>  0/1:12:28:28:24:4:14.29%:5.5746E-2:37:36:14:10:3:1
# * Sometimes INFO doesn't have an "=".
#   1000g2015aug_afr=.;1000g2015aug_eas=.;1000g2015aug_eur=.;ALLELE_END
# * INFO
#   WT/HET/HOM/NC  number of samples WT/HET/HOM/NC (not called).
# * FORMAT
#   SDP            Raw read depth
#   DP             Quality adjusted read depth
#   ADP            Average depth across samples.
#
# For some rows, samtools skips the values.  If this happens, try
# to guess the values based on the INFO column.
# INFO    ADP=16;WT=1;HET=0;HOM=0;NC=0
# FORMAT  GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
# VALUES  ./.:.:16
#
#                   samtools    Platypus          GATK
# num_ref           GENO:RD                       GENO:AD
# num_alt           GENO:AD     GENO:NV  INFO:TR  GENO:AD
# total_reads       GENO:DP     GENO:NR  INFO:TC  GENO:DP
# vaf               GENO:FREQ
# call              GENO:GT     GENO:GT           GENO:GT
#
#                   bcftools   NextGene              Backfill
# num_ref                      GENO:SGCOUNTREF_F/R   GENO:BFILL_REF
# num_alt                      GENO:SGCOUNTALT_F/R   GENO:BFILL_ALT
# total_reads       INFO:DP    GENO:DP               GENO:BFILL_COV
# vaf                                                GENO:BFILL_VAF
# call              GENO:GT
#
#                   JointSNVMix  SomaticSniper
# num_ref           INFO:TR      GENO:DP4
# num_alt           INFO:TA      GENO:DP4
# total_reads                    GENO:DP
# vaf
# call
#
# * BFILL_ format is used for backfilling reads.
# * For multi VCF files, Platypus INFO:TR and INFO:TC just contains
#   the info for the first sample (with reads).  For positions without
#   calls, GENO:NV and GENO:NR are blank.
# * GENO:AD    For GATK, includes REF.  REF,ALT1,ALT2 etc.
#   GENO:FREQ  For samtools, is a percent, e.g. 100%.
#   GENO:FREQ  May be "." or blank.
# * If DP is missing, then try SDP.
# * Sometimes DP is missing for samtools.  Can use INFO:ADP if
#   only 1 sample.
# * NextGENe
#   Comma means different (simultaneous) interpretations of
#   alternate alignments.  Allele frequency can be either.
#   SGCOUNTREF_F  3277
#   SGCOUNTREF_R  2675
#   SGCOUNTALT_F  1500,522
#   SGCOUNTALT_R  1761,515
#   AF            0.353,0.112
#   SGACOV        3261,1037   Sum of ALT reads.
#   DP            9232
#
#
# JoinSNVMix has a non-standard VCF file.  No "FORMAT" column, and
# does not have information about each sample.  Everything in INFO.



class VCFFile:
    def __init__(self, matrix, caller):
        # matrix is an AnnotationMatrix.
        import copy

        # Should start with headers, and then samples after that.
        vcf_headers = [
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT"]
        # JointSNVMix does not have FORMAT column.  This one is
        # optional.
        x1 = matrix.headers[:len(vcf_headers)]
        x2 = matrix.headers[len(vcf_headers):]
        if len(matrix.headers) == 8:
            x1 = matrix.headers[:8]
            x2 = []
        assert x1 == vcf_headers[:len(x1)], "Non-standard VCF headers: %s" % (
            " ".join(x1))
        samples = x2

        self.matrix = copy.deepcopy(matrix)
        self.samples = samples
        self.caller = caller

    def num_variants(self):
        return self.matrix.num_annots()


class Variant:
    # Holds structured information for each variant.
    #
    # chrom            string
    # pos              int
    # id_              string
    # ref              string
    # alt              list of strings
    # qual             float or None
    # filter_          list of strings
    # info_names       list of strings
    # infodict
    # samples          list of strings
    # genotype_names   list of strings.
    #                  Can be empty list if no FORMAT (jointsnvmix).
    # sample2genodict  Can be empty dict if no FORMAT.
    # caller           Caller
    def __init__(
        self, chrom, pos, id_, ref, alt, qual, filter_,
        info_names, infodict, samples, genotype_names, sample2genodict,
        caller):
        import copy
        assert type(chrom) is type("")
        assert type(pos) is type(0)
        assert type(id_) is type("")
        assert type(ref) is type("")
        assert type(alt) is type([])
        assert type(qual) is type(0.0) or qual is None
        assert type(filter_) is type([])
        assert type(info_names) is type([])
        assert type(infodict) is type({})
        assert type(samples) is type([])
        assert type(genotype_names) is type([])
        assert type(sample2genodict) is type({})

        self.chrom = chrom
        self.pos = pos
        self.id_ = id_
        self.ref = ref
        self.alt = alt[:]
        self.qual = qual
        self.filter_ = filter_[:]
        self.info_names = info_names[:]
        self.infodict = infodict.copy()
        self.samples = samples[:]
        self.genotype_names = genotype_names[:]
        self.sample2genodict = copy.deepcopy(sample2genodict)
        self.caller = caller
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [repr(self.chrom),
             repr(self.pos),
             repr(self.id_),
             repr(self.ref),
             repr(self.alt),
             repr(self.qual),
             repr(self.filter_),
             repr(self.info_names),
             repr(self.infodict),
             repr(self.samples),
             repr(self.genotype_names),
             repr(self.sample2genodict),]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x



class Call:
    # Contains call information for a variant.
    #
    # total_reads      int or None
    # num_ref          int or None
    # num_alt          list of (int or None).  Can be empty list.
    # vaf              list of (float or None).  Can be empty list.
    # call             string or None
    def __init__(self, num_ref, num_alt, total_reads, vaf, call):
        assert total_reads is None or type(total_reads) is type(0), \
               "Invalid: %s" % repr(total_reads)
        assert total_reads is None or total_reads >= 0
        assert num_ref is None or type(num_ref) is type(0)
        assert num_ref is None or num_ref >= 0, "Invalid num_ref: %s" % num_ref
        assert type(num_alt) is type([])
        for x in num_alt:
            assert type(x) is type(0) or x is None
            if type(x) is type(0):
                assert x >= 0
        assert type(vaf) is type([])
        for x in vaf:
            if type(x) is type(0):
                x = float(x)
            assert type(x) is type(0.0) or x is None, \
                   "Invalid vaf: %s" % repr(x)
            if type(x) is type(0.0):
                assert x >= 0
        assert len(num_alt) == len(vaf)
        assert type(call) is type("") or call is None
        
        #assert num_alt_alleles >= 1
        #self.num_alt_alleles = num_alt_alleles
        self.total_reads = total_reads
        self.num_ref = num_ref
        self.num_alt = num_alt[:]
        self.vaf = vaf[:]
        self.call = call
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            repr(self.num_ref),
            repr(self.num_alt),
            repr(self.total_reads),
            repr(self.vaf),
            repr(self.call),
            ]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x
    def __cmp__(self, other):
        if not isinstance(other, Call):
            return cmp(id(self), id(other))
        x1 = [
            self.num_ref, self.num_alt,
            self.total_reads, self.vaf, self.call]
        x2 = [
            other.num_alt_alleles, other.num_ref, other.num_alt,
            other.total_reads, other.vaf, other.call]
        return cmp(x1, x2)


class Caller:
    def __init__(self, name):
        self.name = name
    def is_caller(self, header_lines):
        raise NotImplementedError
    def get_call(self, var, sample):
        # Return a Call object.
        raise NotImplementedError
    def get_filter(self, var):
        # Return a string that can be used to filter on.  Usually the
        # FILTER field, but some formats provide this somewhere else.
        raise NotImplementedError
    def is_pass(self, filter_str):
        # Whether the return value from get_filter indicates a passing
        # variant.
        raise NotImplementedError

        


class NextGENeCaller(Caller):
    # FORMAT:SGCOUNTREF_F   ref/forward
    # FORMAT:SGCOUNTREF_R
    # FORMAT:SGCOUNTALT_F
    # FORMAT:SGCOUNTALT_R
    # FORMAT:DP             total reads
    # FORMAT:GT
    def __init__(self):
        Caller.__init__(self, "NextGENe")
    def is_caller(self, header_lines):
        # ##source=NextGENeV2.4.1.2
        x = [x for x in header_lines if x.startswith("##source=NextGENe")]
        if x:
            return True
        return False
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        total_reads = _parse_vcf_value(
            genodict["DP"], to_int=1, len_exactly=1)[0]
        x1 = _parse_vcf_value(
            genodict["SGCOUNTREF_F"], to_int=1, len_exactly=1)[0]
        x2 = _parse_vcf_value(
            genodict["SGCOUNTREF_R"], to_int=1, len_exactly=1)[0]
        x3 = _parse_vcf_value(genodict["SGCOUNTALT_F"], to_int=1)
        x4 = _parse_vcf_value(genodict["SGCOUNTALT_R"], to_int=1)
        num_ref = x1 + x2
        assert len(x3) == len(x4)
        num_alt = []
        for i in range(len(x3)):
            num_alt.append(x3[i] + x4[i])
        vaf = [x/float(total_reads) for x in num_alt]
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Returns PASS.
        
        assert len(var.filter_) == 1
        x = var.filter_[0]
        assert x in ["PASS"]
        return x
    def is_pass(self, filter_str):
        assert filter_str == "PASS"
        return True



class samtoolsCaller(Caller):
    # INFO:DP    total_reads
    # INFO:DP4   high-quality ref-fwd, ref-rev, alt-fwd, alt-rev
    # FORMAT:GT  call
    def __init__(self):
        Caller.__init__(self, "samtools")
    def is_caller(self, header_lines):
        # ##samtoolsVersion=1.2+htslib-1.2.1
        x = [x for x in header_lines if x.startswith("##samtoolsVersion=")]
        if x:
            return True
        return False
    def get_call(self, var, sample):
        infodict = var.infodict
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        total_reads = _parse_vcf_value(
            infodict["DP"], to_int=1, len_exactly=1)[0]
        dp4 = _parse_vcf_value(infodict["DP4"], to_int=1, len_exactly=4)
        num_ref = dp4[0] + dp4[1]
        num_alt = [dp4[2] + dp4[3]]
        vaf = [num_alt[0] / float(total_reads)]
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Return PASS.
        
        # Always ".".
        assert not var.filter_
        #assert len(var.filter_) == 1
        #x = var.filter_[0]
        #assert x == None
        return "PASS"
    def is_pass(self, filter_str):
        assert filter_str == "PASS"
        return True


class GATKCaller(Caller):
    # FORMAT:AD  num_ref, num_alt
    # FORMAT:DP  total_reads
    # FORMAT:GT  call
    #
    # Also see (rare):
    # GT:GQ:PL         (INFO:DP=0  REF=., ALT=C)
    # 1/1:6:49,6,0
    def __init__(self):
        Caller.__init__(self, "GATK")
    def is_caller(self, header_lines):
        # Older versions:
        # ##GATKCommandLine=<ID=HaplotypeCaller
        # Newer versions:
        # ##GATKCommandLine.HaplotypeCaller=<ID=HaplotypeCaller
        x = [x for x in header_lines if x.startswith("##GATKCommandLine")]
        # BUG: does not handle UnifiedGenotypeCalle
        x = [x for x in x if x.find("HaplotypeCaller") >= 0]
        if x:
            return True
        return False
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]

        if "AD" in genodict and "DP" in genodict:
            x = _parse_vcf_value(genodict["AD"], to_int=1, len_at_least=2)
            num_ref = x[0]
            num_alt = x[1:]
            total_reads = _parse_vcf_value(
                genodict["DP"], to_int=1, len_exactly=1)[0]
            if not total_reads:
                vaf = [0.0] * len(num_alt)
            else:
                vaf = [x/float(total_reads) for x in num_alt]
        else:
            num_ref = 0
            num_alt = [0]
            total_reads = 0
            vaf = [0.0]
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Return PASS.
        
        # Always ".".
        assert not var.filter_
        #assert len(var.filter_) == 1, "Unexpected: %s" % str(var.filter_)
        #x = var.filter_[0]
        #assert x == None, "Unexpected: %s" % x
        return "PASS"
    def is_pass(self, filter_str):
        assert filter_str == "PASS"
        return True


class PlatypusCaller(Caller):
    # INFO:TC    total coverage
    # INFO:TCF
    # INFO:TCR
    # INFO:TR    Number of reads with variant.
    # INFO:NR    reverse reads
    # FORMAT:NV  number of reads containing variant.  num_alt
    # FORMAT:NR  number of reads covering variant location
    # FORMAT:GT  call
    def __init__(self):
        Caller.__init__(self, "Platypus")
    def is_caller(self, header_lines):
        # ##source=Platypus_Version_0.8.1
        x = [x for x in header_lines if x.startswith("##source=Platypus")]
        if x:
            return True
        return False
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]

        # FORMAT:NR FORMAT:NV INFO:TC  INFO:TR
        #   12, 12   12, 12     12
        #    6,  6    2,  3      6
        #  879,879  576,576    879
        #    5,  6    2,  6      5       2, 6
        #        0        0      0          0   Q20;badReads;MQ;QD
        # o May have multiple alt alleles that simultaneously overlap.
        #   Will have one FORMAT:NR for each alt allele.
        # o With indels, it's possible for the number of reads to be
        #   different for different alt alleles.
        #     *AA  ALT1
        #     TAA  ALT2   Has more reads than ALT1
        #   >>>    READ
        #   >>>>>  READ
        # o total_reads should be a list.  For simplicity, just use
        #   the largest one.
        # o Only counts good reads.  May be 0 if reads are bad.
        total_list = _parse_vcf_value(genodict["NR"], to_int=1)
        #assert min(total_list) > 0, "Zero total reads: %s %s %d" % (
        #    sample, var.chrom, var.pos)
        total_reads = max(total_list)
        num_alt = _parse_vcf_value(genodict["NV"], to_int=1)
        assert len(total_list) == len(num_alt)
        x = [total_list[i]-num_alt[i] for i in range(len(num_alt))]
        # Minimum num_ref corresponds to max vaf.
        num_ref = min(x)
        vaf = [_safe_div_z(num_alt[i], total_list[i])
               for i in range(len(num_alt))]
        assert num_ref >= 0, "Negative num ref: %s %s %d" % (
            sample, var.chrom, var.pos)
        assert min(vaf) >= 0, "Negative vaf: %s %s %d" % (
            sample, var.chrom, var.pos)
        assert max(vaf) <= 1, "vaf greater than 1: %s %s %d" % (
            sample, var.chrom, var.pos)

        #total_reads = _parse_vcf_value(
        #    var.infodict["TC"], to_int=1, len_exactly=1)[0]
        #x = _parse_vcf_value(genodict["NV"], to_int=1)
        #num_alt = [max(x)]
        #num_ref = total_reads - sum(num_alt)
        #assert num_ref >= 0, "Negative num ref: %s %s %d" % (
        #    sample, var.chrom, var.pos)
        #if total_reads:
        #    vaf = [x/float(total_reads) for x in num_alt]
        #else:
        #    vaf = [None] * len(num_alt)
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Return:
        # PASS, badReads, alleleBias, Q20, MQ
        # Can be comma-separated combination as well.
        return ",".join(var.filter_)
    def is_pass(self, filter_str):
        return filter_str == "PASS"


class MuTectCaller(Caller):
    # FORMAT:DP   Approximate read depth (MQ=255, bad mates filtered)
    # FORMAT:AD   Depth for ref and alt in order listed
    # FORMAT:FA   alternate allele fraction  (variant allele frequency)
    # FORMAT:GT   Genotype
    def __init__(self):
        Caller.__init__(self, "MuTect")
    def is_caller(self, header_lines):
        # MuTect
        # ##GATKCommandLine=<ID=MuTect,Version=3.1-0-g72492bb,Date="Sun
        # Mutect2
        # ##GATKCommandLine.MuTect2=<ID=MuTect2,Version=3.6-0-g89b7209,Date=
        # Make sure we're picking up MuTect and not MuTect2.
        x = [x for x in header_lines if x.startswith("##GATKCommandLine")]
        x = [x for x in x if x.find("MuTect") >= 0]
        x = [x for x in x if x.find("MuTect,") >= 0]
        if x:
            return True
        return False
    #def is_caller_from_vcf(self, vcf):
    #    if not has_all_format(vcf, ["DP", "AD", "FA", "GT"]):
    #        return False
    #    return True
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        total_reads = _parse_vcf_value(
            genodict["DP"], to_int=1, len_exactly=1)[0]
        x = _parse_vcf_value(genodict["AD"], to_int=1, len_at_least=2)
        num_ref = x[0]
        num_alt = x[1:]
        vaf = _parse_vcf_value(
            genodict["FA"], to_float=1, len_exactly=len(num_alt))
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Returns PASS, REJECT.
        
        # Mutect uses "PASS" to indicate a confident somatic variant.
        # INFO:SOMATIC  If present, indicate somatic mutation.
        # INFO:VT       SNP,INS,DEL

        assert len(var.filter_) == 1
        x = var.filter_[0]
        assert x in ["PASS", "REJECT"]
        return x
    def is_pass(self, filter_str):
        assert filter_str in ["PASS", "REJECT"]
        return filter_str == "PASS"


class MuTect2Caller(Caller):
    # FORMAT:DP   Approximate read depth (MQ=255, bad mates filtered)
    #             In header, but doesn't look like it's used.
    # FORMAT:AD   Depth for ref and alt in order listed
    # FORMAT:AF   alternate allele fraction  (variant allele frequency)
    # FORMAT:GT   Genotype
    def __init__(self):
        Caller.__init__(self, "MuTect2")
    def is_caller(self, header_lines):
        # MuTect
        # ##GATKCommandLine=<ID=MuTect,Version=3.1-0-g72492bb,Date="Sun
        # Mutect2
        # ##GATKCommandLine.MuTect2=<ID=MuTect2,Version=3.6-0-g89b7209,Date=
        # Make sure we're picking up MuTect2 and not MuTect.
        x = [x for x in header_lines if x.startswith("##GATKCommandLine")]
        x = [x for x in x if x.find("MuTect2") >= 0]
        if x:
            return True
        return False
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        #total_reads = _parse_vcf_value(
        #    genodict["DP"], to_int=1, len_exactly=1)[0]
        x = _parse_vcf_value(genodict["AD"], to_int=1, len_at_least=2)
        num_ref = x[0]
        num_alt = x[1:]
        vaf = _parse_vcf_value(
            genodict["AF"], to_float=1, len_exactly=len(num_alt))
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        assert len(num_alt) == 1
        total_reads = num_ref + num_alt[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Return:
        # alt_allele_in_normal, clustered_events, germline_risk,
        # homologous_mapping_event, multi_event_alt_allele_in_normal, 
        # panel_of_normals, str_contraction, t_lod_fstar, 
        # triallelic_site
        # PASS
        # Can be comma-separated combination as well.
        return ",".join(var.filter_)
    def is_pass(self, filter_str):
        return filter_str == "PASS"


class VarScan2Caller(Caller):
    # INFO:DP      Depth of quality bases
    # FORMAT:DP    Read depth
    # FORMAT:RD    reference reads
    # FORMAT:AD    variant reads
    # FORMAT:FREQ  variant allele frequency
    # FORMAT:DP4   read counts (ref/fwd, ref/rev, var/fwd, var/rev)
    # FORMAT:GT   Genotype
    def __init__(self):
        Caller.__init__(self, "VarScan2")
    def is_caller(self, header_lines):
        # ##source=VarScan2
        x = [x for x in header_lines if x.startswith("##source=VarScan2")]
        if x:
            return True
        return False
    #def is_caller_from_vcf(self, vcf):
    #    if not has_all_info(vcf, ["DP"]):
    #        return False
    #    if not has_all_format(vcf, ["DP", "RD", "AD", "FREQ", "DP4", "GT"]):
    #        return False
    #    return True
    def get_call(self, var, sample):
        # GENO  GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
        #       ./.:.:1

        num_ref = None
        num_alt = []
        vaf = []
        total_reads = None
        call = None

        # Sometimes these values are blank.
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        x = _parse_vcf_value(genodict["RD"], to_int=1, len_exactly=1)
        if x:
            num_ref = x[0]
        num_alt = _parse_vcf_value(genodict["AD"], to_int=1, len_exactly=1)
        vaf = _parse_vcf_value(genodict["FREQ"], perc_to_dec=1, len_exactly=1)
        x = _parse_vcf_value(genodict["DP"], to_int=1, len_exactly=1)
        if x:
            total_reads = x[0]
        x = _parse_vcf_value(genodict["GT"], len_exactly=1)
        if x:
            call = x[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Returns PASS, REFERENCE, GERMLINE, LOH, or UNKNOWN.

        # For somatic calling.
        if "SS" in var.infodict:
            # Varscan uses "PASS" to mean that the variant could be called
            # correctly, and indicates the somatic variant in the SS info.
            # INFO:SS  0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown
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
        # PASS, indelError, str10
        return ",".join(var.filter_)
        
    def is_pass(self, filter_str):
        # Same for SNP and indel.
        # For non-somatic, just PASS, indelError, or str10 (strand bias).
        assert filter_str in [
            "PASS", "indelError", "str10",
            "REFERENCE", "GERMLINE", "LOH", "UNKNOWN"]
        return filter_str == "PASS"


class SomaticSniperCaller(Caller):
    # FORMAT:DP      total read depth
    # FORMAT:DP4     high quality read (ref/fwd, ref/rev, var/fwd, var/rev)
    # FORMAT:BCOUNT  Occurrence for each base (A,C,G,T)
    # FORMAT:SSC     Somatic score
    # FORMAT:GT      Genotype
    def __init__(self):
        Caller.__init__(self, "SomaticSniper")
    def is_caller(self, header_lines):
        # ##FORMAT=<ID=SSC,Number=1,Type=Integer,Description="Somatic Score">
        x = [x for x in header_lines if x.startswith("##FORMAT=<ID=SSC")]
        x = [x for x in header_lines if x.find("Somatic Score") >= 0]
        if x:
            return True
        return False
    #def is_caller_from_vcf(self, vcf):
    #    if not has_all_format(vcf, ["DP", "DP4", "BCOUNT", "SSC", "GT"]):
    #        return False
    #    return True
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        x = _parse_vcf_value(genodict["DP4"], to_int=1, len_exactly=4)
        num_ref = x[0] + x[1]
        num_alt = [x[2] + x[3]]
        total_reads = _parse_vcf_value(
            genodict["DP"], to_int=1, len_exactly=1)[0]
        assert num_ref + num_alt[0] == total_reads
        vaf = [float(num_alt[0])/total_reads]
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # PASS, REFERENCE, GERMLINE, LOH, UNKNOWN.
        
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
            return "REFERENCE"
        elif ss1 == "1":
            return "GERMLINE"
        elif ss1 == "2":
            return "PASS"
        elif ss1 == "3":
            return "LOH"
        return "UNKNOWN"
    
    def is_pass(self, filter_str):
        assert filter_str in [
            "PASS", "REFERENCE", "GERMLINE", "LOH", "UNKNOWN"]
        return filter_str == "PASS"


class StrelkaCaller(Caller):
    # FORMAT:DP    Read depth for tier1
    # FORMAT:DP2   Read depth for tier2
    # INFO:NT      genotype

    # For SNP calls.
    # FORMAT:AU    A alleles, tiers 1,2
    # FORMAT:CU    C alleles, tiers 1,2
    # FORMAT:GU    G alleles, tiers 1,2
    # FORMAT:TU    T alleles, tiers 1,2

    # For INDEL calls.
    # FORMAT:TAR   reads for alternate alleles in tiers 1,2
    #              means reference and any other conflicting candidate indels
    # FORMAT:TIR   reads for indel allele in tiers 1,2
    # FORMAT:TOR   other reads for tiers 1,2
    #              other are alleles that cannot be distinguished
    
    def __init__(self):
        Caller.__init__(self, "Strelka")
    def is_caller(self, header_lines):
        # ##source=strelka
        x = [x for x in header_lines if x.startswith("##source=strelka")]
        if x:
            return True
        return False
    #def is_caller_from_vcf(self, vcf):
    #    if not has_all_info(vcf, ["NT"]):
    #        return False
    #    if not has_all_format(vcf, ["DP", "AU", "CU", "GU", "TU"]):
    #        return False
    #    return True
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]

        # For SNPs.
        if "AU" in genodict:
            total_reads = _parse_vcf_value(
                genodict["DP"], to_int=1, len_exactly=1)[0]
            # Includes values for two tiers.  Use tier 1 data.
            au = _parse_vcf_value(genodict["AU"], to_int=1, len_exactly=2)
            cu = _parse_vcf_value(genodict["CU"], to_int=1, len_exactly=2)
            gu = _parse_vcf_value(genodict["GU"], to_int=1, len_exactly=2)
            tu = _parse_vcf_value(genodict["TU"], to_int=1, len_exactly=2)
            au, cu, gu, tu = au[0], cu[0], gu[0], tu[0]
            # May have multiple ALT alleles.
            num_alt = []
            for alt in var.alt:
                if alt == "":
                    num_alt.append(0)
                elif alt == "A":
                    num_alt.append(au)
                elif alt == "C":
                    num_alt.append(cu)
                elif alt == "G":
                    num_alt.append(gu)
                elif alt == "T":
                    num_alt.append(tu)
                else:
                    raise AssertionError, "Unknown alt: %s" % alt
        # For INDELs.
        else:
            tar = _parse_vcf_value(genodict["TAR"], to_int=1, len_exactly=2)
            tir = _parse_vcf_value(genodict["TIR"], to_int=1, len_exactly=2)
            tor = _parse_vcf_value(genodict["TOR"], to_int=1, len_exactly=2)
            # Use Tier1 results for indel reads.
            total_reads = tar[0] + tir[0] + tor[0]
            num_alt = [tir[0]]
            
        num_ref = total_reads - sum(num_alt)
        if total_reads > 0:
            vaf = [float(x)/total_reads for x in num_alt]
        else:
            vaf = [None] * len(num_alt)
        call = _parse_vcf_value(var.infodict["NT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Return:
        # PASS, BCNoise, DP, QSS_ref
        # Can be comma-separated combination as well.

        # FILTER seems to be pretty good.
        # ##FILTER=<ID=DP,Description="Greater than 3.0x chromosomal mean depth
        #   in Normal sample">
        # ##FILTER=<ID=BCNoise,Description="Fraction of basecalls filtered at
        #   this site in either sample is at or above 0.4">
        # ##FILTER=<ID=SpanDel,Description="Fraction of reads crossing site
        #   with spanning deletions in either sample exceeeds 0.75">
        # ##FILTER=<ID=QSS_ref,Description="Normal sample is not homozygous ref
        #   or ssnv Q-score < 15, ie calls with NT!=ref or QSS_NT < 15">
        # INFO:SOMATIC   Present in every line
        return ",".join(var.filter_)

    def is_pass(self, filter_str):
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
        # 
        # Same for SNP and indels
        #
        # While VCF file is semicolon separated, SimpleVariantFile is
        # comma separated.
        return filter_str == "PASS"


class JointSNVMixCaller(Caller):
    # INFO:TR     tumor num_ref
    # INFO:TA     tumor num_alt
    # INFO:NR     normal num_ref
    # INFO:NA     normal num_alt
    # FORMAT:RD   Depth for ref in tumor sample  (TR+TA)
    # FORMAT:AD   Depth for alt in tumor sample
    def __init__(self):
        Caller.__init__(self, "mutationSeq")
    def is_caller(self, header_lines):
        # ##source=mutationSeq_4.3.6
        x = [x for x in header_lines if x.startswith("##source=mutationSeq")]
        if x:
            return True
        return False
    #def is_caller_from_vcf(self, vcf):
    #    if not has_all_info(vcf, ["TR", "TA", "NR", "NA"]):
    #        return False
    #    if not has_all_format(vcf, ["RD", "AD"]):
    #        return False
    #    return True
    def get_call(self, var, sample):
        #num_alt_alleles = 1
        num_ref = int(var.infodict["TR"])
        num_alt = [int(var.infodict["TA"])]
        total_reads = num_ref + sum(num_alt)
        vaf = [x/float(total_reads) for x in num_alt]
        call = None
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Return PASS
        
        # PASS  (for snp)
        # INDL  (for indel)   (passing indel)

        # Only PASS calls.  Only shows tumor data.  No information for
        # germline.
        assert len(var.filter_) == 1
        x = var.filter_[0]
        assert x in ["PASS", "INDL"]
        return "PASS"

    def is_pass(self, filter_str):
        assert filter_str == "PASS"
        return True


class MuSECaller(Caller):
    # FORMAT:DP   Read depth
    # FORMAT:AD   Depth for all alleles 0/1/2/3...
    # FORMAT:BQ   Average base quality
    # FORMAT:SS   Variant status
    # FORMAT:GT   Genotype
    def __init__(self):
        Caller.__init__(self, "MuSE")
    def is_caller(self, header_lines):
        # ##MuSE_version="v1.0rc Build Date Aug  6 2015 Build Time 10:11:56"
        x = [x for x in header_lines if x.startswith("##MuSE_version=")]
        if x:
            return True
        return False
    #def is_caller_from_vcf(self, vcf):
    #    if not has_all_format(vcf, ["DP", "AD", "BQ", "SS", "GT"]):
    #        return False
    #    return True
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        x = _parse_vcf_value(genodict["AD"], to_int=1, len_at_least=2)
        num_ref = x[0]
        num_alt = x[1:]
        total_reads = _parse_vcf_value(
            genodict["DP"], to_int=1, len_exactly=1)[0]
        vaf = [float(x)/total_reads for x in num_alt]
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # PASS, Tier1, Tier2, Tier3, Tier4, Tier5
        
        # MuSE uses "PASS", "Tier1", "Tier2", "Tier3", "Tier4", "Tier5".
        # INFO:SOMATIC  If present, indicate somatic mutation.
        #               Seems like always there.
        assert len(var.filter_) == 1
        x = var.filter_[0]
        assert x in ["PASS", "Tier1", "Tier2", "Tier3", "Tier4", "Tier5"]
        return x

    def is_pass(self, filter_str, wgs_or_wes):
        assert wgs_or_wes in ["wgs", "wes"]
        assert filter_str in [
            "PASS", "Tier1", "Tier2", "Tier3", "Tier4", "Tier5"]
        # According to paper, use up to Tier4 for WES, and Tier5 for WGS.
        if filter_str in ["PASS", "Tier1", "Tier2", "Tier3", "Tier4"]:
            return True
        if wgs_or_wes == "wes" and filter_str == "Tier5":
            return True
        return False


class RadiaCaller(Caller):
    # Generates really detailed output that includes a lot of
    # information about bad calls.
    # 
    # INFO:DP     Total read depth for all samples.
    # INFO:AF     Allele frequency for each ALT allele.
    # INFO:FA     Fraction of reads supporting ALT.
    # INFO:VT     Variant type.  SNP, INS, or DEL.
    # FORMAT:DP   Read depth at this position in the sample.
    # FORMAT:AD   Depth of reads supporting alleles
    # FORMAT:AF   Fraction of reads supporting alleles
    # FORMAT:SS   Variant status: 0=WT, 1=germ, 2=somatic, 3=LOH, 4=unk, 5=edit
    #             For SNPs: 0,3,4 not seen.
    def __init__(self):
        Caller.__init__(self, "Radia")
    def is_caller(self, header_lines):
        # ##source="RADIA pipeline v1.1.3"
        x = [x for x in header_lines if x.startswith('##source="RADIA')]
        if x:
            return True
        return False
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        total_reads = None
        x = _parse_vcf_value(genodict["DP"], to_int=1, len_exactly=1)
        if x:
            total_reads = x[0]
        num_ref, num_alt = None, []
        x = _parse_vcf_value(genodict["AD"], to_int=1, len_at_least=2)
        if x:
            num_ref = x[0]
            num_alt = x[1:]
        vaf = []
        if total_reads:
            vaf = [float(x)/total_reads for x in num_alt]
        call = None
        x = _parse_vcf_value(genodict["GT"], len_exactly=1)
        if x:
            call = x[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Return:
        # PASS, REFERENCE, GERMLINE, LOH, UNKNOWN, or one of the
        # error conditions found in FILTER.
        #
        # pbias   fails positional bias
        # multi   Multiple ALT alleles across all samples
        # dnmtb   dna normal total bases is below cutoff
        # ntr     Does not overlap TCGA target region
        # blck    Overlaps 1000 genomes
        # [...]
        assert var.filter_
        if len(var.filter_) != 1 or var.filter_[0] != "PASS":
            return ",".join(var.filter_)
        assert "SS" in var.infodict
        SS = var.infodict["SS"]
        assert SS in ["0", "1", "2", "3", "4", "5"]
        if SS == "2":
            return "PASS"
        elif SS == "0":
            return "REFERENCE"
        elif SS == "1":
            return "GERMLINE"
        elif SS == "3":
            return "LOH"
        elif SS == "5":
            return "EDIT"
        return "UNKNOWN"
        
        return ",".join(var.filter_)
    def is_pass(self, filter_str):
        return filter_str == "PASS"
    

class PindelCaller(Caller):
    # FORMAT:RD    Reference depth, how many reads support the reference
    #              Doesn't appear to be used?
    # FORMAT:AD    Allele depth, how many reads support this allele
    # FORMAT:GT    Genotype
    # INFO:SVTYPE  Type of structural variant
    def __init__(self):
        Caller.__init__(self, "Pindel")
    def is_caller(self, header_lines):
        # ##source=pindel
        x = [x for x in header_lines if x.startswith("##source=pindel")]
        if x:
            return True
        return False
    def get_call(self, var, sample):
        assert sample in var.sample2genodict, "Unknown sample: %s" % sample
        genodict = var.sample2genodict[sample]
        x = _parse_vcf_value(genodict["AD"], to_int=1, len_at_least=2)
        num_ref = x[0]
        num_alt = x[1:]
        assert len(num_alt) == 1
        # See:
        # PASS  GT:AD   0/0:-1,1
        # Don't understand what negative REF frequency means?  Just
        # set to 0.
        if num_ref < 0:
            num_ref = 0
        total_reads = num_ref + num_alt[0]
        if not total_reads:
            vaf = [0] * len(num_alt)
        else:
            vaf = [float(x)/total_reads for x in num_alt]
        call = _parse_vcf_value(genodict["GT"], len_exactly=1)[0]
        return Call(num_ref, num_alt, total_reads, vaf, call)
    def get_filter(self, var):
        # Always "PASS".
        # Return PASS.
        assert len(var.filter_) == 1
        x = var.filter_[0]
        assert x == "PASS"
        return x
    def is_pass(self, filter_str):
        assert filter_str == "PASS"
        return True


# List of classes
CALLERS = [
    samtoolsCaller,
    GATKCaller,
    PlatypusCaller,
    NextGENeCaller,
    # Somatic Callers.
    MuTectCaller,
    MuTect2Caller,
    VarScan2Caller,
    SomaticSniperCaller,
    StrelkaCaller,
    JointSNVMixCaller,
    MuSECaller,
    RadiaCaller,
    PindelCaller,
    ]


def read(filename):
    # Return a VCFFile object.
    import AnnotationMatrix

    caller = get_caller(filename)
    assert caller, "Can't identify caller: %s" % filename
    matrix = AnnotationMatrix.read(filename, header_char="##")
    return VCFFile(matrix, caller)


def write(handle_or_file, vcf):
    import AnnotationMatrix
    AnnotationMatrix.write(handle_or_file, vcf.matrix)


def get_caller(vcf_file):
    # Return a Caller object or None if the file is not recognized.
    from genomicode import filelib

    # Use the header lines to figure out the format.
    header_lines = []
    for line in filelib.openfh(vcf_file):
        if not line.startswith("##"):
            break
        header_lines.append(line)

    for caller in CALLERS:
        obj = caller()
        if obj.is_caller(header_lines):
            return obj
    return None


## def get_caller_from_vcf(vcf):
##     # Return a Caller object or None if the vcf object is not recognized.
##     from genomicode import filelib

##     for name, caller in CALLERS:
##         obj = caller()
##         if obj.is_caller_from_vcf(vcf):
##             return obj
##     return None


def get_variant(vcf, num):
    # Return a Variant object.
    assert num >= 0 and num < vcf.num_variants()

    chrom = vcf.matrix["#CHROM"][num]
    pos = int(vcf.matrix["POS"][num])
    id_ = vcf.matrix["ID"][num]

    # REF and ALT are usually separated by commas.  However, VarScan2
    # splits them by slashes, e.g. C/T .

    # REF should be string.
    ref = vcf.matrix["REF"][num]
    # Occasionally, with Platypus indel calls, there will be multiple
    # REFs.  Just ignore this case.
    # TGG/G
    # GCCC/CC
    if ref.find("/") >= 0:
        x = ref.split("/")
        ref = x[0]
    assert ref.find(",") < 0
    assert ref.find("/") < 0, "%s %s %s" % (chrom, pos, ref)
    assert ref != "."    

    # ALT should be string.
    x = vcf.matrix["ALT"][num]
    if x.find(",") >= 0:
        x = x.split(",")
    elif x.find("/") >= 0:
        x = x.split("/")
    else:
        x = [x]
    for i in range(len(x)):
        if x[i] == ".":
            x[i] = ""
    alt = x
        
    qual = _safe_float(vcf.matrix["QUAL"][num])
    filter_ = _parse_vcf_value(vcf.matrix["FILTER"][num], delim=";")
    x = _parse_info(vcf.matrix["INFO"][num])
    info_names, infodict = x

    # JoinSNVMix can be missing "FORMAT".  Ignore it if it doesn't exist.
    genotype_names = []
    sample2genodict = {}
    if "FORMAT" in vcf.matrix:
        format_str = vcf.matrix["FORMAT"][num]
        for sample in vcf.samples:
            geno_str = vcf.matrix[sample][num]
            try:
                x = _parse_genotype(format_str, geno_str)
            except AssertionError, x:
                # Mismatch: \
                #   GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR .
                if str(x).startswith("Mismatch"):
                    # Make a better error message.
                    x = str(x).split()
                    assert len(x) == 3
                    raise AssertionError, "Mismatch in %s:%d %s: %s %s" % (
                        chrom, pos, sample, x[1], x[2])
                raise
            genotype_names, genodict = x
            sample2genodict[sample] = genodict

    return Variant(
        chrom, pos, id_, ref, alt, qual, filter_, info_names, infodict,
        vcf.samples, genotype_names, sample2genodict, vcf.caller)


def set_variant(vcf, num, var):
    # Update a VCFFile in place with information from Variant.
    assert num >= 0 and num < vcf.num_variants()

    vcf.matrix["#CHROM"][num] = var.chrom
    vcf.matrix["POS"][num] = str(var.pos)
    vcf.matrix["ID"][num] = var.id_

    vcf.matrix["REF"][num] = _format_vcf_value(var.ref, join_char=",")
    vcf.matrix["ALT"][num] = _format_vcf_value(var.alt, join_char=",")
    vcf.matrix["QUAL"][num] = _format_vcf_value(var.qual, join_char=",")
    vcf.matrix["FILTER"][num] = _format_vcf_value(var.filter_, join_char=";")
    
    vcf.matrix["INFO"][num] = _format_info(var.info_names, var.infodict)
    if "FORMAT" in vcf.matrix:
        vcf.matrix["FORMAT"][num] = ":".join(var.genotype_names)
        for sample in var.samples:
            vcf.matrix[sample][num] = _format_genotype(
                var.genotype_names, var.sample2genodict[sample])
    else:
        assert not var.genotype_names, var.genotype_names
        assert not var.samples


def add_variant(vcf, var):
    # Add a variant to a VCFFile in place.
    num = vcf.num_variants()
    # Add a row to the matrix.
    matrix = vcf.matrix
    for key in matrix.header2annots:
        x = matrix.header2annots[key]
        x = x + ["."]
        matrix.header2annots[key] = x
    set_variant(vcf, num, var)


def get_call(var, sample):
    return var.caller.get_call(var, sample)


## def get_call(var, sample):
##     # Return a Call object from a variant.
##     if sample is not None:  # may be missing (e.g. JointSNVMix)
##         assert sample in var.samples, "Unknown sample: %s" % sample

##     info_dict = var.infodict
##     geno_dict = {}
##     if sample is not None:  # may be missing (e.g. JointSNVMix)
##         geno_dict = var.sample2genodict[sample]

##     # Figure out the number of ALT alleles.
##     num_alt_alleles = len(var.alt)
##     # Each is a list of values.
##     num_ref = []
##     num_alt = []
##     total_reads = []
##     vaf = []
##     # Call is None or a string.
##     call = None

##     # For debugging.
##     pos_str = "%s %s" % (var.chrom, var.pos)

##     # This function is a mess.  Should have separate functions for
##     # each caller.
##     if "RD" in geno_dict:
##         x = geno_dict["RD"]
##         # Previously, there was bug in _format_genotype that led to a
##         # list being formatted here.  Detect this bug and account for
##         # it.
##         if x.startswith("[") and x.endswith("]"):
##             x = x[1:-1]
##         num_ref = [_safe_int(x)]
##     if "SGCOUNTREF_F" in geno_dict:
##         x1 = _safe_int(geno_dict["SGCOUNTREF_F"])
##         x2 = _safe_int(geno_dict["SGCOUNTREF_R"])
##         num_ref = [_safe_add(x1, x2)]
##     if "SGCOUNTALT_F" in geno_dict:
##         x1 = geno_dict["SGCOUNTALT_F"].split(",")
##         x2 = geno_dict["SGCOUNTALT_R"].split(",")
##         assert len(x1) == len(x2)
##         x1 = [_safe_int(x) for x in x1]
##         x2 = [_safe_int(x) for x in x2]
##         num_alt = [_safe_add(x, y) for (x, y) in zip(x1, x2)]
##     if "AD" in geno_dict:
##         # ALT   AD
##         # C,T    2
##         # C    2,2   First is REF, second is ALT.
##         x = geno_dict["AD"].split(",")
##         x = [_safe_int(x) for x in x]
##         # Look for situation where first is REF and second is ALT.
##         # No RD means both AD and RD are merged together.
##         if "RD" not in geno_dict and len(x) == num_alt_alleles+1:
##             assert not num_ref
##             num_ref, num_alt = [x[0]], x[1:]

##             # Inconsistent read depths (AD and DP) from MuTect may
##             # lead to negative num_ref:
##             # GT:AD:BQ:DP:FA  0:92,179:.:140:0.661
##             # AD is actual depth, and DP is filtered for quality.
##             # If this is the case, then set total_reads based on AD,
##             # rather than DP.
##             total_reads = [num_ref[0]+x for x in num_alt]
##         else:
##             #assert len(x) == num_alt_alleles, "Bad AD: %s %d %s" % (
##             #    pos_str, num_alt_alleles, geno_dict["AD"])
##             num_alt = x
##     if not total_reads and "DP" in geno_dict:
##         x = geno_dict["DP"].split(",")
##         x = [_safe_int(x) for x in x]
##         total_reads = x

##     # SomaticSniper.
##     if not num_ref and not num_alt and "DP4" in geno_dict:
##         x = geno_dict["DP4"].split(",")
##         x = [_safe_int(x) for x in x]
##         assert len(x) == 4
##         ref_fwd, ref_rev, alt_fwd, alt_rev = x
##         num_ref = [ref_fwd + ref_rev]
##         num_alt = [alt_fwd + alt_rev]
        

##     if not total_reads and "DP" in info_dict:
##         x = info_dict["DP"]
##         assert type(x) is type(0)
##         total_reads = [x]
##     if "FREQ" in geno_dict:
##         x = geno_dict["FREQ"].split(",")
##         x = [_percent_to_decimal(x) for x in x]
##         x = [_safe_float(x) for x in x]
##         vaf = x
##     if "GT" in geno_dict:
##         call = geno_dict["GT"]
##         assert type(call) is type("")

##     # "TR" in info_dict means different things Platypus and
##     # JointSNVMix.  Need to be sure to interpret them correctly.
##     if not num_ref and not num_alt and "TR" in info_dict and "TA" in info_dict:
##         assert not total_reads
##         x = info_dict["TR"].split(",")
##         x = [_safe_int(x) for x in x]
##         assert len(x) == 1
##         num_ref = x
##         x = info_dict["TA"].split(",")
##         x = [_safe_int(x) for x in x]
##         assert len(x) == 1
##         num_alt = x
##         total_reads = [num_ref[0] + num_alt[0]]
##     if not num_alt and "TR" in info_dict:
##         assert "TA" not in info_dict
##         x = info_dict["TR"].split(",")
##         x = [_safe_int(x) for x in x]
##         # Don't bother checking.  Hard to know what Platypus is trying
##         # to do here.  Ex:
##         # ALT  T,TTGTGTGTGTG,TTGTGTGTG,TTGTGTG,TTGTG,TTG
##         # TR   14,14
##         #assert len(x) == 1 or len(x) == num_alt_alleles, \
##         #       "Mismatch alleles: %s %d %s %d %r" % (
##         #    sample, variant_num, alt_alleles, num_alt_alleles,
##         #    info_dict["TR"])
##         num_alt = x
##     if not total_reads and "TC" in info_dict:
##         x = info_dict["TC"].split(",")
##         x = [_safe_int(x) for x in x]
##         #assert len(x) == 1 or len(x) == num_alt_alleles
##         total_reads = x

##     if not num_ref and total_reads and num_alt:
##         # Possibilities:
##         # 1.  1 total reads, 1 alt.
##         # 2.  1 total reads, multiple alts.
##         #     Different alternative alleles.  Sum them up.
##         #     Actually, cannot sum up Platypus.  17 reads, alts: 14,14
##         # 3.  multiple total reads, multiple alts.
##         #     should have same number, and num_ref calculated should
##         #     be the same.

##         # Case 1.
##         if len(total_reads) == 1 and len(num_alt) == 1:
##             if total_reads[0] is not None and num_alt[0] is not None:
##                 num_ref = total_reads[0] - num_alt[0]
##                 assert num_ref >= 0, "%s: %s %s" % (
##                     pos_str, num_alt[0], total_reads[0])
##         # Case 2.
##         elif len(total_reads) == 1 and len(num_alt) > 1:
##             x = [x for x in num_alt if x is not None]
##             #num_ref = total_reads[0] - sum(x)
##             num_ref = max(total_reads) - max(x)
##             assert num_ref >= 0, "%s: %s %s" % (pos_str, num_alt, total_reads)
##         # Case 3.
##         elif len(total_reads) > 1 and len(num_alt) > 1:
##             assert len(total_reads) == len(num_alt), \
##                    "%s: %s %s" % (pos_str, num_alt, total_reads)
##             calc_nr = []
##             for i in range(len(total_reads)):
##                 if total_reads[i] is not None and num_alt[i] is not None:
##                     x = total_reads[i] - num_alt[i]
##                     calc_nr.append(x)
##             calc_nr.sort()
##             if calc_nr:
##                 assert calc_nr[0] == calc_nr[-1]
##                 num_ref = calc_nr[0]
##         else:
##             raise AssertionError, "%s: %s %s" % (pos_str, num_alt, total_reads)

##     if not vaf and num_ref and num_alt:
##         calc_v = [None] * len(num_alt)
##         for i in range(len(num_alt)):
##             if num_alt[i] is None:
##                 continue
##             assert len(num_ref) == 1   # Is this true?
##             total = num_ref[0] + num_alt[i]
##             if total:
##                 x = float(num_alt[i])/total
##                 calc_v[i] = x
##         vaf = calc_v

##     # If no other data available, then use the information from BFILL_
##     # fields.
##     if not num_ref and "BFILL_REF" in geno_dict:
##         num_ref = [_safe_int(geno_dict["BFILL_REF"])]
##     if not num_alt and "BFILL_ALT" in geno_dict:
##         num_alt = [_safe_int(geno_dict["BFILL_ALT"])]
##     if not total_reads and "BFILL_COV" in geno_dict:
##         total_reads = [_safe_int(geno_dict["BFILL_COV"])]
##     if not vaf and "BFILL_VAF" in geno_dict:
##         vaf = [_safe_float(geno_dict["BFILL_VAF"])]

##     assert type(num_ref) is type([])
##     assert type(num_alt) is type([])
##     assert type(total_reads) is type([])
##     assert type(vaf) is type([])

##     # Convert lists to numbers.
##     if type(num_alt) is type([]) and len(num_alt) == 1:
##         num_alt = num_alt[0]
##     if type(total_reads) is type([]) and len(total_reads) == 1:
##         total_reads = total_reads[0]
##     if type(vaf) is type([]) and len(vaf) == 1:
##         vaf = vaf[0]

##     # If total_reads is 0, then sometimes does not include the other ones.
##     if total_reads == 0 and not num_ref and not num_alt:
##         num_ref = 0
##         num_alt = 0
##     if total_reads == 0 and not vaf:
##         vaf = 0.0

##     return Call(
##         num_alt_alleles, num_ref, num_alt, total_reads, vaf, call)


def set_call(variant, sample, call):
    # Update a Variant object in place with the information from the
    # Call object.
    ID = variant.infodict
    GD = variant.sample2genodict[sample]

    # TODO: separate this out into objects.

    # Figure out what kind of file it is.
    if "BFILL_REF" in GD and "BFILL_ALT" in GD and "BFILL_COV" in GD and \
       "BFILL_VAF" in GD:
        if call.num_ref is not None:
            GD["BFILL_REF"] = _format_vcf_value(call.num_ref)
        if call.num_alt is not None:
            GD["BFILL_ALT"] = _format_vcf_value(call.num_alt)
        if call.total_reads is not None:
            GD["BFILL_COV"] = _format_vcf_value(call.total_reads)
        if call.vaf is not None:
            GD["BFILL_VAF"] = _format_vcf_value(call.vaf)
    elif "RD" in GD and "AD" in GD and "DP" in GD and "FREQ" in GD:
        # samtools
        GD["RD"] = _format_vcf_value(call.num_ref)
        GD["AD"] = _format_vcf_value(call.num_alt)
        GD["DP"] = _format_vcf_value(call.total_reads)
        # Convert FREQ to percent.
        x = call.vaf
        if type(x) is not type([]):
            x = [x]
        x = [x*100.0 for x in x]
        x = ["%s%%" % x for x in x]
        GD["FREQ"] = _format_vcf_value(x)
        GD["GT"] = call.call
    elif "TR" in ID and "TC" in ID:
        # Platypus
        ID["TR"] = _format_vcf_value(call.num_ref)
        ID["TC"] = _format_vcf_value(call.total_reads)
        GD["GT"] = call.call
    elif "AD" in GD and "DP" in GD and not "RD" in GD and not "FREQ" in GD:
        # GATK
        GD["RD"] = _format_vcf_value(call.num_ref)
        GD["AD"] = _format_vcf_value(call.num_alt)
        GD["DP"] = _format_vcf_value(call.total_reads)
        GD["GT"] = call.call
    elif "DP" in ID and "GT" in GD and not "AD" in GD:
        # bcftools
        ID["DP"] = _format_vcf_value(call.total_reads)
        GD["GT"] = call.call
    elif "SGCOUNTREF_F" in GD:
        # NextGene
        raise NotImplementedError, "NextGene not implemented yet."
    else:
        raise AssertionError, "Unknown VCF format."


def is_pass_filter(var, QUAL=None, MQ=None, QD=None,
                   min_DP=None, max_DP=None, FILTER_contains=None,
                   FILTER_doesnotcontain=None):
    # Returns True if a variant passes all these filters.
    # <var value> > QUAL    Phred-scale quality score.  Usually 20.
    # <var value> > MQ      Average Mapping Quality.  Usually 40.
    # <var value> > QD      Quality by depth.  Usually 2.0.
    # <var value> > min_DP  Minimum depth.
    # <var value> < max_DP  Maximum depth.
    # FILTER has FILTER_contains (string).
    # FILTER does not have FILTER_doesnotcontain (string).

    ## if FILTER_contains is None:
    ##     FILTER_contains = ""
    ## if FILTER_doesnotcontain is None:
    ##     FILTER_doesnotcontain = ""
    ## if type(FILTER_doesnotcontain) is type(""):
    ##     FILTER_doesnotcontain = [FILTER_doesnotcontain]
    ## # Make sure there's no overlap between FILTER_contains and
    ## # FILTER_doesnotcontain.
    ## x = [x for x in FILTER_contains if x in FILTER_doesnotcontain]
    ## assert not x, "Overlap between FILTER_contains and FILTER_doesnotcontain"

    if QUAL is not None:
        assert var.qual is not None, "Missing QUAL: %s %d" % (
            var.chrom, var.pos)
        if var.qual < QUAL:
            return False
    if MQ is not None:
        assert "MQ" in var.info_names, "Missing MQ: %s %d" % (
            var.chrom, var.pos)
        if var.infodict["MQ"] < MQ:
            return False
    if QD is not None:
        assert "QD" in var.info_names, "Missing QD: %s %d" % (
            var.chrom, var.pos)
        if var.infodict["QD"] < QD:
            return False
    DP = None
    if min_DP is not None or max_DP is not None:
        if "DP" in var.infodict:
            DP = var.infodict["DP"]
    # TODO: Try other ways of getting depth.
    if min_DP is not None:
        assert DP is not None, "Missing DP: %s %d" % (var.chrom, var.pos)
        if DP < min_DP:
            return False
    if max_DP is not None:
        assert DP is not None, "Missing DP: %s %d" % (var.chrom, var.pos)
        if DP > max_DP:
            return False

    if FILTER_contains:
        if FILTER_contains not in var.filter_:
            return False
    if FILTER_doesnotcontain:
        if FILTER_doesnotcontain in var.filter_:
            return False
    return True


def select_variants(vcf, is_selected_fn):
    # is_selected_fn takes a Variant function and return a bool that
    # indicates whether to keep this variant.  Return a VCFFile.
    from genomicode import AnnotationMatrix

    I = []
    for i in range(vcf.num_variants()):
        var = get_variant(vcf, i)
        if is_selected_fn(var):
            I.append(i)
    matrix = AnnotationMatrix.rowslice(vcf.matrix, I)
    return VCFFile(matrix, vcf.caller)


def make_coverage_matrix(vcf, samples=None):
    # Make a num_variants x num_samples matrix where each element is
    # an integer or None.
    assert vcf.num_variants()

    if samples is None:
        samples = vcf.samples
    for s in samples:
        assert s in vcf.samples

    matrix = []
    for i in range(vcf.num_variants()):
        row = []
        var = get_variant(vcf, i)
        for s in samples:
            call = get_call(var, s)
            x = call.total_reads
            assert type(x) in [type(None), type(0)]
            row.append(x)
        matrix.append(row)
    return matrix


def make_vaf_matrix(vcf, samples=None):
    # Each element in the matrix is a float or None.
    assert vcf.num_variants()

    if samples is None:
        samples = vcf.samples
    for s in samples:
        assert s in vcf.samples

    matrix = []
    for i in range(vcf.num_variants()):
        row = []
        var = get_variant(vcf, i)
        for s in samples:
            call = get_call(var, s)
            x = call.vaf
            assert type(x) in [type(None), type(0.0)], x
            row.append(x)
        matrix.append(row)
    return matrix


def _safe_int(x):
    if x in ["", ".", None]:
        return None
    return int(x)


def _safe_float(x):
    if x in ["", ".", None]:
        return None
    return float(x)


def _safe_div_z(x1, x2):
    # Safe divide and ignore zero denominator.  x1/x2.  x1 and x2
    # should be integers, float, or None.  If both are None, return
    # None.  If only 1 is None, interpret it as 0.  If x2 is 0, return
    # 0.
    if x1 is None and x2 is None:
        return None
    if x1 is None:
        x1 = 0
    if x2 is None:
        x2 = 0
    if abs(x2) < 1E-20:
        return 0
    return float(x1)/x2


def _safe_add(x1, x2):
    # Add two numbers.  x1 and x2 should be integers or None.  If both
    # are None, return None.  If only 1 is None, interpret it as 0.
    if x1 is None and x2 is None:
        return None
    if x1 is None:
        x1 = 0
    if x2 is None:
        x2 = 0
    return x1 + x2


def _percent_to_decimal(x):
    # If the VAF is provided as a percent, convert it to decimal.
    # e.g. 14.29% -> 0.1429
    if type(x) != type(""):
        return x
    if not x.endswith("%"):
        return x
    x = x[:-1]
    return float(x) / 100


def _parse_info(info_str):
    # Return a tuple of info_names, info_dict

    # Parse out the INFO line.
    # Examples:
    # BRF=0.89;FR=1.0000;HP=2;HapScore=1;...
    # .
    # DB

    # No.  Don't do any conversion.  Different files may use these
    # differently.
    # Which elements should be int or float.
    #is_int = ["MQ0", "DP", "AN"]
    #is_float = ["FS", "QD", "SOR", "MQ"]
    # MLEAC  2  1,1
    # MLEAF  0.500,0.500
    # AC     1,1
    # AF     0.500,0.500

    if info_str == ".":
        return [], {}

    names = []
    d = {}
    x = info_str.split(";")
    for x in x:
        x = x.split("=")
        assert len(x) in [1, 2]
        if len(x) == 1:
            key = x[0]
            value = None
        else:
            key, value = x
            #if key in is_int:
            #    value = int(value)
            #if key in is_float:
            #    value = float(value)
        d[key] = value
        names.append(key)
    return names, d


def _parse_genotype(format_str, genotype_str):
    # Return tuple of format_names, genotype_dict.
    names = format_str.split(":")
    values = genotype_str.split(":")

    # Some weird conditions to handle:
    # - Varscan sometimes generates lines like:
    #   GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
    #   ./.:.:3
    # - Another Varscan:
    #   GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
    #   .
    # - Radia:
    #   GT:DP:INDEL:START:STOP:AD:AF:BQ:SB
    #   .
    # Detect this and fix it by padding with ".".
    if len(names) > 3 and len(values) in [1, 3] and \
           names[:3] == ["GT", "GQ", "SDP"]:
        x = len(names) - len(values)
        values = values + ["."]*x
    if len(values) == 1 and values[0] == "." and len(names) > 1:
        values = ["."] * len(names)

    assert len(names) == len(values), "Mismatch: %s %s" % (
        format_str, genotype_str)
    d = {}
    for (n, v) in zip(names, values):
        d[n] = v
    return names, d


def _format_info(info_names, info_dict):
    not_found = [x for x in info_dict if x not in info_names]
    assert not not_found, \
           "Unknown names in info_dict: %s" % ", ".join(not_found)

    keyvalues = []
    for name in info_names:
        value = info_dict[name]
        if value is None:
            keyvalues.append(name)
        else:
            keyvalues.append("%s=%s" % (name, value))
    return ";".join(keyvalues)


def _format_genotype(genotype_names, genotype_dict):
    not_found = [x for x in genotype_dict if x not in genotype_names]
    assert not not_found, \
           "Unknown names in genotype_dict: %s" % ", ".join(not_found)

    values = []
    for name in genotype_names:
        if name in genotype_dict:
            x = genotype_dict[name]
            if type(x) in [type([]), type(())]:
                x = ",".join(map(str, x))
            values.append(x)
        else:
            values.append(".")
    values = map(str, values)
    return ":".join(values)


def _parse_vcf_value(
    value, to_int=False, to_float=False, perc_to_dec=False, len_exactly=None,
    len_at_least=None, delim=","):
    # Return a list of strings or None.  Can be an empty list.
    # Examples:
    # .
    # 4
    # 4,4
    # 4,.
    if value == ".":
        return []
    values = value.split(delim)
    for i in range(len(values)):
        if values[i] == ".":
            values[i] = None
    if to_int:
        values = [_safe_int(x) for x in values]
    if to_float:
        values = [_safe_float(x) for x in values]
    if perc_to_dec:
        values = [_percent_to_decimal(x) for x in values]
    if len_exactly is not None:
        assert len(values) == len_exactly, "Expected len %d: %s" % (
            len_exactly, value)
    if len_at_least is not None:
        assert len(values) >= len_at_least
    return values


def _format_vcf_value(value, None_char=".", join_char=","):
    # value can be:
    # - string, int, float, None
    # - list of any of these
    # None_char is used for missing values.  Default is ".".
    # join_char is the character used to join lists.

    # Make a list for consistency.
    if type(value) is not type([]):
        value = [value]
    for i in range(len(value)):
        x = value[i]
        if x is None:
            x = None_char
        elif type(x) in [type(0), type(0.0)]:
            x = str(x)
        value[i] = x
    x = join_char.join(value)
    return x


def has_info(vcf, name, num=None):
    if num is None:
        num = 0
    assert num >= 0 and num < vcf.num_variants()
    var = get_variant(vcf, num)
    if name in var.infodict:
        return True
    return False


def has_all_info(vcf, names, num=None):
    for name in names:
        if not has_info(vcf, name, num=num):
            return False
    return True

def has_format(vcf, name, num=None, sample=None):
    # JointSNVMix sometimes does not have any samples.
    if not vcf.samples:
        return False
    if sample is None:
        sample = vcf.samples[0]
    var = get_variant(vcf, num)
    genodict = var.sample2genodict[sample]
    if name in genodict:
        return True
    return False


def has_all_format(vcf, names, num=None, sample=None):
    for name in names:
        if not has_format(vcf, name, num=num, sample=sample):
            return False
    return True


def simplify_call(call):
    # Return a Call object.  If there are multiple num_alt or vaf,
    # choose the one with the highest vaf.
    assert len(call.num_alt) == len(call.vaf), "Mismatch: %s %s" % (
        call.num_alt, call.vaf)
    num_alt = call.num_alt[:]
    vaf = call.vaf[:]
    if len(num_alt) > 1:
        max_vaf = max_i = None
        for i in range(len(vaf)):
            if max_vaf is None or vaf[i] > max_vaf:
                max_vaf = vaf[i]
                max_i = i
        assert max_i is not None
        num_alt = [num_alt[i]]
        vaf = [vaf[i]]
    x = Call(call.num_ref, num_alt, call.total_reads, vaf, call.call)
    return x
