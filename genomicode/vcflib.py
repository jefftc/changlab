"""
Classes:
VCFFile
Variant
Call

Functions:
read
write

has_info
has_all_info
has_format
has_all_format

get_variant      Get Variant from VCFFile.
set_variant      Set a Variant to a VCFFile.
add_variant      Add a Variant to a VCFFile.
get_call         Get a Call from a Variant.
set_call         Set a Call to a Variant.
is_pass_filter   Whether the variant passes a filter.
select_variants  Select a subset of the variants.

make_coverage_matrix
make_vaf_matrix

"""
# _parse_info
# _parse_genotype
# _format_info
# _format_genotype
# _format_vcf_value
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
# FILTER
# .
# badReads
# Q20;badReads
# PASS
#
#
# JoinSNVMix has a non-standard VCF file.  No "FORMAT" column, and
# does not have information about each sample.  Everything in INFO.


# Callers.
#SAMTOOLS = 0
#BCFTOOLS = 1
#GATK = 2
#PLATYPUS = 3
#JOINTSNVMIX = 4
#SOMATICSNIPER = 5
#STRELKA = 6
#VARSCAN = 7
#NEXTGENE = 8
#BACKFILL = 9

class Caller:
    def is_caller(self, vcf):
        raise NotImplementedError
    def get_call(self, var, sample):
        raise NotImplementedError

class JointSNVMixCaller(Caller):
    # INFO:TR     tumor num_ref
    # INFO:TA     tumor num_alt
    # INFO:NR     normal num_ref
    # INFO:NA     normal num_alt
    # FORMAT:RD   Depth for ref in tumor sample  (TR+TA)
    # FORMAT:AD   Depth for alt in tumor sample
    def is_caller(self, vcf):
        if not has_all_info(vcf, ["TR", "TA", "NR", "NA"]):
            return False
        if not has_all_format(vcf, ["RD", "AD"]):
            return False
        return True
    def get_call(self, var, sample):
        num_alt_alleles = 1
        num_ref = int(var.infodict["TR"])
        num_alt = int(var.infodict["TA"])
        total_reads = num_ref + num_alt
        vaf = num_alt / float(total_reads)
        call = None
        return Call(1, num_ref, num_alt, total_reads, vaf, call)

# TODO: Implement all the other callers.



class VCFFile:
    def __init__(self, matrix):
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

    def num_variants(self):
        return self.matrix.num_annots()




class Variant:
    # Holds structured information for each variant.
    #
    # chrom            string
    # pos              int
    # id_              string
    # ref              string or list of strings
    # alt              string or list of strings
    # qual             float or None
    # filter_          list of strings
    # info_names       list of strings
    # infodict
    # samples          list of strings
    # genotype_names   list of strings.
    #                  Can be empty list if no FORMAT (jointsnvmix).
    # sample2genodict  Can be empty dict if no FORMAT.
    def __init__(
        self, chrom, pos, id_, ref, alt, qual, filter_, info_names, infodict,
        samples, genotype_names, sample2genodict):
        import copy
        self.chrom = chrom
        self.pos = pos
        self.id_ = id_
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter_ = filter_
        self.info_names = info_names[:]
        self.infodict = infodict.copy()
        self.samples = samples[:]
        self.genotype_names = genotype_names[:]
        self.sample2genodict = copy.deepcopy(sample2genodict)
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
    # num_alt_alleles
    # num_ref          None, int.
    # num_alt          Can be list.  Each can be None.
    # total_reads      Can be list.  Each can be None.
    # vaf              Can be list.  Each can be None.
    # call             String or None

    def __init__(self, num_alt_alleles, num_ref, num_alt, total_reads, vaf,
                 call):
        assert num_ref is None or num_ref >= 0, "Invalid num_ref: %s" % num_ref
        assert num_alt_alleles >= 1
        self.num_alt_alleles = num_alt_alleles
        self.num_ref = num_ref
        self.num_alt = num_alt
        self.total_reads = total_reads
        self.vaf = vaf
        self.call = call
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [repr(self.num_alt_alleles),
             repr(self.num_ref),
             repr(self.num_alt),
             repr(self.total_reads),
             repr(self.vaf),
             repr(self.call),]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x
    def __cmp__(self, other):
        if not isinstance(other, Call):
            return cmp(id(self), id(other))
        x1 = [
            self.num_alt_alleles, self.num_ref, self.num_alt,
            self.total_reads, self.vaf, self.call]
        x2 = [
            other.num_alt_alleles, other.num_ref, other.num_alt,
            other.total_reads, other.vaf, other.call]
        return cmp(x1, x2)


def read(filename):
    # Return a VCFFile object.
    import AnnotationMatrix
    matrix = AnnotationMatrix.read(filename, header_char="##")
    return VCFFile(matrix)


def write(handle_or_file, vcf):
    import AnnotationMatrix
    AnnotationMatrix.write(handle_or_file, vcf.matrix)


def get_variant(vcf, num):
    # Return a Variant object.
    assert num >= 0 and num < vcf.num_variants()

    chrom = vcf.matrix["#CHROM"][num]
    pos = int(vcf.matrix["POS"][num])
    id_ = vcf.matrix["ID"][num]
    ref = vcf.matrix["REF"][num]
    # VarScan2 REF and ALT can be separated by slashes, e.g. C/T  .
    if ref.find(",") >= 0:
        ref = ref.split(",")
    elif ref.find("/") >= 0:
        ref = ref.split("/")
    alt = vcf.matrix["ALT"][num]
    if alt.find(",") >= 0:
        alt = alt.split(",")
    elif alt.find("/") >= 0:
        alt = alt.split("/")
    qual = _safe_float(vcf.matrix["QUAL"][num])
    filter_ = vcf.matrix["FILTER"][num].split(";")
    x = _parse_info(vcf.matrix["INFO"][num])
    info_names, infodict = x

    # JoinSNVMix can be missing "FORMAT".
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
        vcf.samples, genotype_names, sample2genodict)


def set_variant(vcf, num, var):
    # Update a VCFFile in place with information from Variant.
    assert num >= 0 and num < vcf.num_variants()

    vcf.matrix["#CHROM"][num] = var.chrom
    vcf.matrix["POS"][num] = str(var.pos)
    vcf.matrix["ID"][num] = var.id_
    ref = var.ref
    if type(ref) is not type(""):
        ref = ",".join(ref)
    vcf.matrix["REF"][num] = ref
    alt = var.alt
    if type(alt) is not type(""):
        alt = ",".join(alt)
    vcf.matrix["ALT"][num] = alt
    qual = var.qual
    if qual is None:
        qual = "."
    vcf.matrix["QUAL"][num] = str(qual)
    assert type(var.filter_) is type([])
    vcf.matrix["FILTER"][num] = ";".join(var.filter_)
    vcf.matrix["INFO"][num] = _format_info(var.info_names, var.infodict)
    if "FORMAT" in vcf.matrix:
        format_str = ":".join(var.genotype_names)
        vcf.matrix["FORMAT"][num] = format_str
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
    # Return a Call object from a variant.
    if sample is not None:  # may be missing (e.g. JointSNVMix)
        assert sample in var.samples, "Unknown sample: %s" % sample

    info_dict = var.infodict
    geno_dict = {}
    if sample is not None:  # may be missing (e.g. JointSNVMix)
        geno_dict = var.sample2genodict[sample]

    # Figure out the number of ALT alleles.
    num_alt_alleles = len(var.alt)
    # Each is a list of values.
    num_ref = []
    num_alt = []
    total_reads = []
    vaf = []
    # Call is None or a string.
    call = None

    # For debugging.
    pos_str = "%s %s" % (var.chrom, var.pos)

    # This function is a mess.  Should have separate functions for
    # each caller.
    if "RD" in geno_dict:
        num_ref = [_safe_int(geno_dict["RD"])]
    if "SGCOUNTREF_F" in geno_dict:
        x1 = _safe_int(geno_dict["SGCOUNTREF_F"])
        x2 = _safe_int(geno_dict["SGCOUNTREF_R"])
        num_ref = [_safe_add(x1, x2)]
    if "SGCOUNTALT_F" in geno_dict:
        x1 = geno_dict["SGCOUNTALT_F"].split(",")
        x2 = geno_dict["SGCOUNTALT_R"].split(",")
        assert len(x1) == len(x2)
        x1 = [_safe_int(x) for x in x1]
        x2 = [_safe_int(x) for x in x2]
        num_alt = [_safe_add(x, y) for (x, y) in zip(x1, x2)]
    if "AD" in geno_dict:
        # ALT   AD
        # C,T    2
        # C    2,2   First is REF, second is ALT.
        x = geno_dict["AD"].split(",")
        x = [_safe_int(x) for x in x]
        # Look for situation where first is REF and second is ALT.
        # No RD means both AD and RD are merged together.
        if "RD" not in geno_dict and len(x) == num_alt_alleles+1:
            assert not num_ref
            num_ref, num_alt = [x[0]], x[1:]

            # Inconsistent read depths (AD and DP) from MuTect may
            # lead to negative num_ref:
            # GT:AD:BQ:DP:FA  0:92,179:.:140:0.661
            # AD is actual depth, and DP is filtered for quality.
            # If this is the case, then set total_reads based on AD,
            # rather than DP.
            total_reads = [num_ref[0]+x for x in num_alt]
        else:
            #assert len(x) == num_alt_alleles, "Bad AD: %s %d %s" % (
            #    pos_str, num_alt_alleles, geno_dict["AD"])
            num_alt = x
    if not total_reads and "DP" in geno_dict:
        x = geno_dict["DP"].split(",")
        x = [_safe_int(x) for x in x]
        total_reads = x

    # SomaticSniper.
    if not num_ref and not num_alt and "DP4" in geno_dict:
        x = geno_dict["DP4"].split(",")
        x = [_safe_int(x) for x in x]
        assert len(x) == 4
        ref_fwd, ref_rev, alt_fwd, alt_rev = x
        num_ref = [ref_fwd + ref_rev]
        num_alt = [alt_fwd + alt_rev]
        

    if not total_reads and "DP" in info_dict:
        x = info_dict["DP"]
        assert type(x) is type(0)
        total_reads = [x]
    if "FREQ" in geno_dict:
        x = geno_dict["FREQ"].split(",")
        x = [_percent_to_decimal(x) for x in x]
        x = [_safe_float(x) for x in x]
        vaf = x
    if "GT" in geno_dict:
        call = geno_dict["GT"]
        assert type(call) is type("")

    # "TR" in info_dict means different things Platypus and
    # JointSNVMix.  Need to be sure to interpret them correctly.
    if not num_ref and not num_alt and "TR" in info_dict and "TA" in info_dict:
        assert not total_reads
        x = info_dict["TR"].split(",")
        x = [_safe_int(x) for x in x]
        assert len(x) == 1
        num_ref = x
        x = info_dict["TA"].split(",")
        x = [_safe_int(x) for x in x]
        assert len(x) == 1
        num_alt = x
        total_reads = [num_ref[0] + num_alt[0]]
    if not num_alt and "TR" in info_dict:
        assert "TA" not in info_dict
        x = info_dict["TR"].split(",")
        x = [_safe_int(x) for x in x]
        # Don't bother checking.  Hard to know what Platypus is trying
        # to do here.  Ex:
        # ALT  T,TTGTGTGTGTG,TTGTGTGTG,TTGTGTG,TTGTG,TTG
        # TR   14,14
        #assert len(x) == 1 or len(x) == num_alt_alleles, \
        #       "Mismatch alleles: %s %d %s %d %r" % (
        #    sample, variant_num, alt_alleles, num_alt_alleles,
        #    info_dict["TR"])
        num_alt = x
    if not total_reads and "TC" in info_dict:
        x = info_dict["TC"].split(",")
        x = [_safe_int(x) for x in x]
        #assert len(x) == 1 or len(x) == num_alt_alleles
        total_reads = x

    if not num_ref and total_reads and num_alt:
        # Possibilities:
        # 1.  1 total reads, 1 alt.
        # 2.  1 total reads, multiple alts.
        #     Different alternative alleles.  Sum them up.
        #     Actually, cannot sum up Platypus.  17 reads, alts: 14,14
        # 3.  multiple total reads, multiple alts.
        #     should have same number, and num_ref calculated should
        #     be the same.

        # Case 1.
        if len(total_reads) == 1 and len(num_alt) == 1:
            if total_reads[0] is not None and num_alt[0] is not None:
                num_ref = total_reads[0] - num_alt[0]
                assert num_ref >= 0, "%s: %s %s" % (
                    pos_str, num_alt[0], total_reads[0])
        # Case 2.
        elif len(total_reads) == 1 and len(num_alt) > 1:
            x = [x for x in num_alt if x is not None]
            #num_ref = total_reads[0] - sum(x)
            num_ref = max(total_reads) - max(x)
            assert num_ref >= 0, "%s: %s %s" % (pos_str, num_alt, total_reads)
        # Case 3.
        elif len(total_reads) > 1 and len(num_alt) > 1:
            assert len(total_reads) == len(num_alt), \
                   "%s: %s %s" % (pos_str, num_alt, total_reads)
            calc_nr = []
            for i in range(len(total_reads)):
                if total_reads[i] is not None and num_alt[i] is not None:
                    x = total_reads[i] - num_alt[i]
                    calc_nr.append(x)
            calc_nr.sort()
            if calc_nr:
                assert calc_nr[0] == calc_nr[-1]
                num_ref = calc_nr[0]
        else:
            raise AssertionError, "%s: %s %s" % (pos_str, num_alt, total_reads)

    if not vaf and num_ref and num_alt:
        calc_v = [None] * len(num_alt)
        for i in range(len(num_alt)):
            if num_alt[i] is None:
                continue
            assert len(num_ref) == 1   # Is this true?
            total = num_ref[0] + num_alt[i]
            if total:
                x = float(num_alt[i])/total
                calc_v[i] = x
        vaf = calc_v

    # If no other data available, then use the information from BFILL_
    # fields.
    if not num_ref and "BFILL_REF" in geno_dict:
        num_ref = [_safe_int(geno_dict["BFILL_REF"])]
    if not num_alt and "BFILL_ALT" in geno_dict:
        num_alt = [_safe_int(geno_dict["BFILL_ALT"])]
    if not total_reads and "BFILL_COV" in geno_dict:
        total_reads = [_safe_int(geno_dict["BFILL_COV"])]
    if not vaf and "BFILL_VAF" in geno_dict:
        vaf = [_safe_float(geno_dict["BFILL_VAF"])]

    assert type(num_ref) is type([])
    assert type(num_alt) is type([])
    assert type(total_reads) is type([])
    assert type(vaf) is type([])

    # Convert lists to numbers.
    if type(num_alt) is type([]) and len(num_alt) == 1:
        num_alt = num_alt[0]
    if type(total_reads) is type([]) and len(total_reads) == 1:
        total_reads = total_reads[0]
    if type(vaf) is type([]) and len(vaf) == 1:
        vaf = vaf[0]

    # If total_reads is 0, then sometimes does not include the other ones.
    if total_reads == 0 and not num_ref and not num_alt:
        num_ref = 0
        num_alt = 0
    if total_reads == 0 and not vaf:
        vaf = 0.0

    return Call(
        num_alt_alleles, num_ref, num_alt, total_reads, vaf, call)


def set_call(variant, sample, call):
    # Update a Variant object in place with the information from the
    # Call object.
    ID = variant.infodict
    GD = variant.sample2genodict[sample]

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
    return VCFFile(matrix)


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
    if x != ".":
        return int(x)
    return None


def _safe_float(x):
    if x == "":
        return None
    if x != ".":
        return float(x)
    return None


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

    # Which elements should be int or float.
    is_int = ["MQ0", "DP", "AN"]
    is_float = ["FS", "QD", "SOR", "MQ"]
    # MLEAC  2  1,1
    # MLEAF  0.500,0.500
    # AC     1,1
    # AF     0.500,0.500

    # Parse out the INFO line.
    # Format: BRF=0.89;FR=1.0000;HP=2;HapScore=1;...
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
        if key in is_int:
            value = int(value)
        if key in is_float:
            value = float(value)
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
    # Detect this and fix it by padding with ".".
    if len(names) > 3 and len(values) in [1, 3] and \
           names[:3] == ["GT", "GQ", "SDP"]:
        x = len(names) - len(values)
        values = values + ["."]*x

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
            values.append(genotype_dict[name])
        else:
            values.append(".")
    values = map(str, values)
    return ":".join(values)


def _format_vcf_value(value, char_for_None="."):
    # value can be:
    # string, int, float, None
    # list of any of these
    # char_for_None is what character for missing values.  Default for
    # VCF format is ".".

    # Make a list for consistency.
    if type(value) is not type([]):
        value = [value]
    for i in range(len(value)):
        x = value[i]
        if x is None:
            x = char_for_None
        elif type(x) in [type(0), type(0.0)]:
            x = str(x)
        value[i] = x
    x = ",".join(value)
    return x


def identify_caller(vcf_file):
    # Returns:
    # "MuTect"
    # "VarScan2"
    # "mutationSeq"  # JointSNVMix
    # "SomaticSniper"
    # "Strelka"
    # None        If I can't figure out the caller.
    from genomicode import filelib

    # Use the header lines to figure out the format.
    header_lines = []
    for line in filelib.openfh(vcf_file):
        if not line.startswith("##"):
            break
        header_lines.append(line)

    # MuTect
    # ##GATKCommandLine=<ID=MuTect,Version=3.1-0-g72492bb,Date="Sun
    x = [x for x in header_lines if x.startswith("##GATKCommandLine")]
    x = [x for x in x if x.find("MuTect") >= 0]
    if x:
        return "MuTect"
    # VarScan2
    # ##source=VarScan2
    x = [x for x in header_lines if x.startswith("##source=VarScan2")]
    if x:
        return "VarScan2"
    # mutationSeq
    # ##source=mutationSeq_4.3.6
    x = [x for x in header_lines if x.startswith("##source=mutationSeq")]
    if x:
        return "mutationSeq"
    # SomaticSniper
    # ##FORMAT=<ID=SSC,Number=1,Type=Integer,Description="Somatic Score">
    x = [x for x in header_lines if x.startswith("##FORMAT=<ID=SSC")]
    x = [x for x in header_lines if x.find("Somatic Score") >= 0]
    if x:
        return "SomaticSniper"
    # Strelka
    # ##source=strelka
    x = [x for x in header_lines if x.startswith("##source=strelka")]
    if x:
        return "Strelka"
    return None



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
