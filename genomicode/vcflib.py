# Classes:
# VCFFile
# Variant
#
# Functions:
# read
# write
# 
# parse_variant         Describes one variant for a specific sample.
# update_variant        Update a VCF object with Variant object.
# make_coverage_matrix  Return a variant x sample matrix of coverage.
# make_vaf_matrix       Return a variant x sample matrix of allele frequencies.
#
# _safe_int
# _safe_float
# _percent_to_decimal
# _fmt_vcf_value



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
#                   samtools   Platypus     GATK  IACS  NextGene
# num_ref           GENO:RD              GENO:AD        GENO:SGCOUNTREF_F/R
# num_alt           GENO:AD     INFO:TR  GENO:AD        GENO:SGCOUNTALT_F/R
# total_reads       GENO:DP     INFO:TC  GENO:DP  calc  GENO:DP
# vaf               GENO:FREQ
# call              GENO:GT     GENO:GT  GENO:GT   ???
# * GENO:AD    For GATK, includes REF.  REF,ALT1,ALT2 etc.
#   GENO:FREQ  For samtools, is a percent, e.g. 100%.
#   GENO:FREQ  May be "." or blank.
# * If DP is missing, then try SDP.
# * Sometimes DP is missing for samtools.  Can use INFO:ADP if
#   only 1 sample.
# * NextGene
#   Comma means different (simulatenous) interpretations of
#   alternate alignments.  Allele frequency can be either.
#   SGCOUNTREF_F  3277
#   SGCOUNTREF_R  2675
#   SGCOUNTALT_F  1500,522
#   SGCOUNTALT_R  1761,515
#   AF            0.353,0.112
#   SGACOV        3261,1037   Sum of ALT reads.
#   DP            9232



class VCFFile:
    def __init__(self, matrix, samples, more_info, genotypes):
        # matrix     AnnotationMatrix
        # samples    list of sample names
        # more_info  list of dicts.  from INFO column.
        # genotypes  dict of sample -> list of dicts.  from genotype columns.
        import copy
        self.matrix = copy.deepcopy(matrix)
        self.samples = samples[:]
        self.more_info = copy.deepcopy(more_info)
        self.genotypes = copy.deepcopy(genotypes)
    def num_variants(self):
        return self.matrix.num_annots()
        


class Variant:
    # Holds structured information for each variant.
    #
    # num_alt_alleles
    # num_ref          Can be None if data missing.
    # num_alt          Can be list.  Each can be None.
    # total_reads      Can be list.  Each can be None.
    # vaf              Can be list.  Each can be None.
    # call             String.
    def __init__(self, num_alt_alleles, num_ref, num_alt, total_reads, vaf,
                 call):
        assert num_alt_alleles >= 1
        #if num_alt_alleles > 1:
        #    assert len(num_alt) == num_alt_alleles
        #    assert len(total_reads) == num_alt_alleles
        #    assert len(vaf) == num_alt_alleles
        self.num_alt_alleles = num_alt_alleles
        self.num_ref = num_ref
        self.num_alt = num_alt
        self.total_reads = total_reads
        self.vaf = vaf
        self.call = call


def read(filename):
    import AnnotationMatrix

    info_header = "INFO"
    format_header = "FORMAT"

    matrix = AnnotationMatrix.read(filename, header_char="##")

    # Should start with headers, and then samples after that.
    vcf_headers = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
        "INFO", "FORMAT"]
    assert matrix.headers[:len(vcf_headers)] == vcf_headers
    samples = matrix.headers[len(vcf_headers):]
    assert samples

    # Parse out the INFO line.
    # Format: BRF=0.89;FR=1.0000;HP=2;HapScore=1;...
    annots = matrix.header2annots[info_header]
    # Make a list of dictionaries.
    info_dicts = [_parse_info_dict(x) for x in annots]

    # Parse out the genotype fields.
    format_strings = matrix.header2annots[format_header]
    genotypes = {}  # sample -> list of dicts
    for sample in samples:
        genotype_strings = matrix[sample]
        geno_dicts = [
            _parse_genotype_dict(fs, gs)
            for (fs, gs) in zip(format_strings, genotype_strings)]
        genotypes[sample] = geno_dicts

    return VCFFile(matrix, samples, info_dicts, genotypes)


def write(handle_or_file, vcf):
    import AnnotationMatrix
    AnnotationMatrix.write(handle_or_file, vcf.matrix)


def parse_variant(vcf, variant_num, sample):
    # Return an Variant object describing a variant for a sample.
    assert sample in vcf.samples, "Unknown sample: %s" % sample
    assert variant_num >= 0
    assert variant_num < vcf.matrix.num_annots()

    info_dict = vcf.more_info[variant_num]
    geno_dict = vcf.genotypes[sample][variant_num]

    # Figure out the number of ALT alleles.
    alt_alleles = vcf.matrix["ALT"][variant_num]
    assert alt_alleles
    x = alt_alleles.split(",")
    num_alt_alleles = len(x)

    num_ref = None
    num_alt = None
    total_reads = None
    vaf = None
    call = None

    # For debugging.
    x1 = vcf.matrix["#CHROM"][variant_num]
    x2 = vcf.matrix["POS"][variant_num]
    pos_str = "%s %s" % (x1, x2)
    
    if "RD" in geno_dict:
        num_ref = _safe_int(geno_dict["RD"])
    if "SGCOUNTREF_F" in geno_dict:
        x1 = _safe_int(geno_dict["SGCOUNTREF_F"])
        x2 = _safe_int(geno_dict["SGCOUNTREF_R"])
        num_ref = _safe_add(x1, x2)
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
        if "RD" not in geno_dict:
            # No RD means both AD and RD are merged together.
            if len(x) == num_alt_alleles+1:
                x = x[1:]
        x = [_safe_int(x) for x in x]
        #assert len(x) == num_alt_alleles, "Bad AD: %s %d %s" % (
        #    pos_str, num_alt_alleles, geno_dict["AD"])
        num_alt = x
    if "DP" in geno_dict:
        x = geno_dict["DP"].split(",")
        x = [_safe_int(x) for x in x]
        total_reads = x
    if "FREQ" in geno_dict:
        x = geno_dict["FREQ"].split(",")
        x = [_percent_to_decimal(x) for x in x]
        x = [_safe_float(x) for x in x]
        vaf = x
    if "GT" in geno_dict:
        call = geno_dict["GT"]

    if num_alt is None and "TR" in info_dict:
        x = info_dict["TR"].split(",")
        x = [_safe_int(x) for x in x]
        assert len(x) == num_alt_alleles
        num_alt = x
    if total_reads is None and "TC" in info_dict:
        x = info_dict["TC"].split(",")
        x = [_safe_int(x) for x in x]
        assert len(x) == num_alt_alleles
        total_reads = x

    if num_ref is None and total_reads is not None and num_alt is not None:
        # Possibilities:
        # 1.  1 total reads, 1 alt.
        # 2.  1 total reads, multiple alts.
        #     Different alternative alleles.  Sum them up.
        # 3.  multiple total reads, multiple alts.
        #     should have same number, and num_ref calculated should
        #     be the same.
        
        # Case 1.
        if len(total_reads) == 1 and len(num_alt) == 1:
            if total_reads[0] is not None and num_alt[0] is not None:
                num_ref = total_reads[0] - num_alt[0]
        # Case 2.
        elif len(total_reads) == 1 and len(num_alt) > 1:
            x = [x for x in num_alt if x is not None]
            num_ref = total_reads[0] - sum(x)
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

    if vaf is None and num_ref is not None and num_alt is not None:
        calc_v = [None] * len(num_alt)
        for i in range(len(num_alt)):
            if num_alt[i] is None:
                continue
            total = num_ref + num_alt[i]
            if total:
                x = float(num_alt[i])/total
                calc_v[i] = x
        vaf = calc_v
        
    # Convert lists to numbers.
    if type(num_alt) is type([]) and len(num_alt) == 1:
        num_alt = num_alt[0]
    if type(total_reads) is type([]) and len(total_reads) == 1:
        total_reads = total_reads[0]
    if type(vaf) is type([]) and len(vaf) == 1:
        vaf = vaf[0]

    return Variant(
        num_alt_alleles, num_ref, num_alt, total_reads, vaf, call)

    
def update_variant(vcf, variant, variant_num, sample):
    # Update the VCF with the information from the Info object.
    # Updates in place.  Will change the vcf variable.
    import copy

    assert variant_num < vcf.matrix.num_annots()
    ID = vcf.more_info[variant_num]
    GD = vcf.genotypes[sample][variant_num]

    # Figure out what kind of file it is.
    if "RD" in GD and "AD" in GD and "DP" in GD and "FREQ" in GD:
        # samtools
        GD["RD"] = _fmt_vcf_value(variant.num_ref)
        GD["AD"] = _fmt_vcf_value(variant.num_alt)
        GD["DP"] = _fmt_vcf_value(variant.total_reads)
        # Convert FREQ to percent.
        x = variant.vaf
        if type(x) is not type([]):
            x = [x]
        x = [x*100.0 for x in x]
        x = ["%s%%" % x for x in x]
        GD["FREQ"] = _fmt_vcf_value(x)
        GD["GT"] = variant.call
    elif "TR" in ID and "TC" in ID:
        # Platypus
        ID["TR"] = _fmt_vcf_value(variant.num_ref)
        ID["TC"] = _fmt_vcf_value(variant.total_reads)
        GD["GT"] = variant.call
    elif "AD" in GD and "DP" in GD and not "RD" in GD and not "FREQ" in GD:
        # GATK
        GD["RD"] = _fmt_vcf_value(variant.num_ref)
        GD["AD"] = _fmt_vcf_value(variant.num_alt)
        GD["DP"] = _fmt_vcf_value(variant.total_reads)
        GD["GT"] = variant.call
    elif "SGCOUNTREF_F" in GD:
        # NextGene
        raise NotImplementedError, "NextGene not implemented yet."
    else:
        raise AssertionError, "Unknown VCF format."

    # Change the values in the MATRIX based on ID and GD.
    annots = vcf.matrix["INFO"]
    x = annots[variant_num]
    keyvalues = x.split(";")
    for i, kv in enumerate(keyvalues):
        x = kv.split("=")
        if len(x) == 1:
            continue
        assert len(x) == 2
        key, value = x
        if key in ID:
            value = ID[key]
        keyvalues[i] = "%s=%s" % (key, value)
    x = ";".join(keyvalues)
    annots[variant_num] = x

    annots = vcf.matrix["FORMAT"]
    x = annots[variant_num]
    keys = x.split(":")
    annots = vcf.matrix[sample]
    x = annots[variant_num]
    values = x.split(":")
    for i, key in enumerate(keys):
        if key in GD:
            values[i] = GD[key]
    x = ":".join(values)
    annots[variant_num] = x

    #return vcf


def make_coverage_matrix(vcf, samples=None):
    # Each element in the matrix is an integer or None.
    if samples is None:
        samples = vcf.samples
    for s in samples:
        assert s in vcf.samples

    num_variants = vcf.num_variants()
    num_samples = len(samples)

    matrix = []
    for i in range(num_variants):
        row = []
        for s in samples:
            var = parse_variant(vcf, i, s)
            x = var.total_reads
            assert type(x) in [type(None), type(0)]
            row.append(x)
        matrix.append(row)
    return matrix
            
    
def make_vaf_matrix(vcf, samples=None):
    # Each element in the matrix is a float or None.
    if samples is None:
        samples = vcf.samples
    for s in samples:
        assert s in vcf.samples

    num_variants = vcf.num_variants()
    num_samples = len(samples)

    matrix = []
    for i in range(num_variants):
        row = []
        for s in samples:
            var = parse_variant(vcf, i, s)
            x = var.vaf
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

def _percent_to_decimal(x):
    # If the VAF is provided as a percent, convert it to decimal.
    # e.g. 14.29% -> 0.1429
    if type(x) != type(""):
        return x
    if not x.endswith("%"):
        return x
    x = x[:-1]
    return float(x) / 100


def _fmt_vcf_value(value):
    # value can be:
    # string, int, float, None
    # list of any of these
    
    # Make a list for consistency.
    if type(value) is not type([]):
        value = [value]
    for i in range(len(value)):
        x = value[i]
        if x is None:
            x = "."
        elif type(x) in [type(0), type(0.0)]:
            x = str(x)
        value[i] = x
    x = ",".join(value)
    return x


def _parse_info_dict(info_str):
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
        d[key] = value
    return d


def _parse_genotype_dict(format_str, genotype_str):
    names = format_str.split(":")
    values = genotype_str.split(":")
    assert len(names) == len(values)
    d = {}
    for (n, v) in zip(names, values):
        d[n] = v
    return d


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
