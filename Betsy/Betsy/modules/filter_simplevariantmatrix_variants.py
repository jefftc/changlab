from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import filelib
        from genomicode import SimpleVariantMatrix
        from Betsy import module_utils as mlib

        simplematrix_file = in_data.identifier
        filelib.assert_exists_nz(simplematrix_file)
        metadata = {}

        x = mlib.get_user_option(
            user_options, "nonsynonymous_and_stopgain_only",
            allowed_values=["no", "yes"])
        nonsynonymous_and_stopgain_only = (x == "yes")

        x = mlib.get_user_option(
            user_options, "sift_polyphen_damaging",
            allowed_values=["no", "yes"])
        sift_polyphen_damaging = (x == "yes")

        min_coverage_in_every_sample = None
        min_callers_in_every_sample = None
        min_callers_in_any_sample = None
        min_gene_expression_in_every_sample = None
        x = mlib.get_user_option(
            user_options, "min_coverage_in_every_sample", type=int)
        if x != "":
            min_coverage_in_every_sample = x
        x = mlib.get_user_option(
            user_options, "min_callers_in_every_sample", type=int)
        if x != "":
            min_callers_in_every_sample = x
        x = mlib.get_user_option(
            user_options, "min_callers_in_any_sample", type=int)
        if x != "":
            min_callers_in_any_sample = x
        x = mlib.get_user_option(
            user_options, "min_gene_expression_in_every_sample", type=float)
        if x != "":
            min_gene_expression_in_every_sample = x
            
        assert not (min_callers_in_every_sample and min_callers_in_any_sample)
        assert nonsynonymous_and_stopgain_only or \
               sift_polyphen_damaging or \
               min_callers_in_every_sample or \
               min_callers_in_any_sample or \
               min_gene_expression_in_every_sample or \
               min_coverage_in_every_sample, \
               "No filters"

        MATRIX = SimpleVariantMatrix.read_as_am(simplematrix_file)

        commands = []
        #in_attrs = in_data.data.attributes
        if nonsynonymous_and_stopgain_only:
            # Actually, just look into the file instead.
            #assert in_attrs["annotated"] == "yes"
            MATRIX = filter_nonsynonymous(MATRIX)
            commands.append("Keep only nonsynonymous and stopgain variants.")
        if sift_polyphen_damaging:
            MATRIX = filter_sift_polyphen_damaging(MATRIX)
            commands.append("Keep only if predicted to be damaging by "
                            "SIFT or Polyphen2.")
        if min_coverage_in_every_sample is not None:
            MATRIX = filter_min_coverage_in_every_sample(
                MATRIX, min_coverage_in_every_sample)
            commands.append(
                "Keep only variants with coverage >= %d "
                "in every sample." % min_coverage_in_every_sample)
        if min_callers_in_every_sample is not None:
            MATRIX = filter_min_callers_in_every_sample(
                MATRIX, min_callers_in_every_sample)
            commands.append(
                "Keep only variants called with >= %d callers "
                "in every sample." % min_callers_in_every_sample)
        if min_callers_in_any_sample is not None:
            MATRIX = filter_min_callers_in_any_sample(
                MATRIX, min_callers_in_any_sample)
            commands.append(
                "Keep only variants called with >= %d callers "
                "in at least one sample." % min_callers_in_any_sample)
        if min_gene_expression_in_every_sample is not None:
            # Actually, just look into the file instead.
            #assert in_attrs["with_gxp"] == "yes"
            MATRIX = filter_min_gene_expression_in_every_sample(
                MATRIX, min_gene_expression_in_every_sample)
            commands.append(
                "Keep only variants with gene expression >= %g "
                "in every sample." % min_gene_expression_in_every_sample)
        metadata["commands"] = commands
        
        SimpleVariantMatrix.write_from_am(out_filename, MATRIX)

        return metadata
            
    
    def name_outfile(self, antecedents, user_options):
        return "calls.txt"


def filter_nonsynonymous(MATRIX):
    # Filter out synonymous variants.
    from genomicode import AnnotationMatrix
    
    # Make sure annotated with Annovar.
    HEADER = "Annovar______ExonicFunc.refGene"
    assert HEADER in MATRIX.headers, "Missing: ExonicFunc.refGene"
    exonic_func = MATRIX[HEADER]
    I_keep = []
    for i, efunc in enumerate(exonic_func):
        assert efunc in [
            "", "nonsynonymous SNV", "synonymous SNV",
            "stopgain", "stoploss",
            "frameshift substitution", "nonframeshift substitution",
            "unknown"], \
            "Unknown exonic_func: %s" % efunc
        if efunc in [
            "nonsynonymous SNV", "stopgain", "stoploss",
            "frameshift substitution"]:
            I_keep.append(i)
    x = AnnotationMatrix.rowslice(MATRIX, I_keep)
    return x


def filter_sift_polyphen_damaging(MATRIX):
    from genomicode import AnnotationMatrix
    
    x = [x for x in MATRIX.headers if x.endswith("SIFT_pred")]
    assert len(x) == 1
    SIFT_pred = MATRIX[x[0]]
    x = [x for x in MATRIX.headers if x.endswith("Polyphen2_HDIV_pred")]
    assert len(x) == 1
    hdiv_pred = MATRIX[x[0]]
    x = [x for x in MATRIX.headers if x.endswith("Polyphen2_HVAR_pred")]
    assert len(x) == 1
    hvar_pred = MATRIX[x[0]]
    
    I_keep = []
    for i, (sift, hdiv, hvar) in enumerate(
        zip(SIFT_pred, hdiv_pred, hvar_pred)):
        if sift == "D" and hdiv in ["D", "P"] and hvar in ["D", "P"]:
            I_keep.append(i)
    x = AnnotationMatrix.rowslice(MATRIX, I_keep)
    return x


def filter_min_gene_expression_in_every_sample(MATRIX, gxp):
    # Gene expression >= 1 in all samples.
    from genomicode import AnnotationMatrix

    assert type(gxp) is type(0.0)
    
    x = MATRIX.headers
    x = [x for x in x if x.startswith("Gene Expression")]
    sample_h = x
    assert sample_h, 'Missing: "Gene Expression" columns'

    I_keep = []
    for i in range(MATRIX.num_annots()):
        keep = True
        for h in sample_h:
            if not MATRIX[h][i]:
                keep = False
                break
            # 5.3
            # 0,0.379
            x = MATRIX[h][i]
            x = x.split(",")
            x = [float(x) for x in x]
            x = max(x)
            exp = x
            if exp < gxp:
                keep = False
                break
        if not keep:
            continue
        I_keep.append(i)
        
    x = AnnotationMatrix.rowslice(MATRIX, I_keep)
    return x



def filter_min_coverage_in_every_sample(MATRIX, coverage):
    from genomicode import AnnotationMatrix

    assert type(coverage) is type(0)
    
    x = MATRIX.headers
    x = [x for x in x if x.startswith("Coverage")]
    sample_h = x
    assert sample_h, 'Missing: "Coverage" columns'

    I_keep = []
    for i in range(MATRIX.num_annots()):
        keep = True
        for h in sample_h:
            if not MATRIX[h][i]:
                keep = False
                break
            # Ref/Alt/VAF
            x = MATRIX[h][i]
            x = x.split("/")
            assert len(x) == 3
            cov = int(x[0]) + int(x[1])
            if cov < coverage:
                keep = False
                break
        if keep:
            I_keep.append(i)
    x = AnnotationMatrix.rowslice(MATRIX, I_keep)
    return x


def filter_min_callers_in_every_sample(MATRIX, num_callers):
    from genomicode import AnnotationMatrix

    assert type(num_callers) is type(0)
    
    x = MATRIX.headers
    x = [x for x in x if x.startswith("Num Callers")]
    callers_h = x
    assert callers_h, 'Missing: "Gene Expression" columns'

    I_keep = []
    for i in range(MATRIX.num_annots()):
        keep = True
        for h in callers_h:
            if not MATRIX[h][i]:
                keep = False
                break
            nc = int(MATRIX[h][i])
            if nc < num_callers:
                keep = False
                break
        if keep:
            I_keep.append(i)
        
    x = AnnotationMatrix.rowslice(MATRIX, I_keep)
    return x


def filter_min_callers_in_any_sample(MATRIX, num_callers):
    from genomicode import AnnotationMatrix

    assert type(num_callers) is type(0)
    
    x = MATRIX.headers
    x = [x for x in x if x.startswith("Num Callers")]
    callers_h = x
    assert callers_h, 'Missing: "Gene Expression" columns'

    I_keep = []
    for i in range(MATRIX.num_annots()):
        keep = False
        for h in callers_h:
            if not MATRIX[h][i]:
                continue
            nc = int(MATRIX[h][i])
            if nc >= num_callers:
                keep = True
                break
        if keep:
            I_keep.append(i)
        
    x = AnnotationMatrix.rowslice(MATRIX, I_keep)
    return x
