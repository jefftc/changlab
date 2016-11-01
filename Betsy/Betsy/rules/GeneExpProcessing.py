# Merging:
# convert_imput_merge_rma        TODO: Don't separate out RMA
# merge_two_classes_rma
# merge_two_classes              Merge class0 and class1 to class0,class1.
# normalize_samples_with_combat
# normalize_samples_with_dwd
# normalize_samples_with_shiftscale

# Need:
# merge_class0_class1_signal_signal
# merge_class0_class1_signal_classlabel
# merge_class0_class1_classlabel_classlabel

from Betsy.bie3 import *
import BasicDataTypes as BDT
YESNO = BDT.YESNO



SimpleLabelFile = DataType(
    "SimpleLabelFile",
    # Should have headers:
    # Sample   Class
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "preprocess", BDT.PREPROCESS, "unknown", "unknown",
        help="preprocess method"),
    help="A simple tab-delimited format with labels for each sample.",
    )

ClassLabelFile = DataType(
    "ClassLabelFile",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "preprocess", BDT.PREPROCESS, "unknown", "unknown",
        help="preprocess method"),
    help="CLS format file with categorical labels.",
    )


IntensityPlot = DataType(
    'IntensityPlot',
    AttributeDef(
        'contents', BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "preprocess", BDT.ANY_PREPROCESS, "unknown", "unknown",
        help="preprocess method"),
    help="Intensity plot file",
    )


UnprocessedSignalFile = DataType(
    "UnprocessedSignalFile",
    AttributeDef(
        "format", ["unknown", "tdf", "pcl", "gct", "res", "jeffs"],
        "unknown", "tdf", help="file format"),
    
    # Properties of the data.
    AttributeDef(
        #"preprocess", PREPROCESS1, "unknown", 'unknown',
        "preprocess", BDT.PREPROCESS, "unknown", "unknown",
        help="preprocess method"),
    AttributeDef(
        "logged", ["unknown", "no", "yes"], "unknown", "yes",
        help="logged or not"),
    #AttributeDef(
    #    "filter_and_threshold", YESNO, "no", "no",
    #    help="Whether the values are filtered and thresholded."),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    help="The input SignalFile which care the format,preprocess,"
    "logged,filter_and_threshold,contents.")

_SignalFile_Impute = DataType(
    "_SignalFile_Impute",
    # UnprocessedSignalFile
    AttributeDef(
        #"preprocess", PREPROCESS1, "unknown", 'unknown',
        "preprocess", BDT.PREPROCESS, "unknown", 'unknown',
        help="preprocess method"),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    # _SignalFile_Impute
    #AttributeDef(
    #    "filter_and_threshold", YESNO, "no", "no",
    #    help="Whether the values are filtered and thresholded."),
    AttributeDef(
        "missing_values", ["unknown", "no", "yes"], "unknown", "no",
        help="missing values unknown, yes or no"),
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill","zero_fill", help="missing algorithm"),
    AttributeDef(
        "filter_missing_values", YESNO, "no", "no",
        help="Whether missing values are filtered."),
    help="The SignalFile after SignalFile_Postprocess, care missing_values, "
    "missing_algorithm and filter.")

_SignalFile_Merge = DataType(
    "_SignalFile_Merge",
    # UnprocessedSignalFile
    AttributeDef(
        #"preprocess", PREPROCESS1, "unknown", 'unknown',
        "preprocess", BDT.PREPROCESS, "unknown", 'unknown',
        help="preprocess method"),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    # _SignalFile_Impute
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef(
        "filter_missing_values", YESNO, "no", "no",
        help="Whether missing values are filtered."),
    # _SignalFile_Merge
    #AttributeDef(
    #    "filter_and_threshold", YESNO, "no", "no",
    #    help="Whether the values are filtered and thresholded."),
    AttributeDef(
        "dwd_norm", YESNO, "no", "no", help="dwd normalization"),
    AttributeDef(
        "bfrm_norm", YESNO, "no", "no", help="bfrm normalization"),
    AttributeDef(
        "quantile_norm", YESNO, "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "shiftscale_norm", YESNO, "no", "no",
        help="shiftscale normalization"),
    AttributeDef(
        "combat_norm", YESNO, "no", "no", help="combat normalization"),
    help="The SignalFile after SignalFile_Impute, care dwd_norm,bfrm_norm,"\
          "quantile_norm,shiftscale_norm,combat_norm.",
    )

_SignalFile_Normalize = DataType(
    "_SignalFile_Normalize",
    # UnprocessedSignalFile
    AttributeDef(
        "format", ["tdf", "pcl"], "tdf", "tdf", help="file format"),
    AttributeDef(
        #"preprocess", PREPROCESS1, "unknown", 'unknown',
        "preprocess", BDT.PREPROCESS, "unknown", 'unknown',
        help="preprocess method"),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    #AttributeDef(
    #    "filter_and_threshold", YESNO, "no", "no",
    #    help="Whether the values are filtered and thresholded."),
    # _SignalFile_Impute
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef(
        "filter_missing_values", YESNO, "no", "no",
        help="Whether missing values are filtered."),
    # _SignalFile_Merge
    AttributeDef(
        "dwd_norm", YESNO, "no", "no", help="dwd normalization"),
    AttributeDef(
        "bfrm_norm", YESNO, "no", "no", help="bfrm normalization"),
    AttributeDef(
        "quantile_norm", YESNO, "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "shiftscale_norm", YESNO, "no", "no",
        help="shiftscale normalization"),
    AttributeDef(
        "combat_norm", YESNO, "no", "no", help="combat normalization"),
    # _SignalFile_Normalize
    AttributeDef(
        "gene_center", ["unknown", "no", "mean", "median"],
        "unknown", "no", help="gene center method"),
    AttributeDef(
        "gene_normalize", ["unknown", "no", "variance", "sum_of_squares"],
        "unknown", "no", help="gene normalize method"),
    help="The SignalFile after SignalFile_Merge, care gene_center,"\
          "gene_normalize.",
    )


_SignalFile_Order = DataType(
    "_SignalFile_Order",
    # UnprocessedSignalFile
    AttributeDef(
        #"preprocess", PREPROCESS1, "unknown", 'unknown',
        "preprocess", BDT.PREPROCESS, "unknown", 'unknown',
        help="preprocess method"),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    # _SignalFile_Impute
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef(
        "filter_missing_values", YESNO, "no", "no",
        help="Whether missing values are filtered."),
    # _SignalFile_Merge
    AttributeDef(
        "dwd_norm", YESNO, "no", "no", help="dwd normalization"),
    AttributeDef(
        "bfrm_norm", YESNO, "no", "no", help="bfrm normalization"),
    AttributeDef(
        "quantile_norm", YESNO, "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "shiftscale_norm", YESNO, "no", "no",
        help="shiftscale normalization"),
    AttributeDef(
        "combat_norm", YESNO, "no", "no", help="combat normalization"),
    # _SignalFile_Normalize
    AttributeDef(
        "gene_center", [ "no", "mean", "median"], "no", "no",
        help="gene center method"),
    AttributeDef(
        "gene_normalize", [ "no", "variance", "sum_of_squares"],
        "no", "no", help="gene normalize method"),
    # _SignalFile_Order
    #AttributeDef(
    #    "filter_and_threshold", YESNO, "no", "no",
    #    help="Whether the values are filtered and thresholded."),
    AttributeDef(
        "gene_order", BDT.GENE_ORDER, "none", "none",
        help="gene order method"),
    help="The SignalFile after SignalFile_Normalize, care gene_order.")


_SignalFile_Annotate = DataType(
    "_SignalFile_Annotate",
    # UnprocessedSignalFile
    AttributeDef(
        #"preprocess", PREPROCESS1, "unknown", 'unknown',
        "preprocess", BDT.PREPROCESS, "unknown", 'unknown',
        help="preprocess method"),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    # _SignalFile_Impute
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef(
        "filter_missing_values", YESNO, "no", "no",
        help="Whether missing values are filtered."),
    # _SignalFile_Merge
    AttributeDef(
        "dwd_norm", YESNO, "no", "no", help="dwd normalization"),
    AttributeDef(
        "bfrm_norm", YESNO, "no", "no", help="bfrm normalization"),
    AttributeDef(
        "quantile_norm", YESNO, "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "shiftscale_norm", YESNO, "no", "no",
        help="shiftscale normalization"),
    AttributeDef(
        "combat_norm", YESNO, "no", "no", help="combat normalization"),
    # _SignalFile_Normalize
    #AttributeDef(
    #    "filter_and_threshold", YESNO, "no", "no",
    #    help="Whether the values are filtered and thresholded."),
    AttributeDef(
        "gene_center", ["no", "mean", "median"], "no", "no",
        help="gene center method"),
    AttributeDef(
        "gene_normalize", [ "no", "variance", "sum_of_squares"], "no", "no",
        help="gene normalize method"),
    # _SignalFile_Order
    AttributeDef(
        "gene_order", BDT.GENE_ORDER, "none", "none",
        help="gene order method"),
    # _SignalFile_Annotate
    AttributeDef(
        "annotate", YESNO, "no", "no", help="annotate file or not"),
    AttributeDef(
        "rename_sample", YESNO, "no", "no",
        help="rename sample or not"),
    AttributeDef(
        # Why is the u133A?
        "platform", ["yes", "no", 'u133A'], "no", "no",
        help="add platform or not"),
    #AttributeDef("has_u133A", ["yes", "no"], "no", "no", help="has hg_u133A platform or not"),
    help="The SignalFile after SignalFile_Order, care annotate,"
    "rename_sample,platform.")


_SignalFile_Filter = DataType(
    "_SignalFile_Filter",
    # UnprocessedSignalFile
    AttributeDef(
        #"preprocess", PREPROCESS1, "unknown", 'unknown',
        "preprocess", BDT.PREPROCESS, "unknown", 'unknown',
        help="preprocess method"),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    # Not sure why these are back now.
    AttributeDef("logged", [ "no", "yes"], "yes", "yes", help="logged or not"),
    AttributeDef("format", [ "tdf", "gct"], "tdf", "tdf", help="file format"),
    # _SignalFile_Impute
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef(
        "filter_missing_values", YESNO, "no", "no",
        help="Whether missing values are filtered."),
    # _SignalFile_Merge
    AttributeDef(
        "dwd_norm", YESNO, "no", "no", help="dwd normalization"),
    AttributeDef(
        "bfrm_norm", YESNO, "no", "no", help="bfrm normalization"),
    AttributeDef(
        "quantile_norm", YESNO, "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "shiftscale_norm", YESNO, "no", "no",
        help="shiftscale normalization"),
    AttributeDef(
        "combat_norm", YESNO, "no", "no", help="combat normalization"),
    # _SignalFile_Normalize
    AttributeDef(
        "gene_center", [ "no", "mean", "median"], "no", "no",
        help="gene center method"),
    AttributeDef(
        "gene_normalize", [ "no", "variance", "sum_of_squares"], "no", "no",
        help="gene normalize method"),
    # _SignalFile_Order
    AttributeDef(
        "gene_order", BDT.GENE_ORDER, "none", "none",
        help="gene order method"),
    # _SignalFile_Annotate
    AttributeDef(
        "annotate", YESNO, "no", "no", help="annotate file or not"),
    AttributeDef(
        "rename_sample", YESNO, "no", "no",
        help="rename sample or not"),
    AttributeDef(
        "platform", ["yes", "no", 'u133A'], "no", "no",
        help="add platform or not"),
    # _SignalFile_Filter
    #AttributeDef("has_u133A", ["yes", "no"], "no", "no", help="has hg_u133A platform or not"),
    AttributeDef(
        "num_features", ["yes", "no"], "no", "no",
        help="select a num of features or not"),
    AttributeDef(
        "filter_and_threshold", YESNO, "no", "no",
        help="Whether the values are filtered and thresholded."),
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no", help="method to get unique genes"),
    AttributeDef(
        "duplicate_probe", ["no", "closest_probe", "high_var_probe"],
        "no", "no", help="method to remove duplicated probes"),
    AttributeDef(
        "group_fc", ["yes", "no"], "no", "no",
        help="group fold change or not"),
    help="The SignalFile after SignalFile_Annotate, care num_features,"\
          "unique_genes,duplicate_probe,group_fc,logged,format.")


SignalFile= DataType(
    "SignalFile",
    # UnprocessedSignalFile
    AttributeDef(
        "preprocess", BDT.ANY_PREPROCESS, "unknown", 'any', 
        #"preprocess", PREPROCESS, "unknown", 'unknown',
        help="preprocess method"),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "logged", YESNO, "yes", "yes", help="logged or not"),
    AttributeDef("format", ["tdf", "gct"], "tdf", "tdf", help="file format"),
    # _SignalFile_Impute
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef(
        "filter_missing_values", YESNO, "no", "no",
        help="Whether missing values are filtered."),
    # _SignalFile_Merge
    AttributeDef(
        "dwd_norm", YESNO, "no", "no", help="dwd normalization"),
    AttributeDef(
        "bfrm_norm", YESNO, "no", "no", help="bfrm normalization"),
    AttributeDef(
        "quantile_norm", YESNO, "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "shiftscale_norm", YESNO, "no", "no",
        help="shiftscale normalization"),
    AttributeDef(
        "combat_norm", YESNO, "no", "no", help="combat normalization"),
    # _SignalFile_Normalize
    AttributeDef(
        "gene_center", ["no", "mean", "median"], "no", "no",
        help="gene center method"),
    AttributeDef(
        "gene_normalize", ["no", "variance", "sum_of_squares"], "no", "no",
        help="gene normalize method"),
    # _SignalFile_Order
    AttributeDef(
        "gene_order", BDT.GENE_ORDER, "none", "none",
        help="gene order method"),
    # _SignalFile_Annotate
    AttributeDef(
        "annotate", YESNO, "no", "no", help="annotate file or not"),
    AttributeDef(
        "rename_sample", YESNO, "no", "no",
        help="rename sample or not"),
    AttributeDef(
        "platform", ["yes", "no", 'u133A'], "no", "no",
        help="add platform or not"),
    #AttributeDef("has_u133A", ["yes", "no"], "no", "no", help="has hg_u133A platform or not"),
    # _SignalFile_Filter
    AttributeDef(
        "num_features", ["yes", "no"], "no", "no",
        help="select a num of features or not"),
    AttributeDef(
        "filter_and_threshold", YESNO, "no", "no",
        help="Whether the values are filtered and thresholded."),
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no", help="method to get unique genes"),
    AttributeDef(
        "duplicate_probe", ["no", "closest_probe", "high_var_probe"],
        "no", "no", help="method to remove duplicated probes"),
    AttributeDef(
        "group_fc", ["yes", "no"], "no", "no",
        help="group fold change or not"),
    help="The SignalFile after SignalFile_Filter, the attributes are the "
    "same as SignalFile_Filter.  Should not be specified by user.  "
    "Data provided by user should be provided as UnprocessedSignalFile.")


all_modules = [
##     ModuleNode('preprocess_tcga_affymetrix',Database.TCGAFile,_SignalFile_Postprocess,
##        Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS),
##        Constraint("tumor_only",MUST_BE, 'yes'),
##        Constraint("preprocess",MUST_BE, 'affymetrix'),
##        Consequence("contents", SAME_AS_CONSTRAINT),
##        Consequence('logged',SET_TO, "unknown"),
##        Consequence('predataset', SET_TO, "no"),
##        Consequence('preprocess',SET_TO, "affymetrix"),
##        Consequence('format',SET_TO, "tdf"),
##        help="preprocess tcga rma file, generate to SignalFile_Postprocess"),
##
##    ModuleNode('preprocess_tcga_agilent',Database.TCGAFile,_SignalFile_Postprocess,
##        Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS),
##        Constraint("tumor_only",MUST_BE, 'yes'),
##        Constraint("preprocess",MUST_BE, 'agilent'),
##        Consequence("contents", SAME_AS_CONSTRAINT),
##        Consequence('logged',SET_TO, "unknown"),
##        Consequence('predataset', SET_TO, "no"),
##        Consequence('preprocess',SET_TO, "agilent"),
##        Consequence('format',SET_TO, "tdf"),
##        help="preprocess tcga agilent file, generate to SignalFile_Postprocess"),
    ####postprocess
    ModuleNode(
        "convert_signal_to_tdf",
        UnprocessedSignalFile, UnprocessedSignalFile,
        Constraint(
            "format", CAN_BE_ANY_OF,
            ["unknown", "pcl", "gct", "res", "jeffs"]),
        Consequence("format", SET_TO, "tdf"),
        help="convert SignalFile_Postprocess to tdf format"
        ),
    
    ModuleNode(
        "check_for_log",
        UnprocessedSignalFile, UnprocessedSignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", BASED_ON_DATA, ["yes", "no"]),
        help="check for log in SignalFile_Postprocess"
        ),

    ModuleNode(
        "filter_and_threshold_genes",
        _SignalFile_Filter, _SignalFile_Filter,
        OptionDef(
            "min_threshold", default="",
            help="Set the minimum value to this number, e.g. 20"),
        OptionDef(
            "max_threshold", default="",
            help="Set the maximum value to this number, e.g. 16000"),
        OptionDef(
            "min_fold_change", default="",
            help="Remove genes that do not have at least this fold change, "
            "e.g. 5 (for 5-fold change)"),
        OptionDef(
            "min_delta", default="",
            help="Remove genes that do not have at least this much different "
            "between the maximum and minimum expressing samples, e.g. 100.0"),
        OptionDef(
            "genes_with_highest_var", default="",
            help="Keep just this number of genes with the highest variance, "
            "e.g. 250"),
        Constraint('format', MUST_BE, "tdf"),
        Consequence('format', SAME_AS_CONSTRAINT),
        Constraint('logged', MUST_BE, "yes"),
        #Constraint('logged', CAN_BE_ANY_OF, ["no", "yes"]),
        Consequence('logged', SAME_AS_CONSTRAINT),
        Constraint('filter_and_threshold', MUST_BE, "no"),
        Consequence('filter_and_threshold', SET_TO, 'yes'),
        help="Either threshold (change the value of) or filter genes.  "
        "All options given here should be non-logged values."
        ),
    ModuleNode(
        "log_signal",
        UnprocessedSignalFile, UnprocessedSignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SET_TO, "yes"),
        ),
    # No.  Causes cycles in the network.
    #ModuleNode(
    #    # Unlog RMA values so it can be filtered and thresholded
    #    # (which requires log=no).  If I allow the module to take
    #    # log=yes, will increase pipelines exponentially.
    #    "unlog_rma_signal",
    #    UnprocessedSignalFile, UnprocessedSignalFile,
    #    Constraint("format", MUST_BE, "tdf"),
    #    Consequence("format", SAME_AS_CONSTRAINT),
    #    Constraint("logged", MUST_BE, "yes"),
    #    Consequence("logged", SET_TO, "no"),
    #    Constraint("preprocess", MUST_BE, "rma"),
    #    Consequence("preprocess", SAME_AS_CONSTRAINT),
    #    ),

    #impute
    ModuleNode(
        "convert_postprocess_impute",
        UnprocessedSignalFile, _SignalFile_Impute,
        Constraint("format", MUST_BE, "tdf"),
        #Constraint("logged", MUST_BE, "yes"),
        #Constraint("logged", CAN_BE_ANY_OF, ["no", "yes", "unknown"]),
        
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
        
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("contents", SAME_AS_CONSTRAINT),
        #Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
        
        Consequence("missing_algorithm", SET_TO, 'zero_fill'),
        Consequence("missing_values", SET_TO, "unknown"),
        Consequence("filter_missing_values", SET_TO, 'no'),
        help="convert SignalFile_Postprocess to SignalFile_Impute"),

    ModuleNode(
        "check_for_missing_values",
        _SignalFile_Impute, _SignalFile_Impute,
        Constraint("missing_values", MUST_BE, "unknown"),
        Consequence("missing_values", BASED_ON_DATA, YESNO),
        help="check missing values in SignalFile_Impute"
        ),

    ModuleNode(
        "filter_genes_by_missing_values",
        _SignalFile_Impute, _SignalFile_Impute,
        OptionDef(
            "filter_value", "0.50",
            help="filter by missing values in percentage, etc.(0-1)"),
        Constraint("missing_values", MUST_BE, "yes"),
        Consequence("missing_values", SAME_AS_CONSTRAINT),
        Constraint("filter_missing_values", MUST_BE, "no"),
        Consequence("filter_missing_values", SET_TO, "yes"),
        help="filter genes by missing values in SignalFile_Impute"
        ),
    ModuleNode(
        "fill_missing_with_zeros",
        _SignalFile_Impute, _SignalFile_Impute,
        Constraint("missing_values", MUST_BE, "yes"),
        Consequence("missing_values", SET_TO, "no"),
        Consequence(
            "missing_algorithm", SET_TO, "zero_fill", side_effect=True),
        help="fill missing values in SignalFile_Impute with zeros"
        ),
    ModuleNode(
        "fill_missing_with_median",
        _SignalFile_Impute, _SignalFile_Impute,
        Constraint('missing_algorithm', MUST_BE, 'zero_fill'),
        Constraint('missing_values', MUST_BE, 'yes'),
        Consequence('missing_algorithm',SET_TO, "median_fill",side_effect=True),
        Consequence('missing_values',SET_TO, 'no'),
        help="fill missing values in SignalFile_Impute with median"),
    
    #merge
    ModuleNode(
        "convert_impute_merge",
        _SignalFile_Impute, _SignalFile_Merge,
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS_WOrma),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        #Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ['none', 'zero_fill', 'median_fill']),
        Constraint("missing_values", MUST_BE, "no"),
        
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        #Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        
        Consequence("quantile_norm", SET_TO, 'no'),
        Consequence("dwd_norm", SET_TO, 'no'),
        Consequence("combat_norm", SET_TO, 'no'),
        Consequence("bfrm_norm", SET_TO, 'no'),
        Consequence("shiftscale_norm", SET_TO, 'no'),
        help="convert SignalFile_Impute to SignalFile_Merge"),

    # TODO: Why is this separate?
    ## ModuleNode(
    ##     "convert_impute_merge_rma",
    ##     _SignalFile_Impute, _SignalFile_Merge,
    ##     Constraint("preprocess", MUST_BE, "rma"),
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
    ##     Constraint("filter", CAN_BE_ANY_OF, ["no", "yes"]),
    ##     Constraint("missing_values", MUST_BE, "no"),
    ##     Constraint(
    ##         "missing_algorithm", CAN_BE_ANY_OF,
    ##         ['none', 'zero_fill', 'median_fill']),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence("preprocess", SAME_AS_CONSTRAINT),
    ##     Consequence("predataset", SAME_AS_CONSTRAINT),
    ##     Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
    ##     Consequence("filter", SAME_AS_CONSTRAINT),
    ##     Consequence("quantile_norm", SET_TO, 'yes'),
    ##     Consequence("dwd_norm", SET_TO, 'no'),
    ##     Consequence("combat_norm", SET_TO, 'no'),
    ##     Consequence("bfrm_norm", SET_TO, 'no'),
    ##     Consequence("shiftscale_norm", SET_TO, 'no'),
    ##     help="convert SignalFile_Impute to SignalFile_Merge with "
    ##     "preprocess=rma"),

    ## ModuleNode(  #did not consider the diff_expr case
    ##     "merge_two_classes_rma", [_SignalFile_Merge, _SignalFile_Merge],
    ##     _SignalFile_Merge,
    ##     Constraint("contents", MUST_BE, "class0", 0),
    ##     Constraint("preprocess", MUST_BE, 'rma', 0),
    ##     Constraint("combat_norm", MUST_BE, 'no', 0),
    ##     Constraint("quantile_norm", MUST_BE, 'yes', 0),
    ##     Constraint("dwd_norm", MUST_BE, "no", 0),
    ##     Constraint("bfrm_norm", MUST_BE, "no", 0),
    ##     Constraint("shiftscale_norm", MUST_BE, "no", 0),
    ##     Constraint("preprocess", MUST_BE, 'rma', 1),
    ##     Constraint("contents", MUST_BE, "class1", 1),
    ##     Constraint("combat_norm", MUST_BE, 'no', 1),
    ##     Constraint("quantile_norm", MUST_BE, "yes", 1),
    ##     Constraint("dwd_norm", MUST_BE, "no", 1),
    ##     Constraint("bfrm_norm", MUST_BE, "no", 1),
    ##     Constraint("shiftscale_norm", MUST_BE, "no", 1),
    ##     Consequence("contents", SET_TO, "class0,class1"),
    ##     Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("combat_norm", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("quantile_norm", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("dwd_norm", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 0),
    ##     DefaultAttributesFrom(0),
    ##     DefaultAttributesFrom(1),
    ##     help="merge two classes SignalFile_Merge with preprocess=rma,generate "
    ##     "SignalFile_Merge",
    ##     ),

    ModuleNode(
        "merge_class0_class1_signal_signal",
        [_SignalFile_Merge, _SignalFile_Merge], _SignalFile_Merge,

        # Interesting attributes.
        Constraint("contents", MUST_BE, "class0", 0),
        Constraint("contents", MUST_BE, "class1", 1),
        Consequence("contents", SET_TO, "class0,class1"),

        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS_WOrma, 0),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        
        Constraint("combat_norm", MUST_BE, 'no', 0),
        Constraint("combat_norm", MUST_BE, 'no', 1),
        Consequence("combat_norm", SAME_AS_CONSTRAINT, 0),
        
        Constraint("quantile_norm", MUST_BE, 'no', 0),
        Constraint("quantile_norm", MUST_BE, "no", 1),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT, 0),
        
        Constraint("dwd_norm", MUST_BE, "no", 0),
        Constraint("dwd_norm", MUST_BE, "no", 1),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT, 0),
        
        Constraint("bfrm_norm", MUST_BE, "no", 0),
        Constraint("bfrm_norm", MUST_BE, "no", 1),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 0),
        
        Constraint("shiftscale_norm", MUST_BE, "no", 0),
        Constraint("shiftscale_norm", MUST_BE, "no", 1),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 0),
        DefaultAttributesFrom(0),
        #DefaultAttributesFrom(1),
        help="merge two classes SignalFile_Merge, generate SignalFile_Merge",
        ),

    ModuleNode(
        "merge_class0_class1_signal_classlabel",
        [_SignalFile_Merge, _SignalFile_Merge], ClassLabelFile,

        # Interesting attributes.
        Constraint("contents", MUST_BE, "class0", 0),
        Constraint("contents", MUST_BE, "class1", 1),
        Consequence("contents", SET_TO, "class0,class1"),

        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS_WOrma, 0),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        
        Constraint("combat_norm", MUST_BE, 'no', 0),
        Constraint("combat_norm", MUST_BE, 'no', 1),
        
        Constraint("quantile_norm", MUST_BE, 'no', 0),
        Constraint("quantile_norm", MUST_BE, 'no', 1),
        
        Constraint("dwd_norm", MUST_BE, "no", 0),
        Constraint("dwd_norm", MUST_BE, "no", 1),
        
        Constraint("bfrm_norm", MUST_BE, "no", 0),
        Constraint("bfrm_norm", MUST_BE, "no", 1),
        
        Constraint("shiftscale_norm", MUST_BE, "no", 0),
        Constraint("shiftscale_norm", MUST_BE, "no", 1),
        help="Generate a ClassLabelFile from two SignalFiles.",
        ),

    ## ModuleNode(  #did not consider the diff_expr case
    ##     "merge_two_classes", [_SignalFile_Merge, _SignalFile_Merge],
    ##     _SignalFile_Merge,
    ##     Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS_WOrma),
    ##     Constraint("contents", MUST_BE, "class0", 0),
    ##     Constraint("combat_norm", MUST_BE, 'no', 0),
    ##     Constraint("quantile_norm", MUST_BE, 'no', 0),
    ##     Constraint("dwd_norm", MUST_BE, "no", 0),
    ##     Constraint("bfrm_norm", MUST_BE, "no", 0),
    ##     Constraint("shiftscale_norm", MUST_BE, "no", 0),
    ##     Constraint("contents", MUST_BE, "class1", 1),
    ##     Constraint("combat_norm", MUST_BE, 'no', 1),
    ##     Constraint("quantile_norm", MUST_BE, "no", 1),
    ##     Constraint("dwd_norm", MUST_BE, "no", 1),
    ##     Constraint("bfrm_norm", MUST_BE, "no", 1),
    ##     Constraint("shiftscale_norm", MUST_BE, "no", 1),
    ##     Consequence("contents", SET_TO, "class0,class1"),
    ##     Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("combat_norm", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("quantile_norm", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("dwd_norm", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 0),
    ##     DefaultAttributesFrom(0),
    ##     DefaultAttributesFrom(1),
    ##     help="merge two classes SignalFile_Merge,generate SignalFile_Merge",
    ##     ),

    ModuleNode(
        "normalize_samples_with_quantile",
        _SignalFile_Merge, _SignalFile_Merge,
        Constraint("combat_norm", MUST_BE, "no"),
        Constraint("shiftscale_norm", MUST_BE, "no"),
        Constraint("bfrm_norm", MUST_BE, "no"),
        Constraint("dwd_norm", MUST_BE, "no"),
        Constraint("quantile_norm", MUST_BE, "no"),
        Consequence("quantile_norm", SET_TO, "yes"),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        help="nommalize SignalFile_Merge with quantile method",
        ),
    ModuleNode(
        "normalize_samples_with_bfrm",
        _SignalFile_Merge, _SignalFile_Merge,
        OptionDef("num_factors", 1, help="num factors for bfrm normalization"),
        Constraint('bfrm_norm', MUST_BE, "no"),
        Constraint('combat_norm', MUST_BE, "no"),
        Constraint('shiftscale_norm', MUST_BE, "no"),
        Constraint('dwd_norm', MUST_BE, "no"),
        Consequence('bfrm_norm',SET_TO, 'yes'),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        help="nommalize SignalFile_Merge with bfrm method",
        ),

    ModuleNode(
        "normalize_samples_with_combat",
        [ClassLabelFile, _SignalFile_Merge], _SignalFile_Merge,
        Constraint("combat_norm", MUST_BE, "no", 1),
        Constraint('shiftscale_norm', MUST_BE, "no", 1),
        Constraint('bfrm_norm', MUST_BE, "no", 1),
        Constraint('dwd_norm', MUST_BE, "no", 1),
        Consequence("combat_norm", SET_TO, "yes"),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 1),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        DefaultAttributesFrom(1),
        help="nommalize SignalFile_Merge with combat method"),

    ModuleNode(
        "normalize_samples_with_dwd",
        [ClassLabelFile, _SignalFile_Merge], _SignalFile_Merge,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("dwd_norm", MUST_BE, "no", 1),
        Consequence("dwd_norm", SET_TO, "yes"),
        Constraint("bfrm_norm", MUST_BE, "no", 1),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 1),
        Constraint("combat_norm", MUST_BE, "no", 1),
        Consequence("combat_norm", SAME_AS_CONSTRAINT, 1),
        Constraint("shiftscale_norm", MUST_BE, "no", 1),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 1),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS, 1),
        #Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
        DefaultAttributesFrom(1),
        help="nommalize SignalFile_Merge with dwd method"),

    ModuleNode(
        "normalize_samples_with_shiftscale",
        [ClassLabelFile, _SignalFile_Merge], _SignalFile_Merge,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("shiftscale_norm", MUST_BE, "no", 1),
        Consequence("shiftscale_norm", SET_TO, "yes"),
        Constraint("bfrm_norm", MUST_BE, "no", 1),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 1),
        Constraint("combat_norm", MUST_BE, "no", 1),
        Consequence("combat_norm", SAME_AS_CONSTRAINT, 1),
        Constraint("dwd_norm", MUST_BE, "no", 1),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT, 1),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS, 1),
        #Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
        DefaultAttributesFrom(1),
        help="nommalize SignalFile_Merge with shiftscale method"),
    
    ###normalize
    ModuleNode(
        "convert_merge_normalize",
        _SignalFile_Merge, _SignalFile_Normalize,
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ['none', 'zero_fill', 'median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        #Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("combat_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("dwd_norm", CAN_BE_ANY_OF, YESNO),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SET_TO, 'unknown'),
        Consequence("gene_normalize", SET_TO, 'unknown'),
        Consequence("format", SET_TO, 'tdf'),
        help="convert SignalFile_Merge to SignalFile_Normalize"
        ),
    ModuleNode(
        "check_gene_center",
        _SignalFile_Normalize, _SignalFile_Normalize,
        Constraint("format", MUST_BE, 'tdf'),
        Constraint("gene_center", MUST_BE, "unknown"),
        Constraint("gene_normalize", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("gene_center",BASED_ON_DATA, ["no", "mean", "median"]),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        help="check SignalFile_Normalize if it is gene center or not"),

    ModuleNode(
        "check_gene_normalize",
        _SignalFile_Normalize, _SignalFile_Normalize,
        Constraint("format", MUST_BE, 'tdf'),
        Constraint("gene_normalize", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence(
            "gene_normalize", BASED_ON_DATA,
            ["no", "variance", "sum_of_squares"]),
        help="check SignalFile_Normalize if gene normalize or nto"),

    ModuleNode(
        "convert_signal_to_pcl",
        _SignalFile_Normalize, _SignalFile_Normalize,
        Constraint("format", MUST_BE, 'tdf'),
        Consequence("format", SET_TO, 'pcl'),
        help="convert SignalFile_Normalize from tdf format to pcl format"),

    ModuleNode(
        "center_genes",
        _SignalFile_Normalize, _SignalFile_Normalize,
        Constraint("format", MUST_BE, "pcl"),
        Constraint("gene_center", MUST_BE, "no"),
        Constraint("gene_normalize", MUST_BE, 'unknown'),
        Consequence("format", SET_TO, "tdf"),
        Consequence("gene_center", SET_TO_ONE_OF, ["mean", "median"]),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        help="center genes in SignalFile_Normalize"),

    ModuleNode(
        "normalize_genes",
        _SignalFile_Normalize, _SignalFile_Normalize,
        Constraint("format", MUST_BE, "pcl"),
        Constraint("gene_normalize", MUST_BE, "no"),
        Consequence("format", SET_TO, "tdf"),
        Consequence(
            "gene_normalize", SET_TO_ONE_OF, ["variance", "sum_of_squares"]),
        help="normalize genes in SignalFile_Normalize"),
    
    ##Order
    ModuleNode(
        "convert_normalize_order",
        _SignalFile_Normalize, _SignalFile_Order,
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ['none', 'zero_fill', 'median_fill']),
        Constraint("format", MUST_BE, 'tdf'),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        #Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("combat_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("dwd_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median", "mean"]),
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ["no", "variance", "sum_of_squares"]),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("gene_order", SET_TO, 'none'),
        help="convert SignalFile_Normalize to SignalFile_Order"
        ),


    ModuleNode(
        "rank_genes_by_class_neighbors",
        [_SignalFile_Order, ClassLabelFile], BDT.GeneListFile,
        OptionDef(
            "cn_num_neighbors", 50,
            help='number of neighbors for class neighbors method'),
        OptionDef(
            "cn_num_perm", 100,
            help='number of permutation for class neighbors method'),
        OptionDef(
            "cn_user_pval", 0.5,
            help='number of user p value for class neighbors method'),
        OptionDef(
            "cn_min_threshold", 10,
            help='min threshold for class neighbors method'),
        OptionDef(
            "cn_max_threshold", 16000,
            help='max threshold for class neighbors method'),
        OptionDef(
            "cn_min_folddiff", 5,
            help='min fold diff for class neighbors method'),
        OptionDef(
            "cn_abs_diff", 50, help='abs diff for class neighbors method'),

        Constraint("gene_order", MUST_BE, "none", 0),
        #Constraint("gene_order", SAME_AS, 0, 1),
        Consequence("gene_order", SET_TO, "class_neighbors"),

        Constraint("contents", MUST_BE, "class0,class1", 0),
        Constraint("contents", SAME_AS, 0, 1),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 0),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        
        #Constraint("gene_order", MUST_BE, "no", 1),
        #Consequence("cn_mean_or_median", SET_TO_ONE_OF, ['mean', 'median']),
        #Consequence("cn_ttest_or_snr", SET_TO_ONE_OF, ['t_test', 'snr']),
        #Consequence("cn_filter_data", SET_TO_ONE_OF, ['yes', 'no']),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        help="rank the genes in SignalFile_Order by class neighbors method"
        ),

    ModuleNode(
        "rank_genes_by_sample_ttest",
        [_SignalFile_Order, ClassLabelFile], BDT.GeneListFile,
        OptionDef(
            "gene_select_threshold", 0.05, help="threshold for sample ttest"),
        
        # Should make this an OptionDef.
        Constraint("gene_order", MUST_BE, "none", 0),
        Consequence("gene_order", SET_TO_ONE_OF, ["ttest_p", "ttest_fdr"]),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        help="rank the genes in SignalFile_Order by ttest method"
        ),

    ModuleNode(
        "reorder_genes",
        [_SignalFile_Order, BDT.GeneListFile], _SignalFile_Order,
        
        Constraint("gene_order", MUST_BE, "none", 0),
        Constraint(
            "gene_order", CAN_BE_ANY_OF, BDT.GENE_ORDER_not_none, 1),
        Consequence("gene_order", SAME_AS_CONSTRAINT, 1),

        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 0),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT),

        Constraint("gene_center", CAN_BE_ANY_OF, ['median', 'mean', 'no'], 0),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ['sum_of_squares', 'variance', 'no'], 0),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        #DefaultAttributesFrom(0),
        help="rank the genes in SignalFile_Order by genes in GeneListFile"
        ),
    
    ## ModuleNode(
    ##     "reorder_genes_with_diff",
    ##     [_SignalFile_Order, BDT.GeneListFile], _SignalFile_Order,
        
    ##     Constraint("gene_order", MUST_BE, "no", 0),
    ##     Constraint(
    ##         "gene_order", CAN_BE_ANY_OF,
    ##         ['diff_ttest', 'diff_sam', 'diff_ebayes', 'diff_fold_change'], 1),
    ##     Consequence("gene_order", SAME_AS_CONSTRAINT, 1),
        
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
    ##     Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 0),
    ##     Constraint("gene_center", CAN_BE_ANY_OF, ['median', 'mean', 'no'], 0),
    ##     Constraint(
    ##         "gene_normalize", CAN_BE_ANY_OF,
    ##         ['sum_of_squares', 'variance', 'no'], 0),
        
    ##     #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 1),
    ##     Consequence("contents", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("gene_center", SAME_AS_CONSTRAINT, 0),
    ##     Consequence("gene_normalize", SAME_AS_CONSTRAINT, 0),
    ##     DefaultAttributesFrom(0),
    ##     help="reorder SignalFile_Order with genes in GeneListFile"
    ##     ),
    
    ##Annotate
    ModuleNode(
        "convert_order_annotate",
        _SignalFile_Order, _SignalFile_Annotate,
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ['none', 'zero_fill', 'median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        #Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
        
        Constraint("quantile_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("combat_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("dwd_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median", "mean"]),
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ["no", "variance", "sum_of_squares"]),
        Constraint("gene_order", CAN_BE_ANY_OF, BDT.GENE_ORDER),
        
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("gene_order", SAME_AS_CONSTRAINT),
        Consequence("annotate", SET_TO, "no"),
        Consequence("rename_sample", SET_TO, "no"),
        Consequence("platform", SET_TO, "no"),
        help="convert SignalFile_Order to SignalFile_Annotate"
        ),
    
    ModuleNode(
        'annotate_probes',
        _SignalFile_Annotate, _SignalFile_Annotate,

        OptionDef(
            "illu_clm", '', help="CLM file for mapping file to sample names"),
        
        Constraint("annotate", MUST_BE, "no"),
        Consequence("annotate", SET_TO, "yes"),
        Constraint("platform", MUST_BE, "no"),
        Consequence("platform", SAME_AS_CONSTRAINT),

        help="annotate SignalFile_Annotate",
        ),

    ModuleNode(
       "relabel_samples",
        [BDT.RenameFile, _SignalFile_Annotate], _SignalFile_Annotate,
        Constraint("rename_sample", MUST_BE, "no", 1),
        Constraint("labels_from", CAN_BE_ANY_OF, ["title", "description"], 0),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("rename_sample", SET_TO, "yes"),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint('annotate', MUST_BE, 'no', 1),
        Consequence("annotate", SAME_AS_CONSTRAINT, 1),
        Constraint('platform', MUST_BE, 'no', 1),
        Consequence("platform", SAME_AS_CONSTRAINT, 1),
        DefaultAttributesFrom(1),
        help="relabel the sample names in SignalFile_Annotate given "
        "RenameFile.",
       ),
    ModuleNode(
         'add_crossplatform_probeid',
         _SignalFile_Annotate, _SignalFile_Annotate,
         OptionDef(
            "platform_name",
            help="given the new platform name to add to the file"),
         Constraint("platform", MUST_BE, "no"),
         Consequence("platform", SET_TO, "yes"),
         help="add a cross platform to SignalFile_Annotate"),

    ModuleNode(
         'add_U133A_probeid',
         _SignalFile_Annotate, _SignalFile_Annotate,
         Constraint("platform", MUST_BE, "no"),
         Consequence("platform", SET_TO, "u133A"),
         help="add a hg_u133A platform to SignalFile_Annotate"),
    #Filter
    ModuleNode(
        "convert_annotate_filter",
        _SignalFile_Annotate, _SignalFile_Filter,
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
        #Constraint("logged", CAN_BE_ANY_OF, ["no", "yes", "unknown"]),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ['none', 'zero_fill', 'median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("filter_and_threshold", SET_TO, "no"),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("combat_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("dwd_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median", "mean"]),
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ["no", "variance", "sum_of_squares"]),
        Constraint("gene_order", CAN_BE_ANY_OF, BDT.GENE_ORDER),
        Constraint("annotate", CAN_BE_ANY_OF, YESNO),
        Constraint("rename_sample", CAN_BE_ANY_OF, YESNO),
        Constraint("platform", CAN_BE_ANY_OF, ["no", "yes", 'u133A']),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("gene_order", SAME_AS_CONSTRAINT),
        Consequence("annotate", SAME_AS_CONSTRAINT),
        Consequence("rename_sample", SAME_AS_CONSTRAINT),
        Consequence("platform", SAME_AS_CONSTRAINT),
        #Consequence("logged", SET_TO, "yes"),
        #Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("format", SET_TO, "tdf"),
        Consequence("num_features", SET_TO, "no"),
        Consequence("unique_genes", SET_TO, "no"),
        Consequence("duplicate_probe", SET_TO, "no"),
        Consequence("group_fc", SET_TO, "no"),
        help="convert SignalFile_Annotate to SignalFile_Filter"
        ),
    ModuleNode(
        'remove_duplicate_genes',
        _SignalFile_Filter, _SignalFile_Filter,
        Constraint("annotate", MUST_BE, "yes"),
        Constraint("num_features", MUST_BE, "no"),
        Constraint("duplicate_probe", MUST_BE, 'no'),
        Constraint("unique_genes", MUST_BE, 'no'),
        Consequence("annotate", SAME_AS_CONSTRAINT),
        Consequence("num_features", SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT),
        Constraint("format", MUST_BE, "tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        #Constraint("logged", MUST_BE, "yes"),
        #Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence(
            "unique_genes", SET_TO_ONE_OF,
            ['average_genes', 'high_var', 'first_gene']),
        help="remove duplicate genes in SignalFile_Filter"),

    ModuleNode(
         'select_first_n_genes',
        _SignalFile_Filter, _SignalFile_Filter,
        OptionDef(
            "num_features_value", 500,
            help="num of features to be selected in the SignalFile"),
        Constraint("format", MUST_BE, "tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Constraint("num_features", MUST_BE, "no"),
        Constraint("duplicate_probe", MUST_BE, 'no'),
        Consequence("num_features", SET_TO, "yes"),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT),
        help="select first n genes in SignalFile_Filter"),

     ModuleNode(
        'remove_duplicate_probes',
        _SignalFile_Filter, _SignalFile_Filter,
        Constraint("format", MUST_BE, "tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Constraint("duplicate_probe", MUST_BE, 'no'),
        Consequence("duplicate_probe", SET_TO, 'high_var_probe'),
        Constraint("platform", CAN_BE_ANY_OF, ["yes", 'u133A']),
        Consequence("platform", SAME_AS_CONSTRAINT),
        help="remove duplciate probes in SignalFile_Filter by "
        "high_var_probe method"),
    ModuleNode(
         'select_probe_by_best_match',
        _SignalFile_Filter, _SignalFile_Filter,
        Constraint("format", MUST_BE, "tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Constraint("duplicate_probe", MUST_BE, 'no'),
        Consequence("duplicate_probe", SET_TO, 'closest_probe'),
        Constraint("platform", CAN_BE_ANY_OF, ["yes", 'u133A']),
        Consequence("platform", SAME_AS_CONSTRAINT),
        help="remove duplciate probes in SignalFile_Filter by closest_probe "
        "method"),

    ModuleNode(
        "filter_genes_by_fold_change_across_classes",
        [ClassLabelFile, _SignalFile_Filter], _SignalFile_Filter,
        OptionDef("group_fc_num", help="group fold change number"),
        Constraint("format", MUST_BE, "tdf", 1),
        Consequence("format", SAME_AS_CONSTRAINT, 1),
        Constraint("logged", MUST_BE, "yes", 1),
        Consequence("logged", SAME_AS_CONSTRAINT, 1),
        Constraint("group_fc", MUST_BE, "no", 1),
        Constraint("num_features", MUST_BE, "no", 1),
        Constraint("duplicate_probe", MUST_BE, "no", 1),
        Constraint("unique_genes", MUST_BE, "no", 1),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("num_features", SAME_AS_CONSTRAINT, 1),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT, 1),
        Consequence("unique_genes", SAME_AS_CONSTRAINT, 1),
        Consequence("group_fc", SET_TO, "yes"),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 1),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS, 1),
        #Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        DefaultAttributesFrom(1),
        help="filter genes in SignalFile_Filter by fold change in "
        "different classes"),
    ModuleNode(
        "convert_signal_to_gct",
        _SignalFile_Filter, _SignalFile_Filter,
        Constraint("format", MUST_BE, 'tdf'),
        Consequence("format", SET_TO, 'gct'),
        help="convert SignalFile_Filter in tdf format to gct format"),
    ModuleNode(
        'unlog_signal',
        _SignalFile_Filter, _SignalFile_Filter,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SET_TO, "no"),
        help="unlog SignalFile_Filter"),

    ModuleNode(
        'convert_signalfile_preprocess',
        SignalFile, SignalFile,
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SET_TO, "any"),
        help='convert preprocess from others to any'),
    
    ModuleNode(
        'convert_filter_final',
        _SignalFile_Filter, SignalFile,
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ['none', 'zero_fill', 'median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("combat_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("dwd_norm", CAN_BE_ANY_OF, YESNO),
        Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median", "mean"]),
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ["no", "variance", "sum_of_squares"]),
        Constraint("gene_order", CAN_BE_ANY_OF, BDT.GENE_ORDER),
        Constraint("annotate", CAN_BE_ANY_OF, YESNO),
        Constraint("rename_sample", CAN_BE_ANY_OF, YESNO),
        Constraint("platform", CAN_BE_ANY_OF, ["no", "yes", 'u133A']),
        Constraint("logged", CAN_BE_ANY_OF, YESNO),
        Constraint("format", CAN_BE_ANY_OF, ['tdf', 'gct']),
        Constraint("num_features", CAN_BE_ANY_OF, ['yes', "no"]),
        Constraint(
            "unique_genes", CAN_BE_ANY_OF,
            ["no", "average_genes", "high_var", "first_gene"]),
        Constraint(
            "duplicate_probe", CAN_BE_ANY_OF,
            ["no", "closest_probe", "high_var_probe"]),
        Constraint("group_fc", CAN_BE_ANY_OF, ['yes', "no"]),

        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("gene_order", SAME_AS_CONSTRAINT),
        Consequence("annotate", SAME_AS_CONSTRAINT),
        Consequence("rename_sample", SAME_AS_CONSTRAINT),
        Consequence("platform", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("num_features", SAME_AS_CONSTRAINT),
        Consequence("unique_genes", SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT),
        Consequence("group_fc", SAME_AS_CONSTRAINT),
        help="transfer SignalFile_Filter to SignalFile"
        ),
    ModuleNode(
        'plot_intensity_boxplot',
        SignalFile, IntensityPlot,
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        help="plot intensity boxplot"
        ),

    ModuleNode(
        "convert_simplelabelfile_to_classlabelfile",
        #SimpleLabelFile, ClassLabelFile,
        [UnprocessedSignalFile, SimpleLabelFile], ClassLabelFile,
        #Constraint("cls_format", MUST_BE, 'label', 0),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 0),
        #Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS1, 1),
        #Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),

        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        #Consequence("cls_format", SET_TO, 'cls'),
        help="Convert a SimpleLabelFile to ClassLabelFile.",
        ),
    ]

all_data_types = [
    SimpleLabelFile,
    ClassLabelFile,
    UnprocessedSignalFile,
    _SignalFile_Impute,
    _SignalFile_Merge,
    _SignalFile_Normalize,
    _SignalFile_Order,
    _SignalFile_Annotate,
    _SignalFile_Filter,
    SignalFile,
    IntensityPlot,
    ]
