# DataTypes:
# SimpleLabelFile
# ClassLabelFile
#
# UnprocessedSignalFile
# _SignalFile_Impute
# _SignalFile_Merge
# _SignalFile_Normalize
# _SignalFile_Order
# _SignalFile_Annotate
# _SignalFile_Filter
# SignalFile
#
# IlluminaControlFile
#
#
# Modules:
#   UnprocessedSignalFile
# convert_signal_to_tdf
# check_for_log
# log_signal
# convert_postprocess_impute
#
#   _SignalFile_Impute
# check_for_missing_values
# filter_genes_by_missing_values
# fill_missing_with_zeros
# fill_missing_with_median
# convert_impute_merge
#
#   _SignalFile_Merge
# merge_class0_class1_signal_signal
# merge_class0_class1_signal_classlabel
# normalize_samples_with_quantile
# normalize_samples_with_bfrm
# normalize_samples_with_combat
# normalize_samples_with_dwd
# normalize_samples_with_shiftscale
# convert_merge_normalize
#
#   _SignalFile_Normalize
# check_gene_center
# check_gene_normalize
# convert_signal_to_pcl
# center_genes
# normalize_genes
# convert_normalize_order
#
#   _SignalFile_Order
# rank_genes_by_class_neighbors
# rank_genes_by_sample_ttest
# reorder_genes
# convert_order_annotate
#
#   _SignalFile_Annotate
# annotate_probes
# relabel_samples
# add_crossplatform_probeid
# add_U133A_probeid           # why is this separate?
# convert_annotate_filter
#
#   _SignalFile_Filter
# remove_duplicate_genes
# select_first_n_genes
# remove_duplicate_probes
# select_probe_by_best_match
# filter_genes_by_fold_change_across_classes
# convert_signal_to_gct
# unlog_signal
# convert_filter_final
#
#   SignalFile
# convert_signalfile_preprocess
#
#
#   Other Stuff
# convert_simplelabelfile_to_classlabelfile




# Merging:
# convert_imput_merge_rma        TODO: Don't separate out RMA
# merge_two_classes_rma
# merge_two_classes              Merge class0 and class1 to class0,class1.

# Need:
# merge_class0_class1_signal_signal
# merge_class0_class1_signal_classlabel
# merge_class0_class1_classlabel_classlabel

from Betsy.bie3 import *
import BasicDataTypes as BDT

YESNO = BDT.YESNO


# Make some variables for common attributes, for convenience.
ATTR_CONTENTS = AttributeDef(
    "contents", BDT.CONTENTS, "unspecified", "unspecified",
    help="contents")
ATTR_PREPROCESS = AttributeDef(
    "preprocess", BDT.PREPROCESS, "unknown", "unknown",
    help="preprocess method")
EXPRESSION_FORMATS = ["tdf", "pcl", "gct", "res", "jeffs"]
ATTR_FORMAT = AttributeDef(
    "format", ["unknown"]+EXPRESSION_FORMATS,
    "unknown", "tdf", help="file format")
ATTR_LOGGED = AttributeDef(
    "logged", ["unknown", "no", "yes"], "unknown", "no",
    help="logged or not")


SimpleClassFile = DataType(
    "SimpleClassFile",
    # Should have headers:
    # Sample   Class
    ATTR_CONTENTS,
    ATTR_PREPROCESS,
    help="A simple tab-delimited format with labels for each sample.  "
    "Can have two or more different classes.",
    )

ClassLabelFile = DataType(
    "ClassLabelFile",
    ATTR_CONTENTS,
    ATTR_PREPROCESS,
    help="CLS format file with categorical labels.  Usually easier to use "
    "SimpleClassFile."
    )

UNPROC_ATTRIBUTES = [
    ATTR_CONTENTS,
    ATTR_PREPROCESS,
    ATTR_LOGGED,
    ]
IMPUTE_ATTRIBUTES = [
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill","zero_fill", help="missing algorithm"),
    AttributeDef(
        "filter_missing_values", YESNO, "no", "no",
        help="Whether missing values are filtered."),
    ]
MERGE_ATTRIBUTES = [
    AttributeDef(
        "quantile_norm", YESNO, "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "combat_norm", YESNO, "no", "no", help="combat normalization"),
    AttributeDef(
        "dwd_norm", YESNO, "no", "no", help="dwd normalization"),
    AttributeDef(
        "shiftscale_norm", YESNO, "no", "no", help="shiftscale normalization"),
    AttributeDef(
        "bfrm_norm", YESNO, "no", "no", help="bfrm normalization"),
    ]
NORMALIZE_ATTRIBUTES = [
    AttributeDef(
        "gene_center", ["unknown", "no", "mean", "median"],
        "unknown", "no", help="gene center method"),
    AttributeDef(
        "gene_normalize", ["unknown", "no", "variance", "sum_of_squares"],
        "unknown", "no", help="gene normalize method"),
    ]
ORDER_ATTRIBUTES = [
    AttributeDef(
        "gene_order", BDT.GENE_ORDER, "none", "none",
        help="gene order method"),
    ]
ANNOTATE_ATTRIBUTES = [
    AttributeDef(
        "relabel_sample", YESNO, "no", "no",
        help="rename sample or not"),
    AttributeDef(
        "annotate", YESNO, "no", "no",
        help="annotate file or not [WHAT DOES THIS MEAN?]"),
    AttributeDef(
        # Why is the u133A?
        "platform", ["yes", "no", 'u133A'], "no", "no",
        help="add platform or not"),
    ]
FILTER_ATTRIBUTES = [
    # TODO: Rename this.  Determines whether to select_first_n_genes.
    # Maybe merge with filter_and_threshold.
    AttributeDef(
        "num_features", YESNO, "no", "no",
        help="select a num of features or not"),
    # Whether to filter and threshold genes.
    AttributeDef(
        "filter_and_threshold", YESNO, "no", "no",
        help="Whether the values are filtered and thresholded."),
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no", help="method to get unique genes"),
    AttributeDef(
        "duplicate_probe", ["no", "closest_probe", "high_var_probe"],
        "no", "no", help="method to remove duplicated probes"),
    # What is this for?
    AttributeDef(
        "group_fc", YESNO, "no", "no", help="group fold change or not"),
    ]

CONVERT_FORMAT = [
    Constraint("format", CAN_BE_ANY_OF, EXPRESSION_FORMATS),
    Consequence("format", SAME_AS_CONSTRAINT),
    ]
CONVERT_UNPROC = [
    Constraint("logged", CAN_BE_ANY_OF, YESNO),
    Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
    Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    Consequence("logged", SAME_AS_CONSTRAINT),
    Consequence("preprocess", SAME_AS_CONSTRAINT),
    Consequence("contents", SAME_AS_CONSTRAINT),
    ]
CONVERT_IMPUTE = [
    Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO),
    Constraint(
        "missing_algorithm", CAN_BE_ANY_OF,
        ['none', 'zero_fill', 'median_fill']),
    Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
    Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
    ]
CONVERT_MERGE = [
    Constraint("quantile_norm", CAN_BE_ANY_OF, YESNO),
    Constraint("combat_norm", CAN_BE_ANY_OF, YESNO),
    Constraint("dwd_norm", CAN_BE_ANY_OF, YESNO),
    Constraint("shiftscale_norm", CAN_BE_ANY_OF, YESNO),
    Constraint("bfrm_norm", CAN_BE_ANY_OF, YESNO),
    Consequence("quantile_norm", SAME_AS_CONSTRAINT),
    Consequence("combat_norm", SAME_AS_CONSTRAINT),
    Consequence("dwd_norm", SAME_AS_CONSTRAINT),
    Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
    Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
    ]
CONVERT_NORMALIZE = [
    Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median", "mean"]),
    Constraint(
        "gene_normalize", CAN_BE_ANY_OF,
        ["no", "variance", "sum_of_squares"]),
    Consequence("gene_center", SAME_AS_CONSTRAINT),
    Consequence("gene_normalize", SAME_AS_CONSTRAINT),
    ]
CONVERT_ORDER = [
    Constraint("gene_order", CAN_BE_ANY_OF, BDT.GENE_ORDER),
    Consequence("gene_order", SAME_AS_CONSTRAINT),
    ]
CONVERT_ANNOTATE = [
    Constraint("annotate", CAN_BE_ANY_OF, YESNO),
    Constraint("relabel_sample", CAN_BE_ANY_OF, YESNO),
    Constraint("platform", CAN_BE_ANY_OF, ["no", "yes", 'u133A']),
    Consequence("annotate", SAME_AS_CONSTRAINT),
    Consequence("relabel_sample", SAME_AS_CONSTRAINT),
    Consequence("platform", SAME_AS_CONSTRAINT),
    ]
CONVERT_FILTER = [
    Constraint("num_features", CAN_BE_ANY_OF, YESNO),
    Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
    Constraint(
        "unique_genes", CAN_BE_ANY_OF,
        ["no", "average_genes", "high_var", "first_gene"]),
    Constraint(
        "duplicate_probe", CAN_BE_ANY_OF,
        ["no", "closest_probe", "high_var_probe"]),
    Constraint("group_fc", CAN_BE_ANY_OF, YESNO),
    Consequence("num_features", SAME_AS_CONSTRAINT),
    Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
    Consequence("unique_genes", SAME_AS_CONSTRAINT),
    Consequence("duplicate_probe", SAME_AS_CONSTRAINT),
    Consequence("group_fc", SAME_AS_CONSTRAINT),
    ]




UnprocessedSignalFile = DataType(
    "UnprocessedSignalFile",
    ATTR_FORMAT,
    *UNPROC_ATTRIBUTES,
    help="Gene expression matrix before any processing.")

_SignalFile_Impute = DataType(
    "_SignalFile_Impute",
    ATTR_FORMAT,
    AttributeDef(
        "missing_values", ["unknown", "no", "yes"], "unknown", "no",
        help="missing values unknown, yes or no"),
    *(UNPROC_ATTRIBUTES+IMPUTE_ATTRIBUTES),
    help="The SignalFile after SignalFile_Postprocess, care missing_values, "
    "missing_algorithm and filter.")

_SignalFile_Merge = DataType(
    "_SignalFile_Merge",
    ATTR_FORMAT,
    *(UNPROC_ATTRIBUTES+IMPUTE_ATTRIBUTES+MERGE_ATTRIBUTES),
    help="Merge signal files and handle of batch effects."
    )

_SignalFile_Normalize = DataType(
    "_SignalFile_Normalize",
    # UnprocessedSignalFile
    #AttributeDef(
    #    "format", ["tdf", "pcl"], "tdf", "tdf", help="file format"),
    ATTR_FORMAT,
    *(UNPROC_ATTRIBUTES+IMPUTE_ATTRIBUTES+MERGE_ATTRIBUTES+
      NORMALIZE_ATTRIBUTES),
    help="Gene level centering and normalization."
    )

_SignalFile_Order = DataType(
    "_SignalFile_Order",
    ATTR_FORMAT,
    *(UNPROC_ATTRIBUTES+IMPUTE_ATTRIBUTES+MERGE_ATTRIBUTES+
      NORMALIZE_ATTRIBUTES+ORDER_ATTRIBUTES),
    help="Changes the order of genes in the expression file."
    )

_SignalFile_Annotate = DataType(
    "_SignalFile_Annotate",
    ATTR_FORMAT,
    *(UNPROC_ATTRIBUTES+IMPUTE_ATTRIBUTES+MERGE_ATTRIBUTES+
      NORMALIZE_ATTRIBUTES+ORDER_ATTRIBUTES+ANNOTATE_ATTRIBUTES),
    help="Annotates the genes with identifiers."
    )

_SignalFile_Filter = DataType(
    "_SignalFile_Filter",
    # Not sure why these are back now.
    #AttributeDef("format", [ "tdf", "gct"], "tdf", "tdf", help="file format"),
    #AttributeDef("logged", [ "no", "yes"], "yes", "yes", help="logged or not"),
    ATTR_FORMAT,
    *(UNPROC_ATTRIBUTES+IMPUTE_ATTRIBUTES+MERGE_ATTRIBUTES+
      NORMALIZE_ATTRIBUTES+ORDER_ATTRIBUTES+ANNOTATE_ATTRIBUTES+
      FILTER_ATTRIBUTES),
    help="Filtering genes."
    )

SignalFile= DataType(
    "SignalFile",
    #AttributeDef("format", ["tdf", "gct"], "tdf", "tdf", help="file format"),
    #AttributeDef("logged", YESNO, "yes", "yes", help="logged or not"),
    # Handle UNPROC_ATTRIBUTES myself, since "preprocess" can also
    # take value of "any".
    ATTR_FORMAT,
    ATTR_CONTENTS,
    AttributeDef(
        "preprocess", BDT.ANY_PREPROCESS, "unknown", "any",
        help="preprocess method"),
    ATTR_LOGGED,
    *(IMPUTE_ATTRIBUTES+MERGE_ATTRIBUTES+
      NORMALIZE_ATTRIBUTES+ORDER_ATTRIBUTES+ANNOTATE_ATTRIBUTES+
      FILTER_ATTRIBUTES),
    help="A processed gene expression matrix."
    )


IlluminaControlFile = DataType(
    "IlluminaControlFile",
    # preprocess=illumina
    # logged=no
    # format=gct
    ATTR_CONTENTS,
    AttributeDef(
        "relabel_sample", YESNO, "no", "no",
        help="rename sample or not"),
    #AttributeDef(
    #    'missing_values', ["unknown", "no", "yes"], "no", "no",
    #    help="missing values yes or not"),
    #AttributeDef(
    #    "missing_algorithm", ["none", "median_fill", "zero_fill"],
    #    "zero_fill", "zero_fill", help="missing algorithm"),
    help="Contains data for Illumina control probes.  "
    "Separate from gene expression data.")


all_data_types = [
    SimpleClassFile,
    ClassLabelFile,
    
    UnprocessedSignalFile,
    _SignalFile_Impute,
    _SignalFile_Merge,
    _SignalFile_Normalize,
    _SignalFile_Order,
    _SignalFile_Annotate,
    _SignalFile_Filter,
    SignalFile,
    IlluminaControlFile,
    ]

all_modules = [
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

        help="Determine whether the gene expression data in this "
        "UnprocessedSignalFile has been logged.  This is done using a "
        "heuristic.  We assume the gene expression is logged if "
        "all values are < 1000.  "
        "If at least 10 or 1% of the values are >= 1000, then we assume the "
        "gene expression has not been logged.  "
        "Otherwise, it is difficult to determine whether the values are "
        "not logged, or if there are just some outliers; and we raise an "
        "error."
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
        Constraint("format", MUST_BE, "tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes"),
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

    # Impute
    ModuleNode(
        "convert_postprocess_impute",
        UnprocessedSignalFile, _SignalFile_Impute,
        Constraint("format", MUST_BE, "tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        *CONVERT_UNPROC,
        # Why are these set?
        #Consequence("missing_algorithm", SET_TO, "zero_fill"),
        #Consequence("missing_values", SET_TO, "unknown"),
        #Consequence("filter_missing_values", SET_TO, "no"),
        help="convert SignalFile_Postprocess to SignalFile_Impute"
        ),

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
        Constraint("missing_values", MUST_BE, "no"),
        *(
            CONVERT_FORMAT +
            CONVERT_UNPROC +
            CONVERT_IMPUTE
            ),
        help="convert SignalFile_Impute to SignalFile_Merge"
        ),

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
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
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
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        DefaultAttributesFrom(1),
        help="nommalize SignalFile_Merge with shiftscale method"),
    
    ###normalize
    ModuleNode(
        "convert_merge_normalize",
        _SignalFile_Merge, _SignalFile_Normalize,
        *(
            CONVERT_FORMAT +
            CONVERT_UNPROC +
            CONVERT_IMPUTE +
            CONVERT_MERGE
            ),
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
        *(
            CONVERT_FORMAT +
            CONVERT_UNPROC +
            CONVERT_IMPUTE +
            CONVERT_MERGE +
            CONVERT_NORMALIZE
            ),
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
        Consequence("gene_order", SET_TO, "class_neighbors"),

        Constraint("contents", MUST_BE, "class0,class1", 0),
        Constraint("contents", SAME_AS, 0, 1),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        
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
        
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT),

        Constraint("gene_center", CAN_BE_ANY_OF, ['median', 'mean', 'no'], 0),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ['sum_of_squares', 'variance', 'no'], 0),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
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
        *(
            CONVERT_FORMAT +
            CONVERT_UNPROC +
            CONVERT_IMPUTE +
            CONVERT_MERGE +
            CONVERT_NORMALIZE +
            CONVERT_ORDER
            ),
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
        [_SignalFile_Annotate, BDT.RenameFile], _SignalFile_Annotate,
        OptionDef(
            "sample_labels_header", 
            help="Header that contains the desired sample labels."),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("relabel_sample", MUST_BE, "no", 0),
        Consequence("relabel_sample", SET_TO, "yes"),
        Constraint('annotate', MUST_BE, "no", 0),
        Consequence("annotate", SAME_AS_CONSTRAINT),
        #Constraint("labels_from", CAN_BE_ANY_OF, ["title", "description"], 0),
        Constraint("platform", MUST_BE, "no", 0),
        Consequence("platform", SAME_AS_CONSTRAINT),
        help="relabel the sample names in SignalFile_Annotate given "
        "RenameFile.",
       ),
    ModuleNode(
        "relabel_samples_illu_control",
        [IlluminaControlFile, BDT.RenameFile], IlluminaControlFile,
        OptionDef(
            "sample_labels_header", 
            help="Header that contains the desired sample labels."),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("relabel_sample", MUST_BE, "no", 0),
        Consequence("relabel_sample", SET_TO, "yes"),
        ),
    ModuleNode(
        'add_crossplatform_probeid',
        _SignalFile_Annotate, _SignalFile_Annotate,
        OptionDef(
            "platform_name",
            help="given the new platform name to add to the file"),
        Constraint("platform", MUST_BE, "no"),
        Consequence("platform", SET_TO, "yes"),
        help="add a cross platform to SignalFile_Annotate",
        ),

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
        *(
            CONVERT_FORMAT +
            CONVERT_UNPROC +
            CONVERT_IMPUTE +
            CONVERT_MERGE +
            CONVERT_NORMALIZE +
            CONVERT_ORDER +
            CONVERT_ANNOTATE
            ),
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
        Consequence(
            "unique_genes", SET_TO_ONE_OF,
            ['average_genes', 'high_var', 'first_gene']),
        help="remove duplicate genes in SignalFile_Filter"),

    ModuleNode(
         'select_first_n_genes',
        _SignalFile_Filter, _SignalFile_Filter,
        OptionDef(
            "num_first_genes", 500,
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
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes"),
        Consequence("logged", SET_TO, "no"),
        help="unlog SignalFile_Filter"),

    ModuleNode(
        'convert_signalfile_preprocess',
        SignalFile, SignalFile,
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SET_TO, "any"),
        help="Convert preprocess from others to any"),
    
    ModuleNode(
        'convert_filter_final',
        _SignalFile_Filter, SignalFile,
        *(
            CONVERT_FORMAT +
            CONVERT_UNPROC +
            CONVERT_IMPUTE +
            CONVERT_MERGE +
            CONVERT_NORMALIZE +
            CONVERT_ORDER +
            CONVERT_ANNOTATE +
            CONVERT_FILTER
            ),
        help="transfer SignalFile_Filter to SignalFile"
        ),

    ModuleNode(
        "convert_simplelabelfile_to_classlabelfile",
        [UnprocessedSignalFile, SimpleClassFile], ClassLabelFile,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        Constraint("format", CAN_BE_ANY_OF, ["unknown"]+EXPRESSION_FORMATS, 0),
        Constraint("logged", CAN_BE_ANY_OF, YESNO, 0),
        help="Convert a SimpleClassFile to ClassLabelFile.",
        ),
    ]
