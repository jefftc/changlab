#PcaAnalysis
from Betsy.bie3 import *
import BasicDataTypes as BDT
import GeneExpProcessing

PcaAnalysis = DataType(
    'PcaAnalysis',
    AttributeDef(
        'contents', BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "preprocess", BDT.PREPROCESS,
        "unknown", "unknown", help="preprocess method"),
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef(
        "logged", ["no", "yes"], "yes", "yes", help="logged or not"),
    # Normalizing the genes.
    AttributeDef(
        "gene_center", ["no", "mean", "median"],
        "no", "no", help="gene center method"),
    AttributeDef(
        "gene_normalize", ["unknown", "no", "variance", "sum_of_squares"],
        "no", "no", help="gene normalize method"),

    # Normalizing the data.  Very difficult to check normalization.
    # If you're not sure if the data is normalized, then the answer is
    # "no".
    AttributeDef(
        "dwd_norm", ["no", "yes"], "no", "no", help="dwd normalization"),
    AttributeDef(
        "bfrm_norm", ["no", "yes"], "no", "no", help="bfrm normalization"),
    AttributeDef(
        "quantile_norm", ["no", "yes"], "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "shiftscale_norm", ["no", "yes"], "no", "no",
        help="shiftscale normalization"),
    AttributeDef(
        "combat_norm", ["no", "yes"], "no", "no", help="combat normalization"),
    # Annotations.
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no", help="method to get unique genes"),
    AttributeDef(
        "duplicate_probe", ["no", "yes", "closest_probe", "high_var_probe"],
        "no", "no", help="method to remove duplicated probes"),
    # Unclassified.
    AttributeDef(
        "num_features", ["yes", "no"], "no", "no",
        help="select a num of features or not"), 
    AttributeDef(
        "filter_and_threshold", ["no", "yes"], "no", "no",
        help="filter_and_threshold or not"),
    AttributeDef(
        "platform", ["yes", "no"], "no", "no", help="add platform or not"),
    AttributeDef(
        "filter_missing_values", ["yes", "no"], "no", "no",
        help="filter missing or not"),
    AttributeDef(
        "group_fc", ["yes", "no"], "no", "no",
        help="group fold change or not"),
    help="Pca Analysis File")

PcaPlot = DataType(
    'PcaPlot',
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    # Properties of the data.
    AttributeDef("preprocess", BDT.PREPROCESS,
        "unknown", "unknown", help="preprocess method"),
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef(
        "logged", ["no", "yes"], "yes", "yes", help="logged or not"),
    # Normalizing the genes.
    AttributeDef(
        "gene_center", ["no", "mean", "median"],
        "no", "no", help="gene center method"),
    AttributeDef(
        "gene_normalize", ["no", "variance", "sum_of_squares"],
        "no", "no", help="gene normalize method"),

    AttributeDef(
        "dwd_norm", ["no", "yes"], "no", "no", help="dwd normalization"),
    AttributeDef(
        "bfrm_norm", ["no", "yes"], "no", "no", help="bfrm normalization"),
    AttributeDef(
        "quantile_norm", ["no", "yes"], "no", "no",
        help="quantile normalization"),
    AttributeDef(
        "shiftscale_norm", ["no", "yes"], "no", "no",
        help="shiftscale normalization"),
    AttributeDef(
        "combat_norm", ["no", "yes"], "no", "no", help="combat normalization"),

    # Annotations.
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no", help="method to get unique genes"),
    AttributeDef(
        "duplicate_probe", ["no", "yes", "closest_probe", "high_var_probe"],
        "no", "no", help="method to remove duplicated probes"),
    
    # Unclassified.
    AttributeDef(
        "num_features", ["yes", "no"], "no", "no",
        help="select a num of features or not"),
    AttributeDef(
        "filter_and_threshold", ["no", "yes"], "no", "no",
        help="filter_and_threshold or not"),
    AttributeDef(
        "platform", ["yes", "no"], "no", "no", help="add platform or not"),
    AttributeDef(
        "filter_missing_values", ["yes", "no"], "no", "no",
        help="filter missing or not"),
    AttributeDef(
        "group_fc", ["yes", "no"], "no", "no",
        help="group fold change or not"),
    help="Pca analysis plot file"
    )
    

all_data_types = [PcaAnalysis,PcaPlot]
all_modules = [
    ModuleNode(
        'analyze_samples_pca',
        GeneExpProcessing.SignalFile, PcaAnalysis,
        OptionDef('pca_gene_num',500, help="number of genes in pca"),
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Constraint("format",MUST_BE,'tdf'),
        Constraint("logged",MUST_BE,'yes'),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, ["yes", "no"]),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("combat_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("gene_center", CAN_BE_ANY_OF, ['mean','median','no']),
        Constraint("gene_normalize", CAN_BE_ANY_OF, ["no", "variance",
                                                  "sum_of_squares"]),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ["none", "median_fill", "zero_fill"]),
        Constraint(
            "unique_genes", CAN_BE_ANY_OF,
            ["no", "average_genes", "high_var", "first_gene"]),
        Constraint("filter_and_threshold", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("group_fc", CAN_BE_ANY_OF, ["yes", "no"]),
        Constraint("num_features", CAN_BE_ANY_OF, ["yes", "no"]),
        Constraint("platform", CAN_BE_ANY_OF, ["yes", "no"]),
        #Constraint("annotate",MUST_BE, "no"),
        #Constraint("annotate",MUST_BE, "yes"),  ## for the big network
        #Constraint("gene_order",MUST_BE,'class_neighbors'),## for the big network
        #Constraint("rename_sample",MUST_BE,'yes'),## for the big network
        #Consequence("annotate", SAME_AS_CONSTRAINT),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),           
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),          
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        Consequence("unique_genes", SAME_AS_CONSTRAINT),
        Consequence("num_features", SAME_AS_CONSTRAINT),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
        Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
        Consequence("platform", SAME_AS_CONSTRAINT),
        Consequence("group_fc", SAME_AS_CONSTRAINT),
        help="analyze samples by pca method"
        ),
    
    ModuleNode(
        'plot_sample_pca_wo_label',
        PcaAnalysis,PcaPlot,
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Constraint("logged",MUST_BE,'yes'),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, ["yes", "no"]),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("combat_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("gene_center", CAN_BE_ANY_OF, ['mean','median','no']),
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ["no", "variance", "sum_of_squares"]),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ["none", "median_fill", "zero_fill"]),
        Constraint(
            "unique_genes", CAN_BE_ANY_OF,
            ["no", "average_genes", "high_var", "first_gene"]),
        Constraint("filter_and_threshold", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("group_fc", CAN_BE_ANY_OF, ["yes", "no"]),
        Constraint("num_features", CAN_BE_ANY_OF, ["yes", "no"]),
        Constraint("platform", CAN_BE_ANY_OF, ["yes", "no"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),           
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),          
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        Consequence("unique_genes", SAME_AS_CONSTRAINT),
        Consequence("num_features", SAME_AS_CONSTRAINT),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT),
        Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
        Consequence("platform", SAME_AS_CONSTRAINT),
        Consequence("group_fc", SAME_AS_CONSTRAINT),
        help="plot pca analysis result"
        ),
    ]
