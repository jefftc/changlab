#PcaAnalysis
from Betsy.bie3 import *
import SignalFile_rule,SignalFile2_rule,SignalFile1_rule
PcaAnalysis = DataType(
    'PcaAnalysis',
    Attribute('contents',["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "unspecified"],"unspecified","unspecified"),
    Attribute(
        "preprocess",["unknown", "illumina",
                                    "agilent", "mas5", "rma", "loess"],
        "unknown","unknown"),
    Attribute(
        "missing_values",["no"],'no','no'),
    Attribute(
        "missing_algorithm",["none", "median_fill", "zero_fill"],
        "none", "none"),
    Attribute(
        "logged",[ "no", "yes"],"yes","yes"),
    # Normalizing the genes.
    Attribute(
        "gene_center",["unknown", "no", "mean", "median"],
        "unknown","unknown"),
    Attribute(
        "gene_normalize",["unknown", "no", "variance", "sum_of_squares"],
        "unknown","unknown"),

    # Normalizing the data.  Very difficult to check normalization.
    # If you're not sure if the data is normalized, then the answer is
    # "no".
    Attribute(
        "dwd_norm",["no", "yes"], "no","no"),
    Attribute(
        "bfrm_norm",["no", "yes"], "no","no"),
    Attribute(
        "quantile_norm",["no", "yes"], "no","no"),
    Attribute(
        "shiftscale_norm",["no", "yes"], "no","no"),
    Attribute(
        "combat_norm",["no", "yes"], "no","no"),
    # Annotations.
    Attribute(
        "unique_genes",["no", "average_genes", "high_var", "first_gene"],
        "no","no"),
    Attribute(
        "duplicate_probe",["no", "yes", "closest_probe", "high_var_probe"],
        "no","no"),
    # Unclassified.
    Attribute("num_features",["yes","no"],"no","no"), 
    Attribute("predataset",["no", "yes"], "no","no"),
    Attribute("platform",["yes","no"],"no","no"),
    Attribute("filter",["yes","no"],"no","no"),
    Attribute("group_fc",["yes","no"],"no","no"))

PcaPlot = DataType(
    'PcaPlot',
    Attribute("pca_gene_num",["yes","no"],"no","no"),
    Attribute("contents",["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "unspecified"],"unspecified","unspecified"),
    # Properties of the data.
    Attribute("preprocess",["unknown", "illumina", "agilent",
                            "mas5", "rma", "loess"],
        "unknown","unknown"),
    Attribute(
        "missing_values",["no"],"no","no"),
    Attribute(
        "missing_algorithm",["none", "median_fill", "zero_fill"],
        "none","none"),
    Attribute(
        "logged",[ "no", "yes"], "yes","yes"),
    # Normalizing the genes.
    Attribute(
        "gene_center",["unknown", "no", "mean", "median"],
        "unknown","unknown"),
    Attribute(
        "gene_normalize",["unknown", "no", "variance", "sum_of_squares"],
        "unknown","unknown"),

    # Normalizing the data.  Very difficult to check normalization.
    # If you're not sure if the data is normalized, then the answer is
    # "no".
    Attribute(
        "dwd_norm",["no", "yes"], "no","no"),
    Attribute(
        "bfrm_norm",["no", "yes"], "no","no"),
    Attribute(
        "quantile_norm",["no", "yes"], "no","no"),
    Attribute(
        "shiftscale_norm",["no", "yes"], "no","no"),
    Attribute(
        "combat_norm",["no", "yes"], "no","no"),

    # Annotations.
    Attribute("annotate",["no", "yes"], "no","no"),
    Attribute(
        "unique_genes",["no", "average_genes", "high_var", "first_gene"],
        "no","no"),
    Attribute(
        "duplicate_probe",["no", "yes", "closest_probe", "high_var_probe"],
        "no","no"),
    Attribute("rename_sample",["no", "yes"], "no","no"),

    # Unclassified.
    Attribute("num_features",["yes","no"], "no","no"),
    Attribute("gene_order",[
            "no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr"],
        "no","no"),
    Attribute("predataset",["no", "yes"], "no","no"),
    Attribute("platform",["yes","no"],"no","no"),
    Attribute("filter",["yes","no"],"no","no"),
    Attribute("group_fc",["yes","no"],"no","no"),
    )
    

list_files = [PcaAnalysis,PcaPlot]
all_modules = [
    Module(
        'analyze_samples_pca',
        SignalFile2_rule.SignalFile2,PcaAnalysis,
        Constraint("contents",CAN_BE_ANY_OF,["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "unspecified"]),
        Constraint("format",MUST_BE,'tdf'),
        Constraint("logged",MUST_BE,'yes'),
        Constraint("filter",CAN_BE_ANY_OF,["yes","no"]),
        Constraint("preprocess",CAN_BE_ANY_OF,["unknown", "illumina",
                                               "agilent", "mas5", "rma", "loess"]),
        Constraint("quantile_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("bfrm_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("combat_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("shiftscale_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("dwd_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("gene_center",CAN_BE_ANY_OF,['mean','median','no','unknown']),
        Constraint("gene_normalize",CAN_BE_ANY_OF,["unknown", "no", "variance",
                                                   "sum_of_squares"]),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,["none", "median_fill", "zero_fill"]),
        Constraint("unique_genes",CAN_BE_ANY_OF,["no", "average_genes", "high_var",
                                                 "first_gene"]),
        Constraint("predataset",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("group_fc",CAN_BE_ANY_OF,["yes","no"]),
        Constraint("num_features",CAN_BE_ANY_OF,["yes","no"]),
        Constraint("platform",CAN_BE_ANY_OF,["yes","no"]),
        Constraint("annotate",MUST_BE,"no"),
        Consequence("contents",SAME_AS_CONSTRAINT),
        Consequence("preprocess",SAME_AS_CONSTRAINT),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("quantile_norm",SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT),
        Consequence("combat_norm",SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm",SAME_AS_CONSTRAINT),
        Consequence("dwd_norm",SAME_AS_CONSTRAINT),
        Consequence("gene_center",SAME_AS_CONSTRAINT),           
        Consequence("gene_normalize",SAME_AS_CONSTRAINT),          
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("unique_genes",SAME_AS_CONSTRAINT),
        Consequence("num_features",SAME_AS_CONSTRAINT),
        Consequence("filter",SAME_AS_CONSTRAINT),
        Consequence("predataset",SAME_AS_CONSTRAINT),
        Consequence("platform",SAME_AS_CONSTRAINT),
        Consequence("group_fc",SAME_AS_CONSTRAINT)),

    Module(
        'plot_sample_pca_wo_label',
        PcaAnalysis,PcaPlot,
        Constraint("contents",CAN_BE_ANY_OF,["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "unspecified"]),
        Constraint("logged",MUST_BE,'yes'),
        Constraint("filter",CAN_BE_ANY_OF,["yes","no"]),
        Constraint("preprocess",CAN_BE_ANY_OF,["unknown", "illumina",
                                               "agilent", "mas5", "rma", "loess"]),
        Constraint("quantile_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("bfrm_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("combat_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("shiftscale_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("dwd_norm",CAN_BE_ANY_OF,['yes','no']),
        Constraint("gene_center",CAN_BE_ANY_OF,['mean','median','no','unknown']),
        Constraint("gene_normalize",CAN_BE_ANY_OF,["unknown", "no", "variance",
                                                   "sum_of_squares"]),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,["none", "median_fill", "zero_fill"]),
        Constraint("unique_genes",CAN_BE_ANY_OF,["no", "average_genes", "high_var",
                                                 "first_gene"]),
        Constraint("predataset",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("group_fc",CAN_BE_ANY_OF,["yes","no"]),
        Constraint("num_features",CAN_BE_ANY_OF,["yes","no"]),
        Constraint("platform",CAN_BE_ANY_OF,["yes","no"]),
        Consequence("contents",SAME_AS_CONSTRAINT),
        Consequence("preprocess",SAME_AS_CONSTRAINT),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("quantile_norm",SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT),
        Consequence("combat_norm",SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm",SAME_AS_CONSTRAINT),
        Consequence("dwd_norm",SAME_AS_CONSTRAINT),
        Consequence("gene_center",SAME_AS_CONSTRAINT),           
        Consequence("gene_normalize",SAME_AS_CONSTRAINT),          
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("unique_genes",SAME_AS_CONSTRAINT),
        Consequence("num_features",SAME_AS_CONSTRAINT),
        Consequence("filter",SAME_AS_CONSTRAINT),
        Consequence("predataset",SAME_AS_CONSTRAINT),
        Consequence("platform",SAME_AS_CONSTRAINT),
        Consequence("group_fc",SAME_AS_CONSTRAINT)
        ),
             
 ]
