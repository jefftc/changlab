#SignalFile1
from Betsy import bie3
from Betsy.bie3 import *
import SignalFile_rule
SignalFile1 = DataType(
    "SignalFile1",
    # Properties of the format.
    Attribute("format",[ "tdf", "pcl"],"tdf","tdf"),
    # Properties of the data.
    Attribute(
        "preprocess",["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        "unknown","unknown"),
    Attribute(
        "missing_values",["no"],"no","no"),
    Attribute(
        "missing_algorithm",["none", "median_fill", "zero_fill"],"none","none"),
    Attribute(
        "logged",["yes"],"yes","yes"),
    Attribute("predataset",["no", "yes"], "no","no"),
    Attribute("filter",["yes","no"], "no","no"),
    Attribute("rename_sample",["no", "yes"], "no","no"),
    Attribute("dwd_norm",["no", "yes"], "no","no"),
    Attribute("bfrm_norm",["no", "yes"], "no","no"),
    Attribute("quantile_norm",["no", "yes"],"no", "no"),
    Attribute("shiftscale_norm",["no", "yes"], "no","no"),
    Attribute("combat_norm",["no", "yes"], "no","no"),
    Attribute(
        "gene_center",["unknown", "no", "mean", "median"],
        "unknown","no"),
    Attribute(
        "gene_normalize",["unknown", "no", "variance", "sum_of_squares"],
        "unknown","no"),
    Attribute("contents",["train0", "train1", "test", 'class0,class1,test',
                        "class0", "class1", "class0,class1", "unspecified"],
                  "unspecified","unspecified")
    )
list_files = [SignalFile1]
all_modules = [
    Module(
        "transfer1",
        SignalFile_rule.SignalFile,SignalFile1,
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("preprocess",CAN_BE_ANY_OF,["unknown", "illumina", "agilent",
                                "mas5", "rma", "loess"]),
        Constraint("rename_sample",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,["none", "median_fill", "zero_fill"]),
        Constraint("predataset",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("filter",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("dwd_norm",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("bfrm_norm",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("quantile_norm",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("shiftscale_norm",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("combat_norm",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("contents",CAN_BE_ANY_OF,["train0", "train1", "test", 'class0,class1,test',
                        "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("format",SAME_AS_CONSTRAINT,0),
        Consequence("logged",SAME_AS_CONSTRAINT,0),
        Consequence("missing_values",SAME_AS_CONSTRAINT,0),
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT,0),           
        Consequence("preprocess",SAME_AS_CONSTRAINT,0),
        Consequence("predataset",SAME_AS_CONSTRAINT,0),
        Consequence("rename_sample",SAME_AS_CONSTRAINT,0),
        Consequence("filter",SAME_AS_CONSTRAINT,0),
        Consequence("dwd_norm",SAME_AS_CONSTRAINT,0),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT,0),
        Consequence("quantile_norm",SAME_AS_CONSTRAINT,0),
        Consequence("shiftscale_norm",SAME_AS_CONSTRAINT,0),
        Consequence("combat_norm",SAME_AS_CONSTRAINT,0),
        Consequence("gene_center",SET_TO,'unknown'),
        Consequence("gene_normalize",SET_TO,'unknown'),
        Consequence("contents",SAME_AS_CONSTRAINT,0)),
    Module(
        "check_gene_center",
        SignalFile1,SignalFile1,
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("gene_center",MUST_BE,"unknown"),
        Constraint("gene_normalize",MUST_BE,"unknown"),
        Consequence("format",SAME_AS_CONSTRAINT,0),
        Consequence("logged",SAME_AS_CONSTRAINT,0),
        Consequence("missing_values",SAME_AS_CONSTRAINT,0),
        Consequence("gene_center",BASED_ON_DATA,["no", "mean", "median"])),
    Module(
        "check_gene_normalize",
        SignalFile1,SignalFile1,
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("gene_normalize",MUST_BE,"unknown"),
        Consequence("format",SAME_AS_CONSTRAINT,0),
        Consequence("logged",SAME_AS_CONSTRAINT,0),
        Consequence("missing_values",SAME_AS_CONSTRAINT,0),
        Consequence("gene_normalize",BASED_ON_DATA,["no", "variance", "sum_of_squares"])),
    Module(   
        "convert_signal_to_pcl",
        SignalFile1,SignalFile1,
        Constraint("format",MUST_BE,'tdf'),
        Consequence("format",SET_TO,'pcl')),
    Module(
        "center_genes",
        SignalFile1,SignalFile1,
        Constraint("format",MUST_BE,"pcl"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("gene_center",MUST_BE,"no"),
        Constraint("gene_normalize",MUST_BE,'unknown'),
        Consequence("format",SAME_AS_CONSTRAINT,0),
        Consequence("logged",SAME_AS_CONSTRAINT,0),
        Consequence("missing_values",SAME_AS_CONSTRAINT,0),
        Consequence("gene_center",SET_TO_ONE_OF,["mean", "median"]),
        Consequence("gene_normalize",SAME_AS_CONSTRAINT,0)),
    Module(
        "normalize_genes",
        SignalFile1,SignalFile1,
        Constraint("format",MUST_BE,"pcl"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("gene_normalize",MUST_BE,"no"),
        Consequence("format",SAME_AS_CONSTRAINT,0),
        Consequence("logged",SAME_AS_CONSTRAINT,0),
        Consequence("missing_values",SAME_AS_CONSTRAINT,0),
        Consequence("gene_normalize",SET_TO_ONE_OF,["variance", "sum_of_squares"]))
    ]
