# SignalFile
#from Betsy import bie3
from Betsy.bie3 import *
import SignalFile_rule

IntensityPlot = DataType(
    'IntensityPlot',
    AttributeDef(
        "contents",
        ["train0", "train1", "test", "class0,class1,test",
         "class0", "class1", "class0,class1", "unspecified"],
        "unspecified", "unspecified")
    )
    
ActbPlot = DataType(
    'ActbPlot',
    AttributeDef(
        "contents",
        ["train0", "train1", "test", "class0,class1,test",
         "class0", "class1", "class0,class1", "unspecified"],
        "unspecified", "unspecified")
    )
Hyb_barPlot = DataType(
    'Hyb_barPlot',
    AttributeDef(
        "contents",
        ["train0", "train1", "test", "class0,class1,test",
         "class0", "class1", "class0,class1", "unspecified"],
        "unspecified", "unspecified"))
ControlPlot = DataType(
    'ControlPlot',
    AttributeDef(
        "contents",
        ["train0", "train1", "test", "class0,class1,test",
         "class0", "class1", "class0,class1", "unspecified"],
        "unspecified", "unspecified"))
BiotinPlot = DataType(
    'BiotinPlot',
    AttributeDef(
        "contents",
        ["train0", "train1", "test", "class0,class1,test",
         "class0", "class1", "class0,class1", "unspecified"],
        "unspecified", "unspecified"))
HousekeepingPlot = DataType(
    'HousekeepingPlot',
    AttributeDef(
        "contents",
        ["train0", "train1", "test", "class0,class1,test",
         "class0", "class1", "class0,class1", "unspecified"],
        "unspecified", "unspecified"))

list_files=[
    IntensityPlot,
    ActbPlot,
    Hyb_barPlot,
    BiotinPlot,
    HousekeepingPlot,
    ControlPlot,
    ]

all_modules = [
    Module(    
        'plot_intensity_boxplot',
        SignalFile_rule.PrettySignalFile, IntensityPlot,
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        ),
                    
    Module(
        'plot_actb_line',
        SignalFile_rule.PrettySignalFile, ActbPlot,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("preprocess",CAN_BE_ANY_OF,['mas5','agilent','loess','unknown','illumina']),
        Constraint("missing_values", MUST_BE,"no"),
        Constraint("quantile_norm",MUST_BE,"no"),
        Constraint("combat_norm", MUST_BE,"no"),
        Constraint("shiftscale_norm", MUST_BE,"no"),
        Constraint("dwd_norm", MUST_BE,"no"),
        Constraint("bfrm_norm", MUST_BE,"no"),
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        ),
    Module(
        'plot_actb_line',
        SignalFile_rule.PrettySignalFile, ActbPlot,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("preprocess",MUST_BE,"rma"),
        Constraint("missing_values", MUST_BE,"no"),
        Constraint("quantile_norm", MUST_BE,"yes",),
        Constraint("combat_norm", MUST_BE,"no"),
        Constraint("shiftscale_norm", MUST_BE,"no"),
        Constraint("dwd_norm", MUST_BE,"no"),
        Constraint("bfrm_norm", MUST_BE,"no"),
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        ),
    Module(
        'plot_affy_affx_line',
        SignalFile_rule.PrettySignalFile,ControlPlot,
        Constraint("annotate", MUST_BE,'no'),
        Constraint("contents",CAN_BE_ANY_OF,["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"]),
        Consequence("contents",SAME_AS_CONSTRAINT)
        ),
    Module(
       'plot_illu_hyb_bar',
       SignalFile_rule.ControlFile,Hyb_barPlot,
       Constraint("preprocess", MUST_BE,"illumina"),
       Constraint("format", MUST_BE,"gct"),
       Constraint("logged", MUST_BE,"no"),
       Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("contents", SAME_AS_CONSTRAINT),),
    Module(
        'plot_illu_biotin_line',
        SignalFile_rule.ControlFile,BiotinPlot,
        Constraint("preprocess", MUST_BE,'illumina'),
        Constraint("format", MUST_BE,"gct"),
        Constraint("logged", MUST_BE,"no"),
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("contents", SAME_AS_CONSTRAINT),),
       
    Module(
        'plot_illu_housekeeping_line',
        SignalFile_rule.ControlFile,HousekeepingPlot,
        Constraint("preprocess", MUST_BE,'illumina'),
        Constraint("format", MUST_BE,"gct"),
        Constraint("logged", MUST_BE,"no"),
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("contents", SAME_AS_CONSTRAINT)),
    ]
