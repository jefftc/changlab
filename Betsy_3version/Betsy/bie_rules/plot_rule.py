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
        "unspecified", "unspecified"),
    AttributeDef("preprocess",SignalFile_rule.PREPROCESS,
                 'unknown','unknown')
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
        SignalFile_rule.SignalFile, IntensityPlot,
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        ),
                    

    Module(
        'plot_actb_line',
        SignalFile_rule.SignalFile_Impute, ActbPlot,
        Constraint("preprocess",CAN_BE_ANY_OF,SignalFile_rule.PREPROCESS),
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        ),
    

    Module(
        'plot_affy_affx_line',
        SignalFile_rule.SignalFile,ControlPlot,
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
