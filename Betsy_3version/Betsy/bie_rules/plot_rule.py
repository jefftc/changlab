# SignalFile
#from Betsy import bie3
from Betsy.bie3 import *
import SignalFile_rule

IntensityPlot = DataType(
    'IntensityPlot',
    AttributeDef(
        "contents",
        SignalFile_rule.CONTENTS,
        "unspecified", "unspecified",help='contents'),
    help="Intensity plot file"
    )
    
ActbPlot = DataType(
    'ActbPlot',
    AttributeDef(
        "contents",
        SignalFile_rule.CONTENTS,
        "unspecified", "unspecified",help='contents'),
    AttributeDef("preprocess",SignalFile_rule.PREPROCESS,
                 'unknown','unknown',help="preprocess method"),
    help="Actb plot file"
    )
Hyb_barPlot = DataType(
    'Hyb_barPlot',
    AttributeDef(
        "contents",
        SignalFile_rule.CONTENTS,
        "unspecified", "unspecified",help='contents'),
    help="Hyb bar plot file")
ControlPlot = DataType(
    'ControlPlot',
    AttributeDef(
        "contents",
        SignalFile_rule.CONTENTS,
        "unspecified", "unspecified",help='contents'),
    help="control plot file")
BiotinPlot = DataType(
    'BiotinPlot',
    AttributeDef(
        "contents",
        SignalFile_rule.CONTENTS,
        "unspecified", "unspecified",help='contents'),
    help="Biotin plot file")
HousekeepingPlot = DataType(
    'HousekeepingPlot',
    AttributeDef(
        "contents",
        SignalFile_rule.CONTENTS,
        "unspecified", "unspecified",help='contents'),
    help="Housekeeping plot file")

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
            "contents", CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot intensity boxplot"
        ),
                    

    Module(
        'plot_actb_line',
        SignalFile_rule.SignalFile_Impute, ActbPlot,
        Constraint("preprocess",CAN_BE_ANY_OF,SignalFile_rule.PREPROCESS),
        Constraint("missing_values",MUST_BE,'no'),
        Constraint(
            "contents", CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        help="plot actb line"
        ),
    

    Module(
        'plot_affy_affx_line',
        SignalFile_rule.SignalFile,ControlPlot,
        Constraint("annotate", MUST_BE,'no'),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="plot affy affx line"
        ),
    Module(
       'plot_illu_hyb_bar',
       SignalFile_rule.ControlFile,Hyb_barPlot,
       Constraint("preprocess", MUST_BE,"illumina"),
       Constraint("format", MUST_BE,"gct"),
       Constraint("logged", MUST_BE,"no"),
       Constraint(
            "contents", CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
       help="plot illumina affx line"),
    Module(
        'plot_illu_biotin_line',
        SignalFile_rule.ControlFile,BiotinPlot,
        Constraint("preprocess", MUST_BE,'illumina'),
        Constraint("format", MUST_BE,"gct"),
        Constraint("logged", MUST_BE,"no"),
        Constraint(
            "contents", CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot illumina biotin line"),
       
    Module(
        'plot_illu_housekeeping_line',
        SignalFile_rule.ControlFile,HousekeepingPlot,
        Constraint("preprocess", MUST_BE,'illumina'),
        Constraint("format", MUST_BE,"gct"),
        Constraint("logged", MUST_BE,"no"),
        Constraint(
            "contents", CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot illumina housekeeping line"),
    ]
