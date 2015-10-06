from Betsy.bie3 import *
#import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS
#import GeneExpProcessing


MACSResults = DataType(
    "MACSResults",
    help="Run a MACS analysis and save the results in this folder.",
    )


all_data_types = [
    MACSResults,
    ]

all_modules = [
    ModuleNode(
        "run_MACS",
        [NGS.BamFolder, NGS.SampleGroupFile], MACSResults,
        OptionDef(
            "treatment_sample", 
            help="Comma-separated names the samples to analyze.",
            ),
        OptionDef(
            "control_sample", default="",
            help="(OPTIONAL) Comma-separated names of the samples for "
            "background.",
            ),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        ),
    ]
