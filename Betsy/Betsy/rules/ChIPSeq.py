from Betsy.bie3 import *
#import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS
#import GeneExpProcessing


MACS14Results = DataType(
    "MACS14Results",
    help="Run a MACS 1.4 analysis and save the results in this folder.",
    )
MACS21Results = DataType(
    "MACS21Results",
    help="Run a MACS 2.1 analysis and save the results in this folder.",
    )
PeakSeqResults = DataType(
    "PeakSeqResults",
    help="Run a peakseq analysis and save the results in this folder.",
    )
SPPResults = DataType(
    "SPPResults",
    help="Run an SPP analysis and save the results in this folder.",
    )
HOMERTagDirectory = DataType(
    "HOMERTagDirectory",
    help="Create a tag directory for use in HOMER analysis.",
    )
HOMERResults = DataType(
    "HOMERResults",
    help="Run a HOMER analysis and save the results in this folder.",
    )

all_data_types = [
    MACS14Results,
    MACS21Results,
    PeakSeqResults,
    SPPResults,
    HOMERTagDirectory,
    HOMERResults,
    ]

all_modules = [
    ModuleNode(
        "run_MACS14",
        [NGS.BamFolder, NGS.SampleGroupFile], MACS14Results,
        OptionDef(
            "treatment_sample", 
            help="Comma-separated names the samples to analyze.",
            ),
        OptionDef(
            "control_sample", default="",
            help="(OPTIONAL) Name of the sample for background.",
            ),
        OptionDef(
            "macs_genome", 
            help="Genome.  Must be: hs, mm, ce, or dm.",
            ),
        OptionDef(
            "macs_shiftsize", default="",
            help="Number of bases to shift peaks.  1/2 of fragment length.  "
            "By default, will try to estimate.  A reasonable one is "
            "73 bp because it is 1/2 the 146 bp that is wrapped in a "
            "nucleosome."
            ),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        ),
    ModuleNode(
        "run_MACS21",
        [NGS.BamFolder, NGS.SampleGroupFile], MACS21Results,
        OptionDef(
            "macs_genome", 
            help="Genome.  Must be: hs, mm, ce, or dm.",
            ),
        OptionDef(
            "treatment_sample", 
            help="Comma-separated names the samples to analyze.",
            ),
        OptionDef(
            "control_sample", default="",
            help="(OPTIONAL) Name of the sample for background.",
            ),
        OptionDef(
            "broad_peaks", default="no",
            help='Set to "yes" to look for Broad peaks.',
            ),
        OptionDef(
            "macs_paired", default="no",
            help='Set to "yes" to do paired analysis.',
            ),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        ),
    ModuleNode(
        "run_peakseq",
        [NGS.BamFolder, NGS.SampleGroupFile], PeakSeqResults,
        OptionDef(
            "treatment_sample", 
            help="Name of the sample to analyze.",
            ),
        OptionDef(
            "control_sample", 
            help="Name of the sample for background.",
            ),
        OptionDef(
            "peakseq_fragment_length", 
            help="Average length of the fragments.",
            ),
        OptionDef(
            "mappability_file", 
            help="Mappability file.  See peakseq documentation.",
            ),
        
        Constraint("sorted", MUST_BE, "coordinate", 0),
        ),
    ModuleNode(
        "run_spp",
        [NGS.BamFolder, NGS.SampleGroupFile], SPPResults,
        OptionDef(
            "treatment_sample", 
            help="Name of the sample to analyze.",
            ),
        OptionDef(
            "control_sample", 
            help="Name of the sample for background.",
            ),
        #Constraint("sorted", MUST_BE, "coordinate", 0),
        ),
    ModuleNode(
        "make_homer_tag_directory",
        NGS.BamFolder, HOMERTagDirectory,
        ),
    ModuleNode(
        "run_homer",
        [HOMERTagDirectory, NGS.SampleGroupFile], HOMERResults,
        OptionDef(
            "treatment_sample", 
            help="Name of the sample to analyze.",
            ),
        OptionDef(
            "control_sample", default="",
            help="(OPTIONAL) Name of the sample for background.",
            ),
        # Not implemented: other styles (broad peaks, etc).
        ),
    ]
