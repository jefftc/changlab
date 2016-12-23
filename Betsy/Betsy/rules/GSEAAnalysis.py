from Betsy.bie3 import *
import BasicDataTypes as BDT
import SignalFile

GSEA_PERMUTATION = ["phenotype", "gene_set"]

GSEA_DATABASE = [
    'computational', 'curated', 'curated:biocarta', 'curated:canonical',
    'curated:kegg', 'curated:reactome', 'gene_ontology',
    'gene_ontology:process', 'motif', 'motif:tfactor', 'positional',
    ]

GSEAResults = DataType(
    "GSEAResults",
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    AttributeDef(
        "permutation_type", GSEA_PERMUTATION, "phenotype", "phenotype",
        help="permutation type"),
    AttributeDef(
        "geneset_database", [
            'computational', 'curated', 'curated:biocarta',
            'curated:canonical', 'curated:kegg', 'curated:reactome',
            'gene_ontology', 'gene_ontology:process', 'motif',
            'motif:tfactor', 'positional'],
        "gene_ontology:process", "gene_ontology:process",
        help="geneset database"),
    AttributeDef(
        "unique_genes", ["average_genes", "high_var", "first_gene"],
        "high_var", "high_var", help="method to get unique genes"),
    help="Folder containing one or more GSEA analyses.",
    )

all_data_types = [
    GSEAResults,
    ]

all_modules=[
    ModuleNode(
        "annotate_genes_with_gsea",
        [SignalFile.SignalFile, SignalFile.ClassLabelFile],
        GSEAResults,

        OptionDef(
            "gsea_fdr_cutoff", default="0.25",
            help="Use this FDR cutoff to determine significant gene sets."),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS,0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("format", MUST_BE, 'gct', 0),
        Constraint("logged", MUST_BE, 'yes', 0),
        Constraint("annotate", MUST_BE, "yes", 0),
        Constraint(
            "unique_genes", CAN_BE_ANY_OF, 
            ["average_genes", "high_var", "first_gene"], 0),
        Consequence("unique_genes", SAME_AS_CONSTRAINT, 0),
        Consequence("permutation_type", SET_TO_ONE_OF, GSEA_PERMUTATION),
        Consequence("geneset_database", SET_TO_ONE_OF, GSEA_DATABASE),
        help="annotate genes with gsea method"),
    ]
