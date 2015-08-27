from Betsy.bie3 import *
import BasicDataTypes as BDT
import GeneExpProcessing

GSEA_PERMUTATION = ["phenotype", "gene_set"]

GSEA_DATABASE = [
    'computational', 'curated', 'curated:biocarta', 'curated:canonical',
    'curated:kegg', 'curated:reactome', 'gene_ontology',
    'gene_ontology:process', 'motif', 'motif:tfactor', 'positional',
    ]


GseaFile = DataType(
    "GseaFile",
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
    help="Gsea file",
    )

all_data_types = [GseaFile]

all_modules=[
    ModuleNode(
        'annotate_genes_with_gsea',
        [GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile],
        GseaFile,
        Constraint("contents",CAN_BE_ANY_OF,BDT.CONTENTS,0),
        Constraint("format",MUST_BE,'gct',1),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("contents",SAME_AS_CONSTRAINT,0),
        Consequence("permutation_type", SET_TO_ONE_OF, GSEA_PERMUTATION),
        Consequence("geneset_database", SET_TO_ONE_OF, GSEA_DATABASE),
        help="annotate genes with gsea method"),
    ]
