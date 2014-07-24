from Betsy.bie3 import *
import SignalFile_rule

GseaFile = DataType('GseaFile',
                    AttributeDef("contents",SignalFile_rule.CONTENTS,
                               'unspecified','unspecified'),
                    AttributeDef("permutation_type",["phenotype", "gene_set"], "phenotype","phenotype"),
                    AttributeDef("geneset_database",['computational', 'curated', 'curated:biocarta',
                        'curated:canonical', 'curated:kegg', 'curated:reactome',
                        'gene_ontology', 'gene_ontology:process', 'motif',
                        'motif:tfactor', 'positional'],
                        "gene_ontology:process","gene_ontology:process"),)
list_files = [GseaFile]
all_modules=[
    Module(
        'annotate_genes_with_gsea',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.SignalFile],GseaFile,
         Constraint("cls_format",MUST_BE,'cls',0),
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
         Constraint("format",MUST_BE,'gct',1),
         Constraint("logged",MUST_BE,'yes',1),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         Consequence("permutation_type",SET_TO_ONE_OF,["phenotype", "gene_set"]),
         Consequence("geneset_database",SET_TO_ONE_OF,['computational', 'curated', 'curated:biocarta',
                        'curated:canonical', 'curated:kegg', 'curated:reactome',
                        'gene_ontology', 'gene_ontology:process', 'motif',
                        'motif:tfactor', 'positional']))]
    


    


