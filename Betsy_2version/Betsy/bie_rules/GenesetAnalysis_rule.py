#GenesetAnalysis
import bie
import SignalFile2_rule
GenesetAnalysis=bie.DataType(
    'GenesetAnalysis',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(geneset=bie.ANYATOM, DEFAULT="no"),
    bie.Attribute(allgenes=['yes','no'], DEFAULT="no"),
    bie.Attribute(automatch=['yes','no'], DEFAULT="no"),
    )

GenesetFile = bie.DataType(
    'GenesetFile',
    bie.Attribute(filename=bie.ANYATOM,DEFAULT="", OPTIONAL=True))
list_files = [GenesetAnalysis,GenesetFile]
all_modules = [
    bie.Module(
        'score_pathway_with_geneset',
        [GenesetFile,
         SignalFile2_rule.SignalFile2(logged='yes',quantile_norm='yes',gene_center='mean',
                     gene_normalize='variance',unique_genes='high_var')],
        GenesetAnalysis(
                    geneset=bie.ANYATOM, allgenes=['yes','no'],
                    automatch=['yes','no']))]
