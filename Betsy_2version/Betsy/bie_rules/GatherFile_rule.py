#GatherFile
import bie
import SignalFile_rule
GatherFile = bie.DataType(
    'GatherFile',
    bie.Attribute(network=['yes','no'],DEFAULT='no'),
    bie.Attribute(homologs=['yes','no'],DEFAULT='no'),
    bie.Attribute(annot_type=['kegg','gene_ontology','transfac'],DEFAULT='kegg'),
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True))
list_files = [GatherFile]
all_modules = [
        bie.Module(
        'annotate_genes_with_gather',
         SignalFile_rule.GeneListFile,
         GatherFile)]
