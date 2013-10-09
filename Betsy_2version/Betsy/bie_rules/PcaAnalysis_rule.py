#PcaAnalysis
import bie
import SignalFile_rule,SignalFile2_rule
PcaAnalysis = bie.DataType(
    'PcaAnalysis',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(pca_gene_num=bie.ANYATOM,DEFAULT='500'),
    bie.Attribute(process = ['before','after'],DEFAULT='before'))

PcaPlot = bie.DataType(
    'PcaPlot',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(pca_gene_num=bie.ANYATOM,DEFAULT='500'),
    bie.Attribute(process = ['before','after'],DEFAULT='before'))

list_files = [PcaAnalysis,PcaPlot]
all_modules = [
    bie.Module(
        'analyze_samples_pca',
        SignalFile2_rule.SignalFile2(format='tdf',logged='yes'),
        PcaAnalysis(process='after')),
    bie.Module(
        'analyze_samples_pca',
        SignalFile_rule.SignalFile(format='tdf',logged='yes'),
        PcaAnalysis(process='before')),
    bie.Module(
        'plot_sample_pca',
        [PcaAnalysis,SignalFile_rule.ClassLabelFile],
        PcaPlot(process=['before','after'])),
    
 ]
