#PcaAnalysis
from Betsy import bie
import SignalFile_rule,SignalFile2_rule
PcaAnalysis = bie.DataType(
    'PcaAnalysis',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(pca_gene_num=bie.ANYATOM,DEFAULT='500'),
    bie.Attribute(process = ['before','after'],DEFAULT='before'),
    bie.Attribute(contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"],DEFAULT='no'))

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
        PcaAnalysis(process='after',contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"])),
    bie.Module(
        'analyze_samples_pca',
        SignalFile_rule.SignalFile(format='tdf',logged='yes',missing_values='no'),
        PcaAnalysis(process='before',contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"])),
##    bie.Module(
##        'plot_sample_pca',
##        [PcaAnalysis(process=['before','after']),
##         SignalFile_rule.ClassLabelFile(cls_format='cls')],
##        PcaPlot(process=['before','after'])),
    bie.Module(
        'plot_sample_pca_wo_label',
        PcaAnalysis(process=['before','after']),
        PcaPlot(process=['before','after'])), 
 ]
