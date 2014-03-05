#ClassifyFile
from Betsy.bie3 import *
import SignalFile_rule
import PcaAnalysis_rule

ClassifyFile=DataType(
    'ClassifyFile',
    AttributeDef("classify_alg",['weighted_voting','svm','random_forest','no'], 'no','no'),
    AttributeDef('wv_feature_stat',['wv_snr', 'wv_ttest', 'wv_snr_median',
                        'wv_ttest_median',
                        'wv_snr_minstd', 'wv_ttest_minstd',
                        'wv_snr_median_minstd',
                        'wv_ttest_median_minstd'],'wv_snr','wv_snr'),
    AttributeDef('svm_kernel', ['linear','polynomial','RBF','sigmoid','precomputed_kernel'],
              'linear','linear'),
    AttributeDef('loocv',['yes','no'],'no','no'),
    AttributeDef('actual_label',['yes','no'],'no','no'))

PredictionPCAPlot = DataType(
    'PredictionPCAPlot',
    AttributeDef('classify_alg',['weighted_voting','svm','random_forest','no'],'no','no'),
    AttributeDef('loocv',['yes','no'],'no','no'),
    AttributeDef('actual_label',['yes','no'],'no','no'))

PredictionPlot = DataType(
    'PredictionPlot',
    AttributeDef('classify_alg',['weighted_voting','svm','random_forest','no'],'no','no'),
    AttributeDef('loocv',['yes','no'],'no','no'),
    AttributeDef('actual_label',['yes','no'],'no','no'))

SvmModel = DataType(
    'SvmModel',
    AttributeDef('classify_alg',['svm','no'], 'no','no'),
    AttributeDef('svm_kernel',['linear','polynomial','RBF','sigmoid','precomputed_kernel'],
              'linear','linear'))

list_files = [ClassifyFile,SvmModel,PredictionPCAPlot,PredictionPlot]

all_modules = [
    Module(
        "merge_files_for_classification",
        [SignalFile_rule.PrettySignalFile,SignalFile_rule.PrettySignalFile],SignalFile_rule.PrettySignalFile,
        Constraint('contents',MUST_BE,"class0,class1",0),
        Constraint('format',MUST_BE,'gct',0),
        Constraint('logged',MUST_BE,"yes",0),
        Constraint('missing_values',MUST_BE,'no',0),
        Constraint('contents',MUST_BE,"test",1),
        Constraint('format',MUST_BE,'gct',1),
        Constraint('logged',MUST_BE,"yes",1),
        Constraint('missing_values',MUST_BE,'no',1),
        Consequence('contents',SET_TO,"class0,class1,test"),
        Consequence('format',SAME_AS_CONSTRAINT,0),
        Consequence('logged',SAME_AS_CONSTRAINT,0),
        Consequence('missing_values',SAME_AS_CONSTRAINT,0),
        Constraint('psf_processing_step',MUST_BE,'processed',0),
        Constraint('psf_processing_step',MUST_BE,'processed',1),
        Consequence('psf_processing_step',SAME_AS_CONSTRAINT,0)
        ),
    Module(
       'classify_with_weighted_voting',
       [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile,
        SignalFile_rule.PrettySignalFile],ClassifyFile,
       UserInputDef('num_features',10),
       UserInputDef('wv_minstd',1),
       Constraint('contents',MUST_BE,'class0,class1',0),
       Constraint('cls_format',MUST_BE,'cls',0),
       Constraint("contents",MUST_BE,'test',1),
       Constraint("format",MUST_BE,'gct',1),
       Constraint("contents",MUST_BE,'class0,class1',2),
       Constraint("format",MUST_BE,'gct',2),
       Consequence("classify_alg",SET_TO,'weighted_voting'),
       Consequence('wv_feature_stat',SET_TO_ONE_OF,['wv_snr','wv_ttest', 'wv_snr_median',
                        'wv_ttest_median',
                        'wv_snr_minstd', 'wv_ttest_minstd',
                        'wv_snr_median_minstd',
                        'wv_ttest_median_minstd'])),
    Module(
       'classify_with_random_forest',
       [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile],ClassifyFile,
       Constraint("contents",MUST_BE,'class0,class1',0),
       Constraint("cls_format",MUST_BE,'cls',0),
       Constraint("contents",MUST_BE,'class0,class1,test',1),
       Constraint("format",MUST_BE,'gct',1),
       Consequence('classify_alg',SET_TO,'random_forest')),
    
    
     Module(
       'train_svm_model',
       [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile],SvmModel,
       Constraint("contents",MUST_BE,'class0,class1',0),
       Constraint("cls_format",MUST_BE,'cls',0),
       Constraint("contents",MUST_BE,'class0,class1,test',1),
       Constraint("format",MUST_BE,'gct',1),
       Constraint("logged",MUST_BE,'yes',1),
       Consequence("classify_alg",SET_TO,'svm'),
       Consequence("svm_kernel",SET_TO_ONE_OF, ['linear','polynomial',
                                                'RBF','sigmoid','precomputed_kernel'])), 
     Module(
       'classify_with_svm',
       [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile,SvmModel],ClassifyFile,
       Constraint("contents",MUST_BE,'class0,class1',0),
       Constraint("cls_format",MUST_BE,'cls',0),
       Constraint("contents",MUST_BE,'class0,class1,test',1),
       Constraint("format",MUST_BE,'gct',1),
       Constraint("classify_alg",MUST_BE,'svm',2),
       Consequence("classify_alg",SAME_AS_CONSTRAINT,2),
       Consequence("svm_kernel",SET_TO_ONE_OF,
                   ['linear','polynomial','RBF','sigmoid','precomputed_kernel'])),
    Module(
       'run_loocv_weighted_voting',
       [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile],ClassifyFile,
       UserInputDef('num_features',10),
       UserInputDef('wv_minstd',1),
       Constraint("contents",MUST_BE,'class0,class1',0),
       Constraint("cls_format",MUST_BE,'cls',0),
       Constraint("contents",MUST_BE,'class0,class1',1),
       Constraint("format",MUST_BE,'gct',1),
       Consequence("classify_alg",SET_TO,'weighted_voting'),
       Consequence('wv_feature_stat',SET_TO_ONE_OF,['wv_snr', 'wv_ttest', 'wv_snr_median',
                        'wv_ttest_median',
                        'wv_snr_minstd', 'wv_ttest_minstd',
                        'wv_snr_median_minstd',
                        'wv_ttest_median_minstd']),
       Consequence('loocv',SET_TO,'yes'),
       Consequence('actual_label',SET_TO,'no')),
    Module(
       'run_loocv_svm',
       [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile],ClassifyFile,
       Constraint('contents',MUST_BE,'class0,class1',0),
       Constraint("cls_format",MUST_BE,'cls',0),
       Constraint("contents",MUST_BE,'class0,class1',1),
       Constraint("format",MUST_BE,'gct',1),
       Consequence('classify_alg',SET_TO,'svm'),
       Consequence('svm_kernel',SET_TO_ONE_OF, ['linear','polynomial',
                                                'RBF','sigmoid','precomputed_kernel']),
       Consequence("loocv",SET_TO,'yes'),
       Consequence("actual_label",SET_TO,'no')),
    Module(
       'run_loocv_random_forest',
       [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile],ClassifyFile,
       Constraint('contents',MUST_BE,'class0,class1',0),
       Constraint("cls_format",MUST_BE,'cls',0),
       Constraint("contents",MUST_BE,'class0,class1',1),
       Constraint("format",MUST_BE,'gct',1),
       Consequence("classify_alg",SET_TO,'random_forest'),
       Consequence("loocv",SET_TO,'yes'),
       Consequence("actual_label",SET_TO,'no')),
                  
    Module(
        'evaluate_prediction',
        [SignalFile_rule.ClassLabelFile,ClassifyFile],ClassifyFile,
        Constraint("contents",MUST_BE,'test',0),
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("loocv",MUST_BE,'no',1),
        Constraint("actual_label",MUST_BE,'no',1),
        Consequence("loocv",SAME_AS_CONSTRAINT,1),
        Consequence("actual_label",SET_TO,'yes'),
        DefaultAttributesFrom(1)),

    Module(
        'plot_prediction',
        ClassifyFile,PredictionPlot,
        Constraint("classify_alg",CAN_BE_ANY_OF,['weighted_voting','svm','random_forest']),
        Constraint("actual_label",CAN_BE_ANY_OF,['yes','no']),
        Constraint("loocv",CAN_BE_ANY_OF,['yes','no']),
        Consequence("classify_alg",SAME_AS_CONSTRAINT),
        Consequence("actual_label",SAME_AS_CONSTRAINT),
        Consequence("loocv",SAME_AS_CONSTRAINT)),

    Module(
        'plot_sample_pca_with_predictions',
         [ClassifyFile,PcaAnalysis_rule.PcaAnalysis],PredictionPCAPlot,
        Constraint('classify_alg',CAN_BE_ANY_OF,['svm','weighted_voting',
                                     'random_forest'],0),
        Constraint("actual_label",CAN_BE_ANY_OF,['yes','no'],0),
        Constraint("loocv",CAN_BE_ANY_OF,['yes','no'],0),
        Constraint("contents",MUST_BE,'test',1),
        Consequence("classify_alg",SAME_AS_CONSTRAINT,0),
        Consequence("actual_label",SAME_AS_CONSTRAINT,0),
        Consequence("loocv",SAME_AS_CONSTRAINT,0)),
    ]
