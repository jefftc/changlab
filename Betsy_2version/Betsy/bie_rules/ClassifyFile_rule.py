#ClassifyFile
from Betsy import bie
import SignalFile_rule
import SignalFile2_rule
import PcaAnalysis_rule
ClassifyFile=bie.DataType(
    'ClassifyFile',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(classify_alg=['weighted_voting','svm','random_forest','no'], DEFAULT='no'),
    bie.Attribute(wv_minstd=bie.ANYATOM,DEFAULT=''),
    bie.Attribute(wv_feature_stat=['wv_snr', 'wv_ttest', 'wv_snr_median',
                        'wv_ttest_median',
                        'wv_snr_minstd', 'wv_ttest_minstd',
                        'wv_snr_median_minstd',
                        'wv_ttest_median_minstd'],DEFAULT='wv_snr'),
    bie.Attribute(num_features=bie.ANYATOM, DEFAULT="all"),
    bie.Attribute(svm_kernel = ['linear','polynomial','RBF','sigmoid','precomputed_kernel'],
              DEFAULT='linear'),
    bie.Attribute(loocv=['yes','no'],DEFAULT='no'),
    bie.Attribute(actual_label=['yes','no'],DEFAULT='no'))
PredictionPCAPlot = bie.DataType(
    'PredictionPCAPlot',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(classify_alg=['weighted_voting','svm','random_forest','no'],DEFAULT='no'),
    bie.Attribute(loocv=['yes','no'],DEFAULT='no'),
    bie.Attribute(actual_label=['yes','no'],DEFAULT='no'))

PredictionPlot = bie.DataType(
    'PredictionPlot',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(classify_alg=['weighted_voting','svm','random_forest','no'],DEFAULT='no'),
    bie.Attribute(loocv=['yes','no'],DEFAULT='no'),
    bie.Attribute(actual_label=['yes','no'],DEFAULT='no'))

SvmModel = bie.DataType(
    'SvmModel',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(classify_alg=['svm','no'], DEFAULT='no'),
    bie.Attribute(svm_kernel = ['linear','polynomial','RBF','sigmoid','precomputed_kernel'],
              DEFAULT='linear'))

list_files = [ClassifyFile,SvmModel,PredictionPCAPlot,PredictionPlot]

all_modules = [
    bie.Module(
       'classify_with_weighted_voting',
       [SignalFile_rule.ClassLabelFile(contents='class0,class1',cls_format='cls'),
        SignalFile2_rule.SignalFile2(contents='test',format='gct'),
        SignalFile2_rule.SignalFile2(contents='class0,class1',format='gct')],
       ClassifyFile(classify_alg='weighted_voting',
                    wv_minstd=bie.ANYATOM,
                    wv_feature_stat=['wv_snr', 'wv_ttest', 'wv_snr_median',
                        'wv_ttest_median',
                        'wv_snr_minstd', 'wv_ttest_minstd',
                        'wv_snr_median_minstd',
                        'wv_ttest_median_minstd'])),
    bie.Module(
       'classify_with_random_forest',
       [SignalFile_rule.ClassLabelFile(contents='class0,class1',cls_format='cls'),
        SignalFile2_rule.SignalFile2(contents='class0,class1,test',format='gct'),
        ],
        ClassifyFile(classify_alg='random_forest')),
    
    
     bie.Module(
       'train_svm_model',
       [SignalFile_rule.ClassLabelFile(contents='class0,class1',cls_format='cls'),
        SignalFile2_rule.SignalFile2(contents='class0,class1,test',format='gct',logged='yes')
        ],
        SvmModel(classify_alg='svm',
                svm_kernel = ['linear','polynomial','RBF','sigmoid','precomputed_kernel']
                    )), 
     bie.Module(
       'classify_with_svm',
       [SignalFile_rule.ClassLabelFile(contents='class0,class1',cls_format='cls'),
        SignalFile2_rule.SignalFile2(contents='class0,class1,test',format='gct'),
        SvmModel(classify_alg='svm')],
        ClassifyFile(classify_alg='svm',
            svm_kernel = ['linear','polynomial','RBF','sigmoid','precomputed_kernel'])),
    bie.Module(
       'run_loocv_weighted_voting',
       [SignalFile_rule.ClassLabelFile(contents='class0,class1',cls_format='cls'),
        SignalFile2_rule.SignalFile2(contents='class0,class1',format='gct')],
        ClassifyFile(classify_alg='weighted_voting',
                    wv_minstd=bie.ANYATOM,
                    wv_feature_stat=['wv_snr', 'wv_ttest', 'wv_snr_median',
                        'wv_ttest_median',
                        'wv_snr_minstd', 'wv_ttest_minstd',
                        'wv_snr_median_minstd',
                        'wv_ttest_median_minstd'],
                     loocv='yes',actual_label='no'
                    )),
    bie.Module(
       'run_loocv_svm',
       [SignalFile_rule.ClassLabelFile(contents='class0,class1',cls_format='cls'),
        SignalFile2_rule.SignalFile2(contents='class0,class1',format='gct')],
        ClassifyFile(classify_alg='svm',
                    svm_kernel = ['linear','polynomial','RBF','sigmoid','precomputed_kernel'],
                     loocv='yes',actual_label='no'
                    )),
    bie.Module(
       'run_loocv_random_forest',
       [SignalFile_rule.ClassLabelFile(contents='class0,class1',cls_format='cls'),
        SignalFile2_rule.SignalFile2(contents='class0,class1',format='gct')],
        ClassifyFile(classify_alg='random_forest',
                     loocv='yes',actual_label='no'
                    )),
    bie.Module(
        'evaluate_prediction',
        [SignalFile_rule.ClassLabelFile(contents='test',cls_format='cls'),
         ClassifyFile(loocv='no',actual_label='no')],
         ClassifyFile(loocv='no',actual_label='yes')),

    bie.Module(
        'plot_prediction',
         ClassifyFile(classify_alg=['weighted_voting','svm','random_forest'],
                        actual_label=['yes','no'],loocv=['yes','no']),
         PredictionPlot(classify_alg=['weighted_voting','svm','random_forest'],
                        actual_label=['yes','no'],loocv=['yes','no'])),

    bie.Module(
        'plot_sample_pca_with_predictions',
         [ClassifyFile(classify_alg=['svm','weighted_voting',
                                     'random_forest'],actual_label=['yes','no'],
                       loocv=['yes','no']),
          PcaAnalysis_rule.PcaAnalysis(process='after',contents='test')],
          PredictionPCAPlot(classify_alg=['svm','weighted_voting',
                                     'random_forest'],
                        actual_label=['yes','no'],loocv=['yes','no'])),
    ]
