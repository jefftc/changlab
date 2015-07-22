#Classification
from Betsy.bie3 import *
import GeneExpProcessing
import PcaAnalysis

ClassifyFile = DataType(
    'ClassifyFile',
    AttributeDef(
        "classify_alg", ['weighted_voting', 'svm', 'random_forest'],
        'svm', 'svm', help="classify algorithm"),
    AttributeDef(
        'wv_feature_stat',
        ['wv_snr', 'wv_ttest', 'wv_snr_median', 'wv_ttest_median',
         'wv_snr_minstd', 'wv_ttest_minstd', 'wv_snr_median_minstd',
         'wv_ttest_median_minstd'], 'wv_snr', 'wv_snr',
        help="weighted voting feature stat"),
    AttributeDef(
        'svm_kernel',
        ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel'],
        'linear', 'linear',
        help="svm kernel"),
    AttributeDef('loocv', ['yes', 'no'], 'no', 'no', help="loocv yes or not"),
    AttributeDef(
        'actual_label', ['yes', 'no'], 'no', 'no',
        help="compared with the actual label or not"),
    help="Classify result file")

PredictionPCAPlot = DataType(
    'PredictionPCAPlot',
    AttributeDef(
        'classify_alg', ['weighted_voting', 'svm', 'random_forest', 'no'],
        'no', 'no', help="classify algorithm"),
    AttributeDef(
        'loocv', ['yes', 'no'], 'no', 'no', help="loocv yes or not"),
    AttributeDef(
        'actual_label', ['yes', 'no'], 'no', 'no',
        help="compared with the actual label or not"),
    help="PCA plot labeled with predication result")

PredictionPlot = DataType(
    'PredictionPlot',
    AttributeDef(
        'classify_alg', ['weighted_voting', 'svm', 'random_forest'],
        'svm', 'svm', help="classify algorithm"),
    AttributeDef(
        'loocv', ['yes', 'no'], 'no', 'no', help="loocv yes or not"),
    AttributeDef(
        'actual_label', ['yes', 'no'], 'no', 'no',
        help="compared with the actual label or not"),
    help='Prediction plot')

SvmModel = DataType(
    'SvmModel',
    AttributeDef(
        'classify_alg', ['svm', 'no'], 'no', 'no', help="classify algorithm"),
    AttributeDef(
        'svm_kernel',
        ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel'],
        'linear', 'linear', help="svm kernel"),
    help='svm model file')

ClassifyReportFile = DataType(
    'ClassifyReportFile', help="Report file for Classify report")

list_files = [
    ClassifyFile,
    SvmModel,
    PredictionPCAPlot,
    PredictionPlot,
    ClassifyReportFile,
    ]

all_modules = [
    Module(
        "merge_files_for_classification",
        [GeneExpProcessing.SignalFile, GeneExpProcessing.SignalFile],
        GeneExpProcessing.SignalFile,
        Constraint('contents', MUST_BE, "class0,class1", 0),
        Constraint('format', MUST_BE, 'gct', 0),
        Constraint('logged', MUST_BE, "yes", 0),
        Constraint('contents', MUST_BE, "test", 1),
        Constraint('format', SAME_AS, 0, 1),
        Constraint('logged', SAME_AS, 0, 1),
        Consequence('contents', SET_TO, "class0,class1,test"),
        Consequence('format', SAME_AS_CONSTRAINT, 0),
        Consequence('logged', SAME_AS_CONSTRAINT, 0),
        Constraint(
            'gene_order', CAN_BE_ANY_OF,
            ["no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr",
             'diff_ttest', 'diff_sam', 'diff_ebayes', 'diff_fold_change'],
            0),
        Constraint('gene_order', MUST_BE, "no", 1),
        Consequence('gene_order', SAME_AS_CONSTRAINT, 0),
        Constraint('num_features', CAN_BE_ANY_OF, ["no", "yes"], 0),
        Constraint('num_features', MUST_BE, "no", 1),
        Consequence('num_features', SAME_AS_CONSTRAINT, 0),
        # Why two DefaultAttributesFrom?
        DefaultAttributesFrom(0),
        DefaultAttributesFrom(1),
        help="merge two files for classification"),
    
    Module(
        'classify_with_weighted_voting',
        [GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile,
         GeneExpProcessing.SignalFile], ClassifyFile,
        OptionDef(
            'wv_num_features', 10,
            help="number of features for weighted voting"),
        OptionDef(
            'wv_minstd', 1,
            help="minstd for weighted voting"),
        Constraint('contents', MUST_BE, 'class0,class1', 0),
        Constraint('cls_format', MUST_BE, 'cls', 0),
        Constraint("contents", MUST_BE, 'test', 1),
        Constraint("format", MUST_BE, 'gct', 1),
        Constraint("contents", MUST_BE, 'class0,class1', 2),
        Constraint("format", MUST_BE, 'gct', 2),
        Consequence("classify_alg", SET_TO, 'weighted_voting'),
        Consequence(
            'wv_feature_stat', SET_TO_ONE_OF,
            ['wv_snr', 'wv_ttest', 'wv_snr_median', 'wv_ttest_median',
             'wv_snr_minstd', 'wv_ttest_minstd', 'wv_snr_median_minstd',
             'wv_ttest_median_minstd']),
        Consequence(
            'svm_kernel', SET_TO_ONE_OF,
            ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel']),
        help="classify with weighted voting method"),
    
    Module(
        'classify_with_random_forest',
        [GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile],
        ClassifyFile,
        
        Constraint("contents", MUST_BE, 'class0,class1', 0),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", MUST_BE, 'class0,class1,test', 1),
        Constraint("format", MUST_BE, 'gct', 1),
        #Constraint("logged",MUST_BE,'yes',1),
        
        Consequence('classify_alg', SET_TO, 'random_forest'),
        Consequence(
            'wv_feature_stat', SET_TO_ONE_OF,
            ['wv_snr', 'wv_ttest', 'wv_snr_median', 'wv_ttest_median',
             'wv_snr_minstd', 'wv_ttest_minstd', 'wv_snr_median_minstd',
             'wv_ttest_median_minstd']),
        Consequence(
            'svm_kernel', SET_TO_ONE_OF,
            ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel']),
        help="classify with random forest method"),
    Module(
        'train_svm_model',
        [GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile],
        SvmModel,
        Constraint("contents", MUST_BE, 'class0,class1', 0),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", MUST_BE, 'class0,class1,test', 1),
        Constraint("format", MUST_BE, 'gct', 1),
        #Constraint("logged",MUST_BE,'yes',1),
        Consequence("classify_alg", SET_TO, 'svm'),
        Consequence(
            "svm_kernel", SET_TO_ONE_OF,
            ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel']),
        help="train data using svm method"),
    Module(
        'classify_with_svm',
        [GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile,
         SvmModel],
        ClassifyFile,
        Constraint("contents", MUST_BE, 'class0,class1', 0),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", MUST_BE, 'class0,class1,test', 1),
        Constraint("format", MUST_BE, 'gct', 1),
        Constraint("classify_alg", MUST_BE, 'svm', 2),
        Constraint(
            "svm_kernel", CAN_BE_ANY_OF,
            ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel'],
            2),
        Consequence("classify_alg", SAME_AS_CONSTRAINT, 2),
        Consequence("svm_kernel", SAME_AS_CONSTRAINT, 2),
        Consequence(
            'wv_feature_stat', SET_TO_ONE_OF,
            ['wv_snr', 'wv_ttest', 'wv_snr_median', 'wv_ttest_median',
             'wv_snr_minstd', 'wv_ttest_minstd', 'wv_snr_median_minstd',
             'wv_ttest_median_minstd']),
        help="classify with svm method"),
    Module(
        'run_loocv_weighted_voting',
        [GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile],
        ClassifyFile,
        OptionDef(
            'wv_num_features', 10,
            help="number of features for weighted voting"),
        OptionDef(
            'wv_minstd', 1, help="minstd for weighted voting"),
        Constraint("contents", MUST_BE, 'class0,class1', 0),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", MUST_BE, 'class0,class1', 1),
        Constraint("format", MUST_BE, 'gct', 1),
        #Constraint("logged",MUST_BE,'yes',1),
        Consequence("classify_alg", SET_TO, 'weighted_voting'),
        Consequence(
            'wv_feature_stat', SET_TO_ONE_OF,
            ['wv_snr', 'wv_ttest', 'wv_snr_median', 'wv_ttest_median',
             'wv_snr_minstd', 'wv_ttest_minstd', 'wv_snr_median_minstd',
             'wv_ttest_median_minstd']),
        Consequence(
            'svm_kernel', SET_TO_ONE_OF,
            ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel']),
        Consequence('loocv', SET_TO, 'yes'),
        Consequence('actual_label', SET_TO, 'no'),
        help="run loocv in weighted voting method"),
    Module(
        'run_loocv_svm',
        [GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile],
        ClassifyFile,
        Constraint('contents', MUST_BE, 'class0,class1', 0),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", MUST_BE, 'class0,class1', 1),
        Constraint("format", MUST_BE, 'gct', 1),
        #Constraint("logged",MUST_BE,'yes',1),
        Consequence('classify_alg', SET_TO, 'svm'),
        Consequence(
            'svm_kernel', SET_TO_ONE_OF,
            ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel']),
        Consequence("loocv", SET_TO, 'yes'),
        Consequence("actual_label", SET_TO, 'no'),
        Consequence(
            'wv_feature_stat', SET_TO_ONE_OF,
            ['wv_snr', 'wv_ttest', 'wv_snr_median', 'wv_ttest_median',
             'wv_snr_minstd', 'wv_ttest_minstd', 'wv_snr_median_minstd',
             'wv_ttest_median_minstd']),
        help="run loocv in svm method"),
    Module(
        'run_loocv_random_forest',
        [GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile],
        ClassifyFile,
        
        Constraint('contents', MUST_BE, 'class0,class1', 0),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", MUST_BE, 'class0,class1', 1),
        Constraint("format", MUST_BE, 'gct', 1),
        #Constraint("logged",MUST_BE,'yes',1),
        
        Consequence("classify_alg", SET_TO, 'random_forest'),
        Consequence("loocv", SET_TO, 'yes'),
        Consequence("actual_label", SET_TO, 'no'),
        Consequence(
            'wv_feature_stat', SET_TO_ONE_OF,
            ['wv_snr', 'wv_ttest', 'wv_snr_median', 'wv_ttest_median',
             'wv_snr_minstd', 'wv_ttest_minstd', 'wv_snr_median_minstd',
             'wv_ttest_median_minstd']),
        Consequence(
            'svm_kernel', SET_TO_ONE_OF,
            ['linear', 'polynomial', 'RBF', 'sigmoid', 'precomputed_kernel']),
        help="run loocv in random forest method"),
    Module(
        'evaluate_prediction',
        [GeneExpProcessing.ClassLabelFile, ClassifyFile], ClassifyFile,
        Constraint("contents", MUST_BE, 'test', 0),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("loocv", MUST_BE, 'no', 1),
        Constraint("actual_label", MUST_BE, 'no', 1),
        Consequence("loocv", SAME_AS_CONSTRAINT, 1),
        Consequence("actual_label", SET_TO, 'yes'),
        DefaultAttributesFrom(1),
        help="evalutate prediction with actual label"),
    Module(
        'plot_prediction', ClassifyFile, PredictionPlot,
        Constraint(
            "classify_alg", CAN_BE_ANY_OF,
            ['weighted_voting', 'svm', 'random_forest']),
        Constraint("actual_label", CAN_BE_ANY_OF, ['yes', 'no']),
        Constraint("loocv", CAN_BE_ANY_OF, ['yes', 'no']),
        Consequence("classify_alg", SAME_AS_CONSTRAINT),
        Consequence("actual_label", SAME_AS_CONSTRAINT),
        Consequence("loocv", SAME_AS_CONSTRAINT),
        help="plot prediction result"),
    Module(
        'plot_sample_pca_with_predictions',
        [ClassifyFile, PcaAnalysis.PcaAnalysis],
        PredictionPCAPlot,
        Constraint(
            'classify_alg', CAN_BE_ANY_OF,
            ['svm', 'weighted_voting', 'random_forest'], 0),
        Constraint("actual_label", CAN_BE_ANY_OF, ['yes', 'no'], 0),
        Constraint("loocv", MUST_BE, 'no', 0),
        Constraint("contents", MUST_BE, 'test', 1),
        Consequence("classify_alg", SAME_AS_CONSTRAINT, 0),
        Consequence("actual_label", SAME_AS_CONSTRAINT, 0),
        Consequence("loocv", SAME_AS_CONSTRAINT, 0),
        help="plot sample pca plot labeled with prediction result"),
    Module(
        'make_classify_report',
        [GeneExpProcessing.SignalFile, ClassifyFile, ClassifyFile,
         PredictionPlot, PredictionPlot, PredictionPCAPlot, ClassifyFile,
         ClassifyFile, PredictionPlot, PredictionPlot, PredictionPCAPlot],
        ClassifyReportFile,
        Constraint('contents', MUST_BE, 'class0,class1,test', 0),
        Constraint('logged', MUST_BE, 'yes', 0),
        Constraint("format", MUST_BE, 'gct', 0),
        Constraint("classify_alg", MUST_BE, 'svm', 1),
        Constraint("actual_label", MUST_BE, 'yes', 1),
        Constraint("classify_alg", MUST_BE, 'svm', 2),
        Constraint("actual_label", MUST_BE, 'no', 2),
        Constraint("loocv", MUST_BE, 'yes', 2),
        Constraint("classify_alg", MUST_BE, 'svm', 3),
        Constraint("actual_label", MUST_BE, 'yes', 3),
        Constraint("loocv", MUST_BE, 'no', 3),
        Constraint("classify_alg", MUST_BE, 'svm', 4),
        Constraint("actual_label", MUST_BE, 'no', 4),
        Constraint("loocv", MUST_BE, 'yes', 4),
        Constraint("classify_alg", MUST_BE, 'svm', 5),
        Constraint("actual_label", MUST_BE, 'yes', 5),
        Constraint("loocv", MUST_BE, 'no', 5),
        Constraint("classify_alg", MUST_BE, 'weighted_voting', 6),
        Constraint("actual_label", MUST_BE, 'yes', 6),
        Constraint("loocv", MUST_BE, 'no', 6),
        Constraint("classify_alg", MUST_BE, 'weighted_voting', 7),
        Constraint("actual_label", MUST_BE, 'no', 7),
        Constraint("loocv", MUST_BE, 'yes', 7),
        Constraint("classify_alg", MUST_BE, 'weighted_voting', 8),
        Constraint("actual_label", MUST_BE, 'yes', 8),
        Constraint("loocv", MUST_BE, 'no', 8),
        Constraint("classify_alg", MUST_BE, 'weighted_voting', 9),
        Constraint("actual_label", MUST_BE, 'no', 9),
        Constraint("loocv", MUST_BE, 'yes', 9),
        Constraint("classify_alg", MUST_BE, 'weighted_voting', 10),
        Constraint("actual_label", MUST_BE, 'yes', 10),
        Constraint("loocv", MUST_BE, 'no', 10),
        help="make classification report"),
    ]
