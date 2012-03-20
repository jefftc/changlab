#classification.py

INPUTS = [
    'gse_dataset',
    'class_label_file',
    'gene_list_file',
    'gse_dataset_and_platform',
    'geo_dataset',
    'input_signal_file']

OUTPUTS = ['svm_predictions','weightedVoting']


PARAMETERS = {
    'preprocess':['rma','mas5','loess','illumina_controls',
                  'illumina','agilent','unknown_preprocess'],
    'gene_center':['mean','median','no_gene_center'],
    'gene_normalize':['variance','sum_of_squares','no_gene_normalize'],
    'quantile':['yes_quantile','no_quantile'],
    'dwd':['yes_dwd','no_dwd'],
    'shiftscale':['yes_shiftscale','no_shiftscale'],
    'combat':['yes_combat','no_combat'],
    'gene_order':['by_sample_ttest','by_gene_list','no_order'],
    'filter_fc':['yes_filter_fc','no_filter_fc'],
    'filter':['no_filter','integer between 0 and 100'],
    'format':['pcl','gct']}

DEFAULT = {
    'preprocess':'rma','gene_center':'no_gene_center',
    'gene_normalize':'no_gene_normalize','quantile':'no_quantile',
    'dwd':'no_dwd','shiftscale':'no_shiftscale',
    'combat':'no_combat','filter':'no_filter',
    'filter_fc':'no_filter_fc',
    'gene_order':'no_order'}

predicate2arguments={
    'gse_dataset':([],'[]'),
    'gse_dataset_and_platform':([],'[]'),
    'geo_dataset':(['version','unknown_version'],'[]'),
    'class_label_file':(['status','given'],'[]'),
    'gene_list_file':([],'[]'),
    'input_signal_file':(['preprocess','rma','format','tdf','is_logged','logged','has_missing_value','no_missing','status','given'],'[]')}
