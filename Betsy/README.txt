Preprocess Usage

1)Betsy can preprocess with rma or mas5 for Affymetrix cel data. 

	Examples:
	A. When given a GSE ID from the geo, and want to do rma preprocess, 
	   the command is

	protocol_engine.py  
      --protocol 'normalize_file' 
      --input 'gse_id:<GSE21947>'  
      --parameters 'preprocess:rma' 
     	------------------------------------------------------------
	B. When given a GSE ID and GPL platform number from the geo, 
	   and want to do rma preprocess, 
	   the command is

	protocol_engine.py  
      --protocol 'normalize_file' 
      --input 'gse_id_and_platform:<GSE17907,GPL570>'  
      --parameters 'preprocess:rma' 
     
	------------------------------------------------------------
	C. When given a folder contains cel file, the command is

	protocol_engine.py  
      --protocol 'normalize_file' 
      --input 'cel_files:<./GSE21947>'     	
      --parameters 'preprocess:rma' 
	------------------------------------------------------------

	The output will be tdf format with annotations, 
      without missing value and logged with base 2.
      Also the result folder will contain pca_plot,intensity_plot,actb_plot,control_plot.

2)Betsy can preprocess with illumina for illumina idat files.

	Example:
	When given a folder contains idat files, the command is

	protocol_engine.py  
      --protocol 'normalize_file' 
      --input 'idat_files:<6991010018>'   
      --parameters 'has_missing_value:zero_fill'
      --parameters 'preprocess:illumina'

	The output will be tdf format with annotations, 
      without missing value and logged with base 2.
      Also the result folder will contain pca_plot,intensity_plot,
	hyb_bar_plot,actb_plot,biotin_plot,control file.


=============================================================================
Process Usage
1) Betsy can do preprocessdataset,log,unlog,gene_filter,quantile,combat, shiftscale,dwd,gene_center,
    gene_normalize,gene_order for signal files
    Example:
    When given a signal_file,do preprocessdataset, quantile, gene_center=mean,   gene_normalize=variance.
    The command is
 
    protocol_engine.py
    --protocol 'normalize_file'
    --input 'input_signal_file:<path of the file>'
    --parameters 'predataset:yes_predataset'
    --parameters 'quantile:yes_quantile'
    --parameters 'gene_center:mean'
    --parameters 'gene_normalize:variance'
    --parameters 'has_missing_value:zero_fill'

     The result folder will contain output signal_file 
     with annotation,pca_plot,intensity_plot,actb_plot,control_plot.
    --------------------------------------------------------------------
    When given a signal_file and gene_list_file, fill the missing value with median, 
    and quantile, gene_center=mean and extract the genes in gene_list_file from signal_file. 
    The output format is tdf. The command is
    protocol_engine.py
    --protocol 'normalize_file'
    --input 'input_signal_file:<path of the file>'
    --input 'gene_list_file:<path of the file>'
    --parameters 'quantile:yes_quantile'
    --parameters 'gene_center:mean'
    --parameters 'has_missing_value:median_fill'
    --parameters 'gene_order:by_gene_list'

    The result folder will contain output signal_file 
    with annotation,pca_plot,intensity_plot,atcb_plot,control_plot.
     --------------------------------------------------------------------
2) Betsy can merge several signal_file or split a signal_file into several signal_file 
    according to the labels.
    Example:
    When given two signal_file, one contains positive sample, one contains negative sample, 
    merge them together. Then command is:
    protocol_engine.py
    --protocol 'normalize_file'
    --input 'input_signal_file:pos:/path of the file'
    --input 'input_signal_file:neg:/path of the file'
    --parameters 'has_missing_value:zero_fill'
    --parameters 'contents:pos,neg'

     The result folder will contain output signal_file
     with annotations,pca_plot,intensity_plot,actb_plot,control_plot
    --------------------------------------------------------------------
    When given a signal_file and its class_label_file, split them into two signal_file 
    according to the labels. The command is:
    protocol_engine.py
    --protocol 'normalize_file'
    --input 'input_signal_file:pos,neg:<path of the file>'
    --input 'class_label_file:pos,neg:<path of the file>'
    --parameters 'has_missing_value:zero_fill'
    --parameters 'contents:pos'
    The result folder will contain output signal_file contains positive samples,
    pca_plot,intensity_plot,actb_plot,control_plot.

=============================================================================
Clustering Usage
Betsy can do clustering for a signal_file and plot the heatmap.
    Example:
    When given a signal_file, gene_center=mean, gene_normalize=variance, 
    cluster the genes by hierarchical method and plot the heatmap, 
    set the heatmap size as x=200,y=1.
    The command is:
    protocol_engine.py
    --protocol 'cluster_genes'
    --input 'input_signal_file:<path of the file>'
    --parameters 'has_missing_value:zero_fill'
    --parameters 'gene_center:mean'
    --parameters 'gene_normalize:variance'
    --parameters 'cluster_alg:hierarchical'
    --parameters 'hm_width:200'
    --parameters 'hm_height:1'
 
    The result folder will contain a clustering file and a png file showing the heatmap.
    --------------------------------------------------------------------
    When given idat files,preprocess=illumina, rename sample names, select genes according to fdr value,gene_center=mean,gene_normalize=variance,clustering and then make heatmap.
    The command is:
    protocol_engine.py
    --protocol 'cluster_genes'
    --input 'idat_files:<path of the file>'
    --input 'class_label_file:<path of the file>'
    --input 'rename_list_file:<path of the file>'
    --parameters 'has_missing_value:zero_fill'
    --parameters 'gene_center:mean'
    --parameters 'gene_normalize:variance'
    --parameters 'cluster_alg:hierarchical'
    --parameters 'hm_width:20'
    --parameters 'hm_height:1'
    --parameters 'rename_sample:yes_rename'
    --parameters 'gene_order:t_test_fdr'

    The result folder will contain a clustering file and a png file showing the heatmap.

=============================================================================
Classification Usage
Betsy can do classification with svm and weighted_voting method for dataset. 
Also it can leave one out cross validation by these two methods.
    Example:
    When given three signal_files, one contains negative samples, one contains positive 
    samples and one contains test samples.Classify the test sample. For the svm training model, 
    set the kernel function as linear. For the weighted-voting method, the number of features 
    is set to 50. The command is:
    protocol_engine.py
    --protocol 'classification'
    --input 'input_signal_file:train:<path of the file>'
    --input 'class_label_file:train:<path of the file>'
    --input 'input_signal_file:test:<path of the file>'
    --parameters 'has_missing_value:zero_fill'
    --parameters 'traincontents:train'
    --parameters 'testcontents:test'
    --parameters 'svm_kernel:linear'
    --parameters 'num_features:50' 
    The result folder contains one predication result for svm, one predication result for 
    weighted_voting,one predication result for leave one out cross validation on svm and one 
    leave one out cross validation on weighted_voting.

=============================================================================
Differential expressed genes analysis usage
Betsy can do the differential expressed genes analysis for signal_files.
Example:
    When given a signal_file and its class_label_file, get the sam result.The delta is set to 1.0
    The command is :
    protocol_engine.py
    --protocol 'differential_expressed_gene_analysis'
    --input 'input_signal_file:<path of the file>'
    --input 'class_label_file:<path of the class_label_file>'
    --parameters 'has_missing_value:zero_fill'
    --parameters 'sam_delta:1.0'

   The result folder will contain a result file for the analysis.
=============================================================================

Heatmap Usage
Betsy can make heatmap for a signal_file without clustering.
    Example:
    When given a signal_file, plot the heatmap and set the heatmap size as x=20,y=20.
    The command is:
    protocol_engine.py
    --protocol 'make_heatmap'
    --input 'input_signal_file:<path of the file>'
    --parameters 'has_missing_value:zero_fill'
    --parameters 'hm_width:200'
    --parameters 'hm_height:10'
    The result folder will contain a png file showing the heatmap.
=============================================================================
Batch effect remove Usage
Example:
When given a signal_file and class_label_file, try five different batch_effect remove methods and plot their pca figure.

     protocol_engine.py
     --protocol 'batch_effect_remove'
     --input 'input_signal_file:<path of the file>'
     --input 'class_label_file:<path of the file>'
     --parameters 'has_missing_value:zero_fill'
     --parameters 'num_factors:1'

