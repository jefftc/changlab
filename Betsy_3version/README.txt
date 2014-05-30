Preprocess Usage

1)Betsy can preprocess with rma or mas5 for Affymetrix cel data. 

	Examples:
	A. When given a GSE ID from the geo, and want to do rma preprocess, 
	   the command is

	python run_rule.py  \
      --intype 'GEOSeries' \
      --user_input 'GSEID=GSE8286' \
      --outtype 'SignalFile' \
      --attr 'SignalFile,preprocess=rma' \
      --attr 'SignalFile,quantile_norm=yes' \
      --png_file 'out.png'
      ----------------------------------------------------------------
      B. When given a GSE ID and GPL platform number from the geo, 
	   and want to do rma preprocess, 
	   the command is 
     
      python run_rule.py  
      --intype 'GEOSeries' 
      --user_input 'GSEID=GSE17907' 
      --user_input 'GPLID=GPL570'
      --outtype 'SignalFile' 
      --attr 'SignalFile,preprocess=rma' 
      --attr 'SignalFile,quantile_norm=yes'
      --png_file 'out.png'
      ----------------------------------------------------------------
     C. When given a folder contains cel file, the command is

	python run_rule.py  
      --intype 'CELFiles'
      --input '<./GSE21947>'
      --outtype 'SignalFile'     	
      --attr 'SignalFile,preprocess=rma'
      --attr 'SignalFile,quantile_norm=yes'
      --png_file 'out.png'
-----------------------------------------------------------------------
2) Betsy can preprocess with illumina for illumina idat files.
     Example:
	When given a folder contains idat files, the command is
     python run_rule.py 
     --intype 'ExpressionFiles'  
     --input '</home/xchen/chencode/betsy_test/6991010018>' 
     --outtype 'SignalFile'  
     --attr 'SignalFile,preprocess=illumina'  
     --png_file 'out.png'
 ----------------------------------------------------------------
3) Betsy can preprocess with agilent for  Agilent files.
     Example:
	When given a folder contains agilent files, the command is
     python run_rule.py 
     --intype 'ExpressionFiles'  
     --input '</home/xchen/chencode/betsy_test/agilent_expression>' 
     --outtype 'SignalFile'  
     --attr 'SignalFile,preprocess=agilent'  
     --png_file 'out.png'
 ----------------------------------------------------------------
4) Betsy can preprocess with gpr for gpr files.
     Example:
	When given a folder contains gpr files, the command is
     python run_rule.py 
     --intype 'ExpressionFiles'  
     --input '</home/xchen/chencode/betsy_test/GSE4189>' 
     --outtype 'SignalFile'  
     --attr 'SignalFile,preprocess=loess'  
     --png_file 'out.png'
=============================================================================
Process Usage
Betsy can do predataset,log,unlog,gene_filter,quantile,combat, shiftscale,dwd,bfrm, predataset,gene_center,gene_normalize,gene_order, annotate, rename_sample, platform,
num_features,unique_genes, duplicate_probe,group_fc,change format,for signal files

The option of the attributes are:
    preprocess:   unknown, illumina, agilent, mas5, rma, loess
    missing_algorithm: none, median_fill, zero_fill             filter:no, yes     dwd_norm: no, yes    bfrm_norm: no, yes    quantile_norm: no, yes    shiftscale_norm: no, yes    combat_norm: no", yes    predataset: no, yes    gene_center: no, mean, median    gene_normalize:  no, variance, sum_of_squares    gene_order:no, class_neighbors, gene_list, t_test_p, t_test_fdr    annotate: no,yes    rename_sample: no, yes    platform: yes,no    num_features: yes,no    unique_genes: no, average_genes, high_var, first_gene,    duplicate_probe:no, closest_probe, high_var_probe        group_fc: yes,no    contents:unspecified, train0, train1, test, class0,class1,test,        class0, class1, class0,class1,          logged: no, yes    format: tdf, gct

 ----------------------------------------------------------------
    Example:
    When given a SignalFile,do predataset, quantile, gene_center=mean,gene_normalize=variance.
    The command is
    python run_rule.py
    --intype 'SignalFile_Postprocess'
    --input '</home/xchen/chencode/betsy_test/all_aml_train_missed.gct>'
    --outtype 'SignalFile'
    --attr 'SignalFile,predataset=yes'
    --attr 'SignalFile,quantile_norm=yes'
    --attr 'SignalFile,gene_center=mean'
    --attr 'SignalFile,gene_normalize=variance'
 ----------------------------------------------------------------
    When given a SignalFile and ClassLabelFile, group_fc=yes and group_fc_num=1
    The command is:
	python run_rule.py 
	--intype 'SignalFile_Postprocess'  
	--input '/home/xchen/chencode/betsy_test/all_aml_train_filt.res' 
	--intype 'ClassLabelFile' 
	--attr 'cls_format=cls' 
	--input '/home/xchen/chencode/betsy_test/all_aml_train.cls' 
	--outtype 'SignalFile'  
	--attr 'SignalFile,group_fc=yes'  
	--user_input 'group_fc_num=1' 
	--png_file 'out.png'
=============================================================================
Heatmap Usage
Betsy can make heatmap for a SignalFile without clustering.
    Example:
    When given a signal_file, plot the heatmap and set the heatmap size as x=20,y=20.
    The command is:
	python run_rule.py \
	--intype 'SignalFile_Postprocess' \
	--input '/home/xchen/chencode/betsy_test/breast_19.mas5'  \
	--outtype 'ReportFile'   \
	--attr 'ReportFile,report_type=heatmap' \
     --user_input 'hm_width=20' \
     --user_input 'hm_height=20' \
     --png_file 'out.png'

   
    The result folder will contain a png file showing the heatmap.
=============================================================================
Clustering Usage
Betsy can do clustering for a SignalFile and plot the heatmap.

    Example:
    When given a signalFile, gene_center=mean, gene_normalize=variance, 
    cluster the genes by hierarchical method and plot the heatmap, 
    set the heatmap size as x=200,y=1.
    The command is:
	python run_rule.py \
	--intype 'SignalFile_Postprocess'  \
	--input '/home/xchen/chencode/betsy_test/breast_19.mas5'  \
	--outtype 'ReportFile'   \
	--attr 'ReportFile,report_type=cluster' \
     --attr 'SignalFile,gene_normalize=variance' \
     --attr 'SignalFile,gene_center=mean' \
	--attr 'ClusterFile,cluster_alg=pca' \
	--attr 'ClusterFile,distance=correlation' \
     --user_input 'hm_width=200' \
     --user_input 'hm_height=1' \
	--png_file 'out.png'

The result folder will contain a clustering file and a png file showing the heatmap.
===============================================================================
Classification Usage
Betsy can do classification with svm and weighted_voting method for dataset. 
Also it can leave one out cross validation by these two methods.
The command is:
	python run_rule.py \
	--intype 'SignalFile_Postprocess' \
	--attr 'contents=class0,class1' \
	--input '/home/xchen/chencode/betsy_test/all_aml_train.res' \
	--intype 'SignalFile_Postprocess' \
	--attr 'contents=test' \
	--input '/home/xchen/chencode/betsy_test/all_aml_test.res' \
	--intype 'ClassLabelFile' \
	--attr 'contents=class0,class1' \
	--attr 'cls_format=cls' \
	--input '/home/xchen/chencode/betsy_test/all_aml_train.cls' \
	--intype 'ClassLabelFile' \
	--attr 'contents=test' \
	--attr 'cls_format=cls' \
	--input '/home/xchen/chencode/betsy_test/all_aml_test.cls' \
	--outtype 'ReportFile' \
	--attr 'ReportFile,report_type=classify'  \
	--png_file 'out.png' \
	--text_file 'out.txt' 

===============================================================================
Differential expressed genes analysis usage
Betsy can do the differential expressed genes analysis for signal_files.
The command is:
	python run_rule.py \
	--intype 'SignalFile_Postprocess'  \
	--input '/home/xchen/chencode/betsy_test/all_aml_train_filt.res' \
	--intype 'ClassLabelFile' \
	--attr 'cls_format=cls' \
	--input '/home/xchen/chencode/betsy_test/all_aml_train.cls' \
	--outtype 'ReportFile'  \
	--attr 'ReportFile,report_type=diffgenes' \
	--png_file 'out.png'
===============================================================================
Geneset Analysis Usage
Example:
When given a signal_file and a gene set file, try do geneset score analysis and plot the result.

	python run_rule.py \
	--outtype 'ReportFile' \
	--attr 'ReportFile,report_type=geneset' \
	--attr 'SignalFile,quantile_norm=yes' \
	--attr 'SignalFile,gene_center=mean' \
      --attr 'SignalFile,gene_normalize=variance' \
	--attr 'SignalFile,unique_genes=high_var' \
	--attr 'SignalFile,annotate=yes'  \
	--intype 'SignalFile_Postprocess'  \
	--input '/home/xchen/chencode/betsy_test/se2fplate6_48.illu.gz' \
	--intype 'GenesetFile' \
	--input '/home/xchen/chencode/betsy_test/genesets.gmt' \
	--user_input 'geneset_value=E2F1n_affy_150_UP' \
	--png_file 'out.png'
===============================================================================
Normalization Usage
Example:
When given a signal_file, try do quantile_norm and gene_center=median,get a normalization report

	python run_rule.py \
	--intype 'SignalFile_Postprocess'  \
	--input '/home/xchen/chencode/betsy_test/all_aml_train_filt.res' \
	--intype 'ClassLabelFile' \
	--attr 'cls_format=cls' \
	--input '/home/xchen/chencode/betsy_test/all_aml_train.cls' \
	--outtype 'ReportFile'  \
	--attr 'ReportFile,report_type=normalize_file' \
	--attr 'SignalFile,quantile_norm=yes'  \
	--attr 'SignalFile,preprocess=unknown' \
	--attr 'SignalFile,gene_center=median' \
	--png_file 'out.png' \
      --json_file 'outjson.txt'
--------------------------------------------------
When given a ExpressionFiles, try rma preprocess, do quantile_norm and gene_center=median,get a normalization report

	python run_rule.py \
	--intype 'ExpressionFiles'  \
	--input '/home/xchen/chencode/betsy_test/GSE8286_folder' \
	--outtype 'ReportFile'  \
	--attr 'ReportFile,report_type=normalize_file' \
	--attr 'SignalFile,quantile_norm=yes'  \
	--attr 'SignalFile,preprocess=rma' \
	--attr 'SignalFile,gene_center=median' \
	--png_file 'out.png'
--------------------------------------------------
When given a ExpressionFiles, try illumina preprocess, do quantile_norm and gene_center=median,get a normalization report

	python run_rule.py \
	--intype 'ExpressionFiles'  \
	--input '/home/xchen/chencode/betsy_test/6991010018' \
	--outtype 'ReportFile'  \
	--attr 'ReportFile,report_type=normalize_file' \
	--attr 'SignalFile,quantile_norm=yes'  \
	--attr 'SignalFile,preprocess=illumina' \
	--attr 'SignalFile,gene_center=median' \
	--png_file 'out.png'