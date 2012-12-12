:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

/*-------------------------------------------------------------------------*/
% Output interface
% The parameters in output can be any length and it will trace to the full length one
% for the parameter which is not provide, will be in no_order

signal_file(Parameters,Modules):-
    % Conditions: the parameters is not full length
    get_desire_parameters_file(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_file,N1),
    N<N1,
    % Input: signal_file with full length parameters
    convert_parameters_file(Parameters,NewParameters),
    signal_file(NewParameters,Modules).
/*-------------------------------------------------------------------------*/
% Input interface
% generate signal_file from signal_norm2 with no_order,unknown_platform,no_unique_genes,yes_missing_probe,yes_dupliate_probe

signal_file(Parameters,Modules):-
   % Conditions: full length and no_order,unknown_platform,no_unique_genes
   get_desire_parameters_file(Parameters,NewParameters1),
   length(NewParameters1,N),
   get_length(n_file,N1),
   N=N1,
   get_value(Parameters,gene_order,no_order,Gene_Order),
   Gene_Order=no_order,
   get_value(Parameters,num_features,0,Num_features),
   Num_features=0,
   get_value(Parameters,platform,unknown_platform,Platform),
   Platform=unknown_platform,
   get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
   Unique_Genes=no_unique_genes, 
   get_value(Parameters,duplicate_probe,yes_duplicate_probe,Duplicate_Probe),
   Duplicate_Probe=yes_duplicate_probe,
   get_value(Parameters,duplicate_data,yes_duplicate_data,Duplicate_Data),
   Duplicate_Data=yes_duplicate_data,
   get_value(Parameters,has_annotation_gene_id,no_gene_id,Has_annotation_gene_id),
   Has_annotation_gene_id=no_gene_id,
    
   /*get_value(Parameters,has_annotations,no_annotations,Has_annotations),
   Has_annotations=no_annotations,*/
   get_value(Parameters,group_fc,no_group_fc,Group_fc),
   Group_fc=no_group_fc,
   % Input: signal_norm2
   get_desire_parameters_norm2(Parameters,NewParameters),
   get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip,num_factors,normalize_order],[],Options),
   append(NewParameters,Options,NewParameters2),
   signal_norm2(NewParameters2,Modules).

/*--------------------------------------------------------------------------*/
% rank genes by sample_ttest
signal_file(Parameters,Modules):-
    % Conditions:Parameters has created,tdf,gene_order in [t_test_p,t_test_fdr]
    %   unknown_platform,no_unique_genes and logged
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format = tdf,
    get_value(Parameters,is_logged,unknown_logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,platform,unknown_platform,Platform),
    Platform=unknown_platform,
    get_value(Parameters,gene_order,no_order,Gene_Order),
    member(Gene_Order,[t_test_p,t_test_fdr]),
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    Unique_Genes = no_unique_genes,
    get_value(Parameters,num_features,0,Num_features),
    Num_features=0,
    % Input: class_label_file and signal_file with no_order,tdf,different status
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    OldGene_Order=no_order,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(OldParameters,Past_Modules_2),
    % Module:ranek_genes_by_sample_ttest,reorder_genes
    % Output Parameters: full length parameters of signal_file
    Newadd=[rank_genes_by_sample_ttest,Parameters,reorder_genes,Parameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
%rank genes by class_neighbors
signal_file(Parameters,Modules):-
    % Conditions: Parameters has created,tdf,by_class_neighbors
    %   unknown_platform,no_unique_genes and logged
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,gene_order,no_order,Gene_Order),
    Gene_Order=by_class_neighbors,
    get_value(Parameters,platform,unknown_platform,Platform),
    Platform=unknown_platform,
    get_value(Parameters,is_logged,unknown_logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    Unique_Genes = no_unique_genes,
    get_value(Parameters,num_features,0,Num_features),
    Num_features=0,
    % Input: class_label_file and signal_file with tdf,no_order and different status
    %   unknown_platform,no_unique_genes and logged
    % Input: class_label_file and signal_file with tdf,no_order and different status
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    OldGene_Order=no_order,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(OldParameters,Past_Modules_2),
    % Module: rank_genes_by_class_neighbors,reorder_genes
    % Output Parameters: full length parameters of signal_file
    Newadd=[rank_genes_by_class_neighbors,Parameters,reorder_genes,Parameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
% rank genes by gene_list_file
% if input is gene_list_file, gene ordering from signal_file with GeneOrder=no_order; 

signal_file(Parameters,Modules):-
    % Conditions:Parameters has created,by_gene_list,tdf
    %   unknown_platform and logged,num_features=0
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,is_logged,unknown_logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,platform,unknown_platform,Platform),
    Platform=unknown_platform,
    get_value(Parameters,num_features,0,Num_features),
    Num_features=0,
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    Unique_Genes = no_unique_genes,
    get_value(Parameters,gene_order,no_order,Gene_Order),
    Gene_Order=by_gene_list,
    % Input: gene_list_file and signal_file with created,pcl,no_order
    %   unknown_platform,no_unique_genes and logged
    get_value(Parameters,contents,[unknown],Contents),
    gene_list_file([contents,Contents],[]),
    OldGene_Order=no_order,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:reorder_genes
    % Output Parameters: full length parameters of signal_file
    Newadd=[reorder_genes,Parameters],
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
% select genes according to num_features
signal_file(Parameters,Modules):-
    % Conditions:Parameters has created,tdf,num_features>0
    %   unknown_platform,no_unique_genes and logged
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format = tdf,
    get_value(Parameters,is_logged,unknown_logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,num_features,0,Num_features),
    Num_features>0,
    get_value(Parameters,platform,unknown_platform,Platform),
    Platform=unknown_platform,
    % Input: signal_file with num_features=0,tdf,different status
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,num_features,0,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:select_first_n_genes
    (member(traincontents,Parameters),
    member(testcontents,Parameters),
    get_value(Parameters,traincontents,[unknown],TrainContents),
    get_value(Parameters,testcontents,[unknown],TestContents),
    convert_parameters_file(Parameters,NewParameters),
    set_value(NewParameters,contents,TrainContents,NewParameters1),
    set_value(NewParameters,contents,TestContents,NewParameters2),
    % Output Parameters: full length parameters of signal_file
    Newadd=[select_first_n_genes,NewParameters1,select_first_n_genes,NewParameters2];
    not(member(traincontents,Parameters)),
    Newadd=[select_first_n_genes,Parameters]),
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
%annotate_gene_id
signal_file(Parameters,Modules):-
    % Conditions: Parameters has created,tdf,unknown_platform and
    %      no_unique_genes
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,is_logged,logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,platform,unknown_platform,Platform),
    Platform=unknown_platform,
    get_value(Parameters,num_features,0,Num_features),
    Num_features=0,
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    Unique_Genes=no_unique_genes,
    get_value(Parameters,has_annotation_gene_id,no_gene_id,Has_annotation_gene_id),
    Has_annotation_gene_id=yes_gene_id,
    % Input: signal_file with no_unique_genes,tdf,unknown_platform
    Pre_has_annotation_gene_id=no_gene_id,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,has_annotation_gene_id,Pre_has_annotation_gene_id,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:annotate_probes
    % Output Parameters: full length parameters of signal_file 
    append(Parameters,[filetype,signal_file,annotate_type,gene_id],NewParameters),
    Newadd=[annotate_probes,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
%get the unique genes
signal_file(Parameters,Modules):-
    % Conditions: Parameters has created,tdf,unknown_platform and
    %     unique_genes in [average_genes,high_var,first_gene]
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,is_logged,logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,platform,unknown_platform,Platform),
    Platform=unknown_platform,
    get_value(Parameters,num_features,0,Num_features),
    Num_features=0,
    get_value(Parameters,has_annotation_gene_id,no_gene_id,Has_annotation_gene_id),
    Has_annotation_gene_id=yes_gene_id,
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    member(Unique_Genes,[average_genes,high_var,first_gene]),
    % Input: signal_file with no_unique_genes,tdf,unknown_platform
    Pre_unique_genes=no_unique_genes,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,unique_genes,Pre_unique_genes,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:remove_duplicate_genes
    % Output Parameters: full length parameters of signal_file 
    Newadd=[remove_duplicate_genes,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
%convert the unknown_platform to desired platform
signal_file(Parameters,Modules):-
    % Conditions: Parameters has created,tdf,logged,not unknown_platform
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,is_logged,logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,duplicate_probe,yes_duplicate_probe,Duplicate_Probe),
    Duplicate_Probe = yes_duplicate_probe,
    get_value(Parameters,platform,unknown_platform,Platform),
    not(Platform=unknown_platform),
    % Input: signal_file with unknown_platform,tdf,logged,different status
    Pre_platform=unknown_platform,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,platform,Pre_platform,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:add_crossplatform_probeid
    % Output Parameters: full length parameters of signal_file
    Newadd=[add_crossplatform_probeid,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% select_unique_probe
signal_file(Parameters,Modules):-
    % Conditions: Parameters has created,tdf,not unknwon_platform and
    % duplicate_probe in  [high_var_probe,closest_probe]
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,is_logged,logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,duplicate_probe,yes_duplicate_probe,Duplicate_Probe),
    (Duplicate_Probe=high_var_probe,
     Module=remove_duplicate_probes;
     Duplicate_Probe=closest_probe,
     Module=select_probe_by_best_match),
    % Input: signal_file with no_unique_genes,tdf,not unknwon_platform,yes_duplicate_probe
    Pre_duplicate_probe=yes_duplicate_probe,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,duplicate_probe,Pre_duplicate_probe,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:remove_duplicate_probes or select_probe_by_best_match
    % Output Parameters: full length parameters of signal_file 
    Newadd=[Module,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% remove_duplicate_data
/*signal_file(Parameters,Modules):-
    % Conditions: Parameters has created,tdf,not unknown_platform and no_duplicate_data
    % no_missing_probe,duplicate_probe in  [high_var_probe,closest_probe]
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,is_logged,logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,platform,unknown_platform,Platform),
    not(Platform=unknown_platform),
    get_value(Parameters,duplicate_probe,yes_duplicate_probe,Duplicate_Probe),
    Duplicate_Probe=high_var_probe,
    get_value(Parameters,duplicate_data,yes_duplicate_data,Duplicate_Data),
    Duplicate_Data=no_duplicate_data,
    % Input: signal_file with no_unique_genes,tdf,not unknown_platform,yes_duplicate_data
    Pre_duplicate_data=yes_duplicate_data,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,duplicate_data,Pre_duplicate_data,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:remove_duplicate_data,
    % Output Parameters: full length parameters of signal_file 
    Newadd=[remove_duplicate_data,Parameters],
    append(Past_Modules, Newadd, Modules).*/

/*-------------------------------------------------------------------------*/
% change tdf format of signal_file to gct.
signal_file(Parameters,Modules):-
    % Conditions: Parameters has created,gct
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=gct,
    % Input: signal_file with tdf and different status
    Pre_format=tdf,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,format,Pre_format,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:convert_signal_to_gct
    % Output Parameters: full length parameters of signal_file
    Newadd=[convert_signal_to_gct,Parameters],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
% unlog the signal file.
signal_file(Parameters,Modules):-
    % Conditions: Parameters has created,pcl,no_logged
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,is_logged,logged,Is_logged),
    Is_logged=no_logged,
    % Input: signal_file with tdf,logged and different status
    Pre_islogged=logged,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,is_logged,Pre_islogged,OldParameters),
    signal_file(OldParameters,Past_Modules),
    % Module:unlog_signal
    % Output Parameters: full length parameters of signal_file
    (member(traincontents,Parameters),
    member(testcontents,Parameters),
    get_value(Parameters,traincontents,[unknown],TrainContents),
    get_value(Parameters,testcontents,[unknown],TestContents),
    convert_parameters_file(Parameters,NewParameters),
    set_value(NewParameters,contents,TrainContents,NewParameters1),
    set_value(NewParameters,contents,TestContents,NewParameters2),
    % Output Parameters: full length parameters of signal_file
    Newadd=[unlog_signal,NewParameters1,unlog_signal,NewParameters2];
    not(member(traincontents,Parameters)),
    Newadd=[unlog_signal,Parameters]),
    append(Past_Modules, Newadd, Modules).
/*---------------------------------------------------------------*/
%algin train signal file and test signal file
signal_file(Parameters,Modules):-
    % Conditions: Parameters has traincontents,testcontents,pcl,
    %     created,logged,unknown_platform
    % Conditions: Parameters has traincontents,testcontents,pcl,created
    member(traincontents,Parameters),
    member(testcontents,Parameters),
    get_value(Parameters,traincontents,[unknown],TrainContents),
    get_value(Parameters,testcontents,[unknown],TestContents),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=tdf,
    get_value(NewParameters,status,created,Status),
    Status=created,
    get_value(Parameters,num_features,0,Num_features),
    Num_features=0,
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    Unique_Genes=no_unique_genes,
    %get_value(NewParameters,is_logged,unknown_logged,Is_logged),
    %Is_logged=logged,
    % Input: signal_file with traincontents,tdf,different status,
    %     logged,unknown_platform and  signal_file with testcontents,pcl,
    %     different status,logged,unknown_platform
    % Input: signal_file with traincontents,pcl,different status and 
    %   signal_file with testcontents,pcl,different status
    member(OldStatus1,[given,jointed,splited,created]),
    member(OldStatus2,[given,jointed,splited,created]),
    set_value(NewParameters,status,OldStatus1,Parameters1),
    set_value(NewParameters,status,OldStatus2,Parameters2),
    set_value(Parameters1,contents,TestContents,NewParameters1),
    signal_file(NewParameters1,Past_Modules_1),
    set_value(Parameters2,contents,TrainContents,NewParameters2),
    signal_file(NewParameters2,Past_Modules_2),
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append([testcontents,TestContents,
         traincontents,TrainContents],NewParameters1,Write_list1),
    append([testcontents,TestContents,
         traincontents,TrainContents],NewParameters2,Write_list2),
    % Module:select_common_genes
    % Output parameters: full length parameters of signal_file and 
    %     [testcontents,TestContents,traincontents,TrainContents]
    Newadd1=[select_common_genes,Write_list2],
    Newadd2=[select_common_genes,Write_list1],
    append(Newadd1,Newadd2,Newadd),
    append(Past_Modules,Newadd,Modules).
/*-------------------------------------------------------------------------*/
/*pca_plot_out(Parameters,Modules):-
    % Input: signal_file with tdf,unknown_platform and class_label_file
    convert_parameters_file(Parameters,NewParameters1),
    set_value(NewParameters1,format,tdf,NewParameters2),
    signal_file(NewParameters2,Past_Modules_2),
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters2,Options,NewParameters3),
    append(NewParameters3,[objecttype,pca_plot_out],NewParameters),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    % Module:plot_pca
    % Output parameters: full length parameters of signal_file 
    Newadd=[plot_pca,NewParameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).*/
/*-------------------------------------------------------------------------*/
/*pca_plot_in(Parameters,Modules):-
    % Input: signal_clean and class_label_file
    convert_parameters_clean_out(Parameters,NewParameters2),
    signal_clean(NewParameters2,Past_Modules_2),
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters2,Options,NewParameters3),
    append(NewParameters3,[objecttype,pca_plot_in],NewParameters),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    % Module:plot_pca
    % Output parameters:full length parameters of signal_clean
    Newadd=[plot_pca,NewParameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).*/
/*-------------------------------------------------------------------------*/
biotin_plot(Parameters,Modules):-
    % Input:control_file with illumina,no_logged,gct
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    get_value(Parameters,contents,[unknown],Contents),
    member(Preprocess,[illumina]),
    convert_parameters_raw([contents,Contents,preprocess,illumina,is_logged,no_logged,
                         format,gct],NewParameters),
    control_file(NewParameters,Past_Modules),
    % Module:plot_illu_biotin_line
    % Output parameters:full length parameters of signal_raw
    Newadd=[plot_illu_biotin_line,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
housekeeping_plot(Parameters,Modules):-
    % Input:control_file with illumina,no_logged,gct
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    get_value(Parameters,contents,[unknown],Contents),
    member(Preprocess,[illumina]),
    convert_parameters_raw([contents,Contents,preprocess,illumina,is_logged,no_logged,
                         format,gct],NewParameters),
    control_file(NewParameters,Past_Modules),
    % Module:plot_illu_housekeeping_line
    % Output parameters:full length parameters of signal_raw
    Newadd=[plot_illu_housekeeping_line,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
control_plot(Parameters,Modules):-
    % Input:control_file with not_illumina
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    not(member(Preprocess,[illumina])),
    convert_parameters_file(Parameters,NewParameters),
    signal_file(NewParameters,Past_Modules),
    % Module:plot_affy_affx_line
    % Output parameters:full length parameters of signal_file
    Newadd=[plot_affy_affx_line,NewParameters],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
actb_plot(Parameters,Modules):-
    % Conditions: Parameters has illumina
    %get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    %Preprocess = illumina,
    % Input: signal_clean 
    convert_parameters_clean_out(Parameters,NewParameters),
    signal_clean(NewParameters,Past_Modules),
    % Module:plot_actb_line
    % Output parameters:full length parameters of signal_clean
    Newadd=[plot_actb_line,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
hyb_bar_plot(Parameters,Modules):-
    % Conditions: Parameters has illumina
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    get_value(Parameters,contents,[unknown],Contents),
    Preprocess=illumina,
    % Input:control_file with illumina,no_logged and gct
    convert_parameters_raw([contents,Contents,preprocess,illumina,is_logged,no_logged,
                         format,gct],NewParameters),
    control_file(NewParameters,Past_Modules),
    % Module:plot_illu_hyb_bar
    % Output parameters:full length parameters of signal_raw
    Newadd=[plot_illu_hyb_bar,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
intensity_plot(Parameters,Modules):-
    % Input:signal_file
    convert_parameters_file(Parameters,NewParameters),
    signal_file(NewParameters,Past_Modules),
    % Module:plot_intensity_boxplot
    % Output parameters:full length parameters of signal_file
    Newadd=[plot_intensity_boxplot,NewParameters],
    append(Past_Modules, Newadd, Modules).


/*--------------------------------------------------------------------------*/
% filter genes by fold change across classes
signal_file(Parameters,Modules):-
    % Conditions:Parameters has created,tdf,gene_order in [t_test_p,t_test_fdr]
    %   unknown_platform,no_unique_genes and logged
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format = tdf,
    get_value(Parameters,is_logged,unknown_logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,platform,unknown_platform,Platform),
    Platform=unknown_platform,
    get_value(Parameters,group_fc,no_group_fc,Group_fc),
    Group_fc=yes_group_fc,
    get_value(Parameters,gene_order,no_order,Gene_Order),
    Gene_Order=no_order,
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    Unique_Genes = no_unique_genes,
    get_value(Parameters,num_features,0,Num_features),
    Num_features=0,
    % Input: class_label_file and signal_file with no_order,tdf,different status
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    OldGroup_fc=no_group_fc,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,group_fc,OldGroup_fc,OldParameters),
    signal_file(OldParameters,Past_Modules_2),
    % Module:filter_genes_by_fold_change_across_classes
    % Output Parameters: full length parameters of signal_file
    Newadd=[filter_genes_by_fold_change_across_classes,Parameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
