:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

/*-------------------------------------------------------------------------*/
% Output interface
% The parameters in output can be any length and it will trace to the full length one
% for the parameter which is not provide, will be in no_order

signal_file(Parameters,Modules):-
    get_desire_parameters_file(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_file,N1),
    N<N1,
    convert_parameters_file(Parameters,NewParameters),
    signal_file(NewParameters,Modules).
/*-------------------------------------------------------------------------*/
% Input interface
% generate signal_file from signal_norm2 with no_order

signal_file(Parameters,Modules):-
   get_desire_parameters_file(Parameters,NewParameters1),
   length(NewParameters1,N),
   get_length(n_file,N1),
   N=N1,
   get_value(Parameters,gene_order,no_order,Gene_Order),
   Gene_Order=no_order,
   get_value(Parameters,platform,unknown_platform,Platform),
   Platform=unknown_platform,
   get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
   Unique_Genes=no_unique_genes,
   get_desire_parameters_norm2(Parameters,NewParameters),
   get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip,num_factors],[],Options),
   append(NewParameters,Options,NewParameters2),
   signal_norm2(NewParameters2,Modules).

/*--------------------------------------------------------------------------*/
% rank genes by sample_ttest
signal_file(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,gene_order,no_order,Gene_Order),
    member(Gene_Order,[t_test_p,t_test_fdr]),
    OldGene_Order=no_order,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(OldParameters,Past_Modules_2),
    Newadd=[rank_gene_by_sample_ttest,Parameters,reorder_genes,Parameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
%rank genes by class_neighbors
signal_file(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,gene_order,no_order,Gene_Order),
    Gene_Order=by_class_neighbors,
    OldGene_Order=no_order,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(OldParameters,Past_Modules_2),
    Newadd=[class_neighbors,Parameters,reorder_genes,Parameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
% rank genes by gene_list_file
% if input is gene_list_file, gene ordering from signal_file with GeneOrder=no_order; 
% the Is_logged should be yes and missing_value has been checked.

signal_file(Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,gene_order,no_order,Gene_Order),
    Gene_Order=by_gene_list,
    get_value(Parameters,contents,[unknown],Contents),
    gene_list_file([contents,Contents],[]),
    OldGene_Order=no_order,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(OldParameters,Past_Modules),
    Newadd=[reorder_genes,Parameters],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
% change the pcl format of signal_file to gct.
signal_file(Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,format,unknown_format,Format),
    Format=gct,
    Pre_format=pcl,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,format,Pre_format,OldParameters),
    signal_file(OldParameters,Past_Modules),
    Newadd=[convert_pcl_gct,Parameters],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
% unlog the signal file.
signal_file(Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,is_logged,logged,Is_logged),
    Is_logged=no_logged,
    Pre_islogged=logged,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,is_logged,Pre_islogged,OldParameters),
    signal_file(OldParameters,Past_Modules),
    Newadd=[unlog_signal_file,Parameters],
    append(Past_Modules, Newadd, Modules).
/*---------------------------------------------------------------*/
signal_file(Parameters,Modules):-
    member(traincontents,Parameters),
    member(testcontents,Parameters),
    get_value(Parameters,traincontents,[unknown],TrainContents),
    get_value(Parameters,testcontents,[unknown],TestContents),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=pcl,
    get_value(NewParameters,status,created,Status),
    Status=created,
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
    Newadd1=[select_common_genes,Write_list1],
    Newadd2=[select_common_genes,Write_list2],
    append(Newadd1,Newadd2,Newadd),
    append(Past_Modules,Newadd,Modules).
/*-------------------------------------------------------------------------*/
%convert the unknown_platform to HG_U133A platform
signal_file(Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,is_logged,logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,platform,unknown_platform,Platform),
    Platform="HG_U133A",
    Pre_platform=unknown_platform,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,platform,Pre_platform,OldParameters),
    signal_file(OldParameters,Past_Modules),
    Newadd=[convert_platform_to_HG_U133A_if_not,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
%get the unique genes
signal_file(Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    member(Unique_Genes,[average_genes,high_var,first_gene]),
    Pre_unique_genes=no_unique_genes,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,unique_genes,Pre_unique_genes,OldParameters),
    signal_file(OldParameters,Past_Modules),
    Newadd=[get_unique_genes,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
pca_plot_out(Parameters,Modules):-
    convert_parameters_file(Parameters,NewParameters1),
    set_value(NewParameters1,format,pcl,NewParameters2),
    set_value(NewParameters1,platform,unknown_platform,NewParameters4),
    signal_file(NewParameters4,Past_Modules_2),
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters1,Options,NewParameters3),
   append(NewParameters3,[objecttype,pca_plot_out],NewParameters),
    Newadd=[pca_plot,NewParameters],
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
pca_plot_in(Parameters,Modules):-
    convert_parameters_clean_out(Parameters,NewParameters2),
    signal_clean(NewParameters2,Past_Modules_2),
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters2,Options,NewParameters3),
    append(NewParameters3,[objecttype,pca_plot_in],NewParameters),
    Newadd=[pca_plot,NewParameters],
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
biotin_plot(Parameters,Modules):-
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    get_value(Parameters,contents,[unknown],Contents),
    member(Preprocess,[illumina]),
    convert_parameters_raw([contents,Contents,preprocess,illumina,is_logged,no_logged,
                         format,gct],NewParameters),
    control_file(NewParameters,Past_Modules),
    Newadd=[plot_biotin,NewParameters],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
actb_plot(Parameters,Modules):-
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    Preprocess = illumina,
    convert_parameters_clean_out(Parameters,NewParameters),
    signal_clean(NewParameters,Past_Modules),
    Newadd=[plot_actb,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
hyb_bar_plot(Parameters,Modules):-
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    get_value(Parameters,contents,[unknown],Contents),
    member(Preprocess,[illumina]),
    convert_parameters_raw([contents,Contents,preprocess,illumina,is_logged,no_logged,
                         format,gct],NewParameters),
    control_file(NewParameters,Past_Modules),
    Newadd=[plot_hyb_bar,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
intensity_plot(Parameters,Modules):-
    convert_parameters_file(Parameters,NewParameters),
    signal_file(NewParameters,Past_Modules),
    Newadd=[plot_intensity,NewParameters],
    append(Past_Modules, Newadd, Modules).

