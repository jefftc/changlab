%signal_file.pl
:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

/*-------------------------------------------------------------------------*/
% Output interface
% The parameters in output can be any length and it will trace to the full length one
% for the parameter which is not provide, will be in no_order

signal_file(DatasetId,Contents,Parameters,Modules):-
    get_desire_parameters_file(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_file,N1),
    N<N1,
    convert_parameters_file(NewParameters1,NewParameters),
    signal_file(DatasetId,Contents,NewParameters,Modules).
/*-------------------------------------------------------------------------*/
% Input interface
% generate signal_file from signal_norm2 with no_order

signal_file(DatasetId,Contents,Parameters,Modules):-
   length(Parameters,N),
    get_length(n_file,N1),
   N=N1,
   get_value(Parameters,gene_order,no_order,Gene_Order),
   Gene_Order=no_order,
   get_desire_parameters_norm2(Parameters,NewParameters),
   signal_norm2(DatasetId,Contents,NewParameters,Modules).

/*--------------------------------------------------------------------------*/
% rank genes by class label
signal_file(DatasetId, Contents,Parameters,Modules):-
    length(Contents,N),
    N is 2,
    class_label_file(DatasetId,Contents,_,Past_Modules_1),
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,gene_order,no_order,Gene_Order),
    Gene_Order=by_sample_ttest,
    OldGene_Order=no_order,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(DatasetId, Contents,OldParameters,Past_Modules_2),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[rank_gene_by_sample_ttest,Write_list],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
% rank genes by gene_list_file
% if input is gene_list_file, gene ordering from signal_file with GeneOrder=no_order; 
% otherwise gene ordering from signal_file with GeneOrder=by_sample_ttest.
% the Is_logged should be yes and missing_value has been checked.

signal_file(DatasetId, Contents,Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,gene_order,no_order,Gene_Order),
    Gene_Order=by_gene_list,
    (gene_list_file(DatasetId,Contents,[],[]),
    OldGene_Order=no_order,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters);
    not(gene_list_file(DatasetId,Contents,[],[])),
    OldGene_Order=by_sample_ttest,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters)),
    signal_file(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[reorder_genes,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% change the pcl format of signal_file to gct.
signal_file(DatasetId, Contents,Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,format,unknown_format,Format),
    Format=gct,
    Pre_format=pcl,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,format,Pre_format,OldParameters),
    signal_file(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=['convert_pcl_gct',Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
pca_plot_out(DatasetId,Contents,Parameters,Modules):-
    signal_file(DatasetId,Contents,Parameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[pca_plot,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
pca_plot_in(DatasetId,Contents,Parameters,Modules):-
    signal_clean(DatasetId,Contents,Parameters,Past_Modules),
    get_desire_parameters_raw(Parameters,NewParameters1),
    convert_parameters_clean_out(NewParameters1,NewParameters),
    append(['DatasetId',
            DatasetId,'Contents',Contents],NewParameters,Write_list),
    Newadd=[pca_plot,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
biotin_plot(DatasetId,Contents,Parameters,Modules):-
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    member(Preprocess,[illumina_controls,illumina]),
    convert_parameters_raw([preprocess,illumina_controls,is_logged,no_logged,
                         format,gct],NewParameters),
    signal_raw(DatasetId,Contents,NewParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],NewParameters,Write_list),
    Newadd=[plot_biotin,Write_list],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
actb_plot(DatasetId,Contents,Parameters,Modules):-
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    Preprocess = illumina,
    signal_clean(DatasetId,Contents,Parameters,Past_Modules),
    get_desire_parameters_raw(Parameters,NewParameters1),
    convert_parameters_clean_out(NewParameters1,NewParameters),
    append(['DatasetId',
            DatasetId,'Contents',Contents],NewParameters,Write_list),
    Newadd=[plot_actb,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
hyb_bar_plot(DatasetId,Contents,Parameters,Modules):-
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    member(Preprocess,[illumina_controls,illumina]),
    convert_parameters_raw([preprocess,illumina_controls,is_logged,no_logged,
                         format,gct],NewParameters),
    signal_raw(DatasetId,Contents,NewParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],NewParameters,Write_list),
    Newadd=[plot_hyb_bar,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
intensity_plot(DatasetId,Contents,Parameters,Modules):-
    signal_file(DatasetId,Contents,Parameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[plot_intensity,Write_list],
    append(Past_Modules, Newadd, Modules).