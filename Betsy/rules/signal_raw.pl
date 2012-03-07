%signal_raw.pl

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

:- dynamic gse_dataset/4.
:- dynamic gse_dataset_and_platform/4.
:- dynamic input_signal_file/4.
:- dynamic geo_dataset/4.
:- dynamic class_label_file/4.
:- dynamic gene_list_file/4.


/*-------------------------------------------------------------------------*/
% Input interface
% given a input_signal_file,generate signal_file that will have a full length of parameters
% the parameter which is not provided is in default value

signal_raw(DatasetId,Contents,NewParameters,Modules):-
    input_signal_file(DatasetId,Contents,OldParameters,Modules),
    convert_parameters_raw(OldParameters,NewParameters).
/*-------------------------------------------------------------------------*/
%Output interface
% The parameters in output can be any length and it will trace to the full length one
% the parameter which is not provided will be a variable

signal_raw(DatasetId,Contents,Parameters,Modules):-
    get_desire_parameters_raw(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_raw,N1),
    N<N1,
    convert_parameters_variable_raw(NewParameters1,NewParameters),
    signal_raw(DatasetId,Contents,NewParameters,Modules).


/*-------------------------------------------------------------------------*/
%% geo_dataset(+DatasetId, +Contents,[version,cc],-Modules)
% generate geo_dataset with version cc from cel_file with version unknown.＝
geo_dataset(DatasetId, Contents,[version,cc], Modules):-
    geo_dataset(DatasetId, Contents,[version,unknown_version], Past_Modules),
    Newadd = ['are_any_cc', ['DatasetId',DatasetId,'Version', cc,
             'Contents',Contents]],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
%% geo_dataset(+DatasetId, +Contents,[version,v3_4], -Modules)
% generate geo_dataset with version v3_4 from cel_file with version cc;
% or from cel_file with version unknown.
geo_dataset(DatasetId, Contents, [version,v3_4],  Modules):-
    geo_dataset(DatasetId,  Contents, [version,cc], Past_Modules),
    Newadd=['convert_v3_4', ['DatasetId',
           DatasetId,'Version', 'v3_4','Contents',Contents]],
    append(Past_Modules, Newadd, Modules);
    geo_dataset(DatasetId, Contents,[version,unknown_version], Past_Modules),
    Newadd=['are_all_v3_4',['DatasetId',
           DatasetId,'Version', v3_4,'Contents',Contents]],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
%% geo_dataset(+DatasetId, +Contents,[version,unknown_version],-Modules)
% generate geo_dataset with version unknown from gse_dataset;
% or from gse_dataset_and_platform.
geo_dataset(DatasetId, Contents,[version,unknown_version], Modules):-
    gse_dataset(DatasetId,Contents,[],[]),
    Newadd=['download_geo_dataset',['DatasetId',
            DatasetId,'Contents',Contents]],
    append([], Newadd, Modules);
    gse_dataset_and_platform(DatasetId,Contents,[],[]),
    Newadd=['download_geo_dataset_GPL',['DatasetId',
            DatasetId,
           'Contents',Contents]],
    append([], Newadd, Modules).

/*-------------------------------------------------------------------------*/
%% geo_dataset(+DatasetId,+Contents,[version,gpr],-Modules)
% generate gpr_file from geo_dataset
geo_dataset(DatasetId,Contents,[version,gpr],Modules):-
    geo_dataset(DatasetId,Contents,[version,unknown_version],Past_Modules),
    Newadd=['are_all_gpr', ['DatasetId',
           DatasetId,'Version', 'gpr','Contents',Contents]],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
%% geo_dataset(+DatasetId,+Contents,[version,illumina],-Modules)
% generate gpr_file from geo_dataset
geo_dataset(DatasetId,Contents,[version,illumina],Modules):-
    geo_dataset(DatasetId,Contents,[version,unknown_version],Past_Modules),
    Newadd=['extract_illumina_idat_files', ['DatasetId',
           DatasetId,'Version', 'illumina','Contents',Contents]],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate gpr_file from geo_dataset
geo_dataset(DatasetId,Contents,[version,agilent],Modules):-
    geo_dataset(DatasetId,Contents,[version,unknown_version],Past_Modules),
    Newadd=['extract_agilent_files', ['DatasetId',
           DatasetId,'Version', 'agilent','Contents',Contents]],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
% Preprocess geo_dataset with gpr format to signal_file by loess,
signal_raw(DatasetId, Contents,Parameters,Modules):-
    convert_parameters_raw([preprocess,loess,is_logged,no_logged,format,pcl],NewParameters),
    Parameters=NewParameters,	
    geo_dataset(DatasetId, Contents,[version,gpr],Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[loess,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% Preprocess geo_dataset with illumine format to signal_file by illumina,
signal_raw(DatasetId, Contents,Parameters,Modules):-
    (convert_parameters_raw([preprocess,illumina,is_logged,no_logged,format,gct],NewParameters);
     convert_parameters_raw([preprocess,illumina_controls,is_logged,no_logged,format,gct],NewParameters)),
    Parameters=NewParameters,	
    geo_dataset(DatasetId, Contents,[version,illumina],Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[illumina,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% Preprocess geo_dataset with agilent format to signal_file by agilent,
signal_raw(DatasetId, Contents,Parameters,Modules):-
    convert_parameters_raw([preprocess,agilent,is_logged,no_logged,format,tdf],NewParameters),
    Parameters=NewParameters,	
    geo_dataset(DatasetId, Contents,[version,agilent],Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[preprocess_agilent,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% Preprocess v3_4 cel_file by MAS5 or RMA to signal_file,
signal_raw(DatasetId, Contents, Parameters,Modules):-
    member((Preprocess, Is_Logged),[(mas5, no_logged),(rma, logged)]),
    convert_parameters_raw([preprocess,Preprocess, is_logged, Is_Logged,format,jeffs,has_missing_value,no_missing],NewParameters),
    Parameters=NewParameters,
    geo_dataset(DatasetId, Contents,[version,v3_4], Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[preprocess,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% log the signal file with Is_Logged is no_logged.
signal_raw(DatasetId, Contents,Parameters,Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    OldIs_Logged=no_logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    set_value(Parameters,is_logged,OldIs_Logged,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=['log_algorithm',Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% change the format of signal_file to pcl.
signal_raw(DatasetId, Contents,Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    (member(OldFormat,[tdf,res,gct,jeffs]),
    set_value(Parameters,format,OldFormat,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=['convert_to_pcl',Write_list],
    append(Past_Modules, Newadd, Modules);
     OldFormat=xls,
     set_value(Parameters,format,OldFormat,OldParameters1),
     set_value(OldParameters1,status,OldStatus,OldParameters),
     signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
     append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
     Newadd=['convert_xls_pcl',Write_list],
     append(Past_Modules, Newadd, Modules)).
/*-------------------------------------------------------------------------*/
% check the logged or not for signal_file with Is_logged is unknown.
signal_raw(DatasetId,Contents,Parameters, Modules):-
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    member((Is_Logged, Module),[(logged, check_logged),(no_logged, check_not_logged)]),
    Old_isLogged=unknown_logged,
    set_value(Parameters,is_logged,Old_isLogged,OldParameters),
    signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[Module,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% check the Format of the signal file with Format unknown,
% when Format is not_xls, Module is check_not_xls,
% when Format is xls, Module is check_xls.
signal_raw(DatasetId,Contents,Parameters, Modules):-
    get_value(Parameters,format,unknown_format,Format),
    member((Format, Module),[(not_xls, check_not_xls), (xls, check_xls)]),
    OldFormat=unknown_format,
    set_value(Parameters,format,OldFormat,OldParameters),
    signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[Module,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% check the Format of the signal file with Format not_xls
signal_raw(DatasetId,Contents,Parameters, Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,format,unknown_format,Format),
    member(Format, [res,tdf,jeffs,gct,pcl]),
    OldFormat=not_xls,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,format,OldFormat,OldParameters),
    signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[check_format,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% check missing value of the signal_file with Has_Missing_Value is unknown,
% The signal_file format is pcl,
% Is_logged is yes.
signal_raw(DatasetId,Contents,Parameters, Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    member((Has_Missing_Value,Module),[(yes_missing, check_missing),(no_missing, check_not_missing)]),
    OldHas_Missing_Value=unknown_missing,
    set_value(Parameters,has_missing_value,OldHas_Missing_Value,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[Module,Write_list],
    append(Past_Modules, Newadd, Modules).
    
/*-------------------------------------------------------------------------*/
% filter genes with the missing value
% The signal_file format is pcl.
% filtering occurs before zero filling and gene ordering,
% Is_logged is logged.
signal_raw(DatasetId, Contents,Parameters,Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    Has_Missing_Value=yes_missing,
    get_value(Parameters,fill,no_fill,Fill),
    Fill=no_fill,
    get_value(Parameters,filter,no_filter,Filter),
    not(atom(Filter)),
    OldFilter=no_filter,
    set_value(Parameters,filter,OldFilter,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[gene_filter,Write_list],
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
% zero filling the missing value,
% zero filling occurs before gene ordering,
% the format is pcl and Is_logged is logged.
signal_raw(DatasetId, Contents,Parameters,Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    Has_Missing_Value=yes_missing,
    get_value(Parameters,fill,no_fill,Fill),
    Fill=yes_fill,
    OldFill=no_fill,
    set_value(Parameters,fill,OldFill,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[zero_fill,Write_list],
    append(Past_Modules, Newadd, Modules).

