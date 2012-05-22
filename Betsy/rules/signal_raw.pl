%signal_raw.pl

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

:- dynamic gse_id/2.
:- dynamic gse_id_and_platform/2.
:- dynamic input_signal_file/2.
:- dynamic gpr_files/2.
:- dynamic cel_files/2.
:- dynamic idat_files/2.
:- dynamic agilent_files/2.
:- dynamic class_label_file/2.
:- dynamic gene_list_file/2.
:- dynamic rename_list_file/2.

/*-------------------------------------------------------------------------*/
% Input interface
% given a input_signal_file,generate signal_file that will have a full length of parameters
% the parameter which is not provided is in default value

signal_raw(NewParameters,Modules):-
    input_signal_file(OldParameters,Modules),
    convert_parameters_raw(OldParameters,NewParameters).
/*-------------------------------------------------------------------------*/
gpr_files([contents,Contents,version,Version],Modules):-
   gpr_files([version,Version,contents,Contents],Modules).

cel_files([contents,Contents,version,Version],Modules):-
   cel_files([version,Version,contents,Contents],Modules).

idat_files([contents,Contents,version,Version],Modules):-
   idat_files([version,Version,contents,Contents],Modules).

agilent_files([contents,Contents,version,Version],Modules):-
   agilent_files([version,Version,contents,Contents],Modules).
/*-------------------------------------------------------------------------*/
%Output interface
% The parameters in output can be any length and it will trace to the full length one

signal_raw(Parameters,Modules):-
    get_desire_parameters_raw(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_raw,N1),
    N<N1,
    convert_parameters_variable_raw(Parameters,NewParameters),
    signal_raw(NewParameters,Modules).

/*-------------------------------------------------------------------------*/
% generate geo_dataset with version v3_4 from cel_file with version cc;
% or from cel_file with version unknown.
cel_files([contents,Contents,version,v3_4], Modules):-
    member(Oldversion,[cc,unknown_version]),
    cel_files([contents,Contents,version,Oldversion], Past_Modules),
    Newadd=[convert_v3_4_if_not_v3_4, [contents,Contents,version,v3_4]],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate gpr_files with version unknown from gse_id;
% or from gse_id_and_platform
gpr_files([contents,Contents,version,unknown_version], Modules):-
    (
    gse_id([contents,Contents],[]),
    Newadd=['download_geo_dataset',[contents,Contents,version,unknown_version,filetype,gpr_files]];
    gse_id_and_platform([contents,Contents],[]),
    Newadd=['download_geo_dataset_GPL',[contents,Contents,version,unknown_version,filetype,gpr_files]]
     ),
    append( [],Newadd,Modules).
/*-------------------------------------------------------------------------*/
% generate cel_files with version unknown from gse_id;
% or from gse_id_and_platform
cel_files([contents,Contents,version,unknown_version], Modules):-
    (
    gse_id([contents,Contents],[]),
    Newadd=['download_geo_dataset',[contents,Contents,version,unknown_version,filetype,cel_files]];
    gse_id_and_platform([contents,Contents],[]),
    Newadd=['download_geo_dataset_GPL',[contents,Contents,version,unknown_version,filetype,cel_files]]
     ),
    append( [],Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate agilent_files with version unknown from gse_id;
% or from gse_id_and_platform
agilent_files([contents,Contents,version,unknown_version], Modules):-
    (
    gse_id([contents,Contents],[]),
    Newadd=['download_geo_dataset',[contents,Contents,version,unknown_version,filetype,agilent_files]];
    gse_id_and_platform([contents,Contents],[]),
    Newadd=['download_geo_dataset_GPL',[contents,Contents,version,unknown_version,filetype,agilent_files]]
     ),
    append([],Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate idat_files with version unknown from gse_id;
% or from gse_id_and_platform
idat_files([contents,Contents,version,unknown_version], Modules):-
    (
    gse_id([contents,Contents],[]),
    Newadd=['download_geo_dataset',[contents,Contents,version,unknown_version,filetype,idat_files]];
    gse_id_and_platform([contents,Contents],[]),
    Newadd=['download_geo_dataset_GPL',[contents,Contents,version,unknown_version,filetype,idat_files]]
     ),
    append([],Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate gpr_file from geo_dataset
gpr_files([contents,Contents,version,gpr],Modules):-
    gpr_files([contents,Contents,version,unknown_version],Past_Modules),
    Newadd=[extract_gpr, [contents,Contents,version,gpr]],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate illumina file from geo_dataset
idat_files([contents,Contents,version,illumina],Modules):-
    idat_files([contents,Contents,version,unknown_version],Past_Modules),
    Newadd=[extract_illumina_idat_files, [contents,Contents,version,illumina]],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate agilent_file from geo_dataset
agilent_files([contents,Contents,version,agilent],Modules):-
    agilent_files([contents,Contents,version,unknown_version],Past_Modules),
    Newadd=[extract_agilent_files, [contents,Contents,version,agilent]],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
% Preprocess geo_dataset with gpr format to signal_file by loess,
signal_raw(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    convert_parameters_raw([contents,Contents,preprocess,loess,is_logged,no_logged,format,pcl],NewParameters),
    Parameters=NewParameters,	
    gpr_files([contents,Contents,version,gpr],Past_Modules),
    Newadd=[loess,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% Preprocess idat_files with illumine format to signal_file by illumina,
signal_raw(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    (convert_parameters_raw([contents,Contents,preprocess,illumina,is_logged,no_logged,format,gct],NewParameters);
     convert_parameters_raw([contents,Contents,preprocess,illumina_controls,is_logged,no_logged,format,gct],NewParameters)),
    get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip,ill_custom_manifest],[],Options),
    append(NewParameters,Options,NewParameters1),
    Parameters = NewParameters1,	
    idat_files([contents,Contents,version,illumina],Past_Modules),
    Newadd=[illumina,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% Preprocess agilent_files with agilent format to signal_file by agilent,
signal_raw(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    convert_parameters_raw([contents,Contents,preprocess,agilent,is_logged,no_logged,format,tdf],NewParameters),
    Parameters=NewParameters,	
    agilent_files([contents,Contents,version,agilent],Past_Modules),
    Newadd=[preprocess_agilent,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% Preprocess v3_4 cel_file by MAS5 or RMA to signal_file,
signal_raw(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    member((Preprocess, Is_Logged),[(mas5, no_logged),(rma, logged)]),
    convert_parameters_raw([contents,Contents,preprocess,Preprocess, is_logged, Is_Logged,format,jeffs,has_missing_value,no_missing],NewParameters),
    Parameters=NewParameters,
    cel_files([contents,Contents,version,v3_4], Past_Modules),
    Newadd=[preprocess,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% log the signal file with Is_Logged is unknown_logged or no_logged.
signal_raw(Parameters,Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    member(OldIs_Logged,[unknown_logged,no_logged]),
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    set_value(Parameters,is_logged,OldIs_Logged,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    Newadd=[log_if_not_log,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% change the format of signal_file to pcl.
signal_raw(Parameters,Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    member(OldFormat,[tdf,res,gct,jeffs]),
    set_value(Parameters,format,OldFormat,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    Newadd=['convert_to_pcl_if_not_pcl',Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% change the Format of the signal file to tdf with Format unknown,
signal_raw(Parameters, Modules):-
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    OldFormat=unknown_format,
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,format,OldFormat,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    Newadd=[convert_to_tdf_if_not_tdf,Parameters],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
% filter genes with the missing value
% The signal_file format is pcl.
% filtering occurs before zero filling and gene ordering,
% Is_logged is logged.
signal_raw(Parameters,Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    Has_Missing_Value=unknown_missing,
    get_value(Parameters,filter,0,Filter),
    not(atom(Filter)),
    Filter>0,
    OldFilter= 0,
    set_value(Parameters,filter,OldFilter,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    Newadd=[gene_filter,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% preprocess_dataset
% The signal_file format is pcl.
% Is_logged is unknown_logged or no_logged.
signal_raw(Parameters,Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    member(Is_Logged,[no_logged,unknown_logged]),
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,predataset,no_predataset,Predataset),
    Predataset=yes_predataset,
    OldPredataset=no_predataset,
    %not(atom(Filter_fc)),
    %OldFilter_fc=no_filter_fc,
    %set_value(Parameters,filter_fc,OldFilter_fc,OldParameters1),
    set_value(Parameters,predataset,OldPredataset,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    Newadd=[preprocessdataset,Parameters],
    append(Past_Modules, Newadd, Modules).

/*--------------------------------------------------------------------------*/
% zero filling the missing value,
% zero filling occurs before gene ordering,
% the format is pcl and Is_logged is logged.
signal_raw(Parameters,Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    Has_Missing_Value = zero_fill,
    Oldmissing = unknown_missing,
    set_value(Parameters,has_missing_value,Oldmissing,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    Newadd=[zero_fill_if_missing,Parameters],
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
% median filling the missing value,
% the format is pcl and Is_logged is logged.
signal_raw(Parameters,Modules):-
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    Has_Missing_Value = median_fill,
    Oldmissing=unknown_missing,
    set_value(Parameters,has_missing_value,Oldmissing,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    Newadd=[median_fill_if_missing,Parameters],
    append(Past_Modules, Newadd, Modules).

/*--------------------------------------------------------------------------*/
%rename the sample name if required
% rename sample occurs after checking missing and check_logged,and check_pcl.
signal_raw(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    rename_list_file([contents,Contents],[]),
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    member(Has_Missing_Value,[median_fill,zero_fill,no_missing]),
    get_value(Parameters,rename_sample,no_rename,Rename_sample),
    Rename_sample=yes_rename,
    OldRename = no_rename,
    set_value(Parameters,rename_sample,OldRename,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    Newadd=[rename_sample,Parameters],
    append(Past_Modules, Newadd, Modules).