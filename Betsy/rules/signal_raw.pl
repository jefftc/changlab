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
:- dynamic geneset_file/2.
/*-------------------------------------------------------------------------*/
% Input interface
% given a input_signal_file,generate signal_file that will have a full length of parameters
% the parameter which is not provided is in default value
signal_raw(NewParameters,Modules):-
    input_signal_file(OldParameters,Modules),
    convert_parameters_raw(OldParameters,NewParameters).
/*-------------------------------------------------------------------------*/
% reorder the parameters from the gpr_files with parameters in unwanted order
gpr_files([contents,Contents,version,Version],Modules):-
   gpr_files([version,Version,contents,Contents],Modules).

% reorder the parameters from the cel_files with parameters in unwanted order
cel_files([contents,Contents,version,Version],Modules):-
   cel_files([version,Version,contents,Contents],Modules).

% reorder the parameters from the idat_files with parameters in unwanted order
idat_files([contents,Contents,version,Version],Modules):-
   idat_files([version,Version,contents,Contents],Modules).

% reorder the parameters from the agilent_files with parameters in unwanted order
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
% generate geo_dataset with version cc or v3_4 from cel_file with version unknown.
cel_files([contents,Contents,version,cc_or_v3_4], Modules):-
    % Input:cel_file with version=unknown_version
    Oldversion=unknown_version,
    cel_files([contents,Contents,version,Oldversion], Past_Modules),
    % Module: extract_CEL_files
    % Output parameters:[contents,Contents,version,v3_4]
    Newadd=[extract_CEL_files, [contents,Contents,version,cc_or_v3_4]],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate geo_dataset with version v3_4 from cel_file with version cc or v3_4.
cel_files([contents,Contents,version,v3_4], Modules):-
    % Input:cel_file with version=cc_or_v3_4
    Oldversion=cc_or_v3_4,
    cel_files([contents,Contents,version,Oldversion], Past_Modules),
    % Module: convert_CEL_to_v3_4
    % Output parameters:[contents,Contents,version,v3_4]
    Newadd=[convert_CEL_to_v3_v4, [contents,Contents,version,v3_4]],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
% generate gpr_files with version unknown from gse_id;
% or from gse_id_and_platform
gpr_files([contents,Contents,version,unknown_version], Modules):-
    (% Input: gse_id
    gse_id([contents,Contents],[]),
    % Module:download_geo_GSEID
    % Output parameters:[contents,Contents,version,unknown_version,filetype,gpr_files]
    Newadd=[download_geo_GSEID,[contents,Contents,version,unknown_version,filetype,gpr_files]];
    % Input: gse_id_and_platform
    gse_id_and_platform([contents,Contents],[]),
    % Module:download_geo_GSEID_GPLID
    % Output parameters:[contents,Contents,version,unknown_version,filetype,gpr_files]
    Newadd=[download_geo_GSEID_GPLID,[contents,Contents,version,unknown_version,filetype,gpr_files]]
     ),
    append( [],Newadd,Modules).
/*-------------------------------------------------------------------------*/
% generate cel_files with version unknown from gse_id;
% or from gse_id_and_platform
cel_files([contents,Contents,version,unknown_version], Modules):-
    (% Input: gse_id
    gse_id([contents,Contents],[]),
    % Module: download_geo_GSEID
    % Output parameters: [contents,Contents,version,unknown_version,filetype,cel_files]
    Newadd=[download_geo_GSEID,[contents,Contents,version,unknown_version,filetype,cel_files]];
    % Input: gse_id_and_platform
    gse_id_and_platform([contents,Contents],[]),
    % Module:download_geo_GSEID_GPLID
    % Output parameters:[contents,Contents,version,unknown_version,filetype,cel_files]
    Newadd=[download_geo_GSEID_GPLID,[contents,Contents,version,unknown_version,filetype,cel_files]]
     ),
    append( [],Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate agilent_files with version unknown from gse_id;
% or from gse_id_and_platform
agilent_files([contents,Contents,version,unknown_version], Modules):-
    (% Input: gse_id
    gse_id([contents,Contents],[]),
     % Module: download_geo_GSEID
     % Output parameters: [contents,Contents,version,unknown_version,filetype,agilent_files]
    Newadd=[download_geo_GSEID,[contents,Contents,version,unknown_version,filetype,agilent_files]];
    % Input: gse_id_and_platform
    gse_id_and_platform([contents,Contents],[]),
    % Module:download_geo_GSEID_GPLID
    % Output parameters: [contents,Contents,version,unknown_version,filetype,agilent_files]
    Newadd=[download_geo_GSEID_GPLID,[contents,Contents,version,unknown_version,filetype,agilent_files]]
     ),
    append([],Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate idat_files with version unknown from gse_id;
% or from gse_id_and_platform
idat_files([contents,Contents,version,unknown_version], Modules):-
    (% Input: gse_id
    gse_id([contents,Contents],[]),
    % Module: download_geo_GSEID
    % Output parameters: [contents,Contents,version,unknown_version,filetype,idat_files]
    Newadd=[download_geo_GSEID,[contents,Contents,version,unknown_version,filetype,idat_files]];
    % Input: gse_id_and_platform
    gse_id_and_platform([contents,Contents],[]),
    % Module:download_geo_GSEID_GPLID
    % Output parameters:[contents,Contents,version,unknown_version,filetype,idat_files]
    Newadd=[download_geo_GSEID_GPLID,[contents,Contents,version,unknown_version,filetype,idat_files]]
     ),
    append([],Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate gpr_files from gpr_files with unknown_version
gpr_files([contents,Contents,version,gpr],Modules):-
    % Input: gpr_files with unknown_version
    gpr_files([contents,Contents,version,unknown_version],Past_Modules),
    % Module: extract_gpr_files
    % Output parameters:[contents,Contents,version,gpr]
    Newadd=[extract_gpr_files, [contents,Contents,version,gpr]],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate illumina file from idat_files with unknown_version
idat_files([contents,Contents,version,illumina],Modules):-
    % Input: idat_files with unknown_version
    idat_files([contents,Contents,version,unknown_version],Past_Modules),
    % Module: extract_illumina_idat_files
    % Output parameters:[contents,Contents,version,illumina]
    Newadd=[extract_illumina_idat_files, [contents,Contents,version,illumina]],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% generate agilent_file from geo_dataset
agilent_files([contents,Contents,version,agilent],Modules):-
    % Input:agilent_files with unknown_version
    agilent_files([contents,Contents,version,unknown_version],Past_Modules),
    % Module:extract_agilent_files
    % Output parameters:[contents,Contents,version,agilent]
    Newadd=[extract_agilent_files, [contents,Contents,version,agilent]],
    append(Past_Modules, Newadd, Modules).

/*-------------------------------------------------------------------------*/
% Preprocess geo_dataset with gpr format to signal_file by loess,
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has loess, no_logged and tdf
    get_value(Parameters,contents,[unknown],Contents),
    convert_parameters_raw([contents,Contents,preprocess,loess,is_logged,no_logged,format,tdf],NewParameters),
    Parameters=NewParameters,	
    %Input: gpr_files with gpr
    gpr_files([contents,Contents,version,gpr],Past_Modules),
    % Module:normalize_with_loess
    % Output parameters:  full length parameters of signal_raw 
    Newadd=[normalize_with_loess,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% Preprocess idat_files with illumina format to illu_folder by illumina,
illu_folder(Parameters,Modules):-
    % Conditions: Parameters has contents,illumina,no_logged,gct
    get_value(Parameters,contents,[unknown],Contents),
    convert_parameters_raw([contents,Contents,preprocess,illumina,is_logged,no_logged,format,gct],NewParameters),
    get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip,ill_custom_manifest],[],Options),
    append(NewParameters,Options,NewParameters1),
    Parameters = NewParameters1,	
    % Input: idat_files with illumina
    idat_files([contents,Contents,version,illumina],Past_Modules),
    % Module: preprocess_illumina
    % Output parameters: full length parameters of signal_raw 
    Newadd=[preprocess_illumina,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% get the control_file with illumina from illu_folder
control_file(Parameters,Modules):-
   convert_parameters_variable_raw(Parameters,NewParameters),
   % Input: illu_folder
   illu_folder(NewParameters,Past_Modules),
   % Module: get_illumina_control
   % Output parameters: full length parameters of signal_raw 
   Newadd=[get_illumina_control,NewParameters],
   append(Past_Modules,Newadd,Modules).
/*-------------------------------------------------------------------------*/
% get the signal_raw with illumina from illu_folder
signal_raw(Parameters,Modules):-
   % Conditions: preprocess is illumina
   get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
   Preprocess=illumina,
   % Input: illu_folder
   illu_folder(Parameters,Past_Modules),
   % Module: get_illumina_signal
   % Output parameters: full length parameters of signal_raw 
   Newadd=[get_illumina_signal,Parameters],
   append(Past_Modules,Newadd,Modules).
/*-------------------------------------------------------------------------*/
% Preprocess agilent_files with agilent format to signal_file by agilent,
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has contents,agilent,no_logged,tdf
    get_value(Parameters,contents,[unknown],Contents),
    convert_parameters_raw([contents,Contents,preprocess,agilent,is_logged,no_logged,format,tdf],NewParameters),
    Parameters=NewParameters,	
    % Input: agilent_files with agilent
    agilent_files([contents,Contents,version,agilent],Past_Modules),
    % Module: preprocess_agilent
    % Output parameters: full length parameters of signal_raw 
    Newadd=[preprocess_agilent,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% Preprocess v3_4 cel_file by MAS5 or RMA to signal_file,
signal_raw(Parameters,Modules):-
    % Conditions: parameters has (mas5 ,no_logged contents,jeffs,no_missing) or (rma,logged,contents,jeffs,no_missing)
    get_value(Parameters,contents,[unknown],Contents),
    member((Preprocess, Is_Logged,Module),[(mas5, no_logged,preprocess_mas5),(rma, logged,preprocess_rma)]),
    convert_parameters_raw([contents,Contents,preprocess,Preprocess, is_logged, Is_Logged,format,jeffs,has_missing_value,no_missing],NewParameters),
    Parameters=NewParameters,
    % Input: cel_files with v3_4
    cel_files([contents,Contents,version,v3_4], Past_Modules),
    % Module:preprocess_rma or preprocess_mas5
    % Output parameters: full length parameters of signal_raw
    Newadd=[Module,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% log the signal file with Is_Logged is unknown_logged or no_logged.
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has created,logged,pcl
    get_value(Parameters,status,given,Status),
    Status=created,
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    % Input: signal_raw has no_logged or unknown_logged,different status,pcl
    member(OldStatus,[given,created,jointed,splited]),
    member(OldIs_Logged,[unknown_logged,no_logged]),
    set_value(Parameters,is_logged,OldIs_Logged,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    % Module: log_signal
    % Output parameters: full length parameters of signal_raw
    Newadd=[log_signal,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% change the format of signal_file to tdf.
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has created,tdf
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    % Input: signal_raw with different status and format is in [pcl,res,gct,jeffs,unknown_format]
    member(OldStatus,[given,created,jointed,splited]),
    member(OldFormat,[pcl,res,gct,jeffs,unknown_format]),
    set_value(Parameters,format,OldFormat,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    % Module:convert_signal_to_tdf
    % Output parameters: full length parameters of signal_raw
    Newadd=[convert_signal_to_tdf,Parameters],
    append(Past_Modules, Newadd, Modules).


/*-------------------------------------------------------------------------*/
% filter genes with the missing value
% filtering occurs before zero_filling and median_fill
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has created,logged,pcl,unknown_missing,filter>0
    get_value(Parameters,status,given,Status),
    Status=created,
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    Has_Missing_Value=unknown_missing,
    get_value(Parameters,filter,0,Filter),
    not(atom(Filter)),
    Filter>0,
    % Input: signal_raw with logged,pcl,unknown_missing,filter=0 and different status
    member(OldStatus,[given,created,jointed,splited]),
    OldFilter= 0,
    set_value(Parameters,filter,OldFilter,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    % Module:filter_genes_by_missing_values
    % Output parameters:full length parameters of signal_raw
    Newadd=[filter_genes_by_missing_values,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
% preprocess_dataset
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has created,unknown_logged or no_logged,pcl,yes_predataset
    get_value(Parameters,status,given,Status),
    Status=created,
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    member(Is_Logged,[no_logged,unknown_logged]),
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,predataset,no_predataset,Predataset),
    Predataset=yes_predataset,
    % Input: signal_raw with no_predataset,unknown_logged or no_logged,pcl and different status
    member(OldStatus,[given,created,jointed,splited]),
    OldPredataset=no_predataset,
    set_value(Parameters,predataset,OldPredataset,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    % Module: filter_and_threshold_genes
    % Output parameters:full length parameters of signal_raw
    Newadd=[filter_and_threshold_genes,Parameters],
    append(Past_Modules, Newadd, Modules).

/*--------------------------------------------------------------------------*/
% zero filling the missing value,
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has created,logged,pcl,zero_fill
    get_value(Parameters,status,given,Status),
    Status=created,
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    Has_Missing_Value = zero_fill,
    % Input: signal_raw with logged,pcl,unknown_missing and different status
    member(OldStatus,[given,created,jointed,splited]),
    Oldmissing = unknown_missing,
    set_value(Parameters,has_missing_value,Oldmissing,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    % Module: fill_missing_with_zeros
    % Output parameters:full length parameters of signal_raw
    Newadd=[fill_missing_with_zeros,Parameters],
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
% median filling the missing value,
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has created,logged,pcl,median_fill
    get_value(Parameters,status,given,Status),
    Status=created,
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    Has_Missing_Value = median_fill,
    % Input: signal_raw with logged,pcl,unknown_missing and different status
    member(OldStatus,[given,created,jointed,splited]),
    Oldmissing=unknown_missing,
    set_value(Parameters,has_missing_value,Oldmissing,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    % Module: fill_missing_with_median
    % Output parameters:full length parameters of signal_raw
    Newadd=[fill_missing_with_median,Parameters],
    append(Past_Modules, Newadd, Modules).

/*--------------------------------------------------------------------------*/
%rename the sample name if required
signal_raw(Parameters,Modules):-
    % Conditions: Parameters has created,logged,pcl,
    %    has_missing_value in [median_fill,zero_fill,no_missing],yes_rename
    get_value(Parameters,status,given,Status),
    Status=created,
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=tdf,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    member(Has_Missing_Value,[median_fill,zero_fill,no_missing]),
    get_value(Parameters,rename_sample,no_rename,Rename_sample),
    Rename_sample=yes_rename,
    % Input: rename_list_file and signal_raw with logged,pcl,
    %    has_missing_value in [median_fill,zero_fill,no_missing],no_rename and different status
    get_value(Parameters,contents,[unknown],Contents),
    rename_list_file([contents,Contents],[]),
    member(OldStatus,[given,created,jointed,splited]),
    OldRename = no_rename,
    set_value(Parameters,rename_sample,OldRename,OldParameters1),
    set_value(OldParameters1,status,OldStatus,OldParameters),
    signal_raw(OldParameters,Past_Modules),
    % Module: relabel_samples
    % Output parameters:full length parameters of signal_raw
    Newadd=[relabel_samples,Parameters],
    append(Past_Modules, Newadd, Modules).

/*--------------------------------------------------------------------------*/
% plot MA figures for agilent_files
ma_plot(Parameters,Modules):-
    % Conditions: Parameters has contents,agilent,no_logged,tdf
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,version,unknown_version,Version),
    Version=agilent,
    % Input: agilent_files with agilent
    agilent_files([contents,Contents,version,agilent],Past_Modules),
    % Module: plot_MA
    % Output parameters: full length parameters of signal_raw 
    Newadd=[plot_MA,Parameters],append(Past_Modules, Newadd, Modules).

