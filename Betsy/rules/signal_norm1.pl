%signal_norm1.pl

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

/*-------------------------------------------------------------------------*/
% Output interface
% The parameters in output can be any length and it will trace to the full length one
% for the parameter which is not provide, will be in default value

signal_norm1(Parameters,Modules):-
    % Conditions: Parameters is not full length 
    get_desire_parameters_norm1(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_norm1,N1),
    N<N1,
    % Input: signal_norm1 with full length parameters
    convert_parameters_norm1(Parameters,NewParameters),
    signal_norm1(NewParameters,Modules).
/*-------------------------------------------------------------------------*/
% Input interface
% generate signal_norm1 from signal_clean with 
%no_combat/no_shiftscale/no_dwd/no_quantile/no_bfrm
signal_norm1(Parameters,Modules):-
   % Conditions:full length and no_combat,no_shiftscale,no_dwd,no_quantile,no_bfrm
   get_desire_parameters_norm1(Parameters,NewParameters1),
   length(NewParameters1,N),
   get_length(n_norm1,N1),
   N=N1,
   get_value(Parameters,quantile,no_quantile,Quantile),
   Quantile=no_quantile,
   get_value(Parameters,combat,no_combat,Combat),
   Combat=no_combat,
   get_value(Parameters,shiftscale,no_shiftscale,Shiftscale),
   Shiftscale=no_shiftscale,
   get_value(Parameters,dwd,no_dwd,Dwd),
   Dwd=no_dwd,
   get_value(Parameters,bfrm,no_bfrm,Bfrm),
   Bfrm=no_bfrm,
   % Input: signal_clean
   get_desire_parameters_raw(Parameters,NewParameters),
   get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip],[],Options),
   append(NewParameters,Options,NewParameters2),
   signal_clean(NewParameters2,Modules).
/*--------------------------------------------------------------------------*/
% Quantile the signal file,
signal_norm1(Parameters, Modules):-
    % Conditions: Parameters has yes_quantile,created
    get_value(Parameters,quantile,no_quantile,Quantile),
    Quantile=yes_quantile,
    get_value(Parameters,status,created,Status),
    Status=created,
    % Input: signal_norm1 with no_quantile and different status
    OldQuantile=no_quantile, 
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,quantile,OldQuantile,OldParameters),
    signal_norm1(OldParameters,Past_Modules),
    % Module: normalize_samples_with_quantile
    % Output parameters: full length parameters of signal_norm1
    Newadd=[normalize_samples_with_quantile,Parameters],
    append(Past_Modules, Newadd, Modules).
    
/*--------------------------------------------------------------------------*/
% Combat the signal file,
signal_norm1(Parameters, Modules):-
    % Conditions: Parameters has yes_combat,created
    get_value(Parameters,combat,no_combat,Combat),
    Combat=yes_combat,
    get_value(Parameters,status,created,Status),
    Status=created,
    % Input: class_label_file and signal_norm1 with no_combat and different status
    OldCombat=no_combat,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,combat,OldCombat,OldParameters),
    signal_norm1(OldParameters,Past_Modules_2),
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    % Module: normalize_samples_with_combat
    % Output parameters: full length parameters of signal_norm1
    Newadd=[normalize_samples_with_combat,Parameters],
    append(Past_Modules, Newadd, Modules).
    
/*--------------------------------------------------------------------------*/
% Dwd the signal file,
signal_norm1(Parameters, Modules):-
    % Conditions: Parameters with yes_dwd,created
    get_value(Parameters,dwd,no_dwd,Dwd),
    Dwd=yes_dwd,
    get_value(Parameters,status,created,Status),
    Status=created,
    % Input:class_label_file and signal_norm1 with no_dwd and different status
    OldDwd=no_dwd,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,dwd,OldDwd,OldParameters),
    signal_norm1(OldParameters,Past_Modules_2),
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    % Module: normalize_samples_with_dwd
    % Output parameters: full length parameters of signal_norm1
    Newadd=[normalize_samples_with_dwd,Parameters],
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
% Bfrm the signal file,
signal_norm1(Parameters, Modules):-
    % Conditions: Parameters with yes_bfrm,created
    get_value(Parameters,bfrm,no_bfrm,Bfrm),
    Bfrm=yes_bfrm,
    get_value(Parameters,status,created,Status),
    Status=created,
    % Input:signal_norm1 with no_bfrm and different status
    OldBfrm=no_bfrm,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,bfrm,OldBfrm,OldParameters),
    signal_norm1(OldParameters,Past_Modules),
    % Module: normalize_samples_with_bfrm
    % Output parameters: full length parameters of signal_norm1
    Newadd=[normalize_samples_with_bfrm,Parameters],
    append(Past_Modules, Newadd, Modules).

/*--------------------------------------------------------------------------*/
% Shiftscale the signal file,
signal_norm1(Parameters, Modules):-
    % Conditions: Parameters with yes_shiftscale,created
    get_value(Parameters,shiftscale,no_shiftscale,Shiftscale),
    Shiftscale=yes_shiftscale,
    get_value(Parameters,status,created,Status),
    Status=created,
    % Input:class_label_file and signal_norm1 with no_shiftscale and different status
    OldShiftscale=no_shiftscale,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,shiftscale,OldShiftscale,OldParameters),
    signal_norm1(OldParameters,Past_Modules_2),
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    % Module: normalize_samples_with_shiftscale
    % Output parameters: full length parameters of signal_norm1
    Newadd=[normalize_samples_with_shiftscale,Parameters],
    append(Past_Modules, Newadd, Modules).