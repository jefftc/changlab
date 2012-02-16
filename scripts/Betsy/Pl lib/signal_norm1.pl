%signal_norm1.pl

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

/*-------------------------------------------------------------------------*/
% Output interface
% The parameters in output can be any length and it will trace to the full length one
% for the parameter which is not provide, will be in 


signal_norm1(DatasetId,Contents,Parameters,Modules):-
    get_desire_parameters_norm1(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_norm1,N1),
    N<N1,
    convert_parameters_norm1(NewParameters1,NewParameters),
    signal_norm1(DatasetId,Contents,NewParameters,Modules),
    write(NewParameters),nl.
/*-------------------------------------------------------------------------*/
% Input interface
% generate signal_norm1 from signal_clean with 
%no_combat/no_shiftscale/no_dwd/no_quantile
signal_norm1(DatasetId,Contents,Parameters,Modules):-
   length(Parameters,N),
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
   get_desire_parameters_raw(Parameters,NewParameters),
   signal_clean(DatasetId,Contents,NewParameters,Modules).
/*--------------------------------------------------------------------------*/
% Quantile the signal file,
signal_norm1(DatasetId,Contents,Parameters, Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,quantile,no_quantile,Quantile),
    Quantile=yes_quantile,
    OldQuantile=no_quantile,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,quantile,OldQuantile,OldParameters),
    signal_norm1(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[quantile,Write_list],
    append(Past_Modules, Newadd, Modules).
    
/*--------------------------------------------------------------------------*/
% Combat the signal file,
signal_norm1(DatasetId,Contents,Parameters, Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,combat,no_combat,Combat),
    Combat=yes_combat,
    OldCombat=no_combat,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,combat,OldCombat,OldParameters),
    signal_norm1(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[combat,Write_list],
    append(Past_Modules, Newadd, Modules).
    
/*--------------------------------------------------------------------------*/
% Dwd the signal file,
signal_norm1(DatasetId,Contents,Parameters, Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,dwd,no_dwd,Dwd),
    Dwd=yes_dwd,
    OldDwd=no_dwd,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,dwd,OldDwd,OldParameters),
    signal_norm1(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[dwd,Write_list],
    append(Past_Modules, Newadd, Modules).
    
/*--------------------------------------------------------------------------*/
% Shiftscale the signal file,
signal_norm1(DatasetId,Contents,Parameters, Modules):-
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,shiftscale,no_shiftscale,Shiftscale),
    Shiftscale=yes_shiftscale,
    OldShiftscale=no_shiftscale,
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,shiftscale,OldShiftscale,OldParameters),
    signal_norm1(DatasetId, Contents,OldParameters,Past_Modules),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[shiftscale,Write_list],
    append(Past_Modules, Newadd, Modules).