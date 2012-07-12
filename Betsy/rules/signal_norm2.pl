%signal_norm2.pl

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).
/*-------------------------------------------------------------------------*/
% Output interface
% The parameters in output can be any length and it will trace to the full length one
% for the parameter which is not provide, will be in no_gene_center/no_gene_normalize

signal_norm2(Parameters,Modules):-
    % Conditions: the Parameters is not full length
    get_desire_parameters_norm2(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_norm2,N1),
    N<N1,
    % Input: signal_norm2 with full length Parameters
    convert_parameters_norm2(Parameters,NewParameters),
    signal_norm2(NewParameters,Modules).
/*-------------------------------------------------------------------------*/
% Input interface
% generate signal_norm2 from signal_norm1 with no_gene_center/no_gene_normalize

signal_norm2(Parameters,Modules):-
   % Conditions: full length and no_gene_center,no_gene_normalize
   get_desire_parameters_norm2(Parameters,NewParameters1),
   length(NewParameters1,N),
   get_length(n_norm2,N1),
   N=N1,
   get_value(Parameters,gene_center,no_gene_center,Gene_Center),
   Gene_Center=no_gene_center,
   get_value(Parameters,gene_normalize,no_gene_normalize,Gene_Normalize),
   Gene_Normalize=no_gene_normalize,
   % Input:signal_norm1
   get_desire_parameters_norm1(Parameters,NewParameters),
   get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip,num_factors],[],Options),
   append(NewParameters,Options,NewParameters2),
   signal_norm1(NewParameters2,Modules).

/*-------------------------------------------------------------------------*/
%gene-center the signal file
signal_norm2(Parameters, Modules):-
    % Conditions: Parameters has created,gene_center in [mean,median],no_gene_normalize
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,gene_center,no_gene_center,Gene_Center),
    member(Gene_Center,[mean,median]),
    get_value(Parameters,gene_normalize,no_gene_normalize,Gene_Normalize),
    Gene_Normalize=no_gene_normalize,
    % Input:signal_norm2 with no_gene_center,no_gene_normalize,different status
    OldGene_center=no_gene_center,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_center,OldGene_center,OldParameters),
    signal_norm2(OldParameters,Past_Modules),
    % Module: centering
    % Output parameters: full length parameters of signal_norm2
    Newadd=[centering,Parameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
%gene_normalize the signal file
signal_norm2(Parameters, Modules):-
    % Conditions:Parameters has created,gene_normalize in [variance,sum_of_squares]
    get_value(Parameters,status,created,Status),
    Status=created,
    get_value(Parameters,gene_normalize,no_gene_normalize,Gene_Normalize),
    member(Gene_Normalize, [variance, sum_of_squares]),
    % Input: signal_norm2 with no_gene_normalize and different status
    OldGene_normalize=no_gene_normalize,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_normalize,OldGene_normalize,OldParameters),
    signal_norm2(OldParameters,Past_Modules),
    % Module:normalize
    % Output parameters: full length parameters of signal_norm2
    Newadd=[normalize,Parameters],
    append(Past_Modules, Newadd, Modules).

