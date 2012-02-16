%prolog_utils.pl

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

use_module(library(lists)).

get_value(Parameters,Key,Default,Value):-
   member(Key,Parameters),
   nth0(N,Parameters,Key),
   X is N+1,
   nth0(X,Parameters,Value);
   not(member(Key,Parameters)),
   Value=Default.

/*-------------------------------------------------------------------------*/
get_value_variable(Parameters,Key,Value):-
   member(Key,Parameters),
   nth0(N,Parameters,Key),
   X is N+1,
   nth0(X,Parameters,Value);
   not(member(Key,Parameters)),
   Value=R.
/*-------------------------------------------------------------------------*/
replace([_|T],0,X,[X|T]). 
replace([H|T],I,X,[H|R]):-I1 is I-1, replace(T,I1,X,R).

/*-------------------------------------------------------------------------*/
set_value(Parameters,Key,Value,NewParameters):-
    member(Key,Parameters),
    nth0(N,Parameters,Key),
    replace(Parameters,N+1,Value,NewParameters).
/*-------------------------------------------------------------------------*/
remove_at(X,[X|Xs],1,Xs).
remove_at(X,[Y|Xs],K,[Y|Ys]) :- K > 1, 
   K1 is K - 1, remove_at(X,Xs,K1,Ys).
    
/*-------------------------------------------------------------------------*/
convert_parameters_raw(Parameters,NewParameters):-
    get_value(Parameters,format,unknown_format,Format),
    member(Format,[tdf,jeffs,res,gct,pcl,unknown_format,xls,not_xls]),
    append([],[format,Format],NewParameters1),
   
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    member(Preprocess,[rma,mas5,loess,unknown_preprocess,illumina,agilent,illumina_controls]),
    append(NewParameters1,[preprocess,Preprocess],NewParameters2),

    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    member(Is_Logged,[logged,no_logged,unknown_logged]),
    append(NewParameters2,[is_logged,Is_Logged],NewParameters3),

    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    member(Has_Missing_Value,[yes_missing,no_missing,unknown_missing]),
    append(NewParameters3,[has_missing_value,Has_Missing_Value],NewParameters4),

    get_value(Parameters,filter,no_filter,Filter),
    (Filter=no_filter;not(atom(Filter))),
    append(NewParameters4,[filter,Filter],NewParameters5),
   
    get_value(Parameters,fill,no_fill,Fill),
    member(Fill,[no_fill,yes_fill]),
    append(NewParameters5,[fill,Fill],NewParameters6),
    
    get_value(Parameters,status,created,Status),
    member(Status,[given,created,jointed,splited]),
    append(NewParameters6,[status,Status],NewParameters).

/*-----------------------------------------------*/
convert_parameters_variable_raw(Parameters,NewParameters):-
    get_value_variable(Parameters,format,Format),
    member(Format,[tdf,jeffs,res,gct,pcl,unknown_format,xls,not_xls]),
    append([],[format,Format],NewParameters1),
    
    get_value_variable(Parameters,preprocess,Preprocess),
    member(Preprocess,[rma,mas5,loess,agilent,illumina,unknown_preprocess,illumina_controls]),
    append(NewParameters1,[preprocess,Preprocess],NewParameters2),

    get_value_variable(Parameters,is_logged,Is_Logged),
    member(Is_Logged,[logged,no_logged,unknown_logged]),
    append(NewParameters2,[is_logged,Is_Logged],NewParameters3),
    
    get_value_variable(Parameters,has_missing_value,Has_Missing_Value),
    member(Has_Missing_Value,[yes_missing,no_missing,unknown_missing]),
    append(NewParameters3,[has_missing_value,Has_Missing_Value],NewParameters4),

    get_value(Parameters,filter,no_filter,Filter),
    (member(Filter,[no_filter,25]);not(atom(Filter))),
    (not(atom(Filter)),Has_Missing_Value=yes_missing;Filter=no_filter,true),
    append(NewParameters4,[filter,Filter],NewParameters5),
    
    get_value_variable(Parameters,fill,Fill),
    member(Fill,[no_fill,yes_fill]),
    (Fill=yes_fill,Has_Missing_Value=yes_missing;not(Fill=yes_fill),true),
    append(NewParameters5,[fill,Fill],NewParameters6),
   
    get_value_variable(Parameters,status,Status),
    member(Status,[given,created,jointed,splited]),
    append(NewParameters6,[status,Status],NewParameters).

/*-------------------------------------------------------------------------*/
%get the parameters list for signal_raw
get_desire_parameters_raw(Parameters,NewParameters):-
    (member(format,Parameters),
    get_value_variable(Parameters,format,Format),
    append([],[format,Format],NewParameters1);
    not(member(format,Parameters)),
    NewParameters1=[]),
    
    (member(preprocess,Parameters),
    get_value_variable(Parameters,preprocess,Preprocess),
    append(NewParameters1,[preprocess,Preprocess],NewParameters2);
    not(member(preprocess,Parameters)),
    NewParameters2=NewParameters1),

    (member(is_logged,Parameters),
    get_value_variable(Parameters,is_logged,Is_Logged),
    append(NewParameters2,[is_logged,Is_Logged],NewParameters3);
    not(member(is_logged,Parameters)),
    NewParameters3=NewParameters2),
    
    (member(has_missing_value,Parameters),
    get_value_variable(Parameters,has_missing_value,Has_Missing_Value),
    append(NewParameters3,[has_missing_value,Has_Missing_Value],NewParameters4);
    not(member(has_missing_value,Parameters)),
    NewParameters4=NewParameters3),

    (member(filter,Parameters),
    get_value_variable(Parameters,filter,Filter),
    append(NewParameters4,[filter,Filter],NewParameters5);
    not(member(filter,Parameters)),
    NewParameters5=NewParameters4),
    
    (member(fill,Parameters),
    get_value_variable(Parameters,fill,Fill),
    append(NewParameters5,[fill,Fill],NewParameters6);
    not(member(fill,Parameters)),
    NewParameters6=NewParameters5),

    (member(status,Parameters),
    get_value_variable(Parameters,status,Status),
    append(NewParameters6,[status,Status],NewParameters);
    not(member(status,Parameters)),
    NewParameters=NewParameters6).

/*-------------------------------------------------------------------------*/
convert_parameters_norm1(Parameters,NewParameters):-
    convert_parameters_variable_raw(Parameters,NewParameters1),

    get_value_variable(Parameters,quantile,Quantile),
    member(Quantile,[yes_quantile,no_quantile]),
    append(NewParameters1,[quantile,Quantile],NewParameters2),

    get_value_variable(Parameters,combat,Combat),
    member(Combat,[yes_combat,no_combat]),
    append(NewParameters2,[combat,Combat],NewParameters3),
     
    get_value_variable(Parameters,shiftscale,Shiftscale),
    member(Shiftscale,[yes_shiftscale,no_shiftscale]),
    append(NewParameters3,[shiftscale,Shiftscale],NewParameters4),

    get_value_variable(Parameters,dwd,Dwd),
    member(Dwd,[yes_dwd,no_dwd]),
    append(NewParameters4,[dwd,Dwd],NewParameters).

/*-------------------------------------------------------------------------*/

%get the parameters list for signal_norm1
get_desire_parameters_norm1(Parameters,NewParameters):-
    get_desire_parameters_raw(Parameters,NewParameters1),
   
    (member(quantile,Parameters),
    get_value_variable(Parameters,quantile,Quantile),
    append(NewParameters1,[quantile,Quantile],NewParameters2);
    not(member(quantile,Parameters)),
    NewParameters2=NewParameters1),

    (member(combat,Parameters),
    get_value_variable(Parameters,combat,Combat),
    append(NewParameters2,[combat,Combat],NewParameters3);
    not(member(combat,Parameters)),
    NewParameters3=NewParameters2),

    (member(shiftscale,Parameters),
    get_value_variable(Parameters,shiftscale,Shiftscale),
    append(NewParameters3,[shiftscale,Shiftscale],NewParameters4);
    not(member(shiftscale,Parameters)),
    NewParameters4=NewParameters3),

    (member(dwd,Parameters),
    get_value_variable(Parameters,dwd,Dwd),
    append(NewParameters4,[dwd,Dwd],NewParameters);
    not(member(dwd,Parameters)),
    NewParameters=NewParameters4).

 /*-------------------------------------------------------------------------*/
convert_parameters_norm2(Parameters,NewParameters):-
    convert_parameters_norm1(Parameters,NewParameters1),
    
    get_value_variable(Parameters,gene_center,Gene_Center),
    member(Gene_Center,[mean,median,no_gene_center]),
    append(NewParameters1,[gene_center,Gene_Center],NewParameters2),
   
    get_value_variable(Parameters,gene_normalize,Gene_Normalize),
    member(Gene_Normalize,[variance,sum_of_squares,no_gene_normalize]),
    append(NewParameters2,[gene_normalize,Gene_Normalize],NewParameters).
 /*-------------------------------------------------------------------------*/
%get the parameters list for signal_norm2
get_desire_parameters_norm2(Parameters,NewParameters):-
    get_desire_parameters_norm1(Parameters,NewParameters1),

    (member(gene_center,Parameters),
    get_value_variable(Parameters,gene_center,Gene_Center),
    append(NewParameters1,[gene_center,Gene_Center],NewParameters2);
    not(member(gene_center,Parameters)),
    NewParameters2=NewParameters1),

    (member(gene_normalize,Parameters),
    get_value_variable(Parameters,gene_normalize,Gene_Normalize),
    append(NewParameters2,[gene_normalize,Gene_Normalize],NewParameters);
    not(member(gene_normalize,Parameters)),
    NewParameters=NewParameters2).
 /*-------------------------------------------------------------------------*/
convert_parameters_file(Parameters,NewParameters):-
    convert_parameters_norm2(Parameters,NewParameters1),
    get_value_variable(Parameters,gene_order,Gene_Order),
    member(Gene_Order,[by_sample_ttest,by_gene_list,no_order]),
    append(NewParameters1,[gene_order,Gene_Order],NewParameters).
 /*-------------------------------------------------------------------------*/
%get the parameters list for signal_file
get_desire_parameters_file(Parameters,NewParameters):-
    get_desire_parameters_norm2(Parameters,NewParameters1),
    (member(gene_order,Parameters),
    get_value_variable(Parameters,gene_order,Gene_Order),
    append(NewParameters1,[gene_order,Gene_Order],NewParameters);
    not(member(gene_order,Parameters)),
    NewParameters=NewParameters1).

/*----------------------------------------------------------------------*/
/*Whether [X|Y] is a subset of Z*/
is_subset([X|Y], Z) :- member(X, Z), is_subset(Y, Z).
is_subset([], _).

/*make a subset and put it into Z.*/
make_subset([_|Y],Z):- make_subset(Y,Z).
make_subset([X|Y],Z):- make_subset(Y,W),append([X],W,Z).
make_subset([],[]).

/*make a proper subset and put it into Z.*/
make_psubset(X,Z):- make_subset(X,Z), not(perm(X,Z)).

/*put the intersection of the first two into the Z.*/
intersection([X|Y],M,[X|Z]):- member(X,M),intersection(Y,M,Z).
intersection([X|Y],M,Z):- \+member(X,M),intersection(Y,M,Z).
intersection([],_,[]).

/*Atoms that are in list [X|Y] but not list M. */
difference([X|Y],M,Z):- member(X,M),difference(Y,M,Z).
difference([X|Y],M,[X|Z]):- \+member(X,M),difference(Y,M,Z).
difference([],_,[]).

takeout(X,[X|R],R).
takeout(X,[F|R],[F|S]):- takeout(X,R,S).
perm([X|Y],Z):- perm(Y,W),takeout(X,Z,W).
perm([],[]).
/*-------------------------------------------------------------------------*/    
%%getdatasetid(+DatasetId,+Contents,+C,+W,-Y)
% given DatasetId and Contents both List_N,
% C is a subset of Contents
% find the corresponding DatasetId of C,
% output is Y, W is the initial output.

getdatasetid(DatasetId,Contents,C,W,Y):-
    length(C,N),
    N>0,
    nth0(0,C,X),
    nth0(Z,Contents,X),
    delete(C,X,R),
    nth0(Z,DatasetId,U),
    append(W,[U],T),
    getdatasetid(DatasetId,Contents,R,T,Y).
getdatasetid(DatasetId,Contents,[],R,R).

get_length(A,B):-
    A=n_raw, B=14;
    A=n_norm1,B=22;
    A=n_norm2,B=26;
    A=n_file,B=28.