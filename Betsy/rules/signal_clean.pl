%signal_clean.pl

/*-------------------------------------------------------------------------*/
% Output interface
% The parameters in output can be any length and it will trace to the full length one
% for the parameter which is not provide, will be a variable

signal_clean(Parameters,Modules):-
    get_desire_parameters_raw(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_raw,N1),
    N<N1,
    convert_parameters_variable_raw(Parameters,NewParameters),
    signal_clean(NewParameters,Modules).

/*-------------------------------------------------------------------------*/
% Input interface
% generate signal_clean from signal_raw with logged,pcl,
%has_missing_value=no_missing,zero_fill,median_fill
signal_clean(Parameters,Modules):-
    get_desire_parameters_raw(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_raw,N1),
    N=N1,
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    member(Has_Missing_Value,[no_missing,zero_fill,median_fill]),
    signal_raw(Parameters,Modules).

/*-------------------------------------------------------------------------*/
%merge the signal file
signal_clean(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    length(Contents,N),N>1,
    make_psubset(Contents,C1),
    difference(Contents,C1,C2),
    length(C1,N2),N2>0,
    length(C2,N3),N3>0,
    get_value(Parameters,status,given,Status),
    Status=jointed,
    member(OldStatus1,[given,jointed,splited,created]),
    member(OldStatus2,[given,jointed,splited,created]),
    set_value(Parameters,status,OldStatus1,Parameters1),
    set_value(Parameters,status,OldStatus2,Parameters2),
    set_value(Parameters1,contents,C1,Parameters3),
    set_value(Parameters2,contents,C2,Parameters4),
    signal_clean(Parameters3,Past_Modules_1),
    signal_clean(Parameters4,Past_Modules_2),
    append(C1,C2,C3),
    C3=Contents,
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append([merge1,C1,merge2,C2],Parameters,Write_list),
    Newadd=[merge_data,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
%split the signal_file
signal_clean(Parameters, Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    length(Contents,N),N>0,
    get_value(Parameters,status,unknown_status,Status),
    Status=splited,
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,C,preprocess,Preprocess,status,given],Past_Modules_1),
    is_subset(Contents,C),
    not(Contents=C),
    member(Oldstatus,[given,splited,created]),
    set_value(Parameters,status,Oldstatus,NewParameters2),
    set_value(NewParameters2,contents,C,Parameters3),
    signal_clean(Parameters3,Past_Modules_2),
    append([precontents,C],Parameters,Write_list),
    Newadd=[split_data,Write_list],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd,Modules).
/*-------------------------------------------------------------------------*/
% the signal content is unknown,
%pull out the data set required
signal_clean(Parameters,Modules):-
    get_value(Parameters,contents,[unknown],Contents),
    not(member(unknown,Contents)),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,given],Past_Modules_1),
    get_value(Parameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,contents,[unknown],OldParameters),
    signal_clean(OldParameters,Past_Modules_2),
    Newadd=[pull_out_dataset,Parameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules,Newadd,Modules).
