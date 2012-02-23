%signal_clean.pl

/*-------------------------------------------------------------------------*/
% Output interface
% The parameters in output can be any length and it will trace to the full length one
% for the parameter which is not provide, will be a variable

signal_clean(DatasetId,Contents,Parameters,Modules):-
    get_desire_parameters_raw(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_raw,N1),
    N<N1,
    convert_parameters_variable_raw(NewParameters1,NewParameters),
    signal_clean(DatasetId,Contents,NewParameters,Modules).
/*-------------------------------------------------------------------------*/
% Input interface
% generate signal_clean from signal_raw with logged,pcl,has_missing/fill or no_missing/no_fill

signal_clean(DatasetId,Contents,Parameters,Modules):-
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(Parameters,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    get_value(Parameters,fill,no_fill,Fill),
    member((Has_Missing_Value,Fill),[(yes_missing,yes_fill),(no_missing,no_fill)]),
    signal_raw(DatasetId, Contents,Parameters,Modules).

/*-------------------------------------------------------------------------*/
%merge the signal file
signal_clean(DatasetId,Contents,Parameters,Modules):-
    length(Contents,N),N>1,
    make_psubset(Contents,C1),
    difference(Contents,C1,C2),
    length(C1,N2),N2>0,
    length(C2,N3),N3>0,
    length(DatasetId,N1),
    (  N1=N,       %consider both DatasetId and Contents are list_N and list_N
       getdatasetid(DatasetId,Contents,C1,[],Y1),
       getdatasetid(DatasetId,Contents,C2,[],Y2);
       N1=1,          %consider DatasetId is length 1 and Contents is list_N
       Y1=DatasetId,
       Y2=DatasetId
    ),
    get_value(Parameters,status,given,Status),
    Status=jointed,
    member(OldStatus1,[given,jointed,splited,created]),
    member(OldStatus2,[given,jointed,splited,created]),
    set_value(Parameters,status,OldStatus1,Parameters1),
    set_value(Parameters,status,OldStatus2,Parameters2),
    signal_clean(Y1,C1,Parameters1,Past_Modules_1),
    signal_clean(Y2,C2,Parameters2,Past_Modules_2),
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(['DatasetId',DatasetId,'Contents',Contents,merge1,C1,merge2,C2,dataset1,
    Y1,dataset2,Y2],Parameters,Write_list),
    Newadd=[merge_data,Write_list],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
%split the signal_file
signal_clean(DatasetId, Contents,Parameters, Modules):-
    length(Contents,N),N>0,
    get_value(Parameters,status,unknown_status,Status),
    Status=splited,
    class_label_file(DatasetId1,C,given,Past_Modules_1),
    member(OldFormat,[res,tdf,jeffs,gct,pcl,not_xls]),
    set_value(Parameters,format,OldFormat,NewParameters1),
    set_value(NewParameters1,status,given,Parameters2),
    signal_clean(DatasetId1,C,Parameters2,Past_Modules_2),
    is_subset(Contents,C),
    not(Contents=C),
    is_subset(DatasetId,DatasetId1),
    getdatasetid(DatasetId1,C,Contents,[],Y),
    Y=DatasetId,
    append(['DatasetId',DatasetId,'Contents',Contents,'PreContents',
     C,'PreDatasetid',DatasetId1],Parameters,Write_list),
    Newadd=[split_data,Write_list],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd,Modules).
/*-------------------------------------------------------------------------*/
%given the DatasetId, but the signal content is unknown,
%pull out the data set required
signal_clean(DatasetId,Contents,Parameters,Modules):-
    not(member(unknown,Contents)),
    class_label_file(DatasetId,Contents,given,Past_Modules_1),
    get_value(Parameters,status,given,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    set_value(Parameters,status,OldStatus,OldParameters),
    signal_clean(DatasetId, [unknown],OldParameters,Past_Modules_2),
    append(['DatasetId',
            DatasetId,'Contents',Contents],Parameters,Write_list),
    Newadd=[pull_out_dataset,Write_list],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules,Newadd,Modules).

