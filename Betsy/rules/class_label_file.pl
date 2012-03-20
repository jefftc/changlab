/*class_label_file.pl*/

/*----------------------------------------------------------------------*/
%%class_label_file(+DatasetId,+[X],[status,created], -Modules)
%make class_label_file from the signal file with the same content.
class_label_file(DatasetId,[X],[status,created], Modules):-  
    input_signal_file(DatasetId,[X],_,Past_Modules),
    Newadd=['make_class_label_file',['DatasetId',
            DatasetId,'Contents',[X]]],
    append(Past_Modules,Newadd,Modules),!.
/*----------------------------------------------------------------------*/
%%class_label_file(+DatasetId,+Contents, [status,joined],-Modules)
%merge different class_label_files to generate one class_label_file.
class_label_file(DatasetId,Contents,[status,joined],Modules):-
    length(Contents,N),N>1,
    make_psubset(Contents,C1),
    difference(Contents,C1,C2),
    length(DatasetId,N1),    
    (  N1=N,       %consider both DatasetId and Contents are list_N and list_N
       getdatasetid(DatasetId,Contents,C1,[],Y1),
       getdatasetid(DatasetId,Contents,C2,[],Y2);
       N1=1,          %consider DatasetId is length 1 and Contents is list_N 
       Y1=DatasetId,
       Y2=DatasetId
    ),
    class_label_file(Y1,C1,_,Past_Modules_1),
    class_label_file(Y2,C2,_,Past_Modules_2),
    append(C1,C2,C3),
    C3=Contents,
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    Newadd=['join_class_label_files',['DatasetId',
  DatasetId,'Contents',Contents,merge1,C1,merge2,C2,dataset1,Y1,dataset2,Y2]],
    append(Past_Modules, Newadd, Modules).
/*----------------------------------------------------------------------*/
%%class_label_file(+DatasetId,+Contents,split, -Modules)
% split class_label_file with Contents from a class 
%label file 
class_label_file(DatasetId,Contents,[status,split],Modules):-
    length(Contents,N),N>0,
    class_label_file(DatasetId1,C,[status,given],Past_Modules),
    is_subset(Contents,C),
    not(Contents=C),
    is_subset(DatasetId,DatasetId1),
    getdatasetid(DatasetId1,C,Contents,[],Y),
    Y=DatasetId,
    Newadd=['split_class_label',['DatasetId',DatasetId,'PreContents',C,
           'Contents',Contents,'PreDatasetid',DatasetId1]],
    append(Past_Modules,Newadd,Modules).
