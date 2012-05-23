/*class_label_file.pl*/
:- dynamic class_label_file/2.
class_label_file([contents,Contents,status,Status],Modules):-
   class_label_file([status,Status,contents,Contents],Modules).
/*----------------------------------------------------------------------*/


class_label_file([contents,[X],preprocess,Preprocess,status,created],Modules):-
    not(class_label_file([contents,[X],status,given],[])),
    signal_raw([contents,[X],preprocess,Preprocess,status,created],Past_Modules),
    Newadd=[make_class_label_file,[contents,[X]]],
    append(Past_Modules,Newadd,Modules),!.

class_label_file([contents,Contents,preprocess,Preprocess,status,given],[]):-
    class_label_file([contents,Contents,status,given],[]),!.
/*----------------------------------------------------------------------*/
%%class_label_file([contents,+Contents,status,joined],-Modules)
%merge different class_label_files to generate one class_label_file.
class_label_file([contents,Contents,preprocess,Preprocess,status,joined],Modules):-
    length(Contents,N),N>1,
    make_psubset(Contents,C1),
    difference(Contents,C1,C2),
    class_label_file([contents,C1,preprocess,Preprocess,status,_],Past_Modules_1),
    class_label_file([contents,C2,preprocess,Preprocess,status,_],Past_Modules_2),
    append(C1,C2,C3),
    C3=Contents,
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    Newadd=[join_class_label_files,[contents,Contents,merge1,C1,merge2,C2]],
    append(Past_Modules, Newadd, Modules).
/*----------------------------------------------------------------------*/
%%class_label_file(split, -Modules)
% split class_label_file with Contents from a class 
%label file 
class_label_file([contents,Contents,preprocess,Preprocess,status,split],Modules):-
    length(Contents,N),N>0,
    class_label_file([contents,C,preprocess,Preprocess,status,given],Past_Modules),
    is_subset(Contents,C),
    not(Contents=C),
    Newadd=[split_class_label,[contents,Contents,precontents,C]],
    append(Past_Modules,Newadd,Modules).
/*----------------------------------------------------------------------*/
