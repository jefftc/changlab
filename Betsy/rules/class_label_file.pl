/*class_label_file.pl*/
:- dynamic class_label_file/2.

class_label_file([contents,Contents,status,Status],Modules):-
   % reorder the parameters from the class_label_file with parameters in unwanted order
   class_label_file([status,Status,contents,Contents],Modules).
/*----------------------------------------------------------------------*/
class_label_file([contents,[X],preprocess,Preprocess,status,created],Modules):-
    % Conditions: the class_label_file is not given
    not(class_label_file([contents,[X],status,given],[])),
    % Input:the class_label_file is generated from signal_raw 
    signal_raw([contents,[X],preprocess,Preprocess,status,created],Past_Modules),
    % Module: make_class_label_file 
    % Output parameters:[[contents,[X]]
    Newadd=[make_class_label_file,[contents,[X]]],
    append(Past_Modules,Newadd,Modules),!.

class_label_file([contents,Contents,preprocess,Preprocess,status,given],[]):-
    % Generate the class_label_file with preprocess from the given class_label_file
    class_label_file([contents,Contents,status,given],[]),!.
/*----------------------------------------------------------------------*/
%% class_label_file([contents,+Contents,status,jointed],-Modules)
% merge different class_label_files to generate one class_label_file.
class_label_file([contents,Contents,preprocess,Preprocess,status,jointed],Modules):-
    % Conditions: the contents length should be bigger than one
    length(Contents,N),N>1,
    % Get the input file contents C1 and C2
    make_psubset(Contents,C1),
    difference(Contents,C1,C2),
    % Input: two class_label_file 
    class_label_file([contents,C1,preprocess,Preprocess,status,_],Past_Modules_1),
    class_label_file([contents,C2,preprocess,Preprocess,status,_],Past_Modules_2),
    append(C1,C2,C3),
    C3=Contents,
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    % Module: join_class_label_file
    % Output parameters:[contents,Contents,merge1,C1,merge2,C2]
    Newadd=[join_class_label_files,[contents,Contents,merge1,C1,merge2,C2]],
    append(Past_Modules, Newadd, Modules).
/*----------------------------------------------------------------------*/
% split class_label_file with Contents from a class label file 
class_label_file([contents,Contents,preprocess,Preprocess,status,splited],Modules):-
    % Conditions: the contents length should be bigger than zero
    length(Contents,N),N>0,
    % Input: a class_label_file which is given and its contents contains the output contents
    class_label_file([contents,C,preprocess,Preprocess,status,given],Past_Modules),
    is_subset(Contents,C),
    not(Contents=C),
    % Module: split_class_label
    % Output parameters: [contents,Contents,precontents,C]
    Newadd=[split_class_label,[contents,Contents,precontents,C]],
    append(Past_Modules,Newadd,Modules).
/*----------------------------------------------------------------------*/
