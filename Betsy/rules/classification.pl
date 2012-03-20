/*classification.pl*/

%% svm_model(DatasetId, +Contents,+Parameters,-Modules)
% generate a svm model file 
svm_model(DatasetId,Contents,Parameters,Modules):-
    class_label_file(DatasetId,Contents,_,Past_Modules_2),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=pcl,
    signal_file(DatasetId,Contents,NewParameters,Past_Modules_1),
    
    append(Past_Modules_2,Past_Modules_1,Past_Modules),
    append(['DatasetId',
              DatasetId,'Contents',Contents],NewParameters,Write_list),
    Newadd=[train_svm_model,Write_list],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/    
%% svm_predictions(DatasetId,+TrainContents,+TestContents,Parameters,-Modules)
% classification the test data set using model

svm_predictions(DatasetId,TrainContents,TestContents,Parameters,Modules):-
    (class_label_file(DatasetId,TestContents,[status,given],Past_Modules_4);
     not(class_label_file(DatasetId,TestContents,[status,given],_)),
     class_label_file(DatasetId,TestContents,[status,created],Past_Modules_4)),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=pcl,
    signal_file(DatasetId,TestContents,NewParameters,Past_Modules_1),
    svm_model(DatasetId, TrainContents,Parameters,Past_Modules_2),
    append(Past_Modules_4,Past_Modules_2,Past_Modules_3),
    append(Past_Modules_3,Past_Modules_1,Past_Modules),
    append(['DatasetId',DatasetId,'TestContents',TestContents,
         'TrainContents',TrainContents],NewParameters,Write_list),
    Newadd=[apply_svm_model,Write_list],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/
%%WeightedVotingXValidation(DatasetId,+Contents,+Parameters,-Modules)
%classification the data set using Weighted votes method

weightedVotingXValidation(DatasetId, Contents,Parameters,Modules):-
    class_label_file(DatasetId,Contents,_,Past_Modules_2),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=gct,
    signal_file(DatasetId,Contents,NewParameters,Past_Modules_1),
    append(Past_Modules_2,Past_Modules_1,Past_Modules),
    append(['DatasetId',
              DatasetId,'Contents',Contents],NewParameters,Write_list),
    Newadd=[weighted_voting_x_validation,Write_list],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/
weightedVoting(DatasetId,TrainContents,TestContents,Parameters,Modules):-
    (class_label_file(DatasetId,TestContents,[status,given],Past_Modules_6);
     not(class_label_file(DatasetId,TestContents,[status,given],_)),
     class_label_file(DatasetId,TestContents,[status,created],Past_Modules_6)),
    class_label_file(DatasetId,TrainContents,_,Past_Modules_4),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=gct,
    signal_file(DatasetId,TestContents,NewParameters,Past_Modules_1),
    signal_file(DatasetId, TrainContents,NewParameters,Past_Modules_2),
    
    append(Past_Modules_6,Past_Modules_4,Past_Modules_3),
    append(Past_Modules_3,Past_Modules_2,Past_Modules_5),
    append(Past_Modules_5,Past_Modules_1,Past_Modules),
    append(['DatasetId',DatasetId,'TestContents',TestContents,'TrainContents',TrainContents],
           NewParameters,Write_list),
    Newadd=[weighted_voting,Write_list],
    append(Past_Modules,Newadd,Modules).


