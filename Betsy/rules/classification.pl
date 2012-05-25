/*classification.pl*/

%% svm_model(+Parameters,-Modules)
% generate a svm model file 

svm_model(Parameters,Modules):-
    convert_parameters_svm(Parameters,NewParameters1),
    (member(traincontents,Parameters),
    get_value(NewParameters1,traincontents,[unknown],TrainContents),
    signal_file(NewParameters1,Past_Modules_1),
    set_value(NewParameters1,contents,TrainContents,NewParameters),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,TrainContents,preprocess,Preprocess,status,_],Past_Modules_2);
    not(member(traincontents,Parameters)),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_2),
    get_value(NewParameters1,format,unknown_format,Format),
    Format=pcl,
    get_value(NewParameters1,status,created,Status),
    Status=created,
    member(OldStatus,[given,jointed,splited,created]),
    set_value(NewParameters1,status,OldStatus,NewParameters),
    signal_file(NewParameters,Past_Modules_1)),
    append(Past_Modules_2,Past_Modules_1,Past_Modules),
    Newadd=[train_svm_model,NewParameters],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/    
%% svm_predictions(Parameters,-Modules)
% classification the test data set using model
svm_predictions(Parameters,Modules):-
    convert_parameters_svm(Parameters,NewParameters),
    svm_model(NewParameters,Past_Modules),
    get_value(Parameters,testcontents,[unknown],TestContents),
    set_value(NewParameters,contents,TestContents,NewParameters1),
    Newadd=[apply_svm_model,NewParameters1],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/
weightedVoting(Parameters,Modules):-
    get_value(Parameters,traincontents,[unknown],TrainContents),
    get_value(Parameters,testcontents,[unknown],TestContents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,TestContents,preprocess,Preprocess,status,_],Past_Modules_6),
    class_label_file([contents,TrainContents,preprocess,Preprocess,status,_],Past_Modules_4),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=gct,
    set_value(NewParameters,contents,TestContents,NewParameters1),
    signal_file(NewParameters1,Past_Modules_1),
    set_value(NewParameters,contents,TrainContents,NewParameters2),
    signal_file(NewParameters2,Past_Modules_2),
    append(Past_Modules_6,Past_Modules_4,Past_Modules_3),
    append(Past_Modules_3,Past_Modules_2,Past_Modules_5),
    append(Past_Modules_5,Past_Modules_1,Past_Modules),
    append([testcontents,TestContents,traincontents,TrainContents],
           NewParameters,Write_list),
    get_options(Parameters,[wv_num_features,wv_minstd,wv_feature_stat],[],Options),
    append(Write_list,Options,Write_list1),
    Newadd=[weighted_voting,Write_list1],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/
loocv(Parameters,Modules):-
    get_value_variable(Parameters,classification,Classification),
    member(Classification,[svm,weightedvoting]),
    get_CV_parameter(Classification,Required_format,List),
    get_value(Parameters,traincontents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_2),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=Required_format,
    set_value(NewParameters,contents,Contents,NewParameters1),
    signal_file(NewParameters1,Past_Modules_1),
    append(Past_Modules_2,Past_Modules_1,Past_Modules),
    append(List,NewParameters1,Write_list),
    Newadd=[loocv,Write_list],
    append(Past_Modules,Newadd,Modules).

get_CV_parameter(A,Format,List):-
    A=svm,
    Format=pcl,
    List=[train_model,train_svm_model,
          predict,apply_svm_model];
    A=weightedvoting,
    Format=gct,
    List=[predict,weighted_voting].
/*---------------------------------------------------------------*/
prediction_plot(Parameters,Modules):-
    get_value_variable(Parameters,class_plot,Class_plot),
    member(Class_plot,[svm,weightedvoting,loocv]),
    (Class_plot=loocv,
    loocv(Parameters,Past_Modules);
    Class_plot=weightedvoting,
    weightedVoting(Parameters,Past_Modules);
    Class_plot=svm,
    svm_predictions(Parameters,Past_Modules)),
    append(Parameters,[class_plot,Class_plot],NewParameters),
    Newadd=[plot_prediction,NewParameters],
    append(Past_Modules,Newadd,Modules).

