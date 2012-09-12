/*classification.pl*/

%% svm_model(+Parameters,-Modules)
% generate a svm model file 

svm_model(Parameters,Modules):-
    % get the full length parameters
    convert_parameters_svm(Parameters,NewParameters1),
    % Conditions: if traincontents in Parameters, 
    % this is for the case there is testdata to test using this svm_model, 
    % the traindata and testdata need to be aligned
    (member(traincontents,Parameters),
    get_value(NewParameters1,traincontents,[unknown],TrainContents),
    % Input: signal_file and class_label_file 
    signal_file(NewParameters1,Past_Modules_1),
    set_value(NewParameters1,contents,TrainContents,NewParameters),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,TrainContents,preprocess,Preprocess,status,_],Past_Modules_2);
    % Conditions: if traincontents not in Parameters,
    not(member(traincontents,Parameters)),
    get_value(NewParameters1,format,unknown_format,Format),
    Format=tdf,
    get_value(NewParameters1,status,created,Status),
    Status=created,
    get_value(NewParameters1,is_logged,unknown_logged,Is_logged),
    Is_logged=logged,
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    % Input: signal_file and class_label_file
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_2),
    member(OldStatus,[given,jointed,splited,created]),
    set_value(NewParameters1,status,OldStatus,NewParameters),
    signal_file(NewParameters,Past_Modules_1)),
    append(Past_Modules_2,Past_Modules_1,Past_Modules),
    % Module: train_svm_model
    % Output parameters: update the input Parameters to a full length and the contents
    Newadd=[train_svm_model,NewParameters],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/    
%% svm_predictions(Parameters,-Modules)
% classification the test data set using model
svm_predictions(Parameters,Modules):-
    convert_parameters_svm(Parameters,NewParameters),
    % Input: svm_model
    svm_model(NewParameters,Past_Modules),
    % Module:apply_svm_model
    % Output parameters: update the Parameters to full length 
    Newadd=[apply_svm_model,NewParameters],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/
weightedVoting(Parameters,Modules):-
    % Conditions: the format should be gct
    convert_parameters_file(Parameters,NewParameters3),
    set_value(NewParameters3,num_features,0,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=gct,
    % Input: train signal_file,test_signal_file,train class_lable_file,test_class_label_file
    get_value(Parameters,traincontents,[unknown],TrainContents),
    get_value(Parameters,testcontents,[unknown],TestContents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,TestContents,preprocess,Preprocess,status,_],Past_Modules_6),
    class_label_file([contents,TrainContents,preprocess,Preprocess,status,_],Past_Modules_4),
    set_value(NewParameters,contents,TestContents,NewParameters1),
    signal_file(NewParameters1,Past_Modules_1),
    set_value(NewParameters,contents,TrainContents,NewParameters2),
    signal_file(NewParameters2,Past_Modules_2),
    append(Past_Modules_6,Past_Modules_4,Past_Modules_3),
    append(Past_Modules_3,Past_Modules_2,Past_Modules_5),
    append(Past_Modules_5,Past_Modules_1,Past_Modules),
    % Output parameters: update the parameters to a full length
    append([testcontents,TestContents,traincontents,TrainContents],
           NewParameters3,Write_list),
    get_options(Parameters,[wv_minstd,wv_feature_stat],[],Options),
    append(Write_list,Options,Write_list1),
    % Module: classify_with_weighted_voting
    Newadd=[classify_with_weighted_voting,Write_list1],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/
loocv(Parameters,Modules):-
    % Conditions: classification should be svm or weightedvoting 
    %             and format is the required format
    get_value_variable(Parameters,classification,Classification),
    member(Classification,[svm,weightedvoting]),
    get_CV_parameter(Classification,Required_format,List),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=Required_format,
    % Input: class_label_file and signal_file
    get_value(Parameters,traincontents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_2),
    set_value(NewParameters,contents,Contents,NewParameters1),
    signal_file(NewParameters1,Past_Modules_1),
    append(Past_Modules_2,Past_Modules_1,Past_Modules),
    % Output parameters: update the parameters to full length
    append(List,NewParameters1,Write_list),
    % Module:run_loocv
    Newadd=[run_loocv,Write_list],
    append(Past_Modules,Newadd,Modules).

get_CV_parameter(A,Format,List):-
    A=svm,
    Format=tdf,
    List=[train_model,train_svm_model,
          predict,apply_svm_model];
    A=weightedvoting,
    Format=gct,
    List=[predict,classify_with_weighted_voting].
/*---------------------------------------------------------------*/
prediction_plot(Parameters,Modules):-
    % Conditions: class_plot should be svm,weightedvoting or loocv
    get_value_variable(Parameters,class_plot,Class_plot),
    member(Class_plot,[svm,weightedvoting,loocv]),
    % Input: loocv,weightedVoting or svm_predictions
    (Class_plot=loocv,
    loocv(Parameters,Past_Modules);
    Class_plot=weightedvoting,
    weightedVoting(Parameters,Past_Modules);
    Class_plot=svm,
    svm_predictions(Parameters,Past_Modules)),
    % Output parameters:update the parameters include class_plot
    append(Parameters,[class_plot,Class_plot],NewParameters),
    % Module:plot_prediction
    Newadd=[plot_prediction,NewParameters],
    append(Past_Modules,Newadd,Modules).

/*---------------------------------------------------------------*/
prediction_pca_plot(Parameters,Modules):-
    % Conditions: class_plot should be svm,weightedvoting or loocv
    get_value_variable(Parameters,class_plot,Class_plot),
    member(Class_plot,[svm,weightedvoting,loocv]),
    % Input: loocv,weightedVoting or svm_predictions
    (Class_plot=loocv,
    loocv(Parameters,Past_Modules);
    Class_plot=weightedvoting,
    weightedVoting(Parameters,Past_Modules);
    Class_plot=svm,
    svm_predictions(Parameters,Past_Modules)),
    % Output parameters:update the parameters include class_plot
    append(Parameters,[class_plot,Class_plot],NewParameters3),
    % Module:plot_prediction_pca
    Newadd=[plot_prediction_pca,NewParameters3],
    append(Past_Modules,Newadd,Modules).

