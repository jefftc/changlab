/*comparison_report.pl*/

% Determine if the First is before the Second in Modules
comes_before(Modules,First,Second):-
    not(Second=none),
    nth0(N,Modules,First),
    nth0(N1,Modules,Second),
    N<N1;
    Second=none,
    member(First,Modules).

/*-------------------------------------------------------------------*/
% 
make_batch_report(Parameters,Modules):-
    % Conditions: desire batch methods and order
    member((Options,First,Second),[([quantile,yes_quantile,dwd,no_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],normalize_samples_with_quantile,none),
([bfrm,yes_bfrm,quantile,no_quantile,shiftscale,no_shiftscale,combat,no_combat,dwd,no_dwd],normalize_samples_with_bfrm,none),
([quantile,yes_quantile,combat,yes_combat,shiftscale,no_shiftscale,dwd,no_dwd,bfrm,no_bfrm],normalize_samples_with_quantile,normalize_samples_with_combat),
([quantile,yes_quantile,dwd,yes_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],normalize_samples_with_quantile,normalize_samples_with_dwd),
([quantile,yes_quantile,shiftscale,yes_shiftscale,combat,no_combat,bfrm,no_bfrm,dwd,no_dwd],normalize_samples_with_quantile,normalize_samples_with_shiftscale)]),
     %Input1: signal_norm1 with created and desired batch methods
     convert_parameters_clean_out(Parameters,NewParameters1),

     (First=bfrm_normalize,
     get_options(Parameters,[num_factors],[],Options3),
     append(NewParameters1,Options3,NewParameters4);
     not(First=bfrm_normalize),
     NewParameters4=NewParameters1),

     get_value(Parameters,status,created,Status),
     Status=created,
     append(NewParameters4,Options,NewParameters),
     signal_norm1(NewParameters,Modules1),
     comes_before(Modules1,First,Second),
   
     %Input2: pca_plot_out 
     get_options(Parameters,[pca_gene_num],[],Options2),
     append(NewParameters,Options2,NewParameters2),
     append(NewParameters2,[objecttype,pca_plot_out],NewParameters3),
     Newadd=[plot_sample_pca,NewParameters3],
     append(Modules1, Newadd, Modules2),

     % Input3: pca_plot_in
     pca_plot_in(NewParameters2,Modules3),

     Modules=[Modules1,Modules2,Modules3].
/*-------------------------------------------------------------------*/
make_diffgenes_report(Parameters,Modules):-
   % Conditions: Parameters has tdf,logged,no_order,created
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=tdf,
    get_value(NewParameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(NewParameters,gene_order,no_order,Gene_Order),
    Gene_Order=no_order,
    get_value(NewParameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,jointed,splited,created]),
    set_value(NewParameters,status,OldStatus,NewParameters1),
    % Input1: differential_expressed_genes with t_test
    append(NewParameters1,[diff_expr,t_test],NewParameters2),
    differential_expressed_genes(NewParameters2,Modules1),
    % Input2: differential_expressed_genes with sam
    append(NewParameters1,[diff_expr,sam],NewParameters3),
    differential_expressed_genes(NewParameters3,Modules2),
    % Input3: cluster_heatmap with t_test_p gene_order and threshold 0.05
    set_value(NewParameters1,gene_order,t_test_p,NewParameters4),
    append(NewParameters4,[cluster_alg,no_cluster_alg,hm_width,50,hm_height,1],NewParameters5),
    cluster_heatmap(NewParameters5,Modules3),
    % Input4: gather with t_test_p
    gather(NewParameters4,Modules4),
    % Input5: gsea
    set_value(NewParameters1,format,gct,NewParameters6),
    gsea(NewParameters6,Modules5),
    Modules = [Modules1,Modules2,Modules3,Modules4,Modules5].

/*-------------------------------------------------------------------*/
make_cluster_report(Parameters,Modules):-
    % Input1: cluster_file 
    convert_cluster_parameters(Parameters,NewParameters1),
    cluster_file(NewParameters1,Past_Modules1),
    append(NewParameters1,[filetype,cluster_file,annotate_type,all],NewParameters2),
    % Module: annot_probe
    % Output parameters:the full length parameters of cluster_file and [filetype,cluster_file]
    append(Past_Modules1,[annotate_probe,NewParameters2],Past_Modules2),

    % Input2: cluster_heatmap
    get_options(Parameters,[hm_width,hm_height,color],[],Options),
    append(NewParameters2,Options,NewParameters3),
    cluster_heatmap(NewParameters3,Past_Modules3),

    Modules = [Past_Modules2,Past_Modules3].

/*-------------------------------------------------------------------*/
/*make_classify_report(Parameters,Modules):-
    % Input1: svm_predictions
    convert_parameters_classify(Parameters,NewParameters),
    append(NewParameters,[class_plot,svm],NewParameters1),
    prediction_plot(NewParameters1,Past_Modules1),

    % Input2:loocv with svm
    append(NewParameters,[classification,svm,class_plot,loocv],NewParameters2),
    prediction_plot(NewParameters2,Past_Modules2),

    % Input3: weightedVoting
    set_value(NewParameters,format,gct,NewParameters3),
    append(NewParameters3,[class_plot,weightedvoting],NewParameters4),
    prediction_plot(NewParameters4,Past_Modules3),

    % loocv with weightedvoting
    append(NewParameters3,[classification,weightedvoting,class_plot,loocv],NewParameters5),
    prediction_plot(NewParameters5,Past_Modules4),

    Modules=[Past_Modules1,Past_Modules2,Past_Modules3,Past_Modules4].*/
/*-------------------------------------------------------------------*/
/*make_classify_report(Parameters,Modules):-
    convert_parameters_classify(Parameters,NewParameters),
    % Input0: test file after common genes aglin
    convert_parameters_svm(NewParameters,NewParameters1),
    %get_value(NewParameters1,testcontents,[],TestContents),
    %set_value(NewParameters1,contents,TestContents,NewParameters2),
    signal_file(NewParameters1,Past_Modules1),

    % Input2: loocv with svm
    append(NewParameters,[classification,svm],NewParameters3),
    loocv(NewParameters3,Past_Modules2),
    
    % Input2: svm_predictions
    svm_predictions(NewParameters,Past_Modules3),

   % Input 3: predication_pca_plot
    append(NewParameters,[class_plot,svm],NewParameters6),
    prediction_pca_plot(NewParameters6,Past_Modules4),

    set_value(NewParameters,format,gct,NewParameters4),
    % Input4:loocv with weightedvoting
    append(NewParameters4,[classification,weightedvoting],NewParameters5),
    loocv(NewParameters5,Past_Modules5),

    % Input5: weightedVoting
    weightedVoting(NewParameters4,Past_Modules6),
   
    % Input6: predication_pca_plot
    append(NewParameters4,[class_plot,weightedvoting],NewParameters7),
    prediction_pca_plot(NewParameters7,Past_Modules7),
    Modules=[Past_Modules1,Past_Modules2,Past_Modules3,Past_Modules4,Past_Modules5,Past_Modules6,Past_Modules7].*/
/*-------------------------------------------------------------------*/
make_classify_report(Parameters,Modules):-
    convert_parameters_classify(Parameters,NewParameters),
    % Input1: test file after common genes aglin
    convert_parameters_svm(NewParameters,NewParameters1),
    %get_value(NewParameters1,testcontents,[],TestContents),
    %set_value(NewParameters1,contents,TestContents,NewParameters2),
    signal_file(NewParameters1,Past_Modules1),

    % Input2: svm_predictions
    svm_predictions(NewParameters,Past_Modules2),

    % Input 3:prediction_plot with svm
    append(NewParameters,[class_plot,svm],NewParameters3),
    prediction_plot(NewParameters3,Past_Modules3),

    % Input4: loocv with svm
    append(NewParameters,[classification,svm],NewParameters4),
    loocv(NewParameters4,Past_Modules4),
    
    % Input5:prediction_plot for loocv with svm
    append(NewParameters,[classification,svm,class_plot,loocv],NewParameters5),
    prediction_plot(NewParameters5,Past_Modules5),

   % Input 6: predication_pca_plot
    append(NewParameters,[class_plot,svm],NewParameters6),
    prediction_pca_plot(NewParameters6,Past_Modules6),

    set_value(NewParameters,format,gct,NewParameters7),
    % Input7:loocv with weightedvoting
    append(NewParameters7,[classification,weightedvoting],NewParameters8),
    loocv(NewParameters8,Past_Modules7),

    % Input 8:prediction_plot for loocv with weightedvoting
    append(NewParameters7,[classification,weightedvoting,class_plot,loocv],NewParameters9),
    prediction_plot(NewParameters9,Past_Modules8),

    % Input9: weightedVoting
    weightedVoting(NewParameters7,Past_Modules9),
   
    % Input10: prediction_plot for weightedVoting
    append(NewParameters7,[class_plot,weightedvoting,classification,weightedvoting],NewParameters10),
    prediction_plot(NewParameters10,Past_Modules10),

    % Input11: predication_pca_plot
    prediction_pca_plot(NewParameters10,Past_Modules11),

    Modules=[Past_Modules1,Past_Modules2,Past_Modules3,Past_Modules4,Past_Modules5,Past_Modules6,Past_Modules7,
Past_Modules8,Past_Modules9,Past_Modules10,Past_Modules11].
/*-------------------------------------------------------------------*/
make_normalize_report(Parameters,Modules):-
    % Input1: signal_file
    convert_parameters_file(Parameters,NewParameters),
    signal_file(NewParameters,Past_Modules1),
    % Module: annot_probes
    % Output parameters:full length parameters of signal_file and [filetype,signal_file]
    append(NewParameters,[filetype,signal_file,annotate_type,all],NewParameters0),
    append(Past_Modules1,[annotate_probes,NewParameters0],Modules1),
    % Input2: pca_plot_in
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters,Options,NewParameters1),
    pca_plot_in(NewParameters1,Modules2),
    % Input3: pca_plot_out
    set_value(NewParameters,format,pcl,NewParameters2),
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters2,Options,NewParameters3),
    pca_plot_out(NewParameters3,Modules3),
    % Input4: intensity_plot
    intensity_plot(NewParameters,Modules4),
    % Input 5:actb_plot
    actb_plot(NewParameters,Modules5),
    % Conditions: if preprocess is illumina
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    (Preprocess=illumina, 
    convert_parameters_raw(Parameters,NewParameters4), 
    set_value(NewParameters4,is_logged,no_logged,NewParameters5),
    set_value(NewParameters5,format,gct,NewParameters6),
    set_value(NewParameters6,has_missing_value,unknown_missing,NewParameters7),
    % Input6: control_file
    control_file(NewParameters7,Modules6),
    % Input7: biotin_plot
    biotin_plot(NewParameters,Modules7),
    % Input8: housekeeping_plot
    housekeeping_plot(NewParameters,Modules8),
    % Input9: hyb_bar_plot
    hyb_bar_plot(NewParameters,Modules9),
    Modules=[Modules1,Modules2,Modules3,Modules4,Modules5,Modules7,Modules8,Modules9,Modules6];
    not(member(Preprocess,[illumina])),
    control_plot(NewParameters,Modules9),
    Modules=[Modules1,Modules2,Modules3,Modules4,Modules5,Modules9]).

/*-------------------------------------------------------------------*/
make_heatmap_report(Parameters,Modules):-
    cluster_heatmap(Parameters,Modules).

/*-------------------------------------------------------------------*/
make_geneset_report(Parameters,Modules):-
    % Input1:geneset_analysis
    geneset_analysis(Parameters,Modules1),
    % Input2:geneset_plot
    geneset_plot(Parameters,Modules2),
    Modules=[Modules1,Modules2].
