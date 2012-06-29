/*comparison_report.pl*/

comes_before(Modules,First,Second):-
    not(Second=none),
    nth0(N,Modules,First),
    nth0(N1,Modules,Second),
    N<N1;
    Second=none,
    member(First,Modules).
/*-------------------------------------------------------------------*/
make_batch_report(Parameters,Modules):-
    member((Options,First,Second),[([quantile,yes_quantile,dwd,no_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],quantile,none),
([quantile,yes_quantile,dwd,yes_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],quantile,dwd),
([quantile,yes_quantile,shiftscale,yes_shiftscale,combat,no_combat,bfrm,no_bfrm,dwd,no_dwd],quantile,shiftscale),
([bfrm,yes_bfrm,quantile,no_quantile,shiftscale,no_shiftscale,combat,no_combat,dwd,no_dwd],bfrm_normalize,none),
([quantile,yes_quantile,combat,yes_combat,shiftscale,no_shiftscale,dwd,no_dwd,bfrm,no_bfrm],quantile,combat)]),
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

     get_options(Parameters,[pca_gene_num],[],Options2),
     append(NewParameters,Options2,NewParameters2),
     append(NewParameters2,[objecttype,pca_plot_out],NewParameters3),
     Newadd=[pca_plot,NewParameters3],
     append(Modules1, Newadd, Modules2),

     pca_plot_in(NewParameters2,Modules3),

     Modules=[Modules1,Modules2,Modules3].
/*-------------------------------------------------------------------*/

make_diffgenes_report(Parameters,Modules):-
    get_value_variable(Parameters,diff_expr,Diff_expr),
    member(Diff_expr,[t_test,sam]),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=pcl,
    get_value(NewParameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(NewParameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,jointed,splited,created]),
    set_value(NewParameters,status,OldStatus,NewParameters1),
    append(NewParameters1,[diff_expr,Diff_expr],NewParameters2),
    differential_expressed_genes(NewParameters2,Modules1),
    append(NewParameters1,[filetype,signal_file],NewParameters3),
    signal_file(NewParameters1,Modules2),
    append(Modules2,[annot_file,NewParameters3],Modules3),
    Modules = [Modules1,Modules3].

/*-------------------------------------------------------------------*/
make_cluster_report(Parameters,Modules):-
    convert_cluster_parameters(Parameters,NewParameters1),
    cluster_file(NewParameters1,Past_Modules1),
    append(NewParameters1,[filetype,cluster_file],NewParameters2),
    append(Past_Modules1,[annot_file,NewParameters2],Past_Modules2),
    get_options(Parameters,[hm_width,hm_height,color],[],Options),
    append(NewParameters2,Options,NewParameters3),
    cluster_heatmap(NewParameters3,Past_Modules3),

    Modules = [Past_Modules2,Past_Modules3].
/*-------------------------------------------------------------------*/
make_classify_report(Parameters,Modules):-
    convert_parameters_classify(Parameters,NewParameters),
    svm_predictions(NewParameters,Past_Modules1),
    set_value(NewParameters,format,gct,NewParameters3),
    weightedVoting(NewParameters3,Past_Modules2),
    append(NewParameters,[classification,svm],NewParameters1),
    append(NewParameters3,[classification,weightedvoting],NewParameters2),
    loocv(NewParameters1,Past_Modules3),
    loocv(NewParameters2,Past_Modules4),
    Modules=[Past_Modules1,Past_Modules2,Past_Modules3,Past_Modules4].
/*-------------------------------------------------------------------*/
make_normalize_report(Parameters,Modules):-
    convert_parameters_file(Parameters,NewParameters),
    signal_file(NewParameters,Past_Modules1),
    append(NewParameters,[filetype,signal_file],NewParameters0),
    append(Past_Modules1,[annot_file,NewParameters0],Modules1),
    
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters,Options,NewParameters1),
    pca_plot_in(NewParameters1,Modules2),

    set_value(NewParameters,format,pcl,NewParameters2),
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters2,Options,NewParameters3),
    pca_plot_out(NewParameters3,Modules3),

    intensity_plot(NewParameters,Modules4),
    
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    (Preprocess=illumina, 
    convert_parameters_raw(Parameters,NewParameters4), 
    set_value(NewParameters4,is_logged,no_logged,NewParameters5),
    set_value(NewParameters5,format,gct,NewParameters6),
    set_value(NewParameters6,has_missing_value,unknown_missing,NewParameters7),
    control_file(NewParameters7,Modules8),
    biotin_plot(NewParameters,Modules5),
    actb_plot(NewParameters,Modules6),
    hyb_bar_plot(NewParameters,Modules7),
    Modules=[Modules1,Modules2,Modules3,Modules4,Modules8,Modules5,Modules6,Modules7];
    not(member(Preprocess,[illumina])),
    Modules=[Modules1,Modules2,Modules3,Modules4]).