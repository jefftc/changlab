/*comparison_report.pl*/

% Determine if the First is before the Second in Modules
comes_before(Modules,First,Second,A,B,C):-
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
    member((Options,First,Second,Label,Pca_label1,Pca_label2),[([quantile,yes_quantile,dwd,no_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],normalize_samples_with_quantile,none,'Normalizing_the_samples_with_quantile_method','Plotting_PCA_plot_for_samples_after_quantile_normalization','Plotting_PCA_plot_for_samples_before_quantile_normalization'),
([bfrm,yes_bfrm,quantile,no_quantile,shiftscale,no_shiftscale,combat,no_combat,dwd,no_dwd],normalize_samples_with_bfrm,none,'Normalizing_the_samples_with_bfrm_method','Plotting_PCA_plot_for_samples_after_bfrm_normalization',
'Plotting_PCA_plot_for_samples_before_bfrm_normalization'),
([quantile,yes_quantile,combat,yes_combat,shiftscale,no_shiftscale,dwd,no_dwd,bfrm,no_bfrm],normalize_samples_with_quantile,normalize_samples_with_combat,'Normalizing_the_samples_with_quantile_and_combat_method',
'Plotting_PCA_plot_for_samples_after_quantile_and_combat_normalization',
'Plotting_PCA_plot_for_samples_before_quantile_and_combat_normalization'),
([quantile,yes_quantile,dwd,yes_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],normalize_samples_with_quantile,normalize_samples_with_dwd,'Normalizing_the_samples_with_quantile_and_dwd_method',
'Plotting_PCA_plot_for_samples_after_quantile_and_dwd_normalization',
'Plotting_PCA_plot_for_samples_before_quantile_and_dwd_normalization'),
([quantile,yes_quantile,shiftscale,yes_shiftscale,combat,no_combat,bfrm,no_bfrm,dwd,no_dwd],normalize_samples_with_quantile,normalize_samples_with_shiftscale,'Normalizing_the_samples_with_quantile_and_shiftscale_method',
'Plotting_PCA_plot_for_samples_after_quantile_and_shiftscale_normalization',
'Plotting_PCA_plot_for_samples_before_quantile_and_shiftscale_normalization')]),
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
     comes_before(Modules1,First,Second,Label,Pca_label1,Pca_label2),
     append(Modules1,[pipeline_label,[label,Label]],Modules1_1),

     %Input2: pca_plot_out 
     get_options(Parameters,[pca_gene_num],[],Options2),
     append(NewParameters,Options2,NewParameters2),
     append(NewParameters2,[objecttype,pca_plot_out],NewParameters3),
     Newadd=[analyze_samples_pca,NewParameters3,plot_sample_pca,NewParameters3],
     append(Modules1, Newadd, Modules2),
     append(Modules2,[pipeline_label,[label,Pca_label1]],Modules2_1),

     % Input3: pca_plot_in
     pca_plot_in(NewParameters2,Modules3),
     append(Modules3,[pipeline_label,[label,Pca_label2]],Modules3_1),

     Modules=[Modules1_1,Modules2_1,Modules3_1].
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
    append(Modules1,[pipeline_label,[label,'Doing_differential_expressed_genes_analysis_with_t_test']],Modules1_1),

    % Input2: differential_expressed_genes with sam
    append(NewParameters1,[diff_expr,sam],NewParameters3),
    differential_expressed_genes(NewParameters3,Modules2),
    append(Modules2,[pipeline_label,[label,'Doing_differential_expressed_genes_analysis_with_sam']],Modules2_1),

    % Input3: cluster_heatmap with t_test_p gene_order and threshold 0.05
    set_value(NewParameters1,gene_order,t_test_p,NewParameters4),
    append(NewParameters4,[cluster_alg,no_cluster_alg,hm_width,50,hm_height,1],NewParameters5),
    cluster_heatmap(NewParameters5,Modules3),
    append(Modules3,[pipeline_label,[label,'Plotting_the_heatmap_of_gene_ordered_by_t_test_p_value']],Modules3_1),

    % Input4: gather with t_test_p
    gather(NewParameters4,Modules4),
    append(Modules4,[pipeline_label,[label,'Doing_gather_analysis_on_the_important_genes']],Modules4_1),

    % Input5: gsea
    set_value(NewParameters1,format,tdf,NewParameters6),
    gsea(NewParameters6,Modules5),
    append(Modules5,[pipeline_label,[label,'Doing_GSEA_on_the_signal_file']],Modules5_1),

    Modules = [Modules1_1,Modules2_1,Modules3_1,Modules4_1,Modules5_1].

/*-------------------------------------------------------------------*/
make_cluster_report(Parameters,Modules):-
    % Input1: cluster_file 
    convert_cluster_parameters(Parameters,NewParameters1),
    cluster_file(NewParameters1,Past_Modules1),
    append(NewParameters1,[filetype,cluster_file,annotate_type,all],NewParameters2),
    % Module: annot_probe
    % Output parameters:the full length parameters of cluster_file and [filetype,cluster_file]
    append(Past_Modules1,[annotate_probes,NewParameters2],Past_Modules2),
    append(Past_Modules2,[pipeline_label,[label,'Clustering_the_samples']],Past_Modules2_1),

    % Input2: cluster_heatmap
    get_options(Parameters,[hm_width,hm_height,color],[],Options),
    append(NewParameters2,Options,NewParameters3),
    cluster_heatmap(NewParameters3,Past_Modules3),
    append(Past_Modules3,[pipeline_label,[label,'Plotting_a_heatmap_of_the_samples']],Past_Modules3_1),

    Modules = [Past_Modules2_1,Past_Modules3_1].

/*-------------------------------------------------------------------*/
make_classify_report(Parameters,Modules):-
    convert_parameters_classify(Parameters,NewParameters),
    % Input1: test file after common genes aglin
    convert_parameters_svm(NewParameters,NewParameters1),
    %get_value(NewParameters1,testcontents,[],TestContents),
    %set_value(NewParameters1,contents,TestContents,NewParameters2),
    signal_file(NewParameters1,Past_Modules1),
    append(Past_Modules1,[pipeline_label,[label,'Preprocessing_the_signal_file']],Past_Modules1_1),

    % Input2: svm_predictions
    svm_predictions(NewParameters,Past_Modules2),
    append(Past_Modules2,[pipeline_label,[label,'Classifying_by_SVM']],Past_Modules2_1),

    % Input 3:prediction_plot with svm
    append(NewParameters,[class_plot,svm],NewParameters3),
    prediction_plot(NewParameters3,Past_Modules3),
    append(Past_Modules3,[pipeline_label,[label,'Plotting_prediction_results_with_SVM']],Past_Modules3_1),

    % Input4: loocv with svm
    append(NewParameters,[classification,svm],NewParameters4),
    loocv(NewParameters4,Past_Modules4),
    append(Past_Modules4,[pipeline_label,[label,'LOOCV_classifying_by_SVM']],Past_Modules4_1),

    % Input5:prediction_plot for loocv with svm
    append(NewParameters,[classification,svm,class_plot,loocv],NewParameters5),
    prediction_plot(NewParameters5,Past_Modules5),
    append(Past_Modules5,[pipeline_label,[label,'Plotting_prediction_results_with_LOOCV_SVM']],Past_Modules5_1),

   % Input 6: predication_pca_plot
    append(NewParameters,[class_plot,svm],NewParameters6),
    prediction_pca_plot(NewParameters6,Past_Modules6),
    append(Past_Modules6,[pipeline_label,[label,'Plot_pca_plot_colored_by_SVM_prediction_results']],Past_Modules6_1),

    set_value(NewParameters,format,gct,NewParameters7),
    % Input7: weightedVoting
    weightedVoting(NewParameters7,Past_Modules7),
    append(Past_Modules7,[pipeline_label,[label,'Classifying_by_weighted_voting']],Past_Modules7_1),

    % Input8: prediction_plot for weightedVoting
    append(NewParameters7,[class_plot,weightedvoting,classification,weightedvoting],NewParameters10),
    prediction_plot(NewParameters10,Past_Modules8),
    append(Past_Modules8,[pipeline_label,[label,'Plotting_prediction_results_with_weighted_voting']],Past_Modules8_1),

    % Input9:loocv with weightedvoting
    append(NewParameters7,[classification,weightedvoting],NewParameters8),
    loocv(NewParameters8,Past_Modules9),
    append(Past_Modules9,[pipeline_label,[label,'LOOCV_classifying_by_weighted_voting']],Past_Modules9_1),

    % Input 10:prediction_plot for loocv with weightedvoting
    append(NewParameters7,[classification,weightedvoting,class_plot,loocv],NewParameters9),
    prediction_plot(NewParameters9,Past_Modules10),
    append(Past_Modules10,[pipeline_label,[label,'Plotting_prediction_results_with_LOOCV_weighted_voting']],Past_Modules10_1),

    % Input11: predication_pca_plot
    prediction_pca_plot(NewParameters10,Past_Modules11),
    append(Past_Modules11,[pipeline_label,[label,'Plot_pca_plot_colored_by_weighted_voting_prediction_results']],Past_Modules11_1),

    Modules=[Past_Modules1_1,Past_Modules2_1,Past_Modules3_1,Past_Modules4_1,Past_Modules5_1,Past_Modules6_1,Past_Modules7_1,
Past_Modules8_1,Past_Modules9_1,Past_Modules10_1,Past_Modules11_1].
/*-------------------------------------------------------------------*/
make_normalize_report(Parameters,Modules):-
    % Input1: signal_file
    convert_parameters_file(Parameters,NewParameters),
    signal_file(NewParameters,Past_Modules1),
    % Module: annot_probes
    % Output parameters:full length parameters of signal_file and [filetype,signal_file]
    append(NewParameters,[filetype,signal_file,annotate_type,all],NewParameters0),
    append(Past_Modules1,[annotate_probes,NewParameters0],Modules1_1),
    append(Modules1_1,[pipeline_label,[label,'Preprocessing_the_signal_file']],Modules1),

    % Input2: pca_plot_in
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters,Options,NewParameters1),
    pca_plot_in(NewParameters1,Modules2_1),
    append(Modules2_1,[pipeline_label,[label,'Making_a_PCA_plot_of_the_samples_before_normalizeion']],Modules2),

    % Input3: pca_plot_out
    set_value(NewParameters,format,pcl,NewParameters2),
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters2,Options,NewParameters3),
    pca_plot_out(NewParameters3,Modules3_1),
    append(Modules3_1,[pipeline_label,[label,'Making_a_PCA_plot_of_the_samples_after_normalizeion']],Modules3),

    % Input4: intensity_plot
    intensity_plot(NewParameters,Modules4_1),
    append(Modules4_1,[pipeline_label,[label,'Making_a_boxplot_of_the_signal_intensities']],Modules4),

    % Input 5:actb_plot
    actb_plot(NewParameters,Modules5_1),
    append(Modules5_1,[pipeline_label,[label,'Plotting_the_ACTB_housekeeping_gene']],Modules5),

    % Conditions: if preprocess is illumina
    get_value(NewParameters,preprocess,unknown_preprocess,Preprocess),
    (Preprocess=illumina, 
    convert_parameters_raw(NewParameters,NewParameters4), 
    set_value(NewParameters4,is_logged,no_logged,NewParameters5),
    set_value(NewParameters5,format,gct,NewParameters6),
    set_value(NewParameters6,has_missing_value,unknown_missing,NewParameters7),
    % Input6: control_file
    control_file(NewParameters7,Modules6_1),
    append(Modules6_1,[pipeline_label,[label,'Plotting_the_control_probes']],Modules6),

    % Input7: biotin_plot
    biotin_plot(NewParameters,Modules7_1),
    append(Modules7_1,[pipeline_label,[label,'Plotting_the_biotin_signal']],Modules7),

    % Input8: housekeeping_plot
    housekeeping_plot(NewParameters,Modules8_1),
    append(Modules8_1,[pipeline_label,[label,'Plotting_the_Illumina_housekeeping_genes']],Modules8),

    % Input9: hyb_bar_plot
    hyb_bar_plot(NewParameters,Modules9_1),
    append(Modules9_1,[pipeline_label,[label,'Plotting_the_hyb_genes']],Modules9),

    Modules=[Modules1,Modules2,Modules3,Modules4,Modules5,Modules7,Modules8,Modules9,Modules6];
    not(member(Preprocess,[illumina])),
    control_plot(NewParameters,Modules9_1),
    append(Modules9_1,[pipeline_label,[label,'Extracting_the_signal_values_of_the_control_probes']],Modules9),
    Modules=[Modules1,Modules2,Modules3,Modules4,Modules5,Modules9]).

/*-------------------------------------------------------------------*/
make_heatmap_report(Parameters,Modules_1):-
    cluster_heatmap(Parameters,Modules),
    append(Modules,[pipeline_label,[label,'Plotting_a_heatmap_of_the_samples']],Modules_1).
/*-------------------------------------------------------------------*/
make_geneset_report(Parameters,Modules):-
    % Input1:geneset_analysis
    geneset_analysis(Parameters,Modules1),
    append(Modules1,[pipeline_label,[label,'Doing_geneset_analysis_on_the_samples']],Modules1_1),
    % Input2:geneset_plot
    geneset_plot(Parameters,Modules2),
    append(Modules2,[pipeline_label,[label,'Plotting_geneset_analysis_result']],Modules2_1),
    Modules=[Modules1_1,Modules2_1].
/*-------------------------------------------------------------------*/

make_call_variants_report(Parameters,Modules):-
    convert_vcf_parameters(Parameters,NewParameters1),
    % Input1: vcf_file with standard call variants
    append(NewParameters1,[reheader,standard],NewParameters2),
    vcf_file(NewParameters2,Past_Modules1),
    append(Past_Modules1,[pipeline_label,[label,'call_variants_by_GATK']],Past_Modules2),
    % Input2: call_variants_by_mpileup
    append(NewParameters1,[reheader,bcftool,filter,yes_filter],NewParameters3),
    vcf_file(NewParameters3,Past_Modules3),
    append(Past_Modules3,[pipeline_label,[label,'call_variants_by_mpileup']],Past_Modules4),
    Modules = [Past_Modules2,Past_Modules4].


/*-------------------------------------------------------------------------*/
convert_vcf_parameters(Parameters,NewParameters):-
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,read,single,Read),
    get_value(Parameters,ref,hg19,Ref),
    get_value(Parameters,recalibration,no_recalibration,Recalibration),
    NewParameters=[contents,Contents,read,Read,ref,Ref,recalibration,Recalibration].
