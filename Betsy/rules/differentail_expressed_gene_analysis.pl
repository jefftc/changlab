/*differential_expressed_gene_analysis.pl*/

%% differential_expressed_genes(+Parameters,-Modules)
differential_expressed_genes(Parameters,Modules):-
    get_value_variable(Parameters,diff_expr,Diff_expr),
    member(Diff_expr,[t_test,sam]),
    member((Diff_expr,Module),[(t_test,t_test),(sam,sam)]),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_2),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=pcl,
    get_value(NewParameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(NewParameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,jointed,splited,created]),
    set_value(NewParameters,status,OldStatus,NewParameters1),
    signal_file(NewParameters1,Past_Modules_1),
    append(Past_Modules_2,Past_Modules_1,Past_Modules),
    (Diff_expr = sam,
    get_options(Parameters,[sam_delta,sam_foldchange],[],Options);
    Diff_expr = t_test,
    Options = []),
    append(NewParameters,Options,NewParameters2),
    append([diff_expr,Diff_expr],NewParameters2,Write_list),
    Newadd=[Module,Write_list],
    append(Past_Modules,Newadd,Modules).

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
    signal_file(NewParameters1,Modules2),
    append(Modules2,[annot_file,NewParameters1],Modules3),
    Modules = [Modules1,Modules3].