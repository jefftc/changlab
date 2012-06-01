/*differential_expressed_gene_analysis.pl*/

%% differential_expressed_genes(+Parameters,-Modules)
differential_expressed_genes(Parameters,Modules):-
    get_value(Parameters,annot,no_annot,Annot),
    Annot=no_annot,
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


differential_expressed_genes(Parameters,Modules):-
    get_value(Parameters,annot,no_annot,Annot),
    Annot=yes_annot,
    convert_parameters_file(Parameters,NewParameters1),
    get_value_variable(Parameters,diff_expr,Diff_expr),
    member(Diff_expr,[t_test,sam]),
    append(NewParameters1,[annot,no_annot,diff_expr,Diff_expr],NewParameters2),
    differential_expressed_genes(NewParameters2,Past_Modules),
    set_value(NewParameters2,annot,yes_annot,NewParameters),
    Newadd=[annot_file,NewParameters],
    append(Past_Modules, Newadd, Modules).