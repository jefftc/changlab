/*differential_expressed_gene_analysis.pl*/

%% differential_expressed_genes(+Parameters,-Modules)
differential_expressed_genes(Parameters,Modules):-
    % Conditions: pcl,logged,created, diff_expr = t_test or sam
    get_value_variable(Parameters,diff_expr,Diff_expr),
    member(Diff_expr,[t_test,sam]),
    member((Diff_expr,Module),[(t_test,t_test),(sam,sam)]),
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=pcl,
    get_value(NewParameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    get_value(NewParameters,status,created,Status),
    Status=created,
    % Input: signal_file and class_label_file
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_2),
    member(OldStatus,[given,jointed,splited,created]),
    set_value(NewParameters,status,OldStatus,NewParameters1),
    signal_file(NewParameters1,Past_Modules_1),
    append(Past_Modules_2,Past_Modules_1,Past_Modules),
    (Diff_expr = sam,
    get_options(Parameters,[sam_delta,sam_foldchange],[],Options);
    Diff_expr = t_test,
    Options = []),
    append(NewParameters,Options,NewParameters2),
    % Module: t_test or sam
    % Output parameters: update the output parameters to full_length
    append([diff_expr,Diff_expr],NewParameters2,Write_list),
    Newadd=[Module,Write_list],
    append(Past_Modules,Newadd,Modules).

