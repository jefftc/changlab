/*pca_analysis*/
pca_analysis(Parameters,Modules):-
   get_value(Parameters,objecttype,pca_plot_in,Objecttype),
   (Objecttype=pca_plot_out,
    convert_parameters_file(Parameters,NewParameters1),
    set_value(NewParameters1,format,tdf,NewParameters2),
    signal_file(NewParameters2,Past_Modules);
    Objecttype=pca_plot_in,
    convert_parameters_clean_out(Parameters,NewParameters2),
    signal_clean(NewParameters2,Past_Modules);
    Objecttype=prediction,
    convert_parameters_file(Parameters,NewParameters1),
    get_value(Parameters,testcontents,[unknown],TestContents),
    set_value(NewParameters1,contents,TestContents,NewParameters2),
    signal_file(NewParameters2,Past_Modules)),
    get_options(Parameters,[pca_gene_num],[],Options),
    append(NewParameters2,Options,NewParameters3),
    append(NewParameters3,[objecttype,Objecttype],NewParameters),
    % Module:analyze_pca
    % Output parameters:full length parameters of signal_file
    Newadd=[analyze_pca,NewParameters],
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
pca_plot_out(Parameters,Modules):-
    convert_parameters_file(Parameters,NewParameters1),
    append(NewParameters1,[objecttype,pca_plot_out],NewParameters),
    pca_analysis(NewParameters,Past_Modules_2),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    % Module:plot_pca
    % Output parameters: full length parameters of signal_file 
    Newadd=[plot_pca,NewParameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*-------------------------------------------------------------------------*/
pca_plot_in(Parameters,Modules):-
    convert_parameters_clean_out(Parameters,NewParameters1),
    append(NewParameters1,[objecttype,pca_plot_in],NewParameters),
    pca_analysis(NewParameters,Past_Modules_2),
    get_value(Parameters,contents,[unknown],Contents),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    % Module:plot_pca
    % Output parameters:full length parameters of signal_clean
    Newadd=[plot_pca,NewParameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).