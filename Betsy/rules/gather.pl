/*gather enrichment*/
% rank genes by sample_ttest

gene_list_file(Parameters,Modules):-
    length(Parameters,N),
    N>2,
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,contents,[unknown],Contents),
    get_value(NewParameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    get_value(NewParameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(NewParameters,gene_order,no_order,Gene_Order),
    member(Gene_Order,[t_test_p,t_test_fdr]),
    OldGene_Order=no_order,
    set_value(NewParameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(OldParameters,Past_Modules_2),
    Newadd=[rank_gene_by_sample_ttest,NewParameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).
/*--------------------------------------------------------------------------*/
%rank genes by class_neighbors
gene_list_file(Parameters,Modules):-
    length(Parameters,N),
    N>2,
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,contents,[unknown],Contents),
    get_value(NewParameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    get_value(NewParameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(NewParameters,gene_order,no_order,Gene_Order),
    Gene_Order=by_class_neighbors,
    OldGene_Order=no_order,
    set_value(NewParameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    signal_file(OldParameters,Past_Modules_2),
    Newadd=[class_neighbors,NewParameters],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd, Modules).

/*--------------------------------------------------------------------------*/
convert_parameters_gather(Parameters,NewParameters):-
   get_value(Parameters,network,no_network,Network),
   member(Network,[yes_network,no_network]),
   get_value(Parameters,homologs,yes_homologs,Homologs),
   member(Homologs,[yes_homologs,no_homologs]),
   get_value(Parameters,annot_type,gene_ontology,Annot_type),
   member(Annot_type,[kegg,gene_ontology,transfac]),
   NewParameters = [network,Network,homologs,Homologs,annot_type,Annot_type].
/*--------------------------------------------------------------------------*/
%run gather website
gather(Parameters,Modules):-
   convert_parameters_gather(Parameters,NewParameters),
   get_value(Parameters,contents,[unknown],Contents),
   (gene_list_file([contents,Contents],[]),
   Newadd=[contents,Contents,NewParameters],
   Modules=[run_gather,Newadd];
   not(gene_list_file([contents,Contents],[])),
   convert_parameters_file(Parameters,NewParameters1),
   gene_list_file(NewParameters1,Past_Modules),
   append(NewParameters1,NewParameters,NewParameters2),
   Newadd=[run_gather,NewParameters2],
   append(Past_Modules,Newadd,Modules)).