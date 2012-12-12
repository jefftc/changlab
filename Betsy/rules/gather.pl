/*gather enrichment*/

% rank genes by sample_ttest
gene_list_file(Parameters,Modules):-
    % Conditions: not a given gene_list_file,
    %    status=created,gene_order=t_test_p or t_test_fdr
    length(Parameters,N),
    N>2,
    convert_parameters_norm2(Parameters,NewParameters),
    get_value(NewParameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(Parameters,gene_order,no_order,Gene_Order),
    member(Gene_Order,[t_test_p,t_test_fdr]),
    % Input: signal_norm2 and class_label_file
    set_value(NewParameters,status,OldStatus,OldParameters),
    get_value(NewParameters,contents,[unknown],Contents),
    get_value(NewParameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    signal_norm2(OldParameters,Past_Modules_2),
    % Module: rank_genes_by_sample_ttest
    % Output parameters: update the parameters to full length
    get_value(Parameters,group_fc,no_group_fc,Group_fc),
    (Group_fc=yes_group_fc,
    append(NewParameters,[group_fc,Group_fc],NewParameters1),
    Newadd1=[filter_genes_by_fold_change_across_classes,NewParameters1];
    Group_fc=no_group_fc,
    NewParameters1=NewParameters,
    Newadd1=[]),
    append(NewParameters1,[gene_order,Gene_Order],NewParameters2),
    Newadd=[rank_genes_by_sample_ttest,NewParameters2],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd1, Modules1),
    append(Modules1,Newadd,Modules).
   
/*--------------------------------------------------------------------------*/
%rank genes by class_neighbors
gene_list_file(Parameters,Modules):-
    % Conditions: not a given gene_list_file,
    %    status=created,gene_order=by_class_neighbors
    length(Parameters,N),
    N>2,
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,status,created,Status),
    Status=created,
    member(OldStatus,[given,created,jointed,splited]),
    get_value(NewParameters,unique_genes,no_unique_genes,Unique_Genes),
    Unique_Genes = no_unique_genes,
    get_value(NewParameters,gene_order,no_order,Gene_Order),
    Gene_Order=by_class_neighbors,
    OldGene_Order=no_order,
    % Input: signal_file and class_label_file
    set_value(NewParameters,status,OldStatus,OldParameters1),
    set_value(OldParameters1,gene_order,OldGene_Order,OldParameters),
    get_value(NewParameters,contents,[unknown],Contents),
    get_value(NewParameters,preprocess,unknown_preprocess,Preprocess),
    class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
    signal_file(OldParameters,Past_Modules_2),
    % Module: rank_genes_by_class_neighbors
    % Output parameters: update the parameters to full length
    get_value(Parameters,group_fc,no_group_fc,Group_fc),
    (Group_fc=yes_group_fc,
    append(NewParameters,[group_fc,Group_fc],NewParameters1),
    Newadd1=[filter_genes_by_fold_change_across_classes,NewParameters1];
    Group_fc=no_group_fc,
    NewParameters1=NewParameters,
    Newadd1=[]),
    append(NewParameters1,[gene_order,Gene_Order],NewParameters2),
    Newadd=[rank_genes_by_class_neighbors,NewParameters2],
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(Past_Modules, Newadd1, Modules1),
    append(Modules1,Newadd,Modules).

/*--------------------------------------------------------------------------*/
% convert the parameters of gather into full length of version
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
   % Input: from the given gene_list_file
   (gene_list_file([contents,Contents],[]),
   % Output parameters:update the output parameters to full length
   Newadd=[contents,Contents,NewParameters],
   % Module: annotate_genes_with_gather
   Modules=[annotate_genes_with_gather,Newadd];
   % Conditions: not given a gene_list_file
   not(gene_list_file([contents,Contents],[])),
   convert_parameters_file(Parameters,NewParameters1),
   % Input: gene_list_file not given
   gene_list_file(NewParameters1,Past_Modules),
   append(NewParameters1,NewParameters,NewParameters2),
   % Module:annotate_genes_with_gather
   % Output parameters:update the output parameters to full length
   Newadd=[annotate_genes_with_gather,NewParameters2],
   append(Past_Modules,Newadd,Modules)).