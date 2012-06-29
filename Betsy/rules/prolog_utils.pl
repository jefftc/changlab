%prolog_utils.pl

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

use_module(library(lists)).

get_value(Parameters,Key,Default,Value):-
   member(Key,Parameters),
   nth0(N,Parameters,Key),
   X is N+1,
   nth0(X,Parameters,Value);
   not(member(Key,Parameters)),
   Value=Default.

/*-------------------------------------------------------------------------*/
get_value_variable(Parameters,Key,Value):-
   member(Key,Parameters),
   nth0(N,Parameters,Key),
   X is N+1,
   nth0(X,Parameters,Value);
   not(member(Key,Parameters)),
   Value=R.
/*-------------------------------------------------------------------------*/
replace([_|T],0,X,[X|T]). 
replace([H|T],I,X,[H|R]):-I1 is I-1, replace(T,I1,X,R).

/*-------------------------------------------------------------------------*/
set_value(Parameters,Key,Value,NewParameters):-
    member(Key,Parameters),
    nth0(N,Parameters,Key),
    replace(Parameters,N+1,Value,NewParameters).
/*-------------------------------------------------------------------------*/
remove_at(X,[X|Xs],1,Xs).
remove_at(X,[Y|Xs],K,[Y|Ys]) :- K > 1, 
   K1 is K - 1, remove_at(X,Xs,K1,Ys).

/*-------------------------------------------------------------------------*/
convert_parameters_raw(Parameters,NewParameters):-
   
    get_value(Parameters,contents,[unknown],Contents),
    NewParameters0=[contents,Contents],
 
    get_value(Parameters,format,unknown_format,Format),
    member(Format,[tdf,jeffs,res,gct,pcl,not_pcl,unknown_format,xls,not_xls]),
    append(NewParameters0,[format,Format],NewParameters1),
     
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    member(Preprocess,[rma,mas5,loess,unknown_preprocess,illumina,agilent]),
    append(NewParameters1,[preprocess,Preprocess],NewParameters2),
    
    get_value(Parameters,is_logged,unknown_logged,Is_Logged),
    member(Is_Logged,[no_logged,logged,unknown_logged]),
    append(NewParameters2,[is_logged,Is_Logged],NewParameters3),

    get_value(Parameters,has_missing_value,unknown_missing,Has_Missing_Value),
    member(Has_Missing_Value,[no_missing,median_fill,zero_fill,unknown_missing]),
    append(NewParameters3,[has_missing_value,Has_Missing_Value],NewParameters4),

    get_value(Parameters,filter,0,Filter),
    not(atom(Filter)),
    append(NewParameters4,[filter,Filter],NewParameters5),
    
    %get_value(Parameters,filter_fc,no_filter_fc,Filter_fc),
    %(Filter_fc=no_filter_fc;not(atom(Filter_fc))),
    %append(NewParameters5,[filter_fc,Filter_fc],NewParameters6),
    
    get_value(Parameters,predataset,no_predataset,Predataset),
    member(Predataset,[no_predataset,yes_predataset]),
    append(NewParameters5,[predataset,Predataset],NewParameters6),

    get_value(Parameters,status,created,Status),
    member(Status,[given,created,jointed,splited]),
    append(NewParameters6,[status,Status],NewParameters7),

    (member(Preprocess,[illumina]),
    get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip,ill_custom_manifest],[],Options),
    append(NewParameters7,Options,NewParameters8);
    not(member(Preprocess,[illumina])),
    NewParameters8=NewParameters7),
    
    get_value(Parameters,rename_sample,no_rename,Rename_sample),
    member(Rename_sample,[no_rename,yes_rename]),
    append(NewParameters8,[rename_sample,Rename_sample],NewParameters).

/*-------------------------------------------------------------------------*/
convert_parameters_clean_out(Parameters,NewParameters):-
    get_value(Parameters,contents,[unknown],Contents),
    append([],[contents,Contents],NewParameters0),

    append(NewParameters0,[format,pcl],NewParameters1),
    get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
    member(Preprocess,[rma,mas5,loess,unknown_preprocess,illumina,agilent]),
    append(NewParameters1,[preprocess,Preprocess],NewParameters2),

    append(NewParameters2,[is_logged,logged],NewParameters3),

    get_value_variable(Parameters,has_missing_value,Has_Missing_Value),
    member(Has_Missing_Value,[median_fill,zero_fill,no_missing]),
    append(NewParameters3,[has_missing_value,Has_Missing_Value],NewParameters4),

    get_value(Parameters,filter,0,Filter),
    not(atom(Filter)),
    append(NewParameters4,[filter,Filter],NewParameters5),

    get_value(Parameters,predataset,no_predataset,Predataset),
    member(Predataset,[yes_predataset,no_predataset]),
    append(NewParameters5,[predataset,Predataset],NewParameters6),

    %get_value(Parameters,filter_fc,no_filter_fc,Filter_fc),
    %(Filter_fc=no_filter_fc;not(atom(Filter_fc))),
    %append(NewParameters5,[filter_fc,Filter_fc],NewParameters6),

    get_value_variable(Parameters,status,Status),
    member(Status,[given,created,jointed,splited]),
    append(NewParameters6,[status,Status],NewParameters7),

    (member(Preprocess,[illumina]),
    get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip],[],Options),
    append(NewParameters7,Options,NewParameters8);
    not(member(Preprocess,[illumina])),
    NewParameters8=NewParameters7),

    get_value(Parameters,rename_sample,no_rename,Rename_sample),
    member(Rename_sample,[yes_rename,no_rename]),
    append(NewParameters8,[rename_sample,Rename_sample],NewParameters).
/*-----------------------------------------------*/
convert_parameters_variable_raw(Parameters,NewParameters):-
    get_value(Parameters,contents,[unknown],Contents),
    append([],[contents,Contents],NewParameters0),

    get_value_variable(Parameters,format,Format),
    member(Format,[tdf,jeffs,res,gct,pcl,unknown_format,xls,not_xls]),
    append(NewParameters0,[format,Format],NewParameters1),
    
    get_value_variable(Parameters,preprocess,Preprocess),
    member(Preprocess,[rma,mas5,loess,agilent,illumina,unknown_preprocess]),
    append(NewParameters1,[preprocess,Preprocess],NewParameters2),

    get_value_variable(Parameters,is_logged,Is_Logged),
    member(Is_Logged,[logged,no_logged,unknown_logged]),
    append(NewParameters2,[is_logged,Is_Logged],NewParameters3),
    
    get_value_variable(Parameters,has_missing_value,Has_Missing_Value),
    member(Has_Missing_Value,[median_fill,zero_fill,no_missing,unknown_missing]),
    append(NewParameters3,[has_missing_value,Has_Missing_Value],NewParameters4),

    get_value(Parameters,filter,0,Filter),
    (member(Filter,[0,25]);not(Filter=25),not(Filter=0),not(atom(Filter))),
    append(NewParameters4,[filter,Filter],NewParameters5),
    
    %get_value_variable(Parameters,filter_fc,Filter_fc),
    %(member(Filter_fc,[no_filter_fc,2]);not(Filter_fc=2),not(atom(Filter_fc))),
    %append(NewParameters5,[filter_fc,Filter_fc],NewParameters6),

    get_value_variable(Parameters,predataset,Predataset),
    member(Predataset,[yes_predataset,no_predataset]),
    append(NewParameters5,[predataset,Predataset],NewParameters6),

    get_value_variable(Parameters,status,Status),
    member(Status,[given,created,jointed,splited]),
    append(NewParameters6,[status,Status],NewParameters7),

    (member(Preprocess,[illumina]),
    get_options(Parameters,[ill_manifest,ill_chip,ill_bg_mode,ill_coll_mode,ill_clm,ill_custom_chip],[],Options),
    append(NewParameters7,Options,NewParameters8);
    not(member(Preprocess,[illumina])),
    NewParameters8=NewParameters7),
    
    get_value_variable(Parameters,rename_sample,Rename_sample),
    member(Rename_sample,[yes_rename,no_rename]),
    append(NewParameters8,[rename_sample,Rename_sample],NewParameters).

/*-------------------------------------------------------------------------*/
%get the parameters list for signal_raw
get_desire_parameters_raw(Parameters,NewParameters):-

    get_value(Parameters,contents,[unknown],Contents),
    append([],[contents,Contents],NewParameters0),

    (member(format,Parameters),
    get_value_variable(Parameters,format,Format),
    append(NewParameters0,[format,Format],NewParameters1);
    not(member(format,Parameters)),
    NewParameters1=NewParameters0),
 
    (member(preprocess,Parameters),
    get_value_variable(Parameters,preprocess,Preprocess),
    append(NewParameters1,[preprocess,Preprocess],NewParameters2);
    not(member(preprocess,Parameters)),
    NewParameters2=NewParameters1),

    (member(is_logged,Parameters),
    get_value_variable(Parameters,is_logged,Is_Logged),
    append(NewParameters2,[is_logged,Is_Logged],NewParameters3);
    not(member(is_logged,Parameters)),
    NewParameters3=NewParameters2),

    (member(has_missing_value,Parameters),
    get_value_variable(Parameters,has_missing_value,Has_Missing_Value),
    append(NewParameters3,[has_missing_value,Has_Missing_Value],NewParameters4);
    not(member(has_missing_value,Parameters)),
    NewParameters4=NewParameters3),

    (member(filter,Parameters),
    get_value_variable(Parameters,filter,Filter),
    append(NewParameters4,[filter,Filter],NewParameters5);
    not(member(filter,Parameters)),
    NewParameters5=NewParameters4),
    
    %(member(filter_fc,Parameters),
    %get_value_variable(Parameters,filter_fc,Filter_fc),
    %append(NewParameters5,[filter_fc,Filter_fc],NewParameters6);
    %not(member(filter_fc,Parameters)),
    %NewParameters6=NewParameters5),

    (member(predataset,Parameters),
    get_value_variable(Parameters,predataset,Predataset),
    append(NewParameters5,[predataset,Predataset],NewParameters6);
    not(member(predataset,Parameters)),
    NewParameters6=NewParameters5),

    (member(status,Parameters),
    get_value_variable(Parameters,status,Status),
    member(Status,[created,jointed,splited,given]),
    append(NewParameters6,[status,Status],NewParameters7);
    not(member(status,Parameters)),
    NewParameters7=NewParameters6),
   
    (member(rename_sample,Parameters),
    get_value_variable(Parameters,rename_sample,Rename_sample),
    member(Rename_sample,[yes_rename,no_rename]),
    append(NewParameters7,[rename_sample,Rename_sample],NewParameters);
    not(member(rename_sample,Parameters)),
    NewParameters=NewParameters7).

/*-------------------------------------------------------------------------*/
convert_parameters_norm1(Parameters,NewParameters):-
    convert_parameters_variable_raw(Parameters,NewParameters1),

    get_value_variable(Parameters,quantile,Quantile),
    member(Quantile,[yes_quantile,no_quantile]),
    append(NewParameters1,[quantile,Quantile],NewParameters2),

    get_value_variable(Parameters,combat,Combat),
    member(Combat,[yes_combat,no_combat]),
    append(NewParameters2,[combat,Combat],NewParameters3),
     
    get_value_variable(Parameters,shiftscale,Shiftscale),
    member(Shiftscale,[yes_shiftscale,no_shiftscale]),
    append(NewParameters3,[shiftscale,Shiftscale],NewParameters4),

    get_value_variable(Parameters,dwd,Dwd),
    member(Dwd,[yes_dwd,no_dwd]),
    append(NewParameters4,[dwd,Dwd],NewParameters5),
 
    get_value_variable(Parameters,bfrm,Bfrm),
    member(Bfrm,[yes_bfrm,no_bfrm]),
    append(NewParameters5,[bfrm,Bfrm],NewParameters6),
    (Bfrm=yes_bfrm,
    get_options(Parameters,[num_factors],[],Options),
    append(NewParameters6,Options,NewParameters);
    Bfrm=no_bfrm,
    NewParameters=NewParameters6).

/*-------------------------------------------------------------------------*/

%get the parameters list for signal_norm1
get_desire_parameters_norm1(Parameters,NewParameters):-
    get_desire_parameters_raw(Parameters,NewParameters1),
   
    (member(quantile,Parameters),
    get_value_variable(Parameters,quantile,Quantile),
    append(NewParameters1,[quantile,Quantile],NewParameters2);
    not(member(quantile,Parameters)),
    NewParameters2=NewParameters1),

    (member(combat,Parameters),
    get_value_variable(Parameters,combat,Combat),
    append(NewParameters2,[combat,Combat],NewParameters3);
    not(member(combat,Parameters)),
    NewParameters3=NewParameters2),

    (member(shiftscale,Parameters),
    get_value_variable(Parameters,shiftscale,Shiftscale),
    append(NewParameters3,[shiftscale,Shiftscale],NewParameters4);
    not(member(shiftscale,Parameters)),
    NewParameters4=NewParameters3),

    (member(dwd,Parameters),
    get_value_variable(Parameters,dwd,Dwd),
    append(NewParameters4,[dwd,Dwd],NewParameters5);
    not(member(dwd,Parameters)),
    NewParameters5=NewParameters4),

    (member(bfrm,Parameters),
    get_value_variable(Parameters,bfrm,Bfrm),
    append(NewParameters5,[bfrm,Bfrm],NewParameters);
    not(member(bfrm,Parameters)),
    NewParameters=NewParameters5).


 /*-------------------------------------------------------------------------*/
convert_parameters_norm2(Parameters,NewParameters):-
    convert_parameters_norm1(Parameters,NewParameters1),
    
    get_value_variable(Parameters,gene_center,Gene_Center),
    member(Gene_Center,[mean,median,no_gene_center]),
    append(NewParameters1,[gene_center,Gene_Center],NewParameters2),
   
    get_value_variable(Parameters,gene_normalize,Gene_Normalize),
    member(Gene_Normalize,[variance,sum_of_squares,no_gene_normalize]),
    append(NewParameters2,[gene_normalize,Gene_Normalize],NewParameters).
 /*-------------------------------------------------------------------------*/
%get the parameters list for signal_norm2
get_desire_parameters_norm2(Parameters,NewParameters):-
    get_desire_parameters_norm1(Parameters,NewParameters1),

    (member(gene_center,Parameters),
    get_value_variable(Parameters,gene_center,Gene_Center),
    append(NewParameters1,[gene_center,Gene_Center],NewParameters2);
    not(member(gene_center,Parameters)),
    NewParameters2=NewParameters1),

    (member(gene_normalize,Parameters),
    get_value_variable(Parameters,gene_normalize,Gene_Normalize),
    append(NewParameters2,[gene_normalize,Gene_Normalize],NewParameters);
    not(member(gene_normalize,Parameters)),
    NewParameters=NewParameters2).
 /*-------------------------------------------------------------------------*/
convert_parameters_file(Parameters,NewParameters):-
    convert_parameters_norm2(Parameters,NewParameters1),
    get_value_variable(Parameters,gene_order,Gene_Order),
    member(Gene_Order,[t_test_p,t_test_fdr,by_gene_list,by_class_neighbors,no_order]),
    append(NewParameters1,[gene_order,Gene_Order],NewParameters2),
    (Gene_Order=by_class_neighbors,
    get_options(Parameters,[cn_num_neighbors,cn_num_perm,cn_user_pval,cn_mean_or_median,cn_ttest_or_snr,cn_filter_data,cn_min_threshold,cn_max_threshold,cn_min_folddiff,cn_abs_diff],[],Options),
    append(NewParameters2,Options,NewParameters3);
    member(Gene_Order,[t_test_p,t_test_fdr]),
    get_options(Parameters,[gene_select_threshold],[],Options1),
    append(NewParameters2,Options1,NewParameters);
    not(member(Gene_Order,[by_class_neighbors,t_test_p,t_test_fdr])),
    NewParameters3=NewParameters2),

    get_value(Parameters,platform,unknown_platform,Platform),
    member(Platform,["HG_U133A",unknown_platform]),
    append(NewParameters3,[platform,Platform],NewParameters4),
    
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    member(Unique_Genes,[no_unique_genes,average_genes,high_var,first_gene]),
    append(NewParameters4,[unique_genes,Unique_Genes],NewParameters).
    

 /*-------------------------------------------------------------------------*/
%get the parameters list for signal_file
get_desire_parameters_file(Parameters,NewParameters):-
    get_desire_parameters_norm2(Parameters,NewParameters1),
    (member(gene_order,Parameters),
    get_value_variable(Parameters,gene_order,Gene_Order),
    append(NewParameters1,[gene_order,Gene_Order],NewParameters2);
    not(member(gene_order,Parameters)),
    NewParameters2=NewParameters1),
    get_value(Parameters,platform,unknown_platform,Platform),
    append(NewParameters2,[platform,Platform],NewParameters3),
    get_value(Parameters,unique_genes,no_unique_genes,Unique_Genes),
    append(NewParameters3,[unique_genes,Unique_Genes],NewParameters).
 /*-------------------------------------------------------------------------*/
convert_parameters_svm(Parameters,NewParameters):-
    convert_parameters_file(Parameters,NewParameters1),
    get_value(NewParameters1,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,svm_kernel,linear,Svm_Kernel),
    member(Svm_Kernel,[linear,polynomial,rbf,sigmoid,precomputed_kernel]),
    append(NewParameters1,[svm_kernel,Svm_Kernel],NewParameters2),
    (member(traincontents,Parameters),
    member(testcontents,Parameters),
    get_value_variable(Parameters,traincontents,TrainContents),
    get_value_variable(Parameters,testcontents,TestContents),
    append(NewParameters2,[traincontents,TrainContents,testcontents,TestContents],NewParameters);
    not(member(traincontents,Parameters)),
    not(member(testcontents,Parameters)),
    NewParameters=NewParameters2).
 /*-------------------------------------------------------------------------*/

convert_parameters_classify(Parameters,NewParameters):-
    convert_parameters_file(Parameters,NewParameters1),
    get_value(NewParameters1,status,created,Status),
    Status=created,
    get_value(NewParameters1,format,unknown_format,Format),
    Format=pcl,
    get_value(Parameters,svm_kernel,linear,Svm_Kernel),
    member(Svm_Kernel,[linear,polynomial,rbf,sigmoid,precomputed_kernel]),
    append(NewParameters1,[svm_kernel,Svm_Kernel],NewParameters2),
    get_options(Parameters,[wv_num_features,wv_minstd,wv_feature_stat],[],Options),
    append(NewParameters2,Options,NewParameters3),
    get_value_variable(Parameters,traincontents,TrainContents),
    get_value_variable(Parameters,testcontents,TestContents),
    append(NewParameters3,[traincontents,TrainContents,testcontents,TestContents],NewParameters).
    

/*----------------------------------------------------------------------*/
/*Whether [X|Y] is a subset of Z*/
is_subset([X|Y], Z) :- member(X, Z), is_subset(Y, Z).
is_subset([], _).

/*make a subset and put it into Z.*/
make_subset([_|Y],Z):- make_subset(Y,Z).
make_subset([X|Y],Z):- make_subset(Y,W),append([X],W,Z).
make_subset([],[]).

/*make a proper subset and put it into Z.*/
make_psubset(X,Z):- make_subset(X,Z), not(perm(X,Z)).

/*put the intersection of the first two into the Z.*/
intersection([X|Y],M,[X|Z]):- member(X,M),intersection(Y,M,Z).
intersection([X|Y],M,Z):- \+member(X,M),intersection(Y,M,Z).
intersection([],_,[]).

/*Atoms that are in list [X|Y] but not list M. */
difference([X|Y],M,Z):- member(X,M),difference(Y,M,Z).
difference([X|Y],M,[X|Z]):- \+member(X,M),difference(Y,M,Z).
difference([],_,[]).

takeout(X,[X|R],R).
takeout(X,[F|R],[F|S]):- takeout(X,R,S).
perm([X|Y],Z):- perm(Y,W),takeout(X,Z,W).
perm([],[]).
/*-------------------------------------------------------------------------*/    
get_length(A,B):-
    A=n_raw, B=18;
    A=n_norm1,B=28;
    A=n_norm2,B=32;
    A=n_file,B=38.
/*-------------------------------------------------------------------------*/  
get_options(Parameters,Keys,S,Options):-
   member(Key,Keys),
   (member(Key,Parameters),
   nth0(N,Parameters,Key),
   X is N+1,
   nth0(X,Parameters,Value),
   append(S,[Key,Value],Option);
   not(member(Key,Parameters)),
   append(S,[],Option)),
   remove_at(Key,Keys,1,Subkeys),
   get_options(Parameters,Subkeys,Option,Options).

get_options(Parameters,[],S,Options):-
   Options=S.
       
    

