/*comparison_report.pl*/

normalization_comparison_report(Parameters,Modules):-
    member((Options,First,Second),[([quantile,yes_quantile,dwd,no_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],quantile,none),
([quantile,yes_quantile,dwd,yes_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],quantile,dwd),
([quantile,yes_quantile,shiftscale,yes_shiftscale,combat,no_combat,bfrm,no_bfrm,dwd,no_dwd],quantile,shiftscale),
([bfrm,yes_bfrm,quantile,no_quantile,shiftscale,no_shiftscale,combat,no_combat,dwd,no_dwd],bfrm_normalize,none),
([quantile,yes_quantile,combat,yes_combat,shiftscale,no_shiftscale,dwd,no_dwd,bfrm,no_bfrm],quantile,combat)]),
     convert_parameters_clean_out(Parameters,NewParameters1),
     get_value(Parameters,status,created,Status),
     Status=created,
     append(NewParameters1,Options,NewParameters),
     signal_norm1(NewParameters,Past_Modules),
     comes_before(Past_Modules,First,Second),
     Modules=Past_Modules.
 

pca_plot_normalization(Parameters,Modules):-
    member((Options,First,Second),[([quantile,yes_quantile,dwd,no_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],quantile,none),
([quantile,yes_quantile,dwd,yes_dwd,combat,no_combat,bfrm,no_bfrm,shiftscale,no_shiftscale],quantile,dwd),
([quantile,yes_quantile,shiftscale,yes_shiftscale,combat,no_combat,bfrm,no_bfrm,dwd,no_dwd],quantile,shiftscale),
([bfrm,yes_bfrm,quantile,no_quantile,shiftscale,no_shiftscale,combat,no_combat,dwd,no_dwd],bfrm_normalize,none),
([quantile,yes_quantile,combat,yes_combat,shiftscale,no_shiftscale,dwd,no_dwd,bfrm,no_bfrm],quantile,combat)]),
     convert_parameters_clean_out(Parameters,NewParameters1),
     get_value(Parameters,status,created,Status),
     Status=created,
     append(NewParameters1,Options,NewParameters),
     signal_norm1(NewParameters,Past_Modules),
     comes_before(Past_Modules,First,Second),
     get_options(Parameters,[pca_gene_num],[],Options2),
     append(NewParameters,Options2,NewParameters2),
     Newadd=[pca_plot,NewParameters2],
     append(Past_Modules, Newadd, Modules).

comes_before(Modules,First,Second):-
    not(Second=none),
    nth0(N,Modules,First),
    nth0(N1,Modules,Second),
    N<N1;
    Second=none,
    member(First,Modules).
/*-------------------------------------------------------------------*/

