/*signature_analysis*/

gsea(Parameters,Modules):-
   convert_parameters_file(Parameters,NewParameters),
   get_value(NewParameters,contents,[unknown],Contents),
   get_value(NewParameters,preprocess,unknown_preprocess,Preprocess),
   class_label_file([contents,Contents,preprocess,Preprocess,status,_],Past_Modules_1),
   get_value(NewParameters,format,unknown_format,Format),
   Format=gct,
   get_value(NewParameters,is_logged,unknown_logged,Is_Logged),
   Is_Logged=logged,
   get_value(NewParameters,status,created,Status),
   Status=created,
   member(OldStatus,[given,jointed,splited,created]),
   set_value(NewParameters,status,OldStatus,NewParameters1),
   signal_file(NewParameters1,Past_Modules_2),
   append(Past_Modules_1,Past_Modules_2,Past_Modules),
   Newadd=[gsea,NewParameters],
   append(Past_Modules,Newadd,Modules).


signature_score(Parameters,Modules):-
   convert_parameters_file(Parameters,NewParameters),
   get_value(NewParameters,format,unknown_format,Format),
   Format=pcl,
   get_value(NewParameters,is_logged,unknown_logged,Is_Logged),
   Is_Logged=logged,
   get_value(NewParameters,status,created,Status),
   Status=created,
   get_value(Parameters,platform,unknown_platform,Platform),
   Platform="HG_U133A",
   member(OldStatus,[given,jointed,splited,created]),
   set_value(NewParameters,status,OldStatus,NewParameters2),
   get_value(Parameters,preprocess,unknown_preprocess,Preprocess),
   (member(Preprocess,[rma,mas5]),
   set_value(NewParameters2,preprocess,rma,NewParameters3),
   set_value(NewParameters2,preprocess,mas5,NewParameters4),
   signal_file(NewParameters3,Past_Modules_1),
   signal_file(NewParameters4,Past_Modules_2),
   append(Past_Modules_2,Past_Modules_1,Past_Modules),
   append([pre1,rma,pre2,mas5],NewParameters,Write_list);
   Preprocess=illumina,
   signal_file(NewParameters2,Past_Modules),
   append([pre1,Preprocess,pre2,Preprocess],NewParameters,Write_list)),
   Newadd=[run_scoresig,Write_list],
   append(Past_Modules,Newadd,Modules).

