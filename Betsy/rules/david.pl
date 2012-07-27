/*David enrichment*/

/*--------------------------------------------------------------------------*/
%run david website
david(Parameters,Modules):-
   get_value(Parameters,contents,[unknown],Contents),
   % Input: from the given gene_list_file
   (gene_list_file([contents,Contents],[]),
   % Output parameters:update the output parameters to full length
   Newadd=[contents,Contents],
   % Module: run_david
   Modules=[run_david,Newadd];
   % Conditions: not given a gene_list_file
   not(gene_list_file([contents,Contents],[])),
   convert_parameters_file(Parameters,NewParameters1),
   % Input: gene_list_file not given
   gene_list_file(NewParameters1,Past_Modules),
   %append(NewParameters1,NewParameters,NewParameters2),
   % Module:run_david
   % Output parameters:update the output parameters to full length
   Newadd=[run_david,NewParameters1],
   append(Past_Modules,Newadd,Modules)).