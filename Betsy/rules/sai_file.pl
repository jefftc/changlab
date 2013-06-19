% sai_file

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

:-dynamic bam_file/2.
:-dynamic sam_file/2.
:-dynamic fastq_file/2.
:-dynamic input_sam_file/2.
/*-------------------------------------------------------------------------*/
% reorder the parameters from the sai_files with parameters in unwanted order
sai_files([contents,Contents,format,Format,read,Read,ref,Ref],Modules):-
   sai_files([format,Format,contents,Contents,read,Read,ref,Ref],Modules).

fastq_file([contents,Contents,read,Read,ref,Ref],[]):-
   fastq_file([read,Read,ref,Ref,contents,Contents],[]).
/*-------------------------------------------------------------------------*/
%generate sai_file from fastq_file.
sai_file([contents,Contents,format,sai,read,Read,ref,Ref],Modules):-
	% Input:fastq_file
	fastq_file([contents,Contents,read,Read,ref,Ref],[]),
      % Module: align_sequence  
      %Output parameters:[contents,Contents,format,sai,read,Read,ref,Ref]
      Newadd=[align_sequence,[contents,Contents,format,sai,read,Read,ref,Ref]],
      append([],Newadd,Modules).


convert_parameters_sai(Parameters,NewParameters):-
    get_value(Parameters,contents,[unknown],Contents),
    NewParameters0=[contents,Contents],
 
    get_value(Parameters,format,unknown_format,Format),
    member(Format,[sai]),
    append(NewParameters0,[format,Format],NewParameters1),
     
    get_value(Parameters,read,single,Read),
    member(Read,[single,pair]),
    append(NewParameters1,[read,Read],NewParameters2),
    
    get_value(Parameters,ref,human,Ref),
    member(Ref,[hg18,hg19,mm9,dm3]),
    append(NewParameters2,[ref,Ref],NewParameters).