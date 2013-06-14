% rna_seq.pl

:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

:-dynamic rna_fastq_file/2.
:-dynamic rna_sam_file/2.
:-dynamic rna_signal_raw/2.
/*-------------------------------------------------------------------------*/
% reorder the parameters from the sai_files with parameters in unwanted order

rna_fastq_file([contents,Contents,format,fa,read,Read,ref,Ref],[]):-
   rna_fastq_file([read,Read,ref,Ref,contents,Contents,format,fa],[]).

rna_sam_file([contents,Contents,format,Format,read,Read,ref,Ref],[]):-
   member(Format,[sam,bam]),
   rna_sam_file([read,Read,ref,Ref,contents,Contents,format,Format],[]).

/*-------------------------------------------------------------------------*/
%generate signal_file from rna_fastq_file 
rna_signal_raw(Parameters,Modules):-
      get_value(Parameters,contents,[unknown],Contents),
	get_value(Parameters,read,single,Read),
	get_value(Parameters,ref,hg19,Ref),
      get_value(Parameters,format,unknown_format,Format),
      Format=fa,
      (Read=single,
	% Input:rna_fastq_file
	rna_fastq_file([contents,Contents,format,fa,read,single,ref,Ref],Past_Modules);
      Read=pair,
      rna_fastq_file([contents,Contents,format,fa,read,pair1,ref,Ref],Past_Modules1),
      rna_fastq_file([contents,Contents,format,fa,read,pair2,ref,Ref],Past_Modules2),
      append(Past_Modules1,Past_Modules2,Past_Modules)),
      % Module: convert_fastq_to_signal_file
      %Output parameters:full length sam_file parameters
      Newadd=[convert_fastq_to_signal_file,Parameters],
      append(Past_Modules,Newadd,Modules).
      
/*-------------------------------------------------------------------------*/
%generate signal_file from rna_sam_file 
rna_signal_raw(Parameters,Modules):-
      get_value(Parameters,contents,[unknown],Contents),
	get_value(Parameters,read,single,Read),
	get_value(Parameters,ref,hg19,Ref),
	get_value(Parameters,format,unknown_format,Format),
      member(Format,[bam,sam]),
      (Read=single,
	% Input:rna_sam_file
	rna_sam_file([contents,Contents,format,Format,read,single,ref,Ref],Past_Modules);
      Read=pair,
      rna_sam_file([contents,Contents,format,Format,read,pair1,ref,Ref],Past_Modules1),
      rna_sam_file([contents,Contents,format,Format,read,pair2,ref,Ref],Past_Modules2),
      append(Past_Modules1,Past_Modules2,Past_Modules)),
      % Module: convert_sam_to_signal_file
      %Output parameters:full length sam_file parameters
      Newadd=[convert_sam_to_signal_file,Parameters],
      append(Past_Modules,Newadd,Modules).
/*-------------------------------------------------------------------------*/
signal_raw(Parameters,Modules):-
      get_value(Parameters,contents,[unknown],Contents),
      member(Format,[sam,bam,fa]),
      member(Read,[single,pair]),
      member(Ref,[hg18,hg19,mm9,dm3]),
      NewParameters=[contents,Contents,format,Format,read,Read,ref,Ref],
      rna_signal_raw(NewParameters,Modules).
      

