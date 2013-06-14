%vcf_file
%Output interface
:- set_prolog_flag(toplevel_print_options,
             [quoted(true), portray(true), max_depth(0)]).

:-dynamic bam_file/2.
:-dynamic sam_file/2.
:-dynamic fastq_file/2.
:-dynamic input_sam_file/2.

vcf_file(Parameters,Modules):-
    % Conditions: the Parameters is not full length
    get_desire_parameters_vcf(Parameters,NewParameters1),
    length(NewParameters1,N),
    get_length(n_vcf,N1),
    N<N1,
    write(yes),
    % Input:vcf_file with full length parameters
    convert_parameters_variable_vcf(Parameters,NewParameters),
    write(NewParameters),nl,
    vcf_file(NewParameters,Modules).
/*-------------------------------------------------------------------------*/
%vcf_file through GATK standard experimental
/*vcf_file(Parameters,Modules):-

      get_value(Parameters,filter,no_filter,Filter),
      Filter=no_filter,

	get_value(Parameters,sorted,no_sorted,Sorted),
      Sorted=yes_sorted,
      
      get_value(Parameters,duplicates_marked,no_marked,Duplicates_marked),
      Duplicates_marked=yes_marked,

      get_value(Parameters,ref,human,Ref),
      get_value(Parameters,recalibration,no_recalibration,Recalibration),
      (Ref=human,
      Recalibration=yes_recalibration;
      Ref=fly,
      Recalibration=no_recalibration),

	get_value(Parameters,has_header,no_fixed,Has_header),
      Has_header=yes_fixed,

      get_value(Parameters,reheader,standard,Reheader),
      (Reheader=standard,
      
      % Input:sam_file
      
      sam_file(Parameters,Past_Modules),

	% Module: call_variants_GATK
      % Output parameters:full length parameters of vcf_file
      Newadd=[call_variants_GATK,Parameters];
      Reheader=bcftool,
      % Input:mpileup_file
      mpileup_file(Parameters,Past_Modules),

	% Module: reheader_bcftool
      % Output parameters:full length parameters of vcf_file
      Newadd=[reheader_bcftool,Parameters]),
      append(Past_Modules,Newadd,Modules).
*/
%-----------------------------------------------------------------------
vcf_file(Parameters,Modules):-
      length(Parameters,N),
      get_length(n_vcf,N1),
      N=N1,
      
      get_value(Parameters,filter,no_filter,Filter),

      get_value(Parameters,sorted,no_sorted,Sorted),
      Sorted=yes_sorted,
      
      get_value(Parameters,duplicates_marked,no_marked,Duplicates_marked),
      Duplicates_marked=yes_marked,
     
      get_value(Parameters,ref,hg19,Ref),
      get_value(Parameters,recalibration,no_recalibration,Recalibration),
      (member(Ref,[hg18,hg19,mm9,dm3]),
      member(Recalibration,[yes_recalibration,no_recalibration]);
      %Recalibration=yes_recalibration;
      member(Ref,[dm3,mm9]),
      Recalibration=no_recalibration),
      
	get_value(Parameters,has_header,no_fixed,Has_header),
      Has_header=yes_fixed,

      get_value(Parameters,reheader,standard,Reheader),
      
      % Input:sam_file
      convert_parameters_sam(Parameters,NewParameters),
      sam_file(NewParameters,Past_Modules),
    
      (Reheader=standard,
      Filter=yes_filter,
	% Module: call_variants_GATK
      % Output parameters:full length parameters of vcf_file
      Newadd=[call_variants_GATK,Parameters];
      Reheader=bcftool,
      Filter=no_filter,
	% Module: call_variants_mpileup
      Newadd=[call_variants_mpileup,Parameters]),
      append(Past_Modules,Newadd,Modules).
/*-------------------------------------------------------------------------*/
% filter the vcf_file
vcf_file(Parameters,Modules):-
      length(Parameters,N),
      get_length(n_vcf,N1),
      N=N1,
      get_value(Parameters,filter,no_filter,Filter),
      Filter=yes_filter,

      % Input:vcf_file
      OldFilter = no_filter,
      set_value(Parameters,filter,OldFilter,OldParameters),

      vcf_file(OldParameters,Past_Modules),

      % Module: filter_vcf_file
      %Output parameters:full length parameters of vcf_file
      Newadd=[filter_vcf_file,Parameters],
      append(Past_Modules,Newadd,Modules).

/*-------------------------------------------------------------------------*/
% mpileup_file

/*mpileup_file(Parameters,Modules):-
       
      convert_parameters_variable_sam(Parameters,NewParameters),
    	get_value(NewParameters,sorted,no_sorted,Sorted),
      Sorted=yes_sorted,
      
      get_value(NewParameters,duplicates_marked,no_marked,Duplicates_marked),
      Duplicates_marked=yes_marked,
      
      get_value(NewParameters,recalibration,no_recalibration,Recalibration),
      Recalibration=yes_recalibration,
      
	get_value(NewParameters,has_header,no_fixed,Has_header),
      Has_header=yes_fixed,

      % Input:sam_file
      sam_file(NewParameters,Past_Modules),
     
	% Module: call_variants_mpileup
      % Output parameters:full length parameters of mpileup_file
      Newadd=[call_variants_mpileup,NewParameters],
      append(Past_Modules,Newadd,Modules).*/
/*-------------------------------------------------------------------------*/
% qc_file

qc_file(Parameters,Modules):-
     get_value(Parameters,format,unknown_format,Format),
     member(Format,[bam,sam]),
     convert_parameters_variable_sam(Parameters,NewParameters),
     sam_file(NewParameters,Past_Modules),
     %Module:check_quality
     %Output parameters:full length parameters of vcf_file
     Newadd=[check_quality,NewParameters],
     append(Past_Modules,Newadd,Modules).
/*-------------------------------------------------------------------------*/
convert_parameters_variable_vcf(Parameters,NewParameters):-
    convert_parameters_variable_sam(Parameters,NewParameters1),

    get_value_variable(Parameters,filter,Filter),
    member(Filter,[yes_filter,no_filter]),
    append(NewParameters1,[filter,Filter],NewParameters2),
    
    get_value_variable(Parameters,reheader,Reheader),
    member(Reheader,[standard,bcftool]),
    append(NewParameters2,[reheader,Reheader],NewParameters).
/*-------------------------------------------------------------------------*/
get_desire_parameters_vcf(Parameters,NewParameters):-
    get_desire_parameters_sam(Parameters,NewParameters1),
    (member(filter,Parameters),
    get_value_variable(Parameters,filter,Filter),
    append(NewParameters1,[filter,Filter],NewParameters2);
    not(member(filter,Parameters)),
    NewParameters2=NewParameters1),

    (member(reheader,Parameters),
    get_value_variable(Parameters,reheader,Reheader),
    append(NewParameters2,[reheader,Reheader],NewParameters);
    not(member(reheader,Parameters)),
    NewParameters=NewParameters2).
