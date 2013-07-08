#call_variants.py

from Betsy.protocol_utils import Parameter
PRETTY_NAME="Calling Variants for Next Generation Sequence Data."
COMMON = 'Common Parameters'
OPTIONAL = 'Optional Parameters'
CATEGORIES=[COMMON,OPTIONAL]

#input predicates
INPUTS=['input_sam_file',
        'fastq_file']

#output predicates
OUTPUTS = 'make_call_variants_report'

# <input predicate name> : information for converting to fact
predicate2arguments = {
    'input_sam_file': (['read','single','ref','hg19'], '[]'),
    'fastq_file': (['read','single','ref','hg19'], '[]')}

PARAMETERS=[Parameter('contents',pretty_name='Contents',
                      type='list',category=OPTIONAL,description='output group information'),
            Parameter('read',pretty_name='Single or Paired',
                      default_value='single',choices=['single','pair'],
                      category=COMMON,description='single or paired'),
            Parameter('ref',pretty_name='Reference genome',default_value='hg19',
                      choices=['hg18','hg19','dm3','mm9'],category=COMMON,
                      description='select a specie of your data'),
            Parameter('recalibration',pretty_name='Base Quality Score Recalibration',
                      default_value='no_recalibration',choices=['no_recalibration','yes_recalibration'],
                      category=COMMON,description='do the recalibration for human species'),
            Parameter('reheader',pretty_name='Reheader',default_value='standard',
                      choices=['standard','bcftool'],category=COMMON,
                      description='method of calling variants')]

   #Parameter('format',pretty_name='File Format',default_value='unknown_format',
            #            choices=['unknown_format','sai', 'sam','bam'],category=COMMON),
            #Parameter('sorted',pretty_name='Sorted or Not',default_value='no_sorted',
            #          choices=['no_sorted','yes_sorted'],category=COMMON),
            #Parameter('duplicates_marked',pretty_name='Mark Duplicates',default_value='no_marked',
            #          choices=['no_marked','yes_marked'],category=COMMON),
            
            #Parameter('has_header',pretty_name='Fix Header',default_value='no_fixed',
            #          choices=['no_fixed','yes_fixed'],category=COMMON),
            #Parameter('filter',pretty_name='Filter the Vcf File',default_value='yes_filter',
            #          choices=['no_filter','yes_filter'],category=COMMON)
