#dna_seq.py
from Betsy.protocol_utils import Parameter
import normalize_file
PRETTY_NAME="Variant Calling."
OPTIONAL = 'Optional Parameters'
DNA = 'DNA Sequence analysis parameters'


CATEGORIES=[OPTIONAL,DNA]

#input  datatype
INPUTS = [
    'FastqFile',
    'SaiFile',
    'SamFile',
    'BamFile',
    'VcfFile',
    ]

#output datatype
OUTPUTS = 'VcfFile'

#parameter objects
PARAMETERS=[Parameter('contents',pretty_name='Contents',category=OPTIONAL,
                      choices=normalize_file.CONTENTS, description='output group information',
                      datatype='VcfFile'),
            Parameter('recalibration',pretty_name='recalibration',
                      choices=["yes", "no"],category=DNA,
                      description='recalibration or not',datatype='VcfFile'),
            Parameter('read',pretty_name='single or pair read',default_value='single',
                        choices=['single', 'pair'],category=DNA,
                      description='single or pair read',datatype='VcfFile'),
            Parameter('ref',pretty_name='Reference',
                         choices=['hg18','hg19','mm9','dm3'],default_value='hg19',
                      category=DNA, description='reference species',datatype='VcfFile'),
            Parameter('vcf_filter',pretty_name='filter vcf file',
                         choices=['yes','no'],default_value='yes',
                      category=DNA, description='filter VcfFile or not',datatype='VcfFile'),
            Parameter('reheader',pretty_name='reheader method',
                         choices=['standard','bcftool'],default_value='standard',
                      category=DNA, description='method to convert to VcfFile',datatype='VcfFile'),
            Parameter('vcf_annotate',pretty_name='annotate vcf file',
                         choices=['yes','no'],default_value='no',
                      category=DNA, description='annotate VcfFile or not',datatype='VcfFile'),

            ]

