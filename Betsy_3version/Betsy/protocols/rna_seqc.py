#rna_seqc.py
from Betsy.protocol_utils import Parameter
import normalize_file
PRETTY_NAME="RNA Sequence Quality Control analysis."

OPTIONAL = 'Optional Parameters'
RNA = 'RNA Sequence analysis parameters'

CATEGORIES=[OPTIONAL,RNA]



#input  datatype
INPUTS = [
    'BamFolder',
    'SampleGroupFile',
    ]

#output datatype
OUTPUTS = 'RNASeQCFile'

PARAMETERS=[Parameter('contents',pretty_name='Contents',category=OPTIONAL,
                      choices=normalize_file.CONTENTS, description='output group information',
                      datatype='RNASeQCFile'),
            Parameter('ref',pretty_name='reference species',default_value='hg19',
                        choices=['hg18', 'hg19'],category=RNA,
                      description='reference species for BamFolder',datatype='BamFolder'),
            Parameter('duplicates_marked',pretty_name='mark duplicates',default_value='yes',
                        choices=['yes'],category=RNA,
                      description='mark duplicates in BamFolder',datatype='RNASeQCFile'),
            Parameter('indexed',pretty_name='indexed file',default_value='yes',
                        choices=['yes'],category=RNA,
                      description='indexed file in BamFolder',datatype='RNASeQCFile'),
            Parameter('sorted',pretty_name='sort bam file',default_value='yes',
                        choices=['yes'],category=RNA,
                      description='sort file in BamFolder',datatype='RNASeQCFile'),
            Parameter('sample_type',pretty_name='RNA or DNA ',default_value='RNA',
                        choices=['RNA'],category=RNA,
                      description='must be RNA',datatype='RNASeQCFile'),
            Parameter('RNA_ref',pretty_name='ref file for RNA_SeQC',
                           category=OPTIONAL,description='ref file for RNA_SeQC',
                      datatype='UserUpload'),
            Parameter('RNA_gtf',pretty_name='gtf file for RNA_SeQC',
                           category=OPTIONAL,description='gtf file for RNA_SeQC',
                      datatype='UserUpload'),
            ]

