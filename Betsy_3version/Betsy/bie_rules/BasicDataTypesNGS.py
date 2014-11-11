#BasicDataTypesNGS
from Betsy.bie3 import *
import Database
FastqFile=DataType(
    'FastqFile',
    AttributeDef("read",['single','pair','pair1','pair2'],
                 "single","single",help='single or pair read'),
    AttributeDef("ref",['hg18','hg19','mm9','dm3'], "hg19","hg19",
                 help='ref species'),
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    help="Fastq file"
    )
SaiFile=DataType(
    'SaiFile',
    AttributeDef("read",['single','pair','pair1','pair2'],
                 "single","single",help='single or pair read'),
    AttributeDef("ref",['hg18','hg19','mm9','dm3'],
                 "hg19","hg19",help='ref species'),
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    help="Sai file"
    )

SamFile=DataType(
    'SamFile',
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',
                 help="contents"),
    AttributeDef("sorted",["yes","no"],"no","no",
                 help='sorted or not'),
    AttributeDef("duplicates_marked",["yes","no"],"no","no",
                 help='mark duplicate or not'),
    AttributeDef("recalibration",["yes","no"],"no","no",
                 help='recalibration or not'),
    AttributeDef("has_header",["yes","no"],"no","no",
                 help='fix header or not'),
    AttributeDef("read",['single','pair'],
                 "single","single",
                 help='single or pair read'),
    AttributeDef("ref",['hg18','hg19','mm9','dm3'],
                 "hg19","hg19",
                 help='ref species'),
    help="Sam file"
    )
BamFile=DataType(
    'BamFile',
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',
                 help="contents"),
    AttributeDef("sorted",["yes","no"],"no","no",
                 help='sorted or not'),
    AttributeDef("duplicates_marked",["yes","no"],"no","no",
                 help='mark duplicate or not'),
    AttributeDef("recalibration",["yes","no"],"no","no",
                 help='recalibration or not'),
    AttributeDef("has_header",["yes","no"],"no","no",
                 help='fix header or not'),
    AttributeDef("read",['single','pair'],
                 "single","single",
                 help='single or pair read'),
    AttributeDef("ref",['hg18','hg19','mm9','dm3'],
                 "hg19","hg19",
                 help='ref species'),
    help='Bam file'
    )
VcfFile=DataType(
    'VcfFile',
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',
                 help="contents"),
    AttributeDef("recalibration",["yes","no"],"no","no",
                 help='recalibration or not'),
    AttributeDef("read",['single','pair'],
                 "single","single",
                 help='single or pair read'),
    AttributeDef("ref",['hg18','hg19','mm9','dm3'],
                 "hg19","hg19",
                 help='ref species'),
    AttributeDef("vcf_filter",['yes','no'], "no","no",
                 help='filter VcfFile or not'),
    AttributeDef("reheader",['standard','bcftool'],
                 "standard","standard",
                 help='method to convert to VcfFile'),
    AttributeDef("vcf_annotate",['yes','no'], "no","no",
                 help='annotate VcfFile or not'),
    help='Vcf file'
    )
RNASeqFile=DataType(
    'RNASeqFile',
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    AttributeDef("format_type",['unknown','not_fastqfolder','not_samfolder','not_bamfolder',
                                'samfolder','bamfolder','fastqfolder'],
                 'unknown','unknown',help="format type"),
    help="RNA Seq File"
    )
FastqFolder=DataType(
    'FastqFolder',
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    help="RNA seq Fastq folder"
    )
SamFolder=DataType(
    'SamFolder',
    AttributeDef("ref",['human','mouse','hg18','hg19'], "human","human",
                 help='ref species'),
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    AttributeDef("sample_type",['RNA','DNA'],
                 'RNA','RNA',help="RNA or DNA type"),
    help="RNA seq Sam folder"
    )

BamFolder=DataType(
    'BamFolder',
    AttributeDef("ref",['human','mouse','hg18','hg19'], "human","human",
                 help='ref species'),
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    AttributeDef("sample_type",['RNA','DNA'],
                 'RNA','RNA',help="RNA or DNA type"),
    AttributeDef("duplicates_marked",["yes","no"],"no","no",
                 help='mark duplicate or not'),
    AttributeDef("sorted",["yes","no",'unknown'],"unknown","unknown",
                 help='sorted or not'),
    AttributeDef("indexed",["yes","no"],"no","no",
                 help='indexed or not'),
    help="RNA seq Bam folder"
    )
SampleGroupFile=DataType(
    'SampleGroupFile',
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    help="File contains sample group infomation"
    )
RNASeQCFile=DataType(
    'RNASeQCFile',
    AttributeDef("contents",Database.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    help="File contains sample group infomation"
    )
list_files = [FastqFile, SaiFile, SamFile, BamFile, VcfFile,
              RNASeqFile, FastqFolder, SamFolder, BamFolder,SampleGroupFile,RNASeQCFile]

all_modules = []
