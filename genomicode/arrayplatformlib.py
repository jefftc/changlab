#arrayplatformlib.py
"""

Functions:
find_platform_by_name    Return a Platform object with a given name.
get_bm_attribute
get_bm_organism
get_priority
prioritize_platforms

find_header

score_platforms                 Given annotations, score for each platform.
score_platform_of_annotations   Return the best match only.
score_all_platforms_of_matrix
score_platform_of_matrix

identify_platform_of_annotations
identify_all_platforms_of_matrix
identify_platform_of_matrix

chipname2filename
chipname2filename_illu
chipname2filename_affy


Classes:
Platform           A description for a probe on a microarray.

Variables:
PLATFORMS          List of all known platforms.

Constants:
PROBE_ID
GENE_ID
GENE_SYMBOL
DESCRIPTION


TODO:
- Need a function like:
  find_header(MATRIX, GENE_SYMBOL) -> "Gene Symbol"
- What is a chipname?
  What is a platform?  How are the names determined?
  score_all_platforms_of_matrix returns entrez_ID_symbol_human.
  Looks like a chipname.  Doesn't match the name of any Platform.
- What are the headers at the bottom?

"""
import os
import re
from genomicode import filelib, config


PROBE_ID, GENE_ID, GENE_SYMBOL, DESCRIPTION = range(4)


class Platform:
    def __init__(self, name, bm_attribute, bm_organism, category, priority):
        # our description of platform, corresponds to files
        # in genemodata/pid2platform, also matches annotation files in
        # /data/genomidata/affymetrix and /data/genomidata/illumina
        self.name = name
        # bm_organism   Which BioMart organism to use.  XXX how to get list?
        # bm_attribute  Which attribute from this mart.
        #
        # To list valid attributes:
        # library(biomaRt)
        # x <- useMart("ensembl", "hsapiens_gene_ensembl")
        # y <- listAttributes(x)
        # y$name
        self.bm_attribute = bm_attribute
        self.bm_organism = bm_organism
        self.category = category  # PROBE_ID, GENE_ID, GENE_SYMBOL, DESCRIPTION
        self.priority = priority # order of platform priority, lower number means higher priority


def find_platform_by_name(name):
    for platform in PLATFORMS:
        if platform.name == name:
            return platform
    return None

    
def get_bm_attribute(platform_name):
    platform = find_platform_by_name(platform_name)
    if platform is None:
        return None
    return platform.bm_attribute


def get_bm_organism(platform_name):
    platform = find_platform_by_name(platform_name)
    if platform is None:
        return None
    return platform.bm_organism


def get_priority(platform_name):
    platform = find_platform_by_name(platform_name)
    if platform is None:
        return None
    return platform.priority


def prioritize_platforms(platform_names):
    """highest priority to lowest, lower number is higher priority"""
    order_prioritize = [
        (get_priority(name), name) for name in platform_names]
    order_prioritize.sort()
    prioritized_names = [i[1] for i in order_prioritize]
    return prioritized_names


def find_header(MATRIX, category):
    # Returns a header from MATRIX that matches this category or None.

    # XXX NEED TO FIX THIS
    name_fix = {
        "entrez_ID_human" : "Entrez_ID_human",
        "entrez_ID_symbol_human" : "Entrez_Symbol_human",
        }
    
    for x in score_all_platforms_of_matrix(MATRIX):
        header, platform_name, score = x
        platform_name = name_fix.get(platform_name, platform_name)
        
        if score < 0.5:
            continue
        platform = find_platform_by_name(platform_name)
        if not platform:
            continue
        if platform.category != category:
            continue
        return header
    return None

def _hash_chipname(filename):
    # RG_U34A_annot.csv.gz           RG_U34A
    # Rat230_2.na26.annot.csv.gz     Rat230_2
    # MoEx_1_0_st_v1.1.mm9.probeset  MoEx_1_0_st_v1
    REMOVE = [
        ".gz", ".csv", ".annot", "_annot",
        ".mm9", ".probeset", ".1",
        ]
    x = os.path.split(filename)[1]
    for r in REMOVE:
        x = x.replace(r, "")
    x = x.replace("-", "_")
    version = 0
    m = re.search(r".na(\d+)",x)
    if m:
        version = int(m.group(1))
        x = x.replace(m.group(0), '')
    return x, version


def chipname2filename(chipname):
    filename = chipname2filename_affy(chipname)
    if not filename:
        filename = chipname2filename_illu(chipname)
    return filename


def chipname2filename_illu(chipname):
    filename = None
    path = config.annot_data_illu
    assert os.path.exists(path), '%s does not exist' % path
    chipname = chipname.replace('_', '-')
    for f in os.listdir(path):
        if chipname in f:
            filename = os.path.join(path, f)
    return filename


def chipname2filename_affy(chipname):
    path = config.annot_data_affy
    assert os.path.exists(path), '%s does not exist' % path

    # Read the files in config.annot_data_affy path.
    chip2file = {}  # chip -> (version, filename)
    for f in os.listdir(path):
        chip, version = _hash_chipname(f)

        x = chip2file.get(chip, (None, None))
        old_version, old_filename = x
        if old_version is None or version > old_version:
            filename = os.path.join(path, f)
            chip2file[chip] = (version, filename)
            
    if chipname not in chip2file:
        return None
    version, filename = chip2file[chipname]
    return filename


def _read_annotations_h():
    #/data/genomidata/pid2platform has 2 folders, case_sensitive and
    # case_insensitive. Platforms in case_sensitive need to match the ids
    # in exactly uppercase and lowercase, while platforms in case_insensitive
    # do not care the uppercase and lowercase.
    
    paths = []
    result = []
    
    root = config.psid2platform
    assert os.path.exists(root), "path not exist: %s" % root
    for subfolder in os.listdir(root):
        if '.DS_Store' in subfolder:
            continue
        assert os.path.isdir(os.path.join(root, subfolder))
        for platform in os.listdir(os.path.join(root, subfolder)):
            paths.append((root, subfolder, platform))
    all_platform = []
    for platform in PLATFORMS:
        all_platform.append(platform.name)
    for x in paths:
        root, subfolder, platform = x
        platform_name = platform[:-4] #get rid of .txt
        assert platform_name in all_platform,('%s is not a platform in PLATFORMS'
                                         %platform_name)
        assert subfolder in ['case_sensitive','case_insensitive']
        f = file(os.path.join(root,subfolder,platform),'r')
        text = f.readlines()
        text = [i.strip() for i in text if len(i.strip())>0]
        f.close()
        chipname = os.path.splitext(platform)[-2]  # remove the '.txt'
        if subfolder == 'case_insensitive':
            result.append((chipname, False, text))
        else:
            result.append((chipname, True, text))
    return result


ALL_PLATFORMS = None
def _read_annotations():
    """Return list of (chip_name, case_sensitive, list of IDs.)"""
    global ALL_PLATFORMS
    if not ALL_PLATFORMS:
        ALL_PLATFORMS = _read_annotations_h()
    return ALL_PLATFORMS


def score_platforms(annotations):
    all_platforms = _read_annotations()
    results = []
    for x in all_platforms:
        chipname, case_sensitive, ids = x
        y = _compare_annotations(annotations, ids, case_sensitive)
        number_shared_annots, only1, only2, match = y
        results.append((chipname, number_shared_annots, match))
    results.sort(key=lambda x: (x[2], x[1]),reverse=True)
    return results


def _compare_annotations(annots1, annots2, case_sensitive):
    # Return tuple of: num_shared, num_annots1_only, num_annots2_only,
    # perc_shared.
    if not case_sensitive:
        annots1 = [psid.upper() for psid in annots1]
        annots2 = [psid.upper() for psid in annots2]
    # Ignore differences in whitespace.
    annots1 = [x.strip() for x in annots1]
    annots2 = [x.strip() for x in annots2]
    # No blank annotations.
    annots1 = [x for x in annots1 if x]
    annots2 = [x for x in annots2 if x]
    # No duplicate annotations.
    annots1 = {}.fromkeys(annots1)
    annots2 = {}.fromkeys(annots2)

    x = [x for x in annots1 if x in annots2]
    num_shared_annots = len(x)
    num_annots1_only = len(annots1) - num_shared_annots
    num_annots2_only = len(annots2) - num_shared_annots
    perc_shared = 0
    if annots1 and annots2:
        perc_shared = float(num_shared_annots)/min(len(annots1), len(annots2))
    x = num_shared_annots, num_annots1_only, num_annots2_only, perc_shared
    return x


def score_platform_of_annotations(annotations):
    # Tuple of (platform, match score) or None if nothing matches.
    possible_platforms = score_platforms(annotations)
    assert possible_platforms
    if possible_platforms[0][1] <= 0:
        return None
    platform, number_shared_annots, match = possible_platforms[0]
    return platform, match


def _parse_matrix_annotations(annots, delim):
    if delim is None:
        return annots
    parsed = []
    for annot in annots:
        x = annot.split(delim)
        parsed.extend(x)
    return parsed


def score_all_platforms_of_matrix(DATA, annot_delim=None):
    """Return a list of (header, platform, score).  score is the
    percentage of the IDs from header that matches the IDs in this
    platform.  List is ordered by increasing platform priority.

    """
    platform2score = {}  # platform_name -> header, match
    for header in DATA.row_names():
        x = DATA.row_names(header)
        annots = _parse_matrix_annotations(x, annot_delim)
        x = score_platform_of_annotations(annots)
        if x is None:
            continue
        platform, score = x
        platform2score[platform] = (header, score)
    order_platforms = prioritize_platforms(platform2score.keys())
    # if chips is empty, will return an empty list
    x = [(platform2score[platform][0], platform, platform2score[platform][1])
         for platform in order_platforms]
    return x


def score_platform_of_matrix(DATA, annot_delim=None):
    # Return the best scoring platform for this matrix.
    platforms = score_all_platforms_of_matrix(DATA, annot_delim=annot_delim)
    if not platforms:
        return None, 0

    # Sort by decreasing score.
    schwartz = [(-x[2], x) for x in platforms]
    schwartz.sort()
    platforms = [x[-1] for x in schwartz]
    
    out_platform = platforms[0][1]
    out_platform_match = platforms[0][2]
    return out_platform, out_platform_match


def identify_platform_of_annotations(annotations):
    x = score_platform_of_annotations(annotations)
    if not x:
        return None
    platform, match = x
    if match == 1:
        return platform
    return None


def identify_all_platforms_of_matrix(DATA, annot_delim=None):
    """return a list of (header, platform) we can identify"""
    platforms = score_all_platforms_of_matrix(DATA, annot_delim=annot_delim)
    result = []
    for x in platforms:
        header, platform, match = x
        if match == 1:
            result.append((header, platform))
    return result

    
def identify_platform_of_matrix(DATA, annot_delim=None):
    x = score_platform_of_matrix(DATA, annot_delim=annot_delim)
    platform_name, match = x
    if match == 1:
        return platform_name
    return None


PLATFORMS = [
    Platform('HG_U95A', "affy_hg_u95a", "hsapiens_gene_ensembl", PROBE_ID, 1),
    Platform(
        'HG_U95Av2', "affy_hg_u95av2", "hsapiens_gene_ensembl", PROBE_ID, 2),
    Platform(
        'HG_U133_Plus_2', "affy_hg_u133_plus_2", "hsapiens_gene_ensembl",
        PROBE_ID, 3),
    Platform(
        'HG_U133A_2', "affy_hg_u133a_2", "hsapiens_gene_ensembl", PROBE_ID, 4),
    Platform(
        'HG_U133A', "affy_hg_u133a", "hsapiens_gene_ensembl", PROBE_ID, 5),
    Platform(
        'HG_U133B', "affy_hg_u133b", "hsapiens_gene_ensembl", PROBE_ID, 6),
    Platform('Hu35KsubA', None, None, PROBE_ID, 7),
    Platform('Hu35KsubB', None, None, PROBE_ID, 8),
    Platform('Hu35KsubC', None, None, PROBE_ID, 9),
    Platform('Hu35KsubD', None, None, PROBE_ID, 10),
    Platform('Hu6800', "affy_hugenefl", "hsapiens_gene_ensembl", PROBE_ID, 11),
    Platform(
        'MG_U74Av2', "affy_mg_u74av2", "mmusculus_gene_ensembl", PROBE_ID, 12),
    Platform(
        'MG_U74Bv2','affy_mg_u74bv2', "mmusculus_gene_ensembl", PROBE_ID, 13),
    Platform(
        'MG_U74Cv2','affy_mg_u74cv2', "mmusculus_gene_ensembl", PROBE_ID, 14),
    Platform(
        'Mouse430_2','affy_mouse430_2', "mmusculus_gene_ensembl", PROBE_ID,
        15),
    Platform(
        'Mouse430A_2','affy_mouse430a_2', "mmusculus_gene_ensembl", PROBE_ID,
        16),
    Platform(
        'Mu11KsubA','affy_mu11ksuba', "mmusculus_gene_ensembl", PROBE_ID, 17),
    Platform(
        'Mu11KsubB', "affy_mu11ksubb", "mmusculus_gene_ensembl", PROBE_ID, 18),
    Platform(
        'RG_U34A','affy_rg_u34a', "rnorvegicus_gene_ensembl", PROBE_ID, 19),
    Platform(
        'RAE230A','affy_rae230a', "rnorvegicus_gene_ensembl", PROBE_ID, 20),
    Platform(
        'HumanHT_12', "illumina_humanht_12_v4", "hsapiens_gene_ensembl",
        PROBE_ID, 21),
    Platform(
        'HumanWG_6', "illumina_humanwg_6_v3", "hsapiens_gene_ensembl",
        PROBE_ID, 22),
    Platform(
        'MouseRef_8', "illumina_mousewg_6_v2", "mmusculus_gene_ensembl",
        PROBE_ID, 23),
    Platform(
        'Entrez_ID_human', "entrezgene", "hsapiens_gene_ensembl", GENE_ID, 24),
    Platform(
        'Entrez_ID_mouse', "entrezgene", "mmusculus_gene_ensembl", GENE_ID,
        25),
    Platform(
        'Entrez_Symbol_human', "hgnc_symbol", "hsapiens_gene_ensembl",
        GENE_SYMBOL, 26),
    Platform(
        'Entrez_Symbol_mouse', "mgi_symbol", "mmusculus_gene_ensembl",
        GENE_SYMBOL, 27),
    Platform(
        'Ensembl_human', "ensembl_gene_id", "hsapiens_gene_ensembl",
        PROBE_ID, 28),
    Platform(
        'Ensembl_mouse', "ensembl_gene_id", "mmusculus_gene_ensembl",
        PROBE_ID, 29),
    Platform(
        'UCSC_human_hg19_kg7', None, None, GENE_ID, 30), 
    Platform(
        'UCSC_human_hg19_kg6', None, None, GENE_ID, 31),  
    Platform(
        'UCSC_human_hg19_kg5', None, None, GENE_ID, 32),
    Platform(
        'UCSC_human_hg38_kg8', None, None, GENE_ID, 33), 
    Platform(
        'Agilent_Human1A', None, None, GENE_ID, 34), 
    Platform(
        'HumanHT_12_control', None, None, GENE_ID, 35),  
    Platform(
        'MouseRef_8_control', None, None, GENE_ID, 36),
    Platform(
        'UCSC_human_hg38_kg8', None, None, GENE_ID, 37),
    Platform(
        'RefSeq_protein_mouse', "refseq_peptide", "mmusculus_gene_ensembl",   
        GENE_SYMBOL, 38),
    Platform(
        'RefSeq_predicted_protein_mouse', "refseq_peptide_predicted", "mmusculus_gene_ensembl",   
        GENE_SYMBOL, 39),

    Platform(
        'RefSeq_protein_human', "refseq_peptide", "hsapiens_gene_ensembl",
        PROBE_ID, 40),
    Platform(
        'RefSeq_predicted_protein_human', "refseq_peptide_predicted", "hsapiens_gene_ensembl",
        PROBE_ID, 41),
    Platform(
        'RefSeq_mRNA_mouse', "refseq_mrna", "mmusculus_gene_ensembl",   
        GENE_SYMBOL, 42),
    Platform(
        'RefSeq_predicted_mRNA_mouse', "refseq_mrna_predicted", "mmusculus_gene_ensembl",   
        GENE_SYMBOL, 43),
    Platform(
        'RefSeq_ncRNA_mouse', "refseq_ncrna", "mmusculus_gene_ensembl",   
        GENE_SYMBOL, 44),
    Platform(
        'RefSeq_predicted_ncRNA_mouse', "refseq_ncrna_predicted", "mmusculus_gene_ensembl",   
        GENE_SYMBOL, 45),
    Platform(
        'RefSeq_mRNA_human', "refseq_mrna", "hsapiens_gene_ensembl",
        PROBE_ID, 46),
    Platform(
        'RefSeq_predicted_mRNA_human', "refseq_mrna_predicted", "hsapiens_gene_ensembl",
        PROBE_ID, 47),
    Platform(
        'RefSeq_ncRNA_human', "refseq_ncrna", "hsapiens_gene_ensembl",
        PROBE_ID, 48),
    Platform(
        'RefSeq_predicted_ncRNA_human', "refseq_ncrna_predicted", "hsapiens_gene_ensembl",
        PROBE_ID, 49),
    Platform(
        'UCSC_mouse_mm10_kg7', None, None, GENE_ID, 50),
    ]



annotate_header = [#The added header names when annotate a file
    'Description', 'Gene Symbol', 'Gene ID', 'Swiss-Prot ID']

affy_headers = {# Header names in affymetrix annotation
                # files in /data/genomidata/affymetrix
    'Description' : ['Target Description'],
    'Gene Symbol' : ['Gene Symbol'],
    'Gene ID' : ['Entrez Gene'],
    'Swiss-Prot ID' : ['SwissProt']
    }
illu_headers = {# Header names in illumina annotation
                # files in /data/genomidata/illumina
    'Description' : ['Definition'],
    'Gene Symbol' : ['Symbol'],
    'Gene ID' : ['Entrez_Gene_ID'],
    'Swiss-Prot ID' : ['swissport_id']
    }
biomart_headers = {# Header names in biomart attributesL
    'Description' : ['description'],
    'Gene Symbol' : ['hgnc_symbol','mgi_symbol'],
    'Gene ID' : ['entrezgene'],
    'Swiss-Prot ID' : ['uniprot_swissprot']
    }

platform_to_GSEA_chipname={# corresponds to chipname in genepattern GSEA module
    'Agilent_Human1A':'Agilent_Human1A',
    'HG_U133A_2':'HG_U133A_2',
    'HG_U133A':'HG_U133A',
    'HG_U133B':'HG_U133B',
    'HG_U133_Plus_2':'HG_U133_Plus_2',
    'HG_U95Av2':'HG_U95Av2',
    'Hu35KsubA':'Hu35KsubA',
    'Hu35KsubB':'Hu35KsubB',
    'Hu35KsubC':'Hu35KsubC',
    'Hu35KsubD':'Hu35KsubD',
    'Hu6800':'Hu6800',
    'MG_U74Av2':'MG_U74Av2',
    'MG_U74Bv2':'MG_U74Bv2',
    'MG_U74Cv2':'MG_U74Cv2',
    'Mouse430_2':'Mouse430_2',
    'Mouse430A_2':'Mouse430A_2',
    'Mu11KsubA':'Mu11KsubA',
    'Mu11KsubB':'Mu11KsubB',
    'RAE230A':'RAE230A',
    'RG_U34A':'RG_U34A',
    'HumanWG_6':'ilmn_HumanWG_6_V3_0_R3_11282955_A',
    'HumanHT_12':'ilmn_HumanHT_12_V3_0_R3_11283641_A',
    'MouseRef_8':'ilmn_MouseRef_8_V2_0_R3_11278551_A'
    } 
