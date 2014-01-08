#arrayplatformlib.py
"""

Functions:
find_platform_by_name
get_bm_attribute
get_bm_organism
get_priority
prioritize_platforms

find_header

score_platforms
score_platform_of_annotations
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
        self.name = name
        self.bm_attribute = bm_attribute
        self.bm_organism = bm_organism
        self.category = category  # PROBE_ID, GENE_ID, GENE_SYMBOL, DESCRIPTION
        self.priority = priority


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
    order_prioritize = [
        (get_priority(name), name) for name in platform_names]
    order_prioritize.sort()
    out_list = [i[1] for i in order_prioritize]
    return out_list


def find_header(MATRIX, category):
    # Returns a header from MATRIX that matches this category or None.

    # XXX NEED TO FIX THIS
    name_fix = {
        "entrez_ID_human" : "Entrez_ID_human",
        "entrez_ID_symbol_human" : "Entrez_symbol_human",
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
    x = os.path.split(filename)[1]
    x = x.replace(".gz", "")
    x = x.replace(".csv", "")
    x = x.replace(".annot", "")
    x = x.replace("_annot", "")
    x = x.replace("-", "_")
    version = 0
    m = re.search(r".na(\d+)",x)
    if m:
        version = m.group(1)
        x = x.replace(m.group(0),'')
    return x,version


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
    filename = None
    path = config.annot_data_affy
    assert os.path.exists(path), '%s does not exist' % path
    chip2file = {}
    for f in os.listdir(path):
        filename1 = os.path.join(path, f)
        chipname1, version = _hash_chipname(filename1)
        if chipname1 in chip2file.keys():
            if version > chip2file[chipname1][0]:
                chip2file[chipname1] = (version,filename1)
        else:
            chip2file[chipname1] = (version, filename1)
    if chipname in chip2file.keys():
        version, filename = chip2file[chipname]
    return filename


def _read_annotations_h():
    paths = []
    result = []
    root = config.psid2platform
    assert os.path.exists(root), "path %s not exist: %s" % root
    for subfolder in os.listdir(root):
        if '.DS_Store' in subfolder:
            continue
        assert os.path.isdir(os.path.join(root, subfolder))
        for platform in os.listdir(os.path.join(root, subfolder)):
            paths.append((root, subfolder, platform))
    for x in paths:
        root, subfolder, platform = x
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
        chipname, case_sensitive, gene = x
        y = _compare_annotations(annotations,gene,case_sensitive)
        number_shared_annots, only1, only2, match = y
        results.append((chipname,number_shared_annots,match))
    results.sort(key=lambda x: (x[2],x[1]),reverse=True)
    return results


def _compare_annotations(annot1, annot2, case_sensitive):
    if not case_sensitive:
        annot1 = [psid.upper() for psid in annot1]
        annot2 = [psid.upper() for psid in annot2]
    num_shared_annots = len(set(annot1).intersection(set(annot2)))
    num_annot1_only = len(set(annot1)) - num_shared_annots
    num_annot2_only = len(set(annot2)) - num_shared_annots
    share_percentage = float(num_shared_annots)/min(len(set(annot1)),len(set(annot2)))
    return num_shared_annots, num_annot1_only, num_annot2_only, share_percentage



def score_platform_of_annotations(annotations):
    platform = None
    match = None
    possible_platforms = score_platforms(annotations)
    if possible_platforms[0][1] > 0:
        platform, number_shared_annots,match = possible_platforms[0]
    return platform, match


def score_all_platforms_of_matrix(DATA):
    """return a list of (header, platform, match) we can guess"""
    # Order of list is arbitrary.
    # XXX what is match?
    chips = dict()
    for name in DATA.row_names():
        x = DATA.row_names(name)
        possible_chip, match = score_platform_of_annotations(x)
        if possible_chip:
            chips[possible_chip] = (name, match)
    order_platforms = prioritize_platforms(chips.keys())
    #if chips is empty, will return an empty list
    x = [(chips[platform][0], platform, chips[platform][1])
         for platform in order_platforms]
    return x


def score_platform_of_matrix(DATA):
    # Return the best scoring platform for this matrix.
    platform_list = score_all_platforms_of_matrix(DATA)
    if not platform_list:
        return None, 0

    # Sort by decreasing score.
    schwartz = [(-x[2], x) for x in platform_list]
    schwartz.sort()
    platform_list = [x[-1] for x in schwartz]
    
    out_platform = platform_list[0][1]
    out_platform_match = platform_list[0][2]
    return out_platform, out_platform_match


def identify_platform_of_annotations(annotations):
    platform, match = score_platform_of_annotations(annotations)
    if match == 1:
        return platform
    return None


def identify_all_platforms_of_matrix(DATA):
    """return a list of (header, platform) we can identify"""
    platform_list = score_all_platforms_of_matrix(DATA)
    result = []
    for x in platform_list:
        header, platform, match = x
        if match == 1:
            result.append((header, platform))
    return result

    
def identify_platform_of_matrix(DATA):
    platform_name, match = score_platform_of_matrix(DATA)
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
        'HumanHT_12', "illumina_humanht_12", "hsapiens_gene_ensembl", PROBE_ID,
        21),
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
        'Entrez_symbol_human', "hgnc_symbol", "hsapiens_gene_ensembl",
        GENE_SYMBOL, 26),
    Platform(
        'Entrez_symbol_mouse', "mgi_symbol", "mmusculus_gene_ensembl",
        GENE_SYMBOL, 27),
    Platform(
        'ensembl_human',"ensembl_gene_id","hsapiens_gene_ensembl",PROBE_ID, 28),
    Platform(
        'ensembl_mouse',"ensembl_gene_id","mmusculus_gene_ensembl",PROBE_ID, 29),
    ]


annotate_header = [
    'Description', 'Gene Symbol', 'Gene ID', 'Swiss-Prot ID']

affy_headers = {
    'Description' : ['Target Description'],
    'Gene Symbol' : ['Gene Symbol'],
    'Gene ID' : ['Entrez Gene'],
    'Swiss-Prot ID' : ['SwissProt']
    }
illu_headers = {
    'Description' : ['Definition'],
    'Gene Symbol' : ['Symbol'],
    'Gene ID' : ['Entrez_Gene_ID'],
    'Swiss-Prot ID' : ['swissport_id']
    }
biomart_headers = {
    'Description' : ['description'],
    'Gene Symbol' : ['hgnc_symbol','mgi_symbol'],
    'Gene ID' : ['entrezgene'],
    'Swiss-Prot ID' : ['uniprot_swissprot']
    }
