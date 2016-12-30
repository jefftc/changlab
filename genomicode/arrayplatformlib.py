"""

Functions:
score_annotations        Given annotations, score for each platform.
score_matrix

find_header              Given a category, try to find the header of a matrix.
categorize_headers
chipname2filename
normalize_id             Normalize a gene name for comparison.
normalize_ids            Normalize a list of IDs.

prioritize_platforms
find_platform_by_name    Return a Platform object with a given name.
get_bm_attribute
get_bm_organism
get_priority


Classes:
Platform      A description for a probe on a microarray.
Score         Describes the similarity between annotations and platform.

Variables:
PLATFORMS     List of all known platforms.

Constants:
PROBE_ID
GENE_ID
GENE_SYMBOL
DESCRIPTION

"""
# _read_annotations
# _normalize_annot
#
# _hash_chipname
# _chipname2filename_illu
# _chipname2filename_affy

# TODO:
# - What is a chipname?
#   What is a platform?  How are the names determined?
#   score_all_platforms_of_matrix returns entrez_ID_symbol_human.
#   Looks like a chipname.  Doesn't match the name of any Platform.
# - What are the headers at the bottom?
# - Merge with arrayannot.py



import os
import re
from genomicode import filelib, config


PROBE_ID, GENE_ID, GENE_SYMBOL, DESCRIPTION = range(4)


class Platform:
    def __init__(self, name, bm_attribute, bm_organism, category, priority):
        # our description of platform, corresponds to files
        # in genomidata/psid2platform, also matches annotation files in
        # /data/genomidata/affymetrix and /data/genomidata/illumina
        self.name = name
        # bm_organism   Which BioMart organism to use.  XXX how to get list?
        # bm_attribute  Which attribute from this mart.
        #
        # To list valid attributes:
        # library(biomaRt)
        # x <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl",
        #   host="www.ensembl.org")
        # y <- listAttributes(x)
        # y$name
        self.bm_attribute = bm_attribute
        self.bm_organism = bm_organism
        # PROBE_ID, GENE_ID, GENE_SYMBOL, DESCRIPTION
        self.category = category
        # order of platform priority, lower number means higher priority
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
    """highest priority to lowest"""
    # lower number is higher priority    
    order_prioritize = [
        (get_priority(name), name) for name in platform_names]
    order_prioritize.sort()
    prioritized_names = [i[1] for i in order_prioritize]
    return prioritized_names


def categorize_headers(MATRIX, remove_version=True):
    # Return a dictionary of category -> header.
    
    #name_fix = {
    #    "entrez_ID_human" : "Entrez_ID_human",
    #    "entrez_ID_symbol_human" : "Entrez_Symbol_human",
    #    }
    scores = score_matrix(MATRIX, remove_version=remove_version)

    # Filter based on a minimum score.
    scores = [x for x in scores if x.max_score >= 0.5]
    
    # Order by increasing platform priority, decreasing score.
    # This is not quite right.  Really want decreasing score, as long
    # as category is not too bad.
    x = [x.platform_name for x in scores]
    all_platforms = {}.fromkeys(x)
    platform_order = prioritize_platforms(all_platforms)
    priority = [platform_order.index(x.platform_name) for x in scores]
    x = [(priority[i], -scores[i].max_score, scores[i])
         for i in range(len(scores))]
    x.sort()
    scores = [x[-1] for x in x]

    cat2header = {}
    for score in scores:
        platform = find_platform_by_name(score.platform_name)
        if not platform:
            continue
        if platform.category in cat2header:
            continue
        cat2header[platform.category] = score.header
    return cat2header


def find_header(MATRIX, category, remove_version=True):
    # Returns a header from MATRIX that matches this category or None.

    cat2header = categorize_headers(MATRIX, remove_version=remove_version)
    return cat2header.get(category)


def _hash_chipname(filename):
    # RG_U34A_annot.csv.gz                   RG_U34A
    # Rat230_2.na26.annot.csv.gz             Rat230_2
    # MoEx_1_0_st_v1.1.mm9.probeset          MoEx_1_0_st_v1  [151015 no v1?]
    # HuGene-2_0-st-v1.na35.hg19.probeset    HuGene_2_0_st
    # HuGene-2_0-st-v1.na35.hg19.transcript  HuGene_2_0_st
    REMOVE = [
        ".gz", ".csv", ".annot", "_annot",
        ".mm9", ".probeset", ".transcript", ".1", ".hg19",
        "-v1",
        ]
    x = os.path.split(filename)[1]
    for r in REMOVE:
        x = x.replace(r, "")
    x = x.replace("-", "_")
    version = 0
    m = re.search(r".na(\d+)",x)
    if m:
        version = int(m.group(1))
        x = x.replace(m.group(0), "")
    return x, version


def chipname2filename(chipname):
    filename = _chipname2filename_affy(chipname)
    if not filename:
        filename = _chipname2filename_illu(chipname)
    return filename


def _chipname2filename_illu(chipname):
    filename = None
    path = config.annot_data_illu
    assert os.path.exists(path), "%s does not exist" % path
    chipname = chipname.replace("_", "-")
    for f in os.listdir(path):
        if chipname in f:
            filename = os.path.join(path, f)
    return filename


def _chipname2filename_affy(chipname):
    path = config.annot_data_affy
    assert os.path.exists(path), "%s does not exist" % path

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
    assert os.path.exists(root), "Path not found: %s" % root
    for subfolder in os.listdir(root):
        if ".DS_Store" in subfolder:
            continue
        x = os.path.join(root, subfolder)
        assert os.path.isdir(x), "Not directory: %s" % x
        for platform in os.listdir(x):
            paths.append((root, subfolder, platform))
    all_platform = []
    for platform in PLATFORMS:
        all_platform.append(platform.name)
    for x in paths:
        root, subfolder, platform = x
        assert platform.lower().endswith(".txt")
        platform_name = platform[:-4] # get rid of .txt
        assert platform_name in all_platform, (
            "%s is not a platform in PLATFORMS" % platform_name)
        assert subfolder in ["case_sensitive","case_insensitive"]
        f = file(os.path.join(root,subfolder,platform),"r")
        text = f.readlines()
        text = [i.strip() for i in text if len(i.strip())>0]
        f.close()
        chipname = os.path.splitext(platform)[-2]  # remove the ".txt"
        if subfolder == "case_insensitive":
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


class Score:
    def __init__(
        self, platform_name, header,
        min_score, max_score, mine_only, platform_only, shared):
        # header may be None.  Only used if scoring a Matrix.
        self.platform_name = platform_name
        self.header = header
        self.min_score = min_score
        self.max_score = max_score
        self.mine_only = list(mine_only)
        self.platform_only = list(platform_only)
        self.shared = list(shared)


# list of (chip_name, case_sensitive, dict of IDs)
ALL_PLATFORMS_CLEAN = None
def score_annotations(annotations, min_score=None, remove_version=False):
    # Return list of Scores.  Sorted by decreasing score.
    # Optimized version of the function.
    global ALL_PLATFORMS_CLEAN

    if ALL_PLATFORMS_CLEAN is None:
        ALL_PLATFORMS_CLEAN = []
        all_platforms = _read_annotations()
        for x in all_platforms:
            chipname, case_sensitive, ids = x
            ids_clean = normalize_ids(
                ids, case_insensitive=(not case_sensitive),
                remove_version_number=remove_version)
            ids_clean = frozenset(ids_clean)
            x = chipname, case_sensitive, ids_clean
            ALL_PLATFORMS_CLEAN.append(x)

    # HACK: Remove version numbers from ENSEMBL IDs.
    for i in range(len(annotations)):
        x = annotations[i]
        if x.startswith("ENS"):
            j = x.find(".")
            if j >= 0:
                x = x[:j]
        annotations[i] = x

    annots_clean_cs = normalize_ids(annotations, case_insensitive=False)
    annots_clean_ci = normalize_ids(annotations, case_insensitive=True)
    annots_clean_cs = frozenset(annots_clean_cs)
    annots_clean_ci = frozenset(annots_clean_ci)

    results = []
    for x in ALL_PLATFORMS_CLEAN:
        chipname, case_sensitive, ids = x

        annots1 = annots_clean_cs
        if not case_sensitive:
            annots1 = annots_clean_ci
        annots2 = ids

        shared = [x for x in annots1 if x in annots2]
        annots1_only = [x for x in annots1 if x not in annots2]
        annots2_only = [x for x in annots2 if x not in annots1]
        num_shared_annots = len(shared)
        num_annots1_only = len(annots1_only)
        num_annots2_only = len(annots2_only)
        
        perc_shared_max = 0
        perc_shared_min = 0
        if annots1 and annots2:
            x = float(num_shared_annots)/min(len(annots1), len(annots2))
            perc_shared_max = x
            x = float(num_shared_annots)/max(len(annots1), len(annots2))
            perc_shared_min = x

        x = Score(
            chipname, None,
            perc_shared_min, perc_shared_max, annots1, annots2, shared)
        results.append(x)

    if min_score is not None:
        results = [x for x in results if x.max_score >= min_score]
        
    # Sort by decreasing score.
    schwartz = [(-x.max_score, -x.min_score, x.platform_name, x)
                for x in results]
    schwartz.sort()
    results = [x[-1] for x in schwartz]
    return results


#score_platforms = _score_platforms_fast
## def score_platforms(annotations):
##     all_platforms = _read_annotations()
##     results = []
##     for x in all_platforms:
##         chipname, case_sensitive, ids = x
##         y = _compare_annotations(annotations, ids, case_sensitive)
##         number_shared_annots, only1, only2, match = y
##         results.append((chipname, number_shared_annots, match))
##     results.sort(key=lambda x: (x[2], x[1]),reverse=True)
##     return results


def normalize_id(gene_id, case_insensitive=False, delimiter=None,
                 remove_version_number=False):
    # Return a list of this ID cleaned up.  gene_id might be split
    # into multiple IDs by the delimiter.
    x = [gene_id]
    if delimiter:
        x = gene_id.split(delimiter)
    if case_insensitive:
        x = [x.upper() for x in x]
    # Clean up whitespace.
    x = [x.strip() for x in x]
    # No empty IDs.
    x = [x for x in x if x]
    # Ignore "---".
    x = [x for x in x if x != "---"]
    if remove_version_number:
        # Hack: Remove version numbers from RefSeq or ENTREZ IDs.
        x = [_remove_version(x) for x in x]
    # No duplicates.
    x = {}.fromkeys(x).keys()
    return x


def normalize_ids(
    gene_ids, case_insensitive=False, delimiter=None,
    remove_version_number=False):
    norm_ids = []
    for x in gene_ids:
        x = normalize_id(
            x, case_insensitive=case_insensitive, delimiter=delimiter,
            remove_version_number=remove_version_number)
        norm_ids.extend(x)
    return norm_ids


def _remove_version(refseq_id):
    i = refseq_id.find(".")
    if i >= 0:
        refseq_id = refseq_id[:i]
    return refseq_id


## def _normalize_annot(annot, case_insensitive=False, delim=None):
##     x = annots
##     if delim:
##         parsed = []
##         for x in x:
##             parsed.extend(x.split(delim))
##         x = parsed
##     if case_insensitive:
##         x = [psid.upper() for psid in x]
##     x = [x.strip() for x in x]   # Ignore differences in whitespace.
##     x = [x for x in x if x]      # No blank annotations.
##     x = {}.fromkeys(x)           # No duplicates
##     return x


## def _compare_annotations(annots1, annots2, case_sensitive):
##     # Return tuple of: num_shared, num_annots1_only, num_annots2_only,
##     # perc_shared.
##     if not case_sensitive:
##         annots1 = [psid.upper() for psid in annots1]
##         annots2 = [psid.upper() for psid in annots2]
##     # Ignore differences in whitespace.
##     annots1 = [x.strip() for x in annots1]
##     annots2 = [x.strip() for x in annots2]
##     # No blank annotations.
##     annots1 = [x for x in annots1 if x]
##     annots2 = [x for x in annots2 if x]
##     # No duplicate annotations.
##     annots1 = {}.fromkeys(annots1)
##     annots2 = {}.fromkeys(annots2)

##     x = [x for x in annots1 if x in annots2]
##     num_shared_annots = len(x)
##     num_annots1_only = len(annots1) - num_shared_annots
##     num_annots2_only = len(annots2) - num_shared_annots
##     perc_shared = 0
##     if annots1 and annots2:
##         perc_shared = float(num_shared_annots)/min(len(annots1), len(annots2))
##     x = num_shared_annots, num_annots1_only, num_annots2_only, perc_shared
##     return x


## def score_platform_of_annotations(annotations):
##     # Tuple of (platform, match score) or None if nothing matches.
##     scores = score_platforms(annotations)
##     assert scores
##     (platform, num_my_annots, num_platform_annots, num_shared_annots,
##      score_min, score_max) = scores[0]
##     if score_max <= 0:
##         return None
##     return platform, score_max


## def _parse_matrix_annotations(annots, delim):
##     # No blank annotations.
##     x = annots
##     if delim:
##         parsed = []
##         for x in x:
##             x = x.split(delim)
##             parsed.extend(x)
##         x = parsed
##     x = [x.strip() for x in x]
##     x = [x for x in x if x]
##     return x


def _score_matrix_h(DATA, annot_delim, min_score, remove_version):
    assert min_score is None or (min_score >= 0 and min_score <= 1.0)
    
    results = []  # List of Score objects.
    for header in DATA.row_names():
        x = DATA.row_names(header)
        annots = normalize_ids(
            x, delimiter=annot_delim, remove_version_number=remove_version)
        x = score_annotations(annots, min_score=min_score)
        for i in range(len(x)):
            x[i].header = header
        results.extend(x)

    # Sort by decreasing score.
    schwartz = [(-x.max_score, -x.min_score, x.platform_name, x)
                for x in results]
    schwartz.sort()
    results = [x[-1] for x in schwartz]
    return results


SCORE_MATRIX_CACHE = {}
def score_matrix(
    DATA, annot_delim=None, min_score=None, remove_version=False):
    """Return a list of Score objects.  List is ordered by decreasing
    score.

    """
    global SCORE_MATRIX_CACHE
    if len(SCORE_MATRIX_CACHE) > 5:
        SCORE_MATRIX_CACHE = {}
    key = id(DATA), annot_delim, min_score, remove_version
    if key not in SCORE_MATRIX_CACHE:
        value = _score_matrix_h(DATA, annot_delim, min_score, remove_version)
        SCORE_MATRIX_CACHE[key] = value
    return SCORE_MATRIX_CACHE[key]


## def score_platform_of_matrix(DATA, annot_delim=None):
##     # Return the best scoring platform for this matrix.
##     platforms = score_all_platforms_of_matrix(DATA, annot_delim=annot_delim)
##     if not platforms:
##         return None, 0

##     # Sort by decreasing score.
##     schwartz = [(-x[2], x) for x in platforms]
##     schwartz.sort()
##     platforms = [x[-1] for x in schwartz]

##     header, platform, score = platforms[0]
##     return platform, score


## def identify_platform_of_annotations(annotations):
##     x = score_platform_of_annotations(annotations)
##     if not x:
##         return None
##     platform, match = x
##     if match == 1:
##         return platform
##     return None


## def identify_all_platforms_of_matrix(DATA, annot_delim=None):
##     """return a list of (header, platform name) we can identify"""
##     platforms = score_all_platforms_of_matrix(DATA, annot_delim=annot_delim)
##     result = []
##     for x in platforms:
##         header, platform, match = x
##         if match == 1:
##             result.append((header, platform))
##     return result


## def identify_platform_of_matrix(DATA, annot_delim=None):
##     x = score_platform_of_matrix(DATA, annot_delim=annot_delim)
##     platform_name, match = x
##     if match == 1:
##         return platform_name
##     return None


PLATFORMS = [
    Platform("HG_U95A", "affy_hg_u95a", "hsapiens_gene_ensembl", PROBE_ID, 1),
    Platform(
        "HG_U95Av2", "affy_hg_u95av2", "hsapiens_gene_ensembl", PROBE_ID, 2),
    Platform(
        "HG_U133_Plus_2", "affy_hg_u133_plus_2", "hsapiens_gene_ensembl",
        PROBE_ID, 3),
    Platform(
        "HG_U133A_2", "affy_hg_u133a_2", "hsapiens_gene_ensembl", PROBE_ID, 4),
    Platform(
        "HG_U133A", "affy_hg_u133a", "hsapiens_gene_ensembl", PROBE_ID, 5),
    Platform(
        "HG_U133B", "affy_hg_u133b", "hsapiens_gene_ensembl", PROBE_ID, 6),
    Platform("Hu35KsubA", None, "hsapiens_gene_ensembl", PROBE_ID, 7),
    Platform("Hu35KsubB", None, "hsapiens_gene_ensembl", PROBE_ID, 8),
    Platform("Hu35KsubC", None, "hsapiens_gene_ensembl", PROBE_ID, 9),
    Platform("Hu35KsubD", None, "hsapiens_gene_ensembl", PROBE_ID, 10),
    Platform("Hu6800", "affy_hugenefl", "hsapiens_gene_ensembl", PROBE_ID, 11),
    Platform(
        "MG_U74Av2", "affy_mg_u74av2", "mmusculus_gene_ensembl", PROBE_ID, 12),
    Platform(
        "MG_U74Bv2","affy_mg_u74bv2", "mmusculus_gene_ensembl", PROBE_ID, 13),
    Platform(
        "MG_U74Cv2","affy_mg_u74cv2", "mmusculus_gene_ensembl", PROBE_ID, 14),
    Platform(
        "Mouse430_2","affy_mouse430_2", "mmusculus_gene_ensembl", PROBE_ID,
        15),
    Platform(
        "Mouse430A_2","affy_mouse430a_2", "mmusculus_gene_ensembl", PROBE_ID,
        16),
    Platform(
        "Mu11KsubA","affy_mu11ksuba", "mmusculus_gene_ensembl", PROBE_ID, 17),
    Platform(
        "Mu11KsubB", "affy_mu11ksubb", "mmusculus_gene_ensembl", PROBE_ID, 18),
    Platform(
        "RG_U34A","affy_rg_u34a", "rnorvegicus_gene_ensembl", PROBE_ID, 19),
    Platform(
        "RAE230A","affy_rae230a", "rnorvegicus_gene_ensembl", PROBE_ID, 20),
    Platform(
        "HumanHT_12", "illumina_humanht_12_v4", "hsapiens_gene_ensembl",
        PROBE_ID, 21),
    Platform(
        "HumanWG_6", "illumina_humanwg_6_v3", "hsapiens_gene_ensembl",
        PROBE_ID, 22),
    Platform(
        "MouseRef_8", "illumina_mouseref_8_v2", "mmusculus_gene_ensembl",
        PROBE_ID, 23),
    Platform(
        "Entrez_ID_human", "entrezgene", "hsapiens_gene_ensembl", GENE_ID, 24),
    Platform(
        "Entrez_ID_mouse", "entrezgene", "mmusculus_gene_ensembl", GENE_ID,
        25),
    Platform(
        "Entrez_Symbol_human", "hgnc_symbol", "hsapiens_gene_ensembl",
        GENE_SYMBOL, 26),
    Platform(
        "Entrez_Symbol_mouse", "mgi_symbol", "mmusculus_gene_ensembl",
        GENE_SYMBOL, 27),
    Platform(
        "Ensembl_human", "ensembl_gene_id", "hsapiens_gene_ensembl",
        PROBE_ID, 28),
    Platform(
        "Ensembl_mouse", "ensembl_gene_id", "mmusculus_gene_ensembl",
        PROBE_ID, 29),
    Platform(
        "UCSC_human_hg19_kg7", None, "hsapiens_gene_ensembl", GENE_ID, 30),
    Platform(
        "UCSC_human_hg19_kg6", None, "hsapiens_gene_ensembl", GENE_ID, 31),
    Platform(
        "UCSC_human_hg19_kg5", None, "hsapiens_gene_ensembl", GENE_ID, 32),
    Platform(
        "UCSC_human_hg38_kg8", None, "hsapiens_gene_ensembl", GENE_ID, 33),
    Platform(
        "Agilent_Human1A", None, "hsapiens_gene_ensembl", GENE_ID, 34),
    Platform(
        "HumanHT_12_control", None, "hsapiens_gene_ensembl", GENE_ID, 35),
    Platform(
        "MouseRef_8_control", None, "hsapiens_gene_ensembl", GENE_ID, 36),
    Platform(
        "UCSC_human_hg38_kg8", None, "hsapiens_gene_ensembl", GENE_ID, 37),
    Platform(
        "RefSeq_protein_mouse", "refseq_peptide", "mmusculus_gene_ensembl",
        PROBE_ID, 38),
    Platform(
        "RefSeq_predicted_protein_mouse", "refseq_peptide_predicted",
        "mmusculus_gene_ensembl", PROBE_ID, 39),
    Platform(
        "RefSeq_protein_human", "refseq_peptide", "hsapiens_gene_ensembl",
        PROBE_ID, 40),
    Platform(
        "RefSeq_predicted_protein_human", "refseq_peptide_predicted",
        "hsapiens_gene_ensembl", PROBE_ID, 41),
    Platform(
        "RefSeq_mRNA_mouse", "refseq_mrna", "mmusculus_gene_ensembl",
        PROBE_ID, 42),
    Platform(
        "RefSeq_predicted_mRNA_mouse", "refseq_mrna_predicted",
        "mmusculus_gene_ensembl", PROBE_ID, 43),
    Platform(
        "RefSeq_ncRNA_mouse", "refseq_ncrna", "mmusculus_gene_ensembl",
        PROBE_ID, 44),
    Platform(
        "RefSeq_predicted_ncRNA_mouse", "refseq_ncrna_predicted",
        "mmusculus_gene_ensembl", PROBE_ID, 45),
    Platform(
        "RefSeq_mRNA_human", "refseq_mrna", "hsapiens_gene_ensembl",
        PROBE_ID, 46),
    Platform(
        "RefSeq_predicted_mRNA_human", "refseq_mrna_predicted",
        "hsapiens_gene_ensembl", PROBE_ID, 47),
    Platform(
        "RefSeq_ncRNA_human", "refseq_ncrna", "hsapiens_gene_ensembl",
        PROBE_ID, 48),
    Platform(
        "RefSeq_predicted_ncRNA_human", "refseq_ncrna_predicted",
        "hsapiens_gene_ensembl", PROBE_ID, 49),
    Platform(
        "UCSC_mouse_mm10_kg7", None, "mmusculus_gene_ensembl", GENE_ID, 50),
    Platform(
        "Entrez_Description_human", None, "hsapiens_gene_ensembl",
        DESCRIPTION, 51),
    Platform(
        "Entrez_Description_mouse", None, "mmusculus_gene_ensembl",
        DESCRIPTION, 52),
    # BioMart only has v3, not v2.
    Platform("HumanRef_8_v2", None, "hsapiens_gene_ensembl", PROBE_ID, 53),
    Platform("HumanRef_8_v3", None, "hsapiens_gene_ensembl", PROBE_ID, 54),

    Platform(
        "MouseWG_6", "illumina_mousewg_6_v2", "mmusculus_gene_ensembl",
        PROBE_ID, 55),
    ]



## annotate_header = [#The added header names when annotate a file
##     "Description", "Gene Symbol", "Gene ID", "Swiss-Prot ID"]

## affy_headers = {# Header names in affymetrix annotation
##                 # files in /data/genomidata/affymetrix
##     "Description" : ["Target Description"],
##     "Gene Symbol" : ["Gene Symbol"],
##     "Gene ID" : ["Entrez Gene"],
##     "Swiss-Prot ID" : ["SwissProt"]
##     }
## illu_headers = {# Header names in illumina annotation
##                 # files in /data/genomidata/illumina
##     "Description" : ["Definition"],
##     "Gene Symbol" : ["Symbol"],
##     "Gene ID" : ["Entrez_Gene_ID"],
##     "Swiss-Prot ID" : ["swissport_id"]
##     }
## biomart_headers = {# Header names in biomart attributesL
##     "Description" : ["description"],
##     "Gene Symbol" : ["hgnc_symbol","mgi_symbol"],
##     "Gene ID" : ["entrezgene"],
##     "Swiss-Prot ID" : ["uniprot_swissprot"]
##     }

platform_to_GSEA_chipname={# corresponds to chipname in genepattern GSEA module
    "Agilent_Human1A" : "Agilent_Human1A",
    "HG_U133A_2" : "HG_U133A_2",
    "HG_U133A" : "HG_U133A",
    "HG_U133B" : "HG_U133B",
    "HG_U133_Plus_2" : "HG_U133_Plus_2",
    "HG_U95Av2" : "HG_U95Av2",
    "Hu35KsubA" : "Hu35KsubA",
    "Hu35KsubB" : "Hu35KsubB",
    "Hu35KsubC" : "Hu35KsubC",
    "Hu35KsubD" : "Hu35KsubD",
    "Hu6800" : "Hu6800",
    "MG_U74Av2" : "MG_U74Av2",
    "MG_U74Bv2" : "MG_U74Bv2",
    "MG_U74Cv2" : "MG_U74Cv2",
    "Mouse430_2" : "Mouse430_2",
    "Mouse430A_2" : "Mouse430A_2",
    "Mu11KsubA" : "Mu11KsubA",
    "Mu11KsubB" : "Mu11KsubB",
    "RAE230A" : "RAE230A",
    "RG_U34A" : "RG_U34A",
    "HumanWG_6" : "ilmn_HumanWG_6_V3_0_R3_11282955_A",
    "HumanHT_12" : "ilmn_HumanHT_12_V3_0_R3_11283641_A",
    "MouseRef_8" : "ilmn_MouseRef_8_V2_0_R3_11278551_A"
    }
