#arrayplatformlib.py
import os
import re
import arrayio
from genomicode import filelib, config

all_platforms = None

class Platform:
    def __init__(self,name,bm_attribute,bm_organism,priority):
        self.name = name
        self.bm_attribute = bm_attribute
        self.bm_organism = bm_organism
        self.priority = priority

    
def get_bm_attribute(platform):
    for one_platform in platforms:
        if one_platform.name == platform:
            bm_attribute = one_platform.bm_attribute
            return bm_attribute
    return None

def get_bm_organism(platform):
    for one_platform in platforms:
        if one_platform.name == platform:
            bm_organism = one_platform.bm_organism
            return bm_organism
    return None


def get_priority(platform):
    for one_platform in platforms:
        if one_platform.name == platform:
            priority = one_platform.priority
            return priority
    return None


def prioritize_platforms(platforms_list):
    order_prioritize = [(get_priority(platform),
                         platform) for platform in platforms_list]
    order_prioritize.sort()
    out_list = [i[1] for i in order_prioritize]
    return out_list


def hash_chipname(filename):
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
    assert os.path.exists(path),'%s does not exist'%path
    chipname=chipname.replace('_','-')
    for file in os.listdir(path):
        if chipname in file:
            filename = os.path.join(path, file)
    return filename


def chipname2filename_affy(chipname):
    filename = None
    path = config.annot_data_affy
    assert os.path.exists(path),'%s does not exist'%path
    chip2file = {}
    for file in os.listdir(path):
        filename1 = os.path.join(path, file)
        chipname1,version = hash_chipname(filename1)
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
    assert os.path.exists(root),'the %s does not exisits'%root
    for subfolder in os.listdir(root):
        if '.DS_Store' in subfolder:
            continue
        assert os.path.isdir(os.path.join(root,subfolder))
        for platform in os.listdir(os.path.join(root,subfolder)):
            paths.append((root,subfolder,platform))
    for x in paths:
        root,subfolder,platform = x
        assert subfolder in ['case_sensitive','case_insensitive']
        f = file(os.path.join(root,subfolder,platform),'r')
        text = f.readlines()
        text = [i.strip() for i in text if len(i.strip())>0]
        f.close()
        chipname = os.path.splitext(platform)[-2] #remove the '.txt'
        if subfolder == 'case_insensitive':
            result.append((chipname,False,text))
        else:
            result.append((chipname,True,text))
    return result


def read_annotations():
    global all_platforms
    if not all_platforms:
        all_platforms = _read_annotations_h()
    return all_platforms


def score_platforms(annotations):
    all_platforms = read_annotations()
    results = []
    for x in all_platforms:
        chipname, case_sensitive, gene = x
        y = compare_annotations(annotations,gene,case_sensitive)
        number_shared_annots, only1, only2, match = y
        results.append((chipname,number_shared_annots,match))
    results.sort(key=lambda x: (x[2],x[1]),reverse=True)
    return results


def compare_annotations(annot1, annot2, case_sensitive):
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
    """return a list of (header,platform,match) we can guess"""
    ids = DATA.row_names()
    chips = dict()
    for id in ids:
        x = DATA.row_names(id)
        possible_chip, match = score_platform_of_annotations(x)
        if possible_chip:
            chips[possible_chip]=(id, match)
    order_platforms = prioritize_platforms(chips.keys())
    #if chips is empty, will return an empty list
    return [(chips[platform][0],platform,chips[platform][1]) for platform in order_platforms]


def score_platform_of_matrix(DATA):
    platform_list = score_all_platforms_of_matrix(DATA)
    if not platform_list:
        return None, 0
    out_platform = platform_list[0][1]
    out_platform_match = platform_list[0][2]
    return out_platform, out_platform_match


def identify_platform_of_annotations(annotations):
    platform, match = score_platform_of_annotations(annotations)
    if match == 1:
        return platform
    return None


def identify_all_platforms_of_matrix(DATA):
    """return a list of (header,platform) we can identify"""
    platform_list = score_all_platforms_of_matrix(DATA)
    result = []
    for x in platform_list:
        header, platform, match = x
        if match == 1:
            result.append((header, platform))
    return result

    
def identify_platform_of_matrix(DATA):
    out_platform, match = score_platform_of_matrix(DATA)
    if match == 1:
        return out_platform
    return None


platforms = [Platform('HG_U95A',"affy_hg_u95a","hsapiens_gene_ensembl",1),
             Platform('HG_U95Av2',"affy_hg_u95av2","hsapiens_gene_ensembl",2),
             Platform('HG_U133_Plus_2',"affy_hg_u133_plus_2","hsapiens_gene_ensembl",3),
             Platform('HG_U133A_2',"affy_hg_u133a_2","hsapiens_gene_ensembl",4),
             Platform('HG_U133A',"affy_hg_u133a","hsapiens_gene_ensembl",5),
             Platform('HG_U133B',"affy_hg_u133b","hsapiens_gene_ensembl",6),
             Platform('Hu35KsubA',None,None,7),
             Platform('Hu35KsubB',None,None,8),
             Platform('Hu35KsubC',None,None,9),
             Platform('Hu35KsubD',None,None,10),
             Platform('Hu6800',"affy_hugenefl","hsapiens_gene_ensembl",11),
             Platform('MG_U74Av2',"affy_mg_u74av2","mmusculus_gene_ensembl",12),
             Platform('MG_U74Bv2','affy_mg_u74bv2',"mmusculus_gene_ensembl",13),
             Platform('MG_U74Cv2','affy_mg_u74cv2',"mmusculus_gene_ensembl",14),
             Platform('Mouse430_2','affy_mouse430_2',"mmusculus_gene_ensembl",15),
             Platform('Mouse430A_2','affy_mouse430a_2',"mmusculus_gene_ensembl",16),
             Platform('Mu11KsubA','affy_mu11ksuba',"mmusculus_gene_ensembl",17),
             Platform('Mu11KsubB',"affy_mu11ksubb","mmusculus_gene_ensembl",18),
             Platform('RG_U34A','affy_rg_u34a',"rnorvegicus_gene_ensembl",19),
             Platform('RAE230A','affy_rae230a',"rnorvegicus_gene_ensembl",20),
             Platform('HumanHT_12',"illumina_humanht_12","hsapiens_gene_ensembl",21),
             Platform('HumanWG_6',"illumina_humanwg_6_v3","hsapiens_gene_ensembl",22),
             Platform('MouseRef_8',"illumina_mousewg_6_v2","mmusculus_gene_ensembl",23),
             Platform('Entrez_ID_human',"entrezgene","hsapiens_gene_ensembl",24),
             Platform('Entrez_ID_mouse',"entrezgene","mmusculus_gene_ensembl",25),
             Platform('Entrez_symbol_human',"hgnc_symbol","hsapiens_gene_ensembl",26),
             Platform('Entrez_symbol_mouse',"mgi_symbol","mmusculus_gene_ensembl",27),
             ]


annotate_header = ['Description', 'Gene Symbol', 'Gene ID', 'Swiss-Prot ID']
affy_headers = {'Description': 'Target Description','Gene Symbol':'Gene Symbol',
               'Gene ID':'Entrez Gene','Swiss-Prot ID': 'SwissProt'}
illu_headers = {'Description': 'Definition', 'Gene Symbol':'Symbol',
               'Gene ID':'Entrez_Gene_ID','Swiss-Prot ID':'swissport_id'}
biomart_headers = {'Description': 'description','Gene Symbol':'hgnc_symbol',
                  'Gene ID':'entrezgene', 'Swiss-Prot ID': 'uniprot_swissprot'}
