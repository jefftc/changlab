#arrayplatformlib.py
import os
import re
import arrayio
from genomicode import filelib, config


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


def identify_platform_of_matrix(DATA):
    platform_list = identify_all_platforms_of_matrix(DATA)
    if not platform_list:
        return None
    out_platform = platform_list[0][1]
    return out_platform


def identify_all_platforms_of_matrix(DATA):
    """return a list of (header,platform) we can identify"""
    ids = DATA.row_names()
    chips = dict()
    for id in ids:
        x = DATA.row_names(id)
        possible_chip = identify_platform_of_annotations(x)
        if possible_chip:
            chips[possible_chip]=id
    order_platforms = prioritize_platforms(chips.keys())
    #if chips is empty, will return an empty list
    return [(chips[platform],platform) for platform in order_platforms]

    
def identify_platform_of_annotations(annotations):
    platform2annot_cs = {}  # chip -> psid -> 1
    platform2annot_ci = {}  # chip -> psid -> 1
    paths = []
    possible_chips = []
    annotations = [i for i in annotations if len(i)>0]
    upannotations = [i.upper() for i in annotations]
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
            for psid in text:
                platform2annot_cs.setdefault(chipname, {})[psid.upper()] = 1
        else:
             for psid in text:
                platform2annot_ci.setdefault(chipname, {})[psid] = 1
                
    for chip in platform2annot_cs:
        for psid in upannotations:
            if psid not in platform2annot_cs[chip]:
                break
        else:
            possible_chips.append(chip)

    for chip in platform2annot_ci:
        for psid in annotations:
            if psid not in platform2annot_ci[chip]:
                break
        else:
            possible_chips.append(chip)  
    
    if not possible_chips:
        return None
    #combine the dict for both case_sensitive and case_insensitive into one dict
    platform2annot = platform2annot_cs.copy()
    for chip in platform2annot_ci:
        platform2annot[chip]=platform2annot_ci[chip]

    # Sort the chips by size, from smallest to largest.
    schwartz = [(len(platform2annot[chip]), chip) for chip in possible_chips]
    schwartz.sort()
    possible_chips = [x[-1] for x in schwartz]
    # Choose the smallest chip that contains all these probe sets.
    return possible_chips[0]


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
             Platform('MouseRef_8',"illumina_mousewg_6_v2","mmusculus_gene_ensembl",22),
             Platform('Entrez_ID_human',"entrezgene","hsapiens_gene_ensembl",23),
             Platform('Entrez_ID_mouse',"entrezgene","mmusculus_gene_ensembl",24),
             Platform('Entrez_symbol_human',"hgnc_symbol","hsapiens_gene_ensembl",25),
             Platform('Entrez_symbol_mouse',"mgi_symbol","mmusculus_gene_ensembl",26),
             ]


##annotate_header = ['Description', 'Gene_symbol', 'Gene_ID', 'Swissprot_ID']
##affy_headers = {'Description': 'Target Description','Gene_symbol':'Gene Symbol',
##               'Gene_ID':'Entrez Gene','Swissprot_ID': 'SwissProt'}
##illu_headers = {'Description': 'Definition','Gene_symbol':'Symbol',
##               'Gene_ID':'Entrez_Gene_ID','Swissprot_ID':'swissport_id'}
##biomart_headers = {'Description': 'description', 'Gene_symbol':'hgnc_symbol',
##                  'Gene_ID':'entrezgene', 'Swissprot_ID': 'uniprot_swissprot'}

annotate_header = ['Gene Symbol', 'Gene ID', 'Swiss-Prot ID']
affy_headers = {'Gene Symbol':'Gene Symbol',
               'Gene ID':'Entrez Gene','Swiss-Prot ID': 'SwissProt'}
illu_headers = {'Gene Symbol':'Symbol',
               'Gene ID':'Entrez_Gene_ID','Swiss-Prot ID':'swissport_id'}
biomart_headers = {'Gene Symbol':'hgnc_symbol',
                  'Gene ID':'entrezgene', 'Swiss-Prot ID': 'uniprot_swissprot'}
