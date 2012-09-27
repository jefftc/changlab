#arrayplatformlib.py

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
    order_prioritize = [(get_priority(platform),platform) for platform in platforms_list]
    order_prioritize.sort()
    out_list = [i[1] for i in order_prioritize]
    return out_list


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

