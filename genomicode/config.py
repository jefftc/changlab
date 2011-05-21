### CONFIGURATION

import os

opj = os.path.join


HOME = os.environ["HOME"]
PROJECTS = opj(HOME, "projects")

SEARCH_PATH = [
    opj(HOME, "bin"),
    "/usr/local/bin",
    "/usr/bin",
    ]

# patserfns
MOTIFSEARCH = opj(PROJECTS, "motifsearch")
patserfns_PATSER_BIN = opj(MOTIFSEARCH, "patser/patser-v3b")
patserfns_ALPHABET_FILE = opj(MOTIFSEARCH, "alphabet.prom_500_100")

motiffns_JASPAR_DB = opj(MOTIFSEARCH, "data/JASPAR_CORE_2008")
motiffns_TRANSFAC_DB = opj(MOTIFSEARCH, "data/transfac.12_1.pfm")
motiffns_JASPAR_INFO = opj(MOTIFSEARCH, "jasparid2info")
motiffns_TRANSFAC_INFO = opj(MOTIFSEARCH, "matid2info")
motiffns_MATID2LENGTH = opj(MOTIFSEARCH, "matid2length")



# genomefns
SCGENOME = opj(PROJECTS, "scgenome")
genomefns_RA_CHROM_HG18 = opj(SCGENOME, "data/data.hg18/ra_chrom")
genomefns_RA_CHROM_MM9 = opj(SCGENOME, "data/data.mm9/ra_chrom")


# tss
TSS = opj(PROJECTS, "tss")
gene_HG18 = opj(TSS, "data/hg18.knownGene.110312.txt")


# primer3fns
primer3fns_PRIMER3_BIN = "primer3_core"
