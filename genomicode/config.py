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

# patser
MOTIFSEARCH = opj(PROJECTS, "motifsearch")
patser_PATSER_BIN = opj(MOTIFSEARCH, "patser/patser-v3b")
patser_ALPHABET_FILE = opj(MOTIFSEARCH, "alphabet.prom_500_100")

motiflib_JASPAR_DB = opj(MOTIFSEARCH, "data/JASPAR_CORE_2008")
motiflib_TRANSFAC_DB = opj(MOTIFSEARCH, "data/transfac.12_1.pfm")
motiflib_JASPAR_INFO = opj(MOTIFSEARCH, "jasparid2info")
motiflib_TRANSFAC_INFO = opj(MOTIFSEARCH, "matid2info")
motiflib_MATID2LENGTH = opj(MOTIFSEARCH, "matid2length")



# genomelib
SCGENOME = opj(PROJECTS, "scgenome")
genomelib_RA_CHROM_HG18 = opj(SCGENOME, "data/data.hg18/ra_chrom")
genomelib_RA_CHROM_MM9 = opj(SCGENOME, "data/data.mm9/ra_chrom")


# tss
TSS = opj(PROJECTS, "tss")
gene_HG18 = opj(TSS, "data/hg18.knownGene.110312.txt")


# primer3
primer3_PRIMER3_BIN = "primer3_core"
