"""

Functions:
convert_gene_ids

"""
# _convert_gene_ids_biomart
# _convert_gene_ids_local
# _convert_entrez_symbol_to_entrez
# _find_entrez_gene
#
# _clean_genes_for_biomart
# _remove_dups
# _start_R

GLOBAL_R = None
def _start_R():
    global GLOBAL_R
    from genomicode import jmath

    if GLOBAL_R is None:
        R = jmath.start_R()
        R('library(biomaRt)')
        GLOBAL_R = R
    return GLOBAL_R


def convert_gene_ids(
    gene_ids, in_platform, out_platform, in_delim=" /// ", out_delim=" /// ",
    keep_dups=True, keep_emptys=True, no_na=True):
    # Return a list of the output IDs, parallel to gene_ids.  If a
    # gene ID to multiple output IDs, the output IDs will be separated
    # by out_delim.  If it is missing, then the output will be an
    # empty string.
    # in_platform and out_platform are the names of the platforms.
    from genomicode import arrayplatformlib as apl

    # Make a cleaned up version of the gene_ids to convert.
    remove_version = False
    if in_platform.lower().find("refseq") >= 0 or \
       in_platform.lower().find("ensembl") >= 0:
        remove_version = True
    x = apl.normalize_ids(
        gene_ids, delimiter=in_delim, remove_version_number=remove_version)
    # No duplicates.
    x = {}.fromkeys(x).keys()
    gene_ids_c = x

    in2out = None
    if in_platform == "Entrez_Symbol_human" and \
       out_platform in ("Entrez_Symbol_human", "Entrez_ID_human"):
        in2out = _convert_entrez_symbol_to_entrez(
            gene_ids_c, out_platform, "9606")
    if in2out is None:
        in2out = _convert_gene_ids_local(in_platform, out_platform)
    if in2out is None:
        in2out = _convert_gene_ids_biomart(
            gene_ids_c, in_platform, out_platform, no_na)
    assert in2out, "I could not convert %s to %s" % (in_platform, out_platform)

    # Make a parallel list of the output IDs.
    output_ids = []
    for gene_id in gene_ids:
        in_ids = apl.normalize_id(
            gene_id, delimiter=in_delim, remove_version_number=remove_version)
        out_ids = []
        for x in in_ids:
            y = in2out.get(x, [""])
            out_ids.extend(y)
        if not keep_emptys:
            out_ids = [x for x in out_ids if x]
        if not keep_dups:
            out_ids = _remove_dups(out_ids)
        x = out_delim.join(out_ids)
        output_ids.append(x)
    return output_ids


def _clean_genes_for_biomart(gene_ids):
    gene_ids = [i.replace("'",'') for i in gene_ids] #modify "2'-PDE"
    return gene_ids


def _convert_gene_ids_biomart(gene_ids, in_platform, out_platform, no_na):
    # Maximum number of genes to request at a time.
    MAX_GENES = 25000

    in2out = {}
    while gene_ids:
        batch = gene_ids[:MAX_GENES]
        gene_ids = gene_ids[MAX_GENES:]
        
        x = _convert_gene_ids_biomart_h(
            batch, in_platform, out_platform, no_na)
        assert x is not None, "Cannot convert: %s %s" % (
            in_platform, out_platform)
        in2out.update(x)
    return in2out


def _convert_gene_ids_biomart_h(gene_ids, in_platform, out_platform, no_na):
    # Return a dictionary of gene_id -> list of converted_ids, or None
    # if these platforms cannot be converted.
    from genomicode import jmath
    from genomicode import arrayplatformlib

    if not gene_ids:
        return {}
    R_fn, R_var = jmath.R_fn, jmath.R_var

    # An attribute is the biomart name for the platform.
    in_attribute = arrayplatformlib.get_bm_attribute(in_platform)
    out_attribute = arrayplatformlib.get_bm_attribute(out_platform)
    #assert in_attribute, "Bad platform: %s" % in_platform
    #assert out_attribute, "Bad platform: %s" % out_platform
    if not in_attribute or not out_attribute:
        return None
    in_mart = arrayplatformlib.get_bm_organism(in_platform)
    out_mart = arrayplatformlib.get_bm_organism(out_platform)
    assert in_mart, "No bm organism for platform: %s" % in_platform

    gene_ids = _clean_genes_for_biomart(gene_ids)

    R = _start_R()
    # Select the BioMart dataset to use.
    #mart = "ensembl"
    mart = "ENSEMBL_MART_ENSEMBL"  # Changed 151120.
    host = "www.ensembl.org"
    R_fn("useMart", mart, in_mart, host=host, RETVAL="in.dataset")
    R_fn("useMart", mart, out_mart, host=host, RETVAL="out.dataset")

    # Link two data sets and retrieve information from the linked datasets.
    jmath.R_equals_vector(gene_ids, 'gene.ids')


    # ERROR:
    #   Error in getLDS(attributes = "ensembl_gene_id", filters =
    #   "ensembl_gene_id", : The query to the BioMart webservice
    #   returned an invalid result: the number of columns in the result
    #   table does not equal the number of attributes in the
    #   query. Please report this to the mailing list.
    # Can mean that the gene IDs are bad.  E.g. version numbers still
    # on entrez IDs.
    R_fn(
        "getLDS", attributes=in_attribute, filters=in_attribute,
        values=R_var("gene.ids"), mart=R_var("in.dataset"),
        attributesL=out_attribute, martL=R_var("out.dataset"),
        RETVAL="homolog")
    homolog = R["homolog"]
    # homolog is DataFrame with two parallel rows:
    # <in_ids>
    # <out_ids>
    assert len(homolog) == 2, \
           "BioMart returned no results mapping from %s:%s to %s:%s." % (
        in_mart, in_attribute, out_mart, out_attribute)

    in_ids = [str(x) for x in homolog[0]]
    out_ids = [str(x) for x in homolog[1]]

    # Sometimes BioMart will generate "NA" if something is missing.
    if no_na:
        for i in range(len(out_ids)):
            if out_ids[i].upper() == "NA":
                out_ids[i] = ""
    
    in2out = {}
    for x, y in zip(in_ids, out_ids):
        if not y.strip():
            continue
        val = in2out.get(x, [])
        val.append(y)
        in2out[x] = sorted(val)

    return in2out


def _convert_gene_ids_local(in_platform, out_platform):
    # Return a dictionary of gene_id -> list of converted_ids, or None
    # if these platforms cannot be converted.
    import os
    from genomicode import config
    from genomicode import filelib

    filelib.assert_exists_nz(config.convert_platform)
    x = "%s___%s.txt" % (in_platform, out_platform)
    filename = os.path.join(config.convert_platform, x)
    if not os.path.exists(filename):
        return None

    in2out = {}
    for cols in filelib.read_cols(filename):
        # <in_id>  <out_id1> ... <out_idn>
        assert len(cols) >= 2
        in_id = cols[0]
        out_ids = cols[1:]
        in2out[in_id] = out_ids
    return in2out
        

FOUND_ID2ENTREZ = {}

def _convert_entrez_symbol_to_entrez(gene_ids, out_platform, db_tax_id):
    # Return a dictionary of gene_id -> list of converted_ids.
    # Need db_tax_id to limit search space.
    global FOUND_ID2ENTREZ

    FOUND = FOUND_ID2ENTREZ

    # Not sure if other organisms are implemented.
    assert db_tax_id in ["9606", "10090"]
    #db_tax_id = "9606"
    #db_tax_id = "10090"

    in2out = {}

    for in_id in gene_ids:
        gene_id = gene_symbol = None

        # Try to find the gene ID from the symbol.
        # First, look to see if it's already been found.
        if in_id in FOUND:
            gene_id, gene_symbol = FOUND[in_id]
        if not gene_id:
            x = _find_entrez_gene(in_id, db_tax_id)
            if x:
                gene_id, gene_symbol = x

        if out_platform == "Entrez_ID_human" and gene_id:
            in2out[in_id] = [str(gene_id)]
        elif out_platform == "Entrez_Symbol_human" and gene_symbol:
            in2out[in_id] = [gene_symbol]
    return in2out


FIND_GENE_ERROR = None
def _find_entrez_gene(gene_symbol, tax_id):
    # Return tuple of (gene_id, gene_symbol) or None.
    global FIND_GENE_ERROR
    from genomicode import genefinder

    FIND_GENE_ERROR = None

    # Try to find the gene.
    try:
        x = genefinder.find_gene(gene_symbol, tax_id=tax_id)
    except AssertionError, x:
        FIND_GENE_ERROR = str(x)
        return None
    if x:
        gene_id = x[0]
        gene_symbol = x[1]
        return gene_id, gene_symbol

    # Could not find gene.  Look for a discontinued gene.
    try:
        x = genefinder.find_discontinued(gene_symbol, tax_id=tax_id)
        if x:
            x = genefinder.find_gene(x, tax_id=tax_id)
    except AssertionError, x:
        FIND_GENE_ERROR = str(x)
        return None
    if x:
        gene_id = x[0]
        gene_symbol = x[1]
        return gene_id, gene_symbol

    return None
    

def _remove_dups(ids):
    # Preserve order.
    ids_c = []
    for x in ids:
        if x not in ids_c:
            ids_c.append(x)
    return ids_c
