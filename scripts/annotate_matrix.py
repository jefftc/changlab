#!/usr/bin/env python

# start_R
# convert_gene_ids
# _convert_gene_ids_biomart
# _convert_gene_ids_local
# _convert_entrez_symbol_to_entrez

# convert_geneset
# convert_matrix
#
# _remove_upds
# _remove_refseq_version
# _clean_id

import os
import sys


GLOBAL_R = None
def start_R():
    global GLOBAL_R
    from genomicode import jmath

    if GLOBAL_R is None:
        R = jmath.start_R()
        R('library(biomaRt)')
        GLOBAL_R = R
    return GLOBAL_R


def _remove_dups(ids):
    ids_c = []
    for x in ids:
        if x not in ids_c:
            ids_c.append(x)
    return ids_c

def _remove_refseq_version(refseq_id):
    i = refseq_id.find(".")
    if i < 0:
        return refseq_id
    return refseq_id[:i]


def _clean_id(gene_id, delimiter, in_platform):
    # Return a list of this ID cleaned up.  gene_id might be split
    # into multiple IDs by the delimiter.
    x = [gene_id]
    if delimiter:
        x = gene_id.split(delimiter)
    # Clean up whitespace.
    x = [x.strip() for x in x]
    # No empty IDs.
    x = [x for x in x if x]
    # Hack: Remove version numbers from RefSeq IDs.
    if in_platform.lower().find("refseq") >= 0:
        x = [_remove_refseq_version(x) for x in x]
    # No duplicates.
    x = {}.fromkeys(x).keys()
    return x


def convert_gene_ids(
    gene_ids, in_platform, out_platform, in_delim, out_delim,
    keep_dups, keep_emptys, no_na):

    # Make a cleaned up version of the gene_ids to convert.
    x = []
    for gene_id in gene_ids:
        x.extend(_clean_id(gene_id, in_delim, in_platform))
    # No duplicates.
    x = {}.fromkeys(x).keys()
    gene_ids_c = x

    in2out = None
    if in_platform == "Entrez_Symbol_human" and \
       out_platform in ("Entrez_Symbol_human", "Entrez_ID_human"):
        in2out = _convert_entrez_symbol_to_entrez(gene_ids_c, out_platform)
    if in2out is None:
        in2out = _convert_gene_ids_local(in_platform, out_platform)
    if in2out is None:
        in2out = _convert_gene_ids_biomart(
            gene_ids_c, in_platform, out_platform, no_na)
    assert in2out, "I could not convert %s to %s" % (in_platform, out_platform)

    # Make a parallel list of the output IDs.
    output_ids = []
    for gene_id in gene_ids:
        in_ids = _clean_id(gene_id, in_delim, in_platform)
        out_ids = []
        for x in in_ids:
            x = in2out.get(x, [""])
            out_ids.extend(x)
        if not keep_emptys:
            out_ids = [x for x in out_ids if x]
        if not keep_dups:
            out_ids = _remove_dups(out_ids)
        x = out_delim.join(out_ids)
        output_ids.append(x)
    return output_ids

def clean_genes_for_biomart(gene_ids):
    gene_ids = [i.replace("'",'') for i in gene_ids] #modify "2'-PDE"
    return gene_ids

def _convert_gene_ids_biomart(gene_ids, in_platform, out_platform, no_na):
    # Return a dictionary of gene_id -> list of converted_ids, or None
    # if these platforms cannot be converted.
    from genomicode import jmath
    from genomicode import arrayplatformlib

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

    R = start_R()
    gene_ids = clean_genes_for_biomart(gene_ids)
    jmath.R_equals_vector(gene_ids, 'gene.ids')
    # Select the BioMart dataset to use.
    R_fn("useMart", "ensembl", in_mart, RETVAL="in.dataset")
    R_fn("useMart", "ensembl", out_mart, RETVAL="out.dataset")

    # Link two data sets and retrieve information from the linked datasets.
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
    from genomicode import config
    from genomicode import filelib
    
    assert os.path.exists(config.convert_platform)
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

def _convert_entrez_symbol_to_entrez(gene_ids, out_platform):
    # Return a dictionary of gene_id -> list of converted_ids.
    global FOUND_ID2ENTREZ

    FOUND = FOUND_ID2ENTREZ

    db_tax_id = "9606"
    #db_tax_id = "10090"

    in2out = {}

    for in_id in gene_ids:
        gene_id = gene_symbol = None

        # Try to find the gene ID from the symbol.
        # First, look to see if it's already been found.
        if in_id in FOUND:
            gene_id, gene_symbol = found[in_id]
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
    


## def convert_gene_ids(
##     gene_ids, in_platform, out_platform, in_delim, out_delim,
##     keep_dups, keep_emptys, no_na):
##     # gene_ids is a list of the gene IDs, one per row of the matrix,
##     # from the in_platform.  Each of the gene_ids may contain multiple
##     # IDs separated by in_delim.  in_delim is either None (no multiple
##     # IDs) or a character indicating the delimiter used to separate
##     # the IDs.
##     #
##     # Return a list parallel to gene_ids that is the IDs in the
##     # out_platform.
##     from genomicode import jmath
##     from genomicode import arrayplatformlib

##     R_fn, R_var = jmath.R_fn, jmath.R_var

##     # Make a cleaned up version of the gene_ids to convert.
##     x = []
##     for gene_id in gene_ids:
##         x.extend(_clean_id(gene_id, in_delim, in_platform))
##     # No duplicates.
##     x = {}.fromkeys(x).keys()
##     gene_ids_c = x

##     # An attribute is the biomart name for the platform.
##     in_attribute = arrayplatformlib.get_bm_attribute(in_platform)
##     out_attribute = arrayplatformlib.get_bm_attribute(out_platform)
##     assert in_attribute, "Bad platform: %s" % in_platform
##     assert out_attribute, "Bad platform: %s" % out_platform
##     in_mart = arrayplatformlib.get_bm_organism(in_platform)
##     out_mart = arrayplatformlib.get_bm_organism(out_platform)

##     R = start_R()
##     jmath.R_equals_vector(gene_ids_c, 'gene_ids')

##     # Select the BioMart dataset to use.
##     R_fn("useMart", "ensembl", in_mart, RETVAL="in_dataset")
##     R_fn("useMart", "ensembl", out_mart, RETVAL="out_dataset")

##     # Link two data sets and retrieve information from the linked datasets.
##     R_fn(
##         "getLDS", attributes=in_attribute, filters=in_attribute,
##         values=R_var("gene_ids"), mart=R_var("in_dataset"),
##         attributesL=out_attribute, martL=R_var("out_dataset"),
##         RETVAL="homolog")
    
##     homolog = R['homolog']
##     # homolog is DataFrame with two parallel rows:
##     # <in_ids>
##     # <out_ids>
##     assert len(homolog) == 2, \
##            "BioMart returned no results mapping from %s to %s." % (
##         in_mart, out_mart)

##     in_ids = [str(x) for x in homolog[0]]
##     out_ids = [str(x) for x in homolog[1]]

##     # Sometimes BioMart will generate "NA" if something is missing.
##     if no_na:
##         for i in range(len(out_ids)):
##             if out_ids[i].upper() == "NA":
##                 out_ids[i] = ""
    
##     in2out = {}
##     for x, y in zip(in_ids, out_ids):
##         if not y.strip():
##             continue
##         val = in2out.get(x, [])
##         val.append(y)
##         in2out[x] = sorted(val)

##     # Make a parallel list of the output IDs.
##     output_ids = []
##     for gene_id in gene_ids:
##         in_ids = _clean_id(gene_id, in_delim, in_platform)
##         out_ids = []
##         for x in in_ids:
##             x = in2out.get(x, [""])
##             out_ids.extend(x)
##         if not keep_emptys:
##             out_ids = [x for x in out_ids if x]
##         if not keep_dups:
##             out_ids = _remove_dups(out_ids)
##         x = out_delim.join(out_ids)
##         output_ids.append(x)
##     return output_ids


def convert_geneset(
    filename, in_delim, out_delim, keep_dups, keep_emptys, no_na, 
    in_genesets, out_platforms, out_geneset, out_format):
    from genomicode import genesetlib
    from genomicode import arrayplatformlib

    # NOT IMPLEMENTED:
    # min_match_score
    # delim -> in_delim, out_delim
    raise NotImplementedError, "Not tested after refactoring."

    assert len(in_genesets) == 1
    geneset = in_genesets[0]

    gene_ids = genesetlib.read_genes(filename, geneset)
    
    x = arrayplatformlib.score_platform_of_annotations(gene_ids)
    assert x, "I could not figure out the platform of the infile."
    in_platform, score = x

    # Convert each of the platforms.
    platform2geneid2outids = {}
    for out_platform in out_platforms:
        x = convert_gene_ids(
            gene_ids, in_platform, out_platform, in_delim, out_delim,
            keep_dups, keep_emptys, no_na)
        platform2geneid2outids[out_platform] = x

    # Write out the gene set.
    genesets = []
    for out_platform in out_platforms:
        geneid2outids = platform2geneid2outids[out_platform]

        name = out_platform
        if out_geneset:
            # BUG: If multiple gene sets, will generate duplicate names.
            name = out_geneset

        description = "na"
        genes = []
        for gene_id in gene_ids:
            genes.extend(geneid2outids.get(gene_id, []))
        genes = sorted({}.fromkeys(genes))
        x = genesetlib.GeneSet(name, description, genes)
        genesets.append(x)

    if out_format == "gmt":
        genesetlib.write_gmt(sys.stdout, genesets)
    elif out_format == "gmx":
        genesetlib.write_gmx(sys.stdout, genesets)
    else:
        raise AssertionError
    
    
def convert_matrix(
    filename, header, header_and_platform, in_delim, out_delim,
    keep_dups, keep_emptys, no_na, 
    out_platforms, min_match_score):
    import arrayio
    from genomicode import Matrix
    from genomicode import arrayplatformlib

    assert not (header and header_and_platform)

    DATA = arrayio.read(filename)

    if header:
        gene_ids = DATA.row_names(header)
        x = arrayplatformlib.score_platform_of_annotations(gene_ids)
        assert x, "I could not identify the platform for %s." % header
        in_platform, score = x
    elif header_and_platform:
        x = header_and_platform.split(",", 1)
        assert len(x) == 2
        header, in_platform = x
        score = 1.0
        gene_ids = DATA.row_names(header)
        assert arrayplatformlib.find_platform_by_name(in_platform), \
               "Unknown platform: %s" % in_platform
    else:
        # Take the platform with the highest match score.
        platforms = arrayplatformlib.score_all_platforms_of_matrix(
            DATA, annot_delim=in_delim)
        assert platforms, "No platforms found"
        schwartz = [(-x[-1], x) for x in platforms]
        schwartz.sort()
        platforms = [x[-1] for x in schwartz]
        header, in_platform, score = platforms[0]
    err = "I could not find any platforms.  The best was %s (%g)." % (
        in_platform, score)
    assert score >= min_match_score, err
    gene_ids = DATA.row_names(header)

    # Convert each of the platforms.
    output_ids_list = []
    for out_platform in out_platforms:
        x = convert_gene_ids(
            gene_ids, in_platform, out_platform, in_delim, out_delim,
            keep_dups, keep_emptys, no_na)
        output_ids_list.append(x)

    # Make a matrix with the new IDs.
    X = DATA._X
    row_names = DATA._row_names.copy()
    row_order = DATA._row_order[:]
    col_names = DATA._col_names.copy()
    col_order = DATA._col_order[:]
    synonyms = DATA._synonyms.copy()

    for (out_platform, output_ids) in zip(out_platforms, output_ids_list):
        header = out_platform
        i = 1
        while header in row_order:
            header = "%s_%d" % (out_platform, i)
            i += 1
        row_order.append(header)
        row_names[header] = output_ids

    # Write the outfile.
    x = Matrix.InMemoryMatrix(
        X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    arrayio.tab_delimited_format.write(x, sys.stdout)
    

def main():
    import argparse
    
    from genomicode import arrayplatformlib

    parser = argparse.ArgumentParser(
        description='Annotate a matrix or a geneset.')
    parser.add_argument("infile")
    
    all_platforms = [platform.name for platform in arrayplatformlib.PLATFORMS]
    x = ", ".join(all_platforms)
    parser.add_argument(
        '--platform', default=[], action='append',
        help="Which platform to add to the matrix.  Options: %s" % x)

    parser.add_argument(
        '--min_match_score', default=0.80, type=float,
        help="When trying to identify the rows of a matrix or geneset, "
        "require at least this portion of the IDs to be recognized.")
    parser.add_argument(
        '--in_delim', 
        help="If a row contains multiple annotations (or gene names), they "
        "are separated by this delimiter, e.g. E2F1,E2F3")
    parser.add_argument(
        '--out_delim', default=" /// ",
        help="Delimiter to use for the converted gene IDs.")
    parser.add_argument(
        '--keep_dups', default=False, action="store_true",
        help="Keep duplicate IDs (e.g. to preserve alignment).")
    parser.add_argument(
        '--keep_emptys', default=False, action="store_true",
        help="Keep empty IDs (e.g. to preserve alignment).")
    parser.add_argument(
        '--no_na', default=False, action="store_true",
        help="If any annotations are NA (e.g. missing), convert to empty "
        "string.")
                
    group = parser.add_argument_group(title="Matrix")
    group.add_argument(
        '--header', 
        help='Which header contains the gene IDs to convert from.  '
        'If not provided, will try to guess')
    group.add_argument(
        "--header_and_platform",
        help="Provide a header and the name of platform.  "
        "Format: <header>,<platform>.")
    
    group = parser.add_argument_group(title="Gene Set")
    group.add_argument(
        '--geneset', default=[], action='append',
        help='Which gene set to annotate (if infile is a gene set file).  '
        'Required for geneset files.')
    group.add_argument("--out_geneset_name")
    group.add_argument(
        "--out_geneset_format", default="gmt", choices=["gmt", "gmx"],
        help="For output geneset file.")
    
    args = parser.parse_args()
    assert os.path.exists(args.infile), "File not found: %s" % args.infile
    assert args.platform, 'Please give at least one platform to add.'
    for x in args.platform:
        assert arrayplatformlib.get_bm_organism(x), "Unknown platform: %s" % x
    assert len(args.geneset) <= 1, "Not implemented."

    assert not (args.header and args.geneset)
    assert not (args.header_and_platform and args.geneset)
    assert not (args.header and args.header_and_platform)
    
    assert type(args.min_match_score) is type(0.0)
    assert args.min_match_score > 0.2, "min_match_score too low"
    assert args.min_match_score <= 1.0, "min_match_score too high"

    if args.geneset:
        convert_geneset(
            args.infile, args.in_delim, args.out_delim,
            args.keep_dups, args.keep_emptys, args.no_na, 
            args.geneset, args.platform,
            args.out_geneset_name, args.out_geneset_format)
    else:
        convert_matrix(
            args.infile, args.header, args.header_and_platform,
            args.in_delim, args.out_delim,
            args.keep_dups, args.keep_emptys, args.no_na, args.platform,
            args.min_match_score)
            
            
if __name__=='__main__':
    main()
