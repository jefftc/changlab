#!/usr/bin/env python

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


def convert_platform(gene_ids, in_platform, out_platform):
    # gene_ids is a list.  Return a dictionary where each key is a
    # gene_id from the in_platform, and the value is a list of the
    # gene_ids from the out_platform.
    from genomicode import jmath
    from genomicode import arrayplatformlib

    R_fn, R_var = jmath.R_fn, jmath.R_var

    # Clean up the gene_ids.
    gene_ids = sorted({}.fromkeys(gene_ids))
    gene_ids = [x for x in gene_ids if x]

    # An attribute is the biomart name for the platform.
    in_attribute = arrayplatformlib.get_bm_attribute(in_platform)
    out_attribute = arrayplatformlib.get_bm_attribute(out_platform)
    assert in_attribute
    assert out_attribute
    in_mart = arrayplatformlib.get_bm_organism(in_platform)
    out_mart = arrayplatformlib.get_bm_organism(out_platform)

    R = start_R()
    jmath.R_equals_vector(gene_ids, 'gene_ids')

    # Select the BioMart dataset to use.
    R_fn("useMart", "ensembl", in_mart, RETVAL="in_dataset")
    R_fn("useMart", "ensembl", out_mart, RETVAL="out_dataset")

    # Link two data sets and retrieve information from the linked datasets.
    R_fn(
        "getLDS", attributes=in_attribute, filters=in_attribute,
        values=R_var("gene_ids"), mart=R_var("in_dataset"),
        attributesL=out_attribute, martL=R_var("out_dataset"),
        RETVAL="homolog")
    
    homolog = R['homolog']
    in_ids = [str(i) for i in homolog[0]]
    out_ids = [str(i) for i in homolog[1]]
    in2out = {}
    for x, y in zip(in_ids, out_ids):
        if not y.strip():
            continue
        val = in2out.get(x, [])
        val.append(y)
        in2out[x] = sorted(val)
    return in2out
    

def main():
    import argparse
    
    import arrayio
    from genomicode import Matrix
    from genomicode import arrayplatformlib
    from genomicode import genesetlib

    parser = argparse.ArgumentParser(
        description='Annotate a matrix or a geneset.')
    parser.add_argument("infile")
    all_platforms = [platform.name for platform in arrayplatformlib.PLATFORMS]
    parser.add_argument(
        '--platform', default=[], action='append',
        help='specify the platform to add:'+ str(all_platforms))

    
    parser.add_argument(
        '--geneset', default=[], action='append',
        help='Annotate a gene set.  infile must be a gene set.')
    parser.add_argument("--out_geneset_name")
    parser.add_argument(
        "--geneset_format", default="gmt", choices=["gmt", "gmx"])
    
    args = parser.parse_args()
    assert os.path.exists(args.infile), "File not found: %s" % args.infile
    assert args.platform, 'please give at least one platform for convert'
    for x in args.platform:
        assert arrayplatformlib.get_bm_organism(x), \
               'we cannot convert to the platform %s' % x
    assert len(args.geneset) <= 1, "Not implemented."

    is_matrix = not args.geneset

    if is_matrix:
        DATA = arrayio.read(args.infile)
        platform_list = arrayplatformlib.identify_all_platforms_of_matrix(DATA)
        assert platform_list, 'we cannot guess the platform for the input file'
        # Why is the first one chosen?
        header = platform_list[0][0]
        gene_ids = DATA.row_names(header)
        in_platform = platform_list[0][1]
    else:
        assert args.geneset
        assert len(args.geneset) == 1
        DATA = genesetlib.read_genes(args.infile, args.geneset[0])
        gene_ids = DATA
        x = arrayplatformlib.score_platform_of_annotations(gene_ids)
        in_platform = x[0]
    assert in_platform, "I could not figure out the platform of the infile."


    # Convert each of the platforms.
    platform2geneid2outids = {}
    for out_platform in args.platform:
        x = convert_platform(gene_ids, in_platform, out_platform)
        platform2geneid2outids[out_platform] = x


    if is_matrix:
        # Make a matrix with the new IDs.
        X = DATA._X
        row_names = DATA._row_names.copy()
        row_order = DATA._row_order[:]
        col_names = DATA._col_names.copy()
        col_order = DATA._col_order[:]
        synonyms = DATA._synonyms.copy()

        for out_platform in args.platform:
            geneid2outids = platform2geneid2outids[out_platform]
            x = [geneid2outids.get(x, []) for x in gene_ids]
            x = [" /// ".join(x) for x in x]
            converted_ids = x

            header = out_platform
            i = 1
            while header in row_order:
                header = "%s_%d" % (out_platform, i)
                i += 1
            row_order.append(header)
            row_names[header] = converted_ids

        # Write the outfile.
        x = Matrix.InMemoryMatrix(
            X, row_names=row_names, col_names=col_names,
            row_order=row_order, col_order=col_order, synonyms=synonyms)
        arrayio.tab_delimited_format.write(x, sys.stdout)
    else:
        # Write out the gene set.
        genesets = []
        for out_platform in args.platform:
            geneid2outids = platform2geneid2outids[out_platform]

            name = out_platform
            if args.out_geneset_name:
                # BUG: If multiple gene sets, will generate duplicate names.
                name = args.out_geneset_name
            
            description = "na"
            genes = []
            for gene_id in gene_ids:
                genes.extend(geneid2outids.get(gene_id, []))
            genes = sorted({}.fromkeys(genes))
            x = genesetlib.GeneSet(name, description, genes)
            genesets.append(x)

        if args.geneset_format == "gmt":
            genesetlib.write_gmt(sys.stdout, genesets)
        elif args.geneset_format == "gmx":
            genesetlib.write_gmx(sys.stdout, genesets)
        else:
            raise AssertionError
            
    
            
if __name__=='__main__':
    main()
