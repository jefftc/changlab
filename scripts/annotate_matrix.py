#!/usr/bin/env python

# convert_geneset
# convert_matrix

import os
import sys


def convert_geneset(
    filename, in_delim, keep_dups, keep_emptys, no_na, 
    in_genesets, all_genesets, out_platforms, out_format):
    from genomicode import genesetlib
    from genomicode import arrayplatformlib as apl
    from genomicode import arrayannot

    # NOT IMPLEMENTED:
    # min_match_score

    assert in_genesets or all_genesets

    # Search for the gene sets to convert.
    name2geneset = {}   # name -> (name, description, gene_ids)
    print "HERE 2", filename
    for x in genesetlib.read_genesets(filename):
        name, description, gene_ids = x
        print "HERE 6", name, in_genesets
        if not all_genesets and name not in in_genesets:
            continue
        assert name not in name2geneset, "Duplicate gene set: %s" % name
        name2geneset[name] = x
        print "HERE 1", name

    # Make sure all the desired gene sets were found.
    for name in in_genesets:
        assert name in name2geneset, "Missing: %s" % name

    # Convert each of the gene sets.
    # list of name, description, in_platform, out_platform, in_ids, out_ids
    converted = []
    for x in sorted(name2geneset.values()):
        name, description, in_ids = x
        
        x = apl.score_annotations(in_ids, min_score=0.50)
        assert x, "Unknown platform: %s" % name
        best_score = x[0]
        in_platform = best_score.platform_name

        # Convert to each of the out platformss.
        for out_platform in out_platforms:
            out_delim = " /// "
            out_ids = arrayannot.convert_gene_ids(
                in_ids, in_platform, out_platform, in_delim, out_delim,
                keep_dups, keep_emptys, no_na)
            # Since gene sets are unordered and unaligned lists of genes,
            # separate any genes that are separated by out_delim.
            x = []
            for out_id in out_ids:
                x.extend(out_id.split(out_delim))
            x = [x.strip() for x in x]
            x = [x for x in x if x]
            x = sorted({}.fromkeys(x))
            out_ids = x

            x = name, description, in_platform, out_platform, in_ids, out_ids
            converted.append(x)

    # Convert to gene sets.
    genesets = []
    for x in converted:
        name, description, in_platform, out_platform, in_ids, out_ids = x
        out_name = name
        if len(out_platforms) > 1:
            out_name = "%s_%s" % (name, out_platform)
        x = genesetlib.GeneSet(out_name, description, out_ids)
        genesets.append(x)
        
    # Save the gene sets.
    if out_format == "gmt":
        genesetlib.write_gmt(sys.stdout, genesets)
    elif out_format == "gmx":
        genesetlib.write_gmx(sys.stdout, genesets)
    else:
        raise AssertionError
    
    
def convert_matrix(
    filename, header, header_and_platform, in_delim, out_delim,
    keep_dups, keep_emptys, no_na, out_platforms, min_match_score, debug):
    import arrayio
    from genomicode import Matrix
    from genomicode import arrayplatformlib as apl
    from genomicode import arrayannot

    MIN_SCORE = 0.80
    REMOVE_VERSION = True

    assert not (header and header_and_platform)

    DATA = arrayio.read(filename)

    if header:
        x = DATA.row_names(header)
        gene_ids = apl.normalize_ids(
            x, delimiter=in_delim, remove_version_number=REMOVE_VERSION)
        x = apl.score_annotations(gene_ids, min_score=0.5)
        assert x, "I could not identify the platform for %s." % header
        best_score = x[0]
        in_platform, score = best_score.platform_name, best_score.max_score
    elif header_and_platform:
        x = header_and_platform.split(",", 1)
        assert len(x) == 2
        header, in_platform = x
        score = 1.0
        x = DATA.row_names(header)
        gene_ids = apl.normalize_ids(
            x, delimiter=in_delim, remove_version_number=REMOVE_VERSION)
        assert apl.find_platform_by_name(in_platform), \
               "Unknown platform: %s" % in_platform
    else:
        # Take the platform with the highest match score.
        scores = apl.score_matrix(
            DATA, annot_delim=in_delim, min_score=None,
            remove_version=REMOVE_VERSION)
        best_score = 0
        if scores:
            best_score = scores[0].max_score
        if best_score < MIN_SCORE and debug and scores:
            header = (
                "Header", "Platform", "Score",
                "Matrix Only", "Plat Only", "Shared",
                "Matrix Only", "Plat Only", "Shared")
            print "\t".join(header)
            for s in scores:
                x1 = sorted(s.mine_only)[:3]
                x2 = sorted(s.platform_only)[:3]
                x3 = sorted(s.shared)[:3]
                x1 = ", ".join(x1)
                x2 = ", ".join(x2)
                x3 = ", ".join(x3)
                x = (
                    s.header, s.platform_name, s.max_score, 
                    len(s.mine_only), len(s.platform_only), len(s.shared),
                    x1, x2, x3)
                assert len(x) == len(header)
                print "\t".join(map(str, x))
        assert best_score >= MIN_SCORE, "No platforms found"
        best_score = scores[0]
        header = best_score.header
        in_platform = best_score.platform_name
        score = best_score = best_score.max_score
    err = "I could not find any platforms.  The best was %s (%g)." % (
        in_platform, score)
    assert score >= min_match_score, err
    gene_ids = DATA.row_names(header)

    # Convert each of the platforms.
    output_ids_list = []
    for out_platform in out_platforms:
        x = arrayannot.convert_gene_ids(
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
    
    from genomicode import arrayplatformlib as apl

    parser = argparse.ArgumentParser(
        description='Annotate a matrix or a geneset.')
    parser.add_argument("infile")
    
    all_platforms = [platform.name for platform in apl.PLATFORMS]
    x = ", ".join(all_platforms)
    parser.add_argument(
        '--platform', default=[], action='append',
        help="Which platform to add to the matrix.  Options: %s" % x)
    parser.add_argument(
        "--debug", action="store_true",
        help="If having trouble identifying the platform for the matrix, "
        "this will print out extra information.")

    group = parser.add_argument_group(title="Platform Match")
    group.add_argument(
        '--min_match_score', default=0.80, type=float,
        help="When trying to identify the rows of a matrix or geneset, "
        "require at least this portion of the IDs to be recognized.")
    group.add_argument(
        '--in_delim', 
        help="If a row contains multiple annotations (or gene names), they "
        "are separated by this delimiter, e.g. E2F1,E2F3")
    group.add_argument(
        '--out_delim', default=" /// ",
        help="Delimiter to use for the converted gene IDs.")
    group.add_argument(
        '--keep_dups', default=False, action="store_true",
        help="Keep duplicate IDs (e.g. to preserve alignment).")
    group.add_argument(
        '--keep_emptys', default=False, action="store_true",
        help="Keep empty IDs (e.g. to preserve alignment).")
    group.add_argument(
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
        'Required for geneset files.  (MULTI)')
    group.add_argument(
        "--all_genesets", action="store_true", help="Convert all genesets.")
    #group.add_argument("--out_geneset_name")
    group.add_argument(
        "--out_geneset_format", default="gmt", choices=["gmt", "gmx"],
        help="For output geneset file (default: gmt).")
    
    args = parser.parse_args()
    assert os.path.exists(args.infile), "File not found: %s" % args.infile
    assert args.platform, 'Please give at least one platform to add.'
    for x in args.platform:
        assert x, "Unknown platform: %s" % x
        #assert arrayplatformlib.get_bm_organism(x), "Unknown platform: %s" % x
    assert len(args.geneset) <= 1, "Not implemented."

    annotate_matrix = False
    annotate_geneset = False
    if args.header or args.header_and_platform:
        annotate_matrix = True
    if args.geneset or args.all_genesets:
        annotate_geneset = True
    #assert annotate_matrix or annotate_geneset
    assert not (annotate_matrix and annotate_geneset)
    assert not (args.header and args.header_and_platform)
    assert not (args.geneset and args.all_genesets)
    
    assert type(args.min_match_score) is type(0.0)
    assert args.min_match_score > 0.2, "min_match_score too low"
    assert args.min_match_score <= 1.0, "min_match_score too high"

    if annotate_geneset:
        print "HERE 3"
        convert_geneset(
            args.infile, args.in_delim,  args.keep_dups, args.keep_emptys,
            args.no_na,  args.geneset, args.all_genesets, args.platform,
            args.out_geneset_format)
    else:
        print "HERE 4"
        convert_matrix(
            args.infile, args.header, args.header_and_platform,
            args.in_delim, args.out_delim,
            args.keep_dups, args.keep_emptys, args.no_na, args.platform,
            args.min_match_score, args.debug)
            
            
if __name__=='__main__':
    main()
