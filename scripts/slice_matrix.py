#!/usr/bin/env python

# TODO: Find some way to specify order of operations.

import os, sys

def parse_names(s):
    return s.split(",")

def parse_indexes(s, max_index):
    from genomicode import parsefns
    
    I = []
    for s, e in parsefns.parse_ranges(s):
        assert s >= 1
        s, e = s-1, e
        e = min(e, max_index)
        I.extend(range(s, e))
    return I

def parse_indexes_or_names(MATRIX, is_row, s):
    # Examples:
    # 5
    # 1,5,10
    # 1-99,215-300
    # E2F1,E2F3,PCNA
    import re

    # Figure out whether this contains indexes or names.
    if not s:
        return []
    if re.match("[a-zA-Z]", s):
        names = parse_names(s)
        params = { "row" : names }
        if not is_row:
            params = { "col" : names }
        I_row, I_col = MATRIX._index(**params)
        I = I_row
        if not is_row:
            I = I_col
    else:
        max_index = MATRIX.nrow()
        if not is_row:
            max_index = MATRIX.ncol()
        I = parse_indexes(s, max_index)
    return I

def parse_geneset(MATRIX, is_row, filename, genesets):
    # Return a list of genes in the desired gene sets.  If genesets is
    # None, return the genes in all the genesets.
    from genomicode import genesetlib
    
    if filename is None:
        return []
    genes = genesetlib.read_genes(filename, *genesets)
    params = { "row" : genes }
    if not is_row:
        params = { "col" : genes }
    I_row, I_col = MATRIX._index(**params)
    I = I_row
    if not is_row:
        I = I_col
    return I

def find_row_indexes(MATRIX, rows, rowfile, row_geneset):
    if not rows and not rowfile and not row_geneset:
        return None
    I1 = parse_indexes_or_names(MATRIX, True, rows)
    I2 = parse_geneset(MATRIX, True, rowfile, row_geneset)
    I = I1 + [i for i in I2 if i not in I1]
    return I

def find_col_indexes(MATRIX, cols, colfile, col_geneset):
    if not cols and not colfile and not col_geneset:
        return None
    I1 = parse_indexes_or_names(MATRIX, False, cols)
    I2 = parse_geneset(MATRIX, False, colfile, col_geneset)
    I = I1 + [i for i in I2 if i not in I1]
    return I

def remove_row_name(MATRIX, name):
    assert name in MATRIX.row_names(), "I could not find name: %s" % name
    MATRIX_clean = MATRIX.matrix()
    Mc = MATRIX_clean
    assert name in Mc._row_names
    del Mc._row_names[name]
    if Mc._row_order and name in Mc._row_order:
        i = Mc._row_order.index(name)
        Mc._row_order.pop(i)
    # Bug: Does not remove from the synonyms list.
    return MATRIX_clean

def _align_geneset_to_matrix(matrix_annots, geneset_genes):
    # Return a list of the indexes required to align the geneset to
    # the rows of the matrix.  If it cannot be completely aligned,
    # then return None.
    from genomicode import jmath

    if len(matrix_annots) > len(geneset_genes):
        return None
    I = jmath.match(matrix_annots, geneset_genes)
    if None in I:
        return None
    return I

def _align_matrix_to_geneset(matrix_annots, geneset_genes):
    # Return a list of the indexes required to align the matrix to the
    # geneset.  If it cannot be completely aligned, then return None.
    from genomicode import jmath

    if len(geneset_genes) > len(matrix_annots):
        return None
    I = jmath.match(geneset_genes, matrix_annots)
    if None in I:
        return None
    return I

def match_rownames_to_geneset(MATRIX, all_genesets, geneset2genes):
    # Return tuple of (I_matrix, I_geneset) or None if no match can be
    # found.
    geneset_aligns = []  # list of (I_geneset, rowname, geneset)
    matrix_aligns = []   # list of (I_matrix, rowname, geneset)
    for name in MATRIX.row_names():
        annots = MATRIX.row_names(name)
        for gs in all_genesets:
            genes = geneset2genes[gs]
            I_geneset = _align_geneset_to_matrix(annots, genes)
            I_matrix = _align_matrix_to_geneset(annots, genes)
            if I_geneset:
                x = I_geneset, name, gs
                geneset_aligns.append(x)
            if I_matrix:
                x = I_matrix, name, gs
                matrix_aligns.append(x)
                
    # First, try to find a geneset that matches the matrix.  Favor
    # geneset_aligns over matrix_aligns to avoid changing the matrix.
    for x in geneset_aligns:
        I_geneset, rowname, geneset = x
        I_matrix = range(MATRIX.nrow())
        assert len(I_matrix) == len(I_geneset)
        return I_matrix, I_geneset

    # Otherwise, choose the match that generates the largest matrix.
    I_matrix = None
    for x in matrix_aligns:
        I, rowname, geneset = x
        if I_matrix is None or len(I) > len(I_matrix):
            I_matrix = I
    if I_matrix:
        I_geneset = range(len(I_matrix))
        return I_matrix, I_geneset

    return None

def add_genesets(MATRIX, geneset_file, genesets):
    # The gene sets in geneset_file should all be aligned in the same
    # order.
    from genomicode import jmath
    from genomicode import genesetlib
    if not genesets:
        return MATRIX
    geneset2genes = {}
    all_genesets = []  # preserve the order of the genesets
    for x in genesetlib.read_genesets(geneset_file):
        geneset, description, genes = x
        geneset2genes[geneset] = genes
        all_genesets.append(geneset)

    # Find an alignment between one of the matrix row_names and the
    # genesets.
    x = match_rownames_to_geneset(MATRIX, all_genesets, geneset2genes)
    assert x, "I could not match the matrix to a geneset."
    I_matrix, I_geneset = x
    MATRIX_new = MATRIX.matrix(I_matrix, None)

    # Add the new annotations to the MATRIX.
    for gs in genesets:
        assert gs in geneset2genes, "Missing geneset: %s" % gs
        assert gs not in MATRIX.row_names(), "duplicate name: %s" % gs
        
        x = geneset2genes[gs]
        genes = [x[i] for i in I_geneset]

        assert len(genes) == MATRIX_new.nrow()
        MATRIX_new._row_order.append(gs)
        MATRIX_new._row_names[gs] = genes
        
    return MATRIX_new
        
def main():
    import argparse
    import arrayio

    parser = argparse.ArgumentParser(
        description="Slice the rows or columns of a matrix.")
    parser.add_argument("filename", help="Matrix to slice.")
    parser.add_argument(
        "-o", default=None, metavar="OUTFILE", dest="outfile",
        help="Save to this file.  By default, writes output to STDOUT.")

    group = parser.add_argument_group(title="Slice options")
    group.add_argument(
        "--rows", default=None, 
        help="Which rows to include e.g. 1-50,75 (1-based, inclusive).  "
        "Can also be comma-separate list of the IDs.")
    group.add_argument(
        "--cols", default=None, 
        help="Which columns to include.")
    group.add_argument(
        "--rowfile", default=None,
        help="File containing the indexes or names of the rows "
        "to include.  This file should be in GMX or GMT format.  By "
        "default, I will select all rows that match any of the gene sets.  "
        "If you wish to select from a limited number of gene sets, you can "
        "specify them using the --row_geneset argument.")
    group.add_argument(
        "--row_geneset", action="append", default=[],
        help="Which gene set to use out of --rowfile.")
    group.add_argument(
        "--colfile", default=None,
        help="File containing the columns to include.")
    group.add_argument(
        "--col_geneset", action="append", default=[],
        help="Which gene set to use out of --colfile.")
        
    group = parser.add_argument_group(title="Other options")
    group.add_argument(
        "--add_geneset", action="append", default=[],
        help="Add this geneset (from rowfile) as an annotation.")
    group.add_argument(
        "--rm_annot", action="append", default=[],
        help="Remove this annotations from the matrix.")

    args = parser.parse_args()
    if not os.path.exists(args.filename):
        parser.error("File not found: %s" % args.filename)
    if args.row_geneset and not args.rowfile:
        parser.error("Row geneset provided without a rowfile.")
    if args.col_geneset and not args.colfile:
        parser.error("Col geneset provided without a colfile.")
    if args.add_geneset and not args.rowfile:
        parser.error("Geneset provided without a rowfile.")
        
    # Read the file.
    fmt_module = arrayio.guess_format(args.filename)
    assert fmt_module, "I could not figure out the format of the matrix file."
    MATRIX = fmt_module.read(args.filename)

    # Slice to a submatrix.
    I_row = find_row_indexes(MATRIX, args.rows, args.rowfile, args.row_geneset)
    I_col = find_col_indexes(MATRIX, args.cols, args.colfile, args.col_geneset)
    MATRIX = MATRIX.matrix(I_row, I_col)

    # Remove unwanted annotations.
    for name in args.rm_annot:
        MATRIX = remove_row_name(MATRIX, name)

    # Add any extra genesets.
    MATRIX = add_genesets(MATRIX, args.rowfile, args.add_geneset)

    # Write the outfile (in the same format).
    handle = sys.stdout
    if args.outfile:
        handle = open(args.outfile, 'w')
    fmt_module.write(MATRIX, handle)
    

if __name__ == '__main__':
    main()
