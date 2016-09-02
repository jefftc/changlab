#!/usr/bin/env python


def parse_indexes(MATRIX, indexes, indexes_include_headers):
    from genomicode import parselib
    
    max_index = MATRIX.ncol()
    num_headers = len(MATRIX._row_names)
    assert max_index, "empty matrix"

    I = []
    for s, e in parselib.parse_ranges(indexes):
        if indexes_include_headers:
            s, e = s-num_headers, e-num_headers
        assert s >= 1, "Index out of range: %s" % s
        assert e <= max_index, "Index out of range: %s" % e
        s, e = s - 1, min(e, max_index)
        I.extend(range(s, e))
    return I


def main():
    import os
    import sys
    import argparse
    
    import arrayio
    from genomicode import dwdnorm

    usage = "dwdnorm.py [options] expression_file"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("expression_file", help="Matrix to normalize.")
    parser.add_argument(
        "--indexes1", 
        help="Which columns in batch 1, E.g. 1-5,8 (1-based, "
        "inclusive).  Columns leftover will be in the remaining batches.")
    parser.add_argument(
        "--indexes_include_headers", "--iih", action="store_true",
        help="If not given (default), then column 1 is the first column "
        "with data.  If given, then column 1 is the very first column in "
        "the file, including the headers.")
    
    args = parser.parse_args()
    assert os.path.exists(args.expression_file), \
        "File not found: %s" % args.expression_file
    assert args.indexes1

    MATRIX = arrayio.read(args.expression_file)

    # Figure out the batch for each sample.
    assert MATRIX.ncol(), "empty matrix"
    # Should be -1 or 1.
    batches = [-1] * MATRIX.ncol()
    I = parse_indexes(MATRIX, args.indexes1, args.indexes_include_headers)
    for i in I:
        batches[i] = 1
    # Make sure there are two batches.
    count = {}
    for b in batches:
        count[b] = count.get(b, 0) + 1
    assert len(count) == 2

    MATRIX_norm = dwdnorm.normalize(MATRIX, batches)

    arrayio.tab_delimited_format.write(MATRIX_norm, sys.stdout)

 
if __name__=='__main__':
    main()
