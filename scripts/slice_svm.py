#!/usr/bin/env python

# Functions:
# filter_min_callers


def filter_min_callers(MATRIX, args):
    if args is None:
        return MATRIX
    from genomicode import AnnotationMatrix
    
    num_callers = args
    assert num_callers >= 1 and num_callers < 20

    I_nc = [i for (i, x) in enumerate(MATRIX.headers)
         if x.startswith("Num Callers")]
    headers_nc = [MATRIX.headers_h[i] for i in I_nc]

    I_remove = []
    for i in range(MATRIX.num_annots()):
        has_sample = False
        for h in headers_nc:
            x = MATRIX.header2annots[h][i]
            if not x.strip():
                continue
            nc = int(x)
            if nc >= num_callers:
                has_sample = True
                break
        if not has_sample:
            I_remove.append(i)

    x = {}.fromkeys(I_remove)
    I_keep = [i for i in range(MATRIX.num_annots()) if i not in x]
    filtered_matrix = AnnotationMatrix.rowslice(MATRIX, I_keep)
    return filtered_matrix


def filter_min_total_reads(MATRIX, args):
    if args is None:
        return MATRIX
    raise NotImplementedError


def exonic_only(MATRIX, args):
    if not args:
        return MATRIX
    from genomicode import AnnotationMatrix

    header = "Annovar______Func.refGene"
    assert header in MATRIX.headers_h

    I_keep = []
    func = MATRIX.header2annots[header]
    for i in range(len(func)):
        # exonic
        # ncRNA_exonic;splicing
        # exonic;splicing
        x = func[i]
        x = x.split(";")
        if "exonic" not in x:
            continue
        I_keep.append(i)
    MATRIX = AnnotationMatrix.rowslice(MATRIX, I_keep)
    return MATRIX


def annotate_linked_variants(MATRIX, args):
    if not args:
        return MATRIX
    from genomicode import filelib
    from genomicode import AnnotationMatrix

    link_file = args
    filelib.assert_exists_nz(link_file)
    coord2perc = {}
    for d in filelib.read_row(link_file, header=1):
        chrom = d.Chrom
        pos = int(d.Pos)
        perc = float(d.Perc_Linked)
        coord2perc[(chrom, pos)] = perc

    chrom = MATRIX.header2annots["______Chrom"]
    pos = MATRIX.header2annots["______Pos"]
    pos = [int(x) for x in pos]

    link_score = [""] * len(chrom)
    for i in range(len(chrom)):
        link_score[i] = coord2perc.get((chrom[i], pos[i]), "")
    
    # Add after:
    # Chrom, Pos, Ref, Alt
    header = "Linkage______Score"
    assert header not in MATRIX.headers
    headers = MATRIX.headers[:4] + [header] + MATRIX.headers[4:]
    all_annots = []
    for h in headers:
        if h != header:
            x = MATRIX[h]
        else:
            x = link_score
        all_annots.append(x)
    return AnnotationMatrix.create_from_annotations(
        headers, all_annots, MATRIX.headerlines)


def main():
    import sys
    import argparse
    from genomicode import SimpleVariantMatrix

    parser = argparse.ArgumentParser(
        description="Perform operations on a SimpleVariantMatrix file.")
    parser.add_argument("filename", nargs=1, help="Annotation file.")

    group = parser.add_argument_group(title="Matrix operations")

    group.add_argument(
        "--filter_min_callers", type=int, 
        help="Discard variants if no samples have at least this many "
        "callers.")
    group.add_argument(
        "--filter_min_total_reads", type=int, 
        help="Discard calls if no samples have at least this many "
        "callers.")
    group.add_argument(
        "--exonic_only", action="store_true",
        help="Keep variants only if they are exonic.")
    group.add_argument(
        "--annotate_linked_variants", 
        help="Add a column that shows the linkage score for each variant.  "
        "Format: <linkage file>.")
    
    args = parser.parse_args()
    assert len(args.filename) == 1
    FILENAME = args.filename[0]

    # Read the matrix.
    MATRIX = SimpleVariantMatrix.read_as_am(FILENAME)

    MATRIX = filter_min_callers(MATRIX, args.filter_min_callers)
    MATRIX = filter_min_total_reads(MATRIX, args.filter_min_total_reads)
    MATRIX = exonic_only(MATRIX, args.exonic_only)
    MATRIX = annotate_linked_variants(MATRIX, args.annotate_linked_variants)

    # Write the matrix back out.
    SimpleVariantMatrix.write_from_am(sys.stdout, MATRIX)

if __name__ == '__main__':
    main()
