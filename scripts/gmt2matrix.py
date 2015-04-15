#!/usr/bin/env python

def main():
    import os
    import argparse
    import math

    from genomicode import genesetlib

    parser = argparse.ArgumentParser()
    parser.add_argument("gmt_file")

    args = parser.parse_args()
    assert os.path.exists(args.gmt_file)

    genesets = []
    for x in genesetlib.read_gmt(args.gmt_file):
        genesets.append(x)
    if not genesets:
        return

    # Figure out the number of columns needed for the file.
    num_cols = None
    for x in genesets:
        name, desc, genes = x
        nc = 1 + 1 + len(genes)
        if num_cols is None or nc > num_cols:
            num_cols = nc
        
    num_digits = int(math.floor(math.log(num_cols, 10))) + 1

    header = ["COL%0*d" % (num_digits, i+1) for i in range(num_cols)]
    print "\t".join(header)
    for x in genesets:
        name, desc, genes = x
        row = [name, desc] + genes
        row = row + [""]*(num_cols-len(row))
        print "\t".join(row)
        
    
    
    
    
    
        
    

if __name__ == '__main__':
    main()
