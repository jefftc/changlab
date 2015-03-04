#!/usr/bin/env python

import os
import sys


# Merge this with align_matrices?
class AnnotationMatrix:
    def __init__(self, name2annots, name_order):
        self.name2annots = name2annots.copy()
        self.name_order = name_order[:]


def read_infile(filename):
    from genomicode import genesetlib

    name_order = []
    name2annots = {}
    num_annots = None
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True):
        name, description, annots = x

        if num_annots is None:
            num_annots = len(annots)
        assert len(annots) == num_annots
        
        name_order.append(name)
        name2annots[name] = annots
    return AnnotationMatrix(name2annots, name_order)


def main():
    import argparse
    
    from genomicode import genomelib

    parser = argparse.ArgumentParser(
        description='Annotate the genes that are close to a genomic location.')
    parser.add_argument("chrom_header")
    parser.add_argument("start_header")
    parser.add_argument("end_header")
    parser.add_argument("infile")
    
    parser.add_argument(
        "--knowngene_file", help="e.g. tss/data/hg38.knownGene.141007.txt")
    parser.add_argument(
        "--distance", default=30000, type=int,
        help="Maximum distance away to look (in base pairs)")
    
    args = parser.parse_args()
    assert os.path.exists(args.infile), "File not found: %s" % args.infile
    assert os.path.exists(args.knowngene_file), "File not found: %s" % \
           args.knowngene_file
    assert args.distance >= 0 and args.distance < 1E7

    matrix = read_infile(args.infile)
    assert args.chrom_header in matrix.name2annots, "Missing header: %s" % \
           args.chrom_header
    assert args.start_header in matrix.name2annots, "Missing header: %s" % \
           args.start_header
    assert args.end_header in matrix.name2annots, "Missing header: %s" % \
           args.end_header

    header = matrix.name_order + [
        "Closest KG ID", "Closest Gene Symbol",
        "All KG ID (< %d bp)" % args.distance, "All Gene ID", \
        "All Gene Symbol", "All Distance"]
    print "\t".join(header)

    num_annots = len(matrix.name2annots[args.chrom_header])
    for zzz in range(num_annots):
        chrom = matrix.name2annots[args.chrom_header][zzz]
        start = matrix.name2annots[args.start_header][zzz]
        end = matrix.name2annots[args.end_header][zzz]
        start, end = int(start), int(end)
        if start > end:
            start, end = end, start

        #if d.chr.startswith("chrUn"):
        #    continue
        #elif d.chr.endswith("_random"):
        #    continue
        #elif d.chr.find("_hap") >= 0:
        #    continue

        # Calculate the middle of the region, and adjust the distance
        # accordingly.
        middle = start + (end-start)/2
        distance = args.distance + middle-start

        # Find the nearest gene.
        chrom = chrom.replace("chr", "")
        # genes sorted by increasing distance
        genes = genomelib.find_nearby_tss(
            chrom, middle, distance, gene_file=args.knowngene_file)

        # Calculate the distance for each gene.
        for gene in genes:
            dist1 = genomelib.calc_tss_base_dist(gene.tss, gene.strand, start)
            dist2 = genomelib.calc_tss_base_dist(gene.tss, gene.strand, end)
            dist = min(dist1, dist2)
            gene.distance = dist

        ## Get rid of blank genes.  Some genes may have known gene IDs
        ## (kg_id), but no gene symbols.  I only want the ones with
        ## gene symbols.
        #y = zip(gene_ids, gene_symbols, gene_distances)
        #y = [x for x in y if x[0]]
        #gene_ids = [x[0] for x in y]
        #gene_symbols = [x[1] for x in y]
        #gene_distances = [x[2] for x in y]

        # Get rid of duplicates.
        i = 0
        while i < len(genes)-1:
            j = i+1
            x1 = genes[i].kg_id, genes[i].gene_id, genes[i].gene_symbol
            x2 = genes[i].kg_id, genes[j].gene_id, genes[j].gene_symbol
            if x1 == x2:
                # Delete the one that is further away.
                del genes[j]
            else:
                i += 1

        closest_kg_id = ""
        closest_gene_symbol = ""
        if genes:
            closest_kg_id = genes[0].kg_id
            closest_gene_symbol = genes[0].gene_symbol

        # Format the information to print out.
        x1 = [x.kg_id for x in genes]
        x2 = [x.gene_id for x in genes]
        x3 = [x.gene_symbol for x in genes]
        x4 = [x.distance for x in genes]
        x1 = ",".join(x1)
        x2 = ",".join(x2)
        x3 = ",".join(x3)
        x4 = ",".join(map(str, x4))

        y1 = [matrix.name2annots[x][zzz] for x in matrix.name_order]
        y2 = closest_kg_id, closest_gene_symbol, x1, x2, x3, x4
        y = y1 + list(y2)
        assert len(y) == len(header)
        print "\t".join(map(str, y))
            
            
if __name__=='__main__':
    main()
