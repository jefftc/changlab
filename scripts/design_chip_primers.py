#!/usr/bin/env python

import sys

def main():
    from optparse import OptionParser, OptionGroup
    from genomicode import genomefns
    from genomicode import primer3fns
    from genomicode import parsefns

    # Defines the targeted region.
    # gene    E2F1 ENSA,0 (<gene>,<tss>).
    # chrom   1 2 X
    # offset  If gene, then relative to TSS.  Otherwise, relative to chrom.
    # length  How long.
    #
    # negative offset is "n100" for -100, otherwize Python will
    # interpret it as a flag.
    usage = "usage: %prog [options] <gene_or_chrom> <offset> <length>"
    parser = OptionParser(usage=usage, version="%prog 01")

    parser.add_option(
        "--product_size", dest="product_size", default=[], action="append",
        help="Add a product size to search, e.g. 75-100.")
    parser.add_option(
        "-n", dest="num_primers", type="int", default=None,
        help="Number of primer pairs to pick.")
    parser.add_option(
        "-v", dest="verbose", action="store_true", default=False,
        help="Make output verbose.")

    options, args = parser.parse_args()
    if len(args) != 3:
        print usage
        sys.exit(-1)
    gene_or_chrom, offset, length = args
    neg = False
    if offset.startswith("n"):
        neg = True
        offset = offset[1:]
    offset = int(offset)
    if neg:
        offset = -offset
    length = int(length)
    assert length > 0

    assert options.num_primers is None or options.num_primers > 0

    if not options.product_size:
        # Default Product Size Range for web is:
        # 150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000
        options.product_size = [
            (50, 100),
            (101, 150),
            (151, 200), 
            (201, 300),
            ]
    for mn, mx in options.product_size:
        assert mn > 0
        assert mn < mx

    # gene, tss_offset, gene_length
    # chrom, pos, chrom_length

    # Figure out if a gene or chrom was requested.
    gene_symbol = chrom = None
    gene_or_chrom = gene_or_chrom.upper()
    try:
        gene_or_chrom = int(gene_or_chrom)
    except ValueError, x:
        pass
    if type(gene_or_chrom) is type(0) or gene_or_chrom in "XYM":
        chrom = gene_or_chrom
        chrom_pos = offset
        chrom_length = length
    else:
        gene_symbol = gene_or_chrom
        gene_offset = offset   # negative is upstream, positive is downstream
        gene_length = length
    assert chrom or gene_symbol
    assert not (chrom and gene_symbol)

    # If a gene was specified, convert to chromosomal coordinates.
    if gene_symbol:
        # See if a transcript num was specified.
        transcript_num = None
        if "," in gene_symbol:
            gene_symbol, x = gene_symbol.split(",", 1)
            transcript_num = int(x)
            assert transcript_num >= 0
        
        x = genomefns.get_gene_coords(gene_symbol)
        genes = genomefns.filter_unique_tss(x)
        assert genes

        # If there is only 1 gene, then use that one.
        if len(genes) == 1 and transcript_num is None:
            transcript_num = 0

        if len(genes) > 1 and transcript_num is None:
            print "Multiple transcripts for %s." % gene_symbol
            print "Please specify a specific transcript."
            x = "Index", "Chrom", "Strand", "TSS"
            print "\t".join(map(str, x))
            for i in range(len(genes)):
                tss = genomefns.transcript2tss(
                    genes[i].txn_start, genes[i].txn_length,genes[i].strand)
                tss = parsefns.pretty_int(tss)
                x = i, genes[i].chrom, genes[i].strand, tss
                print "\t".join(map(str, x))
            return
        
        assert transcript_num is not None
        assert transcript_num >= 0 and transcript_num < len(genes), \
               "Invalid transcript: %s" % transcript_num
        gene = genes[transcript_num]

        x = genomefns.transcript2promoter(
            gene.txn_start, gene.txn_length, gene.strand, gene_offset,
            gene_length)
        base, length, strand = x
        chrom = gene.chrom
        chrom_base = base
        chrom_length = length

        if options.verbose:
            tss = genomefns.transcript2tss(
                gene.txn_start, gene.txn_length, gene.strand)
            tss_str = parsefns.pretty_int(tss)
            print "%s: TSS at chr%s:%s:%s." % (
                gene_symbol, gene.chrom, tss_str, gene.strand)

    base_str = parsefns.pretty_int(chrom_base)
    if options.verbose:
        print "Target at chr%s:%s:%s." % (chrom, base_str, "+")
        print

    # Retrieve the target sequence, plus some flanking region for the
    # primers.
    max_size = 300
    if options.product_size:
        max_size = max([x[1] for x in options.product_size])
    BUFFER = 20  # provide a few bases for the primer to bind
    FLANK_SIZE = max_size + BUFFER
    start = chrom_base - FLANK_SIZE
    length = chrom_length + FLANK_SIZE * 2
    seq = genomefns.get_sequence(chrom, start, length)
    seq = seq.upper()
    target = FLANK_SIZE+1, chrom_length  # primer3 target is 1-based

    if options.verbose:
        # Print out the sequence.
        x1 = parsefns.pretty_int(start)
        x2 = parsefns.pretty_int(start+length)
        title = "chr%s %s-%s" % (chrom, x1, x2)
        s = seq[:FLANK_SIZE].lower() + \
            seq[FLANK_SIZE:FLANK_SIZE+chrom_length] + \
            seq[FLANK_SIZE+chrom_length:].lower()
        genomefns.write_fasta(title, s)
        print

    primers = primer3fns.primer3(
        seq, product_size=options.product_size, target=target,
        num_return=options.num_primers)
    if not primers:
        print "No primers found."
        return

    if options.verbose:
        x = ["Index",
             "L_Seq", "L_Pos", "L_Length", "L_Tm", "L_GC",
             "R_Seq", "R_Pos", "R_Length", "R_Tm", "R_GC",
             "Product Size"]
        print "\t".join(x)

    for i, x in enumerate(primers):
        d1, d2, size = x
        d2.seq_rc = genomefns.revcomp(d2.seq)

        F_i = seq.find(d1.seq) - FLANK_SIZE
        R_i = seq.find(d2.seq_rc) - FLANK_SIZE

        if options.verbose:
            x = (i+1,
                 d1.seq, F_i, d1.length, d1.tm, d1.gc_percent, 
                 d2.seq, R_i, d2.length, d2.tm, d2.gc_percent,
                 size)
            print "\t".join(map(str, x))
        else:
            print d1.seq
            print d2.seq
            print "Product Size=%d" % size
            print 

if __name__ == '__main__':
    main()
