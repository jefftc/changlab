#!/usr/bin/env python

import sys


class ExonSequence:
    def __init__(self, full_seq):
        # seq should be a strings.  Exons in upper case, all other
        # bases in lower case.

        # seq-length vector of exon ids.
        full_exonid = [0] * len(full_seq)  # 0 is intron, 1+ is exon
        id = 0
        for i in range(len(full_seq)):
            if not full_seq[i].isupper():
                continue
            if not i or not full_seq[i-1].isupper():
                id += 1
            full_exonid[i] = id
        assert len(full_seq) == len(full_exonid)

        sub_indexes = []   # exon-length array of indexes into full array
        for i, id in enumerate(full_exonid):
            if id > 0:
                sub_indexes.append(i)

        sub_seq = "".join([full_seq[i] for i in sub_indexes])
        sub_exonid = [full_exonid[i] for i in sub_indexes]

        self.full_seq = full_seq
        self.full_exonid = full_exonid
        self.sub_indexes = sub_indexes
        self.sub_seq = sub_seq
        self.sub_exonid = sub_exonid
        
    def num_exons(self):
        return max(self.full_exonid)
    def get_seq(self):
        return self.full_seq
    def get_exon(self, exon_id):
        # exon_id starts from 1.
        assert exon_id >= 1
        x = [self.full_seq[i] for i in self.sub_indexes
             if self.full_exonid[i] == exon_id]
        x = "".join(x)
        return x
    def get_index(self, exon_id):
        # 0-based index where exon_id is in the sequence.
        assert exon_id >= 1
        x = [i for i, id in enumerate(self.full_exonid) if id == exon_id]
        assert x, "missing exon: %s" % exon_id
        return x[0]
    def get_all_exons(self):
        return self.sub_seq

def format_exon_ids(exon_ids):
    count = {}
    for id in exon_ids:
        count[id] = count.get(id, 0) + 1

    if len(count) == 1:
        id_str = str(count.keys()[0])
    else:
        x = []
        for id in sorted(count):
            x.append("%s(x%d)" % (id, count[id]))
        id_str = "-".join(x)
    return id_str

def main():
    from optparse import OptionParser, OptionGroup
    from genomicode import genomefns
    from genomicode import primer3fns
    from genomicode import parsefns

    # gene    E2F1 ENSA,0 (<gene>,<tss>).
    usage = "usage: %prog [options] <gene>"
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
    if len(args) != 1:
        print usage
        sys.exit(-1)
    gene, = args

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

    gene_symbol, transcript_num = gene.upper(), None
    # See if a transcript num was specified.
    if "," in gene_symbol:
        gene_symbol, x = gene_symbol.split(",", 1)
        transcript_num = int(x)
        assert transcript_num >= 0
        
    genes = genomefns.get_gene_coords(gene_symbol)
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

    seq = genomefns.get_transcript(
        gene.chrom, gene.strand, gene.txn_start, gene.txn_length,
        gene.exon_starts, gene.exon_lengths)
    #print gene.txn_start, gene.txn_length
    #genomefns.write_fasta("HELLO", seq)
    
    # Search for primers on just the exons.
    exon_seq = ExonSequence(seq)
    primers = primer3fns.primer3(
        exon_seq.sub_seq, product_size=options.product_size,
        num_return=options.num_primers)
    if not primers:
        print "No primers found."
        return

    # Revcomp the right primer.
    for d1, d2, size in primers:
        d2.seq_rc = genomefns.revcomp(d2.seq)

    if options.verbose:
        x = ["Index",
             "L_Seq", "L_Pos", "L_Length", "L_Tm", "L_GC", "L_Exon",
             "R_Seq", "R_Pos", "R_Length", "R_Tm", "R_GC", "R_Exon", 
             "Genome Size", "Product Size"]
        print "\t".join(x)

    for zzz, x in enumerate(primers):
        d1, d2, size = x

        # Calculate the size of the genomic product.
        gen_left = exon_seq.full_seq.find(d1.seq.upper())
        gen_right = exon_seq.full_seq.find(d2.seq_rc.upper())
        gen_size = 0
        if gen_left >= 0 and gen_right >= 0:
            gen_size = gen_right - gen_left + len(d2.seq_rc)

        # Calculate the size of the PCR product.
        pcr_left = exon_seq.sub_seq.find(d1.seq.upper())
        pcr_right = exon_seq.sub_seq.find(d2.seq_rc.upper())
        pcr_size = pcr_right - pcr_left + len(d2.seq)
        assert pcr_size == size
        
        # Figure out the exons of the primers.
        L_exon_ID = [exon_seq.sub_exonid[pcr_left+i]
                     for i in range(len(d1.seq))]
        R_exon_ID = [exon_seq.sub_exonid[pcr_right+i]
                     for i in range(len(d2.seq))]
        L_exon_ID_str = format_exon_ids(L_exon_ID)
        R_exon_ID_str = format_exon_ids(R_exon_ID)

        # Ignore primers on the same exon.
        if L_exon_ID == R_exon_ID:
            continue
        
        if options.verbose:
            x = (zzz+1,
                 d1.seq, pcr_left, d1.length, d1.tm, d1.gc_percent,
                 L_exon_ID_str, 
                 d2.seq, pcr_right, d2.length, d2.tm, d2.gc_percent,
                 R_exon_ID_str, 
                 gen_size, size)
            print "\t".join(map(str, x))
        else:
            print d1.seq
            print d2.seq
            print "Product Size=%d" % size
            print 

if __name__ == '__main__':
    main()
