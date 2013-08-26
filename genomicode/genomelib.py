"""

Functions:
get_sequence
get_transcript
len_chrom
revcomp

load_genes
get_gene_coords
get_promoters
filter_unique_tss

find_nearby_tss
find_overlapping_genes

transcript2tss
transcript2promoter
tss2promoter
genbase2tssbase
tssbase2genbase

calc_tss_base_dist       
calc_tss_seq_dist        
calc_seq_seq_dist        WAS calc_seq_distance WAS calc_5p_offset
calc_seq_seq_sep         WAS calc_seq_separation WAS calc_raw_offset

read_fasta
read_fasta_many
write_fasta
read_ra
write_ra

"""
# _get_default_ra_chrom
# _get_default_knowngene
# _load_gene_file
# _overlaps
# _safe_split_int
# _assert_chrom

import os, sys

def _get_default_ra_chrom():
    import config
    return config.RA_CHROM_HG19

def _get_default_knowngene():
    import config
    return config.knowngene_hg19

def get_sequence(chrom, start, length, ra_path=None):
    # Low complexity sequences will be in lower case.
    _assert_chrom(chrom)
    ra_path = ra_path or _get_default_ra_chrom()
    filename = os.path.join(ra_path, "chr%s.ra" % chrom)
    x = read_ra(filename, start, length, "c")
    return "".join(x)

def get_transcript(
    chrom, strand, txn_start, txn_length, exon_starts, exon_lengths,
    mask_repeats=False, ra_path=None):
    assert len(exon_starts) == len(exon_lengths)

    # Pull out the sequence.
    seq = get_sequence(chrom, txn_start, txn_length, ra_path=ra_path)

    # Mark the exons as upper case, introns as lower case.
    seq = seq.lower()
    for start, length in zip(exon_starts, exon_lengths):
        assert start >= txn_start
        assert length > 0
        start -= txn_start
        end = start + length
        seq = seq[:start] + seq[start:end].upper() + seq[end:]

    if mask_repeats:
        # Mask the repeats to N.
        seq = list(seq)
        for i in range(len(seq)):
            if seq[i].islower():
                seq[i] = "N"
        seq = "".join(seq)

    if strand == "-":
        seq = revcomp(seq)
    return seq

def len_chrom(chrom, ra_path=None):
    import stat

    ra_path = ra_path or _get_default_ra_chrom()
    filename = os.path.join(ra_path, "%s.ra" % chrom)
    size = os.stat(filename)[stat.ST_SIZE]
    return size

def revcomp(seq):
    from Bio import Seq
    return Seq.Seq(seq).reverse_complement().tostring()

def get_gene_coords(gene_symbol, gene_file=None):
    # Return a list of objects with members:
    # kg_id
    # genbank_id
    # refseq_id
    # gene_id
    # gene_symbol
    # chrom          1, 2, 3, ..., X, Y, M, (string)
    # strand
    # txn_start
    # txn_length
    # tss
    # cds_start
    # cds_length
    # exon_starts
    # exon_lengths
    #
    # Will return multiple objects if this gene has different
    # transcripts.  Objects will be sorted so that first object is the
    # one whose transcription start is most upstream.

    # Find the transcript info for this gene.
    gene_file = gene_file or _get_default_knowngene()
    genes = load_genes(gene_file)
    gene_symbol = gene_symbol.upper()
    genes = [x for x in genes if x.gene_symbol.upper() == gene_symbol]
    if not genes:
        return []

    # Sort by transcription start that is most upstream.
    # Coordinates are 0-based, exclusive end.
    tss = [x.tss for x in genes]
    if genes[0].strand == "-":
        tss = [-x for x in tss]
    schwartz = zip(tss, genes)
    schwartz.sort()
    genes = [x[-1] for x in schwartz]
    
    return genes

def get_promoters(gene_symbol, prom_offset, prom_length, gene_file=None,
                  ra_path=None):
    # Return list of (chrom, tss, strand, prom_base, prom_length, prom_seq).
    # Will revcomp promoters on the - strand.
    genes = get_gene_coords(gene_symbol, gene_file=gene_file)
    transcripts = filter_unique_tss(genes)
    if not transcripts:
        return []
    promoters = []
    for g in transcripts:
        x = transcript2promoter(
            g.txn_start, g.txn_length, g.strand, prom_offset, prom_length)
        prom_base, prom_length, prom_strand = x

        seq = get_sequence(g.chrom, prom_base, prom_length, ra_path=ra_path)
        if g.strand == "-":
            seq = revcomp(seq)
        x = g.chrom, g.tss, g.strand, prom_base, prom_length, seq
        promoters.append(x)
    return promoters

def filter_unique_tss(genes):
    # Return a list of objects with members:
    # kg_id
    # genbank_id
    # refseq_id
    # gene_id
    # gene_symbol
    # chrom          1, 2, 3, ..., X, Y, M, (string)
    # strand
    # txn_start
    # txn_length
    # tss
    #
    # Objects will be sorted so that first object is the one whose
    # transcript is most upstream.
    import filelib

    # Filter out duplicate entries.
    i = 0
    while i < len(genes)-1:
        if genes[i].tss == genes[i+1].tss:
            del genes[i+1]
        else:
            i += 1

    # Clean up the objects.
    genes = [
        filelib.GenericObject(
            kg_id=x.kg_id, genbank_id=x.genbank_id, 
            refseq_id=x.refseq_id, gene_id=x.gene_id,
            gene_symbol=x.gene_symbol, chrom=x.chrom, strand=x.strand,
            txn_start=x.txn_start, txn_length=x.txn_length, tss=x.tss)
        for x in genes]

    # Sort by transcription start that is most upstream.
    # Coordinates are 0-based, exclusive end.
    tss = [x.tss for x in genes]
    if genes and genes[0].strand == "-":
        tss = [-x for x in tss]
    schwartz = zip(tss, genes)
    schwartz.sort()
    genes = [x[-1] for x in schwartz]
    
    return genes

def find_nearby_tss(chrom, base, max_bases, gene_file=None):
    # Return a list of genes based on distance to the transcription
    # start site.  max_bases is the maximum number of bases to the
    # tss.  base should be 0-based.
    _assert_chrom(chrom)

    gene_file = gene_file or _get_default_knowngene()
    genes = load_genes(gene_file, chrom)
    assert genes, "No genes on chromosome %s." % chrom

    # Sort by increasing distance.
    dist = [calc_tss_base_dist(x.tss, x.strand, base) for x in genes]
    adist = [abs(x) for x in dist]
    schwartz = zip(adist, genes)
    schwartz.sort()
    adist = [x[0] for x in schwartz]
    genes = [x[-1] for x in schwartz]
    genes = [g for (g, dist) in zip(genes, adist) if dist < max_bases]
    return genes

def find_overlapping_genes(chrom, base, length, gene_file=None):
    # Format of chrom: "1", "8", "X", "Y", etc...
    # Return a list of objects (see load_genes for description).
    _assert_chrom(chrom)

    gene_file = gene_file or _get_default_knowngene()
    genes = load_genes(gene_file, chrom)
    #genes = [x for x in genes if x.chrom == chrom]
    assert genes, "No genes on chromosome %s." % chrom

    # Maybe can speed up with bisect?
    genes = [
        x for x in genes if _overlaps(base, length, x.txn_start, x.txn_length)]
    return genes

def transcript2tss(start, length, strand):
    # Given a stretch of transcribed sequence, determine the position
    # at which transcription starts.
    #assert left <= right
    assert strand in "+-"
    assert length > 0
    if strand == "+":
        return start
    return start + length - 1

def transcript2promoter(
    txn_start, txn_length, txn_strand, prom_offset, prom_length):
    # Return base, length, strand.  prom_offset of -1000 means 1000
    # bases uptream of the transcription start.
    assert txn_start >= 0
    assert txn_length > 0
    assert txn_strand in "+-"

    tss_base = transcript2tss(txn_start, txn_length, txn_strand)
    x = tss2promoter(tss_base, txn_strand, prom_offset, prom_length)
    return x

def tss2promoter(tss_base, txn_strand, prom_offset, prom_length):
    # Calculate the left-most base of the promoter.
    assert prom_length > 0
    prom_base = tss_base + prom_offset
    if txn_strand == "-":
        # Should be +1?  If so, will be off-by-1 compared to genome
        # browser.  If not, then math doesn't work out.
        prom_base = tss_base - prom_offset
        prom_base = prom_base - prom_length + 1
    return prom_base, prom_length, txn_strand
    
def genbase2tssbase(genome_base, tss, txn_strand):
    # Convert a genome coordinate to a coordinate relative to the TSS.
    # genpos is a position on the genome.  tss is the position of the
    # TSS on the genome.  strand indicates the strand of the
    # transcript.
    assert txn_strand in "+-"
    if txn_strand == "+":
        return genome_base - tss
    return tss - genome_base

def tssbase2genbase(tss_base, tss, txn_strand):
    assert txn_strand in "+-"
    if txn_strand == "+":
        return tss + tss_base
    return tss - tss_base

def calc_tss_base_dist(tss, txn_strand, base):
    # Calculate the distance between a base and a TSS.
    assert txn_strand in "+-"
    if txn_strand == "+":
        return base - tss
    return tss - base

def calc_tss_seq_dist(tss, txn_strand, base, length):
    # Calculate the distance of a piece of sequence (e.g. TFBS) from
    # the TSS.  The distance is the minimum distance between 
    d1 = calc_tss_base_dist(tss, txn_strand, base)
    d2 = calc_tss_base_dist(tss, txn_strand, base+length-1)
    return min(d1, d2)

def calc_seq_seq_dist(base1, length1, strand1, base2, length2):
    # Calculate the distance of the second stretch of sequence
    # relative to the 5p end of the first stretch.
    #
    # ---01234----43210------
    # ---012-345---1234------
    # ---210----------567----
    #
    #  SEQ1  SEQ2  DISTANCE
    # 01234   012     0
    # 01234   210     0
    # 01234   345     4
    # 43210  1234     0
    # 43210   567    -2
    # 
    # The strand of the second sequence is irrelevant.
    # Semantics defined in journal:070513.
    assert strand1 in ["-", "+"]
    offset = base2 - base1
    if strand1 == "-":
        base1 = base1 + length1 - 1
        base2 = base2 + length2 - 1
        offset = base1 - base2
    return offset

def calc_seq_seq_sep(base1, length1, base2, length2):
    # Calculate the number of bases between two stretches of sequence,
    # e.g. transcription factor binding sites.  A negative distance
    # indicates the number of bases that overlap (up to the length of
    # the shortest tfbs).
    # 
    # ---------12345------
    # --123-----45-6789--0
    #
    #  SEQ1  SEQ2  SEPARATION
    # 12345   123      4
    # 12345    45     -2
    # 12345  6789     -1
    # 12345     0      5
    if base1 > base2:
        base1, length1, base2, length2 = base2, length2, base1, length1
    offset = base2 - (base1+length1)
    offset = max(offset, -length1, -length2)
    return offset

def read_fasta(fh):
    # Return a tuple of (title, sequence)
    for x in read_fasta_many(fh):
        yield x

def read_fasta_many(fh):
    # Yield tuples of (title, sequence)
    import filelib

    handle = filelib.openfh(fh)
    title, sequence = "", []
    for line in handle:
        if line.startswith(">"):
            if title or sequence:
                yield title, "".join(sequence)
            title = line[1:].strip()
            sequence = []
        else:
            sequence.append(line.strip())
    if title or sequence:
        yield title, "".join(sequence)

def write_fasta(title, sequence, width=60, handle=None):
    handle = handle or sys.stdout
    w = handle.write
    w(">%s\n" % title)
    i = 0
    while i < len(sequence):
        s = sequence[i:i+width]
        i += width
        w("%s\n" % s)

def read_ra(filename, start, length, typecode):
    # length is number of items to read (not number of bytes)
    import array
    import stat

    ar = array.array(typecode)
    bytes_per_pos = ar.itemsize
    pos = start * bytes_per_pos

    size = os.stat(filename)[stat.ST_SIZE]
    assert size >= pos+length, "%s (%d): %d+%d" % (filename, size, pos, length)
    handle = open(filename, 'rb')
    handle.seek(pos)
    
    ar.fromfile(handle, length)
    return ar.tolist()

def write_ra(filename, start, data, typecode):
    import array

    ar = array.array(typecode)
    bytes_per_pos = ar.itemsize
    pos = start * bytes_per_pos

    mode = "r+b"
    if not os.path.exists(filename):
        mode = "w+b"
    handle = open(filename, mode)
    handle.seek(pos)
    
    ar.fromlist(data)
    ar.tofile(handle)

GENES_KEY = GENES_ALL = GENES_CHROM = None
def load_genes(gene_file=None, chrom=None):
    # Return a list of objects with members:
    # kg_id          UCSC KnownGene ID
    # genbank_id
    # refseq_id
    # gene_id
    # gene_symbol
    # chrom          1, 2, 3, ..., X, Y, M, (string)
    # strand
    # txn_start
    # txn_length
    # tss
    # cds_start
    # cds_length
    # exon_starts
    # exon_lengths
    # 
    # Genes are sorted by chrom, txn_start.
    global GENES_KEY, GENES_ALL, GENES_CHROM

    gene_file = gene_file or _get_default_knowngene()
    
    if GENES_KEY != gene_file:
        data = _load_genes_h(gene_file)
        # Hash by chromosome.
        GENES_CHROM = {}  # chrom -> list of genes
        for d in data:
            if d.chrom not in GENES_CHROM:
                GENES_CHROM[d.chrom] = []
            _assert_chrom(d.chrom)
            GENES_CHROM[d.chrom].append(d)
        GENES_KEY, GENES_ALL = gene_file, data
    if chrom is not None:
        _assert_chrom(chrom)
        return GENES_CHROM.get(chrom, [])
    return GENES_ALL
    
def _load_genes_h(gene_file):
    # This function is really slow for some reason.
    import filelib

    gene_file = gene_file or _get_default_knowngene()
    gene_data = _load_gene_file(gene_file)

    genes = []
    for d in gene_data:
        # Ignore weird chromosomes, e.g.:
        # chr6_apd_hap1
        # chr17_ctg5_hap1
        # chr1_gl00191_random
        if d.chrom.find("_") >= 0:
            continue
        chrom = d.chrom
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        _assert_chrom(chrom)

        assert d.txn_end >= d.txn_start
        assert d.cds_end >= d.cds_start
        txn_length = d.txn_end - d.txn_start
        tss = transcript2tss(d.txn_start, txn_length, d.strand)
        cds_length = d.cds_end - d.cds_start
        exon_lengths = []
        assert len(d.exon_starts) == len(d.exon_ends)
        for start, end in zip(d.exon_starts, d.exon_ends):
            assert end >= start
            length = end - start
            exon_lengths.append(length)
        assert len(exon_lengths) == len(d.exon_starts)
        
        params = dict(
            refseq_id=d.refseq_id,
            gene_id=d.gene_id, gene_symbol=d.gene_symbol,
            chrom=chrom, strand=d.strand,
            txn_start=d.txn_start, txn_length=txn_length, tss=tss,
            cds_start=d.cds_start, cds_length=cds_length,
            exon_starts=d.exon_starts, exon_lengths=exon_lengths)
        if hasattr(d, "kg_id"):
            params["kg_id"] = d.kg_id
        if hasattr(d, "genbank_id"):
            params["genbank_id"] = d.genbank_id
        if hasattr(d, "ensembl_id"):
            params["ensembl_id"] = d.ensembl_id
        if hasattr(d, "hgnc_symbol"):
            params["hgnc_symbol"] = d.hgnc_symbol
        obj = filelib.GenericObject(**params)
        genes.append(obj)

    # Sort by (chrom, txn_start, kg_id) just to make sure this
    # function always returns by a unique sorting order.
    schwartz = []
    for g in genes:
        try:
            chrom = int(g.chrom)
        except ValueError:
            chrom = g.chrom
        if hasattr(g, "kg_id"):
            x = g.kg_id
        else:
            x = g.refseq_id
        x = chrom, g.txn_start, x, g
        schwartz.append(x)
    schwartz.sort()
    genes = [x[-1] for x in schwartz]
        
    return genes

def _load_gene_file_h(gene_file):
    import filelib

    # kg_id is unique.  Others may be duplicated or missing.
    data = []
    for d in filelib.read_row(gene_file, header=1):
        d.txn_start, d.txn_end = int(d.txn_start), int(d.txn_end)
        d.cds_start, d.cds_end = int(d.cds_start), int(d.cds_end)
        d.exon_starts = _safe_split_int(d.exon_starts)
        d.exon_ends = _safe_split_int(d.exon_ends)
        
        assert d.strand in ["+", "-"]
        assert len(d.exon_starts) == len(d.exon_ends)
        data.append(d)
    return data

TSS_KEY = TSS_VALUE = None
def _load_gene_file(gene_file):
    global TSS_KEY, TSS_VALUE

    # Coordinates are 0-based, exclusive end.
    if TSS_KEY != gene_file:
        data = _load_gene_file_h(gene_file)
        TSS_KEY, TSS_VALUE = gene_file, data
    return TSS_VALUE

def _safe_split_int(s):
    return [int(x) for x in s.split(",") if x]

def _overlaps(s1, l1, s2, l2):
    if s1 > s2:
        s1, l1, s2, l2 = s2, l2, s1, l2
    return s1+l1 > s2

GOOD_CHROM = {}  # chrom -> 1
def _assert_chrom(chrom):
    # There's a limited number of acceptable chromosomes, so cache the
    # good ones so I don't have to check it again.
    global GOOD_CHROM

    if chrom in GOOD_CHROM:
        return
    chrom = chrom.replace("_random", "")
    chrom = chrom.replace("_cox_hap1", "")
    chrom = chrom.replace("_qbl_hap2", "")
    chrom = chrom.replace("_h2_hap1", "")
    try:
        x = int(chrom)
    except ValueError, x:
        assert chrom in "XYM" or chrom == "MT", "Invalid (1) chrom: %s" % chrom
    else:
        assert x >= 1 and x <= 25, "Invalid (2) chrom: %s" % chrom
        # danRer has up to 25.
    GOOD_CHROM[chrom] = 1
