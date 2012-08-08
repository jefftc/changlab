"""

Functions:
score_tfbs_genome    Score the TFBS in a genomic location (from filesystem).

load_matrices
list_matrices        Return a list of (matrix ID, gene symbol).
matid2matrix         Return the matrix for a matrix ID.
gene2matrices        Return all matrices for a gene symbol.
is_matrix_id

find_matrix_file

OBSOLETE
get_tfbs_genome_db   Load the TFBS in a genomic location (from database).
get_tfbs_tss_db      Load the TFBS around a TSS.

"""
# _query_db




import os, sys

def score_tfbs_genome(chrom, start, length, matrices=None, nlp=None,
                      num_procs=1):
    # Return list of matrix, chrom, strand, position, NLP.
    # NLP is given in log_e.
    import genomelib
    import patser

    seq = genomelib.get_sequence(chrom, start, length)
    if not matrices:
        matrices = [x[0] for x in list_matrices()]
    x = [find_matrix_file(x) for x in matrices]
    matrix_files = [x for x in x if x]  # Filter out missing ones.

    # list of sequence_num, matrix_num, (0-based) position, strand, score, nlp
    patser_data = patser.score_tfbs(seq, matrix_files, num_jobs=num_procs)

    # list of matrix, chrom, strand, pos, nlp
    data = []
    for x in patser_data:
        sequence_num, matrix_num, position, strand, score, tf_nlp = x
        assert sequence_num == 0
        matrix = matrices[matrix_num]
        pos = position + start

        if nlp is not None and tf_nlp < nlp:
            continue
        x = matrix, chrom, strand, pos, tf_nlp
        data.append(x)

    # Sort by chrom, position, strand, matrix
    x = [(x[1], x[3], x[2], x[0], x) for x in data]
    x.sort()
    data = [x[-1] for x in x]

    return data

def _load_matrices_h():
    import config
    import filelib

    # MATRIX_ID (all upper case) -> object with members:
    #   matid
    #   accession
    #   gene_id
    #   gene_symbol
    #   organism     (for JASPAR)
    #   length

    # Load the lengths.
    matid2length = {}
    for d in filelib.read_row(
        config.motiflib_MATID2LENGTH, "matid:s length:d"):
        matid2length[d.matid] = d.length

    # Load the matrix information.
    matrices = []
    for d in filelib.read_row(config.motiflib_JASPAR_INFO, header=1):
        assert d.xID in matid2length
        length = matid2length.get(d.xID, 0)
        x = filelib.GenericObject(
            matid=d.xID, accession="", gene_id=d.LocusLink,
            gene_symbol=d.Gene_Symbol, organism=d.Organism, length=length)
        matrices.append(x)
    for d in filelib.read_row(config.motiflib_TRANSFAC_INFO, header=1):
        #assert d.Accession in matid2length, "Missing: %s" % d.Accession
        # Some are missing, e.g. M00316.
        length = matid2length.get(d.Accession, 0)
        x = filelib.GenericObject(
            matid=d.Accession, accession=d.xID, gene_id=d.LocusLink,
            gene_symbol=d.Gene_Symbol, organism="", length=length)
        matrices.append(x)
    return matrices

MATRIX_CACHE = None
def load_matrices():
    # Return a list of objects with members:
    # matid        matrix ID
    # accession    for TRANSFAC only
    # gene_id
    # gene_symbol
    # length
    global MATRIX_CACHE
    if not MATRIX_CACHE:
        MATRIX_CACHE = _load_matrices_h()
    return MATRIX_CACHE
    
def list_matrices():
    # Return a list of matrix_id, gene_symbol
    matrix_db = load_matrices()
    x = [(d.matid, d.gene_symbol) for d in matrix_db]
    return x
    
def matid2matrix(matid):
    # Return a matrix object (see load_matrices).
    matrix_db = load_matrices()
    x = [x for x in matrix_db if x.matid.upper() == matid.upper()]
    assert x, "I couldn't find matrix: %s." % matid
    assert len(x) == 1, "Multiple matches for matid %s" % matid
    return x[0]

def gene2matrices(gene_symbol):
    ugene_symbol = gene_symbol.upper()
    matrices = load_matrices()
    x = [x for x in matrices if x.gene_symbol.upper() == ugene_symbol]
    return x

def is_matrix_id(matid):
    matrix_db = load_matrices()
    x = [x for x in matrix_db if x.matid.upper() == matid.upper()]
    return len(x) > 0

def find_matrix_file(matrix_file_or_id):
    # Take a matrix ID and return the name of the file or None.
    import config

    if os.path.exists(matrix_file_or_id):
        return matrix_file_or_id
    # If it's not a file, it must be an ID.
    matrix_id = matrix_file_or_id
    
    opj = os.path.join
    files = [
        opj(config.motiflib_JASPAR_DB, "%s.pfm" % matrix_id),
        opj(config.motiflib_TRANSFAC_DB, "%s.pfm" % matrix_id),
        #opj(config.MOTIFSEARCH, "matrices/%s.matrix" % matrix_id),
        ]
    for file in files:
        if os.path.exists(file):
            return file
    return None
    #raise AssertionError, "Cannot find matrix file for %s" % matrix_id

def get_tfbs_genome_db(chrom, start, end, matrices=None, nlp=None):
    # Return list of matrix, chrom, strand, position, NLP.
    # NLP is given in log_e.

    import bisect
    assert start >= 0 and end >= 0
    assert start <= end

    # Get a list of the genpos records that I should consider.
    # Should bind variable for chrom, but when I do that, I get 0 rows
    # back?  Also reported on the internet, but no fix.
    x = """SELECT gps_id, gps_position FROM genpos
        WHERE RTRIM(gps_chrom) = :chrom AND
          (gps_position = (
            SELECT MAX(gps_position) FROM genpos WHERE
            RTRIM(gps_chrom) = :chrom AND gps_position <= :start1) OR
          gps_position = (
            SELECT MAX(gps_position) FROM genpos WHERE
            RTRIM(gps_chrom) = :chrom AND gps_position < :end1))
        ORDER BY gps_position
        """
    results = _query_db(x, chrom=chrom, start1=start, end1=end)
    gps_ids = [x[0] for x in results]
    assert gps_ids, "Could not find genome position: %s [%d-%d]" % (
        chrom, start, end)

    x = ",".join(map(str, gps_ids))
    gpsid_query = "AND gps_id IN (%s)" % x
    matrix_query = ""
    if matrices is not None:
        x = ",".join(["'%s'" % x for x in matrices])
        matrix_query = "AND RTRIM(mat_name) IN (%s)" % x
    nlp_query = ""
    if nlp is not None:
        nlp_query = "AND tbs_nlp >= %g" % nlp
    x = """SELECT mat_name, gps_chrom, gps_strand, gps_position, tbs_offset,
            tbs_nlp
        FROM tfbs, genpos, matrix
        WHERE tbs_gpsid=gps_id AND tbs_matid=mat_id
        %(gpsid_query)s
        %(matrix_query)s
        %(nlp_query)s
        AND gps_position+tbs_offset >= :start1
        AND gps_position+tbs_offset < :end1
        ORDER BY gps_chrom, gps_position+tbs_offset, gps_strand, mat_name
        """ % locals()
    results = _query_db(x, start1=start, end1=end)
    BATCH_SIZE = 10000
    data, num = [], 0
    for x in results:
        matrix, chrom, strand, position, offset, nlp = x
        chrom = chrom.strip()
        pos = position + offset
        if num >= len(data):
            data += [None] * BATCH_SIZE
        data[num] = matrix, chrom, strand, pos, nlp
        num += 1
    return data[:num]

def get_tfbs_tss_db(chrom, txn_start, txn_end, txn_strand, bp_upstream, length,
                    matrices=None, nlp=None):
    # Return list of matrix, chrom, strand, gen_pos, tss, tss_pos, NLP.
    # gen_pos is the position of the matrix on the chromosome.
    # tss_pos is the position of the matrix relative to the TSS.
    # NLP is given in log_e.
    import genomelib
    
    x = genomelib.transcript2genome(
        txn_start, txn_end, txn_strand, bp_upstream, length)
    gen_start, gen_end = x
    tss = genomelib.determine_tss(txn_start, txn_end, txn_strand)

    data = get_tfbs_genome_db(
        chrom, gen_start, gen_end, matrices=matrices, nlp=nlp)
    results = []
    for x in data:
        matrix, chrom, mat_strand, mat_pos, nlp = x
        tss_pos = genomelib.genpos2tsspos(mat_pos, tss, txn_strand)
        x = matrix, chrom, mat_strand, mat_pos, tss, tss_pos, nlp
        results.append(x)
    return results

def _query_db(*args, **keywds):
    import cx_Oracle

    # Should put this in a configuration file.
    connection = cx_Oracle.connect("jefftc", "imaoraclem", "student1")
    cursor = connection.cursor()
    cursor.arraysize = 50
    cursor.execute(*args, **keywds)
    return cursor.fetchall()
