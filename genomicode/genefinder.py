"""
Functions:
find_gene               Look up the one best gene for a name.
find_gene_detailed      Return many possible matches, lots of info.

find_many_genes
find_many_genes_detailed

find_discontinued       Look for a gene name in the discontinued list.

"""
# _lookup_gene
# _find_entrez_gene_table
# _list_mapping_tables

# TODO:
# fuzzy match of gene names

def find_gene(name, tax_id=None, ignore_ambiguous_names=False):
    # Return gene_id, symbol, name, tax_id, organism, name_from_query
    # or None.  If name is a RefSeq ID, be sure there is no version
    # number.  E.g. "NP_000680.2" won't be found, but "NP_000680"
    # will.
    #
    # tax_id   9606  human   (should be given as an integer)
    #         10090  mouse
    #
    # If the name is ambiguous, will raise an exception.  If
    # ignore_ambiguous_names is a true value, will return None instead.
    
    # Optimization: an empty string will be ambiguous.
    name = str(name)
    orig_name = name
    x = find_many_genes([name], tax_id=tax_id)
    assert len(x) == 1
    num_matches, gene_id, symbol, name, tax_id, organism, name_from_query =x[0]
    if num_matches == 0:
        return None
    if num_matches > 1 and ignore_ambiguous_names:
        return None
    if num_matches > 1:
        msg = 'Multiple matches (%d) found for "%s".' % (
            num_matches, orig_name)
        x = find_many_genes_detailed([orig_name])
        all_llid = [x[0] for x in x]
        llid_str = ",".join(map(str, all_llid))
        if len(all_llid) <= 5:
            msg = 'Multiple matches (%s) found for "%s".' % (
                llid_str, orig_name)
        assert num_matches == 1, msg
    x = gene_id, symbol, name, tax_id, organism, name_from_query
    return x

## def find_all_genes(name):
##     # Return list of gene_id, symbol, name, tax_id, organism, name_from_query.
##     x = find_gene_detailed(name)
##     results = [(x[0], x[1], x[2], x[3], x[4], x[5]) for x in x]
##     results.sort()
##     # No duplicates.
##     i = 0
##     while i < len(results)-1:
##         if results[i] == results[i+1]:
##             del results[i]
##         else:
##             results += 1
##     return results
    
def find_gene_detailed(name, tax_id=None):
    # Return list of gene_id, symbol, name, tax_id, organism,
    # name_from_query, source_db, name_in_db.
    return find_many_genes_detailed([name], tax_id=tax_id)

def find_many_genes(genes, tax_id=None):
    # Return list of (num_matches, gene_id, symbol, name, tax_id,
    # organism, name_from_query) that is parallel to genes.  If there
    # are no matches or multiple matches, values (other than
    # num_matches and name_from_query) are set to None.  gene_id and
    # tax_id are integers.
    
    x = find_many_genes_detailed(genes, tax_id=tax_id)
    results = [(x[0], x[1], x[2], x[3], x[4], x[5]) for x in x]
    results.sort()
    # No duplicates.
    i = 0
    while i < len(results)-1:
        if results[i] == results[i+1]:
            del results[i]
        else:
            i += 1
    query2hits = {}
    for x in results:
        gene_id, symbol, name, tax_id, organism, name_from_query = x
        if name_from_query not in query2hits:
            query2hits[name_from_query] = []
        query2hits[name_from_query].append(x)
    clean = []
    for name_from_query in genes:
        hits = query2hits.get(name_from_query, [])
        num_matches = len(hits)
        gene_id = symbol = name = tax_id = organism = None
        if num_matches == 1:
            gene_id, symbol, name, tax_id, organism, x = hits[0]
            assert x == name_from_query
        x = num_matches, gene_id, symbol, name, tax_id, organism, \
            name_from_query
        clean.append(x)
    return clean

def find_many_genes_detailed(genes, tax_id=None):
    """Return list of (gene_id, symbol, name, tax_id, organism,
    name_from_query, source_db, name_in_db).  May not be in the same
    order as the query.

    gene_id and tax_id are integers.  Everything else is a string.

    """
    import itertools
    import config
    import dblib

    genes_str = dblib.format_list(genes)
    
    # gene_id, source_db, name_in_db, is_official
    results = []
    if 1:
        # Use the master mapping table.
        columns = "gene_id, symbol, official, source_db"
        q = "SELECT %s FROM %s WHERE symbol in (%s);" % (
            columns, "MASTERMAP", genes_str)
        x = dblib.query(
            q, config.gm_user, config.gm_passwd, config.gm_db, config.gm_host,
            config.gm_port, config.gm_socket)
        for x in x:
            gene_id, name_in_db, is_official, source_db = x
            gene_id = int(gene_id)
            x = gene_id, source_db, name_in_db, is_official
            results.append(x)
    else:
        # Query each mapping table separately.
        table_names = _list_mapping_tables()
        for table_name in table_names:
            columns = "gene_id, symbol, official"
            q = "SELECT %s FROM %s WHERE symbol in (%s);" % (
                columns, table_name, genes_str)
            x = dblib.query(
                q, config.gm_user, config.gm_passwd, config.gm_db,
                config.gm_host, config.gm_port, config.gm_socket)
            for x in x:
                gene_id, name_in_db, is_official = x
                gene_id = int(gene_id)
                x = gene_id, table_name, name_in_db, is_official
                results.append(x)

    # Collect all possible hits, organized by the names in the query.
    # name_from_query -> (gene_id, source_db, name_in_db, is_official)
    hits_by_query = {}
    for (name_from_query, x) in itertools.product(genes, results):
        gene_id, table_name, name_in_db, is_official = x
        if name_from_query.upper() != name_in_db.upper():
            continue
        if name_from_query not in hits_by_query:
            hits_by_query[name_from_query] = []
        hits_by_query[name_from_query].append(x)

    # Clean up the hits a bit.
    for name_from_query in hits_by_query:
        hits = hits_by_query[name_from_query]
        
        # If some of the cases match exactly, then use the exact
        # matches.  Do this before the official hits, because this can
        # give clues to the organism.
        exact_matches = [int(x[2]==name_from_query) for x in hits]
        if sum(exact_matches) > 0:
            hits = [x for x in hits if x[2] == name_from_query]

        # If some of these hits are official, and others aren't, then
        # only keep the official ones.
        is_official = [x[3] for x in hits]
        if sum(is_official) > 0:
            hits = [x for x in hits if x[3]]

        hits_by_query[name_from_query] = hits
        
    clean = {}
    for (name_from_query, hits) in hits_by_query.iteritems():
        for x in hits:
            gene_id, source_db, name_in_db, is_official = x
        
            x = _lookup_gene(gene_id)
            id, symbol, name, tax_id_, organism = x
            x = gene_id, symbol, name, tax_id_, organism, name_from_query, \
                source_db, name_in_db
            clean[x] = 1  # no duplicates
    clean = sorted(clean)

    # If tax_id is given, then only return hits from that tax_id.
    if tax_id:
        tax_id = int(tax_id)
        clean = [x for x in clean if x[3] == tax_id]
            
    return clean


def find_discontinued(name, tax_id=None):
    # Return current gene_id or None if not found.
    import config
    import dblib

    # If name looks like a gene_id, see if it's a discontinued gene ID.
    is_gene_id = True
    try:
        int(name)
    except ValueError, x:
        is_gene_id = False

    # Need to clean up name for security reasons.
    columns = "old_gene_id, old_symbol, gene_id, tax_id"
    q = "SELECT %s FROM %s WHERE old_symbol='%s';" % (
        columns, "DISCONTINUED", name)
    if is_gene_id:
        q = "SELECT %s FROM %s WHERE old_gene_id=%s;" % (
            columns, "DISCONTINUED", name)
    x = dblib.query(
        q, config.gm_user, config.gm_passwd, config.gm_db,
        config.gm_host, config.gm_port, config.gm_socket)
    gene_ids = []
    for x in x:
        old_gene_id, old_symbol, gene_id, tax_id_ = x
        if tax_id is not None and int(tax_id) != int(tax_id_):
            continue
        if gene_id not in gene_ids:
            gene_ids.append(gene_id)
    if not gene_ids:
        return None
    assert len(gene_ids) == 1, "Multiple discontinued for: %s" % name
    gene_id = gene_ids[0]
    return int(gene_id)


def _lookup_gene(id):
    """id should be the Entez ID of a gene, given as an integer.
    
    Return id, symbol, name, tax_id, organism.  id and tax_id will be
    integers.  Everything else is a string.  If not found, then
    everything (except for the id) will be empty strings.
    
    """
    import config
    import dblib

    id = int(id)
    symbol = name = tax_id = organism = ""
    table = _find_entrez_gene_table()
    q = "SELECT symbol, name, tax_id, organism FROM %s WHERE gene_id=%s" % (
        table, id)
    x = dblib.query(
        q, config.gm_user, config.gm_passwd, config.gm_db, config.gm_host,
        config.gm_port, config.gm_socket)
    for x in x:
        symbol, name, tax_id, organism = x
    return id, symbol, name, tax_id, organism

def _find_entrez_gene_table():
    # Return the name of the most recent one.
    import config
    import dblib

    x = dblib.query(
        "SHOW TABLES", config.gm_user, config.gm_passwd, config.gm_db,
        config.gm_host, config.gm_port, config.gm_socket)
    x = [x[0] for x in x]
    x = [x for x in x if x.startswith("entrez_gene_")]
    x.sort()
    assert x
    return x[-1]

def _list_mapping_tables():
    import config
    import dblib

    x = dblib.query(
        "SHOW TABLES", config.gm_user, config.gm_passwd, config.gm_db,
        config.gm_host, config.gm_port, config.gm_socket)
    x = [x[0] for x in x]
    x = [x for x in x if x.startswith("MAP_")]
    x = sorted(x)   # will sort by table name, then date.
    mapping_tables = x

    # If multiple versions of the table exists, use the one with the
    # latest date.
    name2tablename = {}  # ENTREZ_ID -> MAP_ENTREZ_ID_111223
    for tablename in mapping_tables:
        # MAP_ENTREZ_ID_111223
        # MAP_GENBANK_111223
        # MAP_AFFY_HG_U133A_2___na32
        assert tablename.startswith("MAP_")
        x = tablename[4:]

        # If last 6 digits is an integer, then remove it.
        try:
            int(x[-6:])
        except ValueError, y:
            pass
        else:
            assert x[-7] == "_"
            x = x[:-7]
        name = x
        
        # later ones will overwrite earlier ones
        name2tablename[name] = tablename
    mapping_tables = sorted(name2tablename.values())
    
    return mapping_tables
