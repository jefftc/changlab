"""
Functions:
find_gene               Look up the one best gene for a name.
find_gene_detailed      Return many possible matches, lots of info.

find_many_genes
find_many_genes_detailed

"""
# _lookup_gene
# _find_entrez_gene_table
# _list_mapping_tables

# TODO:
# fuzzy match of gene names

def find_gene(name):
    # Return gene_id, symbol, name, tax_id, organism, name_from_query
    # or None.  If name is a RefSeq ID, be sure there is no version
    # number.  E.g. "NP_000680.2" won't be found, but "NP_000680"
    # will.
    orig_name = name
    x = find_many_genes([name])
    assert len(x) == 1
    num_matches, gene_id, symbol, name, tax_id, organism, name_from_query =x[0]
    if num_matches == 0:
        return None
    assert num_matches == 1, 'Multiple matches (%d) found for "%s".' % (
        num_matches, orig_name)
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
    
def find_gene_detailed(name):
    # Return list of gene_id, symbol, name, tax_id, organism,
    # name_from_query, source_db, name_in_db.
    return find_many_genes_detailed([name])

def find_many_genes(genes):
    # Return list of (num_matches, gene_id, symbol, name, tax_id,
    # organism, name_from_query) that is parallel to genes.  If there
    # are no matches or multiple matches, values (other than
    # num_matches and name_from_query) are set to None.
    
    x = find_many_genes_detailed(genes)
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

def find_many_genes_detailed(genes):
    # Return list of gene_id, symbol, name, tax_id, organism,
    #   name_from_query, source_db, name_in_db.
    # May not be in the same order as the query.
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
            q, config.gm_user, config.gm_passwd, config.gm_db, config.gm_host)
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
                config.gm_host)
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
        
    clean = []
    for (name_from_query, hits) in hits_by_query.iteritems():
        for x in hits:
            gene_id, source_db, name_in_db, is_official = x
        
            x = _lookup_gene(gene_id)
            id, symbol, name, tax_id, organism = x
            x = gene_id, symbol, name, tax_id, organism, name_from_query, \
                source_db, name_in_db
            clean.append(x)
    return clean

def _lookup_gene(id):
    # Return id, symbol, name, tax_id, organism.  If not found, then
    # everything (except for the id) will be empty strings.
    import config
    import dblib

    symbol = name = tax_id = organism = ""
    table = _find_entrez_gene_table()
    q = "SELECT symbol, name, tax_id, organism FROM %s WHERE gene_id=%d" % (
        table, id)
    x = dblib.query(
        q, config.gm_user, config.gm_passwd, config.gm_db, config.gm_host)
    for x in x:
        symbol, name, tax_id, organism = x
    return id, symbol, name, tax_id, organism

def _find_entrez_gene_table():
    # Return the name of the most recent one.
    import config
    import dblib

    x = dblib.query(
        "SHOW TABLES", config.gm_user, config.gm_passwd, config.gm_db,
        config.gm_host)
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
        config.gm_host)
    x = [x[0] for x in x]
    x = [x for x in x if x.startswith("MAP_")]
    return x
