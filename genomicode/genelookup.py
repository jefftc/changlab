"""
Functions:
find_gene               Look up the one best gene for a name.
find_all_genes          Return all possible matches for one name.
find_all_genes_verbose  Return all possible matches, lots of info. 

"""
# _lookup_gene
# _find_entrez_gene_table
# _list_mapping_tables

def find_gene(name):
    # Return gene_id, symbol, name, tax_id, organism, name_from_query or None
    results = find_all_genes(name)
    if not results:
        return None
    x = [x[1] for x in results]
    x = ", ".join(x)
    assert len(results) == 1, "Multiple matches found: %s" % x
    return results[0]

def find_all_genes(name):
    # Return list of gene_id, symbol, name, tax_id, organism, name_from_query.
    x = find_all_genes_verbose(name)
    results = [(x[0], x[1], x[2], x[3], x[4], x[5]) for x in x]
    results.sort()
    # No duplicates.
    i = 0
    while i < len(results)-1:
        if results[i] == results[i+1]:
            del results[i]
        else:
            results += 1
    return results
    
def find_all_genes_verbose(name):
    # Return list of gene_id, symbol, name, tax_id, organism,
    # name_from_query, source_db, name_in_db.
    from genomicode import config
    from genomicode import dblib

    case_sensitive = False  # sort out matches later

    table_names = _list_mapping_tables()
    
    results = [] # gene_id, name_in_db, is_official, name_from_query, source_db
    for table_name in table_names:
        columns = "gene_id, symbol, official"
        cs = ""
        if case_sensitive:
            cs = "BINARY"
        q = "SELECT %s FROM %s WHERE %s symbol='%s';" % (
            columns, table_name, cs, name)
        x = dblib.query_db(
            q, config.gm_user, config.gm_passwd, config.gm_db, config.gm_host)
        for x in x:
            id, symbol, official = x
            id = int(id)
            x = id, symbol, official, name, table_name
            results.append(x)

    # If some of the cases match exactly, then use the exact matches.
    # Do this before the official hits, because this can give clues to
    # the organism.
    exact_matches = [int(x[1]==name) for x in results]
    if sum(exact_matches) > 0:
        results = [x for x in results if x[1] == name]

    # If some of these hits are official, and others aren't, then only
    # keep the official ones.
    is_official = [x[2] for x in results]
    if sum(is_official) > 0:
        results = [x for x in results if x[2]]

    clean = []
    for x in results:
        gene_id, name_in_db, is_official, name_from_query, source_db = x
        x = _lookup_gene(gene_id)
        id, symbol, name, tax_id, organism = x
        x = gene_id, symbol, name, tax_id, organism, name_from_query, \
            source_db, name_in_db
        clean.append(x)
    return clean

def _lookup_gene(id):
    # Return id, symbol, name, tax_id, organism.  If not found, then
    # everything (except for the id) will be empty strings.
    from genomicode import config
    from genomicode import dblib

    symbol = name = tax_id = organism = ""
    table = _find_entrez_gene_table()
    q = "SELECT symbol, name, tax_id, organism FROM %s WHERE gene_id=%d" % (
        table, id)
    x = dblib.query_db(
        q, config.gm_user, config.gm_passwd, config.gm_db, config.gm_host)
    for x in x:
        symbol, name, tax_id, organism = x
    return id, symbol, name, tax_id, organism

def _find_entrez_gene_table():
    # Return the name of the most recent one.
    from genomicode import config
    from genomicode import dblib

    x = dblib.query_db(
        "SHOW TABLES", config.gm_user, config.gm_passwd, config.gm_db,
        config.gm_host)
    x = [x[0] for x in x]
    x = [x for x in x if x.startswith("entrez_gene_")]
    x.sort()
    assert x
    return x[-1]

def _list_mapping_tables():
    from genomicode import config
    from genomicode import dblib

    x = dblib.query_db(
        "SHOW TABLES", config.gm_user, config.gm_passwd, config.gm_db,
        config.gm_host)
    x = [x[0] for x in x]
    x = [x for x in x if x.startswith("MAP_")]
    return x
