# Functions:
# query
# format_list

# _dnslookup
# _connect_mysql
# _disconnect_mysql

def query(query, user, passwd, db, host):
    import time
    import MySQLdb

    start = time.time()
    num_tries = 0
    while 1:
        DB = _connect_mysql(user, passwd, db, host)
        cursor = DB.cursor()
        num_tries += 1
        try:
            cursor.execute(query)
        except MySQLdb.OperationalError, x:
            if (str(x).find("MySQL server has gone away") >= 0 and
                num_tries < 5):
                _disconnect_mysql(user, db, host)
                continue
            elif str(x).find("packet bigger than 'max_allowed_packet'") >= 0:
                x = "%s\n%s\n%d bytes" % (str(x), query, len(query))
                raise MySQLdb.OperationalError, x
            else:
                raise
        results = cursor.fetchall()
        cursor.close()
        break
    #DB.close()   # don't disconnect cached connection
    #_disconnect_mysql()
    end = time.time()
    #debuglog(query)
    #print "QUERY %g: %s" % (end-start, query)
    return results

def format_list(x, quote=True):
    assert x, "empty list!"
    if quote:
        x = ['"%s"' % x for x in x]
    else:
        x = map(str, x)
    return ", ".join(x)

def _dnslookup(hostname):
    import socket
    return socket.gethostbyname(hostname)

DB_CACHE = {}  # (user, db, host) -> DB object
def _connect_mysql(user, passwd, db, host):
    # Return a DB object.
    global DB_CACHE
    import time
    import re
    import MySQLdb

    if host is None:
        host = ""

    # If host is not an IP address, convert it to one to prevent
    # repeated DNS lookups.
    if host and not re.match("^[.\d]+$", host):
        host = _dnslookup(host)

    key = user, db, host
    if key in DB_CACHE:
        return DB_CACHE[key]

    known_db_errors = [
        "Lost connection to MySQL server",
        # Don't ignore this.  Can mask if server really isn't running.
        #"Can't connect to MySQL server",
        ]
    # Optimization, try to reuse connections.  Broken: I eventually
    # will get an exception.  Doesn't matter -- mysql connections are
    # really fast.
    # MySQLdb.OperationalError: (1040, 'Too many connections')
    # Optimization: Using the raw _mysql calls not much faster.
    start = time.time()
    while 1:
        try:
            DB = MySQLdb.connect(user=user, passwd=passwd, db=db, host=host)
        except MySQLdb.OperationalError, x:
            for err in known_db_errors:
                if str(x).find(err) >= 0:
                    break
            else:
                raise
        else:
            break
        if time.time() >= start+300:   # try for 5 minutes
            raise AssertionError, "Database timed out"
    DB_CACHE[key] = DB
    return DB

def _disconnect_mysql(user, db, host):
    global DB_CACHE
    key = user, db, host
    if key not in DB_CACHE:
        return
    DB_CACHE[key].close()
    del DB_CACHE[key]
