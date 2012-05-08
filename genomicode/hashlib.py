"""

Functions:
hash_R               Hash a string using the R algorithm for list names.
hash_var             Hash a string to an acceptable variable name.
hash_geneid          Lowercase, no spaces at the end.
hash_sampleid        Lowercase, no punctuation (except _), no initial X.

hash_many_geneids
hash_many_sampleids

"""

def hash_R(s):
    # Hash a string using the R algorithm for list names.
    import re

    #s_orig = s
    # If the string starts with a non word character, prepend an "X".
    if re.match(r"[^a-zA-Z]", s):
        s = "X%s" % s
    # Convert all punctuation (except for _) to ".".
    s = re.sub(r"\W", ".", s)
    return s

def hash_var(name):
    import re
    # Fix the header to be a python variable.
    x = str(name)
    # Replace all non-word character with _.
    x = re.sub(r"\W", "_", x)
    # Replace initial numbers with Xnumber.
    x = re.sub(r"^(\d)", r"X\1", x)
    return x

def hash_geneid(id_):
    return id_.strip().lower()

def hash_sampleid(id_):
    # Hash the sample names so that small differences are ignored.  R
    # will change the sample names, so if one data set has been
    # through R but not the other, the names will be different.
    # 2554_6933_32492_Mock1_HG-U133A+2
    # X2554_6933_32492_Mock1_HG.U133A.2
    import re

    x = id_

    # If there are alphanumeric characters, then assume that
    # punctuation isn't meaningful.
    if re.search(r"\w", x):
        # Change all non-words to '.', like R does.  (This does not
        # change underscores.)
        x = re.subn(r"\W", ".", x)[0]

    # Ignore initial X.
    if re.match(r"X[\d\w]", x):
        x = x[1:]

    # Make case insensitive.
    x = x.lower()

    return x

def hash_many_geneids(ids):
    #return [_hash_geneid(x) for x in ids]
    # Optimization: do this without a function call.
    return [x.strip().lower() for x in ids]

def hash_many_sampleids(ids):
    return [hash_sampleid(x) for x in ids]
