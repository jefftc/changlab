"""

Functions:
hash_R    Hash a string using the R algorithm for list names.

"""

def hash_R(s):
    # Hash a string using the R algorithm for list names.
    import re

    s_orig = s
    # If the string starts with a non word character, prepend an "X".
    if re.match(r"[^a-zA-Z]", s):
        s = "X%s" % s
    # Convert all punctuation (except for _) to ".".
    s = re.sub(r"\W", ".", s)
    return s
