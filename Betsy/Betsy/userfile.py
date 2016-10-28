"""

Functions:
set            Cache a file for a user, returning the cache filename.
has_file

"""
# _make_user_path
# _make_cache_path


# File path:
# <config.OUTPATH>/userfile/<username>/           user_path
#   <uid>___<filename>/                           cache_path
#     <filename>     # can be file or path            (always there)
#     hash_smart.txt # hash based on contents or size (always there)
#     hash_full.txt  # full hash of file contents     (can be missing)

def set(username, filename):
    import os
    import shutil
    import tempfile
    from Betsy import bhashlib

    assert os.path.exists(filename), "File not found: %s" % filename
    filename = os.path.realpath(filename)  # don't cache symlink
    cache_filename = has_file(username, filename)
    if cache_filename:
        return cache_filename
    
    # Make the cache path.
    user_path = _make_user_path(username)
    cache_path = _make_cache_path(username, filename)
    assert not os.path.exists(cache_path)  # should be unique
    
    # Make the name of the file that will exist.
    p, file_ = os.path.split(filename)
    assert file_ not in ["hash_smart.txt", "hash_full.txt"]
    data_filename = os.path.join(cache_path, file_)
    hash_filename = os.path.join(cache_path, "hash_smart.txt")

    # Copy to a temporary path and then move it once it's done
    # copying.  This prevents partial copies.
    temp_path = None
    try:
        x, temp_path = tempfile.mkstemp(dir=user_path); os.close(x)
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        os.mkdir(temp_path)
        
        temp_filename = os.path.join(temp_path, file_)
        if os.path.isdir(filename):
            shutil.copytree(filename, temp_filename)
        else:
            shutil.copy(filename, temp_filename)

        # Make the smart hash file.
        smart_hash = bhashlib.checksum_file_or_path_smart(filename)
        x = os.path.join(temp_path, "hash_smart.txt")
        open(x, 'w').write(smart_hash)
            
        # If the file exists now, then assume it's correct.  This can
        # happen if two processes are trying to generate this file at
        # the same time.
        if not os.path.exists(cache_path):
            os.rename(temp_path, cache_path)
    finally:
        if temp_path is not None and os.path.exists(temp_path):
            assert os.path.isdir(temp_path)
            shutil.rmtree(temp_path)

    assert os.path.exists(data_filename)
    assert os.path.exists(hash_filename)
    return data_filename


def has_file(username, filename):
    # Return the name of the file.
    import os
    from Betsy import bhashlib
    
    user_path = _make_user_path(username)
    assert os.path.exists(filename), "Not found: %s" % filename
    filename = os.path.realpath(filename)  # don't cache symlink
    file_smart_hash = bhashlib.checksum_file_or_path_smart(filename)

    # Look for files that match the smart hash.  If I find any, then
    # do a hash of the full contents of the file.
    found = []  # list of "<uid>___<filename>"
    paths = os.listdir(user_path)
    for path in paths:
        # If a previous "set" operation was interrupted, there might
        # be temp directories leftover with no hash files.  Ignore
        # directories with no hash files.
        cache_file = os.path.join(user_path, path, "hash_smart.txt")
        if not os.path.exists(cache_file):
            # delete this directory?
            continue
        #assert os.path.exists(cache_file), "File not found: %s" % cache_file
        cache_smart_hash = open(cache_file).read()
        if file_smart_hash == cache_smart_hash:
            found.append(path)
    if not found:
        return None

    # For everything found, see if it matches hash.txt.
    # Not implemented.  Just use hash_fast.txt.  Can result in
    # collisions.

    # BUG:
    # There might be multiple files that match.  This can happen if
    # there is a collision, or if two processes generated this file at
    # the same time.  Assume no collisions, and use the first one.
    #assert len(found) == 1, "collision"
    assert len(found) >= 1
    cache_path = os.path.join(user_path, found[0])
    # Look for <filename> in cache_path.
    x = os.listdir(cache_path)
    x = [x for x in x if x not in ["hash_smart.txt", "hash_full.txt"]]
    assert len(x) == 1
    cache_filename = os.path.join(cache_path, x[0])
    return cache_filename
    

## def get_by_checksum(username, checksum, file_length):
##     import os
    
##     user_path = _make_path(username)
##     storefiles = os.listdir(user_path)
##     for storefile in storefiles:
##         x = _unhash_storefile(storefile)
##         username, in_filename, in_checksum, in_file_length = x
##         if in_checksum == checksum and in_file_length == file_length:
##             return os.path.join(user_path, storefile)
##     raise ValueError(
##         'cannot find the file with checksum %s and file_length %s' %
##         (checksum, file_length))


## def dir(username):
##     import os
    
##     user_path = _make_path(username)
##     filenames = os.listdir(user_path)
##     filenames = [i for i in filenames if not i.startswith('.')]
##     result = []
##     for storefile in filenames:
##         x = _unhash_storefile(storefile)
##         result.append(list(x))
##     return result


def _make_user_path(username):
    # Return the full path to store files for this user.
    # <config.CACHE_PATH>/userfile/<username>
    import os
    import config

    p1 = os.path.join(config.CACHE_PATH, 'userfile')
    p2 = os.path.join(p1, username)
    if not os.path.exists(p1):
        os.mkdir(p1)
    if not os.path.exists(p2):
        os.mkdir(p2)
    return p2


## def _hash_storefile(username, filename):
##     # Return a unique name or path to save this file.
##     import os
##     from genomicode import filelib
##     from Betsy import bhashlib
    
##     assert '___' not in username
##     assert '___' not in filename
##     assert '__' not in filename

##     #file_length = os.path.getsize(filename)
##     file_length = filelib.get_file_or_path_size(filename)
##     checksum = bhashlib.checksum_file_or_path_smart(filename)

##     x = filename.replace("/", "__")
##     # Put the file at the end to preserve the file extension.
##     x = [username, checksum, file_length, x]
##     x = "___".join(map(str, x))
##     return x
##     #newfilename = filename.replace('/', '__')
##     #store_name = str(username + '___' + newfilename + '___' + str(checksum) +
##     #                 '___' + str(file_length))
##     #return store_name


## def _unhash_storefile(storefilename):
##     if '___' not in storefilename:
##         return storefilename
##     x = storefilename.split('___')
##     username, checksum, file_length, filename = x
##     real_filename = filename.replace('__', '/')
##     return username, real_filename, checksum, file_length

def _make_cache_path(username, filename):
    import tempfile
    # <uid>___<filename>
    import os
    import tempfile

    user_path = _make_user_path(username)
    uid = None
    try:
        x, uid = tempfile.mkstemp(dir=user_path, prefix=""); os.close(x)
    finally:
        if uid:
            x = os.path.join(user_path, uid)
            if os.path.exists(x):
                os.unlink(x)
    
    p, file_ = os.path.split(filename)
    x = "%s___%s" % (uid, file_)
    x = os.path.join(user_path, x)
    return x
