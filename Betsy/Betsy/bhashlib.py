# Functions:
# checksum_file    Checksum a file.  If fast, just use size and creation date.
# checksum_path    Checksum each file in the path.  BUG: what about path names?
# checksum_file_or_path         Checksum either a file or path.
# checksum_file_or_path_smart   Will do fast checksum if files are too big.

CHUNK_SIZE = 1024*1024

def checksum_file(filename, fast=False):
    import os
    import stat
    from hashlib import md5
    
    assert os.path.exists(filename)
    # Make sure I don't hash a symlink.
    filename = os.path.realpath(filename)
    
    hasher = md5()

    # If fast, just checksum based on size and creation date.
    # Otherwise, checksum the entire contents.
    if fast:
        x = os.stat(filename)
        x = "%d %d" % (x[stat.ST_SIZE], x[stat.ST_MTIME])
        hasher.update(x)
        return hasher.hexdigest()
    
    handle = open(filename)
    while True:
        x = handle.read(CHUNK_SIZE)
        if not x:
            break
        hasher.update(x)
    return hasher.hexdigest()


def checksum_path(path, fast=False):
    import os
    from hashlib import md5
    from genomicode import filelib
    
    hasher = md5()

    # Checksum each file.
    filenames = filelib.list_files_in_path(path)
    for filename in filenames:
        x = checksum_file(filename, fast=fast)
        hasher.update(x)
        
    return hasher.hexdigest()


def checksum_file_or_path(file_or_path, fast=False):
    import os

    #size = get_file_or_path_size(file_or_path)
    if os.path.isdir(file_or_path):
        return checksum_path(file_or_path, fast=fast)
    return checksum_file(file_or_path, fast=fast)


def checksum_file_or_path_smart(file_or_path):
    # Returns a checksum.  If the files are too big, does a fast
    # checksum.
    from genomicode import filelib

    size = filelib.get_file_or_path_size(file_or_path)
    # Do a fast checksum if files are over 128 Mb.
    fast = size > 1024*1024*128
    return checksum_file_or_path(file_or_path, fast=fast)
