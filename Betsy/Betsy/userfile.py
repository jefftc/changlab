"""

Functions:
set               Save a file for a user.
get               Retrieve a file for a user.
get_by_checksum   Retrieve a file based on checksum.
dir               List the files owned by the user.

"""
# _make_path
# _hash_storefile
# _unhash_storefile


def set(username, filename):
    import os
    import shutil
    import tempfile

    assert os.path.exists(filename), '%s does not exists' % filename
    
    out_path = _make_path(username)
    out_file = _hash_storefile(username, filename)
    out_filename = os.path.join(out_path, out_file)
    
    # If this file already exists, don't copy it again.
    if os.path.exists(out_filename):
        return out_filename

    # Copy to a temporary path and then move it once it's done
    # copying.  This prevents partial copies.
    temp_filename = None
    try:
        x, temp_filename = tempfile.mkstemp(dir=out_path); os.close(x)
        if os.path.exists(temp_filename):
            os.unlink(temp_filename)
        
        if os.path.isdir(filename):
            shutil.copytree(filename, temp_filename)
        else:
            shutil.copy(filename, temp_filename)
        os.rename(temp_filename, out_filename)
    finally:
        if temp_filename is not None and os.path.exists(temp_filename):
            if os.path.isdir(temp_filename):
                shutil.rmtree(temp_filename)
            else:
                os.unlink(temp_filename)
    return out_filename


def get(username, filename, file_length=None):
    import os
    
    user_path = _make_path(username)
    storefiles = os.listdir(user_path)
    for storefile in storefiles:
        x = _unhash_storefile(storefile)
        username, in_filename, checksum, in_file_length = x
        if in_filename != filename:
            continue
        if file_length is not None and in_file_length != file_length:
            continue
        return os.path.join(user_path, storefile)
    raise KeyError, filename


def get_by_checksum(username, checksum, file_length):
    import os
    
    user_path = _make_path(username)
    storefiles = os.listdir(user_path)
    for storefile in storefiles:
        x = _unhash_storefile(storefile)
        username, in_filename, in_checksum, in_file_length = x
        if in_checksum == checksum and in_file_length == file_length:
            return os.path.join(user_path, storefile)
    raise ValueError(
        'cannot find the file with checksum %s and file_length %s' %
        (checksum, file_length))


def dir(username):
    import os
    
    user_path = _make_path(username)
    filenames = os.listdir(user_path)
    filenames = [i for i in filenames if not i.startswith('.')]
    result = []
    for storefile in filenames:
        x = _unhash_storefile(storefile)
        result.append(list(x))
    return result


def _make_path(username):
    # Return the full path to store files for this user.
    # <config.OUTPUTPATH>/userfile/<username>
    import os
    import config

    p1 = os.path.join(config.OUTPUTPATH, 'userfile')
    p2 = os.path.join(p1, username)
    if not os.path.exists(p1):
        os.mkdir(p1)
    if not os.path.exists(p2):
        os.mkdir(p2)
    return p2


def _hash_storefile(username, filename):
    # Return a unique name or path to save this file.
    import os
    from genomicode import filelib
    from Betsy import bhashlib
    
    assert '___' not in username
    assert '___' not in filename
    assert '__' not in filename

    #file_length = os.path.getsize(filename)
    file_length = filelib.get_file_or_path_size(filename)
    checksum = bhashlib.checksum_file_or_path_smart(filename)

    x = filename.replace("/", "__")
    # Put the file at the end to preserve the file extension.
    x = [username, checksum, file_length, x]
    x = "___".join(map(str, x))
    return x
    #newfilename = filename.replace('/', '__')
    #store_name = str(username + '___' + newfilename + '___' + str(checksum) +
    #                 '___' + str(file_length))
    #return store_name


def _unhash_storefile(storefilename):
    if '___' not in storefilename:
        return storefilename
    x = storefilename.split('___')
    username, checksum, file_length, filename = x
    real_filename = filename.replace('__', '/')
    return username, real_filename, checksum, file_length
