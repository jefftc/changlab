"""

Functions:
set(username, filename)
get(username, filename, file_length=None)
get_by_checksum(username, checksum, file_length)
dir(user_name)

"""
# _make_user_path(username)
# _hash_storefile(username,filename, checksum,file_length)
# _unhash_storefile(storefilename)



def _make_user_path(username):
    # Return the full path to store files for this user.
    # <config.OUTPUTPATH>/userfile/<username>
    import os
    import config

    p1 = os.path.join(config.OUTPUTPATH, 'userfile')
    if not os.path.exists(p1):
        os.mkdir(p1)
    # TODO: need to hash username
    p2 = os.path.join(p1, username)
    if not os.path.exists(p2):
        os.mkdir(p2)
    return p2


def _hash_storefile(username, filename):
    import os
    import hash_method
    
    assert '___' not in username
    assert '___' not in filename
    assert '__' not in filename

    file_length = os.path.getsize(filename)
    checksum = hash_method.get_input_checksum(filename)

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


def set(username, filename):
    import os
    import shutil

    assert os.path.exists(filename), '%s does not exists' % filename
    
    user_path = _make_user_path(username)
    
    store_name = _hash_storefile(username, filename)
    new_file_path = os.path.join(user_path, store_name)
    if os.path.exists(new_file_path):
        return new_file_path
    if os.path.isdir(filename):
        shutil.copytree(filename, new_file_path)
    else:
        shutil.copy(filename, new_file_path)
    return new_file_path


def get(username, filename, file_length=None):
    import os
    
    user_path = _make_user_path(username)
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
    
    user_path = _make_user_path(username)
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
    
    user_path = _make_user_path(username)
    filenames = os.listdir(user_path)
    filenames = [i for i in filenames if not i.startswith('.')]
    result = []
    for storefile in filenames:
        x = _unhash_storefile(storefile)
        result.append(list(x))
    return result
