#! /usr/bin/env python
#userfile.py
import sys
import os
import shutil
import config
import hash_method


"""
Functions:
set(username,filename)
get(username,filename,file_length=None)
get_by_checksum(username,checksum,file_length)
dir(user_name)
_make_user_path(username)
_hash_storefile(username,filename, checksum,file_length)
_unhash_storefile(storefilename)
"""



def _make_user_path(username):
    output_path = os.path.join(config.OUTPUTPATH, 'userfile')
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    user_path = os.path.join(output_path, username)
    if not os.path.exists(user_path):
        os.mkdir(user_path)
    return user_path


def _hash_storefile(username,filename, checksum,file_length):
    assert '___' not in username
    assert '___' not in filename
    assert '__' not in filename
    newfilename = filename.replace('/', '__')
    store_name = str(username + '___' + newfilename + '___' +
                   str(checksum) + '___' + str(file_length))
    return store_name

def _unhash_storefile(storefilename):
    x = storefilename.split('___')
    username, filename, checksum, file_length = x
    real_filename = filename.replace('__','/')
    return username, real_filename, checksum, file_length


def set(username, filename):
    user_path = _make_user_path(username)
    size = os.path.getsize(filename)
    checksum = hash_method.get_input_checksum(filename)
    store_name = _hash_storefile(username, filename, checksum, size)
    new_file_path = os.path.join(user_path, store_name)
    if os.path.exists(new_file_path):
        return new_file_path
    if os.path.isdir(filename):
        shutil.copytree(filename, new_file_path)
    else:
        shutil.copy(filename, new_file_path)
    return new_file_path


def get(username, filename, file_length=None):
    user_path = _make_user_path(username)
    storefiles = os.listdir(user_path)
    for storefile in storefiles:
        x = _unhash_storefile(storefile)
        username, in_filename, checksum, in_file_length = x
        if in_filename == filename:
            if file_length and in_file_length != file_length:
                continue
            return os.path.join(user_path, storefile) 
    raise ValueError('cannot find the file %s' % filename)


def get_by_checksum(username, checksum, file_length):
    user_path = _make_user_path(username)
    storefiles = os.listdir(user_path)
    for storefile in storefiles:
        x = _unhash_storefile(storefile)
        username, in_filename, in_checksum, in_file_length = x
        if in_checksum == checksum and in_file_length == file_length:
            return os.path.join(user_path, storefile)
    raise ValueError('cannot find the file with checksum %s and file_length %s'
                     % (checksum, file_length))


def dir(username):
    user_path = _make_user_path(username)
    filenames = os.listdir(user_path)
    filenames = [i for i in filenames if not i.startswith('.')]
    result = []
    for storefile in filenames:
        x = _unhash_storefile(storefile)
        result.append(list(x))
    return result


