#! /usr/bin/env python
#run_rule.py
import sys
import os
import shutil
import config
import hash_method


"""
Functions:
store(username,filename)
retrieve(username,filename)
retrieve_by_checksum(username,checksum,file_length)
list(user_name)

"""

output_path = os.path.join(config.OUTPUTPATH, 'userfile')


def store(username, filename):
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    user_path = os.path.join(output_path, username)
    if not os.path.exists(user_path):
        os.mkdir(user_path)
    size = os.path.getsize(filename)
    checksum = hash_method.get_input_checksum(filename)
    newfilename = filename.replace('/', '_')
    new_name = str(username + '_' + newfilename + '_' +
                   str(checksum) + '_' + str(size))
    new_file_path = os.path.join(user_path, new_name)
    if os.path.exists(new_file_path):
        return new_file_path
    if os.path.isdir(filename):
        shutil.copytree(filename, new_file_path)
    else:
        shutil.copy(filename, new_file_path)
    return new_file_path


def retrieve(username, filename):
    user_path = os.path.join(output_path, username)
    assert os.path.exists(user_path), 'no file for username %s' % username
    storefiles = os.listdir(user_path)
    newfilename = filename.replace('/', '_')
    for storefile in storefiles:
        if storefile.startswith(username + '_' + newfilename):
            return os.path.join(user_path, storefile)
    raise ValueError('cannot find the file %s' % filename)


def retrieve_by_checksum(username, checksum, file_length):
    user_path = os.path.join(output_path, username)
    assert os.path.exists(user_path), 'no file for username %s' % username
    storefiles = os.listdir(user_path)
    for storefile in storefiles:
        endstring = str(checksum) + '_' + str(file_length)
        if storefile.startswith(username) and storefile.endswith(endstring):
            return os.path.join(user_path, storefile)
    raise ValueError('cannot find the file with checksum %s and file_length %s'
                     % (checksum, file_length))


def list(username):
    user_path = os.path.join(output_path, username)
    filenames = os.listdir(user_path)
    filenames = [i for i in filenames if not i.startswith('.')]
    result = []
    for filename in filenames:
        items = filename.split('_')
        result.append([items[0], '_'.join(items[1:-2]), items[-2], items[-1]])
    return result
