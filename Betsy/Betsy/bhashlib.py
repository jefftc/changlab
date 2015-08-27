# Functions:
# checksum_file_or_path
# checksum_file
# checksum_path

CHUNK_SIZE = 1024*1024

def checksum_path(path):
    import os
    from hashlib import md5
    
    hasher = md5()

    all_filenames = []
    for x in os.walk(path):
        dirpath, dirnames, files = x
        x = [os.path.join(dirpath, x) for x in files]
        all_filenames.extend(x)
    for filename in all_filenames:
        x = checksum_file(filename)
        hasher.update(x)
    return hasher.hexdigest()


def checksum_file(filename):
    from hashlib import md5
    
    hasher = md5()
    handle = open(filename)
    while True:
        x = handle.read(CHUNK_SIZE)
        if not x:
            break
        hasher.update(x)
    return hasher.hexdigest()


def checksum_file_or_path(file_or_path):
    import os
    if os.path.isdir(file_or_path):
        return checksum_path(file_or_path)
    return checksum_file(file_or_path)
