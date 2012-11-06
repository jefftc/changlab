#hash_method.py
import os


def hash_parameters(inputid, pipeline, **parameters):
    """given a file parameters,generate a hash string"""
    from hashlib import md5
    hashstring = inputid + ' '.join(pipeline)
    for key in sorted(parameters):
        hashstring += str(parameters[key])
    hash = md5()
    hash.update(hashstring)
    hash_result = hash.hexdigest()
    return hash_result


def get_file_checksum(identifier):
    from hashlib import md5, sha1
    chunk_size = 1048576  # 1024 B * 1024 B = 1048576 B = 1 MB
    file_md5_checksum = md5()
    file_sha1_checksum = sha1()
    with open(identifier, "rb") as f:
        byte = f.read(chunk_size)
        byte_size = len(byte)
        while byte:
            file_md5_checksum.update(byte)
            file_sha1_checksum.update(byte)
            byte = f.read(chunk_size)
            byte_size += len(byte)
    byte_size = str(byte_size)
    md5_checksum = file_md5_checksum.hexdigest()
    sha1_checksum = file_sha1_checksum.hexdigest()
    return  byte_size, md5_checksum, sha1_checksum


def get_input_checksum(identifier):
    from hashlib import md5
    chunk_size = 1048576  # 1024 B * 1024 B = 1048576 B = 1 MB
    file_md5_checksum = md5()
    if os.path.isdir(identifier):
        files = os.listdir(identifier)
        byte_size = 0
        for filename in files:
            with open(os.path.join(identifier, filename), "rb") as f:
                byte = f.read(chunk_size)
                byte_size += len(byte)
                while byte:
                    file_md5_checksum.update(byte)
                    byte = f.read(chunk_size)
                    byte_size += len(byte)
    else:
        with open(identifier, "rb") as f:
            byte = f.read(chunk_size)
            byte_size = len(byte)
            while byte:
                file_md5_checksum.update(byte)
                byte = f.read(chunk_size)
                byte_size += len(byte)
    byte_size = str(byte_size)
    md5_checksum = file_md5_checksum.hexdigest()
    return md5_checksum
