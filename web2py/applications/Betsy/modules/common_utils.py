#!/usr/bin/env python
# coding: utf8
from gluon import *
processing_info='/home/xchen/chencode/tmp/processing_info/'
def hash_command(time, command_line):
    from hashlib import md5
    hashstring=time+command_line
    hash = md5()
    hash.update(hashstring)
    hash_result = hash.hexdigest()
    return hash_result
