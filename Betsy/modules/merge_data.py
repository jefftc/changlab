#merge_data.py

import module_utils
import os
import arrayio
import hash_method
import rule_engine
from genomicode import binreg,Matrix

def run(parameters,objects,pipeline):
    """merge two signal file to generate a joined signal file"""
    merge_file1,obj1=module_utils.find_object(parameters,objects,'signal_file','merge1,dataset1')
    merge_file2,obj2=module_utils.find_object(parameters,objects,'signal_file','merge2,dataset2')
    assert os.path.exists(merge_file1)
    assert os.path.exists(merge_file2)
    outfile,new_objects = get_outfile(parameters,objects)
    f = file(outfile,'w')
    module_utils.merge_two_files(merge_file1,merge_file2,f)
    f.close()
    assert os.path.exists(outfile)
    module_utils.write_Betsy_parameters_file(parameters,[obj1,obj2],pipeline)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','merge1,dataset1')


def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,
            objects,'signal_file','merge1,dataset1','signal_file')
    


