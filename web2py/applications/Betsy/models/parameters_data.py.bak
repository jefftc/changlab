# coding: utf8
import os
from Betsy import protocol_utils
protocol_path='/home/xchen/chencode/Betsy/Betsy/protocols/'
labels=[]
label_showname=dict()
protocol_files = os.listdir(protocol_path)
for protocol_file in protocol_files:
     if protocol_file.endswith('.py') and not protocol_file == '__init__.py':
         label = os.path.splitext(protocol_file)[0]
         module = protocol_utils.import_protocol(label)
         labels.append(label)
         label_showname[label]=module.PRETTY_NAME
