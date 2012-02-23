#Betsy_config.py
import os
import sys
import ConfigParser

#write the config file
##config = ConfigParser.RawConfigParser()
##config.add_section('Section1')
##config.set('Section1', 'PREPROCESS', '/home/xchen/chencode/scripts/preprocess.py')
##config.set('Section1', 'ARRAYPLOT', '/home/xchen/chencode/scripts/arrayplot.py')
##with open('Betsy.cfg','wb') as configfile:
##    config.write(configfile)

DEFAULTS={'PREPROCESS': 'preprocess.py', 'ARRAYPLOT': 'arrayplot.py','XLS2TXT':'xls2txt'}
SEARCH_PATH = [
        '/home/xchen/chencode/scripts',
        "/usr/local/bin",
        "/usr/bin",
        '/opt/local/bin',
        os.environ['HOME']+'/bin'
        '/home/changlab/changlab/scripts'
        ]

def read_config():
    #config_file='Betsy.cfg'
    config_file = os.environ["HOME"]+'/.betsyrc'
    var_dict = dict()
    if os.path.exists(config_file):
        config = ConfigParser.ConfigParser()
        config.read(config_file)
        sections = config.sections()
        for i in sections:
            section_content = config.items(i)
            for j in range(len(section_content)):
                var_dict[section_content[j][0].upper()] = section_content[j][1]
    for keys in DEFAULTS.keys():
        if keys not in var_dict.keys():
            var_dict[keys] = DEFAULTS[keys]
            for path in SEARCH_PATH:
                filepath=os.path.join(path,DEFAULTS[keys])
                if os.path.exists(filepath):
                    var_dict[keys] = filepath
                    break
             
    return var_dict


var_dict = read_config()
for keys in var_dict.keys():
    vars()[keys]=var_dict[keys]

del var_dict



