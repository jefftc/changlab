### CONFIGURATION

import os

opj = os.path.join

HOME = os.environ["HOME"]
PROJECTS = opj(HOME, "projects")
CLAB_HOME = "/home/changlab"

SEARCH_PATH = [
    opj(HOME, "bin"),
    "/usr/local/bin",
    "/usr/bin",
    ]

def read_config():
    import ConfigParser
    
    config_file = os.environ["HOME"]+'/.genomicoderc'

    var_dict = dict()
    if not os.path.exists(config_file):
        return var_dict

    config = ConfigParser.ConfigParser()
    config.optionxform = str   # use case sensitive option names
    config.read(config_file)
    for section in config.sections():
        for x in config.items(section):
            name, value = x
            var_dict[name] = value
    return var_dict

var_dict = read_config()
for name, value in var_dict.iteritems():
    vars()[name] = value
del var_dict
