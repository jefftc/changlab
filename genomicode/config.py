### CONFIGURATION

#opj = os.path.join
#
#HOME = os.environ["HOME"]
#PROJECTS = opj(HOME, "projects")
#CLAB_HOME = "/home/changlab"
#
#SEARCH_PATH = [
#    opj(HOME, "bin"),
#    "/usr/local/bin",
#    "/usr/bin",
#    ]

def read_config():
    import os
    import ConfigParser
    
    config_file = os.environ["HOME"]+'/.genomicoderc'
    assert os.path.exists(config_file), \
           "I could not find the configuration file %s." % (config_file)

    # Read the configuration.
    config = ConfigParser.ConfigParser()
    config.optionxform = str   # use case sensitive option names
    config.read(config_file)

    # Set a dictionary of name=value from the configuration file,
    # ignoring section headings.
    var_dict = {}
    for section in config.sections():
        for x in config.items(section):
            name, value = x
            var_dict[name] = value
    return var_dict


var_dict = read_config()
for name, value in var_dict.iteritems():
    vars()[name] = value
del var_dict
