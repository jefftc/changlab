

#SEARCH_PATH = [
#    "/usr/local/bin",
#    "/usr/bin",
#    '/opt/local/bin',
#    os.path.join(os.environ['HOME'], 'bin'),
#    '/home/changlab/changlab/scripts',
#    ]


def read_config():
    import os
    import ConfigParser

    config_file = os.path.join(os.environ["HOME"], ".betsyrc")
    assert os.path.exists(config_file), "File not found: %s" % config_file
    
    # Read the configuration.
    config = ConfigParser.ConfigParser()
    config.optionxform = str   # use case sensitive option names
    config.read(config_file)

    # Set a dictionary of name=value from the configuration file,
    # ignoring section headings.
    var_dict = dict()
    for section in config.sections():
        for (name, value) in config.items(section):
            var_dict[name] = value
    
    #config = ConfigParser.ConfigParser()
    #config.read(config_file)
    #sections = config.sections()
    #for i in sections:
    #    section_content = config.items(i)
    #    for j in range(len(section_content)):
    #        var_dict[section_content[j][0].upper()] = section_content[j][1]
                
    #for key in DEFAULTS.keys():
    #    if key not in var_dict.keys():
    #        var_dict[key] = DEFAULTS[key]
    #        for path in SEARCH_PATH:
    #            filepath = os.path.join(path, DEFAULTS[key])
    #            if os.path.exists(filepath):
    #                var_dict[key] = filepath
    #                break
    return var_dict


var_dict = read_config()
for keys in var_dict.keys():
    vars()[keys] = var_dict[keys]
del var_dict
