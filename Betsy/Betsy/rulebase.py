#rulebase.py


def import_rules():
    import os

    RULE_PATH = "rules"
    module_path = os.path.join(os.path.dirname(__file__), RULE_PATH)
    assert os.path.exists(module_path), \
           "I could not find the rules: %s" % module_path
    filenames = os.listdir(module_path)
    # filenames is a list of:
    # __init__.py
    # ClassifyFile.py
    # ClassifyFile.pyc
    # GseaFile.pyc
    # [...]

    # Convert to module names, e.g.:
    # ClassifyFile
    # GseaFile
    x = filenames
    x = [x for x in x if not x.endswith("~")]
    x = [os.path.splitext(x)[0] for x in x]  # No extensions.
    x = {}.fromkeys(x)  # No duplicates

    x = [x for x in x if not x.startswith('._')]
    x = [x for x in x if not x.startswith('__init__')]
    x = [x for x in x if not x.startswith('.DS_Store')]
    x = [x for x in x if not x.startswith(".")]
    x = [x for x in x if not x.startswith("#")]
    #x = [x for x in x if x.endswith("_rule")]        # Must end with _rule.
    #x = [x for x in x if x.endswith(".py")]        # Must end with _rule.
    module_names = x

    # Load each module.
    all_modules = []
    for name in module_names:
        import_name = "rules.%s" % name
        #print "IMPORTING %s" % import_name
        module = __import__(import_name, globals(), locals(), [name], -1)
        all_modules.extend(module.all_modules)

        # Save all the data types in the global namespace.
        assert hasattr(module, "all_data_types"), \
               "Missing all_data_types var in %s" % repr(module)
        for mod in module.all_data_types:
            globals()[mod.name] = mod

    # BUG: Will load modules multiple times.
    #for file_name in filenames:
    #    module_name = os.path.splitext(file_name)[0]
    #    if module_name.endswith('rule'):
    #        module = __import__('bie_rules.'+module_name, globals(),
    #                            locals(), [module_name], -1)
    #        all_modules.extend(module.all_modules)
    #        for i in module.list_files:
    #            vars()[i.name] = i
    return all_modules


all_modules = import_rules()
