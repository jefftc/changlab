#rulebase.py


def import_rules():
    import os

    RULE_PATH = "bie_rules"
    module_path = os.path.join(os.path.dirname( __file__ ), RULE_PATH)
    assert os.path.exists(module_path), "I could not find the rules: %s" % \
           module_path
    filenames = os.listdir(module_path)
    # filenames is a list of:
    # __init__.py
    # ClassifyFile_rule.py
    # ClassifyFile_rule.pyc
    # GseaFile_rule.pyc
    # [...]

    # Convert to module names, e.g.:
    # ClassifyFile_rule
    # GseaFile_rule
    x = [os.path.splitext(x)[0] for x in filenames]  # No extensions.
    x = {}.fromkeys(x)                               # No duplicates
   
    x = [x for x in x if not x.startswith('._')]
    x = [x for x in x if not x.startswith('__init__')]
    x = [x for x in x if not x.startswith('.DS_Store')]
    #x = [x for x in x if x.endswith("_rule")]        # Must end with _rule.
    #x = [x for x in x if x.endswith(".py")]        # Must end with _rule.
    module_names = x
    # Load each module.
    all_modules = []
    for name in module_names:
        import_name = "bie_rules.%s" % name
        module = __import__(import_name, globals(), locals(), [name], -1)
        all_modules.extend(module.all_modules)
        assert hasattr(module, "list_files"), "Missing list_files var in %s" \
               % repr(module)
        for mod in module.list_files:
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
