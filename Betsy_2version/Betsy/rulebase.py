#rulebase.py
import os
module_lib = 'bie_rules'
file_names = os.listdir(module_lib)
all_modules = []
for file_name in file_names:
    module_name = os.path.splitext(file_name)[0]
    if module_name.endswith('rule'):
        module = __import__('bie_rules.'+module_name, globals(),
                        locals(), [module_name], -1)
        all_modules.extend(module.all_modules)
        for i in module.list_files:
            vars()[i.name] = i

