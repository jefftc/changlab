# goal
# rule         another rule is proposition
# consequent   second half of a rule
# antecedent   first half of a rule
# fact


class Data:
    def __init__(self, datatype, attributes):
        self.datatype = datatype
        self.attributes = attributes.copy()
    def copy(self):
        return Data(self.datatype, self.attributes.copy())
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            repr(self.datatype),
            repr(self.attributes),
            ]
        return "Data(%s)" % ", ".join(x)

class Module:
    def __init__(self, name, in_data, out_data):
        self.name = name
        self.in_data = in_data.copy()
        self.out_data = out_data.copy()
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            repr(self.name),
            repr(self.in_data),
            repr(self.out_data),
            ]
        x = "Module(%s)" % ", ".join(x)
        return x


def make_data(datatype, **params):
    return Data(datatype, params)

antecedent = make_data
consequent = make_data


all_modules = [
    Module(
        "preprocess_rma",
        antecedent("cel_files", cel_version="v3_4"),
        consequent("signal_file",
                  logged="yes", preprocess="rma", quantile_norm="yes")),
    Module(
        "preprocess_mas5",
        antecedent("cel_files", cel_version="v3_4"),
        consequent("signal_file", logged="no", preprocess="mas5")),
    Module(
        "quantile_norm",
        antecedent("signal_file", quantile_norm=[None, "no"]),
        consequent("signal_file", quantile_norm="yes")),
    Module(
        "log_signal",
        antecedent("signal_file", logged=[None, "no"]),
        consequent("signal_file", logged="yes")),
    Module(
        "gene_center",
        antecedent("signal_file", gene_center=[None, "no"]),
        consequent("signal_file", gene_center=["mean", "median"])),
    Module(
        "gene_normalize",
        antecedent("signal_file", gene_normalize=[None, "no"]),
        consequent("signal_file", 
            gene_normalize=["variance", "sum_of_squares"])),
    ]


def _item2seq(x):
    import operator
    if type(x) is type(""):
        return [x]
    if not operator.isSequenceType(x):
        x = [x]
    return x


def intersection(x, y):
    return list(set(x).intersection(y))


def find_modules(moduledb, out_data):
    # Return list of modules, sorted in decreasing order of the number
    # of attributes provided.
    matches = []  # list of (module, num attributes provided)
    for module in moduledb:
        if module.out_data.datatype != out_data.datatype:
            continue
        # Count the number of attributes this module provides.  If any
        # of the attributes conflict, then the number if 0.
        attributes_provided = 0
        for key, out_values in out_data.attributes.iteritems():
            out_values = _item2seq(out_values)
            if key not in module.out_data.attributes:
                continue
            module_values = _item2seq(module.out_data.attributes[key])
            if not intersection(out_values, module_values):
                attributes_provided = 0
                break
            attributes_provided += 1
        if not attributes_provided:
            continue
        x = module, attributes_provided
        matches.append(x)

    # Sort by decreasing number of attributes provided.
    schwartz = [(-x[-1], x) for x in matches]
    schwartz.sort()
    matches = [x[-1] for x in schwartz]
    modules = [x[0] for x in matches]
    return modules


def _backchain_data(module, out_data):
    # Return list of data that can be antecedents of module.  A module
    # can have multiple antecedents if its parameters can be different
    # values.
    import itertools

    datatype = module.in_data.datatype
    attributes = out_data.attributes.copy()
    for key in module.out_data.attributes:
        if key in attributes:
            del attributes[key]

    # Make all combinations of attributes.
    back_data = []
    keys = sorted(module.in_data.attributes)
    values_list = [_item2seq(module.in_data.attributes[x]) for x in keys]
    for values in itertools.product(values_list):
        for k, v in zip(keys, values):
            attributes[k] = v
        x = Data(datatype, attributes)
        back_data.append(x)
    return back_data
    

def find_pipelines(moduledb, in_data, out_data):
    # Do the backwards chaining using a depth-first search.  Keep
    # track of the stack as a list of (list of modules, list of data).
    # If the list of modules is length N, the list of data will be
    # length N+1.  The antecedent and consequent of modules[i] is
    # data[i] and data[i+1], respectively.
    pipelines = []
    
    stack = []
    stack.append(([], [out_data]))
    while stack:
        x = stack.pop()  # DFS: add to the end, pop from the end.
        modules, data_list = x
        # Backwards chain to the previous module that can create this
        # datatype and attributes.
        x = find_modules(moduledb, data_list[0])
        # If nothing can create this, then we're done.
        if not x:
            pipelines.append((modules, data_list))
            continue
        # Add each of the modules identified to the possible pipelines.
        for module in x:
            for x in _backchain_data(module, data_list[0]):
                back_data = x
                x = [module] + modules, [back_data] + data_list
                stack.append(x)

    # XXX check for concordance with in_data.
                
    return pipelines


def test_bie():
    pipelines = find_pipelines(
        all_modules, make_data("cel_files"),
        make_data("signal_file", preprocess="rma", quantile_norm="yes",
                  logged="yes"))
    for i, (modules, data) in enumerate(pipelines):
        print "PIPELINE %d" % (i+1)
        for i in range(len(modules)):
            print "STEP %d" % (i+1)
            print "Antecedent: %s" % data[i]
            print "Module: %s" % modules[i]
            print "Consequent: %s" % data[i+1]
            print 


if __name__ == '__main__':
    test_bie()
