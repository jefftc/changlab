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
    #cel_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id", platform="unknown"),
        consequent("cel_files", cel_version = "unknown")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform",platform="GPL"),
        consequent("cel_files", cel_version = "unknown")),
    Module(
        "extract_CEL_files",
        antecedent("cel_files", cel_version="unknown"),
        consequent("cel_files", cel_version = "cc_or_v3_4")),
    Module(
        "convert_CEL_to_v3_4",
        antecedent("cel_files", cel_version="cc_or_v3_4"),
        consequent("cel_files", cel_version = "v3_4")),
    Module(
        "preprocess_rma",
        antecedent("cel_files", cel_version="v3_4"),
        consequent("signal_raw",
                  logged="yes", preprocess="rma",
                   format='jeffs',missing=[None,"no"])),
    Module(
        "preprocess_mas5",
        antecedent("cel_files", cel_version="v3_4"),
        consequent("signal_raw", logged=[None,"no"], preprocess="mas5",format='jeffs',
                   missing=[None,"no"])),
    #-----------------------------------------------------------------------
    #agilent_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id", platform="unknown"),
        consequent("agilent_files", version = "unknown")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform",platform="GPL"),
        consequent("agilent_files", version = "unknown")),
    Module(
        "extract_agilent_files",
        antecedent("agilent_files", version="unknown"),
        consequent("agilent_files", version = "agilent")),
    Module(
        "preprocess_agilent",
        antecedent("agilent_files", version="agilent"),
        consequent("signal_raw", format = "tdf",logged=[None,"no"],preprocess='agilent',missing='unknown')),
    #-----------------------------------------------------------------------
    #idat_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id", platform="unknown"),
        consequent("idat_files", version = "unknown")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform",platform="GPL"),
        consequent("idat_files", version = "unknown")),
    Module(
        "extract_illumina_idat_files",
        antecedent("idat_files", version="unknown"),
        consequent("idat_files", version = "illumina")),
    Module(
        "preprocess_illumina",
        antecedent("idat_files", version="illumina"),
        consequent("illu_folder", format = "gct",logged=[None,"no"],preprocess='illumina',missing='unknown')),
    Module(
        "get_illumina_signal",
        antecedent("illu_folder", preprocess='illumina'),
        consequent("signal_raw", preprocess='illumina',logged=[None,"no"],missing='unknown')),

    #-----------------------------------------------------------------------
    #gpr_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id", platform="unknown"),
        consequent("gpr_files", version = "unknown")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform",platform="GPL"),
        consequent("gpr_files", version = "unknown")),
    Module(
        "extract_gpr_files",
        antecedent("gpr_files", version="unknown"),
        consequent("gpr_files", version = "gpr")),
    Module(
        "normalize_with_loess",
        antecedent("gpr_files", version="gpr"),
        consequent("signal_raw", format = "tdf",logged=[None,"no"],preprocess='loess',missing='unknown')),
    #-----------------------------------------------------------------------
    Module(
        "filter_genes_by_missing_values",
        antecedent("signal_raw", format='tdf',logged="yes",missing=["yes","unknown"],filter=[None,"no"]),
        consequent("signal_raw", format='tdf',logged="yes",missing=["yes","unknown"],filter=20)),###integer value
    Module(
        "fill_missing_with_median",
        antecedent("signal_raw", format='tdf',logged="yes",missing=["yes","unknown"]),
        consequent("signal_raw", format='tdf',logged="yes",missing="median")),
    Module(
        "fill_missing_with_zeros",
        antecedent("signal_raw", format='tdf',logged="yes",missing=["yes","unknown"]),
        consequent("signal_raw", format='tdf',logged="yes",missing="zero")),
    Module(
        "convert_signal_to_tdf",
        antecedent("signal_raw", format=['pcl','res','gct','jeffs','unknown','xls']),
        consequent("signal_raw", format='tdf')),
    Module(
        "log_signal",
        antecedent("signal_raw", logged=[None, "no"],format='tdf'),
        consequent("signal_raw", logged="yes",format='tdf')),
    Module(
        "filter_and_threshold_genes",
        antecedent("signal_raw", logged=[None, "no"],format='tdf',predataset=[None,'no']),
        consequent("signal_raw", logged=[None, "no"],format='tdf',predataset="yes")),
    Module(#require a rename_list_file
        "relabel_samples",
        antecedent("signal_raw", logged="yes",format='tdf',missing=[None,"no","median","zero"],
                   rename_sample = [None,"no"]),
        consequent("signal_raw", logged="yes",format='tdf',missing=[None,"no","median","zero"],
                   rename_sample = "yes")),
    #------------------------------------------------------------------
    Module(
        "quantile_norm",
        antecedent("signal_raw", quantile_norm=[None, "no"],gene_center=[None,"no"],
                   gene_normalize=[None, "no"],format='tdf',missing=[None,"no","median","zero"]),
        consequent("signal_file", quantile_norm="yes",gene_center=[None,"no"],
                   gene_normalize=[None, "no"],format='tdf',missing=[None,"no","median","zero"])),
    #------------------------------------------------------------------
    Module(
        "gene_center",
        antecedent("signal_file", gene_center=[None, "no"],logged="yes",
                   gene_normalize=[None, "no"],format='tdf',missing=[None,"no","median","zero"]),
        consequent("signal_file", gene_center=["mean", "median"],logged="yes",
                   gene_normalize=[None, "no"],format='tdf',missing=[None,"no","median","zero"])),
    Module(
        "gene_normalize",
        antecedent("signal_file", gene_normalize=[None, "no"],logged="yes",format='tdf',
                   missing=[None,"no","median","zero"]),
        consequent("signal_file", 
            gene_normalize=["variance", "sum_of_squares"],logged="yes",format='tdf',
                   missing=[None,"no","median","zero"]))]
                 


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
    for k, v in zip(keys, values_list):
        attributes[k] = v
    x = Data(datatype, attributes)
    back_data.append(x)
##    for values in itertools.product(values_list):        
##        for k, v in zip(keys, values):
##            attributes[k] = v
##        x = Data(datatype, attributes)
##        back_data.append(x)

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
        all_modules, make_data("gse_id",platform='unknown'),
        make_data("signal_file", preprocess='rma',logged='yes',filter='no',missing=None,
                  predataset='no',rename_sample='no',format='tdf',gene_center='mean',
                  quantile_norm='yes'))
    for i, (modules, data) in enumerate(pipelines):
        print "PIPELINE %d" % (i+1)
        for i in range(len(modules)):
            #print "STEP %d" % (i+1)
            #print "Antecedent: %s" % data[i]
            print "Module: %s" % modules[i].name
            #print "Consequent: %s" % data[i+1]
            print 


if __name__ == '__main__':
    test_bie()
