"""

Glossary:
goal
rule            another rule is proposition
consequent      second half of a rule
antecedent      first half of a rule
fact

Modules "convert" one data object to another.


Functions:
make_data        Make a Data object.
antecedent       Synonym for make_data.
consequent       Synonym for make_data.

backchain
prune_network_by_start

test_bie


Classes:
Data
Module
Network

"""

# _backchain_to_modules       Given a consequent, list compatible modules.
# _backchain_to_antecedents   Given a consequent and module, list antecedents.
# _calc_compat_cons_data      Calc the compatability of a consequent and data.
#
# _print_nothing
# _print_string
# _print_line
# _print_network
# _plot_network_gv
# _pretty_attributes
#
# _find_node
# _is_item
# _item2seq
# _intersection


class Data:
    # Members:
    # datatype     string.  Name of data.
    # attributes   Dict of key -> value or list of values.
    def __init__(self, datatype, attributes):
        self.datatype = datatype
        self.attributes = attributes.copy()
    def copy(self):
        return Data(self.datatype, self.attributes.copy())
    def __cmp__(self, other):
        x1 = [self.datatype, self.attributes]
        x2 = [other.datatype, other.attributes]
        return cmp(x1, x2)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            repr(self.datatype),
            _pretty_attributes(self.attributes),
            ]
        x = [x for x in x if x]
        return "Data(%s)" % ", ".join(x)

class Module:
    def __init__(self, name, ante_data, cons_data):
        self.name = name
        self.ante_data = ante_data.copy()
        self.cons_data = cons_data.copy()
    def __cmp__(self, other):
        x1 = [self.name, self.ante_data, self.cons_data]
        x2 = [other.name, other.ante_data, other.cons_data]
        return cmp(x1, x2)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            repr(self.name),
            repr(self.ante_data),
            repr(self.cons_data),
            ]
        x = "Module(%s)" % ", ".join(x)
        return x


# nodes are Data or Module objects.  Data transition to Modules, and
# Modules transition to Data.
class Network:
    def __init__(self, nodes, transitions):
        self.nodes = nodes[:]
        self.transitions = transitions.copy()


def make_data(datatype, **params):
    return Data(datatype, params)


antecedent = make_data
consequent = make_data


def backchain(moduledb, out_data):
    # Return a Network object.
    nodes = []        # list of Data or Module objects.
    transitions = {}  # list of index -> list of indexes

    MAX_NETWORK_SIZE = 128
    nodes.append(out_data)
    stack = [0]
    while stack:
        assert len(nodes) < MAX_NETWORK_SIZE, "network too large"
        node_id = stack.pop()
        assert node_id < len(nodes)
        node = nodes[node_id]
        node_type = node.__class__.__name__

        if node_type == "Data":
            # Backwards chain to the previous module that can create this
            # datatype and attributes.
            modules = _backchain_to_modules(moduledb, node)
            for m in modules:
                # This module should not already occur.  If it does,
                # then that means there's a cycle?
                nodes.append(m)
                module_id = len(nodes) - 1
                stack.append(module_id)
                if module_id not in transitions:
                    transitions[module_id] = []
                transitions[module_id].append(node_id)
        elif node_type == "Module":
            # Each Module should point to exactly one Data.
            x = transitions[node_id]
            assert len(x) == 1
            consequent_id = x[0]

            for x in _backchain_to_antecedents(node, nodes[consequent_id]):
                back_data = x
                #print back_data
                #print node
                #print nodes[consequent_id]
                back_id = _find_node(nodes, back_data)
                if back_id == -1:
                    nodes.append(back_data)
                    back_id = len(nodes)-1
                    if back_id not in stack:
                        stack.append(back_id)
                if back_id not in transitions:
                    transitions[back_id] = []
                if node_id not in transitions[back_id]:
                    transitions[back_id].append(node_id)
        else:
            raise AssertionError, "Unknown node type: %s" % node_type
    return Network(nodes, transitions)


def _is_compatible_with_start(data, start_data):
    if data.datatype != start_data.datatype:
        return False

    data_attr = data.attributes
    start_attr = start_data.attributes
    
    x = data_attr.keys() + start_attr.keys()
    all_attributes = sorted({}.fromkeys(x))

    # Case 1.  Attribute in data but not start_data.
    #          User does not care.  Is still compatible.
    # Case 2.  Attribute in start_data but not data.
    #          Incompatible.
    # Case 3.  Attribute values are items in both.
    #          Compatible if the items are equal.
    # Case 4.  Attribute values is item in data and sequence in start.
    #          User has a specific item.  Start can accept multiple items.
    #          Compatible if item is in start sequence.
    # Case 5.  Attribute values is sequence in data and item in start.
    #          Incompatible.  Value in user is ambiguous, but start
    #          requires a specific value.
    # Case 6.  Attribute values is sequence in data and start.
    #          Compatible only if sequences are the same.

    compatible = True
    for key in all_attributes:
        key_in_data = key in data_attr
        key_in_start = key in start_attr
        
        data_is_item = key_in_data and _is_item(data_attr[key])
        start_is_item = key_in_start and _is_item(start_attr[key])
        data_is_list = key_in_data and not _is_item(data_attr[key])
        start_is_list = key_in_start and not _is_item(start_attr[key])
        
        # Case 1.
        if key_in_data and not key_in_start:
            pass
        # Case 2.
        elif not key_in_data and key_in_start:
            compatible = False
        # Case 3.
        elif data_is_item and start_is_item:
            if start_attr[key] != data_attr[key]:
                compatible = False
        # Case 4.
        elif data_is_item and start_is_list:
            if data_attr[key] not in start_attr[key]:
                compatible = False
        # Case 5.
        elif data_is_list and start_is_item:
            compatible = False
        # Case 6.
        elif data_is_list and start_is_list:
            if sorted(data_attr[key]) != sorted(start_attr[key]):
                compatible = False
        else:
            raise NotImplementedError
    return compatible
            
        
def prune_network_by_start(network, start_data):
    # Look for the nodes that are compatible with start_data.
    node_ids = []  # list of node_ids.
    for i, n in enumerate(network.nodes):
        if n.__class__.__name__ != "Data":
            continue
        if _is_compatible_with_start(n, start_data):
            node_ids.append(i)

    # For each of these node_ids, do forward chaining to find all
    # nodes that these ones can connect to.
    start_ids = {}
    stack = node_ids[:]
    while stack:
        node_id = stack.pop(0)
        start_ids[node_id] = 1
        for next_id in network.transitions.get(node_id, []):
            if next_id not in start_ids:
                stack.append(next_id)
    start_ids = sorted(start_ids)
    
    new_nodes = [network.nodes[i] for i in start_ids]
    new_transitions = {}
    for node_id, next_ids in network.transitions.iteritems():
        if node_id not in start_ids:
            continue
        for next_id in next_ids:
            if next_id not in start_ids:
                continue
            new_node_id = start_ids.index(node_id)
            new_next_id = start_ids.index(next_id)
            if new_node_id not in new_transitions:
                new_transitions[new_node_id] = []
            new_transitions[new_node_id].append(new_next_id)

    return Network(new_nodes, new_transitions)


def test_bie():
    in_data = make_data("gse_id", platform='unknown')
    #in_data = make_data("gse_id")
    #in_data = make_data("signal_file", logged="yes")
    #out_data = make_data(
    #    "signal_file", format='tdf',preprocess='rma', logged='yes',
    #    filter='no', missing=None, predataset='no', rename_sample='no', 
    #    gene_center='mean', quantile_norm='yes')
    out_data = make_data(
        "signal_file", format='tdf',preprocess='rma',logged='yes',
        filter='no', missing='zero', predataset='no', rename_sample='no',
         quantile_norm='yes',dwd_norm='no',gene_order='class_neighbors',
        shiftscale_norm='no',combat_norm='no',bfrm_norm='yes',gene_center='mean',
        gene_normalize='variance',group_fc='int',
        num_features='int',annotate="yes",unique_genes='average_genes',platform='str',
        duplicate_probe='high_var_probe')
    #out_data = make_data(
    #    "signal_file", filter='no', format='tdf', logged='yes',
    #    preprocess='rma')
    #out_data = make_data(
    #    "signal_file", preprocess='rma', logged='yes')
    
    network = backchain(all_modules, out_data)
    network = prune_network_by_start(network, in_data)
    
    print "DONE FINDING PIPELINES"
    _print_network(network)
    _plot_network_gv("out.png", network)
    #_print_many_pipelines(pipelines, verbose=True)
    #_print_many_pipelines(pipelines, verbose=False)
    #for i, (modules, data) in enumerate(pipelines):
    #    print "PIPELINE %d" % (i+1)
    #    for i in range(len(modules)):
    #        print "STEP %d" % (i+1)
    #        print "Antecedent: %s" % data[i]
    #        print "Module: %s" % modules[i]
    #        print "Consequent: %s" % data[i+1]
    #        print 


def _backchain_to_modules(moduledb, data):
    # Return list of modules that that can generate a consequent that
    # is compatible with data.  The modules will be sorted in
    # decreasing order of the number of attributes that are
    # compatible.

    matches = []  # list of (module, num compatible attributes)
    for module in moduledb:
        num_attributes = _calc_compat_cons_data(module, data)
        if not num_attributes:
            continue
        x = module, num_attributes
        matches.append(x)

    # Sort by decreasing number of attributes provided.
    matches.sort(key=lambda x: -x[-1])
    modules = [x[0] for x in matches]
    return modules


def _backchain_to_antecedents(module, data):
    # Return list of data that can be antecedents of module.  A module
    # can have multiple antecedents if its parameters can be different
    # values.
    import itertools

    # The datatype is from the antecedent of the module.
    datatype = module.ante_data.datatype

    # Back chain the attributes.  Possibilities:
    # Case 1.  Attribute not in antecedent or consequent.
    #          Keep attribute in data unchanged.
    # Case 2.  Attribute in antecedent but not consequent.
    #          E.g. cel_version="v3_4"
    #          Set the value of the attribute from the antecedent.
    # Case 3.  Attribute in consequent but not antecedent.
    #          E.g. preprocess="rma"
    #          Ignore the attribute.
    # Case 4.  Value of attribute in antecedent and consequent are
    #          different.  E.g. logged="no" -> logged="yes"
    #          Set the value of the attribute from the antecedent.
    # Case 5.  Value of attribute in antecedent and consequent are the
    #          same.  E.g. format="tdf" -> format="tdf"
    #          Set the value of the attribute from the antecedent.
    # Case 6.  Value of attribute in antecedent is list, while that
    #          from consequent is a single object and the values in
    #          the consequent is a subset of those in the antecedent.
    #          E.g. format=["cdt", "tdf"] -> format="tdf".
    #          Set the value of the attribute from the antecedent.
    # Case 7.  Value of attribute in antecedent is list, while that
    #          from consequent is a single object and the values in
    #          the consequent are different from those in the
    #          antecedent.
    #          E.g. quantile_norm=[None, "no"] -> quantile_norm="yes".
    #          Set the value of the attribute from the antecedent.
    # Case 8.  Value of attribute in antecedent and consequent are
    #          both lists.  E.g. gene_normalize=[None, "no"] ->
    #          gene_normalize=["variance", "sum_of_squares"]
    #          Set the value of the attribute from the antecedent.
    # Anything else is not yet defined or implemented.

    data_attr = data.attributes
    ante_attr = module.ante_data.attributes
    cons_attr = module.cons_data.attributes

    x = data_attr.keys() + ante_attr.keys() + cons_attr.keys()
    all_attributes = sorted({}.fromkeys(x))

    attributes = {}
    for key in all_attributes:
        key_in_ante = key in ante_attr
        key_in_cons = key in cons_attr
        ante_is_item = key_in_ante and _is_item(ante_attr[key])
        cons_is_item = key_in_cons and _is_item(cons_attr[key])
        ante_is_list = key_in_ante and not _is_item(ante_attr[key])
        cons_is_list = key_in_cons and not _is_item(cons_attr[key])
        
        # Case 1.
        if not key_in_ante and not key_in_cons:
            attributes[key] = data_attr[key]
        # Case 2.
        elif key_in_ante and not key_in_cons:
            attributes[key] = ante_attr[key]
        # Case 3.
        elif not key_in_ante and key_in_cons:
            pass
        # Case 4.
        elif ante_is_item and cons_is_item and ante_attr[key] !=cons_attr[key]:
            attributes[key] = ante_attr[key]
        # Case 5.
        elif ante_is_item and cons_is_item and ante_attr[key] ==cons_attr[key]:
            attributes[key] = ante_attr[key]
        # Case 6.
        elif ante_is_list and cons_is_item and \
                 cons_attr[key] in ante_attr[key]:
            attributes[key] = ante_attr[key]
        # Case 7.
        elif ante_is_list and cons_is_item and \
                 cons_attr[key] not in ante_attr[key]:
            attributes[key] = ante_attr[key]
        # Case 8.
        elif ante_is_list and cons_is_list:
            attributes[key] = ante_attr[key]
        else:
            print key
            print ante_attr.get(key, "MISSING")
            print cons_attr.get(key, "MISSING")
            print data_attr.get(key, "MISSING")
            raise NotImplementedError

    x = Data(datatype, attributes)
    return [x]


def _calc_compat_cons_data(module, data):
    # Return the number of attributes in the consequent of the module
    # that is compatible with the attributes of data.  0 means that
    # module can not produce data.

    p = _print_nothing
    #p = _print_string
    
    p("Evaluating module %s." % module.name)

    # If this module doesn't produce the right data type, then
    # ignore it.
    if module.cons_data.datatype != data.datatype:
        p("Not right data type [%s %s]." % (
            module.cons_data.datatype, data.datatype))
        return 0

    ante_attributes = module.ante_data.attributes
    cons_attributes = module.cons_data.attributes

    # Count the number of attributes in the consequent that is
    # compatible with the attributes for data.
    num_attributes = 0
    for key in data.attributes:
        p("Evaluating attribute %s." % key)
        data_values = _item2seq(data.attributes[key])
        ante_values = _item2seq(ante_attributes.get(key, []))
        cons_values = _item2seq(cons_attributes.get(key, []))

        # Checks:
        # 1.  If module does not produce attribute, then ignore it.
        #     This goes first.  If module is indifferent to this
        #     attribute, then don't do anything.
        # 2.  If module produces an attribute that conflicts the one
        #     with the data, then ignore the whole module.
        # 3.  If this module does not change this attribute, then
        #     ignore it.  This check has to go after #2 to make sure
        #     conflicts are not ignored.

        # If this module does not produce this attribute, then
        # ignore it.
        if not cons_values:
            p("Module does not produce this attribute.")
            continue

        # If this module generates an attribute whose value conflicts
        # with one in the data, then ignore the module.
        if not _intersection(cons_values, data_values):
            p("Module produces incompatible values for this attribute.")
            num_attributes = 0
            break

        # If this module does not change this attribute, then
        # ignore it.
        if sorted(ante_values) == sorted(cons_values):
            p("Module does not change this attribute.")
            continue

        p("Attribute is compatible.")
        num_attributes += 1
    if not num_attributes:
        p("No attributes compatible.")
    else:
        p("%d attributes compatible." % num_attributes)
    return num_attributes


def _print_nothing(s):
    pass


def _print_string(s):
    print s


def _print_line(line, prefix0, prefixn, width):
    lines = []
    p = prefix0
    while line:
        x = line[:(width-len(p))]
        line = line[len(x):]
        lines.append(x)
        p = prefixn
        
    for i in range(len(lines)):
        p = prefixn
        if i == 0:
            p = prefix0
        print "%s%s" % (p, lines[i])


## def _print_pipeline_verbose(modules, data_list):
##     assert len(modules)+1 == len(data_list)

##     line_width = 72
##     p_space0 = " " * 4
##     p_spacen = " " * 6
##     for j in range(len(modules)):
##         p_step = "%d.  " % (j+1)
##         p_space0 = "*" + " " * (len(p_step)-1)
##         p_spacen = " " * (len(p_step)+2)
##         _print_line(str(data_list[j]), p_space0, p_spacen, line_width)
##         _print_line(str(modules[j]), p_step, p_spacen, line_width)
##     _print_line(str(data_list[-1]), p_space0, p_spacen, line_width)


## def _print_pipeline_succinct(modules, data_list):
##     assert len(modules)+1 == len(data_list)

##     for j in range(len(modules)):
##         print "%d.  %s" % (j+1, modules[j].name)


## def _print_pipeline(modules, data_list, verbose=True):
##     if verbose:
##         _print_pipeline_verbose(modules, data_list)
##     else:
##         _print_pipeline_succinct(modules, data_list)


## def _print_many_pipelines(stack, verbose=True):
##     # List of (list of modules, list of data).
##     print "MANY PIPELINES"
##     for i in range(len(stack)):
##         print "PIPELINE %d of %d" % (i+1, len(stack))
##         modules, data_list = stack[i]
##         _print_pipeline(modules, data_list, verbose=verbose)
##         print


def _print_network(network):
    line_width = 72
    for i in range(len(network.nodes)):
        p_step = "%d.  " % i
        p_space = " " * (len(p_step)+2)
        _print_line(str(network.nodes[i]), p_step, p_space, line_width)
    print

    for i in sorted(network.transitions):
        x = [i, "->"] + network.transitions[i]
        print "\t".join(map(str, x))
    

def _plot_network_gv(filename, network):
    from genomicode import graphviz
    
    gv_nodes = []
    gv_edges = []
    gv_node2attr = {}

    id2name = {}
    for node_id, n in enumerate(network.nodes):
        node2attr = {}
        node2attr["style"] = "filled"
        if n.__class__.__name__ == "Data":
            x = n.datatype
            node2attr["shape"] = "note"
            node2attr["fillcolor"] = "#EEEEEE"
        else:
            x = n.name
            node2attr["shape"] = "box"
            node2attr["fillcolor"] = "#80E0AA"
            
        node_name = "%s [%d]" % (x, node_id)
        id2name[node_id] = node_name
        gv_nodes.append(node_name)
        gv_node2attr[node_name] = node2attr

    for node_id, next_ids in network.transitions.iteritems():
        for nid in next_ids:
            x1 = id2name[node_id]
            x2 = id2name[nid]
            gv_edges.append((x1, x2))

    graphviz.draw(
        filename, gv_nodes, gv_edges, node2attributes=gv_node2attr, prog="dot",
        directed=True)

    
def _pretty_attributes(attributes):
    import re

    if not attributes:
        return ""

    # Separate the proper python variables from other attributes.
    proper = []
    improper = []
    for key in sorted(attributes):
        if re.match(r"^[a-z_][a-z0-9_]*$", key, re.I):
            proper.append(key)
        else:
            improper.append(key)

    # For proper python variables, format as <key>=<value>.
    fmt_proper = []
    for key in proper:
        x = "%s=%r" % (key, attributes[key])
        fmt_proper.append(x)
    fmt_proper = ", ".join(fmt_proper)

    # For improper python variables, format as a dict.
    x = {}
    for key in improper:
        x[key] = attributes[key]
    fmt_improper = ""
    if x:
        fmt_improper = repr(x)

    if not fmt_improper:
        return fmt_proper
    return "%s, %s" % (fmt_proper, fmt_improper)


def _find_node(nodes, node):
    for i, n in enumerate(nodes):
        if n.__class__.__name__ != node.__class__.__name__:
            continue
        if node == n:
            return i
    return -1


def _is_item(x):
    import operator
    if type(x) is type(""):
        return True
    if not operator.isSequenceType(x):
        return True
    return False


def _item2seq(x):
    import operator
    if _is_item(x):
        return [x]
    return x


def _intersection(x, y):
    return list(set(x).intersection(y))


all_modules = [
    # cel_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id", platform="unknown"),
        consequent("cel_files", cel_version="unknown")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform", platform="GPL"),
        consequent("cel_files", cel_version="unknown")),
    Module(
        "extract_CEL_files",
        antecedent("cel_files", cel_version="unknown"),
        consequent("cel_files", cel_version="cc_or_v3_4")),
    Module(
        "convert_CEL_to_v3_4",
        antecedent("cel_files", cel_version="cc_or_v3_4"),
        consequent("cel_files", cel_version="v3_4")),
    Module(
        "preprocess_rma",
        antecedent("cel_files", cel_version="v3_4"),
        consequent(
            "signal_file", logged="yes", preprocess="rma", format='jeffs',
            missing="no")),
    Module(
        "preprocess_mas5",
        antecedent("cel_files", cel_version="v3_4"),
        consequent(
            "signal_file", logged="no", preprocess="mas5",
            format="jeffs", missing="no")),
    
    #-----------------------------------------------------------------------
    # agilent_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id", platform="unknown"),
        consequent("agilent_files", version="unknown")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform", platform="GPL"),
        consequent("agilent_files", version="unknown")),
    Module(
        "extract_agilent_files",
        antecedent("agilent_files", version="unknown"),
        consequent("agilent_files", version="agilent")),
    Module(
        "preprocess_agilent",
        antecedent("agilent_files", version="agilent"),
        consequent(
            "signal_file", format="tdf", logged="no",
            preprocess="agilent", missing="unknown")),
    
    #-----------------------------------------------------------------------
    # idat_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id", platform="unknown"),
        consequent("idat_files", version="unknown")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform", platform="GPL"),
        consequent("idat_files", version="unknown")),
    Module(
        "extract_illumina_idat_files",
        antecedent("idat_files", version="unknown"),
        consequent("idat_files", version="illumina")),
    Module(
        "preprocess_illumina",
        antecedent("idat_files", version="illumina"),
        consequent(
            "illu_folder", format="gct", logged="no",
            preprocess='illumina', missing='unknown',
            ill_manifest=['HumanHT-12_V3_0_R2_11283641_A.txt',
                    'HumanHT-12_V4_0_R2_15002873_B.txt',
                    'HumanHT-12_V3_0_R3_11283641_A.txt',
                    'HumanHT-12_V4_0_R1_15002873_B.txt',
                    'HumanMI_V1_R2_XS0000122-MAP.txt',
                    'HumanMI_V2_R0_XS0000124-MAP.txt',
                    'HumanRef-8_V2_0_R4_11223162_A.txt',
                    'HumanRef-8_V3_0_R1_11282963_A_WGDASL.txt',
                    'HumanRef-8_V3_0_R2_11282963_A.txt',
                    'HumanRef-8_V3_0_R3_11282963_A.txt',
                    'HumanWG-6_V2_0_R4_11223189_A.txt',
                    'HumanWG-6_V3_0_R2_11282955_A.txt',
                    'HumanWG-6_V3_0_R3_11282955_A.txt',
                    'MouseMI_V1_R2_XS0000127-MAP.txt',
                    'MouseMI_V2_R0_XS0000129-MAP.txt',
                    'MouseRef-8_V1_1_R4_11234312_A.txt',
                    'MouseRef-8_V2_0_R2_11278551_A.txt',
                    'MouseRef-8_V2_0_R3_11278551_A.txt',
                    'MouseWG-6_V1_1_R4_11234304_A.txt',
                    'MouseWG-6_V2_0_R2_11278593_A.txt',
                    'MouseWG-6_V2_0_R3_11278593_A.txt',
                    'RatRef-12_V1_0_R5_11222119_A.txt'],
            ill_chip=['ilmn_HumanHT_12_V3_0_R3_11283641_A.chip',
                'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
                'ilmn_HumanRef_8_V2_0_R4_11223162_A.chip',
                'ilmn_HumanReF_8_V3_0_R1_11282963_A_WGDASL.chip',
                'ilmn_HumanRef_8_V3_0_R3_11282963_A.chip',
                'ilmn_HumanWG_6_V2_0_R4_11223189_A.chip',
                'ilmn_HumanWG_6_V3_0_R3_11282955_A.chip',
                'ilmn_MouseRef_8_V1_1_R4_11234312_A.chip',
                'ilmn_MouseRef_8_V2_0_R3_11278551_A.chip',
                'ilmn_MouseWG_6_V1_1_R4_11234304_A.chip',
                'ilmn_MouseWG_6_V2_0_R3_11278593_A.chip',
                'ilmn_RatRef_12_V1_0_R5_11222119_A.chip'],
            ill_bg_mode=['false','true'],
            ill_coll_mode=['none','max','median'],
            ill_clm='string',
            ill_custom_chip='string',
            ill_custom_manifest='string')),
    Module(
        "get_illumina_signal",
        antecedent("illu_folder", preprocess='illumina'),
        consequent(
            "signal_file", preprocess='illumina', logged="no",
            missing='unknown',format='gct')),
   Module(
        "get_illumina_control",
        antecedent("illu_folder", preprocess='illumina'),
        consequent(
            "control_file", preprocess='illumina', logged="no",
            missing='unknown',format='gct')),
    #-----------------------------------------------------------------------
    # gpr_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id", platform="unknown"),
        consequent("gpr_files", version="unknown")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform", platform="GPL"),
        consequent("gpr_files", version="unknown")),
    Module(
        "extract_gpr_files",
        antecedent("gpr_files", version="unknown"),
        consequent("gpr_files", version="gpr")),
    Module(
        "normalize_with_loess",
        antecedent("gpr_files", version="gpr"),
        consequent(
            "signal_file", format="tdf", logged="no",
            preprocess="loess", missing="unknown")),
    
    #-----------------------------------------------------------------------
    Module(
        "convert_signal_to_tdf",
        antecedent(
            "signal_file",
            format=['pcl', 'res', 'gct', 'jeffs', 'unknown', 'xls'],
            quantile_norm=[None,'no'],bfrm_norm=[None,'no'],combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm=[None,'no'],gene_center=[None,'no'],gene_normalize=[None,'no'],filter=[None,'no'],
            predataset='no',rename_sample='no',group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent("signal_file", format='tdf',quantile_norm=[None,'no'],
            bfrm_norm=[None,'no'],combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm=[None,'no'],gene_center=[None,'no'],gene_normalize=[None,'no'],filter=[None,'no'],
            predataset='no',rename_sample='no',group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "log_signal",
        antecedent("signal_file", logged=[None, "no",'unknown'], format='tdf',
            quantile_norm=[None,'no'],bfrm_norm=[None,'no'],combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm=[None,'no'],gene_center=[None,'no'],gene_normalize=[None,'no'],
            predataset='no',rename_sample='no',filter=[None,'no'],group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent("signal_file", logged="yes", format='tdf',
            quantile_norm=[None,'no'],bfrm_norm=[None,'no'],combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm=[None,'no'],gene_center=[None,'no'],gene_normalize=[None,'no'],
            predataset='no',rename_sample='no',filter=[None,'no'],group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "filter_genes_by_missing_values",
        antecedent(
            "signal_file", format='tdf', logged="yes",missing=["yes", "unknown"], filter=[None, "no"],
            quantile_norm=[None,'no'],bfrm_norm=[None,'no'],combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm=[None,'no'],gene_center=[None,'no'],gene_normalize=[None,'no'],
            predataset='no',rename_sample='no',group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", format='tdf', logged="yes", missing=["yes", "unknown"], filter=20,
            quantile_norm=[None,'no'],bfrm_norm=[None,'no'],combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm=[None,'no'],gene_center=[None,'no'],gene_normalize=[None,'no'],
            predataset='no',rename_sample='no',group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')), ###integer value
    Module(
        "fill_missing_with_median",
        antecedent(
            "signal_file", format='tdf', logged="yes",
            missing=["yes", "unknown"],quantile_norm=[None,'no'],bfrm_norm=[None,'no'],
            combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm=[None,'no'],gene_center=[None,'no'],gene_normalize=[None,'no'],
            predataset='no',rename_sample='no',group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", format='tdf', logged="yes", missing="median",
            quantile_norm=[None,'no'],bfrm_norm=[None,'no'],combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm=[None,'no'],gene_center=[None,'no'],gene_normalize=[None,'no'],
            predataset='no',rename_sample='no',group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "fill_missing_with_zeros",
        antecedent(
            "signal_file", format='tdf', logged="yes",
            missing=["yes", "unknown"],quantile_norm=[None,'no'],bfrm_norm=[None,'no'],
            combat_norm=[None,'no'],shiftscale_norm=[None,'no'],dwd_norm=[None,'no'],gene_center=[None,'no'],
            gene_normalize=[None,'no'],predataset='no',rename_sample='no',group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent("signal_file", format='tdf', logged="yes", missing="zero",
            quantile_norm=[None,'no'],bfrm_norm=[None,'no'],combat_norm=[None,'no'],shiftscale_norm=[None,'no'],
            dwd_norm='no',gene_center=[None,'no'],gene_normalize=[None,'no'],
            predataset='no',rename_sample='no',group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "filter_and_threshold_genes",
        antecedent(
            "signal_file", logged=[None, "no"], format='tdf',predataset=[None, 'no']),
        consequent(
            "signal_file", logged=[None, "no"], format='tdf',predataset="yes")),
    Module( #require a rename_list_file
        "relabel_samples",
        antecedent(
            "signal_file", logged="yes", format='tdf',missing=[None, "no", "median", "zero"],
            rename_sample=[None, "no"],quantile_norm=[None,'no'],bfrm_norm=[None,'no'],
            combat_norm=[None,'no'],shiftscale_norm=[None,'no'],dwd_norm=[None,'no'],
            gene_center=[None,'no'],gene_normalize=[None,'no'],group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", logged="yes", format='tdf',missing=[None, "no", "median", "zero"],
            rename_sample="yes",quantile_norm=[None,'no'],bfrm_norm=[None,'no'],combat_norm=[None,'no'],
            shiftscale_norm=[None,'no'],dwd_norm=[None,'no'],gene_center=[None,'no'],
            gene_normalize=[None,'no'],group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    #------------------------------------------------------------------
    Module(
        "normalize_samples_with_quantile",
        antecedent(
            "signal_file", quantile_norm=[None, "no"], format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", quantile_norm="yes", format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "normalize_samples_with_combat",#requir class label file
        antecedent(
            "signal_file", combat_norm=[None, "no"], format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", combat_norm="yes", format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "normalize_samples_with_dwd",#requir class label file
        antecedent(
            "signal_file", dwd_norm=[None, "no"], format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", dwd_norm="yes", format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "normalize_samples_with_bfrm",
        antecedent(
            "signal_file", bfrm_norm=[None, "no"], format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", bfrm_norm="yes", format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "normalize_samples_with_shiftscale",#require class label file
        antecedent(
            "signal_file", shiftscale_norm=[None, "no"], format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", shiftscale_norm="yes", format='tdf',
            logged="yes", missing=[None, "no", "median", "zero"],
            gene_center=[None,'no'],gene_normalize=[None,'no'],
            group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    #------------------------------------------------------------------
##    Module(   #may cause loop in the network 
##        "convert_signal_to_pcl",
##        antecedent(
##            "signal_file", logged="yes",
##            format='tdf', missing=[None, "no", "median", "zero"]),
##        consequent(
##            "signal_file",logged="yes", format='pcl',
##            missing=[None, "no", "median", "zero"])),
    Module(
        "gene_center",
        antecedent(
            "signal_file", gene_center=[None, "no"], logged="yes",
            gene_normalize=[None, "no"], format='tdf',
            missing=[None, "no", "median", "zeros"],group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", gene_center=["mean", "median"], logged="yes",
            gene_normalize=[None, "no"], format='tdf',
            missing=[None, "no", "median", "zero"],group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "gene_normalize",
        antecedent(
            "signal_file", gene_normalize=[None, "no"], logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],group_fc=[None,'no'],
            gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", gene_normalize=["variance", "sum_of_squares"],
            logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],group_fc=[None,'no'],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    #------------------------------------------------------------------
    Module(
        "filter_genes_by_fold_change_across_classes", #require class_label_file
        antecedent(
            "signal_file", group_fc=[None, "no"], logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", group_fc="int",logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "rank_genes_by_sample_ttest", #require class_label_file,generate gene_list_file and need reorder_genes
        antecedent(
            "signal_file", logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            gene_order=[None,'no'],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            gene_order=["t_test_p","t_test_fdr"],
            num_features='all',annotate=[None,'no'],unique_genes=[None,'no'],
            platform='unknown_platform',duplicate_probe='yes')),
    Module(
        "rank_genes_by_class_neighbors", #require class_label_file,generate gene_list_file and need reorder_genes
        antecedent(
            "signal_file", logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            platform='unknown_platform',gene_order=[None,'no'],
            unique_genes=[None,'no'],num_features='all',annotate=[None,'no'],duplicate_probe='yes'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            platform='unknown_platform',gene_order='class_neighbors',
            unique_genes=[None,'no'],num_features='all',annotate=[None,'no'],duplicate_probe='yes')),
    Module(
        "reorder_genes", #require gene_list_file
        antecedent(
            "signal_file",logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            platform='unknown_platform',gene_order=[None,'no'],
            unique_genes=[None,'no'],num_features='all',annotate=[None,'no'],duplicate_probe='yes'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            platform='unknown_platform',gene_order=['gene_list'],
            unique_genes=[None,'no'],num_features='all',annotate=[None,'no'],duplicate_probe='yes')),
    Module(
        'annotate_probes',
        antecedent(
            "signal_file",logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            platform='unknown_platform',unique_genes=[None,'no'],
            num_features='all',annotate=[None,"no"],duplicate_probe='yes'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            platform='unknown_platform',unique_genes=[None,'no'],
            num_features="all",annotate="yes",duplicate_probe='yes')),
    Module(
        'remove_duplicate_genes',
        antecedent(
            "signal_file",logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            platform='unknown_platform',
            unique_genes=[None,'no'],num_features='all',annotate='yes',duplicate_probe='yes'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            platform='unknown_platform',
            unique_genes=['average_genes','high_var','first_gene'],
            num_features='all',annotate='yes',duplicate_probe='yes')),
    Module(
        'select_first_n_genes',
        antecedent(
            "signal_file",logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            platform='unknown_platform',num_features='all',duplicate_probe='yes'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            platform='unknown_platform',
            num_features='int',duplicate_probe='yes')),
    
    
    Module(
        'add_crossplatform_probeid',
        antecedent(
            "signal_file",logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            platform='unknown_platform',duplicate_probe='yes'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            platform='str',duplicate_probe='yes')),
    Module(
        'remove_duplicate_probes',
        antecedent(
            "signal_file",logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            duplicate_probe='yes',platform='str'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            duplicate_probe='high_var_probe',platform='str')),
    Module(
        'select_probe_by_best_match',
        antecedent(
            "signal_file",logged="yes",
            format='tdf', missing=[None, "no", "median", "zero"],
            duplicate_probe='yes',platform='str'),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[None, "no", "median", "zero"],
            duplicate_probe='closest_probe',platform='str')),
##    Module(#this may cause loop
##        'convert_signal_to_gct',
##        antecedent(
##            "signal_file",
##            format='tdf', missing=[None, "no", "median", "zero"]),
##        consequent(
##            "signal_file",
##            format='gct',missing=[None, "no", "median", "zero"]
##            )),
##    Module(#this may cause loop
##        'unlog_signal',
##        antecedent(
##            "signal_file",
##            format='tdf', logged='yes',
##            missing=[None, "no", "median", "zero"]),
##        consequent(
##            "signal_file",
##            format='tdf',missing=[None, "no", "median", "zero"],
##            logged='no')),
    ]

if __name__ == '__main__':
    test_bie()
