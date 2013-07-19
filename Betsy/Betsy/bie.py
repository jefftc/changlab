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
Attribute
Data
Module
Network

"""

# _backchain_to_modules       Given a consequent, list compatible modules.
# _backchain_to_antecedents   Given a consequent and module, list antecedents.
# _calc_consequent_data_compatibility
#
# _print_nothing
# _print_string
# _print_line
# _print_network
# _plot_network_gv
# _pretty_attributes
#
# _find_node
# _get_attribute_type
# _intersection



# The value of an attribute can be:
# NOVALUE  Value is not specified.
# GENERIC  Value can be any string.  Used for parameters like the
#          number of genes to filter out of a list, which can be any
#          number.  Any value can match a GENERIC.
# ITEM     e.g. "mean"
# LIST     e.g. ["mean", "median"]
# 
# 
# Use of GENERIC and LIST.
# 1.  In the antecedent, LIST means that the module can accept
#     multiple values.
#     EX: filter_missing_with_zeros    missing=["yes", "unknown"]
#     EX: convert_to_tdf               format=["pcl", "res", "gct", "jeffs"]
# 2.  In the consequent, LIST means that the module can generate
#     multiple values to fit the needs of the next Data object.
#     This can be used for user-settable parameters.
#     EX: gene_center     gene_normalize=[None, "no"] -> ["mean", "median"]
# 3.  In the consequent, LIST can also mean that the module doesn't
#     change the value.
#     EX: quantile_norm   missing=["median", "zero"] -> ["median", "zero"]
# 4.  In the consequent, GENERIC means that the module can generate
#     many different values to fit the needs of the next DATA object.
#     EX: filter_genes   filter=[None, "no"] -> GENERIC

NOVALUE = "___BETSY_NOVALUE___"
GENERIC = "___BETSY_GENERIC___"

TYPE_GENERIC, TYPE_NOVALUE, TYPE_ITEM, TYPE_LIST = range(4)




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
    # data is a node in the network.  start_data is the Data that the
    # user wants to start on.  In start_data, a LIST
    if data.datatype != start_data.datatype:
        return False

    data_attr = data.attributes
    strt_attr = start_data.attributes
    
    x = data_attr.keys() + strt_attr.keys()
    all_attributes = sorted({}.fromkeys(x))

    # A LIST in START_TYPE indicates multiple possible starting
    # places.  Does not have to match exactly.  (Because in
    # antecedent, module can accept multiple values.)

    # CASE  DATA_TYPE  START_TYPE  RESULT
    #   1    NOVALUE    NOVALUE    OK.  START does not care.
    #   2    NOVALUE    GENERIC    Incompatible.
    #   3    NOVALUE     ITEM      Incompatible.
    #   4    NOVALUE     LIST      Check if NOVALUE in LIST.
    #   5    GENERIC    NOVALUE    OK.  START does not care.
    #   6    GENERIC    GENERIC    NotImplementedError.
    #   7    GENERIC     ITEM      OK.
    #   8    GENERIC     LIST      NotImplementedError.
    #   9     ITEM      NOVALUE    OK.  START does not care.
    #  10     ITEM      GENERIC    OK.
    #  11     ITEM       ITEM      Check if items are equal.
    #  12     ITEM       LIST      Compatible if item is in start sequence.
    #  13     LIST      NOVALUE    OK.  START does not care.
    #  14     LIST      GENERIC    Incompatible.  Value in network is
    #                              ambiguous, but start requires a
    #                              specific (but GENERIC) value.
    #  15     LIST       ITEM      OK if ITEM in LIST.
    #  16     LIST       LIST      Compatible if intersection in lists.

    compatible = True
    for key in all_attributes:
        DATA_VALUE = data_attr.get(key)
        STRT_VALUE = strt_attr.get(key)
        DATA_TYPE = _get_attribute_type(data_attr, key)
        STRT_TYPE = _get_attribute_type(strt_attr, key)

        # Case 1.
        if DATA_TYPE == TYPE_NOVALUE and STRT_TYPE == TYPE_NOVALUE:
            pass
        # Cases 2-3.
        elif DATA_TYPE == TYPE_NOVALUE:
            compatible = False
        # Case 4.
        elif DATA_TYPE == TYPE_NOVALUE and STRT_TYPE == TYPE_LIST:
            if DATA_VALUE not in STRT_VALUE:
                compatible = False
        # Case 5.
        elif DATA_TYPE == TYPE_GENERIC and STRT_TYPE == TYPE_NOVALUE:
            raise NotImplementedError
        # Case 6.
        elif DATA_TYPE == TYPE_GENERIC and STRT_TYPE == TYPE_GENERIC:
            raise NotImplementedError
        # Cases 7.
        elif DATA_TYPE == TYPE_GENERIC and STRT_TYPE == TYPE_ITEM:
            pass
        # Case 8.
        elif DATA_TYPE == TYPE_GENERIC and STRT_TYPE == TYPE_LIST:
            raise NotImplementedError
        # Cases 9-10.
        elif DATA_TYPE == TYPE_ITEM and STRT_TYPE in [
            TYPE_NOVALUE, TYPE_GENERIC]:
            pass
        # Case 11.
        elif DATA_TYPE == TYPE_ITEM and STRT_TYPE == TYPE_ITEM:
            if DATA_VALUE != STRT_VALUE:
                compatible = False
        # Case 12.
        elif DATA_TYPE == TYPE_ITEM and STRT_TYPE == TYPE_LIST:
            if DATA_VALUE not in STRT_VALUE:
                compatible = False
        # Case 13.
        elif DATA_TYPE == TYPE_LIST and STRT_TYPE == TYPE_NOVALUE:
            pass
        # Case 14.
        elif DATA_TYPE == TYPE_LIST and STRT_TYPE == TYPE_GENERIC:
            compatible = False
        # Case 15.
        elif DATA_TYPE == TYPE_LIST and STRT_TYPE == TYPE_ITEM:
            if STRT_VALUE not in DATA_VALUE:
                compatible = False
        # Case 16.
        elif DATA_TYPE == TYPE_LIST and STRT_TYPE == TYPE_LIST:
            #if sorted() != sorted(STRT_VALUE):
            if not _intersection(DATA_VALUE, STRT_VALUE):
                compatible = False
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
    #in_data = make_data("gse_id", platform='unknown')
    #in_data = make_data("gse_id")
    in_data = make_data("signal_file", logged="yes")
    #out_data = make_data(
    #    "signal_file", preprocess='rma', logged='yes')
    #out_data = make_data(
    #    "signal_file", preprocess='rma', logged='yes', filter='no',
    #    format='tdf')
    out_data = make_data(
        "signal_file", format='tdf', preprocess='rma', logged='yes',
        filter='no', missing='zero', predataset='no', rename_sample='no',
        gene_center='no', gene_normalize='no', quantile_norm='yes')
    
    network = backchain(all_modules, out_data)
    network = prune_network_by_start(network, in_data)
    
    print "DONE FINDING PIPELINES"
    _print_network(network)
    _plot_network_gv("out.png", network)


def _backchain_to_modules(moduledb, data):
    # Return list of modules that that can generate a consequent that
    # is compatible with data.  The modules will be sorted in
    # decreasing order of the number of attributes that are
    # compatible.

    matches = []  # list of (module, num compatible attributes)
    for module in moduledb:
        num_attributes = _calc_consequent_data_compatibility(module, data)
        if not num_attributes:
            continue
        x = module, num_attributes
        matches.append(x)

    # Sort by decreasing number of attributes provided.
    matches.sort(key=lambda x: -x[-1])
    modules = [x[0] for x in matches]
    return modules


def _backchain_to_antecedents(module, data):
    # Return list of Data objects that can be the antecedents of module.

    # The datatype is from the antecedent of the module.
    datatype = module.ante_data.datatype

    # Back chain the attributes.  Possibilities:
    #
    # DATA_VALUE  Value of attribute in data.
    # ANTE_VALUE  Value of attribute in antecedent.
    # CONS_VALUE  Value of attribute in consequent.
    # ANTE_TYPE   Type of attribute in antecedent.
    # CONS_TYPE   Type of attribute in consequent.
    #
    # CASE  ANTE_TYPE  CONS_TYPE  RESULT
    #   1    NOVALUE    NOVALUE   DATA_VALUE.  Attribute not for this module.
    #   2    NOVALUE    GENERIC   Ignore.  Attribute created by module.
    #   3    NOVALUE     ITEM     Ignore.  preprocess="rma"
    #   4    NOVALUE     LIST     Ignore.
    #   5    GENERIC    NOVALUE   NotImplementedError. When does this happen?
    #   6    GENERIC    GENERIC   NotImplementedError. When does this happen?
    #   7    GENERIC     ITEM     NotImplementedError. When does this happen?
    #   8    GENERIC     LIST     NotImplementedError. When does this happen?
    #   9     ITEM      NOVALUE   ANTE_VALUE.  cel_version="v3_4"
    #  10     ITEM      GENERIC   ANTE_VALUE.
    #  11     ITEM       ITEM     ANTE_VALUE.  logged="no"->"yes"
    #  12     ITEM       LIST     ANTE_VALUE.
    #  13     LIST      NOVALUE   ANTE_VALUE.
    #  14     LIST      GENERIC   ANTE_VALUE.
    #  15     LIST       ITEM     ANTE_VALUE.  format=["cdt","tdf"] -> "tdf"
    #  16     LIST       LIST     DATA_VALUE or ANTE_VALUE.
    
    data_attr = data.attributes
    ante_attr = module.ante_data.attributes
    cons_attr = module.cons_data.attributes

    x = data_attr.keys() + ante_attr.keys() + cons_attr.keys()
    all_attributes = sorted({}.fromkeys(x))

    attributes = {}
    for key in all_attributes:
        DATA_VALUE = data_attr.get(key)
        ANTE_VALUE = ante_attr.get(key)
        CONS_VALUE = cons_attr.get(key)
        DATA_TYPE = _get_attribute_type(data_attr, key)
        ANTE_TYPE = _get_attribute_type(ante_attr, key)
        CONS_TYPE = _get_attribute_type(cons_attr, key)

        # Case 1.
        if ANTE_TYPE == TYPE_NOVALUE and CONS_TYPE == TYPE_NOVALUE:
            attributes[key] = DATA_VALUE
        # Cases 2-4.
        elif ANTE_TYPE == TYPE_NOVALUE:
            pass
        # Cases 5-8.
        elif ANTE_TYPE == TYPE_GENERIC:
            raise NotImplementedError
        # Case 14.
        elif ANTE_TYPE == TYPE_LIST and CONS_TYPE == TYPE_GENERIC:
            # ANTE  filter=[None, 'no']
            # CONS  filter='___GENERIC___'
            # DATA  filter='0.50'
            attributes[key] = ANTE_VALUE
        # Case 16.
        elif ANTE_TYPE == TYPE_LIST and CONS_TYPE == TYPE_LIST:
            # This can happen if:
            # 1.  The module doesn't change the value.
            # 2.  The module can generate multiple values to fit the
            #     next Data object.
            # If it's the first case, use DATA_VALUE, otherwise use
            # ANTE_VALUE.

            assert DATA_TYPE == TYPE_ITEM
            assert DATA_VALUE in CONS_VALUE
            # Case 1.
            if sorted(ANTE_VALUE) == sorted(CONS_VALUE):
                attributes[key] = DATA_VALUE
            # Case 2.
            else:
                attributes[key] = ANTE_VALUE
        # Cases 9-13,15.
        else:
            attributes[key] = ANTE_VALUE

    x = Data(datatype, attributes)
    return [x]


def _calc_consequent_data_compatibility(module, data):
    # Return the number of attributes in data that this module can
    # generate.  Look for attributes in consequent that match
    # attributes in data.  0 means that module can not produce data.
    
    p = _print_nothing
    #p = _print_string
    
    p("Evaluating module %s." % module.name)

    # If this module doesn't produce the right data type, then
    # ignore it.
    if module.cons_data.datatype != data.datatype:
        p("Not right data type [%s %s]." % (
            module.cons_data.datatype, data.datatype))
        return 0

    data_attr = data.attributes
    ante_attr = module.ante_data.attributes
    cons_attr = module.cons_data.attributes

    # Count the number of attributes in the consequent that is
    # compatible with the attributes for data.
    #
    # A module is disqualified if it generates something that is
    # incompatible with the data (so there's no way it can generate
    # this data).
    #
    # CASE  CONS_TYPE  DATA_TYPE  RESULT
    #   1    NOVALUE    NOVALUE   +0.  Module does not care.
    #   2    NOVALUE    GENERIC   +0.  Module does not generate this.
    #   3    NOVALUE     ITEM     +0.
    #   4    NOVALUE     LIST     +0.
    #   5    GENERIC    NOVALUE   +0.
    #   6    GENERIC    GENERIC   +1.
    #   7    GENERIC     ITEM     +1 if ANTE_VALUE != DATA_VALUE, otherwise DQ.
    #   8    GENERIC     LIST     Disqualify.  GENERIC can't match LIST.
    #   9     ITEM      NOVALUE   +0.
    #  10     ITEM      GENERIC   NotImplementedError.
    #  11     ITEM       ITEM     +1 if same, otherwise DQ.
    #  12     ITEM       LIST     +1 if ITEM in LIST, otherwise DQ.
    #  13     LIST      NOVALUE   +0.
    #  14     LIST      GENERIC   NotImplementedError.
    #  15     LIST       ITEM     +1 is ITEM in LIST, otherwise DQ.
    #  16     LIST       LIST     +1 if intersect, otherwise DQ.
    #
    # *** The +1 only counts if module changes this attribute.
    # Otherwise, +0.

    num_attributes = 0
    for key in data.attributes:
        p("Evaluating attribute %s." % key)

        DATA_VALUE = data_attr.get(key)
        ANTE_VALUE = ante_attr.get(key)
        CONS_VALUE = cons_attr.get(key)
        ANTE_TYPE = _get_attribute_type(ante_attr, key)
        CONS_TYPE = _get_attribute_type(cons_attr, key)
        DATA_TYPE = _get_attribute_type(data_attr, key)

        module_changes_value = ANTE_VALUE != CONS_VALUE
        if ANTE_TYPE == TYPE_LIST and CONS_TYPE == TYPE_LIST:
            module_changes_value = sorted(ANTE_VALUE) == sorted(CONS_VALUE)

        disqualify = False
        # Cases 1-4.
        if CONS_TYPE == TYPE_NOVALUE:
            pass
        # Case 5.
        elif CONS_TYPE == TYPE_GENERIC and DATA_TYPE == TYPE_NOVALUE:
            pass
        # Case 6.
        elif CONS_TYPE == TYPE_GENERIC and DATA_TYPE == TYPE_GENERIC:
            if module_changes_value:
                p("Attribute is compatible.")
                num_attributes += 1
        # Case 7.
        elif CONS_TYPE == TYPE_GENERIC and DATA_TYPE == TYPE_ITEM:
            if ANTE_TYPE == TYPE_NOVALUE:
                disqualify = True
            elif ANTE_TYPE == TYPE_GENERIC:
                raise NotImplementedError
            elif ANTE_TYPE == TYPE_ITEM:
                if DATA_VALUE == ANTE_VALUE:
                    disqualify = True
                else:
                    p("Attribute is compatible.")
                    num_attributes += 1
            elif ANTE_TYPE == TYPE_LIST:
                if DATA_VALUE in ANTE_VALUE:
                    disqualify = True
                else:
                    p("Attribute is compatible.")
                    num_attributes += 1
        # Case 8.
        elif CONS_TYPE == TYPE_GENERIC and DATA_TYPE == TYPE_LIST:
            disqualify = True
        # Case 9.
        elif CONS_TYPE == TYPE_ITEM and DATA_TYPE == TYPE_NOVALUE:
            pass
        # Case 10.
        elif CONS_TYPE == TYPE_ITEM and DATA_TYPE == TYPE_GENERIC:
            raise NotImplementedError
        # Case 11.
        elif CONS_TYPE == TYPE_ITEM and DATA_TYPE == TYPE_ITEM:
            if CONS_VALUE == DATA_VALUE:
                if module_changes_value:
                    p("Attribute is compatible.")
                    num_attributes += 1
            else:
                disqualify = True
        # Case 12.
        elif CONS_TYPE == TYPE_ITEM and DATA_TYPE == TYPE_LIST:
            # E.g. preprocess_rma generates format="jeffs", data has
            # format ["pcl", "res", "gct", "jeffs"].  Keep if the
            # CONS_VALUE in DATA_VALUE.  Otherwise, DQ.
            #print "HERE 2"
            #print key
            #print CONS_VALUE
            #print DATA_VALUE
            if CONS_VALUE in DATA_VALUE:
                if module_changes_value:
                    p("Attribute is compatible.")
                    num_attributes += 1
            else:
                disqualify = True
        # Case 13.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_NOVALUE:
            pass
        # Case 14.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_GENERIC:
            raise NotImplementedError
        # Case 15.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_ITEM:
            # missing=["yes","unknown"]  missing="no"
            if DATA_VALUE in CONS_VALUE:
                p("Attribute is compatible.")
                num_attributes += 1
            else:
                disqualify = True
        # Case 16.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_LIST:
            if _intersection(CONS_VALUE, DATA_VALUE):
                if module_changes_value:
                    p("Attribute is compatible.")
                    num_attributes += 1
            else:
                disqualify = True
                
        if disqualify:
            p("Attribute is disqualified.")
            num_attributes = 0
            break
        
        
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


## def _is_item(x):
##     import operator
##     if type(x) is type(""):
##         return True
##     if not operator.isSequenceType(x):
##         return True
##     return False


## def _item2list(x):
##     if _is_item(x):
##         return [x]
##     return x


def _get_attribute_type(attributes, name):
    import operator

    x = attributes.get(name)

    if name not in attributes:
        return TYPE_NOVALUE
    elif x == GENERIC:
        return TYPE_GENERIC
    elif type(x) is type(""):
        return TYPE_ITEM
    elif operator.isSequenceType(x):
        return TYPE_LIST
    raise AssertionError, "Unknown attribute type: %s" % str(x)


def _intersection(x, y):
    return list(set(x).intersection(y))


all_modules = [
    # cel_files
    Module(
        "download_geo_GSEID",
        antecedent("gse_id"),
        consequent("cel_files")),
    Module(
        "download_geo_GSEID_GPLID",
        antecedent("gse_id_and_platform"),
        consequent("cel_files")),
    Module(
        "extract_CEL_files",
        antecedent("cel_files", cel_version="unknown"),
        consequent("cel_files", cel_version="cc_or_v3_4")),
    Module(
        "convert_CEL_to_v3_4",
        # XXX FIX cel_version
        antecedent("cel_files", cel_version="cc_or_v3_4"),
        consequent("cel_files", cel_version="v3_4")),
    Module(
        "preprocess_rma",
        antecedent("cel_files", cel_version="v3_4"),
        consequent(
            "signal_file", logged="yes", preprocess="rma", format="jeffs",
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
            "signal_file", format="tdf", logged="yes",
            preprocess="agilent", missing="yes")),
    
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
            preprocess='illumina', missing='unknown')),
    Module(
        "get_illumina_signal",
        antecedent("illu_folder", preprocess='illumina'),
        consequent(
            "signal_file", preprocess='illumina', logged="no",
            missing='unknown')),

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
            "signal_file", format="tdf", logged="yes",
            preprocess="loess", missing="unknown")),
    
    #-----------------------------------------------------------------------
    Module(
        "filter_genes_by_missing_values",
        antecedent(
            "signal_file", format='tdf', logged="yes",
            missing=["yes", "unknown"], filter=[NOVALUE, "no"]),
        consequent(
            "signal_file", format='tdf', logged="yes",
            missing=["yes", "unknown"], filter=GENERIC)),
    Module(
        "fill_missing_with_median",
        antecedent(
            "signal_file", format='tdf', logged="yes",
            missing=["yes", "unknown"]),
        consequent(
            "signal_file", format='tdf', logged="yes", missing="median")),
    Module(
        "fill_missing_with_zeros",
        antecedent(
            "signal_file", format='tdf', logged="yes",
            missing=["yes", "unknown"]),
        consequent("signal_file", format='tdf', logged="yes", missing="zero")),
    Module(
        "convert_signal_to_tdf",
        antecedent(
            "signal_file",
            format=['pcl', 'res', 'gct', 'jeffs', 'unknown', 'xls']),
        consequent("signal_file", format='tdf')),
    Module(
        "log_signal",
        antecedent("signal_file", logged=[NOVALUE, "no"], format='tdf'),
        consequent("signal_file", logged="yes", format='tdf')),
    Module(
        "filter_and_threshold_genes",
        antecedent(
            "signal_file", logged=[NOVALUE, "no"], format='tdf',
            predataset=[NOVALUE, 'no']),
        consequent(
            "signal_file", logged=[NOVALUE, "no"], format='tdf',
            predataset="yes")),
    Module( #require a rename_list_file
        "relabel_samples",
        antecedent(
            "signal_file", logged="yes", format='tdf',
            missing=[NOVALUE, "no", "median", "zero"],
            rename_sample=[NOVALUE, "no"]),
        consequent(
            "signal_file", logged="yes", format='tdf',
            missing=[NOVALUE, "no", "median", "zero"], rename_sample="yes")),
    
    #------------------------------------------------------------------
    Module(
        "quantile_norm",
        antecedent(
            "signal_file", quantile_norm=[NOVALUE, "no"], format='tdf',
            logged="yes", missing=[NOVALUE, "no", "median", "zero"]),
        consequent(
            "signal_file", quantile_norm="yes", format='tdf',
            logged="yes", missing=[NOVALUE, "no", "median", "zero"])),
    
    #------------------------------------------------------------------
    Module(
        "gene_center",
        antecedent(
            "signal_file", gene_center=[NOVALUE, "no"], logged="yes",
            gene_normalize=[NOVALUE, "no"], format='tdf',
            missing=[NOVALUE, "no", "median", "zero"]),
        consequent(
            "signal_file", gene_center=["mean", "median"], logged="yes",
            gene_normalize=[NOVALUE, "no"], format='tdf',
            missing=[NOVALUE, "no", "median", "zero"])),
    Module(
        "gene_normalize",
        antecedent(
            "signal_file", gene_normalize=[NOVALUE, "no"], logged="yes",
            format='tdf', missing=[NOVALUE, "no", "median", "zero"]),
        consequent(
            "signal_file", gene_normalize=["variance", "sum_of_squares"],
            logged="yes", format='tdf',
            missing=[NOVALUE, "no", "median", "zero"])),
    ]
   

if __name__ == '__main__':
    test_bie()


