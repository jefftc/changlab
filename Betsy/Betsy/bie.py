"""

Glossary:
fact         Starting point.
goal         Desired outcome.
rule         Another word is proposition.
consequent   Second half of a rule.
antecedent   First half of a rule.

Modules "convert" one DataType to another.
Modules "produce" a Data node.


Functions:
backchain
prune_network_by_start

summarize_moduledb
print_modules

test_bie


Classes:
DataType
Data
Module
ModuleDbSummary
Network

"""

# _make_goal
# 
# _is_compatible_with_start
#
# _backchain_to_modules       Given a consequent, find compatible Modules.
# _backchain_to_antecedent    Given a consequent and module, find antecedents.
# _can_module_produce_data
# _can_converting_module_produce_data
# _can_nonconverting_module_produce_data
# 
# _find_duplicate_modules
# _merge_duplicate_modules
# 
# _find_data_node
# _find_module_node
# _get_attribute_type
# _intersection
# 
# _print_nothing
# _print_string
# _print_line
# _print_network
# _plot_network_gv
# _pretty_attributes


# The value of an attribute can be:
# ATOM      e.g. "mean"
# LIST      e.g. ["mean", "median"]
# ANYATOM   Value can be any string.  Used for parameters like the
#           number of genes to filter out of a list, which can be any
#           number.  Any string can match ATOM.
# NOVALUE   Value is not specified.  Meaning is ambiguous (see below).
#           Can be default or don't care.
#
# Module Node: Antecedent
# ATOM      The module requires a specific value.
# LIST      The module can take any of these values.
#           E.g. filter_missing_with_zeros    missing=["yes", "unknown"]
# ANYATOM   UNDEFINED.  When is this used?
# NOVALUE   This attribute is irrelevant for the module.
# 
# Module Node: Consequent
# ATOM      The module generates a specific value.
# LIST      The module can generate any one of these values to fit the
#           requirements of the next Data node.
#           e.g. gene_center  gene_normalize=[None, "no"] -> ["mean", "median"]
#           Also is seen if the module does not change this value.
#           e.g. quantile_norm missing=["median", "zero"] -> ["median", "zero"]
# ANYATOM   The module can generate any value to fit the requirements of
#           the next Data node.
# NOVALUE   If the antecedent is the same DataType, then it must also
#           have NOVALUE.  If this is the case, then this attribute is
#           irrelevant to the module.
#           If the antecedent is a different DataType, then this
#           module produces the default value.
#
# Data Node
# ATOM      A specific value.
# LIST      The actual value is unknown, but can be one of these options.
# ANYATOM   UNDEFINED.  Is this needed?
# NOVALUE   Default value.


NOVALUE = "___BETSY_NOVALUE___"
ANYATOM = "___BETSY_ANYATOM___"

TYPE_ANYATOM, TYPE_NOVALUE, TYPE_ATOM, TYPE_LIST = range(4)


class DataType:
    # XXX for at LIST, first value is default.
    def __init__(self, name, **attributes):
        # Make sure the attributes are value.
        for key, values in attributes.iteritems():
            attr_type = _get_attribute_type(attributes, key)
            # Values can be either LIST or ANYATOM.
            assert attr_type in [TYPE_LIST, TYPE_ANYATOM], "%s %s" % (
                name, key)
        self.name = name
        self.attributes = attributes.copy()
    def get_defaults(self):
        # Return a dictionary of the attributes to their values.
        defaults = {}
        for key, values in self.attributes.iteritems():
            attr_type = _get_attribute_type(self.attributes, key)
            if attr_type == TYPE_ANYATOM:
                defaults[key] = ANYATOM
            elif attr_type == TYPE_LIST:
                defaults[key] = values[0]
            else:
                raise AssertionError, "Invalid attribute type: %s %s" % (
                    key, values)
        return defaults
    def __call__(self, **attributes):
        # Create a Data object.
        for key in attributes:
            assert key in self.attributes, "Unknown attribute for %s: %s" % (
                self.name, key)
        return Data(self, attributes)


class Data:
    # Members:
    # datatype     Datatype object.
    # attributes   Dict of key -> value or list of values.

    def __init__(self, datatype, attributes):
        # Check the attributes.
        # CASE  DATA_TYPE  DATATYPE_TYPE  RESULT
        #   1    NOVALUE     NOVALUE      ERROR.  DATA can't be NOVALUE.
        #   2    NOVALUE     ANYATOM      ERROR.  
        #   3    NOVALUE       ATOM       ERROR.
        #   4    NOVALUE       LIST       ERROR.
        #   5    ANYATOM     NOVALUE      ERROR.  DATATYPE can't be NOVALUE.
        #   6    ANYATOM     ANYATOM      OK.
        #   7    ANYATOM       ATOM       Bad.
        #   8    ANYATOM       LIST       Bad.
        #   9      ATOM      NOVALUE      ERROR.
        #  10      ATOM      ANYATOM      OK.
        #  11      ATOM        ATOM       OK if ATOMs equal.
        #  12      ATOM        LIST       OK if ATOM in LIST.
        #  13      LIST      NOVALUE      ERROR.
        #  14      LIST      ANYATOM      Bad.
        #  15      LIST        ATOM       Bad.
        #  16      LIST        LIST       OK if DATA LIST is subset.
        data_attr = attributes
        dtyp_attr = datatype.attributes
        for key, values in attributes.iteritems():
            assert key in datatype.attributes, "Unknown key: %s" % key
            
            DATA_VALUE = data_attr.get(key)
            DTYP_VALUE = dtyp_attr.get(key)
            DATA_TYPE = _get_attribute_type(data_attr, key)
            DTYP_TYPE = _get_attribute_type(dtyp_attr, key)

            # Cases 1-4.
            if DTYP_TYPE == TYPE_NOVALUE:
                raise AssertionError
            # Cases 5,9,13.
            elif DATA_TYPE == TYPE_NOVALUE:
                raise AssertionError
            # Case 6.
            elif DATA_TYPE == TYPE_ANYATOM and DTYP_TYPE == TYPE_ANYATOM:
                # OK
                pass
            # Cases 7, 8.
            elif DATA_TYPE == TYPE_ANYATOM:
                raise AssertionError, "type mismatch"
            # Case 10.
            elif DATA_TYPE == TYPE_ATOM and DTYP_TYPE == TYPE_ANYATOM:
                pass
            # Case 11.
            elif DATA_TYPE == TYPE_ATOM and DTYP_TYPE == TYPE_ATOM:
                assert DATA_VALUE == DTYP_VALUE
            # Case 12.
            elif DATA_TYPE == TYPE_ATOM and DTYP_TYPE == TYPE_ANYATOM:
                assert DATA_VALUE in DTYP_VALUE
            # Case 14.
            elif DATA_TYPE == TYPE_LIST and DTYP_TYPE == TYPE_ANYATOM:
                raise AssertionError
            # Case 15.
            elif DATA_TYPE == TYPE_LIST and DTYP_TYPE == TYPE_ATOM:
                raise AssertionError
            # Case 16.
            elif DATA_TYPE == TYPE_LIST and DTYP_TYPE == TYPE_LIST:
                for x in DATA_VALUE:
                    assert x in DTYP_VALUE
        self.datatype = datatype
        self.attributes = attributes.copy()
    def copy(self):
        return Data(self.datatype, self.attributes)
    def __cmp__(self, other):
        x1 = [self.datatype, self.attributes]
        x2 = [other.datatype, other.attributes]
        return cmp(x1, x2)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            self.datatype.name,
            _pretty_attributes(self.attributes),
            ]
        x = [x for x in x if x]
        return "Data(%s)" % ", ".join(x)


class Module:
    def __init__(self, name, ante_data, cons_data, **parameters):
        self.name = name
        # If a DataType is provided, convert it into a Data object
        # with no attributes.
        if isinstance(ante_data, DataType):
            ante_data = ante_data()
        if isinstance(cons_data, DataType):
            cons_data = cons_data()
        self.ante_data = ante_data.copy()
        self.cons_data = cons_data.copy()
        self.parameters = parameters.copy()
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


class ModuleDbSummary:
    # module_names    List of module names (strings).
    # name2module     Dictionary of module name (string) to Module object.
    # name2datatypes  Dict of module name to (ante Datatype, cons Datatype).
    # datatypes       List of Datatype objects.
    def __init__(
        self, module_names, name2module, name2datatypes, datatypes):
        self.module_names = module_names[:]
        self.name2module = name2module.copy()
        self.name2datatypes = name2datatypes.copy()
        self.datatypes = datatypes[:]
        

class Network:
    def __init__(self, nodes, transitions):
        # nodes should be a list of Data or Module objects.  Data
        # transition to Modules, and Modules transition to Data.
        
        # Make sure nodes are Data or Module objects.
        for n in nodes:
            assert isinstance(n, Data) or isinstance(n, Module)
        # Make sure transitions point to the right types of objects.
        for node_id, next_ids in transitions.iteritems():
            assert node_id >= 0 and node_id < len(nodes)
            for nid in next_ids:
                assert nid >= 0 and nid < len(nodes)
            for nid in next_ids:
                n1 = nodes[node_id]
                n2 = nodes[nid]
                if isinstance(n1, Data):
                    assert isinstance(n2, Module)
                else:
                    assert isinstance(n2, Data)
        
        self.nodes = nodes[:]
        self.transitions = transitions.copy()

    def iterate(self, node_class=None):
        # Yield tuple of (node_id, next_node_ids).  node_class is the
        # class of the node of the network (either Data or Module).  If
        # provided, will only iterate over that kind of nodes.
        assert node_class in [Data, Module]
        for node_id, node in enumerate(self.nodes):
            if node_class and not isinstance(node, node_class):
                continue
            next_ids = self.transitions.get(node_id, [])
            yield node_id, next_ids

    def delete_node(self, node_id):
        """Delete a node and return a new Network object."""
        assert node_id < len(self.nodes)
        nodes = self.nodes[:]
        nodes.pop(node_id)

        transitions = {}
        for nid, next_ids in self.transitions.iteritems():
            if nid == node_id:
                continue
            elif nid > node_id:
                nid = nid - 1
            next_ids = [x for x in next_ids if x != node_id]
            for i in range(len(next_ids)):
                if next_ids[i] > node_id:
                    next_ids[i] = next_ids[i] - 1
            transitions[nid] = next_ids
        return Network(nodes, transitions)

    def merge_nodes(self, node_ids):
        """node_ids is a list of the indexes of nodes.  Replace all
        these nodes with just a single one.  Returns a new Network
        object."""
        node_ids = sorted(node_ids)

        # Make sure nodes are the same type.
        for i in range(1, len(node_ids)):
            n1 = self.nodes[node_ids[0]]
            n2 = self.nodes[node_ids[i]]
            assert n1.__class__ == n2.__class__

        # Keep the first node, and delete the rest.

        # Make pointers to node_ids all point to the first one.
        transitions = {}
        for nid, next_ids in self.transitions.iteritems():
            next_ids = next_ids[:]
            for i in range(len(next_ids)):
                if next_ids[i] in node_ids[1:]:
                    next_ids[i] = node_ids[0]
            next_ids = sorted({}.fromkeys(next_ids))  # no duplicates
            transitions[nid] = next_ids

        # Make the first node point to all the next_nodes of the other
        # node_ids.
        for i in range(1, len(node_ids)):
            nid0 = node_ids[0]
            nidi = node_ids[i]
            x = transitions.get(nid0, []) + transitions.get(nidi, [])
            x = sorted({}.fromkeys(x))
            transitions[nid0] = x

        merged = Network(self.nodes, transitions)
        for i in range(len(node_ids)-1, 0, -1):
            merged = merged.delete_node(node_ids[i])
        return merged


def backchain(moduledb, goal_datatype, goal_attributes):
    # Return a Network object.
    assert isinstance(goal_datatype, DataType)
    assert type(goal_attributes) is type({})

    goal_data = _make_goal(goal_datatype, goal_attributes)
    
    nodes = []        # list of Data or Module objects.
    transitions = {}  # list of index -> list of indexes

    MAX_NETWORK_SIZE = 1024
    nodes.append(goal_data)
    stack = [0]
    seen = {}
    while stack:
        assert len(nodes) < MAX_NETWORK_SIZE, "network too large"

        # Pop the next node off the stack.
        node_id = stack.pop()
        assert node_id < len(nodes)
        node = nodes[node_id]

        # If I've already seen this node, then don't process it again.
        if node_id in seen:
            continue
        seen[node_id] = 1

        if isinstance(node, Data):
            # Backwards chain to the previous module.
            modules = _backchain_to_modules(moduledb, node)
            for m in modules:
                nodes.append(m)
                m_id = len(nodes) - 1
                stack.append(m_id)
                transitions[m_id] = transitions.get(m_id, [])
                transitions[m_id].append(node_id)
        elif isinstance(node, Module):
            cons_id = transitions[node_id][0]
            cons = nodes[cons_id]
            d = _backchain_to_antecedent(node, cons, goal_attributes)
            d_id = _find_data_node(nodes, d)
            if d_id == -1:
                nodes.append(d)
                d_id = len(nodes)-1
            stack.append(d_id)
            transitions[d_id] = transitions.get(d_id, [])
            transitions[d_id].append(node_id)
        else:
            raise AssertionError, "Unknown node type: %s" % node

    # Remove the duplicates from transitions.
    for nid in transitions:
        next_ids = sorted({}.fromkeys(transitions[nid]))
        transitions[nid] = next_ids
        
    network = Network(nodes, transitions)
    network = _merge_duplicate_modules(network)
    return network


def prune_network_by_start(network, start_data):
    if isinstance(start_data, DataType):
        start_data = start_data()  # convert to Data
    
    # Look for the nodes that are compatible with start_data.
    node_ids = []  # list of node_ids.
    for node_id, next_ids in network.iterate(node_class=Data):
        if _is_compatible_with_start(network.nodes[node_id], start_data):
            node_ids.append(node_id)

    # For each of these node_ids, do forward chaining to find all
    # nodes that these ones can connect to.
    good_ids = {}
    stack = node_ids[:]
    while stack:
        node_id = stack.pop(0)
        if node_id in good_ids:
            continue
        good_ids[node_id] = 1
        stack.extend(network.transitions.get(node_id, []))

    # Delete all the IDs that aren't in good_ids.
    for nid in range(len(network.nodes)-1, -1, -1):
        if nid not in good_ids:
            network = network.delete_node(nid)

    return network


def summarize_moduledb(moduledb):
    """Take a list of Modules and return a ModuleDbSummary object."""
    name2module = {}   # module_name -> Module
    for module in moduledb:
        assert module.name not in name2module
        name2module[module.name] = module
    module_names = sorted(name2module)

    # module_name -> (ante Datatype, cons Datatype)
    name2datatypes = {}
    for name, module in name2module.iteritems():
        ante_datatype = module.ante_data.datatype
        cons_datatype = module.cons_data.datatype
        name2datatypes[name] = ante_datatype, cons_datatype

    # list of all Data objects in moduledb.
    all_data = []
    for name, module in name2module.iteritems():
        for data in [module.ante_data, module.cons_data]:
            if data not in all_data:
                all_data.append(data)

    datatype2data = {}   # Datatype -> list of Data
    for data in all_data:
        if data.datatype not in datatype2data:
            datatype2data[data.datatype] = []
        datatype2data[data.datatype].append(data)
    datatypes = sorted(datatype2data)

    x = ModuleDbSummary(module_names, name2module, name2datatypes, datatypes)
    return x


def print_modules(moduledb):
    summary = summarize_moduledb(moduledb)

    for name in summary.module_names:
        ante_datatype, cons_datatype = summary.name2datatypes[name]
        x = name, ante_datatype.name, cons_datatype.name
        print "\t".join(x)

    for datatype in summary.datatypes:
        print datatype.name
        for attr in sorted(datatype.attributes):
            values = datatype.attributes[attr]
            if type(values) is type(""):
                values = [values]
            print "  %s : %s" % (attr, ", ".join(map(str, values)))
        print


def _make_goal(datatype, attributes):
    # Make a Data object that can represent a goal.  Use the
    # attributes provided by the user that are appropriate for this
    # datatype.  If anything is not given, use the defaults.
    defaults = datatype.get_defaults()
    attrs = {}
    for key in defaults:
        attrs[key] = attributes.get(key, defaults[key])
    return Data(datatype, attrs)


def _is_compatible_with_start(data, start_data):
    # data is a Data node in the network.  start_data is the Data that
    # the user wants to start on.
    if data.datatype != start_data.datatype:
        return False

    # Start Data
    # ATOM      Must match a specific value.
    # LIST      UNDEFINED.
    # ANYATOM   UNDEFINED.
    # NOVALUE   Use default value.
    # 
    # CASE  DATA_TYPE  START_TYPE  RESULT
    #   1    NOVALUE    NOVALUE    ERROR.  START should have default values.
    #   2    ANYATOM    NOVALUE    ERROR.  START should have default values.
    #   3     ATOM      NOVALUE    ERROR.  START should have default values.
    #   4     LIST      NOVALUE    ERROR.  START should have default values.
    #   5    NOVALUE    ANYATOM    No.
    #   6    ANYATOM    ANYATOM    OK.
    #   7     ATOM      ANYATOM    OK.
    #   8     LIST      ANYATOM    No.
    #   9    NOVALUE     ATOM      No.
    #  10    ANYATOM     ATOM      OK.
    #  11     ATOM       ATOM      Check if items are equal.
    #  12     LIST       ATOM      Check if ATOM in LIST.
    #  13    NOVALUE     LIST      NotImplementedError.
    #  14    ANYATOM     LIST      NotImplementedError.
    #  15     ATOM       LIST      NotImplementedError.
    #  16     LIST       LIST      NotImplementedError.

    data_attr = data.attributes
    strt_attr = start_data.attributes

    # Fill in default values for the start.
    strt_attr = strt_attr.copy()
    for key, value in start_data.datatype.get_defaults().iteritems():
        if key not in strt_attr:
            strt_attr[key] = value
    
    compatible = True
    for key in strt_attr:
        DATA_VALUE = data_attr.get(key)
        STRT_VALUE = strt_attr.get(key)
        DATA_TYPE = _get_attribute_type(data_attr, key)
        STRT_TYPE = _get_attribute_type(strt_attr, key)

        # Cases 1-4.
        if STRT_TYPE == TYPE_NOVALUE:
            raise AssertionError
        # Cases 5, 8.
        elif STRT_TYPE == TYPE_ANYATOM and \
                 DATA_TYPE in [TYPE_NOVALUE, TYPE_LIST]:
            compatible = False
        # Cases 6, 7.
        elif STRT_TYPE == TYPE_ANYATOM and \
                 DATA_TYPE in [TYPE_ANYATOM, TYPE_ATOM]:
            pass
        # Cases 13-16.
        elif STRT_TYPE == TYPE_LIST:
            raise NotImplementedError
        # Case 9.
        elif DATA_TYPE == TYPE_NOVALUE:
            compatible = False
        # Case 10.
        elif DATA_TYPE == TYPE_ANYATOM:
            pass
        # Case 11.
        elif DATA_TYPE == TYPE_ATOM:
            if DATA_VALUE != STRT_VALUE:
                compatible = False
        # Case 12.
        elif DATA_TYPE == TYPE_LIST:
            if STRT_VALUE not in DATA_VALUE:
                compatible = False
    return compatible
            
        
def _backchain_to_modules(moduledb, data):
    # Return list of modules that that can generate a consequent that
    # is compatible with data.  The modules will be sorted in
    # decreasing order of the number of attributes that are
    # compatible.

    matches = []  # list of (module, num compatible attributes)
    for module in moduledb:
        num_attributes = _can_module_produce_data(module, data)
        if not num_attributes:
            continue
        x = module, num_attributes
        matches.append(x)

    # Sort by decreasing number of attributes provided.
    matches.sort(key=lambda x: -x[-1])
    modules = [x[0] for x in matches]
    return modules


def _backchain_to_antecedent(module, data, goal_attributes):
    # Return the Data object that is the antecedent of the module and
    # the consequent (data).  goal_attributes are attributes provided
    # by the user for the goal.  This is necessary because some of
    # these attributes may be relevant to new data types.  E.g. The
    # goal datatype is a signal_file, but some of the attributes are
    # relevant for illu_folder.

    # Back chain the attributes.  Possibilities:
    #
    # DATA_VALUE  Value of attribute in data.
    # ANTE_VALUE  Value of attribute in antecedent.
    # CONS_VALUE  Value of attribute in consequent.
    # ANTE_TYPE   Type of attribute in antecedent.
    # CONS_TYPE   Type of attribute in consequent.
    #
    # CASE  ANTE_TYPE  CONS_TYPE  RESULT
    #   1    NOVALUE    NOVALUE   DATA_VALUE.  Irrelevant for module.
    #   2    NOVALUE    ANYATOM   Ignore (NOVALUE).
    #   3    NOVALUE      ATOM    Ignore (NOVALUE).
    #   4    NOVALUE      LIST    Ignore (NOVALUE).
    #   5    ANYATOM    NOVALUE   NotImplementedError. When does this happen?
    #   6    ANYATOM    ANYATOM   NotImplementedError. When does this happen?
    #   7    ANYATOM      ATOM    NotImplementedError. When does this happen?
    #   8    ANYATOM      LIST    NotImplementedError. When does this happen?
    #   9     ATOM      NOVALUE   ANTE_VALUE.  cel_version="v3_4"
    #  10     ATOM      ANYATOM   ANTE_VALUE.
    #  11     ATOM        ATOM    ANTE_VALUE.  logged="no"->"yes"
    #  12     ATOM        LIST    ANTE_VALUE.
    #  13     LIST      NOVALUE   ANTE_VALUE.
    #  14     LIST      ANYATOM   ANTE_VALUE.
    #  15     LIST        ATOM    ANTE_VALUE.
    #  16     LIST        LIST    DATA_VALUE or ANTE_VALUE (see below).
    
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
            assert DATA_TYPE is not TYPE_NOVALUE
            attributes[key] = DATA_VALUE
        # Cases 2-4.
        elif ANTE_TYPE == TYPE_NOVALUE:
            pass
        # Cases 5-8.
        elif ANTE_TYPE == TYPE_ANYATOM:
            raise NotImplementedError
        # Case 16.
        elif ANTE_TYPE == TYPE_LIST and CONS_TYPE == TYPE_LIST:
            # This can happen if:
            # 1.  The module doesn't change the value.
            # 2.  The module can generate multiple values to fit the
            #     next Data object.
            # If it's the first case, use DATA_VALUE, otherwise use
            # ANTE_VALUE.
            assert DATA_TYPE in [TYPE_ATOM, TYPE_LIST]
            if DATA_TYPE == TYPE_ATOM:
                assert DATA_VALUE in CONS_VALUE
            else:
                assert sorted(DATA_VALUE) == sorted(CONS_VALUE)
            
            # Case 1.
            if sorted(ANTE_VALUE) == sorted(CONS_VALUE):
                # DATA_VALUE is a subset of CONS_VALUE.
                attributes[key] = DATA_VALUE
            # Case 2.
            else:
                attributes[key] = ANTE_VALUE
        # Cases 9-15.
        else:
            assert ANTE_TYPE is not TYPE_NOVALUE
            attributes[key] = ANTE_VALUE

    # If we are converting to a different data type, then add the
    # relevant attributes from the goal_attributes.
    datatype = module.ante_data.datatype
    if module.ante_data.datatype == module.cons_data.datatype:
        data = Data(datatype, attributes)
    else:
        attrs = goal_attributes
        attrs.update(attributes)
        data = _make_goal(datatype, attrs)

    return data


def _can_module_produce_data(module, data):
    # Return the number of attributes in data that this module can
    # generate.  Look for attributes in consequent that match
    # attributes in data.  0 means that module can not produce this
    # data object.
    
    # A module can produce this data if:
    # - The consequent has an attribute value that matches the value
    #   of the attribute in the data.
    # - The consequents and antecedents have no values, and their data
    #   types are different.
    #   e.g. download_geo_GSEID  gseid -> expression_files  (no attributes)
    #
    # A module is incompatible if:
    # - None of the attribute values are the same.
    # - One or more of the attribute values conflict.

    # If this module doesn't produce the same data type, then it can't
    # produce this data object.
    if module.cons_data.datatype != data.datatype:
        return 0

    if module.ante_data.datatype == module.cons_data.datatype:
        return _can_nonconverting_module_produce_data(module, data)
    return _can_converting_module_produce_data(module, data)
    

def _can_converting_module_produce_data(module, data):
    assert module.ante_data.datatype != module.cons_data.datatype

    p = _print_nothing
    #p = _print_string
    
    p("Checking if module %s can produce data." % module.name)

    data_attr = data.attributes
    ante_attr = module.ante_data.attributes
    cons_attr = module.cons_data.attributes

    # If there are no attributes to match, then this matches by
    # default.
    # E.g. extract_CEL_files converts ExpressionFiles to GSEID.
    if not ante_attr and not cons_attr:
        if not data_attr:
            p("Match by no attributes.")
            return 1

    # Fill in default values for the consequent.
    cons_attr = cons_attr.copy()
    for key, value in module.cons_data.datatype.get_defaults().iteritems():
        if key not in cons_attr:
            cons_attr[key] = value


    # CASE  CONS_TYPE  DATA_TYPE   RESULT
    #   1    NOVALUE    NOVALUE    ERROR.
    #   2    NOVALUE    ANYATOM    ERROR.
    #   3    NOVALUE      ATOM     ERROR.
    #   4    NOVALUE      LIST     ERROR.
    #   5    ANYATOM    NOVALUE    DQ.  Must match now.
    #   6    ANYATOM    ANYATOM    DQ.  Must be provided.
    #   7    ANYATOM      ATOM     +1.
    #   8    ANYATOM      LIST     DQ.  ANYATOM doesn't match LIST.
    #   9      ATOM     NOVALUE    DQ.  Must match now.
    #  10      ATOM     ANYATOM    NotImplementedError.
    #  11      ATOM       ATOM     +1 if same, otherwise DQ.
    #  12      ATOM       LIST     +1 if ATOM in LIST, otherwise DQ.
    #  13      LIST     NOVALUE    DQ.  Must match now.
    #  14      LIST     ANYATOM    NotImplementedError.
    #  15      LIST       ATOM     +1 is ATOM in LIST, otherwise DQ.
    #  16      LIST       LIST     +1 if intersect, otherwise DQ.

    num_attributes = 0
    for key in cons_attr:
        p("  Evaluating attribute %s." % key)

        DATA_VALUE = data_attr.get(key)
        ANTE_VALUE = ante_attr.get(key)
        CONS_VALUE = cons_attr.get(key)
        DATA_TYPE = _get_attribute_type(data_attr, key)
        ANTE_TYPE = _get_attribute_type(ante_attr, key)
        CONS_TYPE = _get_attribute_type(cons_attr, key)

        disqualify = False
        # Cases 1-4.
        if CONS_TYPE == TYPE_NOVALUE:
            raise AssertionError, "Should not have NOVALUES."
        # Case 7.
        elif CONS_TYPE == TYPE_ANYATOM and DATA_TYPE == TYPE_ATOM:
            num_attributes += 1
        # Cases 5, 6, 8.
        elif CONS_TYPE == TYPE_ANYATOM:
            disqualify = True
        # Case 9.
        elif CONS_TYPE == TYPE_ATOM and DATA_TYPE == TYPE_NOVALUE:
            disqualify = True
        # Case 10.
        elif CONS_TYPE == TYPE_ATOM and DATA_TYPE == TYPE_ANYATOM:
            raise NotImplementedError
        # Case 11.
        elif CONS_TYPE == TYPE_ATOM and DATA_TYPE == TYPE_ATOM:
            if CONS_VALUE == DATA_VALUE:
                num_attributes += 1
            else:
                disqualify = True
        # Case 12.
        elif CONS_TYPE == TYPE_ATOM and DATA_TYPE == TYPE_LIST:
            if CONS_VALUE in DATA_VALUE:
                num_attributes += 1
            else:
                disqualify = True
        # Case 13.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_NOVALUE:
            disqualify = True
        # Case 14.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_ANYATOM:
            raise NotImplementedError
        # Case 15.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_ATOM:
            if DATA_VALUE in CONS_VALUE:
                num_attributes += 1
            else:
                disqualify = True
        # Case 16.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_LIST:
            if _intersection(CONS_VALUE, DATA_VALUE):
                num_attributes += 1
            else:
                disqualify = True
                
        if disqualify:
            p("    Attribute is disqualified.")
            num_attributes = 0
            break

    if not num_attributes:
        p("  No attributes compatible.")
    else:
        p("  %d attributes compatible." % num_attributes)
    return num_attributes


def _can_nonconverting_module_produce_data(module, data):
    assert module.ante_data.datatype == module.cons_data.datatype
    
    p = _print_nothing
    #p = _print_string
    
    p("Checking if module %s can produce data." % module.name)

    data_attr = data.attributes
    ante_attr = module.ante_data.attributes
    cons_attr = module.cons_data.attributes


    # CASE  CONS_TYPE  DATA_TYPE   RESULT
    #   1    NOVALUE    NOVALUE    +0.  Irrelevant for this module.
    #   2    NOVALUE    ANYATOM    +0.
    #   3    NOVALUE      ATOM     +0.
    #   4    NOVALUE      LIST     +0.
    #   5    ANYATOM    NOVALUE    +0.
    #   6    ANYATOM    ANYATOM    +1.
    #   7    ANYATOM      ATOM     +1 if ANTE_VALUE != DATA_VALUE, otherwise DQ
    #   8    ANYATOM      LIST     DQ.  ANYATOM doesn't match LIST.
    #   9      ATOM     NOVALUE    +0.
    #  10      ATOM     ANYATOM    NotImplementedError.
    #  11      ATOM       ATOM     +1 if same, otherwise DQ.
    #  12      ATOM       LIST     +1 if ATOM in LIST, otherwise DQ.
    #  13      LIST     NOVALUE    +0.
    #  14      LIST     ANYATOM    NotImplementedError.
    #  15      LIST       ATOM     +1 is ATOM in LIST, otherwise DQ.
    #  16      LIST       LIST     +1 if intersect, otherwise DQ.
    #
    # *** The +1 only counts if module changes this attribute.
    # Otherwise, +0.

    num_attributes = 0
    for key in module.cons_data.datatype.attributes:
        p("  Evaluating attribute %s." % key)

        DATA_VALUE = data_attr.get(key)
        ANTE_VALUE = ante_attr.get(key)
        CONS_VALUE = cons_attr.get(key)
        DATA_TYPE = _get_attribute_type(data_attr, key)
        ANTE_TYPE = _get_attribute_type(ante_attr, key)
        CONS_TYPE = _get_attribute_type(cons_attr, key)

        module_changes_value = ANTE_VALUE != CONS_VALUE
        if ANTE_TYPE == TYPE_LIST and CONS_TYPE == TYPE_LIST:
            module_changes_value = sorted(ANTE_VALUE) != sorted(CONS_VALUE)

        disqualify = False
        # Cases 1-4.
        if CONS_TYPE == TYPE_NOVALUE:
            # Since all NOVALUE, ignore this attribute.
            assert ANTE_TYPE == TYPE_NOVALUE
            pass
        # Case 5.
        elif CONS_TYPE == TYPE_ANYATOM and DATA_TYPE == TYPE_NOVALUE:
            pass
        # Case 6.
        elif CONS_TYPE == TYPE_ANYATOM and DATA_TYPE == TYPE_ANYATOM:
            if module_changes_value:
                p("    Attribute is compatible.")
                num_attributes += 1
        # Case 7.
        elif CONS_TYPE == TYPE_ANYATOM and DATA_TYPE == TYPE_ATOM:
            if ANTE_TYPE == TYPE_NOVALUE:
                disqualify = True
            elif ANTE_TYPE == TYPE_ANYATOM:
                raise NotImplementedError
            elif ANTE_TYPE == TYPE_ATOM:
                if DATA_VALUE == ANTE_VALUE:
                    disqualify = True
                else:
                    p("    Attribute is compatible.")
                    num_attributes += 1
            elif ANTE_TYPE == TYPE_LIST:
                if DATA_VALUE in ANTE_VALUE:
                    disqualify = True
                else:
                    p("    Attribute is compatible.")
                    num_attributes += 1
        # Case 8.
        elif CONS_TYPE == TYPE_ANYATOM and DATA_TYPE == TYPE_LIST:
            disqualify = True
        # Case 9.
        elif CONS_TYPE == TYPE_ATOM and DATA_TYPE == TYPE_NOVALUE:
            pass
        # Case 10.
        elif CONS_TYPE == TYPE_ATOM and DATA_TYPE == TYPE_ANYATOM:
            raise NotImplementedError
        # Case 11.
        elif CONS_TYPE == TYPE_ATOM and DATA_TYPE == TYPE_ATOM:
            if CONS_VALUE == DATA_VALUE:
                if module_changes_value:
                    p("    Attribute is compatible [ATOM %s ATOM %s]." % (
                        CONS_VALUE, DATA_VALUE))
                    num_attributes += 1
            else:
                disqualify = True
        # Case 12.
        elif CONS_TYPE == TYPE_ATOM and DATA_TYPE == TYPE_LIST:
            if CONS_VALUE in DATA_VALUE:
                if module_changes_value:
                    num_attributes += 1
            else:
                disqualify = True
        # Case 13.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_NOVALUE:
            pass
        # Case 14.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_ANYATOM:
            raise NotImplementedError
        # Case 15.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_ATOM:
            if DATA_VALUE in CONS_VALUE:
                if module_changes_value:
                    p("    Attribute is compatible [LIST %s ATOM %s]." % (
                        str(CONS_VALUE), DATA_VALUE))
                    num_attributes += 1
            else:
                disqualify = True
        # Case 16.
        elif CONS_TYPE == TYPE_LIST and DATA_TYPE == TYPE_LIST:
            if _intersection(CONS_VALUE, DATA_VALUE):
                if module_changes_value:
                    p("    Attribute is compatible.")
                    num_attributes += 1
            else:
                disqualify = True
                
        if disqualify:
            p("    Attribute is disqualified.")
            num_attributes = 0
            break

    if not num_attributes:
        p("  No attributes compatible.")
    else:
        p("  %d attributes compatible." % num_attributes)
    return num_attributes


def _find_duplicate_modules(network):
    # Return list of module_ids for modules that are duplicated.  If
    # no duplicates found, return an empty list.

    # If the same data node points to two of the same module nodes,
    # then those modules are duplicated.
    duplicates = []
    for data_id, next_ids in network.iterate(node_class=Data):
        if len(next_ids) < 2:
            continue
        
        all_same = True
        for i in range(1, len(next_ids)):
            if network.nodes[next_ids[0]] != network.nodes[ next_ids[i]]:
                all_same = False
                break
        if all_same:
            duplicates = next_ids
            break
    return duplicates
    

def _merge_duplicate_modules(network):
    while 1:
        duplicates = _find_duplicate_modules(network)
        if not duplicates:
            break
        network = network.merge_nodes(duplicates)
    return network


def _find_data_node(nodes, node):
    assert isinstance(node, Data)

    # CASE   N1_TYPE  N2_TYPE  RESULT
    #   1    NOVALUE  NOVALUE  OK.
    #   2    NOVALUE  ANYATOM  No.
    #   3    NOVALUE    ATOM   No.
    #   4    NOVALUE    LIST   No.
    #   5    ANYATOM  NOVALUE  No.
    #   6    ANYATOM  ANYATOM  OK.
    #   7    ANYATOM    ATOM   No.
    #   8    ANYATOM    LIST   No.
    #   9     ATOM    NOVALUE  No.
    #  10     ATOM    ANYATOM  No.
    #  11     ATOM      ATOM   OK if ATOM equal.
    #  12     ATOM      LIST   No.
    #  13     LIST    NOVALUE  No.
    #  14     LIST    ANYATOM  No.
    #  15     LIST      ATOM   No.
    #  16     LIST      LIST   OK if LIST equal.
    
    for i, n in enumerate(nodes):
        if not isinstance(n, Data):
            continue
        if n.datatype != node.datatype:
            continue

        n1_attr = n.attributes
        n2_attr = node.attributes
        x1 = n1_attr.keys()
        x2 = n2_attr.keys()
        all_attributes = {}.fromkeys(x1 + x2)

        is_equal = True
        for key in all_attributes:
            N1_VALUE = n1_attr.get(key)
            N2_VALUE = n2_attr.get(key)
            N1_TYPE = _get_attribute_type(n1_attr, key)
            N2_TYPE = _get_attribute_type(n2_attr, key)

            attr_equal = False
            # Case 1.
            if N1_TYPE == NOVALUE and N2_TYPE == NOVALUE:
                attr_equal = True
            # Case 6.
            elif N1_TYPE == TYPE_ANYATOM and N2_TYPE == TYPE_ANYATOM:
                attr_equal = True
            # Case 11.
            elif N1_TYPE == TYPE_ATOM and N2_TYPE == TYPE_ATOM:
                if N1_VALUE == N2_VALUE:
                    attr_equal = True
            # Case 16.
            elif N1_TYPE == TYPE_LIST and N2_TYPE == TYPE_LIST:
                if sorted(N1_VALUE) == sorted(N2_VALUE):
                    attr_equal = True
            if not attr_equal:
                is_equal = False
        if is_equal:
            return i
    return -1


def _find_module_node(nodes, transitions, node):
    assert isinstance(node, Module)

    for i, n in enumerate(nodes):
        if not isinstance(n, Module):
            continue
        if node == n:
            return i
    return -1


def _get_attribute_type(attributes, name):
    import operator

    x = attributes.get(name)

    if name not in attributes:
        return TYPE_NOVALUE
    elif x == ANYATOM:
        return TYPE_ANYATOM
    elif type(x) is type(""):
        return TYPE_ATOM
    elif operator.isSequenceType(x):
        return TYPE_LIST
    raise AssertionError, "Unknown attribute type: %s" % str(x)


def _intersection(x, y):
    return list(set(x).intersection(y))


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
            x = n.datatype.name
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


def test_bie():
    #print_modules(all_modules); return
    in_data = GSEID
    #in_data = SignalFile(logged="yes", preprocess="rma")
    #x = dict(preprocess="rma", missing_values="no", format="jeffs")

    #print _make_goal(SignalFile, x)
    #return

    #goal_datatype = CELFiles
    #goal_attributes = {}
    #goal_attributes = dict(platform="GPL1691")
    #goal_attributes = dict(version=["v3", "v4"])
    
    goal_datatype = SignalFile
    goal_attributes = dict(format='tdf', preprocess='rma', logged='yes')
    #goal_attributes = dict(
    #    format='tdf', preprocess='rma', logged='yes',
    #    quantile_norm="yes", combat_norm="yes", dwd_norm="yes",
    #    missing_values="no")

    #out_data = make_data(
    #    "signal_file", format='tdf', preprocess='rma', logged='yes',
    #    filter='no', missing='zero', predataset='no', rename_sample='no',
    #    gene_center='no', gene_normalize='no', quantile_norm='yes',
    #    dwd_norm="no", gene_order="class_neighbors", shiftscale_norm="no",
    #    unique_genes="average_genes")
    #out_data = make_data(
    #    "signal_file", format='tdf',preprocess='rma',logged='yes',
    #    filter='no', missing='zero', predataset='no', rename_sample='no',
    #    quantile_norm='yes', dwd_norm='no', gene_order='class_neighbors',
    #    shiftscale_norm='no', combat_norm='no', bfrm_norm='yes',
    #    gene_center='mean', gene_normalize='variance', group_fc='int',
    #    num_features='int', annotate="yes", unique_genes='average_genes',
    #    platform='str', duplicate_probe='high_var_probe')
    
    network = backchain(all_modules, goal_datatype, goal_attributes)
    network = prune_network_by_start(network, in_data)
    
    _print_network(network)
    _plot_network_gv("out.png", network)


AgilentFiles = DataType("AgilentFiles")
CELFiles = DataType(
    "CELFiles",
    version=["unknown", "cc", "v3", "v4"])
ControlFile = DataType(
    "ControlFile",
    preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
    missing_values=["unknown", "no", "yes", "median_fill", "zero_fill"],
    logged=["unknown", "no", "yes"],
    format=["unknown", "tdf", "gct", "jeffs", "pcl", "res", "xls"],
    )
ExpressionFiles = DataType("ExpressionFiles")
GPRFiles = DataType("GPRFiles")
GSEID = DataType("GSEID", platform=ANYATOM)
IDATFiles = DataType("IDATFiles")
ILLUFolder = DataType("ILLUFolder")
SignalFile = DataType(
    "SignalFile",
    format=["unknown", "tdf", "gct", "jeffs", "pcl", "res"],
    preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
    
    # Properties of the data.
    missing_values=["unknown", "no", "yes", "median_fill", "zero_fill"],
    logged=["unknown", "no", "yes"],

    # Normalizing the genes.
    gene_center=["unknown", "no", "mean", "median"],
    gene_normalize=["unknown", "no", "variance", "sum_of_squares"],

    # Normalizing the data.  Very difficult to check normalization.
    # If you're not sure if the data is normalized, then the answer is
    # "no".
    dwd_norm=["no", "yes"],
    bfrm_norm=["no", "yes"],
    quantile_norm=["no", "yes"],
    shiftscale_norm=["no", "yes"],
    combat_norm=["no", "yes"],

    # Annotations.
    #annotate=["no", "yes"],
    #unique_genes=["no", "average_genes", "high_var", "first_gene"],
    #duplicate_probe=["no", "yes", "closest_probe", "high_var_probe"],
    #rename_sample=["no", "yes"],

    # Unclassified.
    #num_features=["all", ANYATOM],
    #gene_order=[
    #    "no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr"],
    #predataset=["no", "yes"],
    #platform=[ANYATOM],
    #filter=["no", ANYATOM],
    #group_fc=["no", ANYATOM],
    )

all_modules = [
    # GSEID
    Module("download_geo", GSEID, ExpressionFiles),
    Module("extract_CEL_files", ExpressionFiles, CELFiles(version="unknown")),

    # CELFiles
    Module(
        "detect_CEL_version", CELFiles(version="unknown"),
        CELFiles(version=["cc", "v3", "v4"])),
    Module(
        "convert_CEL_cc_to_CEL_v3", CELFiles(version="cc"),
        CELFiles(version="v3")),
    Module(
        "preprocess_rma",
        CELFiles(version=["v3", "v4"]),
        SignalFile(logged="yes", preprocess="rma", format="jeffs")
        ),
    Module(
        "preprocess_mas5",
        CELFiles(version=["v3", "v4"]),
        SignalFile(logged="no", preprocess="mas5", format="jeffs")),

    # SignalFile
    Module(
        "convert_signal_to_tdf",
        SignalFile(format=['pcl', 'res', 'gct', 'jeffs', 'unknown']),
        SignalFile(format='tdf')),
    Module(
        "log_signal",
        SignalFile(logged='no', format='tdf'),
        SignalFile(logged='yes', format='tdf')),
    Module(
        "fill_missing_with_median",
        SignalFile(format="tdf", logged="yes", missing_values="yes"),
        SignalFile(format="tdf", logged="yes", missing_values="median_fill")),
    Module(
        "fill_missing_with_zeros",
        SignalFile(format="tdf", logged="yes", missing_values="yes"),
        SignalFile(format="tdf", logged="yes", missing_values="zero_fill")),
    Module(
        "check_for_missing_values",
        SignalFile(format="tdf", missing_values="unknown"),
        SignalFile(format="tdf", missing_values=["no", "yes"])),

    # XXX fix missing values
    # Sample normalization.
    Module(
        "normalize_samples_with_quantile",
        SignalFile(
            format="tdf", logged="yes",
            missing_values=["no", "zero_fill", "median_fill"],
            quantile_norm="no"),
        SignalFile(
            format="tdf", logged="yes",
            missing_values=["no", "zero_fill", "median_fill"],
            quantile_norm="yes")),
    Module(
        "normalize_samples_with_combat",
        SignalFile(
            format="tdf", logged="yes",
            missing_values=["no", "zero_fill", "median_fill"],
            combat_norm="no"),
        SignalFile(
            format="tdf", logged="yes",
            missing_values=["no", "zero_fill", "median_fill"],
            combat_norm="yes")),
    Module(
        "normalize_samples_with_dwd",
        SignalFile(
            format="tdf", logged="yes",
            missing_values=["no", "zero_fill", "median_fill"],
            dwd_norm="no"),
        SignalFile(
            format="tdf", logged="yes",
            missing_values=["no", "zero_fill", "median_fill"],
            dwd_norm="yes")),
    ]


## all_modules = [
##     #cel_files
##     Module(
##         "download_geo_GSEID",
##         antecedent("gse_id"),
##         consequent("expression_files")),
##     Module(
##         "download_geo_GSEID_GPLID",
##         # XXX can platform be parameter for gse_id?
##         antecedent("gse_id_and_platform"),
##         consequent("expression_files")),
##     Module(
##         "extract_CEL_files",
##         antecedent("expression_files"),
##         consequent("cel_files", cel_version="cc_or_v3_4")),
##     Module(
##         "convert_CEL_to_v3_4",
##         # XXX FIX cel_version
##         antecedent("cel_files", cel_version="cc_or_v3_4"),
##         consequent("cel_files", cel_version="v3_4")),
##     Module(
##         "preprocess_rma",
##         antecedent("cel_files", cel_version="v3_4"),
##         consequent(
##             "signal_file", logged="yes", preprocess="rma",
##             format="jeffs", missing="no")),
##     Module(
##         "preprocess_mas5",
##         antecedent("cel_files", cel_version="v3_4"),
##         consequent(
##             "signal_file", logged="no", preprocess="mas5",
##             format="jeffs", missing="no")),
    
##     #-----------------------------------------------------------------------
##     #agilent_files
##     Module(
##         "extract_agilent_files",
##         antecedent("expression_files"),
##         consequent("agilent_files")),
##     Module(
##         "preprocess_agilent",
##         antecedent("agilent_files"),
##         consequent(
##             "signal_file", format="tdf", logged="no", preprocess="agilent",
##             missing="unknown")),
##     #-----------------------------------------------------------------------
##     #idat_files
##     Module(
##         "extract_illumina_idat_files",
##         antecedent("expression_files"),
##         consequent("idat_files")),
##     Module(
##         "preprocess_illumina",
##         antecedent("idat_files"),
##         consequent(
##             "illu_folder", preprocess='illumina', 
##             ill_manifest=[
##                 'HumanHT-12_V3_0_R2_11283641_A.txt',
##                 'HumanHT-12_V4_0_R2_15002873_B.txt',
##                 'HumanHT-12_V3_0_R3_11283641_A.txt',
##                 'HumanHT-12_V4_0_R1_15002873_B.txt',
##                 'HumanMI_V1_R2_XS0000122-MAP.txt',
##                 'HumanMI_V2_R0_XS0000124-MAP.txt',
##                 'HumanRef-8_V2_0_R4_11223162_A.txt',
##                 'HumanRef-8_V3_0_R1_11282963_A_WGDASL.txt',
##                 'HumanRef-8_V3_0_R2_11282963_A.txt',
##                 'HumanRef-8_V3_0_R3_11282963_A.txt',
##                 'HumanWG-6_V2_0_R4_11223189_A.txt',
##                 'HumanWG-6_V3_0_R2_11282955_A.txt',
##                 'HumanWG-6_V3_0_R3_11282955_A.txt',
##                 'MouseMI_V1_R2_XS0000127-MAP.txt',
##                 'MouseMI_V2_R0_XS0000129-MAP.txt',
##                 'MouseRef-8_V1_1_R4_11234312_A.txt',
##                 'MouseRef-8_V2_0_R2_11278551_A.txt',
##                 'MouseRef-8_V2_0_R3_11278551_A.txt',
##                 'MouseWG-6_V1_1_R4_11234304_A.txt',
##                 'MouseWG-6_V2_0_R2_11278593_A.txt',
##                 'MouseWG-6_V2_0_R3_11278593_A.txt',
##                 'RatRef-12_V1_0_R5_11222119_A.txt'],
##             ill_chip=[
##                 'ilmn_HumanHT_12_V3_0_R3_11283641_A.chip',
##                 'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
##                 'ilmn_HumanRef_8_V2_0_R4_11223162_A.chip',
##                 'ilmn_HumanReF_8_V3_0_R1_11282963_A_WGDASL.chip',
##                 'ilmn_HumanRef_8_V3_0_R3_11282963_A.chip',
##                 'ilmn_HumanWG_6_V2_0_R4_11223189_A.chip',
##                 'ilmn_HumanWG_6_V3_0_R3_11282955_A.chip',
##                 'ilmn_MouseRef_8_V1_1_R4_11234312_A.chip',
##                 'ilmn_MouseRef_8_V2_0_R3_11278551_A.chip',
##                 'ilmn_MouseWG_6_V1_1_R4_11234304_A.chip',
##                 'ilmn_MouseWG_6_V2_0_R3_11278593_A.chip',
##                 'ilmn_RatRef_12_V1_0_R5_11222119_A.chip'],
##             ill_bg_mode=['false', 'true'],
##             ill_coll_mode=['none', 'max', 'median'],
##             ill_clm=ANYATOM,
##             ill_custom_chip=ANYATOM,
##             ill_custom_manifest=ANYATOM)),
##     Module(
##         "get_illumina_signal",
##         antecedent("illu_folder", preprocess='illumina'),
##         consequent(
##             "signal_file", preprocess='illumina', logged="no",
##             missing='unknown', format='gct')),
##     Module(
##         "get_illumina_control",
##         antecedent("illu_folder", preprocess='illumina'),
##         consequent(
##             "control_file", preprocess='illumina',
##             logged="no", missing='unknown', format='gct')),
    
##     #-----------------------------------------------------------------------
##     # gpr_files
##     Module(
##         "extract_gpr_files",
##         antecedent("expression_files"),
##         consequent("gpr_files")),
##     Module(
##         "normalize_with_loess",
##         antecedent("gpr_files"),
##         consequent(
##             "signal_file", format="tdf", logged="no",
##             preprocess="loess", missing="unknown")),
    
##     #-----------------------------------------------------------------------
##     Module(
##         "convert_signal_to_tdf",
##         antecedent(
##             "signal_file",
##             format=['pcl', 'res', 'gct', 'jeffs', 'unknown', 'xls']),
##         consequent(
##             "signal_file", format='tdf')),
##     Module(
##         "log_signal",
##         antecedent(
##             "signal_file", logged=[NOVALUE, 'no', 'unknown'], format='tdf'),
##         consequent(
##             "signal_file", logged="yes", format='tdf')),
##     Module(
##         "filter_genes_by_missing_values",
##         antecedent(
##             "signal_file", format='tdf', logged="yes",
##             missing=["yes", "unknown"], filter=[NOVALUE, "no"]),
##         consequent(
##             "signal_file", format='tdf', logged="yes",
##             missing=["yes", "unknown"], filter=ANYATOM)),
##     Module(
##         "fill_missing_with_median",
##         antecedent(
##             "signal_file", format='tdf', logged="yes",
##             missing=["yes", "unknown", NOVALUE]),
##         consequent(
##             "signal_file", format='tdf', logged="yes", missing="median")),
##     Module(
##         "fill_missing_with_zeros",
##         antecedent(
##             "signal_file", format='tdf', logged="yes",
##             missing=["yes", "unknown", NOVALUE]),
##         consequent("signal_file", format='tdf', logged="yes", missing="zero")),
##     Module(
##         "filter_and_threshold_genes",
##         antecedent(
##             "signal_file", logged=[NOVALUE, "no"], format='tdf',
##             predataset=[NOVALUE, 'no']),
##         consequent(
##             "signal_file", logged=[NOVALUE, "no"], format='tdf',
##             predataset="yes")),
##     Module(   # require a rename_list_file
##         "relabel_samples",
##         antecedent(
##             "signal_file", format='tdf', rename_sample=[NOVALUE, "no"]),
##         consequent(
##             "signal_file",  format='tdf', rename_sample="yes")),
##     #------------------------------------------------------------------
##     Module(
##         "normalize_samples_with_quantile",
##         antecedent(
##             "signal_file", quantile_norm=[NOVALUE, "no"], format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"]),
##         consequent(
##             "signal_file", quantile_norm="yes", format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"])),
##     Module(
##         "normalize_samples_with_combat",  # require class label file
##         antecedent(
##             "signal_file", combat_norm=[NOVALUE, "no"], format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"]),
##         consequent(
##             "signal_file", combat_norm="yes", format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"])),
##     Module(
##         "normalize_samples_with_dwd",  # require class label file
##         antecedent(
##             "signal_file", dwd_norm=[NOVALUE, "no"], format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"]),
##         consequent(
##             "signal_file", dwd_norm="yes", format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"])),
##     Module(
##         "normalize_samples_with_bfrm",
##         antecedent(
##             "signal_file", bfrm_norm=[NOVALUE, "no"], format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"]),
##         consequent(
##             "signal_file", bfrm_norm="yes", format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"])),
##     Module(
##         "normalize_samples_with_shiftscale",  # require class label file
##         antecedent(
##             "signal_file", shiftscale_norm=[NOVALUE, "no"], format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"]),
##         consequent(
##             "signal_file", shiftscale_norm="yes", format='tdf',
##             logged="yes", missing=[NOVALUE, "no", "median", "zero"])),
    
##     #------------------------------------------------------------------
## ##    Module(   # may cause loop in the network
## ##        "convert_signal_to_pcl",
## ##        antecedent(
## ##            "signal_file", logged="yes",
## ##            format='tdf', missing=[None, "no", "median", "zero"]),
## ##        consequent(
## ##            "signal_file",logged="yes", format='pcl',
## ##            missing=[None, "no", "median", "zero"])),
##     Module(
##         "gene_center",
##         antecedent(
##             "signal_file", gene_center=[NOVALUE, "no"], logged="yes",
##             format='tdf', missing=[NOVALUE, "no", "median", "zero"]),
##         consequent(
##             "signal_file", gene_center=["mean", "median"], logged="yes",
##             format='tdf', missing=[NOVALUE, "no", "median", "zero"])),
##     Module(
##         "gene_normalize",
##         antecedent(
##             "signal_file", gene_normalize=[NOVALUE, "no"], logged="yes",
##             format='tdf', missing=[NOVALUE, "no", "median", "zero"]),
##         consequent(
##             "signal_file", gene_normalize=["variance", "sum_of_squares"],
##             logged="yes", format='tdf',
##             missing=[NOVALUE, "no", "median", "zero"])),
    
##     #------------------------------------------------------------------
##     Module(  # require class_label_file
##         "filter_genes_by_fold_change_across_classes",
##         antecedent(
##             "signal_file", group_fc=[NOVALUE, "no"], logged="yes",
##             format='tdf', missing=[NOVALUE, "no", "median", "zero"]),
##         consequent(
##             "signal_file", group_fc=ANYATOM, logged="yes", format='tdf',
##             missing=[NOVALUE, "no", "median", "zero"])),
##     Module(  # require class_label_file,generate gene_list_file
##              # and need reorder_genes
##         "rank_genes_by_sample_ttest",
##         antecedent(
##             "signal_file", logged="yes",
##             format='tdf', missing=[NOVALUE, "no", "median", "zero"],
##             gene_order=[NOVALUE, 'no']),
##         consequent(
##             "signal_file", logged="yes", format='tdf',
##             missing=[NOVALUE, "no", "median", "zero"],
##             gene_order=["t_test_p", "t_test_fdr"])),
##     Module(  # require class_label_file,generate gene_list_file
##              # and need reorder_genes
##         "rank_genes_by_class_neighbors",
##         antecedent(
##             "signal_file", logged="yes", format='tdf',
##             missing=[NOVALUE, "no", "median", "zero"],
##             gene_order=[NOVALUE, 'no']),
##         consequent(
##             "signal_file", logged="yes", format='tdf',
##             missing=[NOVALUE, "no", "median", "zero"],
##             gene_order='class_neighbors')),
##     Module(
##         "reorder_genes",   # require gene_list_file
##         antecedent(
##             "signal_file", format='tdf', gene_order=[NOVALUE, 'no']),
##         consequent(
##             "signal_file",  format='tdf', gene_order=['gene_list'])),
##     Module(
##         'annotate_probes',
##         antecedent(
##             "signal_file", format='tdf', annotate=[NOVALUE, "no"]),
##         consequent(
##             "signal_file", format='tdf', annotate="yes")),
##     Module(
##         'remove_duplicate_genes',
##         antecedent(
##             "signal_file", format='tdf', 
##             unique_genes=[NOVALUE, 'no'], annotate='yes'),
##         consequent(
##             "signal_file", format='tdf',
##             unique_genes=['average_genes', 'high_var', 'first_gene'],
##             annotate='yes')),
##     Module(
##         'select_first_n_genes',
##         antecedent(
##             "signal_file", format='tdf', num_features='all'),
##         consequent(
##             "signal_file",  format='tdf', num_features=ANYATOM)),
##     Module(
##         'add_crossplatform_probeid',
##         antecedent(
##             "signal_file", format='tdf', platform='unknown',
##             duplicate_probe=[NOVALUE, "no"]),
##         consequent(
##             "signal_file", format='tdf', platform=ANYATOM,
##             duplicate_probe='yes')),
##     Module(
##         'remove_duplicate_probes',
##         antecedent(
##             "signal_file", format='tdf', duplicate_probe='yes'),
##         consequent(
##             "signal_file",  format='tdf', duplicate_probe='high_var_probe')),
##     Module(
##         'select_probe_by_best_match',
##         antecedent(
##             "signal_file", format='tdf', duplicate_probe='yes'),
##         consequent(
##             "signal_file", format='tdf', duplicate_probe='closest_probe')),
##     ##    Module(#this may cause loop
## ##        'convert_signal_to_gct',
## ##        antecedent(
## ##            "signal_file",
## ##            format='tdf', missing=[NOVALUE, "no", "median", "zero"]),
## ##        consequent(
## ##            "signal_file",
## ##            format='gct',missing=[NOVALUE, "no", "median", "zero"]
## ##            )),
## ##    Module(#this may cause loop
## ##        'unlog_signal',
## ##        antecedent(
## ##            "signal_file",
## ##            format='tdf', logged='yes',
## ##            missing=[NOVALUE, "no", "median", "zero"]),
## ##        consequent(
## ##            "signal_file",
## ##            format='tdf',missing=[NOVALUE, "no", "median", "zero"],
## ##            logged='no')),
##     ]


if __name__ == '__main__':
    test_bie()
