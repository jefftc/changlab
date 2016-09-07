"""

Glossary:
data              A unit of information.  Typically a file made by a program.
attribute         Describes something about the data.  e.g. logged="yes".
input attribute   An attribute of an input data.
output attribute  An attribute of an output data.
user attribute    An attribute provided by the user that can specify
                  preferences
data type         Describes a set of data that share the same attributes.

module            Takes one or more data objects as input and produces a single
                  data object as output.
                  Modules "produce" data.
                  Modules can "convert" data from one type to another.
option            Used in the modules, but does not affect the inferencing.


Classes:
AttributeDef      Describes something about the Data.
Attribute
OptionDef         Something that modifies how a Module works.
Option
DataType
DataNode
Constraint
Consequence
DefaultAttributesFrom
ModuleNode
ModuleDbSummary
Network

Functions:
make_network

find_paths
find_paths_by_start_ids
prune_paths

get_input_nodes           Show all possible input nodes for the network.
get_input_datatypes       Show all combinations of datatype that can be inputs.
group_nodes_by_datatype

summarize_moduledb
check_moduledb

print_modules
print_network
plot_network_gv
get_node_name

read_network
write_network

debug_print

"""
# Classes:
# Pathway
# PathSubset
#
# Functions:
# _init_network
# _split_network
# _complete_network
#
# _prune_by_custom_attributes     Should organize these into classes.
# _prune_alternate_attributes1
# _list_alternate_attributes1
# _find_alternate_attributes
# _prune_alternate_attributes2
# _list_alternate_attributes2
# _prune_superset_pipelines
# _is_superset_pipeline
# _does_path_go_through
# _prune_parallel_pipelines
# _is_parallel_pipeline1
# _is_parallel_pipeline2
# _is_parallel_pipeline3
# _compare_paths
# _build_subpath
#
#     InData(s) -> Module -> OutData
# _bc_to_modules       Module <- OutData
# _bc_to_inputs        InDatas <- Module <- OutData            INFERENCE
# _bc_to_one_input     InDatas[i] <- Module <- OutData         INFERENCE
#                      DO NOT CALL.  Helper for _bc_to_inputs.
# _bc_to_input_ids     InDatas IDs <- Module ID <- OutData ID
# _bc_to_input_and_module_ids
# _fc_to_outputs       InDatas -> Module -> OutData            INFERENCE
# _fc_to_output_ids    InDatas IDs -> Module ID -> OutData ID
# _resolve_constraint
#
# _is_valid_inputs          InDatas (?) -> Module ID -> OutData
# _is_valid_input_i         InData_i (?) -> Module ID -> OutData
# _is_valid_input_ids       InData IDs (?) -> Module ID -> OutData ID
# _is_valid_output          Module -> OutData (?)              INFERENCE
# _is_valid_outdata_id_path  path -> OutData ID (?)
# _is_valid_outmodule_id_path  path -> OutData ID (?)
# _is_valid_output_from_input_and_module_ids
#
# _find_same_data          Find a DataNode that is an exact match.
# _find_compat_data        Find a DataNode that is compatible.
#                          WAS _find_start_node
# _score_same_data
# _score_compat_data       WAS _score_start_nodes.
# _is_data_same
# _is_data_compatible
# _is_attribute_same
# _is_attribute_compatible
# _merge_data_nodes
# _merge_attribute_values
# _is_data_node_atomic
# _make_data_node_atomic
# _is_network_atomic
# _get_atomic_data_node_from_pathway
#
# _merge_start_ids         Should only be called by find_paths_by_start_ids
# _does_based_on_data_conflict_with_out_data    Only by find_paths_by_start_ids
#
# _get_attribute_type
# _assign_case_by_type
#
# _get_parents_of          WAS _backchain_to_ids
# _get_children_of
# _make_parents_dict       WAS _make_backchain_dict
# _make_ancestor_dict
# _make_descendent_dict
# _can_reach_by_bc
# _can_reach_by_fc
#
# _get_custom_values
#
# _iter_upper_diag
# _intersect
# _is_subset
# _flatten
# _flatten1_intlist
# _uniq
# _uniq_intlist
# _uniq_flatten1_intlist
# _intlist2bits
#
# _print_nothing
# _print_string
# _print_line
# _pretty_attributes
#
# _fix_node_id_pairs_after_merge
# _fix_node_id_dict_after_merge
# _fix_ancestor_dict_after_merge
# _product_network
# _product_and_chain
#
# _object_to_dict           Use for writing and reading json file
# _dict_to_object           Use for writing and reading json file


# Rules:
# - A single DataNode cannot go into a ModuleNode as different inputs
#   simultaneously.  E.g.
#   VCFFolder.caller=mutect   -> call_variants_all
#   VCFFolder.caller=varscan  ->
#   The VCFFolders must be separate DataNodes.



# Problem in the way custom attributes are currently handled.  They
# are now provided as a list of Attribute objects.  There is no way to
# distinguish attributes for multiple objects of the same data type.
# E.g.
# BamFolder
#   aligner=star
#   sorted=coordinate
#   indexed=yes
# BamFolder
#   aligner=bwa_mem
# Not possible to indicate that the start aligned BamFolder is sorted.


DEBUG_BACKCHAIN = False
DEBUG_PRUNE_PATHS = False
DEBUG_FIND_PATHS = False
DEBUG_PRUNE_CUSTOM_ATTRIBUTES = False



# These are also used in cbie3module.c.  If these values change, you
# need to change that file too!
TYPE_ATOM = 100
TYPE_ENUM = 101

# Constraints
MUST_BE = 200
CAN_BE_ANY_OF = 201
SAME_AS = 202

# Consequences
SET_TO = 300
SET_TO_ONE_OF = 301
BASED_ON_DATA = 302
SAME_AS_CONSTRAINT = 303

# CONSTRAINT
# behavior       arg1               input_index
# MUST_BE        value              <optional> (set to 0 by default)
# CAN_BE_ANY_OF  list of values     <optional> (set to 0 by default)
# SAME_AS        index of datatype  index of datatype
#
# input_index is the index of the DataType that this constraint
# applies to.  So for SAME_AS, that means data[input_index] should get
# its value from data[arg1].

# CONSEQUENCE
# behavior            arg1                 arg2
# SET_TO              <string>             None
# SET_TO_ONE_OF       <list>               None
# BASED_ON_DATA       <list>               None
# SAME_AS_CONSTRAINT  index of input data  None
# o arg2 is not used?  Can get rid of it.

CONST2STR = {
    TYPE_ATOM: "TYPE_ATOM",
    TYPE_ENUM: "TYPE_ENUM",
    MUST_BE: "MUST_BE",
    CAN_BE_ANY_OF: "CAN_BE_ANY_OF",
    SAME_AS: "SAME_AS",
    SET_TO: "SET_TO",
    SET_TO_ONE_OF: "SET_TO_ONE_OF",
    BASED_ON_DATA: "BASED_ON_DATA",
    SAME_AS_CONSTRAINT: "SAME_AS_CONSTRAINT",
}


# When backchaining, should we allow the attributes of the input data
# to be all possible values, or fix it to the default?  All possible
# values is correct, but generates a combinatorial explosion that is
# difficult to manage.
DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES = False
#DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES = True

MAX_NETWORK_SIZE = 1024 * 4


class AttributeDef:
    def __init__(self, name, values, default_in, default_out, help=None):
        # Make sure name and values are valid.
        assert type(name) is type("")
        assert type(values) is type([]), "Value must be list: %s" % \
               type(values)
        for x in values:
            assert type(x) is type("")

        # Make sure no duplicated values.
        seen = {}
        for x in values:
            assert x not in seen, "Duplicated value (%s) in %s." % (x, name)
            seen[x] = 1

        # Make sure default_in and default_out are valid values.
        assert type(default_in) is type(""), "default_in must be ATOM"
        assert type(default_out) is type(""), "default_out must be ATOM"
        assert default_in in values, \
               "Invalid value %r for attribute %r." % (default_in, name)
        assert default_out in values, \
               "Invalid value %r for attribute %r." % (default_out, name)

        self.name = name
        self.values = values
        self.default_in = default_in
        self.default_out = default_out
        self.help = help

    def is_valid_value(self, value):
        if type(value) is type(""):
            return value in self.values
        elif type(value) is type([]):
            return _is_subset(value, self.values)
        raise AssertionError

    def __cmp__(self, other):
        if not isinstance(other, AttributeDef):
            return cmp(id(self), id(other))
        x1 = [self.name, self.values, self.default_in, self.default_out,
              self.help]
        x2 = [other.name, other.values, other.default_in, other.default_out,
              self.help]
        return cmp(x1, x2)

    def __hash__(self):
        x = self.name, tuple(self.values), self.default_in, self.default_out, \
            self.help
        return hash(x)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x = [repr(self.name),
             repr(self.values),
             repr(self.default_in),
             repr(self.default_out), ]
        if self.help is not None:
            x.append("help=%r" % self.help)
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x

    @staticmethod
    def __init_from_dict(args):
        #inst = AttributeDef(**args)
        inst = AttributeDef(
            args["name"], args["values"], args["default_in"],
            args["default_out"], help=args.get("help"))
        return inst


class Attribute:
    def __init__(self, datatype, name, value):
        assert isinstance(datatype, DataType)
        assert type(name) is type("")
        assert type(value) is type("")

        # Check if this is a valid attribute name for the datatype.
        x = [x for x in datatype.attribute_defs if x == name]
        #x = [x for x in datatype.attribute_defs if x.name == name]
        assert len(x) == 1, "datatype %r does not have attribute %r." % (
            datatype.name, name)
        attr = datatype.attribute_defs[x[0]]
        assert datatype.is_valid_attribute_value(name, value), \
               "Invalid value %r for attribute %r.  Valid values are: %s" % (
            value, name, ", ".join(attr.values))

        self.datatype = datatype
        self.name = name
        self.value = value

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x = [self.datatype.name, repr(self.name), repr(self.value)]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x


class CustomAttributes:
    def __init__(self, attributes, all_nodes):
        # attributes is a list of Attribute objects.
        # all_nodes must be True or False
        assert attributes, "Missing: %s" % repr(attributes)
        assert all_nodes in [True, False], "Invalid all_nodes: %s" % all_nodes
        datatype = attributes[0].datatype
        for i in range(1, len(attributes)):
            assert attributes[i].datatype == datatype
        self.datatype = datatype
        self.attributes = attributes[:]
        self.all_nodes = all_nodes
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [repr(self.datatype), repr(self.attributes), str(self.all_nodes)]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x


class OptionDef:
    def __init__(self, name, default=None, help=None):
        assert type(name) is type("")
        self.name = name
        self.default = default
        self.help = help

    def __cmp__(self, other):
        if not isinstance(other, OptionDef):
            return cmp(id(self), id(other))
        x1 = [self.name, self.default, self.help]
        x2 = [other.name, other.default, self.help]
        return cmp(x1, x2)

    def __hash__(self):
        x = self.name, self.default, self.help
        return hash(x)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x = [repr(self.name), ]
        if self.default is not None:
            x.append("default=%r" % self.default)
        if self.help is not None:
            x.append("help=%r" % self.help)
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x

    @staticmethod
    def __init_from_dict(args):
        assert 'name' in args
        assert 'default' in args
        assert 'help' in args
        inst = OptionDef(
            args['name'], default=args['default'], help=args['help'])
        return inst


class Option:
    def __init__(self, module, name, value):
        assert isinstance(module, ModuleNode)
        assert type(name) is type("")
        assert type(value) is type("")
        self.module = module
        self.name = name
        self.value = value


class Constraint(object):
    def __init__(self, name, behavior, arg1=None, input_index=None):
        if behavior == MUST_BE:
            # name   Name of attribute.
            # arg1   Value of the attribute.
            assert type(arg1) is type(""), \
                   "argument for MUST_BE should be string"
        elif behavior == CAN_BE_ANY_OF:
            # name   Name of attribute.
            # arg1   List of values of the attribute.
            assert type(arg1) in [type([]), type(())]
            for x in arg1:
                assert type(x) is type("")
        elif behavior == SAME_AS:
            # name   Name of attribute.
            # arg1   Index of the datatype that this must match.
            assert type(arg1) is type(0), (
                "arg1 should be the index of the datatype with the "
                "same attribute")
            assert input_index is not None, (
                "input_index must be given for SAME_AS constraint"
            )
            assert type(input_index) is type(0)
            if input_index is not None:
                assert arg1 != input_index
        else:
            raise AssertionError, "Invalid behavior (%s) for constraint %s." %\
                  (behavior, name)
        assert input_index is None or type(input_index) is type(0)

        if behavior == CAN_BE_ANY_OF and len(arg1) == 1:
            behavior = MUST_BE
            arg1 = arg1[0]

        self.name = name
        self.behavior = behavior
        self.arg1 = arg1
        self.input_index = input_index or 0

    def __cmp__(self, other):
        if not isinstance(other, Constraint):
            return cmp(id(self), id(other))
        x1 = [self.name, self.behavior, self.arg1, self.input_index]
        x2 = [other.name, other.behavior, other.arg1, other.input_index]
        return cmp(x1, x2)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x = [repr(self.name), CONST2STR[self.behavior], ]
        if self.arg1 is not None:
            x.append(repr(self.arg1))
        if self.input_index is not None:
            x = x + ["input_index=%s" % self.input_index]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x

    @staticmethod
    def __init_from_dict(args):
        assert 'name' in args
        assert 'behavior' in args
        assert 'arg1' in args
        assert 'input_index' in args
        inst = Constraint(
            args['name'], args['behavior'], arg1=args['arg1'],
            input_index=args['input_index'])
        return inst


class Consequence(object):
    def __init__(self, name, behavior,
                 arg1=None,
                 arg2=None,
                 side_effect=False):
        if behavior == SET_TO:
            assert type(arg1) is type("")
            assert arg2 is None
        elif behavior in [SET_TO_ONE_OF, BASED_ON_DATA]:
            assert type(arg1) in [type([]), type(())], "arg should be list"
            for x in arg1:
                assert type(x) is type("")
            assert arg2 is None
        elif behavior == SAME_AS_CONSTRAINT:
            if arg1 is None:  # default to datatype 0.
                arg1 = 0
            assert type(arg1) is type(0), (
                "Argument to SAME_AS_CONSTRAINT should be the index "
                "of the input variable.")
            assert arg2 is None
        else:
            raise AssertionError, "Invalid consequence: %s" % behavior
        self.name = name
        self.behavior = behavior
        self.arg1 = arg1
        self.arg2 = arg2
        self.side_effect = side_effect

    def __cmp__(self, other):
        if not isinstance(other, Consequence):
            return cmp(id(self), id(other))
        x1 = [self.name, self.behavior, self.arg1, self.arg2, self.side_effect]
        x2 = [other.name, other.behavior, other.arg1, other.arg2,
              other.side_effect]
        return cmp(x1, x2)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x = [repr(self.name), CONST2STR[self.behavior], ]
        if self.arg1 is not None:
            x.append(repr(self.arg1))
        if self.arg2 is not None:
            assert self.arg1 is not None
            x.append(repr(self.arg2))
        if self.side_effect:
            x = x + ["side_effect=True"]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x

    @staticmethod
    def __init_from_dict(args):
        assert 'name' in args
        assert 'behavior' in args
        assert 'arg1' in args
        assert 'arg2' in args
        assert 'side_effect' in args
        inst = Consequence(
            args['name'], args['behavior'],
            arg1=args['arg1'],
            arg2=args['arg2'],
            side_effect=args['side_effect'])
        return inst


class DefaultAttributesFrom(object):
    def __init__(self, input_index):
        assert type(input_index) is type(0)
        self.input_index = input_index

    def __cmp__(self, other):
        if not isinstance(other, DefaultAttributesFrom):
            return cmp(id(self), id(other))
        x1 = [self.input_index]
        x2 = [other.input_index]
        return cmp(x1, x2)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x = [str(self.input_index), ]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x

    @staticmethod
    def __init_from_dict(args):
        assert 'input_index' in args
        inst = DefaultAttributesFrom(args['input_index'])
        return inst


class DataType:
    def __init__(self, name, *attribute_defs, **keywds):
        for x in attribute_defs:
            assert isinstance(x, AttributeDef), repr(x)
        for x in keywds:
            assert x in ["help"]

        # Optimizations:
        # 1.  Save attribute_defs as a dictionary for fast lookups.
        # 2.  Pre-hash objects for fast comparisons.
        attr_defs_dict = {}  # name -> AttributeDef
        for adef in attribute_defs:
            assert adef.name not in attr_defs_dict, \
                   "Multiple attributes named %s" % adef.name
            attr_defs_dict[adef.name] = adef

        x = [attr_defs_dict[x] for x in sorted(attr_defs_dict)]
        attr_defs_tuple = tuple(x)

        self.name = name
        #self.attribute_defs = attribute_defs   # AttributeDef
        self.attribute_defs = attr_defs_dict
        self._attribute_names = sorted(attr_defs_dict)   # optimize
        self.help = keywds.get("help")
        self.hash_ = hash((name, hash(attr_defs_tuple), self.help))

    def get_attribute_def(self, name):
        #x = [x for x in self.attribute_defs if x.name == name]
        #assert len(x) > 0, "DataType %s has no attribute %s." % (
        #    repr(self.name), repr(name))
        #assert len(x) == 1, "Multiple attributes with same name?"
        #return x[0]
        if name not in self.attribute_defs:
            raise KeyError, "DataType %s has no attribute %s." % (
                repr(self.name), repr(name))
        return self.attribute_defs[name]

    def get_attribute_names(self):
        #return [x.name for x in self.attribute_defs]
        #return sorted(self.attribute_defs)
        return self._attribute_names

    def is_valid_attribute_name(self, name):
        return name in self.get_attribute_names()

    def is_valid_attribute_value(self, name, value):
        attr = self.get_attribute_def(name)
        return attr.is_valid_value(value)

    def assert_valid_attribute_dict(self, attr_dict):
        # attr_dict is dictionary of name -> value.  Check if
        # everything in this dictionary is valid.

        # This function is called frequently, so inline many of the
        # function calls for speed.
        x = self.get_attribute_names()
        all_names = {}.fromkeys(x)

        for name, value in attr_dict.iteritems():
            assert name in all_names, \
                   "'%s' is not a known attribute for datatype %s." % (
                name, self.name)

            if name not in self.attribute_defs:
                raise KeyError, "DataType %s has no attribute %s." % (
                    repr(self.name), repr(name))
            #attr = self.attribute_defs[name]
            # Optimization for:
            #assert attr.is_valid_value(value), \
            #       "In a %s, '%s' is not a valid value for '%s'." % (
            #    self.name, value, name)
            # Makes code more complicated for a minor speedup.
            #attr_values_dict = attr.values_dict
            #is_valid = True
            #if type(value) is type(""):
            #    is_valid = value in attr_values_dict
            #elif type(value) is type([]):
            #    for x in value:
            #        if x not in attr_values_dict:
            #            is_valid = False
            #            break
            #assert is_valid, \
            #       "In a %s, '%s' is not a valid value for '%s'." % (
            #    self.name, value, name)

    def __cmp__(self, other):
        if not isinstance(other, DataType):
            return cmp(id(self), id(other))
        # Bug: should compare attributes without regard to order.
        #x1 = [self.name, self.attribute_defs, self.help]
        #x2 = [other.name, other.attribute_defs, other.help]
        #return cmp(x1, x2)
        return cmp(self.hash_, other.hash_)

    def __hash__(self):
        #x = self.name, tuple(self.attribute_defs), self.help
        #return hash(x)
        return self.hash_

    def _resolve_attributes(self, attribute_dict, is_input):
        # Make a dictionary of all the attributes.  The values given
        # by the caller take precedence.  Anything else should be set
        # to the default attributes.
        attrdict = {}
        # Priority 1: Set to the attribute dict.
        for (name, value) in attribute_dict.iteritems():
            attrdict[name] = value
        # Priority 2: Set to default attributes.
        for attr in self.attribute_defs.itervalues():
            if attr.name in attrdict:
                continue
            value = attr.default_in
            if not is_input:
                value = attr.default_out
            attrdict[attr.name] = value
        return attrdict

    def input(self, **attribute_dict):
        # Create a DataNode object.
        attrdict = self._resolve_attributes(attribute_dict, True)
        # Don't bother checking the attributes here.  The DataNode
        # object will do that.
        return DataNode(self, **attrdict)

    def output(self, **attribute_dict):
        attrdict = self._resolve_attributes(attribute_dict, False)
        return DataNode(self, **attrdict)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x = [self.name]
        x += [repr(x) for x in self.attribute_defs.itervalues()]
        if self.help:
            x.append("help=%r" % self.help)
        return "DataType(%s)" % ", ".join(x)

    @staticmethod
    def __init_from_dict(args):
        assert 'name' in args
        assert 'attribute_defs' in args
        assert 'help' in args
        #inst = DataType(
        #    args['name'], *args['attribute_defs'], help=args['help'])
        dictionary = args['attribute_defs']
        attributes = []
        for i in dictionary:
            attributes.append(dictionary[i])
        inst = DataType(args['name'], *attributes, **{"help" : args['help']})
        return inst


class DataNode(object):
    # Members:
    # datatype      Datatype object.
    # attributes    Dict of attribute name -> value or list of values

    # Should not be called by the user.  Should always be created from
    # a DataType object.
    def __init__(self, datatype, **keywds):
        # keywds is a dictionary of attribute name ->
        # value (or list of values).

        # Make sure values are provided for every attribute.
        attr_names = datatype.get_attribute_names()
        for name in attr_names:
            assert name in keywds, "No value given for %s." % name

        # Make sure the values of the attributes are legal.
        #for name, value in keywds.iteritems():
        #    assert datatype.is_valid_attribute_name(name), \
        #           "'%s' is not a known attribute for datatype %s." % (
        #        name, datatype.name)
        #    assert datatype.is_valid_attribute_value(name, value), \
        #           "In a %s, '%s' is not a valid value for '%s'." % (
        #        datatype.name, value, name)
        datatype.assert_valid_attribute_dict(keywds)

        self.datatype = datatype
        self.attributes = keywds.copy()

    def __cmp__(self, other):
        if not isinstance(other, DataNode):
            return cmp(id(self), id(other))
        x1 = [self.datatype, self.attributes]
        x2 = [other.datatype, other.attributes]
        return cmp(x1, x2)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        #keywds = self.attributes.copy()
        x = [self.datatype.name,
             #_pretty_attributes(keywds),
             _pretty_attributes(self.attributes), ]
        x = [x for x in x if x]
        return "DataNode(%s)" % ", ".join(x)

    @staticmethod
    def __init_from_dict(args):
        assert 'datatype' in args
        assert 'attributes' in args
        return DataNode(args['datatype'], **args['attributes'])


# DataNode + identifier.
class IdentifiedDataNode:
    def __init__(self, data, identifier=""):
        self.data = data
        self.identifier = identifier
    def __repr__(self):
        x = str(self.data) + ' identifier:' + self.identifier
        return x


def _pretty_node_id(network, node_id, attributes_to_show=[]):
    x = network.nodes[node_id]
    x = _pretty_node(x, attributes_to_show=attributes_to_show)
    x = "[%s] %s" % (node_id, x)
    return x


def _pretty_node(node, attributes_to_show=[]):
    node_name = get_node_name(node)
    attrs = {}
    if isinstance(node, DataNode):
        attrs = node.attributes
    elif isinstance(node, IdentifiedDataNode):
        attrs = node.data.attributes
    attr_info = []
    for name in attributes_to_show:
        if name not in attrs:
            continue
        x = "%s=%s" % (name, attrs[name])
        attr_info.append(x)
    attr_info = ", ".join(attr_info)
    if attr_info:
        pretty = "%s (%s)" % (node_name, attr_info)
    else:
        pretty = node_name
    return pretty


def get_node_name(node):
    if isinstance(node, IdentifiedDataNode):
        return get_node_name(node.data)
    elif isinstance(node, DataNode):
        return node.datatype.name
    elif isinstance(node, ModuleNode):
        return node.name
    assert type(node) is not type(0), \
           "Cannot get name of node.  Was passed integer."
    raise AssertionError, "Cannot get name of node."


class ModuleNode:
    def __init__(self, name, in_datatypes, out_datatype, *params, **keywds):
        # params is a list of Constraint, Consequence, and OptionDef.
        # objects.
        assert type(name) is type("")
        for k in keywds:
            assert k in ["help"]

        # The caller can provide either a single DataType object or a
        # list of DataType objects.  Make sure it is always a list of
        # objects.  If it is a single object, make it a list.
        assert type(in_datatypes) is type([]) or \
               isinstance(in_datatypes, DataType)
        if type(in_datatypes) != type([]):
            in_datatypes = [in_datatypes]
        for x in in_datatypes:
            assert isinstance(x, DataType)
        assert isinstance(out_datatype, DataType)

        # Separate the param objects.
        constraints = []
        consequences = []
        option_defs = []  # OptionDef
        default_attributes_from = []
        for x in params:
            if isinstance(x, Constraint):
                constraints.append(x)
            elif isinstance(x, Consequence):
                consequences.append(x)
            elif isinstance(x, OptionDef):
                option_defs.append(x)
            elif isinstance(x, DefaultAttributesFrom):
                default_attributes_from.append(x)
            else:
                print type(x)
                raise AssertionError, "invalid parameter: %s" % repr(x)

        # Convenience hack: If the module doesn't convert the data
        # type, and no DefaultAttributesFrom is given, then set the
        # default_attributes_from.
        if not default_attributes_from and \
               len(in_datatypes) == 1 and in_datatypes[0] == out_datatype:
            default_attributes_from = [DefaultAttributesFrom(0)]
        # Another hack: If no DefaultAttributesFrom is given, and only
        # one of the in_datatypes is the same as the out_datatype,
        # then set the default_attributes_from.
        if not default_attributes_from and len(in_datatypes) > 1:
            i = [i for i in range(len(in_datatypes))
                 if in_datatypes[i] == out_datatype]
            if len(i) == 1:
                i = i[0]
                default_attributes_from = [DefaultAttributesFrom(i)]

        # Check default_attributes_from.  Should be a list of
        # DefaultAttributesFrom objects.
        assert len(default_attributes_from) <= len(in_datatypes)
        seen = {}  # make sure no duplicates
        for daf in default_attributes_from:
            assert daf.input_index < len(in_datatypes)
            assert out_datatype.name == in_datatypes[daf.input_index].name
            assert daf.input_index not in seen
            seen[daf.input_index] = 1

        # default_attributes_from can be an empty list if the
        # attributes of the output object should come from none of the
        # input data types.  Probably the most common case if the
        # module converts the data type.

        #assert len(default_attributes_from) <= 1
        #x = None
        #if default_attributes_from:
        #    x = default_attributes_from[0]
        #default_attributes_from = x
        #if default_attributes_from:
        #    assert len(in_datatypes) > 1
        #    assert default_attributes_from.input_index < len(in_datatypes)
        #    assert out_datatype == \
        #           in_datatypes[default_attributes_from.input_index]

        # Any checking necessary on Option?

        self.name = name
        self.in_datatypes = in_datatypes
        self.out_datatype = out_datatype
        self.constraints = constraints
        self.consequences = consequences
        self.default_attributes_from = default_attributes_from
        self.option_defs = option_defs
        self.help = keywds.get("help")
        # To optimize __cmp__.
        self.hash_ = hash((
            name, tuple(in_datatypes), out_datatype,
            tuple(constraints), tuple(consequences),
            tuple(default_attributes_from), tuple(option_defs), self.help))

        for x in constraints:
            self._assert_constraint(
                name, in_datatypes, out_datatype,
                constraints, consequences, x, default_attributes_from)
        for x in consequences:
            self._assert_consequence(name, in_datatypes, out_datatype,
                                     constraints, x)

    def _assert_constraint(self, name, in_datatypes, out_datatype, constraints,
                           consequences, constraint, default_attributes_from):
        # Get the input datatype that this constraint refers to.
        i = constraint.input_index
        assert i < len(in_datatypes), \
               "Invalid constraint index %d in module %s" % (i, name)
        in_datatype = in_datatypes[i]

        assert constraint.behavior in [MUST_BE, CAN_BE_ANY_OF, SAME_AS]
        if constraint.behavior in [MUST_BE, CAN_BE_ANY_OF]:
            assert in_datatype.is_valid_attribute_name(constraint.name), \
                ("%r: Invalid attribute %r for datatype %r." %
                 (name, constraint.name, in_datatype.name))
            assert in_datatype.is_valid_attribute_value(
                constraint.name, constraint.arg1), \
                ("%r: Invalid value %r for attribute %r." %
                 (name, constraint.arg1, constraint.name))
        elif constraint.behavior == SAME_AS:
            # Make sure the datatype has this attribute.
            dt = in_datatypes[constraint.arg1]
            assert dt.is_valid_attribute_name(constraint.name)
            # Make sure value can be resolved.
            assert len(in_datatypes) > 1, (
                "%r: SAME_AS constraint requires at least two input "
                "datatypes." % name)
            const = _resolve_constraint(constraint, constraints)
            assert const.behavior in [MUST_BE, CAN_BE_ANY_OF]
            #assert constraint.arg1 < len(in_datatypes)
            #assert constraint.arg1 != constraint.input_index
            ## Make sure there is a MUST_BE or CAN_BE_ANY_OF constraint
            ## on constraint.arg1.
            #x = constraints
            #x = [x for x in x if x.name == constraint.name]
            #x = [x for x in x if x.input_index == constraint.arg1]
            #assert len(x) > 0, (
            #    "%r: %r SAME_AS %d, but datatype %d has no constraint on %r."%
            #    (name, constraint.name, constraint.arg1, constraint.arg1,
            #     constraint.name))
            #assert len(x) == 1
            #x = x[0]
            #assert x.behavior in [MUST_BE, CAN_BE_ANY_OF]
        else:
            raise NotImplementedError

        # For every constraint, there must be a consequent given.
        # Need to specify what the module does with the variable.
        # Exception: when DefaultAttributesFrom specifies a different
        # input node.
        if in_datatype.name == out_datatype.name:
            assert default_attributes_from
            daf = [x.input_index for x in default_attributes_from]
            if constraint.input_index in daf:
                x = [x for x in consequences if x.name == constraint.name]
                assert x, "%r: constraint but no consequence for %r." % (
                    name, constraint.name)

    def _assert_consequence(self, name, in_datatypes, out_datatype,
                            constraints, consequence):
        import itertools
        assert consequence.name in out_datatype.get_attribute_names(), \
               "ModuleNode %r refers to an unknown attribute %r." % (
            self.name, consequence.name)

        if consequence.behavior in [SET_TO, SET_TO_ONE_OF, BASED_ON_DATA]:
            assert out_datatype.is_valid_attribute_value(
                consequence.name, consequence.arg1), \
                "'%s' is not a valid value for '%s' in module '%s'" % (
                consequence.arg1, consequence.name, name)
        elif consequence.behavior == SAME_AS_CONSTRAINT:
            # Make sure index on input variable is reasonable.
            index = consequence.arg1  # index of input variable
            assert index < len(in_datatypes), \
                   "Invalid input index (%s) for module %s:%s." % (
                index, name, consequence.name)
            in_datatype = in_datatypes[index]

            # Make sure there is a valid constraint.
            x = constraints
            #x = [x for x in constraints
            #     if x.behavior in [MUST_BE, CAN_BE_ANY_OF, SAME_AS]]
            x = [x for x in x if x.input_index == index]
            x = [x for x in x if x.name == consequence.name]
            assert len(x) > 0, (
                "%r: I could not find a constraint on %r for input "
                "datatype %d." % (name, consequence.name, index))
            assert len(x) == 1, (
                "%r: Multiple constraints for %r for input "
                "datatype %d." % (name, consequence.name, index))
            cons = x[0]

            # Make sure the values of this constraint are allowed in
            # the input and output datatypes.  The values of the
            # consequent should be a subset of the values of the
            # constraint.
            if cons.behavior in [MUST_BE, CAN_BE_ANY_OF]:
                in_attr = in_datatype.get_attribute_def(consequence.name)
                out_attr = out_datatype.get_attribute_def(consequence.name)
                assert in_attr.is_valid_value(cons.arg1), \
                       "Invalid value for %s (%s) in module %s" % (
                    in_attr.name, cons.arg1, self.name)
                assert out_attr.is_valid_value(cons.arg1), \
                       "Invalid value for %s (%s) in module %s" % (
                    out_attr.name, cons.arg1, self.name)
        else:
            raise AssertionError, consequence.behavior

        # Make sure a MUST_BE constraint and SET_TO constraint aren't
        # the same value.  If so, this will confuse the backchainer
        # into thinking this is the objective of the module, even
        # if it isn't.
        if consequence.behavior == SET_TO:
            for x in itertools.product(in_datatypes, constraints):
                in_datatype, constraint = x
                if in_datatype.name != out_datatype.name:
                    continue
                if consequence.name != constraint.name:
                    continue
                if constraint.behavior != MUST_BE:
                    continue
                assert constraint.arg1 != consequence.arg1, \
                       "MUST_BE and SET_TO should not be set to the " \
                       "same value in module %s (%s)." % (
                    self.name, constraint.name)

    def __hash__(self):
        return self.hash_

    def __cmp__(self, other):
        if not isinstance(other, ModuleNode):
            return cmp(id(self), id(other))
        #x1 = [self.name, self.in_datatypes, self.out_datatype,
        #      self.constraints, self.consequences,
        #      self.default_attributes_from, self.option_defs, self.help]
        #x2 = [other.name, other.in_datatypes, other.out_datatype,
        #      other.constraints, other.consequences,
        #      other.default_attributes_from, other.option_defs, other.help]
        #return cmp(x1, x2)
        return cmp(self.hash_, other.hash_)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x1 = repr(self.name)
        if len(self.in_datatypes) == 1:
            x2 = self.in_datatypes[0].name
        else:
            x2 = [x.name for x in self.in_datatypes]
            x2 = "[%s]" % ", ".join(x2)
        x3 = self.out_datatype.name
        x4 = [repr(x) for x in self.constraints]
        x5 = [repr(x) for x in self.consequences]
        x6 = [repr(x) for x in self.default_attributes_from]
        x7 = [repr(x) for x in self.option_defs]
        x = [x1, x2, x3] + x4 + x5 + x6 + x7
        if self.help is not None:
            x.append("help=%r" % self.help)
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x

    @staticmethod
    def __init_from_dict(args):
        assert 'name' in args
        assert 'in_datatypes' in args
        assert 'out_datatype' in args
        assert 'consequences' in args
        assert 'constraints' in args
        assert 'option_defs' in args
        assert 'default_attributes_from' in args
        assert 'help' in args
        name = args['name']
        in_datatypes = args['in_datatypes']
        out_datatype = args['out_datatype']
        help_ = args['help']
        #params = args.copy()
        #del params['name']
        #del params['in_datatypes']
        #del params['out_datatype']
        #del params['help']
        #params = (params['consequences']+params['constraints']+
        #            params['user_inputs']+params['default_attributes_from'])
        params = (args['consequences'] + args['constraints'] +
                  args['option_defs'] + args['default_attributes_from'])
        inst = ModuleNode(
            name, in_datatypes, out_datatype, *params, **{"help" : help_})
        return inst


class ModuleDbSummary:
    # module_names    List of module names (strings).
    # name2module     Dictionary of module name (string) to ModuleNode
    #                 object.
    # name2datatypes  Dict of module name to (in Datatype, out Datatype).
    # datatypes       List of Datatype objects.
    def __init__(self, module_names, name2module, name2datatypes, datatypes):
        self.module_names = module_names[:]
        self.name2module = name2module.copy()
        self.name2datatypes = name2datatypes.copy()
        self.datatypes = datatypes[:]


class Network:
    def __init__(self, nodes, transitions):
        # nodes should be a list of DataNode or ModuleNode objects.
        # DataNode transition to ModuleNodes, and ModuleNodes
        # transition to DataNode.

        # Make sure nodes are DataNode or ModuleNode objects.
        for n in nodes:
            assert isinstance(n, DataNode) or isinstance(n, ModuleNode)
        # Make sure transitions point to the right types of objects.
        for node_id, next_ids in transitions.iteritems():
            assert node_id >= 0 and node_id < len(nodes)
            for nid in next_ids:
                assert nid >= 0 and nid < len(nodes)
            for nid in next_ids:
                n1 = nodes[node_id]
                n2 = nodes[nid]
                if isinstance(n1, DataNode):
                    assert isinstance(n2, ModuleNode)
                else:
                    assert isinstance(n2, DataNode)

        self.nodes = nodes[:]
        self.transitions = transitions.copy()

    def iterate(self, node_class=None):
        # Yield tuple of (node_id, next_node_ids).  node_class is the
        # class of the node of the network (either DataNode or
        # ModuleNode).  If provided, will only iterate over that kind
        # of nodes.
        assert node_class in [DataNode, ModuleNode]
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

    def delete_nodes(self, node_ids):
        # Make sure no duplicates.
        for i in range(len(node_ids)):
            assert node_ids[i] not in node_ids[i+1:]
        network = Network(self.nodes, self.transitions)
        # Delete from high to low so the node_ids don't get messed up.
        node_ids = reversed(sorted(node_ids))
        for nid in node_ids:
            network = network.delete_node(nid)
        return network

    def merge_nodes(self, node_ids, nodeid2parents=None):
        """node_ids is a list of the indexes of nodes.  Replace all
        these nodes with just a single one.  Returns a new Network
        object."""
        if nodeid2parents is None:
            nodeid2parents = _make_parents_dict(self)
        node_ids = sorted(node_ids)

        # Make sure no duplicate node_ids.
        for i in range(len(node_ids) - 1):
            assert node_ids[i] != node_ids[i + 1], "Duplicate node IDs."

        # Make sure nodes are the same type.
        for i in range(1, len(node_ids)):
            n1 = self.nodes[node_ids[0]]
            n2 = self.nodes[node_ids[i]]
            assert n1.__class__ == n2.__class__, "%s %s" % (
                n1.__class__, n2.__class__)

        # Keep the first node, and delete the rest.

        # Make transitions to any node_ids point to the first one.
        node_id = node_ids[0]
        prev_ids = []
        for nid in node_ids[1:]:
            x = nodeid2parents.get(nid, [])
            prev_ids.extend(x)
        transitions = self.transitions.copy()
        for prev_id in prev_ids:
            x = transitions[prev_id]
            if node_id not in x:
                x = x + [node_id]
            transitions[prev_id] = x

        # Make the first node point to all the next_nodes of the other
        # node_ids.
        nid0 = node_ids[0]
        for nidi in node_ids[1:]:
            x = transitions.get(nid0, []) + transitions.get(nidi, [])
            x = sorted({}.fromkeys(x))
            transitions[nid0] = x

        x = Network(self.nodes, transitions)
        x = x.delete_nodes(node_ids[1:])
        return x

    def __cmp__(self, other):
        if not isinstance(other, Network):
            return cmp(id(self), id(other))
        # Optimization.  Do some quick comparisons first.
        if id(self) == id(other):
            return 0
        x = cmp(len(self.nodes), len(other.nodes))
        if x != 0:
            return x
        x = cmp(self.transitions, other.transitions)
        if x != 0:
            return x
        x1 = [self.nodes, self.transitions]
        x2 = [other.nodes, other.transitions]
        return cmp(x1, x2)

    @staticmethod
    def __init_from_dict(args):
        assert 'nodes' in args
        assert 'transitions' in args
        new_transition = dict()
        for key, value in args['transitions'].items():
            new_transition[int(key)] = value
        inst = Network(args['nodes'], new_transition)
        return inst


def make_network(moduledb, out_data, custom_attributes):
    import copy

    # Clean up this code.
    network = _init_network(moduledb, out_data, custom_attributes)
    #plot_network_gv("network1.png", network, verbose=True)

    # The inferencing engine should generate DataNodes whose
    # attributes are all atomic.
    assert _is_network_atomic(network)

    # Split the data nodes so that everything is TYPE_ATOM.  Fixes
    # problems in the inference, and also makes the other steps easier
    # to handle.  Carefully build the network back up.
    #network = _split_network(network)

    optimizers = [
        # There should not be any cycles.  This is not needed.
        #_OptimizeNoCycles(),
        _OptimizeNoInvalidOutputs(),
        _OptimizeNoDuplicateModules(),
        _OptimizeNoDuplicateData(),
        _OptimizeMergeData1(),
        # Don't do this merging.  See below for reason.
        #_OptimizeMergeData2(),
        ]
    it = 0
    old_network = None
    while old_network != network:
        old_network = copy.deepcopy(network)
        for opt in optimizers:
            #old_network2 = copy.deepcopy(network)
            network = opt.optimize(network, custom_attributes)
            #if old_network2 != network:
            #    print "Optimized with %s." % opt.__class__.__name__
        it += 1

    #plot_network_gv("network2.png", network, verbose=True)
    #num = 0
    #plot_network_gv("test-%02d.png" % num, network, verbose=True); num+=1
    # This makes the network really messy.  Might have to be rewritten.
    #network = _complete_network(network, custom_attributes)

    return network


def _init_network(moduledb, out_data, custom_attributes):
    # Return a Network object.
    check_moduledb(moduledb)
    if isinstance(out_data, DataType):
        out_data = out_data.output()
    assert isinstance(out_data, DataNode)

    #nodes = []  # list of DataNode or ModuleNode objects.
    #transitions = {}  # list of index -> list of indexes

    network = Network([], {})
    network.nodes.append(out_data)
    stack = [0]
    seen = {}
    nit = -1
    while stack:
        nit += 1

        if len(network.nodes) >= MAX_NETWORK_SIZE:
            plot_network_gv("network.png", network, verbose=True)
            write_network("network.json", network)
        assert len(network.nodes) < MAX_NETWORK_SIZE, \
               "network [%d] too large" % len(network.nodes)
        #_print_network(Network(nodes, transitions))

        # Pop the next node off the stack.
        node_id = stack.pop()
        assert node_id < len(network.nodes)

        # If I've already seen this node, then don't process it again.
        if node_id in seen:
            continue
        seen[node_id] = 1

        node = network.nodes[node_id]
        if isinstance(node, DataNode):
            # Backwards chain to the previous module.
            modules = _bc_to_modules(moduledb, node)
            for m in modules:
                network.nodes.append(m)
                m_id = len(network.nodes) - 1
                stack.append(m_id)
                network.transitions[m_id] = network.transitions.get(m_id, [])
                network.transitions[m_id].append(node_id)
        elif isinstance(node, ModuleNode):
            # Networks generated may not be identical.  Nodes can be
            # added onto the stack in different orders, leading to
            # differences in how the network grows.
            all_inputs = {}
            for x in network.transitions[node_id]:
                combos = _bc_to_inputs(network, node_id, x, custom_attributes)
                for combo in combos:
                    for x in combo:
                        all_inputs[x] = 1

            for d in all_inputs:
                d_id = _find_same_data(network.nodes, d)
                if d_id == -1:
                    network.nodes.append(d)
                    d_id = len(network.nodes) - 1
                stack.append(d_id)
                network.transitions[d_id] = network.transitions.get(d_id, [])
                network.transitions[d_id].append(node_id)
        else:
            raise AssertionError, "Unknown node type: %s" % node

        # DEBUG: Check for cycles.
        #plot_network_gv("cycle-%02d.png" % nit, network, verbose=True)
        #_make_ancestor_dict(network)

    # Remove the duplicates from transitions.
    for nid, next_ids in network.transitions.iteritems():
        network.transitions[nid] = _uniq_intlist(next_ids)

    # Check for cycles here.  Will not be able to make ancestor dict
    # if there is a cycle.
    _make_ancestor_dict(network)

    #network = Network(nodes, transitions)
    return network


def _split_network(network):
    # Inferencing can lead to a situation where a ModuleNode points to
    # DataNode that it can't generate.  E.g.
    # trim_adapters -> Fastq.trimmed=["no", "yes"]  (should only be "yes")
    #
    # Solution: split Fastq into multiple objects.
    # _OptimizeNoInvalidOutputs will remove the bad links.
    import itertools

    nodeid2parents = _make_parents_dict(network)

    to_delete = []
    for node_id in range(len(network.nodes)):
        node = network.nodes[node_id]
        if not isinstance(node, DataNode):
            continue
        # Look for attributes with multiple values.  Once found, replace
        # with all possible individual values.
        attr_names = []   # list of attribute names
        attr_values = []  # list of list of attribute values
        for name, value in node.attributes.iteritems():
            if _get_attribute_type(value) != TYPE_ENUM:
                continue
            attr_names.append(name)
            attr_values.append(value)
        if not attr_names:
            continue
        # Make a new DataNode.
        for values in itertools.product(*attr_values):
            attrs = node.attributes.copy()
            assert len(values) == len(attr_names)
            for name, value in zip(attr_names, values):
                attrs[name] = value
            x = DataNode(node.datatype, **attrs)
            network.nodes.append(x)
            nid = len(network.nodes)-1
            # Make sure this points to all the children of the
            # previous node.
            network.transitions[nid] = network.transitions[node_id][:]
            # Make sure all the parent nodes point to this one.
            for pid in nodeid2parents.get(node_id, []):
                network.transitions[pid].append(nid)
        # Mark the old node for deletion.
        to_delete.append(node_id)
    network = network.delete_nodes(to_delete)
    return network


def _complete_network(network, custom_attributes):
    # Sometimes, the network generated by backchaining may be missing
    # some links.  This function will search for missing links and add
    # them back into the network.  Returns a new Network object.
    #
    # Example:
    # 1.  PSF (preprocess=unknown) -> rank_genes_by_class_neighbors ->
    #     GeneListFile
    #     preprocess assigned to unknown because it is the default
    #     value for PSF files.
    # 2.  During inferencing, PSF (preprocess=illumina) is created.
    #     It does not point to rank_genes_by_class_neighbors--it
    #     points to another module.
    # 3.  complete_network will add link:
    #     PSF (preprocess=illumina) -> rank_genes_by_class_neighbors
    #
    # This occurs because of the optimization we made where
    # backchaining created antecedents with default values.  If the
    # antecedents countained all possible values, this would not be
    # necessary.
    import copy
    import itertools

    debug_print(DEBUG_BACKCHAIN, "Completing network.")

    network = copy.deepcopy(network)
    nodeid2parents = _make_parents_dict(network)
    ancestors = _make_ancestor_dict(network)
    descendents = _make_descendent_dict(network)

    # For each DataNode object, check to see if it can be the
    # antecedent of any ModuleNode objects.
    data_ids = [x for x in range(len(network.nodes))
                if isinstance(network.nodes[x], DataNode)]
    module_ids = [x for x in range(len(network.nodes))
                  if isinstance(network.nodes[x], ModuleNode)]
    for x in itertools.product(data_ids, module_ids):
        input_id, module_id = x

        # If data_id already points to module_id, then ignore
        # this.
        if module_id in network.transitions.get(input_id, []):
            continue

        # If this node is not a DataType that the module takes, then
        # don't bother checking.
        found = False
        for dt in network.nodes[module_id].in_datatypes:
            if network.nodes[input_id].datatype.name == dt.name:
                found = True
                break
        if not found:
            continue

        # Don't add a link from data_id to module_id if it would
        # create a cycle.
        if module_id in ancestors[input_id]:
            #debug_print("Skipping DataNode %d -> ModuleNode %d (cycle)." % (
            #    input_id, module_id))
            continue

        # Since modules can take multiple inputs, we need to combine
        # input_id with all previous input IDs and try all possible
        # combinations.
        #x = _get_parents_of(network, module_id)
        x = nodeid2parents.get(module_id, [])
        combined_ids = x + [input_id]

        # Find combinations of inputs that are compatible with the
        # network.
        combos = _bc_to_input_ids(
            network, module_id, custom_attributes, all_input_ids=combined_ids,
            nodeid2parents=nodeid2parents)

        # Add the new transitions.
        added = []
        for id_ in itertools.chain.from_iterable(combos):
            # Probably don't need to search through.  All the id_'s,
            # except for input_id, is already in a parent of this
            # node.
            assert id_ in network.transitions
            if module_id in network.transitions[id_]:
                continue
            # Add id_ -> module_id.
            network.transitions[id_].append(module_id)
            added.append(id_)
            debug_print(
                DEBUG_BACKCHAIN,
                "Completing DataNode %s [%d] -> ModuleNode %s [%d]." % (
                    network.nodes[id_].datatype.name, id_,
                    network.nodes[module_id].name, module_id))

        # id_ is now a parent of module_id.
        for id_ in added:
            if module_id not in nodeid2parents:
                nodeid2parents[module_id] = []
            if id_ not in nodeid2parents[module_id]:
                nodeid2parents[module_id].append(id_)

        # Module and all its descendents inherit the ancestors of all
        # the added nodes (including the added nodes).
        all_ids = [module_id] + descendents.get(module_id, [])
        for node_id in all_ids:
            anc = ancestors[node_id]
            anc.extend(added)
            for id_ in added:
                anc.extend(ancestors.get(id_, []))
            ancestors[node_id] = _uniq_intlist(anc)

        # All the added nodes inherit the descendents of the Module.
        for id_ in added:
            desc = descendents[id_]
            desc.append(module_id)
            desc.extend(descendents.get(module_id, []))
            descendents[id_] = _uniq_intlist(desc)

    return network


class _OptimizeNoCycles:
    # Methods:
    # optimize
    #
    # _find_cycle          Find a cycle in a network.
    # _find_cycle_from_one_node
    # _find_cycle_from_all_nodes
    # _choose_breakpoint   Given a cycle, choose the best breakpoint.
    # _break_cycle         Break a cycle
    def __init__(self):
        pass

    def optimize(self, network, custom_attributes):
        global NUM_TIMES
        # Do backwards chaining.  If I encounter a cycle, then break it.

        # Optimization.  Since many cycles are short, break the
        # 5-cycles first to speed up the search.
        # [131104] Is this still necessary, since the algorithm is now
        # smarter about not searching non-cycle nodes?
        # Length of cycles:
        # 5     168,081   67%
        # 7      84,056   33%
        # 9-17      318    0%

        while True:
            cycle = self._find_cycle(network, 5)
            if not cycle:
                break
            transition = self._choose_breakpoint(network, cycle)
            if transition:
                network = self._break_cycle(network, transition)
        while True:
            cycle = self._find_cycle(network, 0)
            if not cycle:
                break
            transition = self._choose_breakpoint(network, cycle)
            if transition:
                network = self._break_cycle(network, transition)

        return network

    def _list_noncycle_node_ids(self, network, nodeid2parents):
        # Return a list of the node_ids that are not in cycles.

        # The nodes at the top of the tree (no prev nodes) are not in
        # cycles.  The nodes at the bottom of the tree (no next nodes)
        # are not in cycles.
        #
        # If all of a node's next nodes are noncycle, then it is
        # noncycle.
        # If all of a node's prev nodes are noncycle, then it is
        # noncycle.
        noncycle = {}

        while True:
            changed = False
            node_ids = [i for i in range(len(network.nodes))
                        if i not in noncycle]
            for node_id in node_ids:
                prev_ids = nodeid2parents.get(node_id, [])
                next_ids = network.transitions.get(node_id, [])
                x1 = [i for i in prev_ids if i not in noncycle]
                x2 = [i for i in next_ids if i not in noncycle]
                if x1 and x2:
                    continue
                # Either all prev_ids are noncycle (or missing), or
                # all next_are noncycle (or missing).
                noncycle[node_id] = 1
                changed = True
            if not changed:
                break
        return noncycle

    def _find_cycle(self, network, max_path_length):
        assert max_path_length >= 0
        nodeid2parents = _make_parents_dict(network)
        noncycle = self._list_noncycle_node_ids(network, nodeid2parents)

        cycle = None
        node_ids = [i for i in range(len(network.nodes)) if i not in noncycle]
        for start_id in node_ids:
            cycle = self._find_cycle_from_one_node(
                network, start_id, max_path_length, noncycle, nodeid2parents)
            if cycle:
                break
        return cycle

    def _find_cycle_from_one_node(self, network, start_id, max_path_length,
                                  noncycle, nodeid2parents):
        # Do a depth-first search and look for cycles.  Return a cycle
        # (list of node_ids) or None.  The cycle will start and end
        # with the same node_id.
        # Previously did a breadth-first search, but stack.pop(0) was
        # running too slowly (see below).  Depth-first search will be
        # faster to find a cycle, if it exists, anyway.
        if not nodeid2parents:
            nodeid2parents = _make_parents_dict(network)

        cycle = None
        # list of (node_id, path (not including node_id))
        stack = [(start_id, [])]
        while stack:
            # Slow.
            #node_id, path = stack.pop(0)
            # Really slow.
            #node_id, path = stack[0]
            #stack = stack[1:]
            # 10x faster than pop(0).
            node_id, path = stack.pop()
            if node_id in noncycle:
                continue

            if max_path_length and len(path) > max_path_length:
                continue
            if node_id in path:
                # If this node_id is already in the path, then this is
                # a cycle.
                i = path.index(node_id)
                cycle = path[i:] + [node_id]
                break
            # Add node to the path.
            path = path + [node_id]
            for prev_id in nodeid2parents.get(node_id, []):
                stack.append((prev_id, path))
        return cycle

    def _find_depth_of_nodes(self, network):
        # Do a breadth-first search to assign the depth of each node.
        assert network.nodes

        nodeid2parents = _make_parents_dict(network)

        # OPTIMIZE: can memoize this function.
        stack = [(0, -1)]  # node_id, next_depth
        nodeid2depth = {}
        while stack:
            node_id, next_depth = stack.pop(0)
            if node_id in nodeid2depth:
                continue
            depth = next_depth + 1
            nodeid2depth[node_id] = depth
            for prev_id in nodeid2parents.get(node_id, []):
                stack.append((prev_id, depth))
        return nodeid2depth

    def _choose_breakpoint(self, network, cycle):
        # Break the cycle at the point furthest from the root of the
        # network (node 0).  If I choose the wrong breakpoint, can
        # leave a whole section dangling.  Return tuple of (node_id,
        # next_id) of the transition to break.
        nodeid2depth = self._find_depth_of_nodes(network)

        # See if this cycle is already broken from the main network.
        x = [x for x in cycle if x not in nodeid2depth]
        if x:
            return None

        # Furthest point is the one with the highest depth.
        depths = [nodeid2depth[x] for x in cycle]
        schwartz = zip(depths, cycle)
        schwartz.sort()
        highest_depth, highest_node_id = schwartz[-1]

        # Find the transition.  cycle is a list of [id_3, id_2, id_1,
        # id_0], where id_0 points to id_1 (etc).
        transition = None
        for i in range(len(cycle) - 1):
            next_id, node_id = cycle[i], cycle[i + 1]
            if next_id == highest_node_id:
                transition = node_id, next_id
                break
        assert transition is not None
        return transition

    def _break_cycle(self, network, bad_transition):
        node_id, next_id = bad_transition

        transitions = network.transitions.copy()
        x = transitions.get(node_id, [])
        assert next_id in x
        i = x.index(next_id)
        assert i >= 0
        x = x[:i] + x[i + 1:]
        transitions[node_id] = x

        return Network(network.nodes, transitions)


class _OptimizeNoInvalidOutputs:
    # Fixing overlapping data can lead to a situation where a
    # ModuleNode points to DataNode that it can't generate.  E.g.
    # convert_signal_to_tdf -> format=["tdf", "pcl"]
    #   will be changed to:
    # convert_signal_to_tdf -> format=["tdf"]
    # convert_signal_to_tdf -> format=["pcl"]
    #
    # Since one of these is now incorrect, remove it.
    def __init__(self):
        pass

    def optimize(self, network, custom_attributes):
        import copy

        bad_transitions = {}  # (node_id, next_id) -> 1
        for (node_id, next_ids) in network.iterate(node_class=ModuleNode):
            module = network.nodes[node_id]
            for next_id in next_ids:
                node = network.nodes[next_id]
                assert isinstance(node, DataNode)
                if not _is_valid_output(module, node):
                    bad_transitions[(node_id, next_id)] = 1

        network = copy.deepcopy(network)
        for node_id, next_id in bad_transitions:
            x = network.transitions.get(node_id, [])
            assert next_id in x
            i = x.index(next_id)
            assert i >= 0
            x.pop(i)
            network.transitions[node_id] = x
        return network


class _OptimizeNoDuplicateModules:
    def __init__(self):
        pass

    def optimize(self, network, custom_attributes):
        ancestors = _make_ancestor_dict(network)
        nodeid2parents = _make_parents_dict(network)
        while True:
            duplicates = self.find_duplicate_modules(network)
            # Merge each of the duplicates.
            changed = False
            while duplicates:
                n1, n2 = duplicates.pop()

                # Don't merge these if it will create a cycle.  This
                # will happen if one node is an ancestor of the other.
                if n2 in ancestors.get(n1, []) or n1 in ancestors.get(n2, []):
                    continue
                network = network.merge_nodes(
                    [n1, n2], nodeid2parents=nodeid2parents)
                duplicates = _fix_node_id_pairs_after_merge(duplicates, n1, n2)
                nodeid2parents = _fix_node_id_dict_after_merge(
                    nodeid2parents, n1, n2)
                ancestors = _fix_ancestor_dict_after_merge(ancestors, n1, n2)
                changed = True
            if not changed:
                # No duplicates merged.  Either no more duplicates, or
                # would create cycles.
                break
        return network

    def find_duplicate_modules(self, network):
        # Return a list of (node_id1, node_id2) for modules that are
        # duplicated.  If no duplicates found, return an empty list.

        node2parents = _make_parents_dict(network)

        # Duplicated modules:

        # CASE 1: If two modules of the same name have the same
        #         parents, then those modules are duplicated.
        # CASE 2: If two modules of the same name point to the same
        #         children, then those modules are duplicated.
        pairs = {}
        for node_id1 in range(len(network.nodes)-1):
            node1 = network.nodes[node_id1]
            if not isinstance(node1, ModuleNode):
                continue
            children1 = sorted(network.transitions.get(node_id1, []))
            parents1 = sorted(node2parents.get(node_id1, []))
            for node_id2 in range(node_id1+1, len(network.nodes)):
                node2 = network.nodes[node_id2]
                if not isinstance(node2, ModuleNode):
                    continue
                if node1.name != node2.name:
                    continue
                children2 = sorted(network.transitions.get(node_id2, []))
                parents2 = node2parents.get(node_id2, [])
                if children1 and children1 == sorted(children2):
                    pairs[(node_id1, node_id2)] = 1
                if parents1 and parents1 == sorted(parents2):
                    pairs[(node_id1, node_id2)] = 1

        ## pairs = {}
        ## for node_id, next_ids in network.iterate(node_class=DataNode):
        ##     next_ids = sorted(next_ids)
        ##     prev_ids = sorted(node2parents.get(node_id, []))

        ##     for (i, j) in _iter_upper_diag(len(next_ids)):
        ##         node_id1, node_id2 = next_ids[i], next_ids[j]
        ##         pairs[(node_id1, node_id2)] = 1

        ##     # Not sure if this is good, so leave out for now.
        ##     # DEFINITION: If two of the same module nodes points to
        ##     # the same data node, then those modules are duplicated.
        ##     #for (i, j) in _iter_upper_diag(len(prev_ids)):
        ##     #    node_id1, node_id2 = prev_ids[i], prev_ids[j]
        ##     #    pairs[(node_id1, node_id2)] = 1

        dups = []
        for (id1, id2) in pairs:
            node_1 = network.nodes[id1]
            node_2 = network.nodes[id2]
            if node_1.name != node_2.name:
                continue
            dups.append((id1, id2))
        return dups


class _OptimizeNoDuplicateData:
    def __init__(self):
        pass

    def optimize(self, network, custom_attributes):
        # This could be made much more efficient with a better way of
        # finding duplicates.
        ancestors = _make_ancestor_dict(network)
        nodeid2parents = _make_parents_dict(network)
        while True:
            duplicates = self.find_duplicate_data(network)
            changed = False
            # Merge each of the duplicates.
            while duplicates:
                n1, n2 = duplicates.pop()

                # Don't merge these if it will create a cycle.  This
                # will happen if one node is an ancestor of the other.
                if n2 in ancestors.get(n1, []) or n1 in ancestors.get(n2, []):
                    continue

                network = network.merge_nodes(
                    [n1, n2], nodeid2parents=nodeid2parents)
                duplicates = _fix_node_id_pairs_after_merge(duplicates, n1, n2)
                nodeid2parents = _fix_node_id_dict_after_merge(
                    nodeid2parents, n1, n2)
                ancestors = _fix_ancestor_dict_after_merge(ancestors, n1, n2)

                changed = True
            if not changed:
                break
        return network

    def find_duplicate_data(self, network):
        # Return list of (node_id1, node_id2) for DataNode objects
        # that are duplicated.  If no duplicates found, return an
        # empty list.

        # Make a list of all pairs of DataNode objects.
        data_node_ids = [
            node_id for (node_id, node) in enumerate(network.nodes)
            if isinstance(node, DataNode)]

        duplicates = []
        for (i, j) in _iter_upper_diag(len(data_node_ids)):
            node_id1, node_id2 = data_node_ids[i], data_node_ids[j]
            node_1, node_2, = network.nodes[node_id1], network.nodes[node_id2]
            if node_1.datatype.name != node_2.datatype.name:
                continue
            if node_1 == node_2:
                duplicates.append((node_id1, node_id2))
        return duplicates


class _OptimizeMergeData1:
    # Sometimes the inference can lead to two nodes that share the
    # same parents and the same children, and almost the same
    # attributes.  For example:
    # Node1                   Node2
    # preprocess="unknown"    preprocess=<everything else>
    #
    # Exception:
    # If separate nodes are needed for the module, then don't merge
    # them.  It makes the pipeline machinery more complicated.
    # VCFFolder.caller=mutect   ->  call_variants_all
    # VCFFolder.caller=varscan  ->
    #
    # Don't merge PileupSum because they go to different modules.  If
    # merge, then cannot specify different PileupSum as inputs.
    # Don't merge sum_cons_mpileup because they have different
    # children and parents.
    # Don't merge Bam because they go to different modules.
    # Bam.aligner=bwa -> sum_cons_mpileup -> PileupSum.aligner=bwa -> dna_cov
    # Bam.aligner=star -> sum_cons_mpileup -> PileupSum.aligner=star -> rna_cov
    #
    # If this happens, merge them to simplify the network.
    def __init__(self):
        pass

    def optimize(self, network, custom_attributes):
        import copy
        network = copy.deepcopy(network)

        while True:
            nodeid2parents = _make_parents_dict(network)
            similar = self._find_similar_nodes(
                network, custom_attributes, nodeid2parents)
            # Make sure these similar nodes are not separate inputs to
            # the same module.

            # Optimization: Cache calls to _bc_to_input_ids in
            # _is_separate_inputs.
            bc_cache = {}
            i = 0
            while i < len(similar):
                n1, n2 = similar[i]
                if self._is_separate_inputs(
                    network, custom_attributes, nodeid2parents, n1, n2,
                    bc_cache):
                    del similar[i]
                else:
                    i += 1
            if not similar:
                break

            # Merge the similar nodes.
            while similar:
                n1, n2 = similar.pop()
                network = self._merge_nodes(network, n1, n2)
                similar = _fix_node_id_pairs_after_merge(similar, n1, n2)
        return network

    def _is_separate_inputs(
        self, network, custom_attributes, nodeid2parents, node_id1, node_id2,
        bc_cache):
        # Find the common modules for both these nodes.
        children1 = _get_children_of(network, node_id1)
        children2 = _get_children_of(network, node_id2)
        module_ids = [x for x in children1 if x in children2]
        assert module_ids
        for module_id in module_ids:
            # If the module only has one input, then these can't be
            # separate.
            module = network.nodes[module_id]
            if len(module.in_datatypes) == 1:
                continue
            if module_id not in bc_cache:
                combos = _bc_to_input_ids(
                    network, module_id, custom_attributes,
                    nodeid2parents=nodeid2parents)
                bc_cache[module_id] = combos
            combos = bc_cache[module_id]
            for combo in combos:
                if node_id1 in combo and node_id2 in combo:
                    return True
        return False


    def _find_similar_nodes(self, network, custom_attributes, nodeid2parents):
        # Return a list of (node_id1, node_id2).  Can be empty.
        data_node_ids = [
            node_id for (node_id, node) in enumerate(network.nodes)
            if isinstance(node, DataNode)]

        # Group the nodes by datatype.
        dt2nodeids = {} # datatype_name -> list of node_ids
        for node_id in data_node_ids:
            dt = network.nodes[node_id].datatype.name
            if dt not in dt2nodeids:
                dt2nodeids[dt] = []
            dt2nodeids[dt].append(node_id)

        # Optimization: The same calls are made to _fc_to_output_ids,
        # which takes up a lot of the compute time.  Cache these
        # calls.
        fc_cache = {}
        similar = []
        for node_ids in dt2nodeids.itervalues():
            for i, node_id1 in enumerate(node_ids):
                for node_id2 in node_ids[i+1:]:
                    if self._are_nodes_similar(
                        network, node_id1, node_id2, custom_attributes,
                        nodeid2parents, fc_cache):
                        similar.append((node_id1, node_id2))
        return similar

    def _are_nodes_similar(
        self, network, node_id1, node_id2, custom_attributes,
        nodeid2parents, fc_cache):
        node_1 = network.nodes[node_id1]
        node_2 = network.nodes[node_id2]

        # The data type must be the same.
        if node_1.datatype.name != node_2.datatype.name:
            return False

        # They must share the same children.
        c1 = network.transitions.get(node_id1, [])
        c2 = network.transitions.get(node_id2, [])
        if len(c1) != len(c2):
            return False
        if sorted(c1) != sorted(c2):
            return False

        # They might not share the same parents.
        # align_bowtie1 -> SamFolder.aligner (bowtie1)
        # align_bowtie2 -> SamFolder.aligner (bowtie2)
        # Merge to:
        # align_bowtie1 -> SamFolder.aligner (bowtie1, bowtie2)
        # align_bowtie2 ->
        ## They must share the same parents.
        ##p1 = nodeid2parents.get(node_id1, [])
        ##p2 = nodeid2parents.get(node_id2, [])
        ##if len(p1) != len(p2):
        ##    return False
        ##if sorted(p1) != sorted(p2):
        ##    return False

        # They must share all but 1 attribute.
        x, x, diff_attrs = _score_same_data(node_1, node_2)
        if len(diff_attrs) != 1:
            return False

        # After merging, these data nodes must be able to generate all
        # the (grand)children that the unmerged data could generate.
        module_ids = c1
        paths = []  # list of (in_data_ids, module_id, out_data_id)
        for module_id in module_ids:
            if module_id not in fc_cache:
                x = _fc_to_output_ids(
                    network, module_id, custom_attributes,
                    nodeid2parents=nodeid2parents)
                fc_cache[module_id] = x
            paths.extend(fc_cache[module_id])
        # Make sure the input data includes node_id1 or node_id2.
        paths = [x for x in paths if node_id1 in x[0] or node_id2 in x[0]]
        # Make sure the input data does include both node_id1 and node_id2.
        paths = [x for x in paths
                 if not (node_id1 in x[0] and node_id2 in x[0])]


        # The combined data node must be able to generate all these
        # out data nodes.
        merged_data = _merge_data_nodes(node_1, node_2)
        for x in paths:
            in_data_ids, module_id, out_data_id = x
            in_datas = [network.nodes[x] for x in in_data_ids]
            for i in range(len(in_data_ids)):
                if in_data_ids[i] in [node_id1, node_id2]:
                    in_datas[i] = merged_data
                    break
            if not _is_valid_inputs(
                network, module_id, in_datas, custom_attributes,
                out_data_ids=[out_data_id]):
                return False
        return True

    def _merge_nodes(self, network, node_id1, node_id2):
        # Delete the one with the higher node_id (node_id2).
        if node_id1 > node_id2:
            node_id1, node_id2 = node_id2, node_id1

        # Merge the attributes of the nodes.
        n1 = network.nodes[node_id1]
        n2 = network.nodes[node_id2]
        network.nodes[node_id1] = _merge_data_nodes(n1, n2)
        # Everything that pointed to node_id2 now goes to node_id1.
        for node_id, next_ids in network.transitions.iteritems():
            if node_id2 in next_ids and node_id1 not in next_ids:
                next_ids.append(node_id1)
        # They share the same children already.  No need to add.
        return network.delete_node(node_id2)


class _OptimizeMergeData2:
    # is_compressed -> Fastq.trimmed=no -> uncompress
    # is_compressed -> Fastq.trimmed=yes -> uncompress
    # Sometimes will be root nodes (no parents).
    #
    # Actually, don't use this.  It can make the inferencing harder.
    # e.g.
    # Fastq.trimmed (no, yes) -> is_compressed -> Fastq.trimmed (no)
    # Hard to reason whether Fastq.trimmed (no, yes) is a valid
    # antecedent of Fastq.trimmed (no).
    def __init__(self):
        pass

    def optimize(self, network, custom_attributes):
        import copy
        network = copy.deepcopy(network)
        while True:
            similar = self._find_similar_nodes(network)
            if not similar:
                break
            # Merge the similar nodes.
            while similar:
                n1, n2 = similar.pop()
                network = self._merge_nodes(network, n1, n2)
                similar = _fix_node_id_pairs_after_merge(similar, n1, n2)
        return network

    def _find_similar_nodes(self, network):
        # Return a list of (node_id1, node_id2).  Can be empty.
        nodeid2parents = _make_parents_dict(network)

        data_node_ids = [
            node_id for (node_id, node) in enumerate(network.nodes)
            if isinstance(node, DataNode)]

        similar = []
        for i, node_id1 in enumerate(data_node_ids):
            for node_id2 in data_node_ids[i+1:]:
                if self._are_nodes_similar(
                    network, node_id1, node_id2, nodeid2parents):
                    similar.append((node_id1, node_id2))
        return similar

    def _are_nodes_similar(self, network, node_id1, node_id2,
                           nodeid2parents):
        node_1 = network.nodes[node_id1]
        node_2 = network.nodes[node_id2]

        # The data type must be the same.
        if node_1.datatype.name != node_2.datatype.name:
            return False

        # They must share the same children.
        c1 = network.transitions.get(node_id1, [])
        c2 = network.transitions.get(node_id2, [])
        if len(c1) != len(c2):
            return False
        if sorted(c1) != sorted(c2):
            return False

        # They must share the same parents.
        p1 = nodeid2parents.get(node_id1, [])
        p2 = nodeid2parents.get(node_id2, [])
        if len(p1) != len(p2):
            return False
        if sorted(p1) != sorted(p2):
            return False

        # They must share all but 1 attribute.
        x, x, diff_attrs = _score_same_data(node_1, node_2)
        if len(diff_attrs) != 1:
            return False
        return True

    def _merge_nodes(self, network, node_id1, node_id2):
        # Delete the one with the higher node_id (node_id2).
        if node_id1 > node_id2:
            node_id1, node_id2 = node_id2, node_id1

        # Merge the attributes of the nodes.
        n1 = network.nodes[node_id1]
        n2 = network.nodes[node_id2]
        network.nodes[node_id1] = _merge_data_nodes(n1, n2)
        # They share the same parents and children, so nothing needs
        # to be rewired.
        return network.delete_node(node_id2)


def _find_paths_h(network, node_id, custom_attributes, nodeid2parents):
    #import itertools
    assert node_id < len(network.nodes)
    node = network.nodes[node_id]
    prev_ids = nodeid2parents.get(node_id)
    if not prev_ids:
        assert isinstance(node, DataNode)
        yield (node_id,)
        return

    if isinstance(node, DataNode):
        combos = []
        for prev_id in prev_ids:
            combos.append((prev_id,))
    elif isinstance(node, ModuleNode):
        combos = _bc_to_input_ids(
            network, node_id, custom_attributes, nodeid2parents=nodeid2parents)
    for combo in combos:
        # Make a list of the possible paths for each branch.
        branch2paths = []
        for prev_id in combo:  # prev_id is node_id for one branch
            paths = []
            for x in _find_paths_h(
                network, prev_id, custom_attributes, nodeid2parents):
                x = tuple(x)
                paths.append(x)
            assert paths
            branch2paths.append(paths)
        # Merge the paths for each branch.
        for x in _product_and_chain(branch2paths, None):
            x = x + (node_id,)
            yield x


def find_paths(network, custom_attributes, max_paths=None):
    # Iterate over all possible paths from the start nodes to the end
    # nodes.  Each path is a list of the node_ids.
    assert network.nodes, "empty network"
    nodeid2parents = _make_parents_dict(network)
    for i, x in enumerate(
        _find_paths_h(network, 0, custom_attributes, nodeid2parents)):
        yield x
        if max_paths is not None and i >= max_paths:
            break


def _find_paths_by_datatypes_h(
    network, node_id, custom_attributes, datatype_names,
    nodeid2parents, depth):
    # Yield tuples of:
    # path         list of node_ids in this path.
    # used_ids     list of node_ids for nodes from datatype_names
    # missing_ids  list of node_ids not in datatype_names
    import itertools

    assert node_id < len(network.nodes), "%s %d" % (
        repr(node_id), len(network.nodes))
    node = network.nodes[node_id]
    prev_ids = _get_parents_of(network, node_id)

    if isinstance(node, DataNode):
        # If this node is one of these datatypes, then this can be an input.
        if node.datatype.name in datatype_names:
            yield [node_id], [node_id], []
        elif not prev_ids:
            # If this is a start node, then this is a missing input.
            yield [node_id], [], [node_id]
        combos = []
        for prev_id in prev_ids:
            combos.append((prev_id,))
    elif isinstance(node, ModuleNode):
        # Find some combination of inputs that works.
        combos = _bc_to_input_ids(network, node_id, custom_attributes)
    for combo in combos:
        # Each branch is a generator to this recursive function.
        branch2info = [
            _find_paths_by_datatypes_h(
                network, x, custom_attributes, datatype_names, nodeid2parents,
                depth+[node_id])
            for x in combo]
        #branch2info = []
        #for x in combo:
        #    x = _find_paths_by_datatypes_h(
        #        network, user_attributes, x, datatype_names, nodeid2parents,
        #        depth+[node_id])
        #    x = list(x)
        #    branch2info.append(x)

        # Try different combinations of paths for each branch.
        for x in itertools.product(*branch2info):
            # Merge the information from each branch.
            path = []
            possible_ids = []  # either used or missing
            for x in x:
                p, uids, mids = x
                path.extend(p)
                possible_ids.extend(uids)
                possible_ids.extend(mids)
            # path may have duplicates if different branches converge
            # upstream.
            path = {}.fromkeys(path).keys()
            possible_ids = {}.fromkeys(possible_ids).keys()
            # For each of the possible_ids, sort out whether it is
            # used or missing.  Do not use the same datatype more than
            # once.
            names = datatype_names[:]
            used_ids = []
            missing_ids = []
            for id_ in possible_ids:
                n = network.nodes[id_].datatype.name
                if n in names:
                    used_ids.append(id_)
                    names.pop(names.index(n))
                else:
                    missing_ids.append(id_)

            # If this is a DataNode, the transition from prev_id (a
            # ModuleNode) to node_id may be invalid.  e.g.
            # FastqFold (trimmed=no)  -> merge_reads -> FastqFold (trimmed=no)
            # FastqFold (trimmed=yes) ->             -> FastqFold (trimmed=yes)
            # If this is an output FastqFolder where trimmed=no, then this
            # is only valid if the input FastqFolder (trimmed=no) is part
            # of the path.  Only follow the valid paths.
            # Actually, this probably isn't necessary anymore given
            # the new reasoning engine.
            if isinstance(network.nodes[node_id], DataNode):
                if not _is_valid_outdata_id_path(
                    network, path, node_id, custom_attributes, nodeid2parents):
                    continue
            path = path + [node_id]
            yield path, used_ids, missing_ids


def find_paths_by_datatypes(network, custom_attributes, datatype_names):
    # Whether this set of datatypes (by name) provides a complete set
    # of inputs.  Yield tuples of:
    # path         list of node_ids in this path.
    # used_ids     list of node_ids for nodes found in datatype_names
    # missing_ids  list of node_ids not in datatype_names

    nodeid2parents = _make_parents_dict(network)

    # Recursively check if this set of datatype_names can be inputs to
    # this network.
    x = _find_paths_by_datatypes_h(
        network, 0, custom_attributes, datatype_names, nodeid2parents, [])
    return x


class Pathway:
    def __init__(self, node_ids, transitions, start_ids, missing_ids):
        # transitions is dictionary of node_id : frozenset of next_ids
        assert type(node_ids) is type(frozenset())
        assert type(missing_ids) is type(frozenset())
        # Debug: Check for duplicates in node_ids.
        #assert sorted(_uniq(node_ids)) == sorted(node_ids)
        #self.node_ids = node_ids[:]
        self.node_ids = node_ids
        self.transitions = transitions.copy()
        self.start_ids = start_ids[:]
        #self.missing_ids = missing_ids[:]
        self.missing_ids = missing_ids
        self._hash = None
    def __hash__(self):
        if self._hash is None:
            # Optimization: pre-hash this pathway for fast comparisons.
            # This function takes a lot of time.  Still need to
            # optimize this further, somehow.
            #t_hash = []
            #for x in sorted(self.transitions):
            #    t_hash.append((x, frozenset(self.transitions[x])))
            #t_hash = tuple(t_hash)
            # Much faster.
            #t_hash = {}
            #for k, v in self.transitions.iteritems():
            #    t_hash[k] = frozenset(v)
            #t_hash = frozenset(t_hash.items())
            #x = (
            #    frozenset(self.node_ids), t_hash, tuple(self.start_ids),
            #    frozenset(self.missing_ids))
            t_hash = frozenset(self.transitions.items())
            x = (
                self.node_ids, t_hash, tuple(self.start_ids), self.missing_ids)
            self._hash = hash(x)
        return self._hash
    def __eq__(self, other):
        if not self._hash:
            hash(self)
        if not other._hash:
            hash(other)
        return self._hash == other._hash
        #return hash(self) == hash(other)
    def __cmp__(self, other):
        if not isinstance(other, Pathway):
            return cmp(id(self), id(other))
        x1 = [
            self.node_ids, self.transitions, self.start_ids, self.missing_ids]
        x2 = [
            other.node_ids, other.transitions, other.start_ids,
            other.missing_ids]
        return cmp(x1, x2)
        # Compare the hashes.  Order of pathways is meaningless.
        #return cmp(self._hash, other._hash)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [repr(self.node_ids),
             repr(self.transitions),
             repr(self.start_ids),
             repr(self.missing_ids),
             ]
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x


def _find_paths_by_start_ids_hh(
    network, node_id, custom_attributes, node2startids, nodeid2parents,
    depth, cache):
    import itertools
    from genomicode import jmath

    PRUNE_PATHS = True
    #PRUNE_PATHS = False

    assert node_id < len(network.nodes), "%s %d" % (
        repr(node_id), len(network.nodes))
    node = network.nodes[node_id]
    # For debugging.
    #if isinstance(node, DataNode):
    #    node_name = node.datatype.name
    #else:
    #    node_name = node.name

    debug_print(
        DEBUG_FIND_PATHS,
        "[%d] %s finding paths." % (node_id, get_node_name(node)))

    all_paths = {}
    #has_missing = True  # Optimization
    if isinstance(node, DataNode):
        # If this node matches one of the node2startids, then this can
        # be an input.
        for i, start_ids in enumerate(node2startids):
            if node_id not in start_ids:
                continue
            sids = [None] * len(node2startids)
            sids[i] = node_id
            p = Pathway(frozenset([node_id]), {}, sids, frozenset([]))
            all_paths[p] = 1
            #has_missing = False
        # If this node doesn't match any start nodes, then this branch
        # may be missing.
        if not all_paths:
            sids = [None] * len(node2startids)
            p = Pathway(frozenset([node_id]), {}, sids, frozenset([node_id]))
            all_paths[p] = 1
        # Search each of the parents for inputs.
        combos = [(x,) for x in nodeid2parents.get(node_id, [])]
    elif isinstance(node, ModuleNode):
        # Find some combination of inputs that works.
        combos = _bc_to_input_ids(
            network, node_id, custom_attributes, nodeid2parents=nodeid2parents)

    checked_pathway = {}
    for cnum, combo_ids in enumerate(combos):
        # combo is a list of node_ids.
        # Some combos will be seen twice, but the vast majority of the
        # time, each combo is evaluated only once.
        x = [get_node_name(network.nodes[x]) for x in combo_ids]
        debug_print(
            DEBUG_FIND_PATHS, "[%d] Backchaining to %s %s." % (
                node_id, combo_ids, x))

        # For each input, get a list of the Pathway objects.
        innum2branch = [
            _find_paths_by_start_ids_h(
                network, x, custom_attributes, node2startids, nodeid2parents,
                depth+[node_id], cache)
            for x in combo_ids]

        # Make sure there aren't too many combinations to search.
        MAX_COMBINATIONS = 5E4
        branchlen = [len(x) for x in innum2branch]
        total = jmath.prod(branchlen)
        assert total < MAX_COMBINATIONS, "Too many paths (%d)" % total

        if not PRUNE_PATHS:
            # No branches indicate an error somewhere.
            assert total > 0

        assert len(combo_ids) == len(branchlen)
        if DEBUG_FIND_PATHS:
            x = []
            for cid, bl in zip(combo_ids, branchlen):
                x.append("%s(%d)" % (cid, bl))
            x = " ".join(x)
            debug_print(
                DEBUG_FIND_PATHS, "[%d] Found %d possible combinations (%s)." %
                (node_id, total, x))

        # Optimization: If there is already a possibility with
        # nomissing, then don't even consider the paths that will
        # introduce a missing node.
        # This happens occasionally, but not very often.  Does not
        # speed things up.
        #if not has_missing:
        #    has_empty_branch = False
        #    for i in range(len(branch2info)):
        #        x = [x for x in branch2info[i] if not x.missing_ids]
        #        branch2info[i] = x
        #        if not x:
        #            has_empty_branch = True
        #            break
        #    if has_empty_branch:
        #        continue

        # To reduce the effect of combinatorial explosion, merge two
        # nodes at a time.
        # This only speeds up by 5%.
        #merge_cache = {}   # may be None if start IDs conflict
        #merge_i1 = merge_i2 = None
        #if len(branch2info) >= 3:
        #    # Pre-merge the two shortest branchlens.
        #    x = [(branchlen[i], i) for i in range(len(branchlen))]
        #    x = sorted(x)
        #    x = [x[-1] for x in x]
        #    merge_i1 = x[0]   # shortest
        #    merge_i2 = x[1]   # second shortest
        #
        #    key = [None] * len(branch2info)
        #    for j in range(branchlen[merge_i1]):
        #        key[merge_i1] = j
        #        for k in range(branchlen[merge_i2]):
        #            key[merge_i2] = k
        #            b1 = branch2info[merge_i1][j]
        #            b2 = branch2info[merge_i2][k]
        #            x = _merge_paths([b1, b2])
        #            node_ids, transitions, missing_ids = x
        #
        #            p = None
        #            conflict = _does_based_on_data_conflict_with_out_data(
        #                network, node_ids, transitions)
        #            # There's a conflict 8% of the time.
        #            sids = None
        #            if not conflict:
        #                sids = _merge_start_ids([b1, b2])
        #            if sids:
        #                p = Pathway(node_ids, transitions, sids, missing_ids)
        #            merge_cache[tuple(key)] = p

        # Try different combinations of paths for each branch.
        for bnum, branches in enumerate(itertools.product(*innum2branch)):
        #x = [range(x) for x in branchlen]
        #for bnum, branches_i in enumerate(itertools.product(*x)):
        #    branches = [
        #        branch2info[i][branches_i[i]] for i in range(len(branches_i))]

            # This _merge_paths is called 99.999% of the time.
            #if merge_cache:
            #    key = [None] * len(branch2info)
            #    key[merge_i1] = branches_i[merge_i1]
            #    key[merge_i2] = branches_i[merge_i2]
            #    x = merge_cache[tuple(key)]
            #    if x is not None:
            #        x1 = [x]
            #        x2 = [branches[i] for i in range(len(branches))
            #              if i not in [merge_i1, merge_i2]]
            #        x = _merge_paths(x1+x2)
            #else:
            #    x = _merge_paths(branches)
            #if x is None:
            #    continue
            #node_ids, transitions, missing_ids = x
            node_ids, transitions, missing_ids = _merge_paths(branches)

            # Optimization: If there's already a nomissing path, then
            # skip all paths with missing_ids.
            # Actually, this never happens.
            #if not has_missing and missing_ids:
            #    continue

            # Add node_id to this set of node_ids.
            assert node_id not in node_ids
            node_ids = node_ids.union([node_id])
            # Add transitions to this node.
            for x in combo_ids:
                if x not in transitions:
                    transitions[x] = frozenset([node_id])
                else:
                    transitions[x] = transitions[x].union([node_id])


            # See if these paths conflict.  Paths may conflict if:
            # 1.  Different data nodes are used for the same start_id.
            # 2.  A module_node has an consequence that is
            #     BASED_ON_DATA, and the subsequent DataNode has
            #     different values for that attribute.
            #     Ignore this check if there are missing_ids.  Is
            #     already an invalid merged pipeline.
            # 3.  If data node, can't transition from previous module
            #     node.
            # 4.  If module node, can't transition from previous data
            #     nodes.
            # 5.  The same data node is seen in different branches,
            #     but it has different parents in each branch.
            
            # If the paths conflict, then skip it.

            # Case 1.  See if there are conflicting start_ids.
            # (As a side effect, also merge the start_ids).
            start_ids = _merge_start_ids(branches)
            if start_ids is None:
                # This almost never happens.
                continue
            
            # Merge each of the branches into one path.
            # Optimization: If this pathway has already been checked,
            # then don't check it again.
            path = Pathway(node_ids, transitions, start_ids, missing_ids)
            if path in checked_pathway:
                # This happens very often.  This can happen if there
                # are multiple equivalent routes through a network
                # (e.g. steps done in different order).  Branches take
                # different combinations of routes, but when you merge
                # them, they end up the same.
                continue
            checked_pathway[path] = 1

            # Case 2.  Look for modules with a BASED_ON_DATA
            # consequence.  Then, make sure the attribute values of
            # the data are the same.

            # If there are missing_ids, there may be lots of conflicts
            # because everything's merged.  So only examine this if
            # there are no missing_ids.
            if not missing_ids:
                conflict = _does_based_on_data_conflict_with_out_data(
                    network, node_ids, transitions)
                if conflict and not missing_ids:
                    # This happens 60% of the time.
                    continue

            # Case 3.  If this is a DataNode, the transition from
            # prev_id (a ModuleNode) to node_id may be invalid.  Only
            # follow the valid paths.
            # Example:
            # FastqFold (trimmed=no)  -> merge_reads -> FastqFold (trimmed=no)
            # FastqFold (trimmed=yes) ->             -> FastqFold (trimmed=yes)
            # If this is an output FastqFolder where trimmed=no, then this
            # is only valid if the input FastqFolder (trimmed=no) is part
            # of the path.
            if isinstance(network.nodes[node_id], DataNode):
                if not _is_valid_outdata_id_path(
                    network, node_ids, node_id, custom_attributes,
                    nodeid2parents):
                    # This happens 3% of the time.
                    continue

            # Case 4: If this is a ModuleNode, then make sure
            # transition from previous nodes is valid.
            # 
            # Example:
            # Bam.mouse=no -> sort -> Bam.mouse=yes,no -> count_with_htseq
            #               ReadStrandedness.mouse=yes ->
            # count_with_htseq has constraint that mouse is the same,
            # but value for mouse conflicts.
            #
            # Another example:
            # BamFolder.trimmed=yes        -> count_with_htseq
            # ReadStrandedness.trimmed=no  ->
            # Because further up the tree, different parents for align
            # get merged:
            # Fastq.trimmed=no -> trim -> Fastq.trimmed=yes -> align
            # Fastq.trimmed=no ->                              align
            if not missing_ids:
                # Don't bother checking if this is a dead pathway
                # (missing start_ids).  The merging below can ccreate
                # a non-atomic DataNode that will lead to an exception
                # from _is_valid_outmodule_id_path, e.g.:
                # BamFold.dup=no -> mrk_dup -> BF.dup=yes -> index -> BF.dup=?
                #                ->                                ->
                if isinstance(network.nodes[node_id], ModuleNode):
                    # This is run multiple times for similar pathways.
                    # May be able to memoize.
                    if not _is_valid_outmodule_id_path(
                        network, branches, combo_ids, node_id,
                        custom_attributes, nodeid2parents):
                        continue

            # Case 5. If the merged branches would generate data nodes
            # with different parents, then this conflicts.
            if not missing_ids:
                multiple_parents = False
                for nid in path.node_ids:
                    if not isinstance(network.nodes[nid], DataNode):
                        continue
                    x = _get_pathway_parents_of(
                        network, path, nid, nodeid2parents=nodeid2parents)
                    if len(x) >= 2:
                        multiple_parents = True
                        break
                if multiple_parents:
                    continue

            #if not path.missing_ids and has_missing:
            #    has_missing = False
            assert path not in all_paths

            debug_print(
                DEBUG_FIND_PATHS, "[%d] Adding path %s." % (node_id, path))

            all_paths[path] = 1

    # If any of the paths have no missing nodes, then remove all paths
    # with missing nodes.
    no_missing = [x for x in all_paths if not x.missing_ids]
    if no_missing:
        paths = no_missing
    else:
        # If there are no working paths, then merge the paths by
        # start_ids.
        # Can not do this with working paths.  Otherwise, this may end
        # up merging alternate routes through the network.
        # Fastq.compress (unknown) -> is_compressed -> Fastq.compress (yes)
        #                                           -> Fastq.compress (no)
        # This can create paths that don't work.
        sids2paths = {}  # start_ids -> list of paths
        for path in all_paths:
            sids = tuple(path.start_ids)
            if sids not in sids2paths:
                sids2paths[sids] = []
            sids2paths[sids].append(path)
        paths = []
        for ps in sids2paths.itervalues():
            p = ps[0]
            if len(ps) > 1:
                node_ids, transitions, missing_ids = _merge_paths(ps)
                p = Pathway(
                    node_ids, transitions, ps[0].start_ids, missing_ids)
            paths.append(p)
        debug_print(DEBUG_FIND_PATHS, "[%d] No working paths." % node_id)

    # Only prune working paths.  These paths start from a node, but
    # may not end up at the bottom node yet.
    if no_missing and PRUNE_PATHS:
        # Can't prune by custom attributes here.  Because paths aren't
        # finished, we don't know which are the bottom-mode nodes to
        # apply the attributes to.
        # Also don't prune by superset pipelines.  Can't be sure if
        # the pipelines aren't finished yet.
        #orig_paths = paths
        debug_print(DEBUG_FIND_PATHS,
                    "[%d] Pruning from %d path(s)." % (node_id, len(paths)))
            
        paths = prune_paths(
            paths, network, custom_attributes,
            prune_custom_attributes=False, prune_superset=False,
            ignore_incomplete_pathway=True, prune_most_inputs=False)
        debug_print(DEBUG_FIND_PATHS,
                    "[%d] Final %d path(s) remain." % (node_id, len(paths)))

    return paths


def _find_paths_by_start_ids_h(
    network, node_id, custom_attributes, node2startids, nodeid2parents,
    depth, cache):
    if node_id not in cache:
        x = _find_paths_by_start_ids_hh(
            network, node_id, custom_attributes, node2startids, nodeid2parents,
            depth, cache)
        cache[node_id] = x
    return cache[node_id]


def find_paths_by_start_ids(network, custom_attributes, node2startids):
    # node2startids should be a list of lists indicating the possible
    # start_ids for each input node.
    #
    # This function will search through the network for pipelines that
    # start from this and return a list of Pathway objects.
    #
    # node_ids     list of node IDs in this path.
    # transitions  node_id -> tuple of next node IDs
    # start_ids    list of node IDs, parallel to node2startids
    # missing_ids  list of node IDs
    #
    # node_ids is only provided if there are no missing_ids.

    # Pre-calculate the parents.
    nodeid2parents = _make_parents_dict(network)

    # Pre-calculate which node IDs have BASED_ON_DATA consequences.
    # No. Does not save much time.
    #moduleid2bod = {}  # module_id -> list of attrs that are BASED_ON_DATA
    #for module_id in range(len(network.nodes)):
    #    if not isinstance(network.nodes[module_id], ModuleNode):
    #        continue
    #    x = network.nodes[module_id].consequences
    #    x = [x for x in x if x.behavior == BASED_ON_DATA]
    #    x = [x.name for x in x]
    #    if x:
    #        moduleid2bod[module_id] = x

    x = _find_paths_by_start_ids_h(
        network, 0, custom_attributes, node2startids, nodeid2parents,
        [], {})
    # One of these may be a trivial pathway [0].  Remove this.
    #x = [x for x in x if x[0] != [0]]
    x = [x for x in x if x.node_ids != [0]]
    return x


def _merge_paths(paths):
    # Merge a list of paths and return a tuple of (node_ids,
    # transitions, missing_ids).  start_ids are handled separately.

    # Originally:
    #   2% of time have 1 path.
    #   3% of time have 2 paths.
    # Now that we pre-cache pairs of paths:
    #   99% of time have 2 paths.
    # Writing code specifically for the 2 path case doesn't speed things up.
    # Caching the previously calculated transitions doesn't speed things up.

    # Merge the node_ids and missing_ids.
    x1 = [x.node_ids for x in paths]
    x2 = [x.missing_ids for x in paths]
    node_ids = frozenset().union(*x1)
    missing_ids = frozenset().union(*x2)
    # 76% of time is spent merging the transitions.
    transitions = _merge_transitions(paths)
    return node_ids, transitions, missing_ids


def _merge_transitions(paths):
    # Merge the transitions.
    transitions = paths[0].transitions.copy()
    for p in paths[1:]:
        for k, v in p.transitions.iteritems():
            if k not in transitions:
                transitions[k] = v
            elif transitions[k] != v:
                # 90% of the time, there are only 2 values.
                # This union called 6x more than the ones above.
                transitions[k] = transitions[k].union(v)
    return transitions


def prune_paths(paths, network, custom_attributes,
                prune_custom_attributes=True, prune_superset=True,
                ignore_incomplete_pathway=False, prune_most_inputs=True):
    # This may be called by _find_paths_by_start_ids when the pathway
    # is incomplete.  If this is the case, then may not be able to
    # resolve TYPE_ENUM.  If we can't resolve because the pathway is
    # incomplete, return None rather than raise an exception.
    nodeid2parents = _make_parents_dict(network)

    # Do the O(N) pruning, then fast O(NN) pruning, then slow O(NN)
    # pruning.
    # Just try different orders until I find the fastest.

    ns = len(paths)
    paths = _prune_by_generated_inputs(network, paths)
    name = "_prune_by_generated_inputs"
    debug_print(
        DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (name, ns, len(paths)))

    ns = len(paths)
    paths = _prune_by_complete_inputs(network, paths)
    name = "_prune_by_complete_inputs"
    debug_print(
        DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (name, ns, len(paths)))

    if prune_most_inputs:
        ns = len(paths)
        paths = _prune_by_most_inputs(network, paths)
        name = "_prune_by_most_inputs"
        debug_print(
            DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (
                name, ns, len(paths)))

    if prune_custom_attributes:
        ns = len(paths)
        paths = _prune_by_custom_attributes(
            network, custom_attributes, paths, nodeid2parents)
        name = "_prune_by_custom_attributes"
        debug_print(
            DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (
                name, ns, len(paths)))
        # If len(paths) == 0, suggests that there's a problem with
        # custom attributes.
 
    #ns = len(paths)
    #paths = _prune_alternate_attributes3(
    #    network, custom_attributes, paths, nodeid2parents)
    #name = "_prune_alternate_attributes3"
    #debug_print(
    #    DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (name, ns, len(paths)))

    ns = len(paths)
    paths = _prune_alternate_attributes2(
        network, custom_attributes, paths, nodeid2parents)
    name = "_prune_alternate_attributes2"
    debug_print(
        DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (name, ns, len(paths)))

    ns = len(paths)
    paths = _prune_alternate_attributes1(
        network, custom_attributes, paths, nodeid2parents,
        ignore_incomplete_pathway)
    name = "_prune_alternate_attributes1"
    debug_print(
        DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (name, ns, len(paths)))

    if prune_superset:
        ns = len(paths)
        paths = _prune_superset_pipelines(network, paths, nodeid2parents)
        name = "_prune_superset_pipelines"
        debug_print(
            DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (
                name, ns, len(paths)))

    ns = len(paths)
    paths = _prune_parallel_pipelines(network, paths, nodeid2parents)
    name = "_prune_parallel_pipelines"
    debug_print(
        DEBUG_PRUNE_PATHS, "Pruned [%s] %d -> %d." % (name, ns, len(paths)))

    return paths


def _prune_by_generated_inputs(network, paths):
    # If a pathway contains an --input specified by the user, then
    # there should be no modules that create this input.

    to_delete = {}
    for i, p in enumerate(paths):
        start_ids = [x for x in p.start_ids if x is not None]
        for next_ids in p.transitions.itervalues():
            # If something transitions to a start_id, then delete this
            # pipeline.
            x = next_ids.intersection(start_ids)
            if x:
                to_delete[i] = 1
                break
    if to_delete:
        paths = [p for (i, p) in enumerate(paths) if i not in to_delete]
    return paths


def _prune_by_complete_inputs(network, paths):
    # If one of the pathways is complete, then choose only the
    # pathways that contain complete inputs.
    if not paths:
        return paths
    complete = {}
    for i, p in enumerate(paths):
        if None not in p.start_ids:
            complete[i] = i
    if complete and len(complete) < len(paths):
        paths = [paths[i] for i in sorted(complete)]
    return paths


def _prune_by_most_inputs(network, paths):
    if not paths:
        return paths
    max_inputs = None
    for i, p in enumerate(paths):
        x = [x for x in p.start_ids if x is not None]
        if max_inputs is None or len(x) > max_inputs:
            max_inputs = len(x)
    assert max_inputs is not None
    to_delete = {}
    for i, p in enumerate(paths):
        x = [x for x in p.start_ids if x is not None]
        if len(x) < max_inputs:
            to_delete[i] = 1
    if to_delete:
        paths = [p for (i, p) in enumerate(paths) if i not in to_delete]
    return paths



def _prune_by_custom_attributes(
    network, custom_attributes, paths, nodeid2parents):
    # Keep only the paths that match the attributes desired by the
    # user.  custom_attributes is a list of CustomAttribute objects.
    if not custom_attributes:
        return paths

    # A module's input node is subject to the user's attributes if:
    # 1.  The input datatype is different from the output datatype.
    # AND
    # 2.  The input datatype has at least one custom attribute that
    #     is not subject to a MUST_BE constraint.  The ones that are
    #     subject to a MUST_BE constraint must match the custom
    #     attribute.
    # AND
    # 3.  The data node has no descendents with the same datatype.
    #     Must be bottom-most of that type.

    # Strategy to generate nodes that should be subject to custom
    # attributes.
    # 1.  Start with pathway node IDs.
    # 2.  Only node IDs with a datatype subject to a custom attribute.
    # 3.  Make a list of the module IDs that these data nodes can go
    #     through.
    # 4.  Calculate the inputs node IDs into each of the modules.
    #     Variable is called sub_pathways.
    # 5.  Only want sub_pathways where the input datatype is different
    #     from the output datatype.
    # 6.  Assign each of the custom attributes to each of the nodes.
    # 7.  Make sure the attribute does not conflict with a MUST_BE
    #     constraint.
    # 8.  Make sure the data node is bottom-most data node of that
    #     type.  (for custom_attributes where all_nodes=False)

    # First, cache some variables for convenience.
    # Cache names of data types with custom attributes.
    all_dnames = {}
    for x in custom_attributes:
        all_dnames[x.datatype.name] = 1
    # Cache node IDs that belong in a path.
    x = [x.node_ids for x in paths]
    path_node_ids = x[0].union(*x[1:])
    # Cache some network related dictionaries.
    nodeid2parents = _make_parents_dict(network)
    descendents = _make_descendent_dict(network)


    # Step 1: Start with a list of all data node IDs.
    x = path_node_ids
    x = [x for x in x if isinstance(network.nodes[x], DataNode)]
    data_node_ids = x
    debug_print(
        DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "%d data_node_ids (all)" %
        len(data_node_ids))

    # Step 2: Only node IDs with a datatype subject to a custom
    # attribute.
    x = data_node_ids
    x = [x for x in x if network.nodes[x].datatype.name in all_dnames]
    data_node_ids = set(x)
    debug_print(
        DEBUG_PRUNE_CUSTOM_ATTRIBUTES,
        "%d data_node_ids (same datatype as custom attribute)" %
        len(data_node_ids))

    # Step 3: Make a list of the module IDs that these data nodes can
    # go through.
    module_ids = []
    for node_id in data_node_ids:
        x = network.transitions.get(node_id, [])
        module_ids.extend(x)
    module_ids = [x for x in module_ids if x in path_node_ids]
    module_ids = set(module_ids)
    debug_print(
        DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "%d module_ids" % len(module_ids))

    # Step 4: Calculate the inputs node IDs into each of the modules.
    # list of (combo_ids, module_id, data_node_id, i_combo)
    #
    # combo_ids (including data_node_id) -> module_id
    #
    # combo_ids     list of node_ids that are the inputs to this module
    # data_node_id  node_id that may be subject to a custom attribute
    # i_combo       index of data_node_id into combo_ids
    sub_pathways = []
    for module_id in module_ids:
        module = network.nodes[module_id]

        x = nodeid2parents[module_id]
        x = [x for x in x if x in path_node_ids]
        combos = _bc_to_input_ids(
            network, module_id, custom_attributes, all_input_ids=x,
            nodeid2parents=nodeid2parents)
        for combo_ids in combos:
            for i, node_id in enumerate(combo_ids):
                if node_id not in data_node_ids:
                    continue
                x = combo_ids, module_id, node_id, i
                sub_pathways.append(x)
    debug_print(
        DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "%d sub_pathways" % len(sub_pathways))

    # Step 5: Only want sub_pathways where the input datatype is
    # different from the output datatype.
    good = []
    for (combo_ids, module_id, data_node_id, i_combo) in sub_pathways:
        module = network.nodes[module_id]

        # If the module takes only one kind of datatype and produces
        # the same kind, then this doesn't apply.
        if len(module.in_datatypes) == 1 and \
           module.out_datatype == module.in_datatypes[0]:
            continue

        # Make sure the input datatype is not DefaultAttributesFrom.
        daf = [x.input_index for x in module.default_attributes_from]
        if i_combo in daf:
            continue
        good.append((combo_ids, module_id, data_node_id, i_combo))
    sub_pathways = good
    debug_print(
        DEBUG_PRUNE_CUSTOM_ATTRIBUTES,
        "%d sub_pathways (input and output different)" % len(sub_pathways))

    # Step 6: Assign each of the custom attributes to each of the
    # nodes.  Make sure the data types match.
    good = []
    for (combo_ids, module_id, data_node_id, i_combo) in sub_pathways:
        node = network.nodes[data_node_id]
        for cattrs in custom_attributes:
            if cattrs.datatype.name != node.datatype.name:
                continue
            x = combo_ids, module_id, data_node_id, i_combo, cattrs
            good.append(x)
    sub_pathways = good
    debug_print(
        DEBUG_PRUNE_CUSTOM_ATTRIBUTES,
        "%d sub_pathways (matched to custom attributes)" % len(sub_pathways))

    # Step 7: Make sure the attribute does not conflict with a MUST_BE
    # constraint.
    good = []
    for (combo_ids, module_id, data_node_id, i_combo, cattrs) in sub_pathways:
        module = network.nodes[module_id]

        # For convenience, list the MUST_BE constraints for this module.
        cons_must_be = {}  # input index -> attr name -> value
        for cons in module.constraints:
            if cons.behavior != MUST_BE:
                continue
            index = cons.input_index
            assert type(index) is type(0)
            if index not in cons_must_be:
                cons_must_be[index] = {}
            cons_must_be[index][cons.name] = cons.arg1

        # If the custom attributes conflict with this constraint,
        # then don't use it.
        cons_attrs = cons_must_be.get(i_combo, {})
        conflict = False
        for attr in cattrs.attributes:
            if attr.name in cons_attrs and attr.value != cons_attrs[attr.name]:
                conflict = True
                break
        if conflict:
            continue
        good.append((combo_ids, module_id, data_node_id, i_combo, cattrs))
    sub_pathways = good
    debug_print(
        DEBUG_PRUNE_CUSTOM_ATTRIBUTES,
        "%d sub_pathways (no MUST_BE conflicts)" % len(sub_pathways))

    # Step 8: Make sure the data node is bottom-most data node of that
    # type.

    # Cache the node IDs that don't have any descendents.
    no_desc = set()
    for node_id in data_node_ids:
        x = descendents.get(node_id, [])
        x = [x for x in x if x in path_node_ids]
        x = [x for x in x if isinstance(network.nodes[x], DataNode)]
        x = [network.nodes[x].datatype.name for x in x]
        if network.nodes[node_id].datatype.name not in x:
            no_desc.add(node_id)

    good = []
    for (combo_ids, module_id, data_node_id, i_combo, cattrs) in sub_pathways:
        node = network.nodes[data_node_id]

        if data_node_id not in no_desc and not cattrs.all_nodes:
            continue
        good.append((combo_ids, module_id, data_node_id, i_combo, cattrs))
    sub_pathways = good
    debug_print(
        DEBUG_PRUNE_CUSTOM_ATTRIBUTES,
        "%d sub_pathways (bottom-most)" % len(sub_pathways))


    # At this point, we have a list of all the nodes that are subject
    # to custom_attributes, as well as an assignment to the proper
    # custom_attribute.  List the node_ids that either conflict with
    # the custom attributes, or are ambiguous, e.g.
    #     Fastq.adapters = [no, yes]; user wants no
    # If the values are all the same, then ignore.
    match_node_ids = {}      # node_id -> CustomAttributes
    mismatch_node_ids = {}   # node_id -> CustomAttributes
    ambiguous_node_ids = {}  # node_id -> CustomAttributes
    for (combo_ids, module_id, data_node_id, i_combo, cattrs) in sub_pathways:
        node = network.nodes[data_node_id]

        # Print out the nodes.
        debug_print(
            DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "[%d] %s (%d) -> [%d] %s" % (
                data_node_id, get_node_name(node), i_combo,
                module_id, get_node_name(network.nodes[module_id])))
        for attr in sorted(cattrs.attributes):
            debug_print(
                DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "  %s=%s (node is %s)" % (
                    attr.name, attr.value, node.attributes[attr.name]))

        # If the custom attributes conflict with this data node, then
        # don't use it.
        conflict = False    # any attribute conflicts
        ambiguous = False   # any attribute is ambiguous
        for attr in cattrs.attributes:
            assert attr.name in node.attributes
            uvalue = attr.value                    # user value
            dvalue = node.attributes[attr.name]    # network value
            utype = _get_attribute_type(uvalue)
            dtype = _get_attribute_type(dvalue)

            # CASE   USER   NETWORK   RESULT
            #   1    ATOM    ATOM     CONFLICT if items aren't equal.
            #   2    ATOM    ENUM     AMBIG if ATOM in ENUM.  else CONFLICT.
            #   3    ENUM    ATOM     error
            #   4    ENUM    ENUM     error
            assert utype is TYPE_ATOM
            if dtype == TYPE_ATOM:
                if uvalue != dvalue:
                    conflict = True
            elif dtype == TYPE_ENUM:
                if uvalue in dvalue:
                    ambiguous = True
                else:
                    conflict = True
            else:
                raise AssertionError
        if conflict:
            mismatch_node_ids[data_node_id] = cattrs
            debug_print(DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "conflict")
        elif ambiguous:
            ambiguous_node_ids[data_node_id] = cattrs
            debug_print(DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "ambiguous")
        else:
            match_node_ids[data_node_id] = cattrs
            debug_print(DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "good")

    # If a node matches a custom attribute (or is ambiguous), then
    # don't make it a mismatch to another.
    x = match_node_ids.keys() + ambiguous_node_ids.keys()
    for node_id in x:
        if node_id in mismatch_node_ids:
            del mismatch_node_ids[node_id]

    if not mismatch_node_ids and not ambiguous_node_ids:
        return paths

    # For each pathway, see if it contains a mismatch or ambiguous node.
    bc_cache = {}
    fc_cache = {}
    path_cache = {}

    delete = {}
    reason = {}  # path_i -> reason why deleted (for DEBUGGING)
    for (i, p) in enumerate(paths):
        debug_print(
            DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "Checking path %d." % i)

        # If this path contains a node that does not match custom
        # attributes, then delete this path.
        mismatch_ids = p.node_ids.intersection(mismatch_node_ids)
        if mismatch_ids:
            delete[i] = 1
            reason[i] = "Nodes (%s) do not match custom attribute." % \
                        ",".join(map(str, mismatch_ids))
            debug_print(
                DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "%d.  %s" % (i, reason[i]))
            continue

        # Make a list of the ambiguous node IDs in this path.  If this
        # path does not contain any ambiguous nodes, then keep this
        # path.
        ambig_ids = p.node_ids.intersection(ambiguous_node_ids)
        if not ambig_ids:
            continue

        # See if this pathway can produce a node that has an attribute
        # that matches the custom_attribute.

        # Make a list of all the sub-pathways that generate this
        # ambiguous node.
        check_paths = []  # list of (in_ids, module_id, out_id)
        for out_data_id in ambig_ids:
            x = _bc_to_input_and_module_ids(
                network, out_data_id, custom_attributes, p.node_ids,
                nodeid2parents, bc_cache)
            check_paths.extend(x)
        assert check_paths

        # Make a list of all the node IDs from ambig_ids that matches
        # the custom attribute.
        match = {}  # node_id -> 1
        for x in check_paths:
            in_data_ids, module_id, out_data_id = x
            if out_data_id in match:  # don't check if already done
                continue

            cattrs = ambiguous_node_ids[out_data_id]
            attrs = {}
            for x in cattrs.attributes:
                attrs[x.name] = x.value

            #name = network.nodes[out_data_id].datatype.name
            #attrs = dname2attrs[name]
            #attrs = nodeid2attrs[out_data_id]
            key = module_id, in_data_ids, out_data_id
            if key not in path_cache:
                # See if this set of attributes can be created from
                # this module and inputs.
                x = _is_valid_output_from_input_and_module_ids(
                    network, in_data_ids, module_id, out_data_id,
                    attrs, fc_cache)
                path_cache[key] = x
            if path_cache[key]:
                match[out_data_id] = 1
            # Optimization: stop checking if we've already resolved
            # all the ambiguous nodes.
            if len(match) == len(ambig_ids):
                break
        # If we cannot generate all ambiguous nodes (based on the
        # custom attributes), then delete this pathway.
        if len(match) != len(ambig_ids):
            delete[i] = 1
            x = [x for x in ambig_ids if x not in match]
            reason[i] = \
                      "Custom attributes can not be generated in nodes (%s)." \
                      % ",".join(map(str, x))
            debug_print(
                DEBUG_PRUNE_CUSTOM_ATTRIBUTES, "%d.  %s" % (i, reason[i]))

    paths = [x for (i, x) in enumerate(paths) if i not in delete]
    return paths


def _prune_alternate_attributes1(
    network, custom_attributes, paths, nodeid2parents,
    ignore_incomplete_pathway):
    # Distinct pipelines may generate the same data node (but with
    # differing attributes).  Only one of these pipelines are needed.
    # Arbitrarily choose one, and delete the other pipelines.
    #
    # align_with_bowtie1 -> SamFolder.align=[bowtie1,bowtie2,bwa_back,bwa_mem]
    # align_with_bowtie2 -> SamFolder.align=[bowtie1,bowtie2,bwa_back,bwa_mem]
    # align_with_bwa_aln -> SamFolder.align=[bowtie1,bowtie2,bwa_back,bwa_mem]
    # align_with_bwa_mem -> SamFolder.align=[bowtie1,bowtie2,bwa_back,bwa_mem]

    # Approach:
    # 1.  Look for DataObject with attribute with multiple values.
    #     - Has attribute is TYPE_ENUM.
    # 2.  Split up the paths based on which value is set.
    #     - Different Modules pointing to that DataObject.
    #     - Modules have Consequence SET_TO to different values in the
    #       same attribute.
    # 3.  Choose the value arbitrarily (which comes first
    #     alphabetically).
    # 4.  Remove the pipelines that don't lead to that value.
    #     - Be careful not to remove pipelines that do not contain the
    #       data object.
    # 5.  Print a message showing which value was chosen.

    # TODO: Make this obsolete.
    ignore_incomplete_pathway = True

    if not paths:
        return []
    paths = paths[:]

    # Make a list of all the node_ids in the paths.
    x = [x.node_ids for x in paths]
    path_node_ids = x[0].union(*x[1:])

    # Make list of all nodes that have alternate values.
    alternates = []  # list of node_id, attr_name
    for node_id, node in enumerate(network.nodes):
        if node_id not in path_node_ids:
            continue
        if not isinstance(node, DataNode):
            continue
        for (name, value) in node.attributes.iteritems():
            t = _get_attribute_type(value)
            if t != TYPE_ENUM:
                continue
            x = node_id, name
            alternates.append(x)
    if not alternates:
        return paths

    # Only want the alternates with multiple parent nodes.
    nodeid2parents = _make_parents_dict(network)
    good = []
    for (node_id, attr_name) in alternates:
        x = nodeid2parents.get(node_id, [])
        x = [x for x in x if x in path_node_ids]
        parent_ids = x
        if len(parent_ids) > 1:
            good.append((node_id, attr_name))
    alternates = good

    # Clean up one alternate at a time to avoid conflicts.  Don't want
    # a situation in which one alternate leads to an impossible
    # situation for another alternate.
    for (node_id, name) in alternates:
        # If this alternate is not in a pathway somewhere (either
        # never in, or was deleted), then ignore.
        x = [i for (i, p) in enumerate(paths) if node_id in p.node_ids]
        if not x:
            continue

        # Sort the pipelines based on the values.
        dont_resolve = False
        value2pathids = {}  # value -> list of pathway indexes
        for index, path in enumerate(paths):
            if node_id not in path.node_ids:
                continue
            # This function gets called multiple times for similar
            # pathways.  Can find a way to optimize this.
            try:
                atomic_node = _get_atomic_data_node_from_pathway(
                    network, path, node_id, custom_attributes,
                    ignore_based_on_data=True, ignore_unchanged=True,
                    nodeid2parents=nodeid2parents)
            except AssertionError, x:
                # This can fail for multiple reasons.
                # 1.  The pathway is incomplete, and the leaf nodes
                #     could not be found.
                #     Can ignore this.
                # 2.  A node is ambiguous because it is BASED_ON_DATA.
                #     check_for_missing_values -> File.missing=[no, yes]
                #     Can ignore this.
                # 3.  A value can be ambiguous, but not changed from
                #     the parents.  Can ignore this.
                #
                # Careful: _get_atomic_data_node_from_pathway may
                # generate multiple errors.  Be sure to handle this
                # possibility.
                if str(x).startswith("No parents") and \
                       ignore_incomplete_pathway:
                    # This can happen if the pathway isn't finished
                    # generating, or the user has given an internal
                    # node as an --input.
                    dont_resolve = True
                else:
                    raise
            if dont_resolve:
                break

            value = atomic_node.attributes[name]
            # Can happen for BASED_ON_DATA or unchanged values.
            if _get_attribute_type(value) == TYPE_ENUM:
                value = tuple(value)
            #assert _get_attribute_type(value) == TYPE_ATOM
            if value not in value2pathids:
                value2pathids[value] = []
            value2pathids[value].append(index)

        # Cannot resolve TYPE_ENUM.  Ignore this alternate.
        if dont_resolve:
            continue

        # If there is only one possible value, then ignore this.
        assert value2pathids
        if len(value2pathids) == 1:
            continue

        # Choose a value.
        chosen_value = None
        # Look in custom_attributes to see if one is favored by the user.
        # Not sure if this is correct.  Need to look carefully to make
        # sure right value goes to right node?
        values = _get_custom_values(
            custom_attributes, network.nodes[node_id].datatype.name, name)
        assert len(values) <= 1, "Multiple custom_attributes."
        if values:
            chosen_value = values[0]
        if chosen_value is None:
            if "no" in value2pathids:
                chosen_value = "no"
        if chosen_value is None:
            # Choose the first one alphabetically.
            chosen_value = sorted(value2pathids)[0]

        # Ignore all the pathway IDs that does not have this value.
        # Be sure not to delete pathways that do not include this
        # node.
        to_delete = []
        for (value, path_ids) in value2pathids.iteritems():
            if value == chosen_value:
                continue
            to_delete.extend(path_ids)
        to_delete = {}.fromkeys(to_delete)

        paths = [p for (i, p) in enumerate(paths) if i not in to_delete]

    #_print_attributes_pruned(deleted)
    return paths


## def _prune_alternate_attributes1(
##     network, custom_attributes, paths, nodeid2parents):
##     # If there an object has an attribute that can take different
##     # values, then pipelines may be generated that can lead to each
##     # value.  This causes unnecessary computation, because only one
##     # value needs to be used.  Arbitrarily choose one, and delete the
##     # other pipelines.
##     #
##     # align_with_bowtie1 -> SamFolder.align=[bowtie1,bowtie2,bwa_back,bwa_mem]
##     # align_with_bowtie2 -> SamFolder.align=[bowtie1,bowtie2,bwa_back,bwa_mem]
##     # align_with_bwa_aln -> SamFolder.align=[bowtie1,bowtie2,bwa_back,bwa_mem]
##     # align_with_bwa_mem -> SamFolder.align=[bowtie1,bowtie2,bwa_back,bwa_mem]

##     # 1.  Look in two pipelines for:
##     #     - Same DataObject.
##     #     - Different Modules pointing to that DataObject.
##     #     - Modules have Consequence SET_TO to different values in the
##     #       same attribute.
##     #     - That attribute is TYPE_ENUM in the DataObject.
##     # 2.  Keep the pipeline with the Module that generates the value
##     #     that comes first alphabetically.
##     # 3.  Print a message showing which value was chosen.

##     if not paths:
##         return []

##     # List of (datatype name, attr name, kept value, deleted value)
##     deleted = []
##     paths = paths[:]

##     # For each data node, list the alternate attributes.
##     dataid2alts = {}
##     for i in range(len(network.nodes)):
##         if not isinstance(network.nodes[i], DataNode):
##             continue
##         x = _list_alternate_attributes1(network, i, nodeid2parents)
##         dataid2alts[i] = x

##     # Only check pathways that has a data node in dataid2alts.
##     has_alts = []
##     for i, p in enumerate(paths):
##         x = [x for x in p.node_ids if dataid2alts.get(x, [])]
##         if x:
##             has_alts.append(i)
##     has_alts = {}.fromkeys(has_alts)

##     ## Convert node_ids to bitwise numbers for faster superset
##     ## comparison.
##     #paths_node_ids_bit = [_intlist2bits(x[0]) for x in paths]

##     seen = {}  # (path_i, path_j) -> 1
##     while True:
##         found = None
##         # Actually more efficient to restart this loop after deleting
##         # paths than to run through it without deleting anything.
##         # After deleting paths, leads to overall fewer calls to
##         # _find_alternate_attributes.  O(N*N) growth.
##         for i in range(len(paths)-1):
##             if not has_alts.get(i):
##                 continue
##             for j in range(i+1, len(paths)):
##                 if not has_alts.get(j):
##                     continue
##                 if (i, j) in seen:
##                     continue
##                 # Almost all time spent in this function.
##                 x = _find_alternate_attributes(
##                     network, paths[i], paths[j], dataid2alts)
##                 seen[(i, j)] = 1
##                 if x:
##                     found = i, j, x
##                     break
##             if found:
##                 break
##         if not found:
##             break

##         path_1, path_2, x = found
##         attr_name, module_id_1, attr_value_1, \
##                    module_id_2, attr_value_2, data_id = x
##         assert attr_value_1 != attr_value_2
##         datatype_name = network.nodes[data_id].datatype.name

##         # Figure out which path to delete.
##         path_to_delete = None
##         # Look in custom_attributes to see if one of them is requested
##         # by the user.
##         user_request = None  # value requested by user
##         for x in custom_attributes:
##             if x.datatype.name != datatype_name:
##                 continue
##             if x.name != attr_name:
##                 continue
##             assert user_request is None, "Multiple custom_attributes."
##             user_request = x.value
##         deleted_based_on_user = True
##         if user_request:
##             if user_request == attr_value_2:
##                 path_to_delete = 1
##             elif user_request == attr_value_1:
##                 path_to_delete = 2
##         if path_to_delete is None:
##             # Keep the one that comes first alphabetically.
##             path_to_delete = 1
##             if attr_value_1 < attr_value_2:
##                 path_to_delete = 2
##             deleted_based_on_user = False

##         # Delete the appropriate path.
##         assert path_to_delete in [1, 2]
##         if path_to_delete == 1:
##             del paths[path_1]
##             #del paths_node_ids_bit[path_1]
##             p1, p2 = path_2, path_1
##             v1, v2 = attr_value_2, attr_value_1
##         else:
##             del paths[path_2]
##             #del paths_node_ids_bit[path_2]
##             p1, p2 = path_1, path_2
##             v1, v2 = attr_value_1, attr_value_2

##         # Fix the indexes of the seen variable.
##         x = seen.keys()
##         x = _fix_node_id_pairs_after_merge(x, p1, p2)
##         seen = {}.fromkeys(x)

##         x = datatype_name, attr_name, v1, v2, deleted_based_on_user
##         deleted.append(x)

##     #_print_attributes_pruned(deleted)
##     return paths


def _list_alternate_attributes1(network, data_id, nodeid2parents):
    # Return list of (data_id, parent_ids, attr_name, attr_values).
    # align_with_bowtie1 -> SamFolder (bowtie1)
    # align_with_bowtie2 -> SamFolder (bowtie2)
    # align_with_bwa_aln -> SamFolder (bwa_backtrack)
    # align_with_bwa_mem -> SamFolder (bwa_mem)
    if data_id not in nodeid2parents:
        return []
    prev_ids = nodeid2parents[data_id]

    # 2 parents must have Consequences SET_TO to the same
    # attribute.
    attr2values = {}     # name -> list of values to set to
    attr2moduleids = {}  # name -> list of module_ids
    for prev_id in prev_ids:
        module_node = network.nodes[prev_id]
        for cons in module_node.consequences:
            if cons.behavior != SET_TO:
                continue
            n, v = cons.name, cons.arg1
            assert v is not None
            if n not in attr2values:
                attr2values[n] = []
                attr2moduleids[n] = []
            if v not in attr2values[n]:
                attr2values[n].append(v)
                attr2moduleids[n].append(prev_id)
    # List of attributes with at least 2 values.
    attrs = [n for (n, v) in attr2values.iteritems() if len(v) >= 2]
    retvals = []
    for attr_name in attrs:
        x = data_id, attr2moduleids, attr_name, attr2values[attr_name]
        retvals.append(x)
    return retvals


def _find_alternate_attributes(network, path_1, path_2, dataid2alts):
    #shared_ids = [x for x in path_1.node_ids if x in path_2.node_ids]
    #unique_ids_1 = [x for x in path_1.node_ids if x not in path_2.node_ids]
    #unique_ids_2 = [x for x in path_2.node_ids if x not in path_1.node_ids]
    shared_ids = path_1.node_ids.intersection(path_2.node_ids)
    unique_ids_1 = path_1.node_ids.difference(path_2.node_ids)
    unique_ids_2 = path_2.node_ids.difference(path_1.node_ids)

    # Look for a DataNode with alternates that is shared in both pathways.
    data_ids = [x for x in shared_ids if dataid2alts.get(x, [])]
    if not data_ids:
        return None

    # Look for ModuleNodes that are unique for each pathway.
    module_ids_1 = [x for x in unique_ids_1
                    if isinstance(network.nodes[x], ModuleNode)]
    module_ids_2 = [x for x in unique_ids_2
                    if isinstance(network.nodes[x], ModuleNode)]
    # This is slower.
    #module_ids_1 = [x for x in node_ids_1 if x not in shared_ids and
    #                isinstance(network.nodes[x], bie3.ModuleNode)]
    #module_ids_2 = [x for x in node_ids_2 if x not in shared_ids and
    #                isinstance(network.nodes[x], bie3.ModuleNode)]
    if not module_ids_1 or not module_ids_2:
        return None

    # Check each data_id carefully.
    for data_id in data_ids:
        for x in dataid2alts[data_id]:
            x, parent_ids, attr_name, attr_values = x

            # Look for parent_ids that are in different pathways.
            good_1 = [x for x in parent_ids if x in module_ids_1]
            good_2 = [x for x in parent_ids if x in module_ids_2]
            if not good_1 or not good_2:
                continue

            # Found one!
            module_id_1 = good_1[0]
            module_id_2 = good_2[0]
            i = parent_ids.index(module_id_1)
            j = parent_ids.index(module_id_2)
            attr_value_1 = attr_values[i]
            attr_value_2 = attr_values[j]
            x = attr_name, module_id_1, attr_value_1, \
                module_id_2, attr_value_2, data_id
            return x
    return None


def _prune_alternate_attributes2(
    network, custom_attributes, paths, nodeid2parents):
    # If a module takes DataNodes that set an attribute to different
    # values, choose one value and don't calculate the other.
    #
    # Fastq.trimmed=no                              -> align
    # Fastq.trimmed=no -> trim -> Fastq.trimmed=yes -> align
    #
    # 1.  Look in two pipelines for:
    #     - Same ModuleNode.
    #     - Different DataNodes can go into that Node.
    #     - One DataNode is upstream of the other DataNode.
    # 2.  Keep the pipeline that is shorter (or found in
    #     custom_attributes).
    # 3.  Print a message showing which value was chosen.

    if not paths:
        return []

    #path_ids = {}
    #transitions = {}
    #for path in paths:
    #    path_ids.update(path.node_ids)
    #    for node_id, next_ids in path.transitions.iteritems():
    #        x = transitions.get(node_id, set()).union(next_ids)
    #        transitions[node_id] = x
    x = _merge_paths(paths)
    path_ids, transitions, x = x
    ancestors = _make_ancestor_dict(network)

    # Find module nodes that can take DataNodes with different
    # attributes.
    alternates = []
    for i in range(len(network.nodes)):
        if not isinstance(network.nodes[i], ModuleNode):
            continue
        x = _list_alternate_attributes2(
            network, i, custom_attributes, path_ids, transitions,
            nodeid2parents)
        alternates.extend(x)

    # For each of the alternates, select the desired transition and
    # rule out the others.
    desired_alternates = [None] * len(alternates)
    for i, x in enumerate(alternates):
        module_id, parent_ids, attr_name, attr_values = x
        # If there is a value specified in the custom attribute, use
        # that one.
        values = _get_custom_values(
            custom_attributes, network.nodes[parent_ids[0]].datatype.name,
            attr_name)
        for value in values:
            if value in attr_values:
                desired_alternates[i] = attr_values.index(value)
                break
        if desired_alternates[i] is not None:
            continue
        # Hack: If the options are "no" and "yes", choose "no".
        if sorted(attr_values) == ["no", "yes"]:
            desired_alternates[i] = attr_values.index("no")
        if desired_alternates[i] is not None:
            continue
        # If one is an ancestor of the other, choose the ancestor.
        if len(parent_ids) == 2:
            if parent_ids[0] in ancestors.get(parent_ids[1], []):
                desired_alternates[i] = 0
            elif parent_ids[1] in ancestors.get(parent_ids[0], []):
                desired_alternates[i] = 1
        if desired_alternates[i] is not None:
            continue
        # Choose the one that comes first in the alphabet.
        x = sorted(attr_values)[0]
        desired_alternates[i] = attr_values.index(x)

    # Make a list of the transitions to avoid in the network.
    to_prune = {}  # transitions to avoid
    for i, x in enumerate(alternates):
        module_id, parent_ids, attr_name, attr_values = x
        for j, id_ in enumerate(parent_ids):
            if j == desired_alternates[i]:
                continue
            to_prune[(parent_ids[j], module_id)] = 1

    # Find the paths to prune.
    delete = {}
    for i, path in enumerate(paths):
        p = False
        for parent_id, child_id in to_prune:
            if child_id in path.transitions.get(parent_id, []):
                p = True
                break
        if p:
            delete[i] = 1

    paths = [x for (i, x) in enumerate(paths) if i not in delete]
    return paths


def _list_alternate_attributes2(
    network, module_id, custom_attributes, path_ids, transitions,
    nodeid2parents):
    # Return list of (module_id, parent_ids, attr_name, attr_values).
    # attr_values can be ATOM or ENUM.
    import itertools

    # Fastq.trimmed=no                              -> align
    # Fastq.trimmed=no -> trim -> Fastq.trimmed=yes -> align

    # At least 2 DataNodes of the same DataType must transition
    # into this ModuleNode.
    x = nodeid2parents[module_id]
    x = [x for x in x if x in path_ids]
    if len(x) < 2:
        return []
    x = [x for x in x if module_id in transitions.get(x, [])]
    if len(x) < 2:
        return []
    parent_ids = x

    # Keep only if there are at least two nodes with the same datatype.
    datatype_names = [network.nodes[x].datatype.name for x in parent_ids]
    counts = {}
    for n in datatype_names:
        counts[n] = counts.get(n, 0) + 1
    x = [x for (i, x) in enumerate(parent_ids)
         if counts[datatype_names[i]] >= 2]
    parent_ids = x
    if len(parent_ids) < 2:
        return []

    # Previous DataNodes must be alternates (i.e. be the same input in
    # different combinations.
    inputnum2parentids = {}  # which input to module -> list of parent IDs
    combos = _bc_to_input_ids(
        network, module_id, custom_attributes, nodeid2parents=nodeid2parents)
    for x in itertools.product(combos, parent_ids):
        combo, parentid = x
        if parentid not in combo:
            continue
        input_num = combo.index(parentid)
        if input_num not in inputnum2parentids:
            inputnum2parentids[input_num] = []
        inputnum2parentids[input_num].append(parentid)

    # Previous DataNodes must have different attribute values.
    alt_attributes = []
    for input_num, parent_ids in inputnum2parentids.iteritems():
        if len(parent_ids) < 2:
            continue
        id_1 = parent_ids[0]
        node_1 = network.nodes[id_1]
        attr_names = []
        for id_2 in parent_ids[1:]:
            node_2 = network.nodes[id_2]
            x, x, diff_attrs = _score_same_data(node_1, node_2)
            # How to handle multiple different attributes?
            if len(diff_attrs) != 1:
                continue
            attr_names.append(diff_attrs[0])
        # Must have all the same different attributes.
        if len(attr_names) != len(parent_ids)-1:
            continue
        all_same = True
        for i in range(1, len(attr_names)):
            if attr_names[i] != attr_names[0]:
                all_same = False
                break
        if not all_same:
            continue
        name = attr_names[0]

        # Ignore this attribute if it is SAME_AS_CONSTRAINT.  Module
        # does not alter the value of the attribute.
        # Example??
        module_passes_attr = False
        for cons in network.nodes[module_id].consequences:
            if cons.name == name and \
                   cons.behavior in [SAME_AS_CONSTRAINT]:
                module_passes_attr = True
                break
        if module_passes_attr:
            continue

        # Do not prune DataNodes whose values are BASED_ON_DATA.  It's
        # not for this module to choose the value.
        # FastQ -> is_compressed -> FastQ.comp=gz -> merge_reads
        #                        -> FastQ.comp=no ->
        based_on_data = False
        for pid in parent_ids:
            # Get the module that generated this DataNode.
            x = nodeid2parents.get(pid, [])
            x = [x for x in x if x in path_ids]
            # No module generated this DataNode.  Must be leaf.  Thus,
            # not BASED_ON_DATA.
            if not x:
                continue
            for mid in x:
                for cons in network.nodes[mid].consequences:
                    if cons.name == name and \
                           cons.behavior in [BASED_ON_DATA]:
                        based_on_data = True
                        break
        if based_on_data:
            continue

        values = [network.nodes[x].attributes[name] for x in parent_ids]
        x = module_id, parent_ids, name, values
        alt_attributes.append(x)
    return alt_attributes


def _prune_superset_pipelines(network, paths, nodeid2parents):
    # Remove pipelines that are just supersets of another pipeline.
    import itertools

    path2length = [len(x.node_ids) for x in paths]

    # Convert node_ids to bitwise numbers for faster superset
    # comparison.
    paths_node_ids_bit = [_intlist2bits(x.node_ids) for x in paths]

    superset = []
    for (i, j) in itertools.product(range(len(paths)), range(len(paths))):
        if i == j:
            continue

        # Optimization: Check for length here to avoid function call.
        if path2length[i] <= path2length[j]:
            continue
        # Optimization: Check for superset here.
        #if not set(node_ids_1).issuperset(node_ids_2):
        b1, b2 = paths_node_ids_bit[i], paths_node_ids_bit[j]
        if (b1 | b2) != b1:
            continue
        if _is_superset_pipeline(network, paths[i], paths[j], nodeid2parents):
            superset.append(i)
    superset = {}.fromkeys(superset)
    paths = [x for (i, x) in enumerate(paths) if i not in superset]
    return paths


def _is_superset_pipeline(network, path_1, path_2, nodeid2parents):
    # Test if path_1 is superset, given path_2.  I.e. path_1 has more
    # processing steps that are not necessary in path_2.

    # If they're superset, then the start nodes must be different.
    # Actually, this is not true.  They can have the same start nodes,
    # but different internal paths.
    #if sorted(start_ids_1) == sorted(start_ids_2):
    #    return False

    # Comparisons now done in calling function.
    ## If they're superset, path_1 must be a superset of path_2.
    #if len(node_ids_1) <= len(node_ids_2):
    #    # Returns from here ~60% of the time.
    #    return False
    ## Takes about 50% of the time in this function.
    #if not set(node_ids_1).issuperset(node_ids_2):
    #    return False

    # At least one of the extra nodes must be a start node.
    # Actually, this is not true.  They can have the same start nodes,
    # but different internal paths.
    #x = list(set(node_ids_1).difference(node_ids_2))
    #x = [x for x in x if x in start_ids_1]
    #super_start_ids_1 = x
    #if not super_start_ids_1:
    #    return False

    # All paths starting from start_ids_1 must go through one of
    # start_ids_2.
    start_ids_1 = [x for x in path_1.start_ids if x is not None]
    start_ids_2 = [x for x in path_2.start_ids if x is not None]
    for start_id in start_ids_1:
        if not _does_path_go_through(
            network, start_id, path_1.node_ids, start_ids_2):
            return False

    # Is a superset:
    #   is_compressed -> Fastq (yes) -> uncompress -> Fastq (no) -> merge
    #   is_compressed -> Fastq (yes)                             -> merge
    # Looks like a superset, but it's not:
    #   is_compressed -> Fastq (yes) -> uncompress -> Fastq (no)
    #   is_compressed                              -> Fastq (no)
    super_ids = [x for x in path_1.node_ids if x not in path_2.node_ids]
    assert super_ids
    # If any of the nodes in the super_ids list is downstream of a
    # Module with a BASED_ON_DATA Consequence, then this is not a
    # superset.
    for node_id in super_ids:
        if not isinstance(network.nodes[node_id], DataNode):
            continue
        x = nodeid2parents.get(node_id, [])
        x = [x for x in x if x in path_1.node_ids]
        x = [x for x in x if node_id in path_1.transitions[x]]
        if not x:
            continue
        for mid in x:
            node = network.nodes[mid]
            x = [x for x in node.consequences if x.behavior == BASED_ON_DATA]
            if x:
                return False

    #for (node_id, node) in enumerate(network.nodes):
    #    if node_id in super_ids: # ignore if this is part of the superset
    #        continue
    #    if not isinstance(node, ModuleNode):
    #        continue
    #    x = [x for x in node.consequences if x.behavior == BASED_ON_DATA]
    #    if not x:
    #        continue
    #    next_ids = network.transitions.get(node_id, [])
    #    x = [x for x in next_ids if x in super_ids]
    #    if x:
    #        return False
    return True


def _does_path_go_through(network, start_id, good_ids, intermediate_ids):
    # Do all paths starting from start_id go through intermediate_ids?
    stack = [start_id]
    while stack:
        node_id = stack.pop(0)
        if node_id in intermediate_ids:
            continue
        x = network.transitions.get(node_id, [])
        next_ids = [x for x in x if x in good_ids]
        if not next_ids:
            # Goes to end of network.  Did not encounter intermediate_ids.
            return False
        stack.extend(next_ids)
    return True

def _prune_parallel_pipelines(network, paths, nodeid2parents):
    # Remove pipelines that are parallel to another pipeline.
    global SUBPATH_CACHE   # for optimization
    SUBPATH_CACHE = {}

    #plot_pipelines(
    #    "pipeline", network, paths, {}, max_pipelines=16, verbose=True)

    path_lengths = [len(x.node_ids) for x in paths]

    # Prune parallel pipelines with the same lengths first to get a
    # smaller set of pipelines.  Then, do the more complex pruning.
    prune = {}
    for i in range(len(paths)-1):
        for j in range(i+1, len(paths)):
            if i in prune and j in prune:
                continue

            # Optimization: Don't check if the lengths are not the
            # same.
            if path_lengths[i] != path_lengths[j]:
                continue
            p = _is_parallel_pipeline1(
                network, paths[i], paths[j], nodeid2parents)
            if not p:
                p = _is_parallel_pipeline2(
                    network, paths[i], paths[j], nodeid2parents)
            #if not p:
            #    p = _is_parallel_pipeline3(
            #        network, paths, i, j, nodeid2parents)
            if not p:
                continue
            if p == 1:
                prune[i] = 1
            else:
                prune[j] = 1
    #orig_path_ids = [
    #    x for (i, x) in enumerate(range(len(paths))) if i not in prune]
    paths = [x for (i, x) in enumerate(paths) if i not in prune]
    #path_lengths = [x for (i, x) in enumerate(path_lengths) if i not in prune]

    # Optimization: Cache the modules with only one child.
    modules_with_one_child = []
    for node_id in range(len(network.nodes)):
        if not isinstance(network.nodes[node_id], ModuleNode):
            continue
        x = network.transitions.get(node_id, [])
        if len(x) == 1:
            modules_with_one_child.append(node_id)
    modules_with_one_child = frozenset(modules_with_one_child)

    # Optimization: Study the pathways to minimize the number of nodes
    # that need to be checked.
    # Make a list of the nodes that are shared across all pathways.
    # These are not going to be parallel.
    node2counts = {}
    for p in paths:
        for nid in p.node_ids:
            node2counts[nid] = node2counts.get(nid, 0) + 1
    shared_node_ids = []
    for nid, count in node2counts.iteritems():
        if count == len(paths):
            shared_node_ids.append(nid)
    shared_node_ids = frozenset(shared_node_ids)
    all_node_ids = frozenset(node2counts.keys())

    # Make a list of the node_ids in the paths that aren't shared.
    path2nodeids = []
    for p in paths:
        x = p.node_ids.difference(shared_node_ids)
        path2nodeids.append(x)

    # Group the node_ids that aren't shared into clusters.  Since I'm
    # looking for parallel pipelines, each parallel portion should be
    # unconnected in the network.  Group together the nodes in each
    # unconnected part.  Only compare two pathways if they differ by
    # at most one unconnected part.
    
    # Start with each node_id in a separate cluster.  Then try to join
    # the clusters together based on transitions.
    transitions = {}
    for p in paths:
        for nid, nextids in p.transitions.iteritems():
            if nid not in transitions:
                transitions[nid] = []
            x = nextids.union(transitions[nid])
            transitions[nid] = x
    clusters = []
    for nid in all_node_ids.difference(shared_node_ids):
        clusters.append([nid])
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(clusters)-1:
            j = i+1
            while j < len(clusters):
                merge = False
                # If any of clusters[i] transitions to clusters[j],
                # then merge.
                trans1 = []
                for nid1 in clusters[i]:
                    trans1.extend(transitions.get(nid1, []))
                trans2 = []
                for nid2 in clusters[j]:
                    trans2.extend(transitions.get(nid2, []))
                    
                # If any of clusters[i] transitions to clusters[j],
                # then merge.
                if set(trans1).intersection(clusters[j]) or \
                   set(trans2).intersection(clusters[i]):
                    clusters[i] = clusters[i] + clusters[j]
                    del clusters[j]
                    changed = True
                else:
                    j += 1
            i += 1

    # Each cluster should have a distinct set of nodes.  Choose an
    # exemplar for each cluster.  Arbitrarily, just use the lowest
    # node ID.
    exemplars = frozenset([min(x) for x in clusters])
                    
    # Now do the more computationally expensive pruning.
    prune = {}
    for i in range(len(paths)-1):
        for j in range(i+1, len(paths)):
            if i in prune and j in prune:
                continue
            # If two paths differ by more than 1 cluster, then they
            # cannot be parallel.
            c1 = exemplars.intersection(path2nodeids[i])
            c2 = exemplars.intersection(path2nodeids[j])
            x1 = c1.difference(c2)
            x2 = c2.difference(c1)
            if len(x1) != 1 or len(x2) != 1:
                continue
            
            p = _is_parallel_pipeline3(
                network, paths, i, j, nodeid2parents, modules_with_one_child,
                path2nodeids)
            if not p:
                continue
            if p == 1:
                prune[i] = 1
            else:
                prune[j] = 1

    paths = [x for (i, x) in enumerate(paths) if i not in prune]
    return paths


def _is_parallel_pipeline1(network, path_1, path_2, nodeid2parents):
    # Test if path_1 is parallel to path_2.  If not parallel, returns
    # False.  Otherwise, returns a 1 or 2 indicating which one should
    # be pruned.
    #
    # Parallel:
    # BAM (no) -> sort by name  -> BAM (name)  -> count_htseq
    # BAM (no) -> sort by coord -> BAM (coord) -> count_htseq
    # Different from:
    # FASTQ (unknown) -> is_compress                              -> FASTQ (no)
    # FASTQ (unknown) -> is_compress -> FASTQ (yes) -> uncompress -> FASTQ (no)

    #node_ids_1, start_ids_1, data_indexes_1 = path_1
    #node_ids_2, start_ids_2, data_indexes_2 = path_2

    # If they're parallel, then they have the same number of nodes.
    if len(path_1.node_ids) != len(path_2.node_ids):
        return False
    # They differ by only two nodes.
    unique_ids_1 = [x for x in path_1.node_ids if x not in path_2.node_ids]
    if len(unique_ids_1) != 2:
        return False
    unique_ids_2 = [x for x in path_2.node_ids if x not in path_1.node_ids]
    if len(unique_ids_2) != 2:
        return False
    # The module node is upstream of the data node.
    module_id_1, data_id_1 = unique_ids_1
    module_id_2, data_id_2 = unique_ids_2
    if isinstance(network.nodes[module_id_1], DataNode):
        module_id_1, data_id_1 = data_id_1, module_id_1
    if isinstance(network.nodes[module_id_2], DataNode):
        module_id_2, data_id_2 = data_id_2, module_id_2
    if not isinstance(network.nodes[module_id_1], ModuleNode):
        return False
    if not isinstance(network.nodes[data_id_1], DataNode):
        return False
    if not isinstance(network.nodes[module_id_2], ModuleNode):
        return False
    if not isinstance(network.nodes[data_id_2], DataNode):
        return False
    if data_id_1 not in network.transitions.get(module_id_1, []):
        return False
    if data_id_2 not in network.transitions.get(module_id_2, []):
        return False
    # The module nodes have the same parent.
    parent_ids_1 = nodeid2parents.get(module_id_1, [])
    parent_ids_2 = nodeid2parents.get(module_id_2, [])
    common_ids = [x for x in parent_ids_1 if x in parent_ids_2]
    if not common_ids:
        return False
    # The data nodes have the same child.
    child_ids_1 = network.transitions.get(data_id_1, [])
    child_ids_2 = network.transitions.get(data_id_2, [])
    common_ids = [x for x in child_ids_1 if x in child_ids_2]
    if not common_ids:
        return False
    # The data nodes have the same type.
    data_1 = network.nodes[data_id_1]
    data_2 = network.nodes[data_id_2]
    if data_1.datatype != data_2.datatype:
        return False
    # The data nodes differ by only one attribute.
    assert sorted(data_1.attributes) == sorted(data_2.attributes)
    diff = [x for x in data_1.attributes
            if data_1.attributes[x] != data_2.attributes[x]]
    if len(diff) != 1:
        return False
    # Prune the one whose value comes later in the alphabet.
    name = diff[0]
    value_1 = data_1.attributes[name]
    value_2 = data_2.attributes[name]
    if value_1 < value_2:
        return 2
    return 1


def _is_parallel_pipeline2(network, path_1, path_2, nodeid2parents):
    # Test if path_1 is parallel to path_2.  If not parallel, returns
    # False.  Otherwise, returns a 1 or 2 indicating which one should
    # be pruned.
    # Does a fast and sloppy check.
    #
    # Parallel:
    # BAM (no) -> sort -> BAM (sort=y) -> addgroup ->  BAM (sort=y, group=y)
    #          -> addgroup -> BAM (group=y) -> sort ->

    # If they're parallel, then they have the same number of nodes.
    if len(path_1.node_ids) != len(path_2.node_ids):
        return False
    # They differ by only three nodes.
    unique_ids_1 = [x for x in path_1.node_ids if x not in path_2.node_ids]
    if len(unique_ids_1) != 3:
        return False
    unique_ids_2 = [x for x in path_2.node_ids if x not in path_1.node_ids]
    if len(unique_ids_2) != 3:
        return False

    # The module node is upstream of the data node.
    top_id_1, data_id_1, bottom_id_1 = unique_ids_1
    top_id_2, data_id_2, bottom_id_2 = unique_ids_2
    if isinstance(network.nodes[top_id_1], DataNode):
        top_id_1, data_id_1 = data_id_1, top_id_1
    if isinstance(network.nodes[bottom_id_1], DataNode):
        bottom_id_1, data_id_1 = data_id_1, bottom_id_1
    if data_id_1 in network.transitions[bottom_id_1]:
        top_id_1, bottom_id_1 = bottom_id_1, top_id_1
    if isinstance(network.nodes[top_id_2], DataNode):
        top_id_2, data_id_2 = data_id_2, top_id_2
    if isinstance(network.nodes[bottom_id_2], DataNode):
        bottom_id_2, data_id_2 = data_id_2, bottom_id_2
    if data_id_2 in network.transitions[bottom_id_2]:
        top_id_2, bottom_id_2 = bottom_id_2, top_id_2

    # Make sure the node types are correct.
    if not isinstance(network.nodes[top_id_1], ModuleNode):
        return False
    if not isinstance(network.nodes[data_id_1], DataNode):
        return False
    if not isinstance(network.nodes[bottom_id_1], ModuleNode):
        return False
    if not isinstance(network.nodes[top_id_2], ModuleNode):
        return False
    if not isinstance(network.nodes[data_id_2], DataNode):
        return False
    if not isinstance(network.nodes[bottom_id_2], ModuleNode):
        return False

    # Make sure the connections are correct.
    if data_id_1 not in path_1.transitions.get(top_id_1, []):
        return False
    if bottom_id_1 not in path_1.transitions.get(data_id_1, []):
        return False
    if data_id_2 not in path_2.transitions.get(top_id_2, []):
        return False
    if bottom_id_2 not in path_2.transitions.get(data_id_2, []):
        return False

    # The top nodes have the same parent.
    parent_ids_1 = nodeid2parents.get(top_id_1, [])
    parent_ids_2 = nodeid2parents.get(top_id_2, [])
    x = [x for x in parent_ids_1 if x in parent_ids_2]
    x = [x for x in x if x in path_1.node_ids and x in path_2.node_ids]
    common_ids = x
    if not common_ids:
        return False

    # The bottom nodes have the same child.
    child_ids_1 = network.transitions.get(bottom_id_1, [])
    child_ids_2 = network.transitions.get(bottom_id_2, [])
    x = [x for x in child_ids_1 if x in child_ids_2]
    x = [x for x in x if x in path_1.node_ids and x in path_2.node_ids]
    common_ids = x
    if not common_ids:
        return False

    # The data nodes have the same type.
    data_1 = network.nodes[data_id_1]
    data_2 = network.nodes[data_id_2]
    if data_1.datatype != data_2.datatype:
        return False

    # The data nodes differ by two attributes.
    assert sorted(data_1.attributes) == sorted(data_2.attributes)
    diff = [x for x in data_1.attributes
            if data_1.attributes[x] != data_2.attributes[x]]
    if len(diff) != 2:
        return False

    # Prune the one whose value comes later in the alphabet.
    name = sorted(diff)[0]
    value_1 = data_1.attributes[name]
    value_2 = data_2.attributes[name]
    if value_1 < value_2:
        return 2
    return 1


def _is_parallel_pipeline3(
    network, paths, path_id_1, path_id_2, nodeid2parents,
    modules_with_one_child, path2nodeids):
    # Test if path_1 is parallel to path_2.  If not parallel, returns
    # False.  Otherwise, returns a 1 or 2 indicating which one should
    # be pruned.
    # More careful and slower version of _is_parallel_pipeline2.
    #
    # Parallel:
    # DataNode -> sort_coord -> mark_dup -> add_read_group -> DataNode
    #          -> add_read_group -> sort_coord -> mark_dup ->
    path_1, path_2 = paths[path_id_1], paths[path_id_2]
    ids_1, ids_2 = path2nodeids[path_id_1], path2nodeids[path_id_2]
    #assert path_1.start_ids == path_2.start_ids, "%s %s" % (
    #    path_1.start_ids, path_2.start_ids)
    if path_1.start_ids != path_2.start_ids:
        return False
    x = _compare_paths(
        network, path_1, path_2, ids_1, ids_2, modules_with_one_child)
    shared_ids, unique_ids_1, unique_ids_2 = x
    # shared_ids is almost always a very long list.

    # Optimization: The same set of unique_ids are tested over and
    # over again.  With one network, this was run 907,625 times on
    # only 767 unique sets of unique_ids (and transitions).  Don't
    # need to re-test the same ones over and over again.
    # Actually, this wouldn't save much because most of the time is
    # spent in _compare_paths.  But as other things optimize, the time
    # for this is relatively increasing.
    
    #for unique_ids in [unique_ids_1, unique_ids_2]:
    #    h = list(unique_ids)
    #    for nid1 in unique_ids:
    #        for nid2 in path_1.transitions.get(nid1, []):
    #            x = (nid1, nid2)
    #            h.append(x)

    if not unique_ids_1 or not unique_ids_2:
        return False

    # For downstream tests, calculate some useful variables describing
    # the networks.
    # _build_subpath is relatively fast.  No use to optimize.
    subpath_1 = _build_subpath(
        network, paths, path_id_1, unique_ids_1, nodeid2parents)
    subpath_2 = _build_subpath(
        network, paths, path_id_2, unique_ids_2, nodeid2parents)

    # There should only be 1 top and 1 bottom node.  This test removes
    # 95% of the possibilities.  The others don't filter as much.
    if len(subpath_1.top_node_ids) != 1 or len(subpath_1.bottom_node_ids) != 1:
        return False
    if len(subpath_2.top_node_ids) != 1 or len(subpath_2.bottom_node_ids) != 1:
        return False

    # The top should be modules.  The bottom might be a DataNode if
    # the last step is the same.
    # path1 -> BamFolder.sorted=name  -> sort_by_contig
    # path2 -> BamFolder.sorted=coord -> sort_by_contig
    x = [subpath_1.top_node_ids[0], subpath_2.top_node_ids[0]]
    for x in x:
        if not isinstance(network.nodes[x], ModuleNode):
            return False

    # The parents of top and children of bottom should be the same.
    if sorted(subpath_1.parents_of_top) != sorted(subpath_2.parents_of_top):
        return False
    if sorted(subpath_1.children_of_bottom) != \
           sorted(subpath_2.children_of_bottom):
        return False

    # Any attributes based on data should be the same.
    bod1 = subpath_1.based_on_data
    bod2 = subpath_2.based_on_data
    if sorted(bod1) != sorted(bod2):
        return False
    for name in bod1:
        if bod1[name] != bod2[name]:
            return False

    # Be sure to always take the same parallel pipeline.  Use the one
    # whose module name comes first in the alphabet.
    id_1, id_2 = subpath_1.top_node_ids[0], subpath_2.top_node_ids[0]
    node_1, node_2 = network.nodes[id_1], network.nodes[id_2]
    value_1, value_2 = node_1.name, node_2.name
    # Prune the one whose value comes later in the alphabet.
    if value_1 < value_2:
        return 2
    return 1


def _compare_paths(network, path_1, path_2, node_ids_1, node_ids_2,
                   modules_with_one_child):
    # This function takes a long time and is called frequently.  Good
    # target for optimization.
    
    # If a node occurs in only one pathway, then it is unique for
    # sure.
    #shared_ids = path_1.node_ids.intersection(path_2.node_ids)
    #unique_ids_1 = path_1.node_ids.difference(path_2.node_ids)
    #unique_ids_2 = path_2.node_ids.difference(path_1.node_ids)
    shared_ids = node_ids_1.intersection(node_ids_2)
    unique_ids_1 = node_ids_1.difference(node_ids_2)
    unique_ids_2 = node_ids_2.difference(node_ids_1)
    shared_ids = set(shared_ids)
    unique_ids_1 = set(unique_ids_1)
    unique_ids_2 = set(unique_ids_2)
    ## To Do (optimization): Should change this to use set operations.
    #unique_ids_1 = _dict_diff(path_1.node_ids, path_2.node_ids, _as_dict=True)
    #unique_ids_2 = _dict_diff(path_2.node_ids, path_1.node_ids, _as_dict=True)
    #shared_ids = _dict_diff(path_1.node_ids, unique_ids_1, _as_dict=True)

    # If a Module occurs in both pathways, it might be unique if it
    # has transitions to different children.
    # sort_by_contig -> BamFolder.readgroups (no)
    # sort_by_contig -> BamFolder.readgroups (yes)
    # i.e same module, but different uses

    # Pre-caching modules_with_one_child leads to 40% speedup.
    #for node_id in list(shared_ids):
    for node_id in shared_ids.intersection(modules_with_one_child):
        if not isinstance(network.nodes[node_id], ModuleNode):
            continue

        # Can happen if the path is not finished yet, and pruning
        # while building pipelines.
        if node_id not in path_1.transitions:
            continue
        if node_id not in path_2.transitions:
            continue

        #assert node_id in path_1.transitions
        #assert node_id in path_2.transitions
        children_1 = path_1.transitions[node_id]
        children_2 = path_2.transitions[node_id]
        # Actually takes more time to compare lengths first.  Sorted
        # pretty fast.
        #if len(children_1) != len(children_2):
        #    continue
        #if sorted(children_1) == sorted(children_2):
        if children_1 == children_2:
            continue
        shared_ids.remove(node_id)
        unique_ids_1.add(node_id)
        unique_ids_2.add(node_id)
        #del shared_ids[node_id]
        #unique_ids_1[node_id] = 1
        #unique_ids_2[node_id] = 1
    shared_ids = frozenset(shared_ids)
    unique_ids_1 = frozenset(unique_ids_1)
    unique_ids_2 = frozenset(unique_ids_2)
    return shared_ids, unique_ids_1, unique_ids_2



class PathSubset:
    # Have a network that includes all nodes.
    # A subset of those nodes comprises a path through the network.
    # A subset of the path comprises a subpath
    #
    # path_node_ids
    # sub_node_ids
    # path_transitions
    #
    # top_node_ids         list of sub node IDs (no parents in subpath)
    # bottom_node_ids      list of sub node IDs (no children in subpath)
    # parents_of_top       list of path node IDs (above top_node_ids)
    # children_of_bottom   list of path node IDs (below top_node_ids)
    #
    # sub_parents      dict of sub_node_id -> ids of parents in subpath
    # sub_children     dict of sub_node_id -> ids of children in subpath
    # path_parents     dict of sub_node_id -> ids of parents not in subpath
    # path_children    dict of sub_node_id -> ids of children not in subpath
    #
    # based_on_data    dict of attrname -> value; from subpath
    #
    # Since this is very expensive computationally, calculate these
    # variables on demand.
    def __init__(
        self, network, path_node_ids, path_transitions, sub_node_ids,
        nodeid2parents):
        # Make sure these are dicts, or will be really slow.
        assert type(path_node_ids) in [type({}), frozenset, set]
        assert type(sub_node_ids) in [type({}), frozenset, set]

        self._network = network
        self._nodeid2parents = nodeid2parents
        self.path_node_ids = path_node_ids
        self.path_transitions = path_transitions
        self.sub_node_ids = sub_node_ids
    def __getattr__(self, name):
        var2info = {
            "sub_parents" : (
                "_sub_parents", self._calc_sub_parents_and_top_nodes),
            "sub_children" : (
                "_sub_children", self._calc_sub_children_and_bottom_nodes),
            "path_parents" : ("_path_parents", self._calc_path_parents),
            "path_children" : ("_path_children", self._calc_path_children),
            "top_node_ids" : (
                "_top_node_ids", self._calc_sub_parents_and_top_nodes),
            "bottom_node_ids" : (
                "_bottom_node_ids", self._calc_sub_children_and_bottom_nodes),
            "parents_of_top" : ("_parents_of_top", self._calc_parents_of_top),
            "children_of_bottom" : (
                "_children_of_bottom", self._calc_children_of_bottom),
            "based_on_data" : ("_based_on_data", self._calc_based_on_data),
            }
        if name not in var2info:
            raise AttributeError, name
        x = var2info[name]
        member_name, fn = x
        if not hasattr(self, member_name):
            fn()
            assert hasattr(self, member_name)
            #setattr(self, member_name, fn())
        return getattr(self, member_name)
    def _calc_sub_parents_and_top_nodes(self):
        pars = {}
        top_node_ids = []
        for node_id in self.sub_node_ids:
            x = self._nodeid2parents.get(node_id, [])
            x = [x for x in x if x in self.sub_node_ids]
            if x:
                pars[node_id] = x
            else:
                top_node_ids.append(node_id)
        self._sub_parents = pars
        self._top_node_ids = top_node_ids
    def _calc_sub_children_and_bottom_nodes(self):
        children = {}
        bottom_node_ids = []
        for node_id in self.sub_node_ids:
            x = self.path_transitions.get(node_id, [])
            x = [x for x in x if x in self.sub_node_ids]
            if x:
                children[node_id] = x
            else:
                bottom_node_ids.append(node_id)
        self._sub_children = children
        self._bottom_node_ids = bottom_node_ids
    def _calc_path_parents(self):
        pars = {}
        for node_id in self.sub_node_ids:
            x = self._nodeid2parents.get(node_id, [])
            x = [x for x in x if x not in self.sub_node_ids]
            if x:
                pars[node_id] = x
        self._path_parents = pars
    def _calc_path_children(self):
        pars = {}
        for node_id in self.sub_node_ids:
            x = self.path_transitions.get(node_id, [])
            x = [x for x in x if x not in self.sub_node_ids]
            if x:
                pars[node_id] = x
        self._path_children = pars
    #def _calc_top_node_ids(self):
    #    sub_parents = self.sub_parents
    #    x = [x for x in self.sub_node_ids if not sub_parents.get(x)]
    #    self._top_node_ids = x
    #def _calc_bottom_node_ids(self):
    #    sub_children = self.sub_children
    #    x = [x for x in self.sub_node_ids if not sub_children.get(x)]
    #    self._bottom_node_ids = x
    def _calc_parents_of_top(self):
        parents_of_top = {}
        for node_id in self.top_node_ids:
            x = [x for x in self._nodeid2parents.get(node_id, [])]
            # Is this right?
            x = [x for x in x if x in self.path_node_ids]
            for x in x:
                parents_of_top[x] = 1
        self._parents_of_top = parents_of_top.keys()
    def _calc_children_of_bottom(self):
        children_of_bottom = {}
        for node_id in self.bottom_node_ids:
            x = [x for x in self.path_transitions.get(node_id, [])]
            x = [x for x in x if x in self.path_node_ids]
            for x in x:
                children_of_bottom[x] = 1
        self._children_of_bottom = children_of_bottom.keys()
    def _calc_based_on_data(self):
        # For all the nodes in the subpath, get the values of the
        # attributes that are set based on data.
        based_on_data = {}
        for node_id in self.sub_node_ids:
            node = self._network.nodes[node_id]
            if not isinstance(node, ModuleNode):
                continue
            # Make a list of the attributes that are BASED_ON_DATA.
            x = [
                x for x in node.consequences if
                x.behavior == BASED_ON_DATA]
            attr_names = [x.name for x in x]
            if not attr_names:
                continue

            # Find the values of the attributes in the children.
            x = self.path_transitions.get(node_id, [])
            # Is this necessary?
            child_ids = [x for x in x if x in self.path_node_ids]
            for name in attr_names:
                for node_id in child_ids:
                    node = self._network.nodes[node_id]
                    if name not in node.attributes:
                        continue
                    value = node.attributes[name]
                    assert name not in based_on_data
                    based_on_data[name] = value
                    #if name not in based_on_data:
                    #    based_on_data[name] = []
                    #based_on_data[name].append(value)
        self._based_on_data = based_on_data


SUBPATH_CACHE = {}  # (path_id, sub_node_ids) -> PathSubset
def _build_subpath(network, paths, path_id, sub_node_ids, nodeid2parents):
    global SUBPATH_CACHE
    key = path_id, tuple(sorted(sub_node_ids))
    if key not in SUBPATH_CACHE:
        x = _build_subpath_h(
            network, paths, path_id, sub_node_ids, nodeid2parents)
        SUBPATH_CACHE[key] = x
    return SUBPATH_CACHE[key]


def _build_subpath_h(network, paths, path_id, sub_node_ids, nodeid2parents):
    path = paths[path_id]
    return PathSubset(
        network, path.node_ids, path.transitions, sub_node_ids, nodeid2parents)


def get_input_nodes(
    network, custom_attributes, skip_datatypes=None,
    skip_private_datatypes=False, max_inputs=None):
    # Return a list of tuples of node ids.  Each tuple contains a set
    # of node IDs that can serve as the inputs to this network.
    #
    # Example return value:
    #   [(1, 5), (8,)]
    # This means that the set of nodes 1 and 5 would make a valid input
    # to this network.  Or, node 8 by itself would also be a valid
    # input.

    # Should be list of names of datatypes to ignore.
    skip_datatypes = skip_datatypes or []

    data_node_ids = [
        node_id for (node_id, node) in enumerate(network.nodes)
        if isinstance(node, DataNode)]

    # Make list of node_ids to skip.
    skip_ids = [
        x for x in data_node_ids
        if network.nodes[x].datatype.name in skip_datatypes]
    if skip_private_datatypes:
        x = [x for x in data_node_ids
             if network.nodes[x].datatype.name.startswith("_")]
        skip_ids.extend(x)
    skip_ids = {}.fromkeys(skip_ids)

    class nodeid2nodeid_fn:
        def __init__(self, skip_ids):
            self.skip_ids = skip_ids
        def __call__(self, node_id):
            if node_id in self.skip_ids:
                return None
            return node_id
    fn = nodeid2nodeid_fn(skip_ids)

    nodeid_combos = _product_network(
        network, custom_attributes, max_nodes=max_inputs, nodeid2id_fn=fn)
    return nodeid_combos


def get_input_datatypes(
    network, custom_attributes,
    skip_datatypes=None, skip_private_datatypes=False, max_inputs=None):
    # Return a list of tuples of datatypes.  Each tuple contains a set
    # of DataType objects that can serve as the inputs to this
    # network.
    #
    # Example return value:
    #   [(UnprocessedSignalFile, ClassLabelFile), (UnprocessedSignalFile,)]

    # Should be list of names of datatypes to ignore.
    skip_datatypes = skip_datatypes or []

    data_node_ids = [
        node_id for (node_id, node) in enumerate(network.nodes)
        if isinstance(node, DataNode)]

    # Make list of node_ids to skip.
    skip_ids = [
        x for x in data_node_ids
        if network.nodes[x].datatype.name in skip_datatypes]
    if skip_private_datatypes:
        x = [x for x in data_node_ids
             if network.nodes[x].datatype.name.startswith("_")]
        skip_ids.extend(x)
    skip_ids = {}.fromkeys(skip_ids)

    class nodeid2dtid_fn:
        def __init__(self, network, skip_ids):
            # Make a list of all the datatypes in the network.
            # Datatype ID is index into this list.
            id2name = []         # id of datatype -> name of datatype
            name2datatype = {}   # name of datatype -> DataType
            name2id = {}         # name of datatype -> id of datatype
            for node_id in data_node_ids:
                dt = network.nodes[node_id].datatype
                name2datatype[dt.name] = dt
            id2name = sorted(name2datatype)
            for i, dt in enumerate(id2name):
                name2id[dt] = i

            self.network = network
            self.skip_ids = skip_ids
            self.name2datatype = name2datatype
            self.id2name = id2name
            self.name2id = name2id
        def __call__(self, node_id):
            if node_id in self.skip_ids:
                return None
            node = self.network.nodes[node_id]
            return self.name2id[node.datatype.name]
    fn = nodeid2dtid_fn(network, skip_ids)

    # Take the product of the data types across the network.
    dtid_combos = _product_network(
        network, custom_attributes, max_nodes=max_inputs, nodeid2id_fn=fn)

    # Convert from ID back to datatypes.
    datatype_combos = []
    for dtid_combo in dtid_combos:
        names = [fn.id2name[x] for x in dtid_combo]
        x = tuple([fn.name2datatype[x] for x in names])
        datatype_combos.append(x)
    return datatype_combos


def group_nodes_by_datatype(network, inputs):
    # Take a network and inputs (e.g. from get_inputs) and return a
    # dictionary where the keys are tuples of DataTypes and the values
    # are the inputs with that pattern of DataTypes.
    #
    # Example return value:
    # {
    #    (SignalFile,) : [(2,), (4,), (8,)],
    #    (SignalFile, ClassLabelFile) : [(5, 10)],
    # }
    # The order of the nodes may not be preserved.

    data_node_ids = [
        node_id for (node_id, node) in enumerate(network.nodes)
        if isinstance(node, DataNode)]

    # Make a list of all the datatypes in the network.  Datatype ID is
    # index into this list.
    #all_datatypes = []  # id of datatype -> name of datatype
    name2datatype = {}   # name of datatype -> DataType
    datatype2id = {}     # name of datatype -> id of datatype
    for node_id in data_node_ids:
        dt = network.nodes[node_id].datatype
        name2datatype[dt.name] = dt
    all_datatypes = sorted(name2datatype)
    for i, dt in enumerate(all_datatypes):
        datatype2id[dt] = i

    # Figure out the datatype for each node.
    # Really fast.
    nodeid2dtid = {}
    for node_id in data_node_ids:
        dtid = datatype2id[network.nodes[node_id].datatype.name]
        nodeid2dtid[node_id] = dtid

    # Figure out the nodes for each data type.  Very slow.  There's
    # lots of different inputs, all unique.  There are many fewer
    # combinations of datatypes.  Thus, iterating over datatypes will
    # be much faster.
    data = {}  # tuple of datatype ids -> tuple of node IDs -> 1
    for inp in inputs:
        # This line takes ~50% of the time in this function.
        x = [nodeid2dtid[x] for x in inp]
        x = sorted(x)     # Takes ~25% of time.
        dtids = tuple(x)
        if dtids not in data:
            data[dtids] = {}
        data[dtids][inp] = 1

    # Convert the datatype IDs to objects.
    # Relatively fast.
    clean = {}  # tuple of datatype objects -> list of tuple of node IDs
    for dtids, value in data.iteritems():
        names = [all_datatypes[x] for x in dtids]
        datatypes = tuple([name2datatype[x] for x in names])
        clean[datatypes] = value.keys()
    return clean


def summarize_moduledb(moduledb):
    """Take a list of ModuleNodes and return a ModuleDbSummary object."""
    name2module = {}  # module_name -> ModuleNode
    for module in moduledb:
        assert module.name not in name2module
        name2module[module.name] = module
    module_names = sorted(name2module)

    # module_name -> (list of in Datatypes, out Datatype)
    name2datatypes = {}
    for name, module in name2module.iteritems():
        name2datatypes[name] = module.in_datatypes, module.out_datatype

    # All DataType objects found in moduledb.
    datatypes = {}  # name -> DataType object
    for (in_datatypes, out_datatype) in name2datatypes.itervalues():
        dts = in_datatypes + [out_datatype]
        for dt in dts:
            if dt.name in datatypes:
                continue
            datatypes[dt.name] = dt
    datatypes = datatypes.values()

    x = ModuleDbSummary(module_names, name2module, name2datatypes, datatypes)
    return x


def check_moduledb(moduledb):
    seen = {}
    for module in moduledb:
        assert isinstance(module, ModuleNode)
        seen[module.name] = seen.get(module.name, 0) + 1

    # Make sure no duplicate modules.
    dups = ["%s (%d times)" % (x, seen[x]) for x in seen if seen[x] > 1]
    msg = ""
    x = dups
    if len(x) > 5:
        x = x[:5] + ["... plus %s more" % (len(dups) - 5)]
        msg = "\n".join(x)
    assert not dups, "Duplicate modules: %s" % msg


def print_modules(moduledb):
    summary = summarize_moduledb(moduledb)

    for name in summary.module_names:
        in_datatypes, out_datatype = summary.name2datatypes[name]
        in_names = [x.name for x in in_datatypes]
        x = ", ".join(in_names)
        x = name, x, out_datatype.name
        print "\t".join(x)

    print
    for datatype in summary.datatypes:
        print datatype.name
        for attr in sorted(datatype.attributes):
            values = datatype.attributes[attr]
            if type(values) is type(""):
                values = [values]
            print "  %s : %s" % (attr, ", ".join(map(str, values)))
        print


def print_network(network, outhandle=None):
    outhandle = outhandle or sys.stdout
    if type(outhandle) is type(""):
        outhandle = open(outhandle, 'w')
    line_width = 72
    for i in range(len(network.nodes)):
        p_step = "%2d.  " % i
        p_space = " " * (len(p_step) + 2)
        _print_line(str(network.nodes[i]), p_step, p_space, line_width,
                    outhandle=outhandle)
    print >>outhandle

    for i in sorted(network.transitions):
        x = [i, "->"] + network.transitions[i]
        print >>outhandle, "\t".join(map(str, x))


def _format_datanode_gv(node, node_id):
    from genomicode import parselib

    title_size = 32

    node_name = "%s [%d]" % (node.datatype.name, node_id)

    LINE_WIDTH = 60

    lines = []
    w = lines.append
    w('<FONT POINT-SIZE="%d"><B>%s</B></FONT>' % (title_size, node_name))

    if node.attributes:
        w('<BR/>')
        w('<BR/>')
        w("<U>Data Attributes</U>:")
        w('<BR ALIGN="LEFT"/>')
    for name in sorted(node.attributes):
        value = node.attributes[name]
        x = "%s = %s" % (name, value)
        for x in parselib.linesplit(x, prefix1=0, prefixn=4, width=LINE_WIDTH):
            w(x)
            w('<BR ALIGN="LEFT"/>')
    return "<%s>" % "".join(lines)


def _format_modulenode_gv(node, node_id, options):
    title_size = 32

    node_name = "%s [%d]" % (node.name, node_id)

    if options is None:
        options = {}

    lines = []
    w = lines.append
    w('<FONT POINT-SIZE="%d"><B>%s</B></FONT>' % (title_size, node_name))
    if node.option_defs:
        w('<BR/>')
        w('<BR/>')
        w("<U>Module Attributes</U>:")
        w('<BR ALIGN="LEFT"/>')
    for option in node.option_defs:
        name = option.name
        value = options.get(name)
        if value is None:
            value = option.default

        if value is None:
            x = ' <B><FONT COLOR="RED">MISSING</FONT></B>'
        elif value:
            x = " = %s" % value
        else:
            x = ""
        w("%s%s" % (name, x))
        w('<BR ALIGN="LEFT"/>')
    return "<%s>" % "".join(lines)


def plot_network_gv(
    filename, network, options=None, bold=[], bold_transitions=[],
    highlight_green=[], highlight_orange=[], highlight_purple=[],
    highlight_yellow=[], nodeid2color={}, show_node_ids=None, verbose=False):
    # bold              List of node IDs to bold.
    # highlight[1-2]    List of node IDs to highlight.
    # bold_transitions  List of tuples (node_id1, node_id2)
    # show_node_ids     If not None, then show only the node ids in this list.
    # nodeid2color      node_id -> <color> (e.g. "#FDAF91")
    from genomicode import graphviz

    if show_node_ids is None:
        show_node_ids = {}.fromkeys(range(len(network.nodes)))
    if bold is None:
        bold = []
    if bold_transitions is None:
        bold_transitions = []
    if highlight_green is None:
        highlight_green = []
    if highlight_orange is None:
        highlight_orange = []
    if highlight_purple is None:
        highlight_purple = []
    if highlight_yellow is None:
        highlight_yellow = []

    gv_nodes = []
    gv_edges = []
    gv_node2attr = {}
    gv_edge2attr = {}

    bold_width = 3
    #bold_color = "#960308"
    bold_color = "#CD1A1C"
    #data_color = "#EEEEEE"
    #module_color = "#80E0AA"
    # Brewer brbg-div.  Looks OK.
    #data_color = "#E0C37E"
    #module_color = "#80CDC2"
    # Brewer bugn-seq.  Looks OK.
    #data_color = "#EEF9FC"
    #module_color = "#99D9CA"
    # Brewer blues-seq.  Best.
    data_color = "#F0F4FF"
    module_color = "#BFD4E7"
    # Brewer qualitative.
    #highlight_color = "#FDAF91"    # red
    #highlight1_color = "#E51517"    # red   # set1
    #highlight2_color = "#50B24E"    # green
    #highlight3_color = "#F97808"    # orange
    #highlight1_color = "#63C2A4"    # green   # set2
    #highlight2_color = "#FB875D"    # orange
    #highlight3_color = "#50B24E"    # purple
    highlight_green_color = "#B4E2CE"    # green   # pastel
    highlight_orange_color = "#FACEAE"    # orange
    highlight_purple_color = "#F5CAE4"    # purple
    highlight_yellow_color = "#FEF2AD"    # yellow
    # Brewer rdbu-div.  Too much contrast.
    #data_color = "#F5A582"
    #module_color = "#92C6DF"
    # Brewer piyg-div.  Too pastel
    #data_color = "#F2B7DB"
    #module_color = "#B8E286"

    id2name = {}
    for node_id, node in enumerate(network.nodes):
        if node_id not in show_node_ids:
            continue
        node2attr = {}
        node2attr["style"] = "filled"
        node2attr["penwidth"] = "1"
        if node.__class__.__name__ == "DataNode":
            name = node.datatype.name
            #node2attr["shape"] = "note"
            node2attr["shape"] = "box"
            node2attr["style"] = "rounded,filled"
            node2attr["fillcolor"] = data_color
            if verbose:
                node2attr["label"] = _format_datanode_gv(node, node_id)

        elif node.__class__.__name__ == "ModuleNode":
            name = node.name
            node2attr["shape"] = "box"
            node2attr["fillcolor"] = module_color
            if verbose:
                node2attr["label"] = _format_modulenode_gv(
                    node, node_id, options)
        else:
            raise AssertionError
        #if node_id in highlight_node1:
        #    node2attr["fillcolor"] = highlight_color
        if node_id in bold:
            node2attr["penwidth"] = bold_width
            node2attr["color"] = bold_color
        if node_id in highlight_yellow:
            node2attr["fillcolor"] = highlight_yellow_color
        elif node_id in highlight_green:
            node2attr["fillcolor"] = highlight_green_color
        elif node_id in highlight_orange:
            node2attr["fillcolor"] = highlight_orange_color
        elif node_id in highlight_purple:
            node2attr["fillcolor"] = highlight_purple_color
        elif node_id in nodeid2color:
            node2attr["fillcolor"] = nodeid2color[node_id]

        node_name = "%s [%d]" % (name, node_id)
        id2name[node_id] = node_name
        gv_nodes.append(node_name)
        gv_node2attr[node_name] = node2attr

    for node_id, next_ids in network.transitions.iteritems():
        if node_id not in show_node_ids:
            continue
        for next_id in next_ids:
            if next_id not in show_node_ids:
                continue
            edge2attr = {}
            x1 = id2name[node_id]
            x2 = id2name[next_id]
            #if node_id in bold and next_id in bold:
            #if (node_id, next_id) in bold_transitions or \
            #       (next_id, node_id) in bold_transitions:
            if (node_id, next_id) in bold_transitions:
                edge2attr["penwidth"] = bold_width
                edge2attr["color"] = bold_color
            gv_edges.append((x1, x2))
            if edge2attr:
                gv_edge2attr[(x1, x2)] = edge2attr

    G = graphviz.make_graph(
        gv_nodes, gv_edges, node2attributes=gv_node2attr,
        edge2attributes=gv_edge2attr, prog="dot", directed=True)
    G.draw(filename)
    #G.write("test.gv")


def read_network(file_or_handle):
    import json

    handle = file_or_handle
    if type(handle) is type(""):
        handle = open(file_or_handle, 'r')
    text = handle.read()
    network = json.loads(text, object_hook=_dict_to_object)
    return network


def write_network(file_or_handle, network):
    import json

    handle = file_or_handle
    if type(handle) is type(""):
        handle = open(file_or_handle, 'w')
    json.dump(network, handle, default=_object_to_dict, indent=2)


def debug_print(print_, s):
    from genomicode import parselib

    if not print_:
        return
    parselib.print_split(s)


def _bc_to_modules(moduledb, out_data):
    # Return list of modules that can generate an output that is
    # compatible with data.

    modules = []  # list of (module, num compatible attributes)
    for module in moduledb:
        if _is_valid_output(module, out_data):
            modules.append(module)
    return modules


def _bc_to_inputs(
    network, module_id, out_data_id, custom_attributes,
    force_default_input_attribute_to_be_all_values=False):
    # Return a list of tuples of input objects.
    import itertools
    # If INPUTs have ENUM, can lead to cycles.
    # subtract_mouse_reads -> Fastq.sub=[yes,no] -> align -> CIGAR ->
    #   subtract_mouse_reads
    # If this is True, will generate inputs that do not have ENUM.
    inputs_have_atomic_attribute_values = True
    #inputs_have_atomic_attribute_values = False

    module = network.nodes[module_id]
    out_data = network.nodes[out_data_id]

    all_attributes = []
    all_attrsource = []
    force = force_default_input_attribute_to_be_all_values
    for in_num in range(len(module.in_datatypes)):
        x = _bc_to_one_input(
            network, module_id, in_num, out_data_id, custom_attributes,
            force_default_input_attribute_to_be_all_values=force)
        attributes, attrsource = x
        all_attributes.append(attributes)
        all_attrsource.append(attrsource)

    # Handle the SAME_AS constraints here.  Can't do in
    # _bc_to_one_input because these constraints require the
    # comparison of multiple inputs.
    for constraint in module.constraints:
        if constraint.behavior != SAME_AS:
            continue

        # BUG: Need to handle user attributes with SAME_AS.

        # If this constraint is the SAME_AS another one, then use the
        # value of the copied constraint.
        i_src = constraint.arg1
        i_dst = constraint.input_index
        assert i_src < len(module.in_datatypes)
        assert i_dst < len(module.in_datatypes)
        #input_src = all_inputs[i_src]
        #input_dst = all_inputs[i_dst]

        attr_src = all_attributes[i_src]
        attr_dst = all_attributes[i_dst]

        name = constraint.name
        assert name in attr_src, \
               "Invalid attribute %s (%s:%s)" % (
            name, module.name, module.in_datatypes[i_src].name)
        assert name in attr_dst, \
               "Invalid attribute %s (%s:%s)" % (
            name, module.name, module.in_datatypes[i_dst].name)
        attr_dst[name] = attr_src[name]
        all_attrsource[i_dst][name] = "SAME_AS,%d" % (i_src)

    # Make the input objects.
    all_inputs = []
    for in_num in range(len(module.in_datatypes)):
        in_datatype = module.in_datatypes[in_num]
        out_datatype = out_data.datatype
        attributes = all_attributes[in_num]
        attrsource = all_attrsource[in_num]

        # TODO: Have some better way of diagnosing this error.
        # To diagnose:
        # AssertionError: 'gene_order' is not a known attribute for datatype
        #   ClassLabelFile.
        #print get_node_name(module)

        if inputs_have_atomic_attribute_values:
            attr_names = sorted(attributes)
            attr_values = [attributes[x] for x in attr_names]
            for i in range(len(attr_values)):
                if type(attr_values[i]) is type(""):
                    attr_values[i] = [attr_values[i]]
            inputs_i = []
            for x in itertools.product(*attr_values):
                d = {}
                for key, value in zip(attr_names, x):
                    d[key] = value
                x = in_datatype.output(**d)
                inputs_i.append(x)
            all_inputs.append(inputs_i)
        else:
            x = in_datatype.output(**attributes)
            all_inputs.append([x])

        # Optimization: Don't call debug_print and sorted.
        if not DEBUG_BACKCHAIN:
            continue
        debug_print(
            DEBUG_BACKCHAIN,
            "Backchaining %s (input=%d) -> %s -> %s." %
            (in_datatype.name, in_num, module.name, out_datatype.name))
        for name in sorted(attributes):
            debug_print(
                DEBUG_BACKCHAIN,
                "  %s=%s (%s)" %
                (name, attributes[name], attrsource[name]))

    combos = [tuple(x) for x in itertools.product(*all_inputs)]
    #return all_inputs
    return combos


## # Sources of attribute values.
## SRC_CONSEQUENCE = "consequence"
## SRC_CONSTRAINT = "constraint"
## SRC_CUSTOM = "custom"
## SRC_OUT = "out"
## SRC_DEFAULT = "default"

def _bc_to_one_input(
    network, module_id, in_num, out_data_id, custom_attributes,
    force_default_input_attribute_to_be_all_values=False):
    # Given a module and output_data, return the input_data object
    # that can generate the output.  This should only be called by
    # _bc_to_inputs.  The SAME_AS constraint is handled
    # there.  If this is called by itself, the attributes might not be
    # right.
    #
    # force_default_input_attribute_to_be_all_values can be useful
    # when checking for the possibility that an input node can go into
    # a module.  However, it should be False when generating the
    # network to prevent combinatorial explosion.
    module = network.nodes[module_id]
    out_data = network.nodes[out_data_id]
    assert in_num < len(module.in_datatypes)

    in_datatype = module.in_datatypes[in_num]

    # Can't generate debug messages here because the SAME_AS
    # constraints aren't handled in this function.

    # The attributes for the input object should come from (in
    # decreasing priority):
    # 1.  Consequence (i.e. SAME_AS_CONSTRAINT).
    # 2.  Constraint.
    # 3.  custom attribute    (only if it doesn't conflict with #1 or #2)
    # 4.  out_data            (default_attributes_from)
    # 5.  default output value of input datatype
    #
    # If the user attribute conflicts with a Constraint, the
    # Constraint is higher priority to make sure no objects are
    # generated that the module cannot handle.  However, if the user
    # attribute or default value is a part of the constraint,
    # (e.g. one of many options), then we can refine it with the user
    # attribute (or default).
    #
    # The user attribute also only applies to the first time this node
    # is generated.  E.g.
    # Bam_1.recal=no -> create_targets -> RealignTarget ->
    #   Bam_2.recal=no -> recal -> Bam_3.recal=yes -> call_variants
    # If the user has requested a specific value for recal, then it
    # only applies to Bam_3.
    #
    # Consequence (SAME_AS_CONSTRAINT) is higher priority than
    # constraint because it indicates a more specific value than the
    # constraint.  e.g.:
    #   Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"])
    #   Consequence("quantile_norm", SAME_AS_CONSTRAINT, 0)

    # Start with empty attributes.
    attributes = {}

    # Keep track of the source of the attribute, for debugging.
    attrsource = {}  # attribute name -> source

    # Case 1.  Set the attributes based on the consequences.  If there
    # is a Consequence that is SAME_AS_CONSTRAINT, then the attribute
    # should be determined by the out_data.  e.g.
    # Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"])
    # Consequence("quantile_norm", SAME_AS_CONSTRAINT)
    #
    # The module takes anything, produces the same value.  So the
    # backchainer needs to preserve the value from the out_data.

    # Set values based on SAME_AS_CONSTRAINT.
    # Nothing to do for SET_TO, SET_TO_ONE_OF, BASED_ON_DATA.  (Does
    # not affect input nodes.)
    x = [x for x in module.consequences if x.behavior == SAME_AS_CONSTRAINT]
    # Get the consequences that are based on this datatype.
    x = [x for x in x if x.arg1 == in_num]
    consequences = x
    for consequence in consequences:
        # Copy the value from the output data.
        attributes[consequence.name] = out_data.attributes[consequence.name]
        attrsource[consequence.name] = "consequence"


    # Case 2.  Set the attributes based on the constraints.
    x = [x for x in module.constraints if x.input_index == in_num]
    constraints = x
    for constraint in constraints:
        name = constraint.name

        if constraint.behavior == SAME_AS:
            constraint = _resolve_constraint(constraint, module.constraints)
        assert constraint.behavior in [MUST_BE, CAN_BE_ANY_OF]

        if name not in attributes:
            attributes[name] = constraint.arg1
            attrsource[name] = "constraint"
            continue

        cons_value = constraint.arg1
        attr_value = attributes[name]
        attr_type = _get_attribute_type(attributes[name]) # type of attr_value

        # Since more values may be allowed in the consequence,
        # further refine based on the constraint.  E.g.:
        # Consequence [A, B, C, D]     Case 1
        # Constraint  [A, B]           Case 2

        if constraint.behavior == MUST_BE:
            # cons_value must be TYPE_ATOM
            if attr_type == TYPE_ATOM:
                assert attr_value == cons_value
            elif attr_type == TYPE_ENUM:
                assert cons_value in attr_value
                attributes[name] = cons_value
                attrsource[name] = "consequence+constraint"
            else:
                raise AssertionError
        elif constraint.behavior == CAN_BE_ANY_OF:
            # cons_value must be TYPE_ENUM
            if attr_type == TYPE_ATOM:
                assert attr_value in cons_value  # don't do anything
            elif attr_type == TYPE_ENUM:
                common = _intersect(attr_value, cons_value)
                assert common
                attributes[name] = cons_value
                attrsource[name] = "consequence+constraint"
            else:
                raise AssertionError

    # Case 3.  If the input data object does not proceed to the output
    # data object, then use the attribute provided by the user.
    x = [x for x in module.default_attributes_from if x.input_index == in_num]
    if not x:
        # Look for relevant custom attributes.
        found = False
        for cattrs in custom_attributes:
            # Ignore attributes for other data types.
            if _does_custom_attribute_apply(
                cattrs, network, module_id, in_datatype, attributes):
                found = True
                break
        if found:
            for custom in cattrs.attributes:
                # If the attribute is already set by a constraint,
                # then refine the value by the custom attribute.
                if attrsource.get(custom.name) == "constraint":
                    attr_value = attributes[custom.name]
                    attr_type = _get_attribute_type(attributes[custom.name])
                    cust_type = _get_attribute_type(custom.value)
                    if attr_type == TYPE_ENUM:
                        if cust_type == TYPE_ATOM:
                            if custom.value in attr_value:
                                attributes[custom.name] = custom.value
                                attrsource[custom.name] = "constraint,custom"
                        elif cust_type == TYPE_ENUM:
                            x = _intersect(attr_value, custom.value)
                            if x:
                                attributes[custom.name] = x
                                attrsource[custom.name] = "constraint,custom"
                        else:
                            raise AssertionError

                if custom.name not in attributes:
                    attributes[custom.name] = custom.value
                    attrsource[custom.name] = "default"


    # Case 4.  If default_attributes_from is the same as in_num, then
    # fill with the same values as the out_data.
    indexes = [x.input_index for x in module.default_attributes_from]
    if in_num in indexes:
        for name, value in out_data.attributes.iteritems():
            if name not in attributes:
                attributes[name] = value
                attrsource[name] = "out"

    # Case 5.  Set values from defaults.
    #
    # Using default attributes makes things a lot simpler and
    # mitigates combinatorial explosion.  However, it can close up
    # some possibilities.  What if the value should be something other
    # than the default?
    for attrdef in in_datatype.attribute_defs.itervalues():
        default = attrdef.default_out
        if DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES or \
               force_default_input_attribute_to_be_all_values:
            default = attrdef.values

        # If the attribute is set by a constraint, and the default is
        # a subset of this, should I refine the constraint by the
        # default?
        # o YES: network may suggest that only the default for an
        #   attribute is acceptable, when other values would work.
        # o NO: more correct, but may generate large network.
        if False and attrsource.get(attrdef.name) == "constraint":
            attr_value = attributes[attrdef.name]
            attr_type = _get_attribute_type(attributes[attrdef.name])
            def_type = _get_attribute_type(default)
            if attr_type == TYPE_ENUM:
                if def_type == TYPE_ATOM:
                    if default in attr_value:
                        attributes[attrdef.name] = default
                        attrsource[attrdef.name] = "constraint,default"
                elif def_type == TYPE_ENUM:
                    x = _intersect(attr_value, default)
                    if x:
                        attributes[attrdef.name] = x
                        attrsource[attrdef.name] = "constraint,default"
                else:
                    raise AssertionError

        if attrdef.name not in attributes:
            attributes[attrdef.name] = default
            attrsource[attrdef.name] = "default"


    ## # Start with empty attributes.
    ## attributes = {}

    ## # Keep track of the source of the attribute, for debugging.
    ## attrsource = {}  # attribute name -> source

    ## # Set the attributes in increasing order of priority.  Higher
    ## # priority overwrites lower priority.

    ## # Case 5.  Set values from defaults.
    ## #
    ## # Using default attributes makes things a lot simpler and
    ## # mitigates combinatorial explosion.  However, it can close up
    ## # some possibilities.  What if the value should be something other
    ## # than the default?
    ## for attrdef in in_datatype.attribute_defs.itervalues():
    ##     default = attrdef.default_out
    ##     if DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES or \
    ##            force_default_input_attribute_to_be_all_values:
    ##         default = attrdef.values
    ##     attributes[attrdef.name] = default
    ##     attrsource[attrdef.name] = "default"

    ## # Case 4.  If default_attributes_from is the same as in_num, then
    ## # fill with the same values as the out_data.
    ## indexes = [x.input_index for x in module.default_attributes_from]
    ## if in_num in indexes:
    ##     for name, value in out_data.attributes.iteritems():
    ##         attributes[name] = value
    ##         attrsource[name] = "out_data"

    ## # Do Case 3 last, so can see if it conflicts with Case 1 or Case 2.

    ## # Case 2.  Set the attributes based on the constraints.
    ## for constraint in module.constraints:
    ##     if constraint.input_index != in_num:
    ##         continue

    ##     if constraint.behavior == MUST_BE:
    ##         attributes[constraint.name] = constraint.arg1
    ##         attrsource[constraint.name] = "constraint"
    ##     elif constraint.behavior == CAN_BE_ANY_OF:
    ##         value = constraint.arg1
    ##         source = "constraint"

    ##         # If the user specified an attribute, then refine by it.
    ##         # Refine by default value?
    ##         # o If YES: network may suggest that only the default for
    ##         #   an attribute is acceptable, when other values would
    ##         #   work.
    ##         # o If NO: may generate large network.

    ##         #if attrsource.get(constraint.name) in ["user", "default"]:
    ##         if attrsource.get(constraint.name) in ["user"]:
    ##             x = _get_attribute_type(attributes[constraint.name])
    ##             if x == TYPE_ATOM:
    ##                 if attributes[constraint.name] in value:
    ##                     value = attributes[constraint.name]
    ##                     source = "constraint,%s" % attrsource[constraint.name]
    ##             elif x == TYPE_ENUM:
    ##                 x = _intersect(attributes[constraint.name], value)
    ##                 if x:
    ##                     value = x
    ##                     source = "constraint,%s" % attrsource[constraint.name]
    ##             else:
    ##                 raise AssertionError
    ##         attributes[constraint.name] = value
    ##         attrsource[constraint.name] = source
    ##     elif constraint.behavior == SAME_AS:
    ##         # Handled in _bc_to_inputs.
    ##         pass
    ##     else:
    ##         raise AssertionError

    ## # Case 1.  Set the attributes based on the consequences.  If there
    ## # is a Consequence that is SAME_AS_CONSTRAINT, then the attribute
    ## # should be determined by the out_data.  e.g.
    ## # Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"])
    ## # Consequence("quantile_norm", SAME_AS_CONSTRAINT)
    ## #
    ## # The module takes anything, produces the same value.  So the
    ## # backchainer needs to preserve the value from the out_data.

    ## # Get SAME_AS_CONSTRAINT.  Nothing to do for SET_TO,
    ## # SET_TO_ONE_OF, BASED_ON_DATA.
    ## x = [x for x in module.consequences if x.behavior == SAME_AS_CONSTRAINT]
    ## # Get the consequences that are based on this datatype.
    ## x = [x for x in x if x.arg1 == in_num]
    ## consequences = x
    ## for consequence in consequences:
    ##     n = consequence.name
    ##     source = "consequence"

    ##     # Copy the value from the output data.
    ##     data_value = out_data.attributes[n]
    ##     data_type = _get_attribute_type(data_value)

    ##     # Since more values may be allowed in the consequence,
    ##     # further refine based on the constraint.  E.g.:
    ##     # Constraint  [A, B]
    ##     # Consequence [A, B, C, D]
    ##     x = [x for x in module.constraints if x.name == n]
    ##     x = [x for x in x if x.input_index == in_num]
    ##     assert len(x) > 0
    ##     assert len(x) == 1
    ##     constraint = x[0]
    ##     if constraint.behavior == SAME_AS:
    ##         constraint = _resolve_constraint(constraint, module.constraints)
    ##     if constraint.behavior == MUST_BE:
    ##         # constraint.arg1  <value>
    ##         if data_type == TYPE_ATOM:
    ##             assert constraint.arg1 == data_value
    ##         elif data_type == TYPE_ENUM:
    ##             assert constraint.arg1 in data_value
    ##             data_value = constraint.arg1
    ##             source = "consequence+constraint"
    ##         else:
    ##             raise AssertionError
    ##     elif constraint.behavior == CAN_BE_ANY_OF:
    ##         # constraint.arg1  list of <values>
    ##         if data_type == TYPE_ATOM:
    ##             assert data_value in constraint.arg1
    ##         elif data_type == TYPE_ENUM:
    ##             common = _intersect(constraint.arg1, data_value)
    ##             assert common
    ##             data_value = common
    ##             source = "consequence+constraint"
    ##         else:
    ##             raise AssertionError
    ##     else:
    ##         raise AssertionError

    ##     attributes[n] = data_value
    ##     attrsource[n] = source

    ## # Case 3.  If the input data object does not proceed to the output
    ## # data object, then use the attribute provided by the user.
    ## # Applies:
    ## # 1.  Only the the first (lowest) time this data type is seen in
    ## #     the network.
    ## # 2.  Every time.
    ## attrs1 = {}
    ## attrs2 = {}
    ## x = [x for x in module.default_attributes_from if x.input_index == in_num]
    ## if not x:
    ##     # A constraint may refine a custom attributes.

    ##     # Look for relevant custom attributes.
    ##     for cattrs in custom_attributes:
    ##         # Ignore attributes for other data types.
    ##         _does_custom_attribute_apply(
    ##             cattrs, network, module_id, in_datatype,
    ##             attributes, attrsource)
    ##         if cattrs.datatype.name != in_datatype.name:
    ##             continue
    ##         attrs = attrs1
    ##         if cattrs.all_nodes:
    ##             attrs = attrs2
    ##         for x in cattrs.attributes:
    ##             attrs[x.name] = x.value

    ## # Found potentially relevant custom attributes.  Apply them if
    ## # there are no descendents with the same data type.
    ## if attrs1 and not \
    ##        _has_descendent_of_datatype(network, module_id, in_datatype.name):
    ##     for name, value in attrs1.iteritems():
    ##         attributes[name] = value
    ##         attrsource[name] = "user"
    ## # Set values that should apply to every node in the network.
    ## for name, value in attrs2.iteritems():
    ##     attributes[name] = value
    ##     attrsource[name] = "user"


    return attributes, attrsource


def _does_custom_attribute_apply(
    cattrs, network, module_id, in_datatype, attributes):
    # Should only be called by:
    #   _bc_to_one_input
    if cattrs.datatype.name != in_datatype.name:
        return False
    if not cattrs.all_nodes and _has_descendent_of_datatype(
        network, module_id, in_datatype.name):
        return False
    # If there are conflicts with any of the higher priority
    # attributes, then it doesn't apply.
    for attr in cattrs.attributes:
        value = attributes.get(attr.name)
        if value is None:
            continue
        if not _is_attribute_compatible(attr.value, value):
            return False
    return True


def _has_descendent_of_datatype(network, node_id, datatype_name):
    stack = [node_id]
    seen = {}
    while stack:
        node_id = stack.pop()
        if node_id in seen:
            continue
        seen[node_id] = 1
        node = network.nodes[node_id]
        if isinstance(node, DataNode) and node.datatype.name == datatype_name:
            return True
        stack.extend(network.transitions.get(node_id, []))
    return False


def _bc_to_input_ids(
    network, module_id, custom_attributes,
    all_input_ids=None, all_output_ids=None, nodeid2parents=None):
    # Return a list of tuples of input_ids that match the input
    # datatypes of the module.  Checks the datatypes and makes sure
    # that this combination of data can generate an output in the
    # network.
    # all_input_ids and all_output_ids can be used to restrict the
    # search to specific input or output IDs.
    import itertools

    if not nodeid2parents:
        nodeid2parents = _make_parents_dict(network)
    assert module_id in nodeid2parents
    if all_input_ids is None:
        all_input_ids = nodeid2parents[module_id]
    if all_output_ids is None:
        all_output_ids = network.transitions.get(module_id, [])

    # In about 97% of the cases, the module only has 1 datatype.  In
    # ~90% of the cases, there are 2 input IDs.
    module = network.nodes[module_id]
    assert isinstance(module, ModuleNode)

    # For each in_datatype, find all the data objects that match this
    # type.
    args = []   # list of list node IDs.  Parallel to module.in_datatypes
    for datatype in module.in_datatypes:
        x = [x for x in all_input_ids
             if network.nodes[x].datatype.name == datatype.name]
        args.append(x)

    # If there are many in_datatypes, then this could cause a
    # combinatorial explosion of possibilities.  Optimize by throwing
    # out data objects that fail a constraint.
    if len(args) > 1:
        for i in range(len(args)):
            x = args[i]
            x = [x for x in x if _is_valid_input_i(
                network, module_id, i, network.nodes[x], custom_attributes)]
            args[i] = x

    valid = []
    # This doesn't work.  The same module can have multiple inputs and
    # outputs.
    ## Optimization: Assume existing inputs in the network are valid
    ## and don't check them again.
    ## Fastq.trimmed (no)  -> is_compressed -> Fastq.trimmed (no)
    ## Fastq.trimmed (yes) -> is_compressed -> Fastq.trimmed (yes)
    #if len(module.in_datatypes) == 1:
    #    ids = nodeid2parents.get(module_id, [])
    #    x = [(x, ) for x in ids]
    #    valid.extend(x)

    # The multiple in_datatypes case is harder, because we don't know
    # the order of the inputs.
    for input_ids in itertools.product(*args):
        assert len(input_ids) == len(module.in_datatypes)

        # Don't check again if already done.
        if input_ids in valid:
            continue

        # Actually, duplicated IDs are OK.  Some modules may accept
        # the same data type twice.  Usually, this is a bug in the
        # rules, but we should allow it.
        # Why allow it?

        # No duplicated IDs.
        x = {}.fromkeys(input_ids)
        if len(x) != len(input_ids):
            continue

        # Make sure the inputs are compatible with the module.
        input_datas = [network.nodes[x] for x in input_ids]
        if not _is_valid_inputs(
            network, module_id, input_datas, custom_attributes):
            continue

        # Can't use _fc_to_output_ids -- recursive.
        num_found = 0
        output_datas = _fc_to_outputs(module, input_datas)
        for output_data in output_datas:
            for out_id in all_output_ids:
                if _is_data_compatible(output_data, network.nodes[out_id]):
                    num_found += 1
                    break
        # Make sure at least one of the outputs generated by this
        # module can be found in the network.
        if not num_found:
            continue
        ## Make sure all of the outputs generated by this module are
        ## found in the network.
        #if num_found < len(output_datas):
        #    continue

        # Passes all the tests.
        valid.append(input_ids)

    return valid


def _bc_to_input_and_module_ids(
    network, out_id, custom_attributes, allowed_ids, nodeid2parents, cache):
    paths = []
    module_ids = nodeid2parents.get(out_id, [])
    module_ids = [x for x in module_ids if x in allowed_ids]
    for module_id in module_ids:
        key = module_id, out_id
        if key not in cache:
            x = _bc_to_input_ids(
                network, module_id, custom_attributes, all_output_ids=[out_id],
                nodeid2parents=nodeid2parents)
            cache[key] = x
        combos = cache[key]
        for in_ids in combos:
            assert module_id in allowed_ids
            assert out_id in allowed_ids
            x = [x for x in in_ids if x in allowed_ids]
            if len(x) != len(in_ids):
                continue
            x = in_ids, module_id, out_id
            paths.append(x)
    return paths


def _fc_to_outputs(module, in_datas, ignore_based_on_data=False,
                   ignore_unchanged=False):
    # Generate a list of DataNode objects that can be generated from
    # module and in_datas.  Multiple objects can be generated because
    # the consequences can vary.  E.g. center_genes can set
    # gene_center to either "mean" or "median".
    #
    # In general, the DataNodes will be atomic unless
    # ignore_based_on_data or ignore_unchanged are used.
    #
    # If ignore_based_on_data is a true value, then will leave
    # attributes that are BASED_ON_DATA as TYPE_ENUM (won't be
    # atomic).
    #
    # ignore_unchanged indicates that we can ignore attributes that
    # are TYPE_ENUM, if the parents have the same attribute.  This can
    # be used to distinguish attributes that are merged from different
    # pathways, from those where the ambiguity doesn't matter.
    # E.g.
    #   preprocess_rma  -> preprocess=[rma, mas5]
    #   preprocess_mas5 ->
    # from:
    #   preprocess=[rma, mas5] -> fill_with_zeros -> preprocess=[rma, mas5]
    import itertools

    # Check the input variables.
    assert len(module.in_datatypes) == len(in_datas), module.name
    for i in range(len(module.in_datatypes)):
        assert in_datas[i].datatype.name == module.in_datatypes[i].name

    # Assume that in_datas fulfill the constraints of the module.

    datatype = module.out_datatype

    attributes = {}      # key -> value (ATOM)
    possibilities = {}   # key -> value (ENUM)

    # Priorities (in increasing order):
    # PRI 1.  Default values
    # PRI 2.  Consequences

    # PRI 1.  Set based on default values.
    # Case 1: default_attributes_from is given.
    #         Use the attributes from this default.
    # Case 2: Fill with the default input values of the output
    #         datatype.
    # Case 1.
    if module.default_attributes_from:
        # If there are multiple default_attributes_from, just use the
        # first one.
        # BUG: Is this always the right thing to do?  There shouldn't
        # be multiple ones.
        input_index = module.default_attributes_from[0].input_index
        assert input_index < len(in_datas)
        data = in_datas[input_index]
        assert data.datatype.name == datatype.name
        for key, value in data.attributes.iteritems():
            x = _get_attribute_type(value)
            if x is TYPE_ATOM:
                attributes[key] = value
            elif ignore_unchanged:
                # TYPE_ENUM, but value is not changed.
                attributes[key] = value
            else:
                possibilities[key] = value
        #attributes.update(data.attributes)
    # Case 2.
    else:
        for attrdef in datatype.attribute_defs.itervalues():
            x = _get_attribute_type(attrdef.default_in)
            assert x is TYPE_ATOM
            attributes[attrdef.name] = attrdef.default_in

    ## Make sure these attributes are atomic.
    #for name, value in attributes.iteritems():
    #    x = _get_attribute_type(value)
    #    assert x is TYPE_ATOM, "Not atomic: %s.%s (%s)" % (
    #        datatype.name, name, value)

    # PRI 2.  Set from the consequences.
    #ignore_possibilities = {}  # don't resolve these TYPE_ENUM
    for cons in module.consequences:
        if cons.behavior == SET_TO:
            attributes[cons.name] = cons.arg1
            if cons.name in possibilities:
                del possibilities[cons.name]
        elif cons.behavior == SET_TO_ONE_OF:
            possibilities[cons.name] = cons.arg1
            if cons.name in attributes:
                del attributes[cons.name]
        elif cons.behavior == BASED_ON_DATA:
            if ignore_based_on_data:
                #possibilities[cons.name] = [cons.arg1[0]]
                attributes[cons.name] = cons.arg1
                if cons.name in possibilities:
                    del possibilities[cons.name]
            else:
                possibilities[cons.name] = cons.arg1
                if cons.name in attributes:
                    del attributes[cons.name]
        elif cons.behavior == SAME_AS_CONSTRAINT:
            input_index = cons.arg1
            data = in_datas[input_index]
            value = data.attributes[cons.name]
            x = _get_attribute_type(value)
            if x == TYPE_ATOM or ignore_unchanged:
                attributes[cons.name] = value
                if cons.name in possibilities:
                    del possibilities[cons.name]
            else:
                possibilities[cons.name] = value
                if cons.name in attributes:
                    del attributes[cons.name]
        else:
            raise AssertionError

    # Make sure nothing is in both possibilities and attributes.
    for x in possibilities:
        assert x not in attributes

    # If no possibilities, then make one output variable.
    if not possibilities:
        x = DataNode.__new__(DataNode)
        x.datatype = datatype
        x.attributes = attributes.copy()
        #return _make_data_node_atomic(x)
        return [x]

    # Otherwise, put them into different output variables.
    names = sorted(possibilities)
    args = [possibilities[x] for x in names]
    outputs = []
    for values in itertools.product(*args):
        for key, value in zip(names, values):
            attributes[key] = value
        # Optimization: DataNode.__init__ is very expensive
        # because of all the checks.  Skip the checks and
        # instantiate the class directly.
        #x = DataNode(datatype, **attributes)
        x = DataNode.__new__(DataNode)
        x.datatype = datatype
        x.attributes = attributes.copy()
        outputs.append(x)

    # TODO: Should make a better check here.
    #if not ignore_based_on_data and not ignore_unchanged:
    #    for x in outputs:
    #        assert _is_data_node_atomic(x)

    return outputs


def _fc_to_output_ids(
    network, module_id, custom_attributes,
    all_input_ids=None, all_output_ids=None, nodeid2parents=None):
    # Return a list of (in_data_ids, module_id, out_data_id) that can be
    # generated by this module.
    combos = _bc_to_input_ids(
        network, module_id, custom_attributes, all_input_ids=all_input_ids,
        all_output_ids=all_output_ids, nodeid2parents=nodeid2parents)
    paths = []
    for in_data_ids in combos:
        # For a given combo and module, find the children it
        # generates.
        in_datas = [network.nodes[x] for x in in_data_ids]
        output_data_ids = network.transitions[module_id]
        if all_output_ids is not None:
            output_data_ids = [
                x for x in output_data_ids if x in all_output_ids]
        # output_datas should be atomic.
        output_datas = _fc_to_outputs(network.nodes[module_id], in_datas)
        for out_data_id in output_data_ids:
            out_node = network.nodes[out_data_id]
            found = False
            for output_data in output_datas:
                if _is_data_compatible(output_data, out_node):
                    found = True
                    break
            if found:
                x = in_data_ids, module_id, out_data_id
                paths.append(x)
                break
    return paths


def _resolve_constraint(constraint, all_constraints):
    # If this should be the same as another constraint, then check the
    # other constraint.
    # CONSEQUENCE   NAME
    # CONSTRAINT 1  NAME  CAN_BE_AN_OF
    # CONSTRAINT 2  NAME  SAME_AS  1
    # CONSTRAINT 3  NAME  SAME_AS  2
    #
    # Given CONSTRAINT_2 or CONSTRAINT_3, return CONSTRAINT_1 (the one
    # that has the actual value.
    const = constraint
    assert const.behavior == SAME_AS
    assert const.arg1 != const.input_index
    #assert const.arg1 < module.in_datatypes

    while const.behavior == SAME_AS:
        #x = [x for x in module.constraints if x.name == const.name]
        x = [x for x in all_constraints if x.name == const.name]
        x = [x for x in x if x.input_index == const.arg1]
        assert len(x) > 0, (
            "%r SAME_AS %d, but datatype %d has no constraint on %r." %
            (const.name, const.arg1, const.arg1, const.name))
        assert len(x) == 1
        const = x[0]
    return const


def _is_valid_inputs(network, module_id, in_datas, custom_attributes,
                     out_data_ids=None):
    # If in_datas is compatible with any of the out_datas, then return
    # True.

    # in_datas is a list of DataNodes.  The DataNodes may contain ENUM
    # attributes, which are hard to compare.  Split them up into
    # atomic values.
    # List of lists of DataNodes.
    atomic_in_datas = [_make_data_node_atomic(x) for x in in_datas]

    module = network.nodes[module_id]
    assert len(in_datas) == len(module.in_datatypes)
    if out_data_ids is None:
        out_data_ids = network.transitions.get(module_id, [])
    for out_data_id in out_data_ids:
        # Make list of all the inputs that can lead to this
        # out_data_id.
        combos = _bc_to_inputs(
            network, module_id, out_data_id, custom_attributes)
        # If each of the in_datas is compatible with the corresponding
        # data nodes in the combos, then, this is compatible.
        for i in range(len(atomic_in_datas)):
            nodes_user = atomic_in_datas[i]
            nodes_bc = [x[i] for x in combos]
            if not _is_atomic_data_compatible(nodes_user, nodes_bc):
                break
        else:
            return True
        #all_inputs = _bc_to_inputs(
        #    network, module_id, out_data_id, custom_attributes,
        #    force_default_input_attribute_to_be_all_values=True)
        #for i in range(len(in_datas)):
        #    if not _is_data_compatible(in_datas[i], all_inputs[i]):
        #        break
        #else:
        #    return True
    return False


def _is_valid_input_i(
    network, module_id, input_num, in_data, custom_attributes):
    module = network.nodes[module_id]
    assert input_num < len(module.in_datatypes)
    # If in_datas is compatible with any of the out_datas, then return
    # True.

    # No.  This doesn't happen with new inference engine.
    ## Subtle corner case.
    ## Fastq.trimmed (yes, no) -> merge_reads -> Fastq.trimmed (no)
    ## Fastq.trimmed (yes, no) -> merge_reads -> Fastq.trimmed (yes)
    ## Will test in_data against each of the outputs.  Since in_data is
    ## a superset of each output, will erroneously say that module
    ## cannot take this input.
    ## Solution: merge the outputs when possible.

    # List of data nodes.
    atomic_in_data = _make_data_node_atomic(in_data)

    out_data_ids = network.transitions.get(module_id, [])
    for out_data_id in out_data_ids:
        combos = _bc_to_inputs(
            network, module_id, out_data_id, custom_attributes)
        nodes_bc = [x[input_num] for x in combos]
        if _is_atomic_data_compatible(atomic_in_data, nodes_bc):
            return True
        #all_inputs = _bc_to_inputs(
        #    network, module_id, out_data_id, custom_attributes,
        #    force_default_input_attribute_to_be_all_values=True)
        #if _is_data_compatible(in_data, all_inputs[input_num]):
        #    return True

        #all_inputs = _bc_to_inputs(
        #    network, module_id, out_data_id, custom_attributes,
        #    force_default_input_attribute_to_be_all_values=True)
        #if _is_data_compatible(in_data, all_inputs[input_num]):
        #    return True
    return False


def _is_valid_input_ids(
    network, module_id, in_data_ids, out_id, custom_attributes,
    nodeid2parents):
    # Optimization: Only check for modules that can generate two
    # or more outputs.
    #next_ids = network.transitions.get(module_id, [])
    #if len(next_ids) <= 1:
    #    return True
    assert type(out_id) is type(0)
    if _bc_to_input_ids(
        network, module_id, custom_attributes, all_input_ids=in_data_ids,
        all_output_ids=[out_id], nodeid2parents=nodeid2parents):
        return True
    return False


def _is_valid_output(module, data):
    # Return whether this module can produce this data object.

    # A module cannot produce this data if:
    # - The module's output data type is not the same as the data.
    # - One or more of the consequences conflict.
    # - An input does not go into the output, and a user_attribute
    #   doesn't match a constraint.  OBSOLETE.  Constraint now takes
    #   precedence over user_attribute.  So user_attribute doesn't
    #   matter.
    # - An input does not go into the output, and the data attribute
    #   doesn't match the (in) defaults of the output data type.
    #
    # A module can produce this data if:
    # - An consequence (SET_TO, SET_TO_ONE_OF, BASED_ON_DATA) that is not
    #   a side effect has a value that matches the value of the data.
    # - The module only converts the datatype.  There are no
    #   consequences (SET_TO, SET_TO_ONE_OF, BASED_ON_DATA), and the
    #   output data type has no attributes.
    #   e.g. download_geo_GSEID  gseid -> expression_files  (no attributes)
    #
    # These rules need to match the policies in _bc_to_one_input.

    # If this module doesn't produce the same data type, then it can't
    # produce this data object.
    if module.out_datatype.name != data.datatype.name:
        #debug_print(
        #    "ModuleNode %s can't generate data type: %s." %
        #    (module.name, data.datatype.name))
        return False

    debug_print(DEBUG_BACKCHAIN, "Testing if module %s can produce data %s." %
                (repr(module.name), str(data)))
    # If any of the consequences conflict, then the module can't produce
    # this data object.
    for consequence in module.consequences:
        assert consequence.name in data.attributes
        data_value = data.attributes[consequence.name]
        data_type = _get_attribute_type(data_value)

        assert data_type in [TYPE_ATOM, TYPE_ENUM]

        if consequence.behavior == SET_TO:
            outc_value = consequence.arg1
            outc_type = TYPE_ATOM
        elif consequence.behavior in [SET_TO_ONE_OF, BASED_ON_DATA]:
            outc_value = consequence.arg1
            outc_type = TYPE_ENUM
        elif consequence.behavior == SAME_AS_CONSTRAINT:
            # Get the value from the constraint.
            datatype_index = consequence.arg1
            assert type(datatype_index) is type(0)
            assert datatype_index < len(module.in_datatypes)
            x = [x for x in module.constraints if x.name == consequence.name]
            x = [x for x in x if x.input_index == datatype_index]
            assert len(x) == 1
            constraint = x[0]

            # If this should be the same as another constraint, then
            # check the other constraint.
            if constraint.behavior == SAME_AS:
                assert constraint.arg1 < len(module.in_datatypes)
                constraint = _resolve_constraint(
                    constraint, module.constraints)
                #x = [x for x in module.constraints
                #     if x.name == consequence.name]
                #x = [x for x in x if x.input_index == constraint.arg1]
                #assert len(x) == 1
                #constraint = x[0]
            if constraint.behavior == MUST_BE:
                outc_value = constraint.arg1
                outc_type = TYPE_ATOM
            elif constraint.behavior == CAN_BE_ANY_OF:
                outc_value = constraint.arg1
                outc_type = TYPE_ENUM
            else:
                raise NotImplementedError, constraint.behavior
        else:
            raise AssertionError

        assert data_type in [TYPE_ATOM, TYPE_ENUM]
        assert outc_type in [TYPE_ATOM, TYPE_ENUM]
        case = _assign_case_by_type(data_type, outc_type)

        if case == 1:
            if data_value != outc_value:
                msg = ("Consequence '%s' requires '%s', "
                       "but data contains '%s'." %
                       (consequence.name, outc_value, data_value))
                debug_print(DEBUG_BACKCHAIN, msg)
                return False
        elif case == 2:
            # ModuleNode can produce any of a list of values.  Check
            # if the data's value can be produced by the module.
            if data_value not in outc_value:
                debug_print(
                    DEBUG_BACKCHAIN, "Consequence %s conflicts." %
                    consequence.name)
                return False
        elif case == 3:
            # ModuleNode produces a specific value.  DataNode could be
            # one of many values.
            if outc_value not in data_value:
                debug_print(
                    DEBUG_BACKCHAIN, "Consequence %s conflicts." %
                    consequence.name)
                return False
        elif case == 4:
            if not _intersect(data_value, outc_value):
                debug_print(
                    DEBUG_BACKCHAIN, "Consequence %s conflicts." %
                    consequence.name)
                return False
        else:
            raise AssertionError

    # Make sure the module's constraints are aligned with the
    # user_attributes.  THIS IS OBSOLETE.  Constraint now takes
    # precedence over user_attribute.  So user_attribute doesn't
    # matter.

    # Get a list of the in_datatypes that don't continue into the
    # out_datatype.
    #indexes = [x.input_index for x in module.default_attributes_from]
    #for i in range(len(module.in_datatypes)):
    #    if i in indexes:
    #        # The values from this datatype should be passed through.
    #        # The user attributes does not apply.
    #        continue
    #
    #    user_attrs = [
    #       x for x in user_attributes
    #        if x.datatype == module.in_datatypes[i]]
    #    for attr in user_attrs:
    #        x = [x for x in module.constraints if x.input_index == i]
    #        x = [x for x in x if x.name == attr.name]
    #        constraints = x
    #
    #        for cons in constraints:
    #            if cons.behavior == MUST_BE:
    #                if attr.value != cons.arg1:
    #                    debug_print(
    #                        "Consequence %s conflicts with user attribute." %(
    #                            cons.name))
    #                    return False
    #            elif cons.behavior == CAN_BE_ANY_OF:
    #                if attr.value not in cons.arg1:
    #                    debug_print(
    #                        "Consequence %s conflicts with user attribute." %(
    #                            cons.name))
    #                    return False
    #            elif cons.behavior == SAME_AS:
    #                # No conflict with user_attribute.
    #                pass
    #            else:
    #                raise AssertionError

    # If the module converts the datatype, and no
    # DefaultAttributesFrom is specified, then the data should match
    # the (in) defaults from the output data type.
    if not module.default_attributes_from:
        debug_print(
            DEBUG_BACKCHAIN,
            "ModuleNode converts datatype.  Checking default attributes.")
        consequence_names = [x.name for x in module.consequences]
        for attrdef in module.out_datatype.attribute_defs.itervalues():
            # Ignore the attributes that have consequences.
            if attrdef.name in consequence_names:
                debug_print(
                    DEBUG_BACKCHAIN,
                    "Attr %r: Skipping--has consequence." % attrdef.name)
                continue
            assert attrdef.name in data.attributes
            data_value = data.attributes[attrdef.name]
            data_type = _get_attribute_type(data_value)
            assert data_type in [TYPE_ATOM, TYPE_ENUM]

            if data_type == TYPE_ATOM:
                if attrdef.default_in != data_value:
                    debug_print(
                        DEBUG_BACKCHAIN,
                        "Attr %r: Conflicts (module %r, data %r)." %
                        (attrdef.name, attrdef.default_in, data_value))
                    return False
            elif data_type == TYPE_ENUM:
                if attrdef.default_in not in data_value:
                    debug_print(
                        DEBUG_BACKCHAIN,
                        "Attr %r: Conflicts (module %r, data %r)." %
                        (attrdef.name, attrdef.default_in, data_value))
                    return False
            else:
                raise AssertionError
            debug_print(
                DEBUG_BACKCHAIN, "Attr %r: matches defaults." % attrdef.name)

    # TESTING.
    # If the module converts the datatype, the consequences don't
    # conflict, and the default attributes don't conflict, then this
    # should match.
    if not module.default_attributes_from:
        debug_print(DEBUG_BACKCHAIN, "Match because of converting datatype.")
        return True

    # At this point, the module produces this datatype and there are
    # no conflicts.  Look for an consequence that is not a side effect
    # that changes the value of an output attribute.
    for consequence in module.consequences:
        if consequence.name not in data.attributes:
            continue
        if consequence.side_effect:
            continue

        if consequence.behavior in [SET_TO, SET_TO_ONE_OF, BASED_ON_DATA]:
            debug_print(
                DEBUG_BACKCHAIN, "Consequence '%s' matches." %
                consequence.name)
            return True

        assert consequence.behavior == SAME_AS_CONSTRAINT

        # If the value of the output attribute is the same as the
        # input attribute, then this does not change the data.
        #
        # If:
        # - this consequence refers to a different object, and
        # - there is a constraint that refers to the same object with
        #   a different value,
        # then this is a match.
        #
        # Example:
        # [GeneListFile, SignalFile] -> SignalFile
        # Constraint("gene_order", CAN_BE_ANY_OF, ["pvalue", "fdr"], 0)
        # Constraint("gene_order", MUST_BE, "no", 1)
        # Consequence("gene_order", SAME_AS_CONSTRAINT, 0)

        # Make sure it refers to a different object.
        indexes = [x.input_index for x in module.default_attributes_from]
        if consequence.arg1 in indexes:
            continue

        # Find the constraint that goes with this consequence.
        x = [
            x for x in module.constraints
            if x.name == consequence.name and x.input_index == consequence.arg1
        ]
        assert len(x) == 1
        const1 = x[0]

        # Find the constraint that refers to the same object.
        x = [x for x in module.constraints
             if x.name == consequence.name and x.input_index in indexes]
        if not x:
            continue
        # If there are multiple constraints, make sure they have the
        # same values.
        if len(x) >= 2:
            raise NotImplementedError
        const2 = x[0]

        # Follow SAME_AS.
        if const1.behavior == SAME_AS:
            const1 = _resolve_constraint(const1, module.constraints)
        if const2.behavior == SAME_AS:
            const2 = _resolve_constraint(const2, module.constraints)
        #while const1.behavior == SAME_AS:
        #    x = [x for x in module.constraints
        #         if x.name == consequence.name and
        #         x.input_index == const1.arg1]
        #    assert len(x) == 1
        #    const1 = x[0]
        #while const2.behavior == SAME_AS:
        #    x = [x for x in module.constraints
        #         if x.name == const2.name and x.input_index == const2.arg1]
        #    assert len(x) == 1
        #    const2 = x[0]

        assert const1.behavior in [MUST_BE, CAN_BE_ANY_OF]
        assert const2.behavior in [MUST_BE, CAN_BE_ANY_OF]
        if (const1.behavior, const2.behavior) == (MUST_BE, MUST_BE):
            if const1.arg1 != const2.arg1:
                debug_print(
                    DEBUG_BACKCHAIN, "Consequence '%s' matches." %
                    consequence.name)
                return True
        elif (const1.behavior, const2.behavior) == (MUST_BE, CAN_BE_ANY_OF):
            if const1.arg1 not in const2.arg1:
                debug_print(
                    DEBUG_BACKCHAIN, "Consequence '%s' matches." %
                    consequence.name)
                return True
        elif (const1.behavior, const2.behavior) == (CAN_BE_ANY_OF, MUST_BE):
            if const2.arg1 not in const1.arg1:
                debug_print(
                    DEBUG_BACKCHAIN, "Consequence '%s' matches." %
                    consequence.name)
                return True
        elif (const1.behavior, const2.behavior) == \
                 (CAN_BE_ANY_OF, CAN_BE_ANY_OF):
            if not _intersect(const1.arg1, const2.arg1):
                debug_print(
                    DEBUG_BACKCHAIN, "Consequence '%s' matches." %
                    consequence.name)
                return True
        else:
            raise AssertionError

    # No conflicts, and the module has no consequences.
    if not module.consequences:
        debug_print(
            DEBUG_BACKCHAIN, "Match because there are no consequences.")
        return True

    # No consequences match.
    debug_print(DEBUG_BACKCHAIN, "No consequences match.")
    return False


def _is_valid_outdata_id_path(network, path, out_id, custom_attributes,
                              nodeid2parents):
    # Can the nodes in this pathway produce out_id.  out_id must be a
    # Data node.

    # Find the parents of out_id.
    x = nodeid2parents.get(out_id, [])
    parent_ids = [x for x in x if x in path]
    assert parent_ids

    assert isinstance(network.nodes[out_id], DataNode)

    for module_id in parent_ids:
        x = nodeid2parents.get(module_id, [])
        prev_ids = [x for x in x if x in path]
        if _is_valid_input_ids(
            network, module_id, prev_ids, out_id, custom_attributes,
            nodeid2parents):
            return True
    return False


def _is_valid_outmodule_id_path(network, paths, parent_ids, module_id,
                                custom_attributes, nodeid2parents):
    # Can the nodes in this pathway produce module_id.  out_id must be a
    # Module node.
    # paths should be parallel to parent_ids.  Each should be a path
    # that creates the parent_id.
    assert isinstance(network.nodes[module_id], ModuleNode)

    # If these are all the same, then no need to check again.
    # module_id, parent_ids, grandparent_ids

    # Look for conflicts in the SAME_AS constraints.
    # Bam.mouse=no -> sort -> Bam.mouse=yes,no -> count_with_htseq
    #               ReadStrandedness.mouse=yes ->
    module = network.nodes[module_id]
    x = [
        x for x in module.constraints if x.behavior == SAME_AS]
    same_as_cons = x

    # If no SAME_AS constraints, then there are no conflicts.
    if not same_as_cons:
        return True

    combo_nodes = [
        _get_atomic_data_node_from_pathway(
            network, paths[i], parent_ids[i], custom_attributes,
            ignore_based_on_data=True, ignore_unchanged=True,
            nodeid2parents=nodeid2parents)
        for i in range(len(parent_ids))]

    # Check if these inputs violate SAME_AS constraints.
    for cons in same_as_cons:
        i_src = cons.arg1
        i_dst = cons.input_index
        assert i_src < len(module.in_datatypes)
        assert i_dst < len(module.in_datatypes)

        attrs_src = combo_nodes[i_src].attributes
        attrs_dst = combo_nodes[i_dst].attributes
        a_src = attrs_src[cons.name]
        a_dst = attrs_dst[cons.name]
        if a_src != a_dst:
            return False
    return True


def _is_valid_output_from_input_and_module_ids(
    network, in_data_ids, module_id, out_data_id, user_attrs, cache):
    # See if user_attrs matches any of the outputs from:
    # in_data_ids -> module_id
    in_datas = [network.nodes[x] for x in in_data_ids]
    module = network.nodes[module_id]

    key = in_data_ids, module_id
    if key not in cache:
        x = _fc_to_outputs(module, in_datas)
        cache[key] = x
    out_datas = cache[key]

    num_outputs_matched = 0
    dt_name = network.nodes[out_data_id].datatype.name
    for out_data in out_datas:
        assert out_data.datatype.name == dt_name, "Bad: %s %s" % (
            out_data.datatype.name, dt_name)
        match = True
        for name, uvalue in user_attrs.iteritems():
            dvalue = out_data.attributes[name]
            utype = _get_attribute_type(uvalue)
            dtype = _get_attribute_type(dvalue)
            assert utype == TYPE_ATOM
            assert dtype == TYPE_ATOM
            if uvalue != dvalue:
                match = False
        if match:
            num_outputs_matched += 1
    return (num_outputs_matched > 0)


def _find_same_data(nodes, node):
    # All values need to be exactly equal.  Return index into nodes.
    # -1 if not found.
    assert isinstance(node, DataNode)

    for i, n in enumerate(nodes):
        if not isinstance(n, DataNode):
            continue
        if n.datatype.name != node.datatype.name:
            continue

        attr1 = n.attributes
        attr2 = node.attributes
        is_equal = True
        for key in attr1:
            assert key in attr2
            value1 = attr1[key]
            value2 = attr2[key]
            if not _is_attribute_same(value1, value2):
                is_equal = False
                break
        if is_equal:
            return i
    return -1


def _find_compat_data(network, data_node, ids_to_score=None):
    # Return a list of node_ids that match the user_data data node exactly.
    scores = _score_compat_data(network, data_node, ids_to_score=ids_to_score)
    x = [x for x in scores if x[0] == 0]
    node_ids = [x[1] for x in x]
    return node_ids


def _score_same_data(node1, node2):
    # Return a tuple (score, has_same_datatype, list of attributes with
    # different values).  If the two nodes have a different datatype,
    # then the attributes will be None.  The score is the number of
    # attributes that are different.  It is 99 if the datatypes are
    # different.
    if node1.datatype.name != node2.datatype.name:
        return (99, False, None)

    attrs = []
    for name in node1.attributes:
        V1 = node1.attributes[name]
        V2 = node2.attributes[name]
        if not _is_attribute_same(V1, V2):
            attrs.append(name)
    return (len(attrs), True, attrs)


def _score_compat_data(network, data_node, ids_to_score=None):
    # Return a list of (score, node_id, data_node, list of (name,
    # value in network node, value in user_data).  Sorted by
    # increasing score.  ids_to_score is a list of node IDs in the
    # network to score.  If None, will score them all.

    # Look for the nodes in the network that are compatible with
    # user_data.
    results = []
    for node_id, next_ids in network.iterate(node_class=DataNode):
        if ids_to_score is not None and node_id not in ids_to_score:
            continue
        netw_node = network.nodes[node_id]
        if netw_node.datatype.name != data_node.datatype.name:
            continue

        # Look for incompatible attributes.
        netw_attr = netw_node.attributes
        data_attr = data_node.attributes
        attrs = []
        for key in netw_attr:
            assert key in data_attr
            netw_value = netw_attr[key]
            data_value = data_attr[key]
            if not _is_attribute_compatible(data_value, netw_value):
                attrs.append(key)

        attr_values = []
        for attr in attrs:
            x = attr, netw_attr[attr], data_attr[attr]
            attr_values.append(x)
        x = len(attrs), node_id, data_node, attr_values
        results.append(x)

    return sorted(results)


def _is_data_same(data_1, data_2):
    x = _score_same_data(data_1, data_2)
    score, has_same_datatype, diff_attributes = x
    return score == 0


def _is_data_compatible(data_specific, data_general):
    # Return boolean indicating whether data_specific is compatible
    # with data_general.
    data_s, data_g = data_specific, data_general
    if data_s.datatype.name != data_g.datatype.name:
        return False
    assert len(data_s.attributes) == len(data_g.attributes)
    assert sorted(data_s.attributes) == sorted(data_g.attributes)
    for name in data_s.attributes:
        s_value = data_s.attributes[name]
        g_value = data_g.attributes[name]
        if not _is_attribute_compatible(s_value, g_value):
            return False
    return True


def _is_atomic_data_compatible(datas_specific, datas_general):
    # Return boolean indicating whether datas_specific is compatible
    # with datas_general.  datas_specific and datas_general are lists
    # of DataNodes that are all atomic.

    # Each of the datas_specific must be found in datas_general.
    for d_s in datas_specific:
        for d_g in datas_general:
            if _is_data_compatible(d_s, d_g):
                break
        else:
            return False
    return True


def _is_attribute_same(values1, values2):
    # CASE   N1_TYPE    N2_TYPE    RESULT
    #   1      ATOM       ATOM     OK if ATOM equal.
    #   2      ATOM       ENUM     No.
    #   3      ENUM       ATOM     No.
    #   4      ENUM       ENUM     OK if ENUM equal.
    type1 = _get_attribute_type(values1)
    type2 = _get_attribute_type(values2)
    case = _assign_case_by_type(type1, type2)
    if case == 1:
        if values1 == values2:
            return True
    elif case == 4:
        if sorted(values1) == sorted(values2):
            return True
    return False


def _is_attribute_compatible(value_specific, value_general):
    # CASE  SPECIFIC   GENERAL   RESULT
    #   1     ATOM      ATOM     Check if items are equal.
    #   2     ATOM      ENUM     Check if ATOM in ENUM.
    #   3     ENUM      ATOM     Not compatible.
    #   4     ENUM      ENUM     Check if SPECIFIC is subset of GENERAL

    s_value = value_specific
    g_value = value_general

    # Optimization.  Code below is too slow due to function calls.
    if type(s_value) is type(""):
        if type(g_value) is type(""):
            if s_value != g_value:
                return False
        else:  # assume g_value is a sequence type
            if s_value not in g_value:
                return False
    else:  # assume s_value is a sequence type
        if type(g_value) is type(""):
            # If specific is ENUM and general is ATOM, can't be more
            # compatible.  specific item takes too many possible
            # values.
            return False
        else:
            if not _is_subset(s_value, g_value):
                return False
    #s_type = _get_attribute_type(s_value)
    #g_type = _get_attribute_type(g_value)
    #case = _assign_case_by_type(s_type, g_type)
    #
    #if case == 1:
    #    if s_value != g_value:
    #        return False
    #elif case == 2:
    #    if s_value not in g_value:
    #        return False
    #elif case == 3:
    #    # If specific is ENUM and general is ATOM, can't be more
    #    # compatible.  specific item takes too many possible
    #    # values.
    #    return False
    #elif case == 4:
    #    if not _is_subset(s_value, g_value):
    #        return False
    #else:
    #    raise AssertionError
    return True


def _merge_data_nodes(node1, node2):
    assert isinstance(node1, DataNode)
    assert isinstance(node2, DataNode)
    assert node1.datatype == node2.datatype

    attrs = {}
    for n, v1 in node1.attributes.iteritems():
        v2 = node2.attributes[n]
        attrs[n] = _merge_attribute_values(v1, v2)
    return DataNode(node1.datatype, **attrs)



def _merge_attribute_values(values1, values2):
    t1 = _get_attribute_type(values1)
    t2 = _get_attribute_type(values2)
    if t1 == TYPE_ATOM and t2 == TYPE_ATOM:
        if values1 == values2:
            return values1
        return [values1, values2]
    elif t1 == TYPE_ATOM and t2 == TYPE_ENUM:
        if values1 in values2:
            return values2
        return [values1] + values2
    elif t1 == TYPE_ENUM and t2 == TYPE_ATOM:
        if values2 in values1:
            return values1
        return values1 + [values2]
    else:
        x = [x for x in values1 if x not in values2]
        return x + values2
    raise AssertionError, "How did I get here?"


def _intersect_attributes(attrs1, attrs2):
    # Given two sets of attributes, return the intersection (make
    # everything more strict).  If two attributes conflict, will raise
    # an exception.
    assert sorted(attrs1) == sorted(attrs2)
    attrs = {}
    for key in attrs1:
        v1 = attrs1[key]
        v2 = attrs2[key]
        t1 = _get_attribute_type(v1)
        t2 = _get_attribute_type(v2)

        if t1 == TYPE_ATOM and t2 == TYPE_ATOM:
            assert v1 == v2
            attrs[key] = v1
        elif t1 == TYPE_ATOM and t2 == TYPE_ENUM:
            assert v1 in v2
            attrs[key] = v1
        elif t1 == TYPE_ENUM and t2 == TYPE_ATOM:
            assert v2 in v1
            attrs[key] = v2
        else:
            x = [x for x in v1 if x in v2]
            assert x
            if len(x) == 1:
                x = x[0]
            attrs[key] = x
    return attrs


def _is_data_node_atomic(node):
    # Return whether the attributes of this data node are all
    # TYPE_ATOM.
    assert isinstance(node, DataNode)
    for name, value in node.attributes.iteritems():
        x = _get_attribute_type(value)
        if x != TYPE_ATOM:
            return False
    return True


def _make_data_node_atomic(node):
    # Return a list of DataNodes that are all atomic.
    import itertools

    assert isinstance(node, DataNode)

    # Optimization
    if _is_data_node_atomic(node):
        return [node]

    # Get all the attributes.
    attr_names = []   # list of attribute names
    attr_values = []  # list of list of attribute values
    for name, values in node.attributes.iteritems():
        if _get_attribute_type(values) == TYPE_ATOM:
            values = [values]
        attr_names.append(name)
        attr_values.append(values)

    # Make new DataNodes.
    atomic_nodes = []
    for values in itertools.product(*attr_values):
        attrs = {}
        for name, value in zip(attr_names, values):
            attrs[name] = value
        x = DataNode.__new__(DataNode)
        x.datatype = node.datatype
        x.attributes = attrs
        #x = DataNode(node.datatype, **attrs)
        atomic_nodes.append(x)
    return atomic_nodes


def _is_network_atomic(network):
    # Return whether all data nodes in this network are atomic.
    for node in network.nodes:
        if not isinstance(node, DataNode):
            continue
        if not _is_data_node_atomic(node):
            return False
    return True


def _find_grandparents_in_path(network, path, node_id, nodeid2parents):
    # Find the parents of this node.  Must be a module ID.
    x = nodeid2parents[node_id]
    x = [x for x in x if x in path.node_ids]
    x = [x for x in x if x in path.transitions]
    x = [x for x in x if node_id in path.transitions[x]]
    parent_ids = x
    # This can happen if there are no parents in the network, or if
    # this is a start node.
    assert parent_ids, "No parents 1"
    assert len(parent_ids) == 1, "Multiple parents of data: %d %s %s %s" % (
        node_id, get_node_name(node), parent_ids,
        [get_node_name(network.nodes[x]) for x in parent_ids])
    module_id = parent_ids[0]

    # Find the parents of the module (the grandparents of node_id).
    x = nodeid2parents[module_id]
    x = [x for x in x if x in path.node_ids]
    x = [x for x in x if x in path.transitions]
    x = [x for x in x if module_id in path.transitions[x]]
    parent_ids = x
    assert parent_ids, "No parents 2"

    return module_id, parent_ids


def _get_atomic_data_node_from_pathway(
    network, path, node_id, custom_attributes, ignore_based_on_data=False,
    ignore_unchanged=False, nodeid2parents=None):
    # Return an atomic data node.  Resolves all the TYPE_ENUM
    # attributes based on the pathway.  node_id may not be atomic.
    # However, since this is a pathway, only one atomic DataNode
    # should be able to be generated.  Create it.
    #
    # The data node returned may not be atomic if:
    # 1.  DataNode is the leaf and contains TYPE_ENUM.  e.g.
    #     BamFolder.duplicates_marked=["no", "yes"]
    
    # OPTIMIZATION: Maybe can memoize this function?  Seems to be
    # called repeatedly for the same nodes.  Actually, doesn't take up
    # that much time.
    assert node_id in path.node_ids
    node = network.nodes[node_id]
    assert isinstance(node, DataNode)
    if _is_data_node_atomic(node):
        return node

    if nodeid2parents is None:
        nodeid2parents = _make_parents_dict(network)
    # If there are no parents, or this is the start node, then we
    # cannot resolve any more attributes.
    if node_id not in nodeid2parents:  # no parents
        return node
    if node_id in path.start_ids:
        return node
 
    module_id, parent_ids = _find_grandparents_in_path(
        network, path, node_id, nodeid2parents)
    nid2atomic = {}  # node_id -> Node (atomic)
    for nid in parent_ids:
        x = _get_atomic_data_node_from_pathway(
            network, path, nid, custom_attributes,
            ignore_based_on_data=ignore_based_on_data,
            ignore_unchanged=ignore_unchanged, nodeid2parents=nodeid2parents)
        nid2atomic[nid] = x

    # Find all combinations of the grandparents.
    combos = _bc_to_input_ids(
        network, module_id, custom_attributes, all_input_ids=parent_ids,
        all_output_ids=[node_id], nodeid2parents=nodeid2parents)
    assert combos, "no input ids"

    # For every combination of input IDs, make the potential data types.
    out_datas = []
    for combo in combos:
        module = network.nodes[module_id]
        #in_datas = [network.nodes[x] for x in combo]
        in_datas = [nid2atomic[x] for x in combo]
        x = _fc_to_outputs(
            module, in_datas, ignore_based_on_data=ignore_based_on_data,
            ignore_unchanged=ignore_unchanged)
        out_datas.extend(x)
    assert out_datas, "no out_datas"

    # Not good check.  May not be atomic due to ignore_based_on_data.
    # Make sure all out_datas are atomic.
    #for x in out_datas:
    #    assert _is_data_node_atomic(x)

    # Should only have one atomic out_data.  If there are multiple
    # ones, try to diagnose.
    if len(out_datas) > 1:
        # Figure out which attributes vary.
        attr2values = {}  # name -> list of values
        for out_data in out_datas:
            for name, value in out_data.attributes.iteritems():
                if name not in attr2values:
                    attr2values[name] = []
                if value not in attr2values[name]:
                    attr2values[name].append(value)

        # Make a list of the BASED_ON_DATA consequences.
        # Cannot resolve if it is BASED_ON_DATA.
        based_on_data = []
        for cons in module.consequences:
            if cons.behavior != BASED_ON_DATA:
                continue
            if cons.name in based_on_data:
                continue
            based_on_data.append(cons.name)

        if False:
            # Plot a network for debugging.
            transitions = []
            for n1, n2s in path.transitions.iteritems():
                for n2 in n2s:
                    transitions.append((n1, n2))
            x = [x for x in path.node_ids if x != node_id]
            plot_network_gv(
                "err.png", network, verbose=True, highlight_orange=x,
                highlight_green=[node_id], bold_transitions=transitions)

        err = []
        err.append("Cannot resolve %s [%d]" % (node.datatype.name, node_id))
        err.append("Path: %s" % " ".join(map(str, sorted(path.node_ids))))
        for name in sorted(attr2values):
            values = attr2values[name]
            if len(values) <= 1:
                continue
            if name in based_on_data:
                err.append("%s is BASED_ON_DATA" % cons.name)
            else:
                err.append("Multiple values for %s: %s" % (name, values))
        x = "\n".join(err)
        raise AssertionError, x
    return out_datas[0]


def _merge_start_ids(paths):
    # Return a list of start IDs.  Return None if there are conflicts.
    start_ids = paths[0].start_ids[:]
    for p in paths[1:]:
        assert len(start_ids) == len(p.start_ids)
        for i in range(len(start_ids)):
            if p.start_ids[i] is None:
                pass  # no start ID.  ignore
            elif start_ids[i] == p.start_ids[i]:
                pass  # Same.  ignore.
            elif start_ids[i] is not None:
                return None  # Conflict.  Not None, not the same.
            else:
                start_ids[i] = p.start_ids[i]
    return start_ids


def _does_based_on_data_conflict_with_out_data(network, node_ids, transitions):
    # Make a list of the module_node_ids in any of these
    # branches with BASED_ON_DATA consequences.

    based_on_data = []  # list of (module_node_id, cons name)
    # Precalculating based_on_data for the whole network
    # does not save much time.
    #for module_id in node_ids:
    #    x = [(module_id, x)
    #         for x in moduleid2bod.get(module_id, [])]
    #    if x:
    #        based_on_data_1.extend(x)
    module_node_ids = [x for x in node_ids
         if isinstance(network.nodes[x], ModuleNode)]
    for module_id in module_node_ids:
        x = network.nodes[module_id].consequences
        x = [x for x in x if x.behavior == BASED_ON_DATA]
        x = [(module_id, x.name) for x in x]
        based_on_data.extend(x)

    # Make sure there are no conflicts.
    #conflict = False
    values = {}  # (module_id, name) -> value
    for (module_id, name) in based_on_data:
        key = module_id, name
        for data_id in transitions.get(module_id, []):
            node = network.nodes[data_id]
            assert name in node.attributes
            value = node.attributes[name]
            if key not in values:
                values[key] = value
            elif values[key] != value:
                return True
    return False


import types
def _get_attribute_type(value):
    t = type(value)
    if t is types.StringType:
        return TYPE_ATOM
    elif t is types.ListType:
        return TYPE_ENUM
    elif t is types.TupleType:
        return TYPE_ENUM
    #if type(name) is type(""):
    #    return TYPE_ATOM
    #elif type(name) in [type([]), type(())]:
    #    return TYPE_ENUM
    raise AssertionError, "Unknown attribute type: %s" % str(value)


def _assign_case_by_type(type1, type2):
    types = [(TYPE_ATOM, TYPE_ATOM),
             (TYPE_ATOM, TYPE_ENUM),
             (TYPE_ENUM, TYPE_ATOM),
             (TYPE_ENUM, TYPE_ENUM), ]
    x = (type1, type2)
    assert x in types, "Unknown types: %s %s" % (type1, type2)
    i = types.index(x)
    return i + 1


def _get_pathway_parents_of(network, pathway, node_id, nodeid2parents=None):
    # Return a list of the IDs of the parents that are also in the
    # pathway.
    if nodeid2parents is None:
        nodeid2parents = _make_parents_dict(network)
    assert node_id in pathway.node_ids
    x = nodeid2parents.get(node_id, [])
    x = [x for x in x if x in pathway.node_ids]
    x = [x for x in x if x in pathway.transitions]
    x = [x for x in x if node_id in pathway.transitions[x]]
    return x


def _get_pathway_grandparents_of(
    network, pathway, node_id, nodeid2parents=None):
    # Return list of tuples (parent_id, list of grandparent IDs).
    parent_ids = _get_pathway_parents_of(
        network, pathway, node_id, nodeid2parents=nodeid2parents)
    grandparent_ids = []
    for parent_id in parent_ids:
        x = _get_pathway_parents_of(
            network, pathway, parent_id, nodeid2parents=nodeid2parents)
        x = parent_id, x
        grandparent_ids.append(x)
    return grandparent_ids


def _get_parents_of(network, node_id):
    # Return a list of IDs that point to this node_id.
    assert node_id < len(network.nodes)
    nodeid2parents = _make_parents_dict(network)
    return nodeid2parents.get(node_id, [])
    #ids = {}
    #for nid, next_ids in network.transitions.iteritems():
    #    if node_id in next_ids:
    #        ids[nid] = 1
    #return ids.keys()


def _get_children_of(network, node_id):
    assert node_id < len(network.nodes)
    return network.transitions.get(node_id, [])


def _make_parents_dict_h(network):
    # Return a dictionary of node_id -> prev_node_ids
    nodeid2parents = {}
    for prev_id, node_ids in network.transitions.iteritems():
        for node_id in node_ids:
            if node_id not in nodeid2parents:
                nodeid2parents[node_id] = []
            nodeid2parents[node_id].append(prev_id)
    return nodeid2parents


BACKCHAIN_CACHE = None  # tuple of (network, nodeid2parents)
def _make_parents_dict(network):
    global BACKCHAIN_CACHE
    import copy

    network_cache = nodeid2parents_cache = None
    if BACKCHAIN_CACHE is not None:
        network_cache, nodeid2parents_cache = BACKCHAIN_CACHE
    # Even though this comparison is slow, caching saves a lot of time.
    if network_cache != network:
        network_cache = copy.deepcopy(network)
        nodeid2parents_cache = _make_parents_dict_h(network)
        x1 = network_cache
        x2 = nodeid2parents_cache
        BACKCHAIN_CACHE = x1, x2
    return nodeid2parents_cache


def _make_ancestor_dict_h(network):
    # Return a dictionary of node_id -> all node_ids that are
    # ancestors of node_id.
    from genomicode import parselib
    if not network.nodes:
        return {}

    node2parents = _make_parents_dict(network)

    # Usually see ~1000 iterations.
    MAX_ITER = 1E5

    # Very inefficient algorithm.
    ancestors = {}  # node id -> list of parent node ids.
    all_nodes = list(range(len(network.nodes)))
    niter = 0
    while all_nodes:
        niter += 1
        if not (niter < MAX_ITER):
            errmsg = (
                "Cycle in network.  This can be caused by a problem in the "
                "rules or a bug in the inferencing code.")
            # Try to diagnose the cycle.
            nid = all_nodes[0]
            cycle = [nid]
            transitions = []
            while True:
                nid = cycle[-1]
                x1 = node2parents.get(nid, [])
                x2 = [x for x in x1 if x in all_nodes]
                x3 = [x for x in x2 if x not in cycle]
                if not x3:
                    # Add last node.
                    assert x2
                    cycle.append(x2[0])
                    transitions.append((x2[0], nid))
                    break
                cycle.append(x3[0])
                transitions.append((x3[0], nid))
            plot_network_gv(
                "cycle.png", network, verbose=True, highlight_green=cycle,
                bold_transitions=transitions, show_node_ids=cycle)
            lines = []
            for nid in cycle:
                node = network.nodes[nid]
                if isinstance(node, DataNode):
                    name = node.datatype.name
                else:
                    name = node.name
                x = "%s[%d]" % (name, nid)
                # DEBUG
                #a = "mouse_reads_subtracted"
                #if isinstance(node, DataNode) and a in node.attributes:
                #    x += " %s=%s" % (a, node.attributes[a])
                lines.append(x)
            x = " <- ".join(lines)
            x = parselib.linesplit(x)
            x = "\n".join(x)
            x = "%s\n%s" % (errmsg, x)
            #assert niter < MAX_ITER, errmsg
            assert niter < MAX_ITER, x
        node_id = all_nodes.pop(0)
        parents = node2parents.get(node_id, [])

        # If this node has no parents, then it has no ancestors.
        if not parents:
            ancestors[node_id] = []
            continue

        # If I haven't found the ancestors of all my parents, try
        # again later.
        # BUG: This will result in an infinite loop if there are
        # cycles.  If there is a cycle, there will never be a node
        # with no parents.
        all_found = True
        for parent_id in parents:
            if parent_id not in ancestors:
                all_found = False
                break
        if not all_found:
            all_nodes.append(node_id)
            continue

        # If all the parents are in ancestors already, then the
        # ancestors of this node are the parents and all their
        # ancestors.
        nodes = parents[:]
        for parent_id in parents:
            nodes.extend(ancestors[parent_id])
        nodes = _uniq_intlist(nodes)
        ancestors[node_id] = nodes

    return ancestors


ANCESTOR_CACHE = None  # tuple of (network, ancestor_dict)
def _make_ancestor_dict(network):
    global ANCESTOR_CACHE
    import copy

    network_cache = ancestor_cache = None
    if ANCESTOR_CACHE is not None:
        network_cache, ancestor_cache = ANCESTOR_CACHE
    if network_cache != network:
        ancestor_cache = _make_ancestor_dict_h(network)
        ANCESTOR_CACHE = copy.deepcopy(network), ancestor_cache
    return ancestor_cache


def _make_descendent_dict(network):
    # Return a dictionary of node_id -> all node_ids that are
    # descendents of node_id.
    if not network.nodes:
        return {}

    # Very inefficient algorithm.
    descendents = {}  # node id -> list of descendent node ids.
    stack = list(range(len(network.nodes)))
    niter = 0
    while stack:
        niter += 1
        assert niter < 1E6, "cycle in network"
        node_id = stack.pop(0)
        children = network.transitions.get(node_id, [])

        # If there are no children, then it has no descendents.
        if not children:
            descendents[node_id] = []
            continue

        # If I haven't found the descendents of all my children, try
        # again later.  This will generate an infinite loop if there
        # are cycles.
        all_found = True
        for child_id in children:
            if child_id not in descendents:
                all_found = False
                break
        if not all_found:
            stack.append(node_id)
            continue

        # If all the children are in descendents already, then the
        # descendents of this node are the children and all their
        # descendents.
        nodes = children[:]
        for child_id in children:
            nodes.extend(descendents[child_id])
        nodes = _uniq_intlist(nodes)
        descendents[node_id] = nodes

    return descendents


def _can_reach_by_bc(network, node_id, good_ids):
    # Return a dictionary of all the node IDs that can be reached by
    # backwards chaining.
    reachable_ids = {}
    stack = [node_id]
    while stack:
        nid = stack.pop(0)
        if nid in reachable_ids:
            continue
        if nid not in good_ids:
            continue
        reachable_ids[nid] = 1
        ids = _get_parents_of(network, nid)
        stack.extend(ids)
    return reachable_ids


def _can_reach_by_fc(network, node_id, good_ids=None):
    # Return a dictionary of all the node IDs that can be reached by
    # forward chaining.
    reachable_ids = {}
    stack = [node_id]
    while stack:
        nid = stack.pop(0)
        if nid in reachable_ids:
            continue
        if good_ids is not None and nid not in good_ids:
            continue
        reachable_ids[nid] = 1
        ids = network.transitions.get(nid, [])
        stack.extend(ids)
    return reachable_ids


def _get_custom_names(custom_attributes, datatype_name):
    # Return a list of the attribute names for this datatype.
    names = []
    for cattr in custom_attributes:
        if cattr.datatype.name != datatype_name:
            continue
        for attr in cattr.attributes:
            # Names might be duplicated, if multiple nodes have
            # custom attributes for the same attribute.
            #assert attr.name not in names, "Duplicate: %s" % attr.name
            if attr.name not in names:
                names.append(attr.name)
    return names


def _get_custom_values(custom_attributes, datatype_name, attribute_name):
    # Return a list of values specified for this datatype and
    # attribute name.
    values = []
    for cattr in custom_attributes:
        if cattr.datatype.name != datatype_name:
            continue
        for attr in cattr.attributes:
            if attr.name != attribute_name:
                continue
            if attr.value not in values:
                values.append(attr.value)
    return values


# THIS FUNCTION IS REALLY SLOW.  DO NOT USE.
ITER_UPPER_DIAG_CACHE = {}
def _iter_upper_diag(n):
    global ITER_UPPER_DIAG_CACHE
    if n >= 1024:
        return _iter_upper_diag_h(n)
    if n not in ITER_UPPER_DIAG_CACHE:
        x = list(_iter_upper_diag_h(n))
        ITER_UPPER_DIAG_CACHE[n] = x
    return ITER_UPPER_DIAG_CACHE[n]


def _iter_upper_diag_h(n):
    for i in range(n - 1):
        for j in range(i + 1, n):
            yield i, j


def _intersect(x, y):
    return list(set(x).intersection(y))


def _is_subset(x, y):
    # Return whether x is a subset of y.
    for i in range(len(x)):
        if x[i] not in y:
            return False
    return True


def _flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                #l[i:i+1] = l[i]  # causes problem in pychecker
                l = l[:i] + l[i] + l[i+1:]
        i += 1
    return ltype(l)


def _flatten1_intlist(l):
    # Flatten integer lists, nested at most 1 level.
    return _flatten(l)


def _uniq(seq):
    return {}.fromkeys(seq).keys()


def _uniq_intlist(seq):
    return _uniq(seq)


def _uniq_flatten1_intlist(seq):
    #return _uniq(_flatten(seq))
    return _uniq_intlist(_flatten1_intlist(seq))


def _intlist2bits(int_list):
    bits = 0
    for i in int_list:
        bits = bits | (1<<i)
    return bits


def _print_nothing(s):
    pass


def _print_string(s):
    print s


def _print_line(line, prefix1, prefixn, width, outhandle=None):
    from genomicode import parselib

    outhandle = outhandle or sys.stdout
    lines = parselib.linesplit(
        line, prefix1=prefix1, prefixn=prefixn, width=width)
    for x in lines:
        print >>outhandle, x


def _pretty_attributes(attributes):
    import re

    if not attributes:
        return ""

    all_attributes = sorted(attributes)

    # Separate the proper python variables from other attributes.
    proper = []
    improper = []
    for key in all_attributes:
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


def _fix_node_id_pairs_after_merge(node_id_pairs, merge_id1, merge_id2):
    # n2 was merged into n1.
    # n2 doesn't exist anymore.
    assert merge_id1 != merge_id2
    n1, n2 = merge_id1, merge_id2
    if n1 > n2:   # make sure n1 < n2, for convenience.
        n1, n2 = n2, n1

    # d  < n1 < n2   Don't change d.
    # n1 <  d < n2   Don't change d.
    # n1 < n2 <  d   Subtract d by 1.
    # n1 < n2 =  d   d is now n1.

    # Watch out for weird situations, e.g. 3-way merges.
    # A-B, A-C, B-C
    # After we merge A and B, we are left with A-C, A-C.
    # Then, after we merge the second A-C, we are left with
    # A-A.  Ignore duplicates.

    pairs = []
    for i in range(len(node_id_pairs)):
        d1, d2 = node_id_pairs[i]
        if d1 == n2:
            d1 = n1
        elif d1 > n2:
            d1 -= 1

        if d2 == n2:
            d2 = n1
        elif d2 > n2:
            d2 -= 1

        if d1 != d2:
            pairs.append((d1, d2))
    return pairs


def _fix_node_id_dict_after_merge(node_id_dict, merge_id1, merge_id2):
    # node_id_dict is node_id -> list of node_ids
    # merge_id2 was merged into merge_id1.
    # merge_id2 doesn't exist anymore.
    assert merge_id1 != merge_id2
    n1, n2 = merge_id1, merge_id2
    if n1 > n2:   # make sure n1 < n2, for convenience.
        n1, n2 = n2, n1

    # d  < n1 < n2   Don't change d.
    # n1 <  d < n2   Don't change d.
    # n1 < n2 <  d   Subtract d by 1.
    # n1 < n2 =  d   d is now n1.

    fixed_dict = {}
    for (key_id, value_ids) in node_id_dict.iteritems():
        if key_id == n2:
            key_id = n1
        elif key_id > n2:
            key_id -= 1

        value_ids = value_ids[:]
        for i in range(len(value_ids)):
            if value_ids[i] == n2:
                value_ids[i] = n1
            elif value_ids[i] > n2:
                value_ids[i] -= 1
        if key_id in fixed_dict:
            # happens when merging n2 into n1
            value_ids = value_ids + fixed_dict[key_id]
        # key_id should not be in value_ids
        assert key_id not in value_ids
        # no duplicates in value_ids
        value_ids = {}.fromkeys(value_ids).keys()
        #if key_id in value_ids:
        #    del value_ids[key]
        fixed_dict[key_id] = value_ids
    return fixed_dict


def _fix_ancestor_dict_after_merge(ancestors, merge_id1, merge_id2):
    # 1.  Merge the ancestors of n1 and n2.
    # 2.  Every node that includes n1 or n2 as an ancestor
    #     now inherits all these ancestors.
    ancestors = ancestors.copy()
    x = ancestors[merge_id1] + ancestors[merge_id2]
    ancs = {}.fromkeys(x).keys()
    ancestors[merge_id1] = ancs
    ancestors[merge_id2] = ancs
    for n in ancestors:
        if merge_id1 in ancestors[n] or merge_id2 in ancestors[n]:
            x = ancestors[n] + ancs
            ancestors[n] = {}.fromkeys(x).keys()
    ancestors = _fix_node_id_dict_after_merge(ancestors, merge_id1, merge_id2)
    return ancestors


def _product_network(network, custom_attributes, max_nodes=None,
                     nodeid2id_fn=None):
    # Perform a product operation (find all combinations of inputs to
    # the node) over every node in the network.  max_nodes is the
    # maximum number of nodes that I should perform a product over.
    # Return a list of tuples.
    if max_nodes is None:
        max_nodes = len(network.nodes)

    nodeid2parents = _make_parents_dict(network)
    cache = {}
    x = _product_network_h(
        network, 0, custom_attributes, nodeid2parents, max_nodes,
        nodeid2id_fn, cache)
    return x.keys()


def _product_network_h(
    network, node_id, custom_attributes, nodeid2parents, max_nodes,
    nodeid2id_fn, cache):
    if node_id not in cache:
        x = _product_network_hh(
            network, node_id, custom_attributes, nodeid2parents, max_nodes,
            nodeid2id_fn, cache)
        cache[node_id] = x
    return cache[node_id]


def _product_network_hh(
    network, node_id, custom_attributes, nodeid2parents, max_nodes,
    nodeid2id_fn, cache):
    # Gets called exactly once per node in the network (due to caching
    # in _product_network_h).
    node = network.nodes[node_id]

    inputs = {}
    if isinstance(node, DataNode):
        #if node_id not in skip_ids:
        #    inputs[(node_id, )] = 1
        x = node_id
        if nodeid2id_fn is not None:
            x = nodeid2id_fn(node_id)
        if x is not None:
            inputs[(x, )] = 1
        for previd in nodeid2parents.get(node_id, []):
            x = _product_network_h(
                network, previd, custom_attributes, nodeid2parents, max_nodes,
                nodeid2id_fn, cache)
            inputs.update(x)
    elif isinstance(node, ModuleNode):
        # Find all potential sets of DataNode notes that can feed into me.
        assert node_id in nodeid2parents
        # list of tuples
        combos = _bc_to_input_ids(
            network, node_id, custom_attributes, nodeid2parents=nodeid2parents)

        # Get the inputs from each of the combinations.
        for combo in combos:
            # No.  Combo is the number of input nodes.  However, if
            # multiple of those input nodes can be created by the same
            # upstream node, then the total number of nodes might be
            # fewer.  So it is inappropriate to check for max_nodes
            # here.
            #if max_nodes is not None and len(combo) > max_nodes:
            #    continue
            # Get the inputs for each branch of this combination.
            # list (for each branch) of list of tuples (possible
            # inputs from this branch).
            branch2inputs = [
                _product_network_h(
                    network, x, custom_attributes, nodeid2parents, max_nodes,
                    nodeid2id_fn, cache)
                for x in combo
            ]
            branch2inputs = [x.keys() for x in branch2inputs]

            # Debug.  See what tuples can make each branch.
            #for i in range(len(combo)):
            #    for x in branch2inputs[i]:
            #        print node_id, combo[i], x

            # Optimization: if only one branch, then no need to find
            # combinations.  Just handle it here.  Actually, this
            # optimization doesn't really save much time.
            if len(branch2inputs) == 1:
                for x in branch2inputs[0]:
                    inputs[x] = 1
                continue

            # If any branches are empty (e.g. the branch consists
            # entirely of skip_ids), then skip this set of
            # combinations.
            empty_branches = False
            for x in branch2inputs:
                if not x:
                    empty_branches = True
                    break
            if empty_branches:
                continue

            # Doing the product and chaining takes ~50% of the time in
            # this function.
            #x = itertools.product(*branch2inputs)
            #x = [itertools.chain(*x) for x in x]
            x = _product_and_chain(branch2inputs, max_nodes)
            inputs.update(x)
    else:
        raise AssertionError

    return inputs


def _product_and_chain(lists_of_tuples, max_length):
    # list of list of tuples of integers
    # Find all permutations for the inputs for each branch.
    # [ [(49, 40)]              Each branch has a set of inputs.
    #   [(38,), (35,)] ]
    # ->                        itertools.product
    # [ ((49, 40), (38,)),
    #   ((49, 40), (35,)) ]
    # ->                        flatten or itertools.chain
    # [ (49, 40, 38), (49, 40, 35) ]
    # Is shallow list, so don't need to flatten arbitrarily
    # deeply.
    # Flatten the tuples.
    #x = [_flatten(x) for x in x]
    assert lists_of_tuples

    if max_length is None:
        max_length = 1E6  # shouldn't be network this big

    product = {}
    for x in lists_of_tuples[0]:
        if len(x) <= max_length:
            product[x] = 1
    for i in range(1, len(lists_of_tuples)):
        list_of_tuples = lists_of_tuples[i]
        # Do a product of everything in product with everything in
        # list_of_tuples.
        # Nearly all the time in this function is spent in this
        # function.
        product = _product_and_chain_h(
            product.keys(), list_of_tuples, max_length)
    return product


def _product_and_chain_h(tup_list1, tup_list2, max_length):
    results = {}
    for x1 in tup_list1:
        for x2 in tup_list2:
            x = x1 + x2
            if len(x) > max_length:
                continue
            # No duplicates.
            x = {}.fromkeys(x)
            x = tuple(sorted(x))
            if len(x) > max_length:
                continue
            results[x] = 1
    return results


def _object_to_dict(obj):
    # Convert objects to a dictionary of their representation
    d = {'__class__': obj.__class__.__name__, '__module__': obj.__module__, }
    d.update(obj.__dict__)
    return d


def _dict_to_object(d):
    assert isinstance(d, dict)
    args = dict((key.encode('ascii'), value) for key, value in d.items())
    for key, value in args.iteritems():
        if isinstance(value, dict):
            args[key] = _dict_to_object(value)
        elif isinstance(value, list):
            if value:
                if isinstance(value[0], unicode):
                    args[key] = [i.encode('ascii') for i in value]
        elif isinstance(value, unicode):
            args[key] = value.encode('ascii')
        else:
            assert 'not expected type %s' % value
    inst = args
    if '__class__' in args:
        class_name = args.pop('__class__')
        module_name = args.pop('__module__')
        if '.' in module_name:
            module = __import__(module_name, globals(), locals(),
                                [module_name.split('.')[-1]], -1)
        else:
            module = __import__(module_name)
        class_ = getattr(module, class_name)
        fn = getattr(class_, '_' + class_name + '__init_from_dict')
        inst = fn(args)
    return inst


import sys
try:
    import cbie3
except ImportError:
    pass
else:
    this_module = sys.modules[__name__]
    for name_ in cbie3.__dict__.keys():
        if name_.startswith("__"):
            continue
        this_module.__dict__[name_] = cbie3.__dict__[name_]
