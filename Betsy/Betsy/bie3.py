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
backchain
complete_network
optimize_network
select_start_node
find_node

get_input_datatypes       Show all combinations of datatype that can be inputs.
get_input_nodes           Show all possible input nodes for the network.
group_nodes_by_datatype

summarize_moduledb
check_moduledb

print_modules
print_network
plot_network_gv

read_network
write_network

diagnose_start_node
debug_print

"""
# Functions:
# _backchain_to_modules     Given an output DataNode, find compatible
#                           ModuleNodes.
# _backchain_to_input       DO NOT CALL.
# _backchain_to_all_inputs  Given an output DataNode and ModuleNode, make
#                           all in DataNodes.
# _forwardchain_to_outputs  Given ModuleNode and DataNodes, make output
#                           DataNode.
#
# _backchain_to_ids         Given a Network and target ID, give source IDs.
# _make_backchain_dict
# _make_ancestor_dict
#
# _can_module_take_data                 Need to supply out_node attributes.
# _can_module_take_one_data
# _can_module_in_network_take_data      Uses out_node attributes from Network.
# _can_module_in_network_take_one_data
# _can_module_produce_data
# _get_valid_input_combinations
#
# _can_reach_by_bc
# _can_reach_by_fc
#
# _find_data_node
# _find_start_nodes
# _score_start_nodes
# _score_one_start_node
# _is_data_compatible_with
#
# _assign_case_by_type
# _get_attribute_type
# _intersection
# _is_subset
# _flatten
#
# _print_nothing
# _print_string
# _print_line
# _pretty_attributes
#
# _compare_two_nodes

# _object_to_dict           Use for writing and reading json file
# _dict_to_object           Use for writing and reading json file

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

DEBUG = False
#DEBUG = True

# When backchaining, should we allow the attributes of the input data
# to be all possible values, or fix it to the default?  All possible
# values is correct, but generates a combinatorial explosion that is
# difficult to manage.
DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES = False
#DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES = True

MAX_NETWORK_SIZE = 1024 * 8


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
        #attr = x[0]
        assert datatype.is_valid_attribute_value(name, value), \
               "Invalid value %r for attribute %r." % (value, name)

        self.datatype = datatype
        self.name = name
        self.value = value

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        x = [self.datatype.name, repr(self.name), repr(self.value)]
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


class Constraint:
    def __init__(self, name, behavior, arg1=None, input_index=None):
        if behavior == MUST_BE:
            # name   Name of attribute.
            # arg1   Value of the attribute.
            assert type(arg1) is type("")
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
            raise AssertionError, "Invalid behavior (%s) for constraint %s." % (
                behavior, name)
        assert input_index is None or type(input_index) is type(0)

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


class Consequence:
    def __init__(self, name, behavior,
                 arg1=None,
                 arg2=None,
                 side_effect=False):
        if behavior == SET_TO:
            assert type(arg1) is type("")
            assert arg2 is None
        elif behavior in [SET_TO_ONE_OF, BASED_ON_DATA]:
            assert type(arg1) in [type([]), type(())]
            for x in arg1:
                assert type(x) is type("")
            assert arg2 is None
        elif behavior == SAME_AS_CONSTRAINT:
            if arg1 is None:  # default to datatype 0.
                arg1 = 0
            assert type(arg1) is type(0)  # index of input variable
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


class DefaultAttributesFrom:
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

        ## # Make sure no overlap between the attributes and user inputs.
        ## attr_names = [x.name for x in attributes]
        ## user_names = [x.name for x in user_inputs]
        ## for x in attr_names:
        ##     assert x not in user_names, "%s overlaps in DataType %s." % (
        ##         x, name)

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
    ## def _resolve_attributes(self, attribute_objs, attribute_dict, is_input):
    ##     # Make a dictionary of all the attributes.  The values given
    ##     # by the caller take precedence.  Anything else should be set
    ##     # to the default attributes.
    ##     attrdict = {}
    ##     # Priority 1: Set to the attribute objects.
    ##     for attr in attribute_objs:
    ##         # Ignore attributes for other data types.
    ##         if attr.datatype.name != self.name:
    ##             continue
    ##         attrdict[attr.name] = attr.value
    ##     # Priority 2: Set to the attribute dict.
    ##     for (name, value) in attribute_dict.iteritems():
    ##         if name in attrdict:
    ##             continue
    ##         # Handle the prioritization for the user.
    ##         #assert name not in attrdict, "Conflict: %s" % name
    ##         attrdict[name] = value
    ##     # Priority 3: Set to default attributes.
    ##     for attr in self.attributes:
    ##         if attr.name in attrdict:
    ##             continue
    ##         value = attr.default_in
    ##         if not is_input:
    ##             value = attr.default_out
    ##         attrdict[attr.name] = value
    ##     return attrdict
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

        ## attributes = {}
        ## for name, value in keywds.iteritems():
        ##     if name in attr_names:
        ##         attributes[name] = value
        ##     elif name in user_names:
        ##         user_inputs[name] = value
        ##     else:
        ##         raise AssertionError, "Unknown: %s" % name

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
        inst = DataNode(args['datatype'], **args['attributes'])
        return inst


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

        for x in constraints:
            self._assert_constraint(name, in_datatypes, out_datatype,
                                    constraints, consequences, x)
        for x in consequences:
            self._assert_consequence(name, in_datatypes, out_datatype,
                                     constraints, x)

    def _assert_constraint(self, name, in_datatypes, out_datatype, constraints,
                           consequences, constraint):
        # Get the input datatype that this constraint refers to.
        i = constraint.input_index
        assert i < len(in_datatypes), \
               "Invalid constraint index %d in module %s" % (i, name)
        in_datatype = in_datatypes[i]

        assert constraint.behavior in [MUST_BE, CAN_BE_ANY_OF, SAME_AS]
        if constraint.behavior in [MUST_BE, CAN_BE_ANY_OF]:
            assert in_datatype.is_valid_attribute_value(
                constraint.name,
                constraint.arg1), ("%r: Invalid value %r for attribute %r." %
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
        if in_datatype.name == out_datatype.name:
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
            assert len(x) == 1
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

    def __cmp__(self, other):
        if not isinstance(other, ModuleNode):
            return cmp(id(self), id(other))
        x1 = [self.name, self.in_datatypes, self.out_datatype,
              self.constraints, self.consequences,
              self.default_attributes_from, self.option_defs, self.help]
        x2 = [other.name, other.in_datatypes, other.out_datatype,
              other.constraints, other.consequences,
              other.default_attributes_from, other.option_defs, other.help]
        return cmp(x1, x2)

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

    def points_to(self, node_id):
        # Return a list of all the ids that point to this node.
        prev_ids = []
        for id_ in self.transitions:
            if node_id in self.transitions[id_]:
                prev_ids.append(id_)
        return prev_ids

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

    def merge_nodes(self, node_ids):
        """node_ids is a list of the indexes of nodes.  Replace all
        these nodes with just a single one.  Returns a new Network
        object."""
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
            x = _backchain_to_ids(self, nid)
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


def make_network(moduledb, out_data, *user_attributes):
    x = backchain(moduledb, out_data, user_attributes)
    x = complete_network(x, user_attributes)
    x = optimize_network(x, user_attributes)
    return x


def backchain(moduledb, out_data, user_attributes):
    # Return a Network object.
    assert type(user_attributes) in [type([]), type(())]
    for x in user_attributes:
        assert isinstance(x, Attribute)

    check_moduledb(moduledb)
    if isinstance(out_data, DataType):
        attrdict = {}
        for attr in user_attributes:
            if attr.datatype.name != out_data.name:
                continue
            attrdict[attr.name] = attr.value
        out_data = out_data.output(**attrdict)
    assert isinstance(out_data, DataNode)

    nodes = []  # list of DataNode or ModuleNode objects.
    transitions = {}  # list of index -> list of indexes

    nodes.append(out_data)
    stack = [0]
    seen = {}
    while stack:
        #if len(nodes) >= 100:
        #    break
        assert len(nodes) < MAX_NETWORK_SIZE, "network too large"
        #_print_network(Network(nodes, transitions))

        # Pop the next node off the stack.
        node_id = stack.pop()
        assert node_id < len(nodes)
        node = nodes[node_id]

        # If I've already seen this node, then don't process it again.
        if node_id in seen:
            continue
        seen[node_id] = 1

        if isinstance(node, DataNode):
            # Backwards chain to the previous module.
            modules = _backchain_to_modules(moduledb, node)
            for m in modules:
                nodes.append(m)
                m_id = len(nodes) - 1
                stack.append(m_id)
                transitions[m_id] = transitions.get(m_id, [])
                transitions[m_id].append(node_id)
        elif isinstance(node, ModuleNode):
            cons_id = transitions[node_id][0]
            cons = nodes[cons_id]
            all_inputs = _backchain_to_all_inputs(node, cons, user_attributes)
            for d in all_inputs:
                d_id = _find_data_node(nodes, d)
                if d_id == -1:
                    nodes.append(d)
                    d_id = len(nodes) - 1
                stack.append(d_id)
                transitions[d_id] = transitions.get(d_id, [])
                transitions[d_id].append(node_id)
        else:
            raise AssertionError, "Unknown node type: %s" % node

    # Remove the duplicates from transitions.
    for nid in transitions:
        #next_ids = sorted({}.fromkeys(transitions[nid]))
        next_ids = {}.fromkeys(transitions[nid]).keys()
        transitions[nid] = next_ids

    network = Network(nodes, transitions)
    return network


def _make_ancestor_dict(network):
    # Return a dictionary of node_id -> all node_ids that are
    # ancestors of node_id.
    if not network.nodes:
        return {}

    node2parents = _make_backchain_dict(network)

    ancestors = {}  # node id -> list of parent node ids.
    all_nodes = list(range(len(network.nodes)))
    niter = 0
    while all_nodes:
        niter += 1
        assert niter < 1E6, "cycle in network"
        node_id = all_nodes.pop(0)
        parents = node2parents.get(node_id, [])

        # If there are no parents, then it has no ancestors.
        if not parents:
            ancestors[node_id] = []
            continue

        # If I haven't found the ancestors of all my parents, try
        # again later.
        # BUG: This will generate an infinite loop if there are
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
        ancestors[node_id] = nodes

    return ancestors


def complete_network(network, user_attributes):
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

    debug_print("Completing network.")

    network = copy.deepcopy(network)
    node2previds = _make_backchain_dict(network)
    ancestors = _make_ancestor_dict(network)

    # For each DataNode object, check to see if it can be the
    # antecedent of any ModuleNode objects.
    data_ids = [x for x in range(len(network.nodes))
                if isinstance(network.nodes[x], DataNode)]
    module_ids = [x for x in range(len(network.nodes))
                  if isinstance(network.nodes[x], ModuleNode)]
    for x in itertools.product(data_ids, module_ids):
        input_id, module_id = x
        ##input_node = network.nodes[input_id]
        ##module_node = network.nodes[module_id]

        # If data_id already points to module_id, then ignore
        # this.
        if module_id in network.transitions.get(input_id, []):
            continue

        # Don't add a link from data_id to module_id if it would
        # create a cycle.
        if module_id in ancestors[input_id]:
            #debug_print("Skipping DataNode %d -> ModuleNode %d (cycle)." % (
            #    input_id, module_id))
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

        # Since modules can take multiple inputs, we need to combine
        # input_id with all previous input IDs and try all possible
        # combinations.
        #x = _backchain_to_ids(network, module_id)
        combined_ids = node2previds.get(module_id, []) + [input_id]

        # Find combinations of inputs that are compatible with the
        # network.
        combos = _get_valid_input_combinations(
            network, module_id, combined_ids, user_attributes,
            node2previds=node2previds)

        # Add the new transitions.
        for id_ in itertools.chain.from_iterable(combos):
            assert id_ in network.transitions
            if module_id in network.transitions[id_]:
                continue
            network.transitions[id_].append(module_id)
            debug_print("Completing DataNode %d -> ModuleNode %d." % (
                id_, module_id))
    return network


def optimize_network(network, user_attributes):
    # Clean up the network by removing cycles, extra nodes, etc.
    # Returns a new Network object.
    optimizers = [
        _OptimizeNoCycles(),
        _OptimizeNoDanglingNodes(),
        _OptimizeNoDuplicateModules(),
        _OptimizeNoDuplicateData(),
        _OptimizeNoOverlappingData(),
        _OptimizeNoInvalidOutputs(),    # Fixes NoOverlappingData.
        _OptimizeMergeNodes(),
        ]

    it = 0
    old_network = None
    while old_network != network:
        old_network = network
        for opt in optimizers:
            network = opt.optimize(network, user_attributes)
        it += 1

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

    def optimize(self, network, user_attributes):
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

    def _list_noncycle_node_ids(self, network, nodeid2previds):
        # Return a list of the node_ids that are not in cycles.

        # The nodes at the top of the tree (no prev nodes) are not in
        # cycles.  The nodes at the bottom of the tree (no next nodes)
        # are not in cycles.
        #
        # If all of a nodes next nodes are noncycle, then it is
        # noncycle.
        # If all of a nodes prev nodes are noncycle, then it is
        # noncycle.
        noncycle = {}

        while True:
            changed = False
            node_ids = [i for i in range(len(network.nodes))
                        if i not in noncycle]
            for node_id in node_ids:
                prev_ids = nodeid2previds.get(node_id, [])
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
        nodeid2previds = _make_backchain_dict(network)
        noncycle = self._list_noncycle_node_ids(network, nodeid2previds)

        cycle = None
        node_ids = [i for i in range(len(network.nodes)) if i not in noncycle]
        for start_id in node_ids:
            cycle = self._find_cycle_from_one_node(
                network, start_id, max_path_length, noncycle, nodeid2previds)
            if cycle:
                break
        return cycle

    def _find_cycle_from_one_node(self, network, start_id, max_path_length,
                                  noncycle, nodeid2previds):
        # Do a depth-first search and look for cycles.  Return a cycle
        # (list of node_ids) or None.  The cycle will start and end
        # with the same node_id.
        # Previously did a breadth-first search, but stack.pop(0) was
        # running too slowly (see below).  Depth-first search will be
        # faster to find a cycle, if it exists, anyway.
        if not nodeid2previds:
            nodeid2previds = _make_backchain_dict(network)

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
            for prev_id in nodeid2previds.get(node_id, []):
                stack.append((prev_id, path))
        return cycle

    def _find_depth_of_nodes(self, network):
        # Do a breadth-first search to assign the depth of each node.
        assert network.nodes

        nodeid2previds = _make_backchain_dict(network)

        # OPTIMIZE: can memoize this function.
        stack = [(0, -1)]  # node_id, next_depth
        nodeid2depth = {}
        while stack:
            node_id, next_depth = stack.pop(0)
            if node_id in nodeid2depth:
                continue
            depth = next_depth + 1
            nodeid2depth[node_id] = depth
            for prev_id in nodeid2previds.get(node_id, []):
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


class _OptimizeNoDuplicateData:
    def __init__(self):
        pass

    def optimize(self, network, user_attributes):
        # This could be made much more efficient with a better way of
        # finding duplicates.
        while True:
            duplicates = self.find_duplicate_data(network)
            if not duplicates:
                break
            # Merge each of the duplicates.
            while duplicates:
                n1, n2 = duplicates.pop()
                network = network.merge_nodes([n1, n2])
                duplicates = _fix_node_ids_after_merge(n1, n2, duplicates)
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
        for i, node_id1 in enumerate(data_node_ids):
            node_1 = network.nodes[node_id1]
            for node_id2 in data_node_ids[i+1:]:
                node_2 = network.nodes[node_id2]
                if node_1.datatype.name != node_2.datatype.name:
                    continue
                if node_1 == node_2:
                    duplicates.append((node_id1, node_id2))
        return duplicates
        
        ## data_pairs = itertools.product(data_nodes, data_nodes)

        ## # Can optimize this by sorting.
        ## for (node_id1, node_id2) in data_pairs:
        ##     if node_id1 == node_id2:
        ##         continue
        ##     data1, data2 = network.nodes[node_id1], network.nodes[node_id2]
        ##     # Everything has to be the same.
        ##     if data1 == data2:
        ##         return [node_id1, node_id2]
        ## return []


class _OptimizeNoDuplicateModules:
    def __init__(self):
        pass

    def optimize(self, network, user_attributes):
        while True:
            duplicates = self.find_duplicate_modules(network)
            if not duplicates:
                break
            # Merge each of the duplicates.
            while duplicates:
                n1, n2 = duplicates.pop()
                network = network.merge_nodes([n1, n2])
                duplicates = _fix_node_ids_after_merge(n1, n2, duplicates)
        return network

    def find_duplicate_modules(self, network):
        # Return a list of (node_id1, node_id2) for modules that are
        # duplicated.  If no duplicates found, return an empty list.

        # DEFINITION: If the same data node points to two of the same
        # module nodes, then those modules are duplicated.
        dups = {}
        for node_id, next_ids in network.iterate(node_class=DataNode):
            if len(next_ids) < 2:
                continue
            for i in range(len(next_ids) - 1):
                node_id1 = next_ids[i]
                node_1 = network.nodes[node_id1]
                for j in range(i+1, len(next_ids)):
                    node_id2 = next_ids[j]
                    node_2 = network.nodes[node_id2]
                    if node_1.name != node_2.name:
                        continue
                    id1, id2 = node_id1, node_id2
                    if id1 > id2:
                        id1, id2 = id2, id1
                    assert id1 != id2
                    dups[(id1, id2)] = 1
        dups = dups.keys()
        return dups


def _fix_node_ids_after_merge(merge_id1, merge_id2, node_id_pairs):
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
        elif d1 >= n2:
            d1 -= 1
            
        if d2 == n2:
            d2 = n1
        elif d2 >= n2:
            d2 -= 1
            
        if d1 != d2:
            pairs.append((d1, d2))
    return pairs


class _OptimizeNoDanglingNodes:
    def __init__(self):
        pass

    def optimize(self, network, user_attributes):
        # Remove nodes that have been made irrelevant due to
        # optimizing.
        while True:
            dangling = self.find_dangling_nodes(network, user_attributes)
            if not dangling:
                break
            # Make sure root not is never deleted.
            assert 0 not in dangling
            network = network.delete_nodes(dangling)
        return network

    def find_dangling_nodes(self, network, user_attributes):
        # Return a list of node_ids.
        dangling = []
        for node_id in range(len(network.nodes)):
            if self.is_dangling_node(network, node_id, user_attributes):
                dangling.append(node_id)
        return dangling

    def is_dangling_node(self, network, node_id, user_attributes):
        # 1.  ModuleNodes with no valid input_combinations.
        # 2.  ModuleNodes with no outputs.
        # 3.  DataNode nodes (except for node 0) that don't point to
        #     any ModuleNodes.
        
        # Case 1.
        node = network.nodes[node_id]
        if isinstance(node, ModuleNode):
            prev_ids = _backchain_to_ids(network, node_id)
            x = _get_valid_input_combinations(
                network, node_id, prev_ids, user_attributes)
            if not x:
                return True
        # Case 2.
        if isinstance(node, ModuleNode):
            if not network.transitions.get(node_id, []):
                return True
        # Case 3.
        if isinstance(node, DataNode):
            if node_id != 0 and not network.transitions.get(node_id, []):
                return True

        return False


class _OptimizeNoOverlappingData:
    # Some data objects may have overlapping attributes.
    # 1.  DataNode(SignalFile, format=['tdf', 'pcl', 'gct'], preprocess='rma')
    # 2.  DataNode(SignalFile, format=['tdf', 'res', 'pcl', 'jeffs', 'unknown',
    #       'xls'], preprocess='rma')
    #
    # ModuleNodes that point to 1 can generate 3 different formats.
    # Those that point to 2 can generate 6.  To resolve this, split
    # the DataNode objects and rewire the ModuleNodes.
    # 3.  DataNode(SignalFile, format=['tdf', 'pcl'], preprocess='rma')
    # 4.  DataNode(SignalFile, format='gct', preprocess='rma')
    # 5.  DataNode(SignalFile, format=['res', 'jeffs', 'unknown', 'xls'],
    #       preprocess='rma')

    # Methods:
    # optimize
    #
    # _find_overlapping_data
    # _find_overlapping_attribute
    #
    # _fix_overlapping_data
    # _remove_atom_from_list    Happens very often.
    # _remove_list_from_list    Rare
    # _split_list               Rare

    def __init__(self):
        pass

    def optimize(self, network, user_attributes):
        # Look for pairs of DataNode objects (DataNode_1, DataNode_2)
        # where an attribute is overlapping.
        import copy

        # Make a copy of the network now.  All changes to the network
        # will be done in place.
        network = copy.deepcopy(network)
        while True:
            overlaps = self._find_overlapping_data(network)
            if not overlaps:
                break
            for x in overlaps:
                node_id1, node_id2, attr_name = x
                # The changes may add new nodes to the network.
                # However, they'll append to the end, so the node IDs
                # should still be valid.
                self._fix_overlapping_data(
                    network, node_id1, node_id2, attr_name, user_attributes)

        return network

    def _find_overlapping_data(self, network):
        # Return list of (node_id1, node_id2, name of overlapping
        # attribute).  Empty list means no overlapping data found.
        data_node_ids = [
            node_id for (node_id, node) in enumerate(network.nodes)
            if isinstance(node, DataNode)]

        dt2nodeids = {}  # name of datatype -> list of node IDs
        for node_id in data_node_ids:
            name = network.nodes[node_id].datatype.name
            if name not in dt2nodeids:
                dt2nodeids[name] = []
            dt2nodeids[name].append(node_id)
            
        overlapping = []
        for dt, node_ids in dt2nodeids.iteritems():
            for i, node_id1 in enumerate(node_ids):
                for node_id2 in node_ids[i+1:]:
                    node_1 = network.nodes[node_id1]
                    node_2 = network.nodes[node_id2]
                    attr = self._find_overlapping_attribute(
                        network, node_id1, node_id2)
                    if not attr:
                        continue
                    x = node_id1, node_id2, attr
                    overlapping.append(x)
        return overlapping

    def _find_overlapping_attribute(self, network, data_id1, data_id2):
        # Return the name of the single overlapping attribute or None.
        # data1 and data2 should be exactly the same, except for one
        # attribute with overlapping values.
        data1 = network.nodes[data_id1]
        data2 = network.nodes[data_id2]
        # Not needed.  Check done in _find_overlapping_data.
        #assert isinstance(data1, DataNode)
        #assert isinstance(data2, DataNode)
        if data1.datatype.name != data2.datatype.name:
            return False

        # CASE    DATA1      DATA2     RESULT
        #   1      ATOM       ATOM     OK if ATOMs are equal.
        #   2      ATOM       ENUM     OVERLAP if ATOM in ENUM; DATA2 not root.
        #   3      ENUM       ATOM     OVERLAP if ATOM in ENUM; DATA1 not root.
        #   4      ENUM       ENUM     OVERLAP if ENUMs share ATOMs; not root.
        data1_attr = data1.attributes
        data2_attr = data2.attributes
        assert sorted(data1_attr) == sorted(data2_attr)

        mismatch = False
        overlapping = []  # list of attribute names

        all_attributes = data1_attr
        for key in all_attributes:
            #assert key in data1_attr
            #assert key in data2_attr
            DATA1_VALUE = data1_attr[key]
            DATA2_VALUE = data2_attr[key]
            DATA1_TYPE = _get_attribute_type(DATA1_VALUE)
            DATA2_TYPE = _get_attribute_type(DATA2_VALUE)
            case = _assign_case_by_type(DATA1_TYPE, DATA2_TYPE)

            if case == 1:
                if DATA1_VALUE != DATA2_VALUE:
                    mismatch = True
            elif case == 2:
                if data_id2 == 0:
                    mismatch = True
                elif len(DATA2_VALUE) == 1 and DATA1_VALUE == DATA2_VALUE[0]:
                    pass
                elif DATA1_VALUE in DATA2_VALUE:
                    overlapping.append(key)
                else:
                    mismatch = True
            elif case == 3:
                if data_id1 == 0:
                    mismatch = True
                elif len(DATA1_VALUE) == 1 and DATA2_VALUE == DATA1_VALUE[0]:
                    pass
                elif DATA2_VALUE in DATA1_VALUE:
                    overlapping.append(key)
                else:
                    mismatch = True
            elif case == 4:
                if data_id1 == 0 or data_id2 == 0:
                    mismatch = True
                elif sorted(DATA1_VALUE) == sorted(DATA2_VALUE):  # OK
                    pass
                elif _intersection(DATA1_VALUE, DATA2_VALUE):
                    overlapping.append(key)
                else:
                    mismatch = True
            else:
                raise AssertionError
            # Optimization.  Short circuit the comparison.  Saves a
            # little bit of time, but not a lot.
            if mismatch or len(overlapping) > 1:
                break

        if mismatch:
            return None
        if len(overlapping) != 1:
            return None
        return overlapping[0]

    def _fix_overlapping_data(self, network, data_id1, data_id2, attr_name,
                              user_attributes):
        # Changes network in place.

        data1 = network.nodes[data_id1]
        data2 = network.nodes[data_id2]
        DATA1_VALUE = data1.attributes[attr_name]
        DATA2_VALUE = data2.attributes[attr_name]
        DATA1_TYPE = _get_attribute_type(DATA1_VALUE)
        DATA2_TYPE = _get_attribute_type(DATA2_VALUE)
        if DATA1_TYPE is TYPE_ENUM and DATA2_TYPE is TYPE_ENUM:
            COMMON_VALUES = _intersection(DATA1_VALUE, DATA2_VALUE)
            s_COMMON_VALUES = sorted(COMMON_VALUES)

        if DATA1_TYPE is TYPE_ATOM:
            assert DATA2_TYPE is TYPE_ENUM
            self._remove_atom_from_list(network, data_id2, data_id1, attr_name)
        elif DATA2_TYPE is TYPE_ATOM:
            assert DATA1_TYPE is TYPE_ENUM
            self._remove_atom_from_list(network, data_id1, data_id2, attr_name)
        elif sorted(DATA1_VALUE) == s_COMMON_VALUES:
            self._remove_list_from_list(network, data_id2, data_id1, attr_name)
        elif sorted(DATA2_VALUE) == s_COMMON_VALUES:
            self._remove_list_from_list(network, data_id1, data_id2, attr_name)
        else:
            self._split_list(
                network, data_id1, data_id2, attr_name, user_attributes)

    def _remove_atom_from_list(self, network, data_id1, data_id2, attr_name):
        # Changes network in place.  data1 is a ENUM and data2 is an
        # ATOM.  Remove the ATOM from data1.

        data1 = network.nodes[data_id1]
        data2 = network.nodes[data_id2]
        DATA1_VALUE = data1.attributes[attr_name]
        DATA2_VALUE = data2.attributes[attr_name]
        DATA1_TYPE = _get_attribute_type(DATA1_VALUE)
        DATA2_TYPE = _get_attribute_type(DATA2_VALUE)
        assert DATA1_TYPE == TYPE_ENUM and DATA2_TYPE == TYPE_ATOM
        assert DATA2_VALUE in DATA1_VALUE

        # Remove the ATOM from the ENUM.
        DATA1_VALUE = DATA1_VALUE[:]
        i = DATA1_VALUE.index(DATA2_VALUE)
        DATA1_VALUE.pop(i)
        if len(DATA1_VALUE) == 1:
            DATA1_VALUE = DATA1_VALUE[0]
        data1.attributes[attr_name] = DATA1_VALUE

        # Every module that pointed to data_id1 should now also point
        # to data_id2.
        module_ids = _backchain_to_ids(network, data_id1)
        for node_id in module_ids:
            if data_id2 not in network.transitions[node_id]:
                network.transitions[node_id].append(data_id2)

        # Since some of the workflow from data1 is being rerouted to
        # data2, data2 should point to the children of data1.
        for node_id in network.transitions[data_id1]:
            if node_id not in network.transitions.get(data_id2, []):
                if data_id2 not in network.transitions:
                    network.transitions[data_id2] = []
                network.transitions[data_id2].append(node_id)

    def _remove_list_from_list(self, network, data_id1, data_id2, attr_name):
        # Changes network in place.  data1 and data2 are both a ENUMs.
        # The values of data2 is a subset of data1.  Remove all values
        # of data2 from data1.

        data1 = network.nodes[data_id1]
        data2 = network.nodes[data_id2]
        DATA1_VALUE = data1.attributes[attr_name]
        DATA2_VALUE = data2.attributes[attr_name]
        DATA1_TYPE = _get_attribute_type(DATA1_VALUE)
        DATA2_TYPE = _get_attribute_type(DATA2_VALUE)
        assert DATA1_TYPE == TYPE_ENUM and DATA2_TYPE == TYPE_ENUM
        assert sorted(DATA1_VALUE) != sorted(DATA2_VALUE)
        assert _is_subset(DATA2_VALUE, DATA1_VALUE)

        # Remove the ATOM from the ENUM.
        DATA1_VALUE = DATA1_VALUE[:]
        for value in DATA2_VALUE:
            i = DATA1_VALUE.index(value)
            DATA1_VALUE.pop(i)
        if len(DATA1_VALUE) == 1:
            DATA1_VALUE = DATA1_VALUE[0]
        data1.attributes[attr_name] = DATA1_VALUE

        # Every module that pointed to data_id1 should now also point
        # to data_id2.
        module_ids = _backchain_to_ids(network, data_id1)
        for node_id in module_ids:
            if data_id2 not in network.transitions[node_id]:
                network.transitions[node_id].append(data_id2)

        # Since some of the workflow from data1 is being rerouted to
        # data2, data2 should point to the children of data1.
        for node_id in network.transitions[data_id1]:
            if node_id not in network.transitions.get(data_id2, []):
                if data_id2 not in network.transitions:
                    network.transitions[data_id2] = []
                network.transitions[data_id2].append(node_id)

    def _split_list(self, network, data_id1, data_id2, attr_name,
                    user_attributes):
        # Changes network in place.  data1 and data2 are both a ENUMs.

        data1 = network.nodes[data_id1]
        data2 = network.nodes[data_id2]
        DATA1_VALUE = data1.attributes[attr_name]
        DATA2_VALUE = data2.attributes[attr_name]
        DATA1_TYPE = _get_attribute_type(DATA1_VALUE)
        DATA2_TYPE = _get_attribute_type(DATA2_VALUE)
        assert DATA1_TYPE == TYPE_ENUM and DATA2_TYPE == TYPE_ENUM
        COMMON_VALUE = _intersection(DATA1_VALUE, DATA2_VALUE)
        assert len(COMMON_VALUE) < len(DATA1_VALUE)
        assert len(COMMON_VALUE) < len(DATA2_VALUE)

        # Remove the common values from the two ENUMs.
        DATA1_VALUE = DATA1_VALUE[:]
        DATA2_VALUE = DATA2_VALUE[:]
        for value in COMMON_VALUE:
            i1 = DATA1_VALUE.index(value)
            i2 = DATA2_VALUE.index(value)
            DATA1_VALUE.pop(i1)
            DATA2_VALUE.pop(i2)
        if len(DATA1_VALUE) == 1:
            DATA1_VALUE = DATA1_VALUE[0]
        if len(DATA2_VALUE) == 1:
            DATA2_VALUE = DATA2_VALUE[0]
        data1.attributes[attr_name] = DATA1_VALUE
        data2.attributes[attr_name] = DATA2_VALUE

        # Make a new DataNode object that contains the common values.
        attributes = data1.attributes.copy()
        attributes[attr_name] = COMMON_VALUE
        data3 = DataNode(data1.datatype, **attributes)
        network.nodes.append(data3)
        data_id3 = len(network.nodes) - 1

        # Every module that pointed to data_id1 or data_id2 should now
        # also point to data_id3.
        x1 = _backchain_to_ids(network, data_id1)
        x2 = _backchain_to_ids(network, data_id2)
        module_ids = x1 + x2
        for node_id in module_ids:
            if data_id3 not in network.transitions[node_id]:
                network.transitions[node_id].append(data_id3)

        # data_id3 should now point to everywhere that data_id1 and
        # data_id2 point to, if data_id3 is compatible with the
        # module.
        x1 = network.transitions.get(data_id1, [])
        x2 = network.transitions.get(data_id2, [])
        #module_ids = sorted({}.fromkeys(x1 + x2))
        module_ids = {}.fromkeys(x1 + x2)
        for module_id in module_ids:
            if _can_module_in_network_take_data(
                network, module_id, [data3], user_attributes):
                if data_id3 not in network.transitions:
                    network.transitions[data_id3] = []
                network.transitions[data_id3].append(module_id)


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

    def optimize(self, network, user_attributes):
        import copy

        bad_transitions = {}  # (node_id, next_id) -> 1
        for (node_id, next_ids) in network.iterate(node_class=ModuleNode):
            module = network.nodes[node_id]
            for next_id in next_ids:
                node = network.nodes[next_id]
                assert isinstance(node, DataNode)
                if not _can_module_produce_data(module, node):
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


class _OptimizeMergeNodes:
    # Sometimes the inference can lead to two nodes that share the
    # same parents and the same children, and almost the same
    # attributes.  For example:
    # Node1                   Node2
    # preprocess="unknown"    preprocess=<everything else>
    #
    # If this happens, merge them to simplify the network.
    def __init__(self):
        pass

    def optimize(self, network, user_attributes):
        while True:
            similar = self._find_similar_nodes(network)
            if not similar:
                break
            # Merge the similar nodes.
            while similar:
                n1, n2 = similar.pop()
                network = self._merge_nodes(network, n1, n2)
                similar = _fix_node_ids_after_merge(n1, n2, similar)
        return network

    def _find_similar_nodes(self, network):
        # Return a list of (node_id1, node_id2).  Can be empty.
        nodeid2previds = _make_backchain_dict(network)
        
        data_node_ids = [
            node_id for (node_id, node) in enumerate(network.nodes)
            if isinstance(node, DataNode)]
        
        similar = []
        for i, node_id1 in enumerate(data_node_ids):
            for node_id2 in data_node_ids[i+1:]:
                if self._are_nodes_similar(
                    network, node_id1, node_id2, nodeid2previds):
                    similar.append((node_id1, node_id2))
        return similar
    
        ## for node_id1 in range(len(network.nodes)):
        ##     node_1 = network.nodes[node_id1]
        ##     if not isinstance(node_1, DataNode):
        ##         continue
        ##     for node_id2 in range(node_id1+1, len(network.nodes)):
        ##         node_2 = network.nodes[node_id2]
        ##         if not isinstance(node_2, DataNode):
        ##             continue
        ##         if self._are_nodes_similar(
        ##             network, node_id1, node_id2, nodeid2previds):
        ##             return (node_id1, node_id2)
        ## return None

    def _are_nodes_similar(self, network, node_id1, node_id2, nodeid2previds):
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
        p1 = nodeid2previds.get(node_id1, [])
        p2 = nodeid2previds.get(node_id2, [])
        if len(p1) != len(p2):
            return False
        if sorted(p1) != sorted(p2):
            return False

        # They must share all but 1 attribute.
        x, diff_attrs = _compare_two_nodes(node_1, node_2)
        if len(diff_attrs) > 1:
            return False
        return True

    def _merge_nodes(self, network, node_id1, node_id2):
        # Since these nodes share the same parents, I can just delete
        # one of these nodes.

        # Delete the one with the higher node_id (node_id2) to disturb the
        # network as little as possible.
        if node_id1 > node_id2:
            node_id1, node_id2 = node_id2, node_id1

        # Merge the attributes of the nodes.
        n1 = network.nodes[node_id1]
        n2 = network.nodes[node_id2]
        
        for name in n1.attributes:
            V1 = n1.attributes[name]
            V2 = n2.attributes[name]
            if V1 == V2:
                continue
            T1 = _get_attribute_type(V1)
            T2 = _get_attribute_type(V2)
            case = _assign_case_by_type(T1, T2)

            if case == 1:
                V = [V1, V2]
            elif case == 2:
                V = V2
                if V1 not in V:
                    V.append(V1)
            elif case == 3:
                V = V1
                if V2 not in V:
                    V.append(V2)
            else:
                V = V1
                for x in V2:
                    if x not in V:
                        V.append(x)
            n1.attributes[name] = V
        x = network.delete_node(node_id2)
        return x


def _compare_two_nodes(node1, node2):
    # Return a tuple (has_same_datatype, list of attributes with
    # different values).  If the two nodes have a different datatype,
    # then the attributes will be None.

    if node1.datatype.name != node2.datatype.name:
        return (False, None)

    attrs = []
    for name in node1.attributes:
        V1 = node1.attributes[name]
        V2 = node2.attributes[name]
        T1 = _get_attribute_type(V1)
        T2 = _get_attribute_type(V2)
        case = _assign_case_by_type(T1, T2)

        if case == 1:
            if V1 != V2:
                attrs.append(name)            
        elif case in [2, 3]:  # different types
            attrs.append(name)
        elif case == 4:
            if sorted(V1) != sorted(V2):
                attrs.append(name)
        else:
            raise AssertionError
    return (True, attrs)



def select_start_node(network, start_data, user_attributes):
    # TODO: Rename this.  Actually prunes network based on start node.
    # start_data may be a single DataNode object or a list of DataNode
    # objects.  DataTypes are also allowed in lieu of DataNode
    # objects.

    # Strategy:
    # 1.  Include all nodes that can reach both a start and end node.
    # 2.  Remove modules that have no inputs.
    # 3.  Repeat steps 1-2 until convergence.
    start_ids = _find_start_nodes(network, start_data)
    good_ids = {}.fromkeys(range(len(network.nodes)))

    while good_ids:
        # If any of the nodes can't reach both the start and the end,
        # then delete it.
        delete_ids = {}
        for node_id in good_ids:
            good_by_bc = _can_reach_by_bc(network, node_id, good_ids)
            good_by_fc = _can_reach_by_fc(network, node_id, good_ids)
            # If it can't reach the goal, then delete it.
            if 0 not in good_by_fc:
                delete_ids[node_id] = 1
            # If it can't reach any starts, then delete it.
            x = [x for x in start_ids if x in good_by_bc]
            if not x:
                delete_ids[node_id] = 1

        # If a module requires multiple inputs, make sure each of the
        # inputs are in the network.  If not, remove the module.
        for node_id in good_ids:
            node = network.nodes[node_id]
            if not isinstance(node, ModuleNode):
                continue
            assert len(node.in_datatypes) > 0
            if len(node.in_datatypes) <= 1:
                continue
            in_ids = _backchain_to_ids(network, node_id)
            combos = _get_valid_input_combinations(
                network, node_id, in_ids, user_attributes)
            # Delete the IDs in in_ids if they aren't a part of a
            # combination that's all good.
            in_good_combo = {}
            for combo in combos:
                x = [x for x in combo if x in good_ids]
                if len(x) == len(combo):
                    # All good.
                    in_good_combo.update({}.fromkeys(x))
            # Figure out which nodes to delete.
            x = in_ids
            x = [x for x in x if x not in in_good_combo]
            x = [x for x in x if x in good_ids]
            delete_ids.update({}.fromkeys(x))
            #nid = [x for x in in_ids if x in good_ids]
            #if len(nid) != len(in_ids):
            #    delete_ids[node_id] = 1

        for node_id in delete_ids:
            if node_id in good_ids:
                del good_ids[node_id]
        if not delete_ids:
            break

    # Delete all the IDs that aren't in good_ids.
    bad_ids = [x for x in range(len(network.nodes)) if x not in good_ids]
    network = network.delete_nodes(bad_ids)

    # This can generate networks that contain internal nodes that
    # can't be generated.
    # DATA1 -> MODULE1 -> DATA2
    # DATA3 -> MODULE1 -> DATA4
    # If DATA1 is deleted, DATA2 will still be in the network, even
    # though it can't be generated anymore.  Should go through and
    # remove all internal nodes that can't be generated anymore.
    #
    # The problem with keeping them in there, is that it's difficult
    # to know whether a network can be created or not.

    # Check each internal DataNode and delete if it can't be
    # generated.
    node2previds = _make_backchain_dict(network)
    bad_ids = []
    for node_id in range(len(network.nodes)):
        if not isinstance(network.nodes[node_id], DataNode):
            continue
        # in_ids -> module_id -> node_id    ID
        # in_datas -> module -> out_data    objects
        module_ids = node2previds.get(node_id, [])
        if not module_ids:  # not internal node
            continue
        # Make sure there is a combination of inputs that can generate
        # this output.
        is_reachable = False
        for module_id in module_ids:
            all_in_ids = node2previds[module_id]
            combos = _get_valid_input_combinations(
                network, module_id, all_in_ids, user_attributes)
            for combo in combos:
                module = network.nodes[module_id]
                in_datas = [network.nodes[x] for x in combo]
                out_data = network.nodes[node_id]
                if _can_module_take_data(
                    module, in_datas, out_data, user_attributes):
                    is_reachable = True
                    break
            if is_reachable:
                break
        if not is_reachable:
            bad_ids.append(node_id)
    network = network.delete_nodes(bad_ids)
    #network = optimize_network(network, user_attributes)
    return network


def find_node(network, node):
    # Return a node_id or None.
    assert isinstance(node, DataNode)

    found = []
    for i in range(len(network.nodes)):
        net_node = network.nodes[i]
        if not isinstance(net_node, DataNode):
            continue
        if _is_data_compatible_with(node, net_node):
            found.append(i)
    assert len(found) <= 1
    if not found:
        return None
    return found[0]


## def remove_data_node(network, *attributes):
##     # Look for the nodes that are compatible with the attributes.  If
##     # an attribute has multiple values, DataNode that matches any of
##     # those values will be removed.
##     node_ids = []  # list of node_ids to remove
##     for node_id, next_ids in network.iterate(node_class=DataNode):
##         node = network.nodes[node_id]
##         if not isinstance(node, DataNode):
##             continue

##         # If compatible, remove or alter this node.
##         remove_attr = {}  # name -> list of values
##         for attr in attributes:
##             assert isinstance(attr, Attribute)
##             if node.datatype != attr.datatype:
##                 break
##             assert attr.name in node.attributes
##             assert attr.name not in remove_attr

##             DATA_VALUE = node.attributes[attr.name]
##             ATTR_VALUE = attr.value
##             DATA_TYPE = _get_attribute_type(DATA_VALUE)
##             ATTR_TYPE = _get_attribute_type(ATTR_VALUE)
##             case = _assign_case_by_type(DATA_TYPE, ATTR_TYPE)

##             if case == 1:
##                 if DATA_VALUE != ATTR_VALUE:
##                     break
##                 remove_attr[attr.name] = [DATA_VALUE]
##             elif case == 2:
##                 # DATA is ATOM, ATTR is ENUM.
##                 if DATA_VALUE not in ATTR_VALUE:
##                     break
##                 remove_attr[attr.name] = [DATA_VALUE]
##             elif case == 3:
##                 # DATA is ENUM, ATTR is ATOM
##                 if ATTR_VALUE not in DATA_VALUE:
##                     break
##                 remove_attr[attr.name] = [ATTR_VALUE]
##             else:
##                 # DATA is ENUM, ATTR is ATOM
##                 x = _intersection(DATA_VALUE, ATTR_VALUE)
##                 if not x:
##                     break
##                 remove_attr[attr.name] = x
##         else:
##             # Remove or alter this node.
##             for name in remove_attr:
##                 assert name in node.attributes
##                 values = node.attributes[name][:]
##                 if type(values) is type(""):
##                     values = [values]
##                 for value in remove_attr[name]:
##                     i = values.index(value)
##                     values.pop(i)
##                 if not values:
##                     # Removed everything, delete this node.
##                     node_ids.append(node_id)
##                     break
##                 if len(values) == 1:
##                     values = values[0]
##                 node.attributes[name] = values

##     # Remove each of these node_ids from the network.
##     network = network.delete_nodes(node_ids)
##     network = _OptimizeNoDanglingNodes().optimize(network)
##     return network


def get_input_nodes(
    network, user_attributes, skip_datatypes=None,
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
        network, user_attributes, max_nodes=max_inputs, nodeid2id_fn=fn)
    return nodeid_combos


def _product_network(
    network, user_attributes, max_nodes=None, nodeid2id_fn=None):
    # Perform a product operation over every node in the network.
    # max_nodes is the maximum number of nodes that I should perform a
    # product over.
    # Return a list of tuples.
    if max_nodes is None:
        max_nodes = len(network.nodes)

    nodeid2previds = _make_backchain_dict(network)
    cache = {}
    x = _product_network_h(
        network, 0, nodeid2previds, user_attributes, max_nodes, nodeid2id_fn,
        cache)
    return x.keys()


def _product_network_h(
    network, node_id, nodeid2previds, user_attributes, max_nodes, nodeid2id_fn,
    cache):
    if node_id not in cache:
        x = _product_network_hh(
            network, node_id, nodeid2previds, user_attributes, max_nodes,
            nodeid2id_fn, cache)
        cache[node_id] = x
    return cache[node_id]


def _product_network_hh(
    network, node_id, nodeid2previds, user_attributes, max_nodes, nodeid2id_fn,
    cache):
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
        for previd in nodeid2previds.get(node_id, []):
            x = _product_network_h(
                network, previd, nodeid2previds, user_attributes, max_nodes,
                nodeid2id_fn, cache)
            inputs.update(x)
    elif isinstance(node, ModuleNode):
        # Find all potential sets of DataNode notes that can feed into me.
        assert node_id in nodeid2previds
        # list of tuples
        combos = _get_valid_input_combinations(
            network, node_id, nodeid2previds[node_id], user_attributes,
            node2previds=nodeid2previds)

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
                    network, x, nodeid2previds, user_attributes, max_nodes,
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


def get_input_datatypes(
    network, user_attributes, skip_datatypes=None,
    skip_private_datatypes=False, max_inputs=None):
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
        network, user_attributes, max_nodes=max_inputs, nodeid2id_fn=fn)

    # Convert from ID back to datatypes.
    datatype_combos = []
    for dtid_combo in dtid_combos:
        names = [fn.id2name[x] for x in dtid_combo]
        x = tuple([fn.name2datatype[x] for x in names])
        datatype_combos.append(x)
    return datatype_combos


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
    import sys

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
    
    node_name = "%s [%d]" % (node.datatype.name, node_id)

    LINE_WIDTH = 60
    
    lines = []
    w = lines.append
    w("<B>%s</B>" % node_name)

    if node.attributes:
        w('<BR/>')
        w('<BR/>')
        w("<U>Data Attributes</U>:")
        w('<BR ALIGN="LEFT"/>')
    for x in node.attributes.iteritems():
        name, value = x
        x = "%s = %s" % (name, value)
        for x in parselib.linesplit(x, prefix1=0, prefixn=4, width=LINE_WIDTH):
            w(x)
            w('<BR ALIGN="LEFT"/>')
    return "<%s>" % "".join(lines)


def _format_modulenode_gv(node, node_id, options):
    node_name = "%s [%d]" % (node.name, node_id)

    if options is None:
        options = {}

    lines = []
    w = lines.append
    w("<B>%s</B>" % node_name)
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
    filename, network, options=None, highlight_node_ids=[], verbose=False):
    from genomicode import graphviz

    gv_nodes = []
    gv_edges = []
    gv_node2attr = {}

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
    highlight_color = "#FDAF91"
    # Brewer rdbu-div.  Too much contrast.
    #data_color = "#F5A582"
    #module_color = "#92C6DF"
    # Brewer piyg-div.  Too pastel
    #data_color = "#F2B7DB"
    #module_color = "#B8E286"

    id2name = {}
    for node_id, node in enumerate(network.nodes):
        node2attr = {}
        node2attr["style"] = "filled"
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
        if node_id in highlight_node_ids:
            node2attr["fillcolor"] = highlight_color

        node_name = "%s [%d]" % (name, node_id)
        id2name[node_id] = node_name
        gv_nodes.append(node_name)
        gv_node2attr[node_name] = node2attr

    for node_id, next_ids in network.transitions.iteritems():
        for next_id in next_ids:
            x1 = id2name[node_id]
            x2 = id2name[next_id]
            gv_edges.append((x1, x2))

    G = graphviz.make_graph(
        gv_nodes, gv_edges, node2attributes=gv_node2attr,
        prog="dot", directed=True)
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


##def write_network(file_or_handle, network):
##    import pickle
##
##    handle = file_or_handle
##    if type(handle) is type(""):
##        handle = open(file_or_handle, 'w')
##    pickle.dump(network, handle)
##

##def read_network(file_or_handle):
##    import pickle
##
##    handle = file_or_handle
##    if type(handle) is type(""):
##        handle = open(file_or_handle, 'r')
##    return pickle.load(handle)


def diagnose_start_node(network, user_data, outhandle=None):
    import sys

    outhandle = outhandle or sys.stdout
    
    results = _score_start_nodes(network, user_data)

    # TODO: Should prettify this output.

    # This can happen if the network doesn't contain any of the same
    # datatypes as user_data.  This can be if the network is really
    # small, e.g. if there's no way at all to generate the desired out
    # data type.
    if not results:
        node_or_nodes = "node"
        if len(network.nodes) > 1:
            node_or_nodes = "nodes"
        print >>outhandle, \
              "The network (%d %s) contains no matching data types." % (
            len(network.nodes), node_or_nodes)
        return
    header = [
        "Mismatches", "Datatype", "Node ID", "Attribute", "Network", "User"]
    print >>outhandle, "\t".join(header)
    for x in results:
        score, node_id, user_data, attr_values = x
        dt_name = network.nodes[node_id].datatype.name
        if not attr_values:
            x = score, dt_name, node_id, "", "", ""
            assert len(x) == len(header)
            print >>outhandle, "\t".join(map(str, x))
        for name, netw_value, user_value in attr_values:
            x = score, dt_name, node_id, name, netw_value, user_value
            assert len(x) == len(header)
            print >>outhandle, "\t".join(map(str, x))


def _backchain_to_modules(moduledb, data):
    # Return list of modules that can generate an output that is
    # compatible with data.

    modules = []  # list of (module, num compatible attributes)
    for module in moduledb:
        if _can_module_produce_data(module, data):
            modules.append(module)
    return modules


def _backchain_to_input(module, in_num, out_data, user_attributes,
                        force_default_input_attribute_to_be_all_values=False):
    # Given a module and output_data, return the input_data object
    # that can generate the output.  This should only be called by
    # _backchain_to_all_inputs.  The SAME_AS constraint is handled
    # there.  If this is called by itself, the attributes might not be
    # right.
    #
    # force_default_input_attribute_to_be_all_values can be useful
    # when checking for the possibility that an input node can go into
    # a module.  However, it should be False when generating the
    # network to prevent combinatorial explosion.
    assert in_num < len(module.in_datatypes)

    in_datatype = module.in_datatypes[in_num]
    #out_datatype = out_data.datatype

    # Can't generate debug messages here because the SAME_AS
    # constraints aren't handled in this function.

    # The attributes for the input object should come from (in
    # decreasing priority):
    # 1.  Consequence (i.e. SAME_AS_CONSTRAINT).
    # 2.  Constraint.
    # 3.  user attribute
    # 4.  out_data              (default_attributes_from)
    # 5.  default output value of input datatype
    #
    # If the user attribute conflicts with a Constraint, the
    # Constraint is higher priority to make sure no objects are
    # generated that the module cannot handle.  However, if the user
    # attribute or default value is a part of the constraint,
    # (e.g. one of many options), then we can refine it with the user
    # attribute (or default).
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

    # Set the attributes in increasing order of priority.  Higher
    # priority overwrites lower priority.

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
        attributes[attrdef.name] = default
        attrsource[attrdef.name] = "default"

    # Case 4.  If default_attributes_from is the same as in_num, then
    # fill with the same values as the out_data.
    indexes = [x.input_index for x in module.default_attributes_from]
    if in_num in indexes:
        for name, value in out_data.attributes.iteritems():
            attributes[name] = value
            attrsource[name] = "out_data"

    # Case 3.  If the input data object does not proceed to the output
    # data object, then set the attribute provided by the user.
    x = [x for x in module.default_attributes_from if x.input_index == in_num]
    if not x:
        # Set values from user attributes.
        for attr in user_attributes:
            # Ignore attributes for other data types.
            if attr.datatype.name != in_datatype.name:
                continue

            ### Make sure this value doesn't conflict with a higher
            ### priority value.
            ##old_value = attributes.get(attr.name, attr.value)
            ##old_type = _get_attribute_type(old_value)
            ##
            ### Make sure the value is compatible with the data being
            ### created.
            ##if old_type == TYPE_ATOM:
            ##    msg = (
            ##        "ModuleNode %s requires a %s with %s=%s, but user requests "
            ##        "it to be %s." % (
            ##            module.name, in_datatype.name, attr.name, old_value,
            ##            attr.value))
            ##    assert attr.value == old_value, msg
            ##elif old_type == TYPE_ENUM:
            ##    assert attr.value in old_value, "incompatible"
            ##else:
            ##    raise AssertionError
            attributes[attr.name] = attr.value
            attrsource[attr.name] = "user"

    # Case 2.  Set the attributes based on the constraints.
    for constraint in module.constraints:
        if constraint.input_index != in_num:
            continue

        if constraint.behavior == MUST_BE:
            attributes[constraint.name] = constraint.arg1
            attrsource[constraint.name] = "constraint"
        elif constraint.behavior == CAN_BE_ANY_OF:
            value = constraint.arg1
            source = "constraint"

            # Refine by the user attribute or default value.
            if attrsource.get(constraint.name) in ["user", "default"]:
                x = _get_attribute_type(attributes[constraint.name])
                if x == TYPE_ATOM:
                    if attributes[constraint.name] in value:
                        value = attributes[constraint.name]
                        source = "constraint,%s" % attrsource[constraint.name]
                elif x == TYPE_ENUM:
                    x = _intersection(attributes[constraint.name], value)
                    if x:
                        value = x
                        source = "constraint,%s" % attrsource[constraint.name]
                else:
                    raise AssertionError
            attributes[constraint.name] = value
            attrsource[constraint.name] = source
        elif constraint.behavior == SAME_AS:
            # Handled in _backchain_to_all_inputs.
            pass
        else:
            raise AssertionError

    # Case 1.  Set the attributes based on the consequences.  If there
    # is a Consequence that is SAME_AS_CONSTRAINT, then the attribute
    # should be determined by the out_data.  e.g.
    # Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"])
    # Consequence("quantile_norm", SAME_AS_CONSTRAINT)
    #
    # The module takes anything, produces the same value.  So the
    # backchainer needs to preserve the value from the out_data.

    # Get SAME_AS_CONSTRAINT.  Nothing to do for SET_TO,
    # SET_TO_ONE_OF, BASED_ON_DATA.
    x = [x for x in module.consequences if x.behavior == SAME_AS_CONSTRAINT]
    # Get the consequences that are based on this datatype.
    x = [x for x in x if x.arg1 == in_num]
    consequences = x
    for consequence in consequences:
        n = consequence.name
        source = "consequence"

        # Copy the value from the output data.
        data_value = out_data.attributes[n]
        data_type = _get_attribute_type(data_value)

        # Since more values may be allowed in the consequence,
        # further refine based on the constraint.  E.g.:
        # Constraint  [A, B]
        # Consequence [A, B, C, D]
        x = [x for x in module.constraints if x.name == n]
        x = [x for x in x if x.input_index == in_num]
        assert len(x) > 0
        assert len(x) == 1
        constraint = x[0]
        if constraint.behavior == SAME_AS:
            constraint = _resolve_constraint(constraint, module.constraints)
        if constraint.behavior == MUST_BE:
            # constraint.arg1  <value>
            if data_type == TYPE_ATOM:
                assert constraint.arg1 == data_value
            elif data_type == TYPE_ENUM:
                assert constraint.arg1 in data_value
                data_value = constraint.arg1
                source = "consequence+constraint"
            else:
                raise AssertionError
        elif constraint.behavior == CAN_BE_ANY_OF:
            # constraint.arg1  list of <values>
            if data_type == TYPE_ATOM:
                assert data_value in constraint.arg1
            elif data_type == TYPE_ENUM:
                common = _intersection(constraint.arg1, data_value)
                assert common
                data_value = common
                source = "consequence+constraint"
            else:
                raise AssertionError
        else:
            raise AssertionError

        attributes[n] = data_value
        attrsource[n] = source

    return attributes, attrsource


def _backchain_to_all_inputs(
    module, out_data, user_attributes,
    force_default_input_attribute_to_be_all_values=False):
    all_attributes = []
    all_attrsource = []
    force = force_default_input_attribute_to_be_all_values
    for in_num in range(len(module.in_datatypes)):
        x = _backchain_to_input(
            module, in_num, out_data, user_attributes,
            force_default_input_attribute_to_be_all_values=force)
        attributes, attrsource = x
        all_attributes.append(attributes)
        all_attrsource.append(attrsource)

    # Handle the SAME_AS constraints here.  Can't do in
    # _backchain_to_input because these constraints require the
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

        x = in_datatype.output(**attributes)
        all_inputs.append(x)

        # Optimization: Don't call debug_print and sorted.
        if not DEBUG:
            continue
        debug_print(
            "Backchaining %s (input=%d) -> %s -> %s." %
            (in_datatype.name, in_num, module.name, out_datatype.name))
        #debug_print("Generating a %s with attributes:" % in_datatype.name)
        for name in sorted(attributes):
            debug_print("  %s=%s (%s)" %
                        (name, attributes[name], attrsource[name]))

    return all_inputs


def _forwardchain_to_outputs(module, in_datas):
    # Generate a list of DataNode objects that can be generated from
    # module and in_datas.  Multiple objects can be generated because
    # the consequences can vary.  E.g. center_genes can set
    # gene_center to either "mean" or "median".  It can be either, but
    # must be one of them.
    import itertools

    # Check the input variables.
    assert len(module.in_datatypes) == len(in_datas), module.name
    for i in range(len(module.in_datatypes)):
        assert in_datas[i].datatype.name == module.in_datatypes[i].name

    # Assume that in_datas fulfill the constraints of the module.

    datatype = module.out_datatype
    attributes = {}

    # Set the attributes based on the default values.
    # Case 1: default_attributes_from is given.
    #         Use the attributes from this default.
    # Case 2: Fill with the default input values of the output
    #         datatype.

    # Case 1.
    if module.default_attributes_from:
        # If there are multiple default_attributes_from, just use the
        # first one.
        # BUG: Is this always the right thing to do?
        input_index = module.default_attributes_from[0].input_index
        assert input_index < len(in_datas)
        data = in_datas[input_index]
        assert data.datatype.name == datatype.name
        attributes.update(data.attributes)
    # Case 2.
    else:
        for attrdef in datatype.attribute_defs.itervalues():
            attributes[attrdef.name] = attrdef.default_in

    # Set the attributes based on the consequences of the module.
    possibilities = {}

    for cons in module.consequences:
        if cons.behavior == SET_TO:
            attributes[cons.name] = cons.arg1
        elif cons.behavior == SET_TO_ONE_OF:
            possibilities[cons.name] = cons.arg1
        elif cons.behavior == BASED_ON_DATA:
            possibilities[cons.name] = cons.arg1
        elif cons.behavior == SAME_AS_CONSTRAINT:
            input_index = cons.arg1
            data = in_datas[input_index]
            attributes[cons.name] = data.attributes[cons.name]
        else:
            raise AssertionError

    # If no possibilities, then make one output variable.
    if not possibilities:
        x = DataNode.__new__(DataNode)
        x.datatype = datatype
        x.attributes = attributes.copy()
        return [x]

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

    return outputs


def _backchain_to_ids(network, node_id):
    # Return a list of IDs that point to this node_id.
    assert node_id < len(network.nodes)
    ids = {}
    for nid, next_ids in network.transitions.iteritems():
        if node_id in next_ids:
            ids[nid] = 1
    #return sorted(ids)
    return ids.keys()


def _make_backchain_dict(network):
    # Return a dictionary of node_id -> prev_node_ids
    nodeid2previds = {}
    for prev_id, node_ids in network.transitions.iteritems():
        for node_id in node_ids:
            if node_id not in nodeid2previds:
                nodeid2previds[node_id] = []
            nodeid2previds[node_id].append(prev_id)
    return nodeid2previds

## def _can_module_take_one_data(module, input_num, data):
##     if data.datatype != module.in_datatypes[input_num]:
##         return False

##     # Make sure the data satisfies each of the module's constraints.
##     for constraint in module.constraints:
##         # If the constraint does not apply to this data object,
##         # then ignore.
##         if constraint.input_index != input_num:
##             continue

##         assert constraint.name in data.attributes
##         data_value = data.attributes.get(constraint.name)
##         data_type = _get_attribute_type(data_value)
##         assert data_type in [TYPE_ATOM, TYPE_ENUM]

##         # If this is a SAME_AS constraint, follow SAME_AS.  This value
##         # must be the same as the value for another input data.  If
##         # there is a constraint on that input data, we might be able
##         # to figure out what it is.
##         # e.g.
##         # input_0   quantile_norm MUST_BE "no"
##         # input_1   quantile_norm SAME_AS input_0
##         #
##         # OR:
##         # input_0   quantile_norm MUST_BE "no"
##         # input_1   quantile_norm SAME_AS input_0
##         # input_2   quantile_norm SAME_AS input_1
##         #
##         while constraint and constraint.behavior == SAME_AS:
##             x_cons = [
##                 x for x in module.constraints
##                 if x.input_index == constraint.arg1 and
##                 x.name == constraint.name]
##             x_same = [x for x in x_cons if x.behavior == SAME_AS]
##             x_value = [x for x in x_cons
##                        if x.behavior in [MUST_BE, CAN_BE_ANY_OF]]
##             assert not (x_same and x_value)
##             if x_same:
##                 constraint = x_same[0]
##             elif x_value:
##                 assert len(x_value) == 1
##                 constraint = x_value[0]
##             else:
##                 # If the SAME_AS chain does not end with MUST_BE or
##                 # CAN_BE_ANY_OF, then we cannot test this.
##                 constraint = None
##         # SAME_AS did not lead to a testable constraint.
##         if not constraint:
##             continue

##         if constraint.behavior == MUST_BE:
##             if data_type == TYPE_ATOM:
##                 if data_value != constraint.arg1:
##                     return False
##             elif data_type == TYPE_ENUM:
##                 return False
##             else:
##                 raise AssertionError
##         elif constraint.behavior == CAN_BE_ANY_OF:
##             if data_type == TYPE_ATOM:
##                 if data_value not in constraint.arg1:
##                     return False
##             elif data_type == TYPE_ENUM:
##                 # data_value contains the possible values of this DataNode
##                 # object.  The values that are acceptable by module is
##                 # in constraint.arg1.  Make sure the module can handle
##                 # all of the possible values.
##                 if not _is_subset(data_value, constraint.arg1):
##                     return False
##             else:
##                 raise AssertionError
##         elif constraint.behavior == SAME_AS:
##             # Should not end up with SAME_AS constraints.
##             raise AssertionError
##         else:
##             raise AssertionError
##     return True

## def _can_module_take_data(module, datas):
##     # Return True/False if a module can take this list of DataNode nodes
##     # as an input.
##     if len(module.in_datatypes) != len(datas):
##         return False
##     for input_num in range(len(module.in_datatypes)):
##         data = datas[input_num]
##         if not _can_module_take_one_data(module, input_num, data):
##             return False

##         # Test SAME_AS constraints for this input_num.
##         x = module.constraints
##         x = [x for x in x if x.behavior == SAME_AS]
##         x = [x for x in x if x.input_index == input_num]
##         constraints = x
##         for constraint in constraints:
##             data_value = data.attributes.get(constraint.name)
##             # Should only be encountered if there are multiple input
##             # data types.
##             target_data = datas[constraint.arg1]
##             if data_value != target_data.attributes[constraint.name]:
##                 return False
##     return True


def _can_module_take_data(module, in_datas, out_data, user_attributes):
    # Return True/False if a module can take in_datas (list of
    # DataNode nodes) as input.

    assert len(in_datas) == len(module.in_datatypes)
    all_inputs = _backchain_to_all_inputs(
        module, out_data, user_attributes,
        force_default_input_attribute_to_be_all_values=True)
    for i in range(len(in_datas)):
        if not _is_data_compatible_with(in_datas[i], all_inputs[i]):
            return False
    return True


def _can_module_take_one_data(
    module, input_num, in_data, out_data, user_attributes):
    assert input_num < len(module.in_datatypes)

    # BUG: should not call _backchain_to_input!!!
    #good_in_data = _backchain_to_input(
    #    module, input_num, out_data, user_attributes)
    all_inputs = _backchain_to_all_inputs(
        module, out_data, user_attributes,
        force_default_input_attribute_to_be_all_values=True)
    return _is_data_compatible_with(in_data, all_inputs[input_num])


def _can_module_in_network_take_data(
    network, module_id, in_datas, user_attributes):
    module = network.nodes[module_id]
    # If in_datas is compatible with any of the out_datas, then return
    # True.
    out_data_ids = network.transitions.get(module_id, [])
    for out_data_id in out_data_ids:
        out_data = network.nodes[out_data_id]
        if _can_module_take_data(module, in_datas, out_data, user_attributes):
            return True
    return False


def _can_module_in_network_take_one_data(
    network, module_id, input_num, in_data, user_attributes):
    module = network.nodes[module_id]
    # If in_datas is compatible with any of the out_datas, then return
    # True.

    out_data_ids = network.transitions.get(module_id, [])
    for out_data_id in out_data_ids:
        out_data = network.nodes[out_data_id]
        if _can_module_take_one_data(
            module, input_num, in_data, out_data, user_attributes):
            return True
    return False


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


def _can_module_produce_data(module, data):
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
    # These rules need to match the policies in _backchain_to_input.

    # If this module doesn't produce the same data type, then it can't
    # produce this data object.
    if module.out_datatype.name != data.datatype.name:
        #debug_print(
        #    "ModuleNode can't generate data type: %s." % data.datatype.name)
        return False

    debug_print("Testing if module %s can produce data %s." %
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
                debug_print(msg)
                return False
        elif case == 2:
            # ModuleNode can produce any of a list of values.  Check
            # if the data's value can be produced by the module.
            if data_value not in outc_value:
                debug_print("Consequence %s conflicts." % consequence.name)
                return False
        elif case == 3:
            # ModuleNode produces a specific value.  DataNode could be
            # one of many values.
            if outc_value not in data_value:
                debug_print("Consequence %s conflicts." % consequence.name)
                return False
        elif case == 4:
            if not _intersection(data_value, outc_value):
                debug_print("Consequence %s conflicts." % consequence.name)
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
            "ModuleNode converts datatype.  Checking default attributes.")
        consequence_names = [x.name for x in module.consequences]
        for attrdef in module.out_datatype.attribute_defs.itervalues():
            # Ignore the attributes that have consequences.
            if attrdef.name in consequence_names:
                debug_print(
                    "Attr %r: Skipping--has consequence." % attrdef.name)
                continue
            assert attrdef.name in data.attributes
            data_value = data.attributes[attrdef.name]
            data_type = _get_attribute_type(data_value)
            assert data_type in [TYPE_ATOM, TYPE_ENUM]

            if data_type == TYPE_ATOM:
                if attrdef.default_in != data_value:
                    debug_print("Attr %r: Conflicts (module %r, data %r)." %
                                (attrdef.name, attrdef.default_in, data_value))
                    return False
            elif data_type == TYPE_ENUM:
                if attrdef.default_in not in data_value:
                    debug_print("Attr %r: Conflicts (module %r, data %r)." %
                                (attrdef.name, attrdef.default_in, data_value))
                    return False
            else:
                raise AssertionError
            debug_print("Attr %r: matches defaults." % attrdef.name)

    # TESTING.
    # If the module converts the datatype, the consequences don't
    # conflict, and the default attributes don't conflict, then this
    # should match.
    if not module.default_attributes_from:
        debug_print("Match because of converting datatype.")
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
            debug_print("Consequence '%s' matches." % consequence.name)
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
                debug_print("Consequence '%s' matches." % consequence.name)
                return True
        elif (const1.behavior, const2.behavior) == (MUST_BE, CAN_BE_ANY_OF):
            if const1.arg1 not in const2.arg1:
                debug_print("Consequence '%s' matches." % consequence.name)
                return True
        elif (const1.behavior, const2.behavior) == (CAN_BE_ANY_OF, MUST_BE):
            if const2.arg1 not in const1.arg1:
                debug_print("Consequence '%s' matches." % consequence.name)
                return True
        elif (const1.behavior, const2.behavior) == \
                 (CAN_BE_ANY_OF, CAN_BE_ANY_OF):
            if not _intersection(const1.arg1, const2.arg1):
                debug_print("Consequence '%s' matches." % consequence.name)
                return True
        else:
            raise AssertionError

    # No conflicts, and the module has no consequences.
    if not module.consequences:
        debug_print("Match because there are no consequences.")
        return True

    # No consequences match.
    debug_print("No consequences match.")
    return False


def _get_valid_input_combinations(
    network, module_id, all_input_ids, user_attributes, node2previds=None):
    # Given a list of all input IDs that point to a module, yield
    # tuples of input_ids that can match the input datatypes of the
    # module.  Only checks the datatypes, and does not do any other
    # checks to make sure the inputs are valid.
    import itertools

    # In about 97% of the cases, the module only has 1 datatype.  In
    # ~90% of the cases, there are 2 input IDs.
    module = network.nodes[module_id]

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
            x = [x for x in x if _can_module_in_network_take_one_data(
                network, module_id, i, network.nodes[x], user_attributes)]
            args[i] = x

    valid = []
    # Optimization: Assume existing inputs in the network are valid
    # and don't check them again.
    if len(module.in_datatypes) == 1:
        if node2previds:
            ids = node2previds.get(module_id, [])
        else:
            ids = _backchain_to_ids(network, module_id)
        x = [(x, ) for x in ids]
        valid.extend(x)
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
        # No duplicated IDs.
        #x = {}.fromkeys(input_ids)
        #if len(x) != len(input_ids):
        #    continue

        # Make sure the inputs are compatible with the module.
        input_datas = [network.nodes[x] for x in input_ids]
        if not _can_module_in_network_take_data(
            network, module_id, input_datas, user_attributes):
            continue

        # Make sure the outputs are compatible with the module.
        output_datas = _forwardchain_to_outputs(module, input_datas)

        output_is_compatible = False
        node_ids = network.transitions.get(module_id, [])
        for x in itertools.product(node_ids, output_datas):
            node_id, output_data = x
            node = network.nodes[node_id]

            if _is_data_compatible_with(output_data, node):
                output_is_compatible = True
                break
        if not output_is_compatible:
            continue
        valid.append(input_ids)
    return valid


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
        ids = _backchain_to_ids(network, nid)
        stack.extend(ids)
    return reachable_ids


def _can_reach_by_fc(network, node_id, good_ids):
    # Return a dictionary of all the node IDs that can be reached by
    # forward chaining.
    reachable_ids = {}
    stack = [node_id]
    while stack:
        nid = stack.pop(0)
        if nid in reachable_ids:
            continue
        if nid not in good_ids:
            continue
        reachable_ids[nid] = 1
        ids = network.transitions.get(nid, [])
        stack.extend(ids)
    return reachable_ids


def _find_data_node(nodes, node):
    assert isinstance(node, DataNode)

    # All values need to be exactly equal.
    #
    # CASE   N1_TYPE    N2_TYPE    RESULT
    #   1      ATOM       ATOM     OK if ATOM equal.
    #   2      ATOM       ENUM     No.
    #   3      ENUM       ATOM     No.
    #   4      ENUM       ENUM     OK if ENUM equal.

    for i, n in enumerate(nodes):
        if not isinstance(n, DataNode):
            continue
        if n.datatype.name != node.datatype.name:
            continue

        n1_attr = n.attributes
        n2_attr = node.attributes
        # Attributes should be the same.
        #x1 = n1_attr.keys()
        #x2 = n2_attr.keys()
        #all_attributes = {}.fromkeys(x1 + x2)

        is_equal = True
        for key in n1_attr:
            assert key in n1_attr
            assert key in n2_attr
            N1_VALUE = n1_attr[key]
            N2_VALUE = n2_attr[key]
            N1_TYPE = _get_attribute_type(N1_VALUE)
            N2_TYPE = _get_attribute_type(N2_VALUE)
            case = _assign_case_by_type(N1_TYPE, N2_TYPE)

            attr_equal = False
            if case == 1:
                if N1_VALUE == N2_VALUE:
                    attr_equal = True
            elif case == 4:
                if sorted(N1_VALUE) == sorted(N2_VALUE):
                    attr_equal = True
            if not attr_equal:
                is_equal = False
        if is_equal:
            return i
    return -1


def _find_start_nodes(network, user_data):
    scores = _score_start_nodes(network, user_data)
    x = [x for x in scores if x[0] == 0]
    node_ids = [x[1] for x in x]
    return node_ids


def _score_start_nodes(network, user_data):
    # Return a list of (score, node_id, user_data, list of (name,
    # value in network node, value in user_data) of incompatible
    # attributes).  Sorted by increasing score.
    import operator

    # Make a list of all the desired start nodes, as DataNode objects.
    user_datas = user_data
    if not operator.isSequenceType(user_data):
        user_datas = [user_data]
    for i, x in enumerate(user_datas):
        if isinstance(x, DataType):
            x = x.input()  # convert to DataNode
        assert isinstance(x, DataNode)
        user_datas[i] = x

    # Look for the nodes in the network that are compatible with
    # user_data.
    results = []
    for node_id, next_ids in network.iterate(node_class=DataNode):
        for user_data in user_datas:
            netw_data = network.nodes[node_id]
            if netw_data.datatype.name != user_data.datatype.name:
                continue

            attrs = _score_one_start_node(netw_data, user_data)
            attr_values = []
            for attr in attrs:
                x1 = network.nodes[node_id].attributes[attr]
                x2 = user_data.attributes[attr]
                x = attr, x1, x2
                attr_values.append(x)
            x = len(attrs), node_id, user_data, attr_values
            results.append(x)
            
    return sorted(results)


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


def _score_one_start_node(network_data, user_data):
    # network_data is a DataNode node in the network.  user_data is
    # the DataNode that the user wants to start on.  Return a list of
    # <attributes> that are incompatible between the network_data and
    # user_data.  An empty list indicates that the two nodes are
    # compatible.
    assert network_data.datatype.name == user_data.datatype.name

    # Start DataNode
    # ATOM      Must match a specific value.
    # ENUM      UNDEFINED.
    #
    # CASE  USER_TYPE  NETW_TYPE   RESULT
    #   1      ATOM       ATOM     Check if items are equal.
    #   2      ATOM       ENUM     Check if ATOM in ENUM.
    #   3      ENUM       ATOM     NotImplementedError.
    #   4      ENUM       ENUM     NotImplementedError.

    netw_attr = network_data.attributes
    user_attr = user_data.attributes

    attributes = []
    for key in user_attr:
        assert key in netw_attr
        assert key in user_attr
        NETW_VALUE = netw_attr[key]
        USER_VALUE = user_attr[key]
        NETW_TYPE = _get_attribute_type(NETW_VALUE)
        USER_TYPE = _get_attribute_type(USER_VALUE)
        CASE = _assign_case_by_type(USER_TYPE, NETW_TYPE)

        if CASE == 1:
            if NETW_VALUE != USER_VALUE:
                attributes.append(key)
        elif CASE == 2:
            if USER_VALUE not in NETW_VALUE:
                attributes.append(key)
        elif CASE in [3, 4]:
            raise NotImplementedError
        else:
            raise AssertionError
    return attributes


def _is_data_compatible_with(data_specific, data_general):
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


def _assign_case_by_type(type1, type2):
    types = [(TYPE_ATOM, TYPE_ATOM),
             (TYPE_ATOM, TYPE_ENUM),
             (TYPE_ENUM, TYPE_ATOM),
             (TYPE_ENUM, TYPE_ENUM), ]
    x = (type1, type2)
    assert x in types, "Unknown types: %s %s" % (type1, type2)
    i = types.index(x)
    return i + 1


import types
def _get_attribute_type(name):
    t = type(name)
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
    raise AssertionError, "Unknown attribute type: %s" % str(name)


def _intersection(x, y):
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
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)


def _print_nothing(s):
    pass


def _print_string(s):
    print s


def _print_line(line, prefix1, prefixn, width, outhandle=None):
    import sys
    from genomicode import parselib

    outhandle = outhandle or sys.stdout
    lines = parselib.linesplit(
        line, prefix1=prefix1, prefixn=prefixn, width=width)
    for x in lines:
        print >>outhandle, x

    #lines = []
    #p = prefix0
    #while line:
    #    x = line[:(width - len(p))]
    #    line = line[len(x):]
    #    lines.append(x)
    #    p = prefixn

    #for i in range(len(lines)):
    #    p = prefixn
    #    if i == 0:
    #        p = prefix0
    #    print >> outhandle, "%s%s" % (p, lines[i])


## def _split_line(line, prefix0, prefixn, width):
##     # Return a list of lines.  prefix0 is a string to put in front of
##     # the first line.  prefixn is what to put in front of the
##     # remaining lines.
##     lines = []
##     p = prefix0
##     while line:
##         x = line[:(width - len(p))]
##         line = line[len(x):]
##         lines.append(x)
##         p = prefixn

##     for i in range(len(lines)):
##         p = prefixn
##         if i == 0:
##             p = prefix0
##         lines[i] = p + lines[i]
##     return lines


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


def debug_print(s):
    from genomicode import parselib
    
    if not DEBUG:
        return
    parselib.print_split(s)


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



try:
    import cbie3
except ImportError:
    pass
else:
    import sys
    this_module = sys.modules[__name__]
    for name in cbie3.__dict__.keys():
        if name.startswith("__"):
            continue
        this_module.__dict__[name] = cbie3.__dict__[name]
