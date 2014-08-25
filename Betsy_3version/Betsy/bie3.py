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
Data
Constraint
Consequence
DefaultAttributesFrom
Module
ModuleDbSummary
Network

Functions:
make_network
backchain
complete_network
optimize_network
select_start_node
remove_data_node

get_inputs
group_inputs_by_datatype

summarize_moduledb
check_moduledb

print_modules
print_network
plot_network_gv

read_network
write_network

diagnose_start_node

"""
# Functions:
# _backchain_to_modules     Given an output Data, find compatible Modules.
# _backchain_to_input       DO NOT CALL.
# _backchain_to_all_inputs  Given an output Data and Module, make all in Datas.
# _forwardchain_to_outputs  Given Module and Datas, make output Data.
# 
# _backchain_to_ids         Given a Network and target ID, give source IDs.
# _make_backchain_dict
# _make_ancestor_dict
# 
# _can_module_take_data
# _can_module_take_data_network    Uses information from Network.
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
# debug_print

# _object_to_dict           Use for writing and reading json file
# _dict_to_object           Use for writing and reading json file



# TODO:
# - Some input data objects flow through the module and becomes the
#   output.  Need a name for this.


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
    TYPE_ATOM : "TYPE_ATOM",
    TYPE_ENUM : "TYPE_ENUM",
    
    MUST_BE : "MUST_BE",
    CAN_BE_ANY_OF : "CAN_BE_ANY_OF",
    SAME_AS : "SAME_AS",
    
    SET_TO : "SET_TO",
    SET_TO_ONE_OF : "SET_TO_ONE_OF",
    BASED_ON_DATA : "BASED_ON_DATA",
    SAME_AS_CONSTRAINT : "SAME_AS_CONSTRAINT",
    }


DEBUG = False
#DEBUG = True

DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES = False
#DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES = True

MAX_NETWORK_SIZE = 1024*8


class AttributeDef:
    def __init__(self, name, values, default_in, default_out, help=None):
        # Make sure name and values are valid.
        assert type(name) is type("")
        assert type(values) is type([])
        for x in values:
            assert type(x) is type("")

        # Make sure no duplicated values.
        seen = {}
        for x in values:
            assert x not in seen, "Duplicated value (%s) in %s." % (
                x, name)
            seen[x] = 1
            
        # Make sure default_in and default_out are valid values.
        assert type(default_in) is type("")
        assert type(default_out) is type("")
        assert default_in in values
        assert default_out in values
        
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
        x1 = [
            self.name, self.values, self.default_in, self.default_out,
            self.help]
        x2 = [
            other.name, other.values, other.default_in, other.default_out,
            self.help]
        return cmp(x1, x2)
    def __hash__(self):
        x = self.name, tuple(self.values), self.default_in, self.default_out, \
            self.help
        return hash(x)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            repr(self.name),
            repr(self.values),
            repr(self.default_in),
            repr(self.default_out),
            ]
        if self.help is not None:
            x.append("help=%r" % self.help)
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x
    @staticmethod
    def __init_from_dict(args):
        inst = AttributeDef(**args)
        return inst


class Attribute:
    def __init__(self, datatype, name, value):
        assert isinstance(datatype, DataType)
        assert type(name) is type("")
        assert type(value) is type("")
        
        # Check if this is a valid attribute name for the datatype.
        x = [x for x in datatype.attribute_defs if x.name == name]
        assert len(x) == 1, "datatype %r does not have attribute %r." % (
            datatype.name, name)
        attr = x[0]
        assert attr.is_valid_value(value), \
               "Invalid value %r for attribute %r." % (value, name)
        
        self.datatype = datatype
        self.name = name
        self.value = value
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            self.datatype.name,
            repr(self.name),
            repr(self.value),
            ]
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
        x = [
            repr(self.name),
            ]
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
        assert isinstance(module, Module)
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
                "input_index must be given for SAME_AS constraint")
            assert type(input_index) is type(0)
            if input_index is not None:
                assert arg1 != input_index
        else:
            raise AssertionError, "Invalid behavior (%s) for constraint %s." %(
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
        x = [
            repr(self.name),
            CONST2STR[self.behavior],
            ]
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
        inst = Constraint(args['name'],args['behavior'],
                          arg1=args['arg1'],input_index=args['input_index'])
        return inst


class Consequence:
    def __init__(self, name, behavior, arg1=None, arg2=None,
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
        x = [
            repr(self.name),
            CONST2STR[self.behavior],
            ]
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
        inst = Consequence(args['name'],args['behavior'],
                          arg1=args['arg1'],arg2=args['arg2'],
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
        x = [
            str(self.input_index),
            ]
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
            assert isinstance(x, AttributeDef)
        for x in keywds:
            assert x in ["help"]

        ## # Make sure no overlap between the attributes and user inputs.
        ## attr_names = [x.name for x in attributes]
        ## user_names = [x.name for x in user_inputs]
        ## for x in attr_names:
        ##     assert x not in user_names, "%s overlaps in DataType %s." % (
        ##         x, name)
            
        self.name = name
        self.attribute_defs = attribute_defs   # AttributeDef
        self.help = keywds.get("help")
    def get_attribute_def(self, name):
        x = [x for x in self.attribute_defs if x.name == name]
        assert len(x) > 0, "DataType %s has no attribute %s." % (
            repr(self.name), repr(name))
        assert len(x) == 1, "Multiple attributes with same name?"
        return x[0]
    def get_attribute_names(self):
        return [x.name for x in self.attribute_defs]
    def is_valid_attribute_name(self, name):
        return name in self.get_attribute_names()
    def is_valid_attribute_value(self, name, value):
        attr = self.get_attribute_def(name)
        return attr.is_valid_value(value)
    def __cmp__(self, other):
        if not isinstance(other, DataType):
            return cmp(id(self), id(other))
        # Bug: should compare attributes without regard to order.
        x1 = [self.name, self.attribute_defs, self.help]
        x2 = [other.name, other.attribute_defs, other.help]
        return cmp(x1, x2)
    def __hash__(self):
        x = self.name, tuple(self.attribute_defs), self.help
        return hash(x)
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
        for attr in self.attribute_defs:
            if attr.name in attrdict:
                continue
            value = attr.default_in
            if not is_input:
                value = attr.default_out
            attrdict[attr.name] = value
        return attrdict
    def input(self, **attribute_dict):
        # Create a Data object.
        attrdict = self._resolve_attributes(attribute_dict, True)
        # Don't bother checking the attributes here.  The Data object
        # will do that.
        return Data(self, **attrdict)
    def output(self, **attribute_dict):
        attrdict = self._resolve_attributes(attribute_dict, False)
        return Data(self, **attrdict)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [self.name]
        x += [repr(x) for x in self.attribute_defs]
        if self.help:
            x.append("help=%r" % self.help)
        return "DataType(%s)" % ", ".join(x)
    @staticmethod
    def __init_from_dict(args):
        assert 'name' in args
        assert 'attribute_defs' in args
        assert 'help' in args
        inst = DataType(
            args['name'], *args['attribute_defs'], help=args['help'])
        return inst


class Data(object):
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
        for name, value in keywds.iteritems():
            assert datatype.is_valid_attribute_name(name), \
                   "'%s' is not a known attribute for datatype %s." % (
                name, datatype.name)
            assert datatype.is_valid_attribute_value(name, value), \
                   "In a %s, '%s' is not a valid value for '%s'." % (
                datatype.name, value, name)

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
        if not isinstance(other, Data):
            return cmp(id(self), id(other))
        x1 = [self.datatype, self.attributes]
        x2 = [other.datatype, other.attributes]
        return cmp(x1, x2)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        #keywds = self.attributes.copy()
        x = [
            self.datatype.name,
            #_pretty_attributes(keywds),
            _pretty_attributes(self.attributes),
            ]
        x = [x for x in x if x]
        return "Data(%s)" % ", ".join(x)
    @staticmethod
    def __init_from_dict(args):
        assert 'datatype' in args
        assert 'attributes' in args
        inst = Data(args['datatype'],**args['attributes'])
        return inst


class Module:
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
        option_defs = []   # OptionDef
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
        if len(in_datatypes) == 1 and in_datatypes[0] == out_datatype and \
           not default_attributes_from:
            default_attributes_from = [DefaultAttributesFrom(0)]

        # Check default_attributes_from.  Should be a list of
        # DefaultAttributesFrom objects.
        assert len(default_attributes_from) <= len(in_datatypes)
        seen = {}  # make sure no duplicates
        for daf in default_attributes_from:
            assert daf.input_index < len(in_datatypes)
            assert out_datatype == in_datatypes[daf.input_index]
            assert daf.input_index not in seen
            seen[daf.input_index] = 1

        # default_attributes_from can be an empty list if the
        # attributes of the output object should come from none of the
        # input data types.  Probably the most common case if if the
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
            self._assert_constraint(
                name, in_datatypes, out_datatype, constraints, consequences, x)
        for x in consequences:
            self._assert_consequence(
                name, in_datatypes, out_datatype, constraints, x)
            
    def _assert_constraint(
        self, name, in_datatypes, out_datatype, constraints, consequences,
        constraint):
        # Get the input datatype that this constraint refers to.
        i = constraint.input_index
        assert i < len(in_datatypes)
        in_datatype = in_datatypes[i]

        assert constraint.behavior in [MUST_BE, CAN_BE_ANY_OF, SAME_AS]
        if constraint.behavior in [MUST_BE, CAN_BE_ANY_OF]:
            assert in_datatype.is_valid_attribute_value(
                constraint.name, constraint.arg1), (
                "%r: Invalid value %r for attribute %r." % (
                    name, constraint.arg1, constraint.name))
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
        if in_datatype == out_datatype:
            x = [x for x in consequences if x.name == constraint.name]
            assert x, "%r: constraint but no consequence for %r." % (
                name, constraint.name)

    def _assert_consequence(
        self, name, in_datatypes, out_datatype, constraints, consequence):
        import itertools
        assert consequence.name in out_datatype.get_attribute_names(), \
               "Module %r refers to an unknown attribute %r." % (
            self.name, consequence.name)
        
        if consequence.behavior in [SET_TO, SET_TO_ONE_OF, BASED_ON_DATA]:
            assert out_datatype.is_valid_attribute_value(
                consequence.name, consequence.arg1), \
                "'%s' is not a valid value for '%s' in module '%s'" % (
                consequence.arg1, consequence.name, name)
        elif consequence.behavior == SAME_AS_CONSTRAINT:
            # Make sure index on input variable is reasonable.
            index = consequence.arg1   # index of input variable
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
                if in_datatype != out_datatype:
                    continue
                if consequence.name != constraint.name:
                    continue
                if constraint.behavior != MUST_BE:
                    continue
                assert constraint.arg1 != consequence.arg1, \
                       "MUST_BE and SET_TO should not be set to the " \
                       "same value in module %s (%s)." % self.name

    def __cmp__(self, other):
        if not isinstance(other, Module):
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
        params = (args['consequences']+args['constraints']+
                    args['option_defs']+args['default_attributes_from'])
        inst = Module(name, in_datatypes, out_datatype, *params, help=help_)
        return inst


class ModuleDbSummary:
    # module_names    List of module names (strings).
    # name2module     Dictionary of module name (string) to Module object.
    # name2datatypes  Dict of module name to (in Datatype, out Datatype).
    # datatypes       List of Datatype objects.
    def __init__(self, module_names, name2module, name2datatypes, datatypes):
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

    def delete_nodes(self, node_ids):
        # Make sure no duplicates.
        for i in range(len(node_ids)):
            assert node_ids[i] not in node_ids[i+1:]
        network = Network(self.nodes, self.transitions)
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
        for i in range(len(node_ids)-1):
            assert node_ids[i] != node_ids[i+1], "Duplicate node IDs."

        # Make sure nodes are the same type.
        for i in range(1, len(node_ids)):
            n1 = self.nodes[node_ids[0]]
            n2 = self.nodes[node_ids[i]]
            assert n1.__class__ == n2.__class__

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

        merged = Network(self.nodes, transitions)
        merged = merged.delete_nodes(node_ids[1:])
        return merged
    
    def __cmp__(self, other):
        if not isinstance(other, Network):
            return cmp(id(self), id(other))
        x1 = [self.nodes, self.transitions]
        x2 = [other.nodes, other.transitions]
        return cmp(x1, x2)


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
            if attr.datatype != out_data:
                continue
            attrdict[attr.name] = attr.value
        out_data = out_data.output(**attrdict)
    assert isinstance(out_data, Data)


    nodes = []        # list of Data or Module objects.
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
            all_inputs = _backchain_to_all_inputs(node, cons, user_attributes)
            for d in all_inputs:
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
    return network


def _make_ancestor_dict(network):
    # Return a dictionary of node_id -> all node_ids that are
    # ancestors of node_id.
    if not network.nodes:
        return {}
    
    node2parents = _make_backchain_dict(network)

    ancestors = {}  # node id -> list of parent node ids.
    all_nodes = list(range(len(network.nodes)))
    while all_nodes:
        node_id = all_nodes.pop(0)
        parents = node2parents.get(node_id, [])

        # If there are no parents, then it has no ancestors.
        if not parents:
            ancestors[node_id] = []
            continue

        # If I haven't found the ancestors of all my parents, try
        # again later.
        # BUG: will have infinite loop if there are cycles.
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

    # For each Data object, check to see if it can be the antecedent
    # of any Module objects.
    data_ids = [x for x in range(len(network.nodes))
                if isinstance(network.nodes[x], Data)]
    module_ids = [x for x in range(len(network.nodes))
                  if isinstance(network.nodes[x], Module)]
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
            #debug_print("Skipping Data %d -> Module %d (cycle)." % (
            #    input_id, module_id))
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
            debug_print("Completing Data %d -> Module %d." % (id_, module_id))
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
        _OptimizeNoInvalidOutputs(),  # Fixes NoOverlappingData.
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
            node_ids = [
                i for i in range(len(network.nodes)) if i not in noncycle]
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
                network, start_id, max_path_length,
                noncycle, nodeid2previds)
            if cycle:
                break
        return cycle

    def _find_cycle_from_one_node(
        self, network, start_id, max_path_length, noncycle, nodeid2previds):
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
        for i in range(len(cycle)-1):
            next_id, node_id = cycle[i], cycle[i+1]
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
        x = x[:i] + x[i+1:]
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
            # Will merge high node_id into low node_id, so id of root
            # node will never be changed.
            network = network.merge_nodes(duplicates)
        return network
    
    def find_duplicate_data(self, network):
        # Return list of node_ids for Data objects that are
        # duplicated.  If no duplicates found, return an empty list.
        import itertools

        # Make a list of all pairs of Data objects.
        data_nodes = [
            node_id for (node_id, node) in enumerate(network.nodes)
            if isinstance(node, Data)]
        data_pairs = itertools.product(data_nodes, data_nodes)

        # Can optimize this by sorting.
        for (node_id1, node_id2) in data_pairs:
            if node_id1 == node_id2:
                continue
            data1, data2 = network.nodes[node_id1], network.nodes[node_id2]
            # Everything has to be the same.
            if data1 == data2:
                return [node_id1, node_id2]
        return []


class _OptimizeNoDuplicateModules:
    def __init__(self):
        pass
    
    def optimize(self, network, user_attributes):
        while 1:
            duplicates = self.find_duplicate_modules(network)
            if not duplicates:
                break
            network = network.merge_nodes(duplicates)
            assert isinstance(network, Network)
        return network
    
    def find_duplicate_modules(self, network):
        # Return list of module_ids for modules that are duplicated.  If
        # no duplicates found, return an empty list.

        # DEFINITION: If the same data node points to two of the same
        # module nodes, then those modules are duplicated.
        for data_id, next_ids in network.iterate(node_class=Data):
            if len(next_ids) < 2:
                continue

            for i in range(len(next_ids)-1):
                node_id1 = next_ids[i]
                node_1 = network.nodes[node_id1]
                dups = [node_id1]
                for j in range(i+1, len(next_ids)):
                    node_id2 = next_ids[j]
                    node_2 = network.nodes[node_id2]
                    if node_1 == node_2:
                        dups.append(node_id2)
                if len(dups) > 1:
                    return dups
        return None


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
        dangling = []
        for node_id in range(len(network.nodes)):
            if self.is_dangling_node(network, node_id, user_attributes):
                dangling.append(node_id)
        return dangling
    def is_dangling_node(self, network, node_id, user_attributes):
        # 1.  Module nodes with no valid input_combinations.
        # 2.  Module nodes with no outputs.
        # 3.  Data nodes (except for node 0) that don't point to any
        #     Modules.

        # Case 1.
        node = network.nodes[node_id]
        if isinstance(node, Module):
            prev_ids = _backchain_to_ids(network, node_id)
            x = _get_valid_input_combinations(
                network, node_id, prev_ids, user_attributes)
            if not x:
                return True
        # Case 2.
        if isinstance(node, Module):
            if not network.transitions.get(node_id, []):
                return True
        # Case 3.
        if isinstance(node, Data):
            if node_id != 0 and not network.transitions.get(node_id, []):
                return True
            
        return False


class _OptimizeNoOverlappingData:
    # Some data objects may have overlapping attributes.
    # 1.  Data(SignalFile, format=['tdf', 'pcl', 'gct'], preprocess='rma')
    # 2.  Data(SignalFile, format=['tdf', 'res', 'pcl', 'jeffs', 'unknown',
    #       'xls'], preprocess='rma')
    #
    # Modules that point to 1 can generate 3 different formats.  Those
    # that point to 2 can generate 6.  To resolve this, split the Data
    # objects and rewire the Modules.
    # 3.  Data(SignalFile, format=['tdf', 'pcl'], preprocess='rma')
    # 4.  Data(SignalFile, format='gct', preprocess='rma')
    # 5.  Data(SignalFile, format=['res', 'jeffs', 'unknown', 'xls'],
    #       preprocess='rma')

    # Methods:
    # optimize
    # 
    # _find_overlapping_data
    # _find_overlapping_attribute
    # 
    # _fix_overlapping_data
    # _remove_atom_from_list
    # _remove_list_from_list
    # _split_list
    
    def __init__(self):
        pass

    def optimize(self, network, user_attributes):
        # Look for pairs of Data objects (Data_1, Data_2) where an
        # attribute is overlapping.
        import copy

        # Make a copy of the network now.  All changes to the network
        # will be done in place.
        network = copy.deepcopy(network)
        while True:
            x = self._find_overlapping_data(network)
            if not x:
                break
            node_id1, node_id2, attr_name = x
            self._fix_overlapping_data(
                network, node_id1, node_id2, attr_name, user_attributes)
            
        return network

    def _find_overlapping_data(self, network):
        # Return (node_id1, node_id2, name of overlapping attribute)
        # or None.
        for (node_id1, node_id2) in _iter_upper_diag(len(network.nodes)):
            node_1 = network.nodes[node_id1]
            node_2 = network.nodes[node_id2]
            if not isinstance(node_1, Data):
                continue
            if not isinstance(node_2, Data):
                continue
            attr = self._find_overlapping_attribute(
                network, node_id1, node_id2)
            if attr:
                return node_id1, node_id2, attr
        return None

    def _find_overlapping_attribute(self, network, data_id1, data_id2):
        # Return the name of the single overlapping attribute or None.
        # data1 and data2 should be exactly the same, except for one
        # attribute with overlapping values.
        data1 = network.nodes[data_id1]
        data2 = network.nodes[data_id2]
        assert isinstance(data1, Data)
        assert isinstance(data2, Data)
        if data1.datatype != data2.datatype:
            return False
        
        # CASE    DATA1      DATA2     RESULT
        #   1      ATOM       ATOM     OK if ATOMs are equal.
        #   2      ATOM       ENUM     OVERLAP if ATOM in ENUM; DATA2 not root.
        #   3      ENUM       ATOM     OVERLAP if ATOM in ENUM; DATA1 not root.
        #   4      ENUM       ENUM     OVERLAP if ENUMs share ATOMs; not root.
        data1_attr = data1.attributes
        data2_attr = data2.attributes

        mismatch = False
        overlapping = []   # list of attribute names

        x = data1_attr.keys() + data2_attr.keys()
        all_attributes = sorted({}.fromkeys(x))
        for key in all_attributes:
            assert key in data1_attr
            assert key in data2_attr
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

        if mismatch:
            return None
        if len(overlapping) != 1:
            return None
        return overlapping[0]

    def _fix_overlapping_data(
        self, network, data_id1, data_id2, attr_name, user_attributes):
        # Changes network in place.

        data1 = network.nodes[data_id1]
        data2 = network.nodes[data_id2]
        DATA1_VALUE = data1.attributes[attr_name]
        DATA2_VALUE = data2.attributes[attr_name]
        DATA1_TYPE = _get_attribute_type(DATA1_VALUE)
        DATA2_TYPE = _get_attribute_type(DATA2_VALUE)
        if DATA1_TYPE is TYPE_ENUM and DATA2_TYPE is TYPE_ENUM:
            COMMON_VALUES = _intersection(DATA1_VALUE, DATA2_VALUE)

        if DATA1_TYPE is TYPE_ATOM:
            assert DATA2_TYPE is TYPE_ENUM
            self._remove_atom_from_list(network, data_id2, data_id1, attr_name)
        elif DATA2_TYPE is TYPE_ATOM:
            assert DATA1_TYPE is TYPE_ENUM
            self._remove_atom_from_list(network, data_id1, data_id2, attr_name)
        elif sorted(DATA1_VALUE) == sorted(COMMON_VALUES):
            self._remove_list_from_list(network, data_id2, data_id1, attr_name)
        elif sorted(DATA2_VALUE) == sorted(COMMON_VALUES):
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

        # Make a new Data object that contains the common values.
        attributes = data1.attributes.copy()
        attributes[attr_name] = COMMON_VALUE
        data3 = Data(data1.datatype, **attributes)
        network.nodes.append(data3)
        data_id3 = len(network.nodes)-1
            
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
        module_ids = sorted({}.fromkeys(x1+x2))
        for module_id in module_ids:
            if _can_module_in_network_take_data(
                network, module_id, [data3], user_attributes):
                if data_id3 not in network.transitions:
                    network.transitions[data_id3] = []
                network.transitions[data_id3].append(module_id)
            # Bug here.  Need to check output and user attributes.
            #module = network.nodes[module_id]
            #if _can_module_take_data(module, [data3]):
            #    if data_id3 not in network.transitions:
            #        network.transitions[data_id3] = []
            #    network.transitions[data_id3].append(module_id)
        

class _OptimizeNoInvalidOutputs:
    # Fixing overlapping data can lead to a situation where a Module
    # points to Data that it can't generate.  E.g.
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
        for (node_id, next_ids) in network.iterate(node_class=Module):
            module = network.nodes[node_id]
            for next_id in next_ids:
                node = network.nodes[next_id]
                assert isinstance(node, Data)
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


def select_start_node(network, start_data):
    # start_data may be a single Data object or a list of Data
    # objects.  DataTypes are also allowed in lieu of Data objects.

    raise NotImplementedError, "Not finished yet."


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
            if not isinstance(node, Module):
                continue
            assert len(node.in_datatypes) > 0
            if len(node.in_datatypes) <= 1:
                continue
            in_ids = _backchain_to_ids(network, node_id)
            nid = [x for x in in_ids if x in good_ids]
            if len(nid) != len(in_ids):
                delete_ids[node_id] = 1

        for node_id in delete_ids:
            del good_ids[node_id]
        if not delete_ids:
            break
        
    # Delete all the IDs that aren't in good_ids.
    bad_ids = [x for x in range(len(network.nodes)) if x not in good_ids]
    network = network.delete_nodes(bad_ids)
    return network


def remove_data_node(network, *attributes):
    # Look for the nodes that are compatible with the attributes.  If
    # an attribute has multiple values, Data that matches any of those
    # values will be removed.
    node_ids = []  # list of node_ids to remove
    for node_id, next_ids in network.iterate(node_class=Data):
        node = network.nodes[node_id]
        if not isinstance(node, Data):
            continue

        # If compatible, remove or alter this node.
        remove_attr = {}  # name -> list of values
        for attr in attributes:
            assert isinstance(attr, Attribute)
            if node.datatype != attr.datatype:
                break
            assert attr.name in node.attributes
            assert attr.name not in remove_attr

            DATA_VALUE = node.attributes[attr.name]
            ATTR_VALUE = attr.value
            DATA_TYPE = _get_attribute_type(DATA_VALUE)
            ATTR_TYPE = _get_attribute_type(ATTR_VALUE)
            case = _assign_case_by_type(DATA_TYPE, ATTR_TYPE)

            if case == 1:
                if DATA_VALUE != ATTR_VALUE:
                    break
                remove_attr[attr.name] = [DATA_VALUE]
            elif case == 2:
                # DATA is ATOM, ATTR is ENUM.
                if DATA_VALUE not in ATTR_VALUE:
                    break
                remove_attr[attr.name] = [DATA_VALUE]
            elif case == 3:
                # DATA is ENUM, ATTR is ATOM
                if ATTR_VALUE not in DATA_VALUE:
                    break
                remove_attr[attr.name] = [ATTR_VALUE]
            else:
                # DATA is ENUM, ATTR is ATOM
                x = _intersection(DATA_VALUE, ATTR_VALUE)
                if not x:
                    break
                remove_attr[attr.name] = x
        else:
            # Remove or alter this node.
            for name in remove_attr:
                assert name in node.attributes
                values = node.attributes[name][:]
                if type(values) is type(""):
                    values = [values]
                for value in remove_attr[name]:
                    i = values.index(value)
                    values.pop(i)
                if not values:
                    # Removed everything, delete this node.
                    node_ids.append(node_id)
                    break
                if len(values) == 1:
                    values = values[0]
                node.attributes[name] = values

    # Remove each of these node_ids from the network.
    network = network.delete_nodes(node_ids)
    network = _OptimizeNoDanglingNodes().optimize(network)
    return network


def _get_inputs_h(network, node_id, nodeid2previds, user_attributes):
    import itertools
    node = network.nodes[node_id]

    inputs = []
    if isinstance(node, Data):
        # I am a valid input, as well as all inputs generated from any
        # of the previous nodes.
        inputs = [(node_id,)]
        for previd in nodeid2previds.get(node_id, []):
            x = _get_inputs_h(network, previd, nodeid2previds, user_attributes)
            inputs.extend(x)
    elif isinstance(node, Module):
        # Find all potential sets of Data notes that can feed into me.
        assert node_id in nodeid2previds
        previds = nodeid2previds[node_id]
        combos = _get_valid_input_combinations(
            network, node_id, previds, user_attributes, 
            node2previds=nodeid2previds)

        # Get the inputs from each of the combinations.
        for combo in combos:
            # Get the inputs for each branch of this combination.
            branch2inputs = [
                _get_inputs_h(network, x, nodeid2previds, user_attributes)
                for x in combo]

            # Find all permutations for the inputs for each branch.
            # [ [(49, 40)]              Each branch has a set of inputs.
            #   [(38,), (35,)] ]
            # ->                        itertools.product
            # [ ((49, 40), (38,)),
            #   ((49, 40), (35,)) ]
            # ->                        flatten
            # [ (49, 40, 38), (49, 40, 35) ]
            x = itertools.product(*branch2inputs)
            # Flatten the tuples.
            x = [_flatten(x) for x in x]
            inputs = x
    else:
        raise AssertionError

    inputs = sorted({}.fromkeys(inputs))
    return inputs
    

def get_inputs(network, user_attributes):
    # Return a list of tuples of node ids.  Each tuple contains a set
    # of node IDs that can serve as the inputs to this network.
    # 
    # Example return value:
    #   [(1, 5), (8,)]
    # This means that the set of nodes 1 and 5 would make a valid input
    # to this network.  Or, node 8 by itself would also be a valid
    # input.
    nodeid2previds = _make_backchain_dict(network)
    inputs = _get_inputs_h(network, 0, nodeid2previds, user_attributes)
    return inputs


def group_inputs_by_datatype(network, inputs):
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

    dt2inputs = {}
    for inp in inputs:
        nodes = [network.nodes[x] for x in inp]
        datatypes = [x.datatype for x in nodes]
        y = sorted(zip(datatypes, inp))
        datatypes = tuple([x[0] for x in y])
        inp = tuple([x[1] for x in y])

        if datatypes not in dt2inputs:
            dt2inputs[datatypes] = []
        dt2inputs[datatypes].append(inp)
    return dt2inputs
    
    
def summarize_moduledb(moduledb):
    """Take a list of Modules and return a ModuleDbSummary object."""
    name2module = {}   # module_name -> Module
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
        assert isinstance(module, Module)
        seen[module.name] = seen.get(module.name, 0) + 1

    # Make sure no duplicate modules.
    dups = ["%s (%d times)" % (x, seen[x]) for x in seen if seen[x] > 1]
    assert not dups, "Duplicate modules: %s" % "\n".join(dups)
        

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
    line_width = 72
    for i in range(len(network.nodes)):
        p_step = "%2d.  " % i
        p_space = " " * (len(p_step)+2)
        _print_line(
            str(network.nodes[i]), p_step, p_space, line_width,
            outhandle=outhandle)
    print >>outhandle

    for i in sorted(network.transitions):
        x = [i, "->"] + network.transitions[i]
        print >>outhandle, "\t".join(map(str, x))
    

def plot_network_gv(filename, network):
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
        elif n.__class__.__name__ == "Module":
            x = n.name
            node2attr["shape"] = "box"
            node2attr["fillcolor"] = "#80E0AA"
        elif n.__class__.__name__ == "QueryModule":
            x = n.name
            node2attr["shape"] = "box"
            node2attr["fillcolor"] = "#60A08A"
        else:
            raise AssertionError
            
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


def read_network(file_or_handle):
    import json
    handle = file_or_handle
    if type(handle) is type(""):
        handle = open(file_or_handle, 'r')
    text =handle.read()
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


def diagnose_start_node(network, user_data):
    results = _score_start_nodes(network, user_data)

    header = ["Score", "Node ID", "Attribute", "Network", "User"]
    print "\t".join(header)
    for x in results:
        score, node_id, user_data, attr_values = x
        for name, netw_value, user_value in attr_values:
            x = score, node_id, name, netw_value, user_value
            assert len(x) == len(header)
            print "\t".join(map(str, x))

    
def _backchain_to_modules(moduledb, data):
    # Return list of modules that can generate an output that is
    # compatible with data.

    modules = []  # list of (module, num compatible attributes)
    for module in moduledb:
        if _can_module_produce_data(module, data):
            modules.append(module)
    return modules


def _backchain_to_input(module, in_num, out_data, user_attributes):
    # Given a module and output_data, return the input_data object
    # that can generate the output.  This should only be called by
    # _backchain_to_all_inputs.  The SAME_AS constraint is handled
    # there.  If this is called by itself, the attributes might not be
    # right.
    assert in_num < len(module.in_datatypes)

    in_datatype = module.in_datatypes[in_num]
    out_datatype = out_data.datatype
    
    # Can't generate debug messages here because the SAME_AS
    # constraints aren't handled in this function.

    # The attributes for the input object should come from (in
    # decreasing priority):
    # 1.  Consequence (i.e. SAME_AS_CONSTRAINT).
    # 2.  Constraint.
    # 3.  user attribute
    # 4.  out_data
    # 5.  default output value of input datatype
    #
    # If the user attribute conflicts with a Constraint, thte
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
    for attrdef in in_datatype.attribute_defs:
        default = attrdef.default_out
        if DEFAULT_INPUT_ATTRIBUTE_IS_ALL_VALUES:
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
            if attr.datatype != in_datatype:
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
            ##        "Module %s requires a %s with %s=%s, but user requests "
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
            else:
                raise AssertionError
        else:
            raise AssertionError
                
        attributes[n] = data_value
        attrsource[n] = "consequence"

    return attributes, attrsource


def _backchain_to_all_inputs(module, out_data, user_attributes):
    all_attributes = []
    all_attrsource = []
    for in_num in range(len(module.in_datatypes)):
        x = _backchain_to_input(module, in_num, out_data, user_attributes)
        attributes, attrsource = x
        all_attributes.append(attributes)
        all_attrsource.append(attrsource)

    # Handle the SAME_AS constraints here.  Difficult to do it in
    # _backchain_to_input because not all inputs have been created
    # yet.
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
        assert name in attr_src
        assert name in attr_dst
        attr_dst[name] = attr_src[name]
        all_attrsource[i_dst][name] = "SAME_AS,%d" % (i_src)

    # Make the input objects.
    all_inputs = []
    for in_num in range(len(module.in_datatypes)):
        in_datatype = module.in_datatypes[in_num]
        out_datatype = out_data.datatype
        attributes = all_attributes[in_num]
        attrsource = all_attrsource[in_num]
        
        debug_print("Backchaining %s (input=%d) -> %s -> %s." % (
            in_datatype.name, in_num, module.name, out_datatype.name))
        #debug_print("Generating a %s with attributes:" % in_datatype.name)
        for name in sorted(attributes):
            debug_print("  %s=%s (%s)" % (
                name, attributes[name], attrsource[name]))
        x = in_datatype.output(**attributes)
        all_inputs.append(x)

    return all_inputs


def _forwardchain_to_outputs(module, in_datas):
    # Generate a list of Data objects that can be generated from
    # module and in_datas.  Multiple objects can be generated because
    # the consequences can vary.  E.g. center_genes can set
    # gene_center to either "mean" or "median".  It can be either, but
    # must be one of them.
    import itertools
    
    # Check the input variables.
    assert len(module.in_datatypes) == len(in_datas)
    for i in range(len(module.in_datatypes)):
        assert in_datas[i].datatype == module.in_datatypes[i]

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
        assert data.datatype == datatype
        attributes.update(data.attributes)
    # Case 2.
    else:
        for attrdef in datatype.attribute_defs:
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
        x = Data.__new__(Data)
        x.datatype = datatype
        x.attributes = attributes.copy()
        return [x]

    names = sorted(possibilities)
    args = [possibilities[x] for x in names]
    outputs = []
    for values in itertools.product(*args):
        for key, value in zip(names, values):
            attributes[key] = value
            # Optimization: Data.__init__ is very expensive because of
            # all the checks.  Skip the checks and instantiate the
            # class directly.
            #x = Data(datatype, **attributes)
            x = Data.__new__(Data)
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
    return sorted(ids)


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
##                 # data_value contains the possible values of this Data
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
##     # Return True/False if a module can take this list of Data nodes
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


def _can_module_take_one_data(
    module, input_num, in_data, out_data, user_attributes):
    assert input_num < len(module.in_datatypes)

    # BUG: should not call _backchain_to_input!!!
    #good_in_data = _backchain_to_input(
    #    module, input_num, out_data, user_attributes)
    all_inputs = _backchain_to_all_inputs(module, out_data, user_attributes)
    return _is_data_compatible_with(in_data, all_inputs[input_num])


def _can_module_take_data(module, in_datas, out_data, user_attributes):
    # Return True/False if a module can take in_datas (list of Data
    # nodes) as input.

    assert len(in_datas) == len(module.in_datatypes)
    all_inputs = _backchain_to_all_inputs(module, out_data, user_attributes)
    for i in range(len(in_datas)):
        if not _is_data_compatible_with(in_datas[i], all_inputs[i]):
            return False
    return True


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
    if module.out_datatype != data.datatype:
        #debug_print(
        #    "Module can't generate data type: %s." % data.datatype.name)
        return False

    debug_print("Testing if module %s can produce data %s." % (
        repr(module.name), str(data)))
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
                       "but data contains '%s'." % (
                           consequence.name, outc_value, data_value))
                debug_print(msg)
                return False
        elif case == 2:
            # Module can produce any of a list of values.  Check if
            # the data's value can be produced by the module.
            if data_value not in outc_value:
                debug_print("Consequence %s conflicts." % consequence.name)
                return False
        elif case == 3:
            # Module produces a specific value.  Data could be one of
            # many values.
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
        debug_print("Module converts datatype.  Checking default attributes.")
        consequence_names = [x.name for x in module.consequences]
        for attrdef in module.out_datatype.attribute_defs:
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
                    debug_print("Attr %r: Conflicts (module %r, data %r)." % (
                        attrdef.name, attrdef.default_in, data_value))
                    return False
            elif data_type == TYPE_ENUM:
                if attrdef.default_in not in data_value:
                    debug_print("Attr %r: Conflicts (module %r, data %r)." % (
                        attrdef.name, attrdef.default_in, data_value))
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
        x = [x for x in module.constraints if
             x.name == consequence.name and 
             x.input_index == consequence.arg1]
        assert len(x) == 1
        const1 = x[0]

        # Find the constraint that refers to the same object.
        x = [x for x in module.constraints if
             x.name == consequence.name and x.input_index in indexes]
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

    module = network.nodes[module_id]

    # For each in_datatype, find all the data objects that match this
    # type.
    args = []
    for datatype in module.in_datatypes:
        x = [x for x in all_input_ids if network.nodes[x].datatype == datatype]
        args.append(x)

    # If there are many in_datatypes, then this could cause a
    # combinatorial explosion of possibilities.  Optimize by throwing
    # out data objects that fail a constraint.
    for i in range(len(args)):
        x = args[i]
        x = [x for x in x
             if _can_module_in_network_take_one_data(
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
        x = [(x,) for x in ids]
        valid.extend(x)
    # The multiple in_datatypes case is harder, because we don't know
    # the order of the inputs.  
    
    for input_ids in itertools.product(*args):
        assert len(input_ids) == len(module.in_datatypes)

        # Don't check again if already done.
        if input_ids in valid:
            continue
        
        # No duplicated IDs.
        x = {}.fromkeys(input_ids)
        if len(x) != len(input_ids):
            continue

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
    assert isinstance(node, Data)

    # All values need to be exactly equal.
    # 
    # CASE   N1_TYPE    N2_TYPE    RESULT
    #   1      ATOM       ATOM     OK if ATOM equal.
    #   2      ATOM       ENUM     No.
    #   3      ENUM       ATOM     No.
    #   4      ENUM       ENUM     OK if ENUM equal.
    
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
    
    # Make a list of all the desired start nodes, as Data objects.
    user_datas = user_data
    if not operator.isSequenceType(user_data):
        user_datas = [user_data]
    for i, x in enumerate(user_datas):
        if isinstance(x, DataType):
            x = x.input()  # convert to Data
        assert isinstance(x, Data)
        user_datas[i] = x

    # Look for the nodes in the network that are compatible with
    # user_data.
    results = []
    for node_id, next_ids in network.iterate(node_class=Data):
        for user_data in user_datas:
            netw_data = network.nodes[node_id]
            if netw_data.datatype != user_data.datatype:
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
    results.sort()

    return results


def _iter_upper_diag(n):
    for i in range(n-1):
        for j in range(i+1, n):
            yield i, j


def _score_one_start_node(network_data, user_data):
    # network_data is a Data node in the network.  user_data is the
    # Data that the user wants to start on.  Return a list of
    # <attributes> that are incompatible between the network_data and
    # user_data.  An empty list indicates that the two nodes are
    # compatible.
    assert network_data.datatype == user_data.datatype

    # Start Data
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
    if data_s.datatype != data_g.datatype:
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
    types = [
        (TYPE_ATOM, TYPE_ATOM),
        (TYPE_ATOM, TYPE_ENUM),
        (TYPE_ENUM, TYPE_ATOM),
        (TYPE_ENUM, TYPE_ENUM),
        ]
    x = (type1, type2)
    assert x in types, "Unknown types: %s %s" % (type1, type2)
    i = types.index(x)
    return i+1


def _get_attribute_type(name):
    if type(name) is type(""):
        return TYPE_ATOM
    elif type(name) in [type([]), type(())]:
        return TYPE_ENUM
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


def _print_line(line, prefix0, prefixn, width, outhandle=None):
    import sys
    
    outhandle = outhandle or sys.stdout

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
        print >>outhandle, "%s%s" % (p, lines[i])


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
    if DEBUG:
        print s


def _object_to_dict(obj):
    # Convert objects to a dictionary of their representation
    d = { '__class__':obj.__class__.__name__, 
          '__module__':obj.__module__,
          }
    d.update(obj.__dict__)
    return d


def _dict_to_object(d):
    assert isinstance(d,dict)
    args = dict((key.encode('ascii'), value) for key, value in d.items())
    for key,value in args.iteritems():
        if isinstance(value,dict):
            args[key]=_dict_to_object(value)
        elif isinstance(value,list):
            if value:
                if isinstance(value[0],unicode):
                    args[key]=[i.encode('ascii') for i in value]
        elif isinstance(value,unicode):
            args[key]=value.encode('ascii')
        else:
            assert 'not expected type %s' %value
    inst = args
    if '__class__' in args:
        class_name = args.pop('__class__')
        module_name = args.pop('__module__')
        if '.' in module_name:
            module = __import__(
                module_name, globals(), locals(), [module_name.split('.')[-1]],
                -1)
        else:
            module = __import__(module_name)
        class_ = getattr(module, class_name)
        fn = getattr(class_,'_' + class_name + '__init_from_dict')
        inst = fn(args)
    return inst
