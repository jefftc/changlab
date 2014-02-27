"""

Glossary:
data              A unit of information.  Typically a file made by a program.
attribute         Describes something about the data.  e.g. logged="yes".
user input        Used in the modules, but does not affect the inferencing.
data type         Describes a set of data that share the same attributes.
module            Takes one or more data objects as input and produces a single
                  data object as output.
                  Modules "produce" data.
                  Modules can "convert" data from one type to another.
input attribute   An attribute of an input data.
output attribute  An attribute of an output data.


Classes:
AttributeDef
Attribute
UserInputDef
UserInput
DataType
Data
Constraint
Consequence
DefaultAttributesFrom
Module
ModuleDbSummary
Network

Functions:
backchain
select_start_node
optimize_network

summarize_moduledb

print_modules
print_network
plot_network_gv
diagnose_start_node

"""
# Functions:
# _backchain_to_modules     Given an output, find compatible Modules.
# _backchain_to_input       Given an output and module, make the input.
# _backchain_to_all_inputs  Given an output and module, make all inputs.
# _backchain_to_ids
# _make_backchain_dict
# 
# _can_module_take_data
# _can_module_produce_data
# 
# _can_reach_by_bc
# _can_reach_by_fc
#
# _find_data_node
# _find_start_nodes
# _score_start_nodes
# _score_one_start_node
#
# _assign_case_by_type
# _get_attribute_type
# _intersection
# _is_subset
# 
# _print_nothing
# _print_string
# _print_line
# _pretty_attributes


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


#DEBUG = False
DEBUG = True


class AttributeDef:
    def __init__(self, name, values, default_in, default_out):
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
            self.name, self.values, self.default_in, self.default_out]
        x2 = [
            other.name, other.values, other.default_in, other.default_out]
        return cmp(x1, x2)
    def __hash__(self):
        x = self.name, tuple(self.values), self.default_in, self.default_out
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
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x


class Attribute:
    def __init__(self, datatype, name, value):
        assert isinstance(datatype, DataType)
        assert type(name) is type("")
        assert type(value) is type("")
        
        # Check if this is a valid attribute name for the datatype.
        x = [x for x in datatype.attributes if x.name == name]
        assert len(x) == 1, "datatype %r does not have attribute %r." % (
            datatype.name, name)
        attr = x[0]
        assert attr.is_valid_value(value), \
               "Invalid value %r for attribute %r." % (value, name)
        
        self.datatype = datatype
        self.name = name
        self.value = value


class UserInputDef:
    def __init__(self, name, default=None):
        assert type(name) is type("")
        self.name = name
        self.default = default
    def __cmp__(self, other):
        if not isinstance(other, UserInputDef):
            return cmp(id(self), id(other))
        x1 = [self.name, self.default]
        x2 = [other.name, other.default]
        return cmp(x1, x2)
    def __hash__(self):
        x = self.name, self.default
        return hash(x)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        x = [
            repr(self.name),
            ]
        if self.default is not None:
            x.append(repr(self.default))
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x


class UserInput:
    def __init__(self, module, name, value):
        assert type(module) is type(Module)
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


class DataType:
    def __init__(self, name, *attributes):
        for x in attributes:
            assert isinstance(x, AttributeDef)

        ## # Make sure no overlap between the attributes and user inputs.
        ## attr_names = [x.name for x in attributes]
        ## user_names = [x.name for x in user_inputs]
        ## for x in attr_names:
        ##     assert x not in user_names, "%s overlaps in DataType %s." % (
        ##         x, name)
            
        self.name = name
        self.attributes = attributes   # AttributeDef
        ## self.user_inputs = user_inputs
    def get_attribute(self, name):
        x = [x for x in self.attributes if x.name == name]
        assert len(x) > 0, "DataType %s has no attribute %s." % (
            repr(self.name), repr(name))
        assert len(x) == 1, "Multiple attributes with same name?"
        return x[0]
    def get_attribute_names(self):
        return [x.name for x in self.attributes]
    #def get_user_input_names(self):
    #    return [x.name for x in self.user_inputs]
    def is_valid_attribute_name(self, name):
        return name in self.get_attribute_names()
    def is_valid_attribute_value(self, name, value):
        attr = self.get_attribute(name)
        return attr.is_valid_value(value)
    def __cmp__(self, other):
        if not isinstance(other, DataType):
            return cmp(id(self), id(other))
        # Bug: should compare attributes and user_inputs without regard
        # to order.
        #x1 = [self.name, self.attributes, self.user_inputs]
        #x2 = [other.name, other.attributes, self.user_inputs]
        x1 = [self.name, self.attributes]
        x2 = [other.name, other.attributes]
        return cmp(x1, x2)
    def __hash__(self):
        #x = self.name, tuple(self.attributes), tuple(self.user_inputs)
        x = self.name, tuple(self.attributes)
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
        for attr in self.attributes:
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
        x += [repr(x) for x in self.attributes]
        #x += [repr(x) for x in self.user_inputs]
        return "DataType(%s)" % ", ".join(x)


class Data:
    # Members:
    # datatype      Datatype object.
    # attributes    Dict of attribute name -> value.

    # Should not be called by the user.  Should always be created from
    # a DataType object.
    def __init__(self, datatype, **keywds):
        # keywds is a dictionary of (attribute or user_input) name ->
        # value (or list of values).
        
        # Make sure values are provided for every attribute.
        attr_names = datatype.get_attribute_names()
        for name in attr_names:
            assert name in keywds, "No value given for %s." % name

        # Make sure the values of the attributes are legal.
        for name, value in keywds.iteritems():
            assert datatype.is_valid_attribute_name(name), \
                   "%s is not a known attribute for datatype %s." % (
                name, datatype.name)
            assert datatype.is_valid_attribute_value(name, value)

        ## attributes = {}
        ## user_inputs = {}
        ## for name, value in keywds.iteritems():
        ##     if name in attr_names:
        ##         attributes[name] = value
        ##     elif name in user_names:
        ##         user_inputs[name] = value
        ##     else:
        ##         raise AssertionError, "Unknown: %s" % name

        self.datatype = datatype
        self.attributes = keywds.copy()
        #self.user_inputs = user_inputs
    def __cmp__(self, other):
        if not isinstance(other, Data):
            return cmp(id(self), id(other))
        #x1 = [self.datatype, self.attributes, self.user_inputs]
        #x2 = [other.datatype, other.attributes, self.user_inputs]
        x1 = [self.datatype, self.attributes]
        x2 = [other.datatype, other.attributes]
        return cmp(x1, x2)
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        #keywds = self.attributes.copy()
        #keywds.update(self.user_inputs)
        x = [
            self.datatype.name,
            #_pretty_attributes(keywds),
            _pretty_attributes(self.attributes),
            ]
        x = [x for x in x if x]
        return "Data(%s)" % ", ".join(x)


class Module:
    def __init__(self, name, in_datatypes, out_datatype, *params):
        # params is a list of Constraint, Consequence, and UserInputDef.
        # objects.
        assert type(name) is type("")

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
        user_inputs = []
        default_attributes_from = []
        for x in params:
            if isinstance(x, Constraint):
                constraints.append(x)
            elif isinstance(x, Consequence):
                consequences.append(x)
            elif isinstance(x, UserInputDef):
                user_inputs.append(x)
            elif isinstance(x, DefaultAttributesFrom):
                default_attributes_from.append(x)
            else:
                raise AssertionError, "invalid parameter: %s" % repr(x)


        # Check default_attributes_from.  Should be None or a valid
        # object.
        assert len(default_attributes_from) <= 1
        x = None
        if default_attributes_from:
            x = default_attributes_from[0]
        default_attributes_from = x
        if default_attributes_from:
            assert len(in_datatypes) > 1
            assert default_attributes_from.input_index < len(in_datatypes)
            assert out_datatype == \
                   in_datatypes[default_attributes_from.input_index]

        # Any checking necessary on UserInput?
            
        self.name = name
        self.in_datatypes = in_datatypes
        self.out_datatype = out_datatype
        self.constraints = constraints
        self.consequences = consequences
        self.default_attributes_from = default_attributes_from
        self.user_inputs = user_inputs

        for x in constraints:
            self._assert_constraint(
                name, in_datatypes, out_datatype, constraints, consequences, x)
        for x in consequences:
            self._assert_consequence(
                name, in_datatypes, out_datatype, constraints, consequences, x)
            
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
            assert len(in_datatypes) > 1, (
                "%r: SAME_AS constraint requires at least two input "
                "datatypes." % name)
            assert constraint.arg1 < len(in_datatypes)
            assert constraint.arg1 != constraint.input_index
            # Make sure there is a MUST_BE or CAN_BE_ANY_OF constraint
            # on constraint.arg1.
            x = constraints
            x = [x for x in x if x.name == constraint.name]
            x = [x for x in x if x.input_index == constraint.arg1]
            assert len(x) > 0, (
                "%r: %r SAME_AS %d, but datatype %d has no constraint on %r." %
                (name, constraint.name, constraint.arg1, constraint.arg1,
                 constraint.name))
            assert len(x) == 1
            x = x[0]
            assert x.behavior in [MUST_BE, CAN_BE_ANY_OF]
            # Make sure the datatype has this attribute.
            dt = in_datatypes[constraint.arg1]
            assert dt.is_valid_attribute_name(constraint.name)
        else:
            raise NotImplementedError

        # For every constraint, there must be a consequent given.
        # Need to specify what the module does with the variable.
        if in_datatype == out_datatype:
            x = [x for x in consequences if x.name == constraint.name]
            assert x, "%r: constraint but no consequence for %r." % (
                name, constraint.name)

    def _assert_consequence(
        self, name, in_datatypes, out_datatype, constraints, consequences,
        consequence):
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
            # the input and output datatypes.
            if cons.behavior in [MUST_BE, CAN_BE_ANY_OF]:
                in_attr = in_datatype.get_attribute(consequence.name)
                out_attr = out_datatype.get_attribute(consequence.name)
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
              self.default_attributes_from, self.user_inputs]
        x2 = [other.name, other.in_datatypes, other.out_datatype,
              other.constraints, other.consequences,
              other.default_attributes_from, other.user_inputs]
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
        x6 = []
        if self.default_attributes_from:
            x6 = [repr(self.default_attributes_from)]
        x7 = [repr(x) for x in self.user_inputs]
        x = [x1, x2, x3] + x4 + x5 + x6 + x7
        x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
        return x


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


def backchain(moduledb, out_data, *user_attributes):
    # Return a Network object.

    if isinstance(out_data, DataType):
        attrdict = {}
        for attr in user_attributes:
            if attr.datatype != out_data:
                continue
            attrdict[attr.name] = attr.value
        out_data = out_data.output(**attrdict)
    assert isinstance(out_data, Data)

    for x in user_attributes:
        assert isinstance(x, Attribute)

    nodes = []        # list of Data or Module objects.
    transitions = {}  # list of index -> list of indexes

    MAX_NETWORK_SIZE = 1024
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


def select_start_node(network, start_data):
    # start_data may be a single Data object or a list of Data
    # objects.  DataTypes are also allowed in lieu of Data objects.


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


def optimize_network(network):
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
            network = opt.optimize(network)
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
    def optimize(self, network):
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
    
    def optimize(self, network):
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
    
    def optimize(self, network):
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
    def optimize(self, network):
        # Remove nodes that have been made irrelevant due to
        # optimizing.
        while True:
            dangling = self.find_dangling_nodes(network)
            if not dangling:
                break
            # Make sure root not is never deleted.
            assert 0 not in dangling
            network = network.delete_nodes(dangling)
        return network
    def find_dangling_nodes(self, network):
        dangling = []
        for (node_id, node) in enumerate(network.nodes):
            if self.is_dangling_node(network, node_id, node):
                dangling.append(node_id)
        return dangling
    def is_dangling_node(self, network, node_id, node):
        # 1.  Module nodes with no inputs.
        # 2.  Module nodes with no outputs.
        # 3.  Data nodes (except for node 0) that don't point to any
        #     Modules.
        if isinstance(node, Module):
            if not network.transitions.get(node_id, []):
                return True
            prev_ids = _backchain_to_ids(network, node_id)
            if not prev_ids:
                return True
        elif isinstance(node, Data):
            if node_id != 0 and not network.transitions.get(node_id, []):
                return True
        else:
            raise AssertionError
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

    def optimize(self, network):
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
            self._fix_overlapping_data(network, node_id1, node_id2, attr_name)
            
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

    def _fix_overlapping_data(self, network, data_id1, data_id2, attr_name):
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
            self._split_list(network, data_id1, data_id2, attr_name)

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


    def _split_list(self, network, data_id1, data_id2, attr_name):
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
            module = network.nodes[module_id]
            if _can_module_take_data(module, [data3]):
                if data_id3 not in network.transitions:
                    network.transitions[data_id3] = []
                network.transitions[data_id3].append(module_id)
        

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
    def optimize(self, network):
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


def _backchain_to_input_old(module, in_num, out_data, user_attributes):
    # Given a module and output_data, return the input_data object
    # that can generate the output.
    assert in_num < len(module.in_datatypes)

    in_datatype = module.in_datatypes[in_num]
    out_datatype = out_data.datatype
    debug_print("Backchaining %s from %s to %s [%d]." % (
        module.name, out_datatype.name, in_datatype.name, in_num))

    # BUG: What if module has two in_datatypes with the same datatype?
    # Which one to copy?
    if in_datatype == out_datatype:
        # Start with the attributes in the out_data.
        attributes = out_data.attributes.copy()
    else:
        # If the datatypes are different, then clear all the values.
        attributes = {}
        
    # Modify the attributes based on the Constraints and Consequences
    # of the module.

    # If there is a Consequence that is SAME_AS_CONSTRAINT, then
    # the attribute should be determined by the out_data.  e.g.
    # Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"])
    # Consequence("quantile_norm", SAME_AS_CONSTRAINT)
    #
    # The module takes anything, produces the same value.  So
    # the backchainer needs to preserve the value from the
    # out_data.
    for consequence in module.consequences:
        n = consequence.name
        if consequence.behavior == SAME_AS_CONSTRAINT:
            # Keep the same attribute.
            if consequence.arg1 == in_num:
                attributes[n] = out_data.attributes[n]
        elif consequence.behavior in [
            SET_TO, SET_TO_ONE_OF, BASED_ON_DATA]:
            # Don't know what it should be.
            if n in attributes:
                del attributes[n]
        else:
            raise AssertionError
        
    debug_print("Taking attributes from out_data %s." % attributes)

    # Now make sure the attributes meet the constraints.
    for constraint in module.constraints:
        if constraint.input_index != in_num:
            continue
        if constraint.behavior == MUST_BE:
            attributes[constraint.name] = constraint.arg1
        elif constraint.behavior == CAN_BE_ANY_OF:
            if constraint.name not in attributes:
                attributes[constraint.name] = constraint.arg1
        elif constraint.behavior == SAME_AS:
            continue
        else:
            raise AssertionError

    # make_out.  This module should take in a "finished" Data object,
    # and use it to generate a new Data object.
    debug_print("Generating a new %s with attributes %s." % (
        in_datatype.name, attributes))
    return in_datatype.output(*user_attributes, **attributes)


def _backchain_to_all_inputs(module, out_data, user_attributes):
    all_inputs = []
    for in_num in range(len(module.in_datatypes)):
        x = _backchain_to_input(module, in_num, out_data, user_attributes)
        all_inputs.append(x)

    # Handle the SAME_AS constraints here.  Difficult to do it in
    # _backchain_to_input because not all inputs have been created
    # yet.
    for constraint in module.constraints:
        if constraint.behavior != SAME_AS:
            continue
        
        # If this constraint is the SAME_AS another one, then use the
        # value of the copied constraint.
        i_src = constraint.arg1
        i_dst = constraint.input_index
        assert i_src < len(all_inputs)
        assert i_dst < len(all_inputs)
        input_src = all_inputs[i_src]
        input_dst = all_inputs[i_dst]

        name = constraint.name
        assert name in input_src.attributes
        assert name in input_dst.attributes
        input_dst.attributes[name] = input_src.attributes[name]

    return all_inputs
        

def _backchain_to_input_new(module, in_num, out_data, user_attributes):
    # Given a module and output_data, return the input_data object
    # that can generate the output.  This should only be called by
    # _backchain_to_all_inputs.
    assert in_num < len(module.in_datatypes)

    in_datatype = module.in_datatypes[in_num]
    out_datatype = out_data.datatype
    debug_print("Backchaining %s <- %s <- %s [%d]." % (
        out_datatype.name, module.name, in_datatype.name, in_num))

    # Start with empty attributes.
    attributes = {}

    # Set the attributes based on the consequences.  If there is a
    # Consequence that is SAME_AS_CONSTRAINT, then the attribute
    # should be determined by the out_data.  e.g.
    # Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"])
    # Consequence("quantile_norm", SAME_AS_CONSTRAINT)
    #
    # The module takes anything, produces the same value.  So
    # the backchainer needs to preserve the value from the
    # out_data.
    for consequence in module.consequences:
        n = consequence.name
        if consequence.behavior == SAME_AS_CONSTRAINT:
            # Keep the same attribute.
            if consequence.arg1 == in_num:
                attributes[n] = out_data.attributes[n]
        elif consequence.behavior in [
            SET_TO, SET_TO_ONE_OF, BASED_ON_DATA]:
            # Don't know what it should be.
            pass
        else:
            raise AssertionError
        
    debug_print("Taking attributes from out_data %s." % attributes)

    # Now set the attributes based on the constraints.
    for constraint in module.constraints:
        if constraint.input_index != in_num:
            continue

        if constraint.behavior == MUST_BE:
            attributes[constraint.name] = constraint.arg1
        elif constraint.behavior == CAN_BE_ANY_OF:
            # If this attribute is not already set by a consequence
            # above (e.g. the output data is already a specific
            # value), then set it here.
            #   Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"])
            #   Consequence("quantile_norm", SAME_AS_CONSTRAINT, 0)
            if constraint.name not in attributes:
                attributes[constraint.name] = constraint.arg1
        elif constraint.behavior == SAME_AS:
            # Handled in _backchain_to_all_inputs.
            pass
        else:
            raise AssertionError

    # Fill in the rest of the attributes.
    # 1.  If the module doesn't convert the datatype, then fill it in
    #     with the same values as the out_data.
    # 2.  If it does convert the datatype, and a
    #     default_attributes_from is the same as in_num, then fill it
    #     in with the same values as the out_data.
    # 3.  If it does convert the datatype, and there is a
    #     user_attribute provided, then use the value from the
    #     user_attribute.
    # 4.  Otherwise, fill it in with the default output values of the
    #     input datatype.
    if len(module.in_datatypes) == 1 and in_datatype == out_datatype:
        def_attributes = out_data.attributes
    elif len(module.in_datatypes) > 1 and \
         module.default_attributes_from and \
         module.default_attributes_from.input_index == in_num:
        assert in_datatype == out_datatype
        def_attributes = out_data.attributes
    else:
        def_attributes = {}
        # Set values from user attributes.
        for attr in user_attributes:
            # Ignore attributes for other data types.
            if attr.datatype != in_datatype:
                continue
            def_attributes[attr.name] = attr.value
        # Set values from defaults.
        for attr in in_datatype.attributes:
            if attr.name not in def_attributes:
                def_attributes[attr.name] = attr.default_out
    for name, value in def_attributes.iteritems():
        if name not in attributes:
            attributes[name] = value

    # make_out.  This module should take in a "finished" Data object,
    # and use it to generate a new Data object.
    debug_print("Generating a %s with attributes %s." % (
        in_datatype.name, attributes))
    return in_datatype.output(**attributes)


#_backchain_to_input = _backchain_to_input_old
_backchain_to_input = _backchain_to_input_new


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


def _can_module_take_data(module, datas):
    # Return True/False if a module can take this list of Data nodes
    # as an input.
    if len(module.in_datatypes) != len(datas):
        return False
    if len(datas) > 1:
        raise NotImplementedError
    data = datas[0]
    if data.datatype != module.in_datatypes[0]:
        return False

    # Make sure the data satisfies each of the module's constraints.
    for constraint in module.constraints:
        assert constraint.name in data.attributes
        data_value = data.attributes.get(constraint.name)
        data_type = _get_attribute_type(data_value)
        assert data_type in [TYPE_ATOM, TYPE_ENUM]
        
        if constraint.behavior == MUST_BE:
            if data_type == TYPE_ATOM:
                if data_value != constraint.arg1:
                    return False
            elif data_type == TYPE_ENUM:
                return False
            else:
                raise AssertionError
        elif constraint.behavior == CAN_BE_ANY_OF:
            if data_type == TYPE_ATOM:
                if data_value not in constraint.arg1:
                    return False
            elif data_type == TYPE_ENUM:
                # data_value contains the possible values of this Data
                # object.  The values that are acceptable by module is
                # in constraint.arg1.  Make sure the module can handle
                # all of the possible values.
                if not _is_subset(data_value, constraint.arg1):
                    return False
            else:
                raise AssertionError
        else:
            raise AssertionError

    return True


def debug_print(s):
    if DEBUG:
        print s


def _can_module_produce_data(module, data):
    # Return whether this module can produce this data object.

    # A module cannot produce this data if:
    # - The module's output data type is not the same as the data.
    # - One or more of the consequences conflict.
    # - The module converts the datatype, and the default values of
    #   the attributes (without consequences) conflict.
    #   THIS RULE IS IN TESTING.
    # 
    # A module can produce this data if:
    # - An consequence (SET_TO, SET_TO_ONE_OF, BASED_ON_DATA) that is not
    #   a side effect has value that matches the value of the data.
    # - The module only converts the datatype.  There are no
    #   consequences (SET_TO, SET_TO_ONE_OF, BASED_ON_DATA), and the
    #   output data type has no attributes.
    #   e.g. download_geo_GSEID  gseid -> expression_files  (no attributes)

    debug_print("Testing if module %s can produce data %s." % (
        repr(module.name), str(data)))

    # If this module doesn't produce the same data type, then it can't
    # produce this data object.
    if module.out_datatype != data.datatype:
        debug_print(
            "Module can't generate data type: %s." % data.datatype.name)
        return False

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
                x = [x for x in module.constraints
                     if x.name == consequence.name]
                x = [x for x in x if x.input_index == constraint.arg1]
                assert len(x) == 1
                constraint = x[0]
            
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
                debug_print("Consequence %s conflicts [%s %s]." % (
                    consequence.name, outc_value, data_value))
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

    # If the module converts the datatype, and no
    # DefaultAttributesFrom is specified, then the data should match
    # the (in) defaults from the output data type.
    if module.in_datatypes != [module.out_datatype] and \
           not module.default_attributes_from:
        debug_print("Module converts datatype.  Checking default attributes.")
        consequence_names = [x.name for x in module.consequences]
        for attr in module.out_datatype.attributes:
            # Ignore the attributes that have consequences.
            if attr.name in consequence_names:
                debug_print(
                    "Attr %r: Skipping--has consequence." % attr.name)
                continue
            assert attr.name in data.attributes
            data_value = data.attributes[attr.name]
            data_type = _get_attribute_type(data_value)
            assert data_type in [TYPE_ATOM, TYPE_ENUM]

            if data_type == TYPE_ATOM:
                if attr.default_in != data_value:
                    debug_print("Attr %r: Conflicts (module %r, data %r)." % (
                        attr.name, attr.default_in, data_value))
                    return False
            elif data_type == TYPE_ENUM:
                if attr.default_in not in data_value:
                    debug_print("Attr %r: Conflicts (module %r, data %r)." % (
                        attr.name, attr.default_in, data_value))
                    return False
            else:
                raise AssertionError
            debug_print("Attr %r: matches defaults." % attr.name)
        

    # TESTING.
    # If the module converts the datatype, the consequences don't
    # conflict, and the default attributes don't conflict, then this
    # should match.
    if module.in_datatypes != [module.out_datatype]:
        debug_print("Match because of converting datatype.")
        return True


    # At this point, the module produces this datatype and there are
    # no conflicts.  Look for an consequence that is not a side effect
    # that matches this data.
    for consequence in module.consequences:
        if consequence.behavior not in [SET_TO, SET_TO_ONE_OF, BASED_ON_DATA]:
            continue
        if consequence.name not in data.attributes:
            continue
        if consequence.side_effect:
            continue
        debug_print("Consequence %s matches." % consequence.name)
        return True

    # No conflicts, and the module has no consequences.
    if not module.consequences:
        debug_print("Match because there are no consequences.")
        return True

    # No consequences match.
    debug_print("No consequences match.")
    return False


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


def _find_start_nodes(network, user_data):
    scores = _score_start_nodes(network, user_data)
    x = [x for x in scores if x[0] == 0]
    node_ids = [x[1] for x in x]
    return node_ids


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


## class Parameter:
##     def __init__(self, module_name, **parameters):
##         # parameters is a dictionary of name -> value.
##         assert type(module_name) is type("")
##         for name, value in parameters.iteritems():
##             assert type(name) is type("")
##         self.module_name = module_name
##         self.parameters = parameters.copy()
##     def __str__(self):
##         return self.__repr__()
##     def __repr__(self):
##         x = [
##             repr(self.module_name),
##             ]
##         for key, value in parameters:
##             x.append("%s=%s" % (key, repr(value)))
##         x = "%s(%s)" % (self.__class__.__name__, ", ".join(x))
##         return x


GEOSeries = DataType("GEOSeries")
ExpressionFiles = DataType("ExpressionFiles")
CELFiles = DataType(
    "CELFiles",
    AttributeDef(
        "version", ["unknown", "cc", "v3", "v4"], "unknown", "v3"),
    )
SignalFile = DataType(
    "SignalFile",
    AttributeDef(
        "format", ["unknown", "tdf", "pcl", "gct", "res", "jeffs"],
        "unknown", "tdf"),

    # Properties of the data.
    AttributeDef(
        "preprocess",
        ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        "unknown", "rma"),
    AttributeDef(
        "missing_values", ["unknown", "no", "yes"], "unknown", "no"),
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "none", "none"),
    AttributeDef("logged", ["unknown", "no", "yes"], "unknown", "yes"),
    
    # This is not necessary.  Remove.
    AttributeDef("filtered", ["no", "yes"], "no", "no"),

    # Normalization of the data.
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no"),

    # Other attributes.
    AttributeDef("predataset", ["no", "yes"], "no", "no"),
    AttributeDef("rename_sample", ["no", "yes"], "no", "no"),
    AttributeDef("contents", [
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],
              "unspecified", "unspecified"),
    )
    

all_modules = [
    Module(
        "download_geo", GEOSeries, ExpressionFiles,
        UserInputDef("GSEID"), UserInputDef("GPLID")),
    Module(
        "extract_CEL_files", ExpressionFiles, CELFiles,
        Consequence("version", SET_TO, "unknown"),
        ),
    Module(
        "detect_CEL_version",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "unknown"),
        Consequence("version", BASED_ON_DATA, ["cc", "v3", "v4"]),
        ),
    Module(
        "convert_CEL_cc_to_CEL_v3",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "cc"),
        Consequence("version", SET_TO, "v3"),
        ),
    Module(
        "preprocess_rma",
        CELFiles, SignalFile,
        Constraint("version", CAN_BE_ANY_OF, ["v3", "v4"]),
        Consequence("logged", SET_TO, "yes"),
        Consequence("preprocess", SET_TO, "rma"),
        Consequence("format", SET_TO, "jeffs"),
        Consequence("missing_values", SET_TO, "no"),
        Consequence("quantile_norm", SET_TO, "yes"),
        ),
    Module(
        "preprocess_mas5",
        CELFiles, SignalFile,
        Constraint("version", CAN_BE_ANY_OF, ["v3", "v4"]),
        Consequence("logged", SET_TO, "no"),
        Consequence("preprocess", SET_TO, "mas5"),
        Consequence("format", SET_TO, "jeffs"),
        Consequence("missing_values", SET_TO, "no"),
        Consequence("quantile_norm", SET_TO, "no"),
        ),
    Module(
        "quantile_normalize",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("quantile_norm", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SET_TO, "yes"),
        ),
    Module(
        "convert_signal_to_tdf",
        SignalFile, SignalFile,
        Constraint("format", CAN_BE_ANY_OF, ['pcl', 'gct', 'res', 'jeffs']),
        Consequence("format", SET_TO, "tdf"),
        ),
    Module(
        "check_for_log",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", BASED_ON_DATA, ["yes", "no"]),
        ),
    Module(
        "log_signal",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SET_TO, "yes"),
        ),
    Module(
        "check_for_missing_values",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("missing_values", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("missing_values", BASED_ON_DATA, ["no", "yes"]),
        ),
    Module(
        "filter_genes_by_missing_values",
        SignalFile, SignalFile,
        UserInputDef("filter_genes_with_missing_values", 0.50),
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("missing_values", MUST_BE, "yes"),
        Constraint("filtered", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("missing_values", SAME_AS_CONSTRAINT),
        Consequence("filtered", SET_TO, "yes")
        ),
    Module(
        "fill_missing_with_zeros",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("missing_values", MUST_BE, "yes"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("missing_values", SET_TO, "no"),
        Consequence(
            "missing_algorithm", SET_TO, "zero_fill", side_effect=True),
        ),

    Module(
        "merge_two_classes", [SignalFile, SignalFile], SignalFile,
        Constraint("contents", MUST_BE, "class0", 0),
        Constraint("format", MUST_BE, "tdf", 0),
        Constraint("logged", MUST_BE, "yes", 0),
        Constraint("preprocess", MUST_BE, "mas5", 0),
        Constraint("contents", MUST_BE, "class1", 1),
        Constraint("format", MUST_BE, "tdf", 1),
        Constraint("logged", MUST_BE, "yes", 1),
        Constraint("preprocess", MUST_BE, "mas5", 1),
        #Constraint("format", SAME_AS, 0, 1),
        #Constraint("logged", SAME_AS, 0, 1),
        #Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("contents", SET_TO, "class0,class1"),
        Consequence("format", SAME_AS_CONSTRAINT, 0),
        Consequence("logged", SAME_AS_CONSTRAINT, 0),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
        DefaultAttributesFrom(0),
        )
    ]

def test_bie():
    #print_modules(all_modules)

    #x = SignalFile(
    #    format="unknown", group_fc=ANYATOM, logged="unknown",
    #    missing_values="unknown", preprocess="unknown")
    #x = SignalFile(preprocess="illumina")
    #in_data = [GEOSeries, ClassLabelFile]
    #in_data = [x, ClassLabelFile]
    #in_data = GEOSeries
    in_data = SignalFile.input(preprocess="rma", format="jeffs")
    #in_data = SignalFile.make_in(
    #    filename="test.txt", logged="yes", preprocess="rma", format="jeffs")
    #in_data = SignalFile.input(
    #    logged="yes", preprocess="rma", format="jeffs", missing_values="no",
    #    missing_algorithm="median_fill", quantile_norm="yes")
        
    #in_data = [
    #    SignalFile(preprocess="rma", format="jeffs", filename='a'),
    #    ClassLabelFile(filename='b')]
    #in_data = GEOSeries
    #x = dict(preprocess="rma", missing_values="no", format="jeffs")
    #in_data = [SignalFile(contents='class0',logged="yes", preprocess="rma"),
    #           SignalFile(contents='class1',logged="yes", preprocess="rma")]
    #in_data = [
    #    SignalFile(contents='class0', preprocess="rma", format="jeffs"),
    #    SignalFile(contents='class1', preprocess="rma", format="jeffs")]

    #print _make_goal(SignalFile, x)
    #return

    #out_data = CELFiles
    #out_data = SignalFile.output(
    #    format="tdf", logged=["no", "yes"], preprocess=["rma", "mas5"],
    #    quantile_norm=["no", "yes"])
    #out_data = SignalFile.output(
    #    format="tdf", preprocess="mas5", logged=["no", "yes"],
    #    missing_values="no", quantile_norm=["no", "yes"])
    out_data = SignalFile.output(
        format="tdf", preprocess="mas5", logged="yes",
        #missing_values="unknown", 
        missing_values="no", quantile_norm="no",
        contents="class0,class1")
    #parameters = [
    #    Parameter("download_geo", GSEID="GSE2034", GPLID="GPL9196"),
    #    Parameter("filter_genes_by_missing_values", filter=0.50),
    #    ]
    #out_data = SignalFile(format='tdf', logged='no', missing_values="no",
    #    missing_algorithm="median_fill")
    
    #goal_datatype = CELFiles
    #goal_attributes = {}
    #goal_attributes = dict(platform="GPL1691")
    #goal_attributes = dict(version=["v3", "v4"])

    #goal_datatype = ILLUFolder
    #goal_attributes = dict(
    #    ill_bg_mode='true',
    #    ill_chip='ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
    #    ill_clm='',
    #    ill_coll_mode='none',
    #    ill_custom_chip='',
    #    ill_custom_manifest='',
    #    illu_manifest='HumanHT-12_V4_0_R2_15002873_B.txt')

    #goal_datatype = ExpressionFiles
    #goal_attributes = {}

    #goal_datatype = SignalFile
    #goal_attributes = dict(format='tdf')
    #goal_attributes = dict(
    #    format='tdf', preprocess="rma", logged='yes', missing_values="no",
    #    missing_algorithm="median_fill")
    #goal_attributes = dict(
    #    format=['jeffs', 'gct'], preprocess='rma', logged='yes',
    #    missing_values="no", missing_algorithm="median_fill")
    #goal_attributes = dict(
    #    format='tdf', preprocess='illumina', logged='no', missing_values="no",
    #    missing_algorithm="median_fill")
    #goal_attributes = dict(contents='class0,class1',
    #    format='tdf', preprocess='rma', logged='yes',
    #    missing_values="no")
    #goal_attributes = dict(
    #    #format=['tdf', 'pcl', 'gct', 'res', 'jeffs', 'unknown', 'xls'],
    #    format='gct',
    #    preprocess='illumina', logged='no',
    #    missing_values="no")
    #goal_attributes = dict(
    #    format='tdf', preprocess='illumina', logged='yes',
    #    missing_values="no", quantile_norm="yes", combat_norm="yes")
    #goal_attributes = dict(
    #    format='tdf', preprocess='illumina', logged='yes',
    #    missing_values="no", group_fc="5")
    #goal_attributes = dict(
    #    format='tdf', preprocess="rma", logged='yes', 
    #    missing_values="no")
    #goal_attributes = dict(format='tdf', preprocess='rma', logged='yes')
    #goal_attributes = dict(
    #    format='tdf', preprocess='rma', logged='yes',
    #    quantile_norm="yes", combat_norm="yes", dwd_norm="yes",
    #    missing_values="no",gene_order='gene_list')
    #goal_attributes = dict(
    #    format='tdf', preprocess='illumina', logged='yes')
    #goal_attributes = dict(
    #    format='tdf', preprocess='illumina', logged='yes',
    #    #quantile_norm="yes", combat_norm="yes", dwd_norm="yes",
    #    missing_values="zero_fill", group_fc="5", ill_bg_mode='true',
    #    illu_coll_mode="max")
    #goal_attributes = dict(
    #    format='tdf', preprocess='illumina', logged='yes',
    #    #quantile_norm="yes", combat_norm="yes", dwd_norm="yes",
    #    missing_values="zero_fill", platform='hg19',
    #    duplicate_probe='closest_probe')
    #goal_attributes = dict(
    #    format='tdf', preprocess='illumina', logged='yes',
    #    quantile_norm="no", combat_norm="no", dwd_norm="no",
    #    missing_values="zero_fill", platform='hg19',
    #    duplicate_probe='closest_probe')
    #goal_attributes = dict(
    #    format='tdf', preprocess='illumina', logged='yes',
    #    quantile_norm="yes", combat_norm="yes", dwd_norm="yes",
    #    missing_values="zero_fill", platform='hg19',
    #    duplicate_probe='closest_probe')

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
    
    network = backchain(all_modules, out_data)
    network = optimize_network(network)
    #network = select_start_node(network, in_data)
    #diagnose_start_node(network, in_data)
    #network = prune_network_by_internal(
    #    network, SignalFile(quantile_norm="yes", combat_norm="no"))
    #network = prune_network_by_internal(
    #    network, SignalFile(combat_norm="yes", dwd_norm="no"))

    print "INPUT:"
    print in_data
    print
    
    print "OUTPUT:"
    print out_data
    print
    
    print_network(network)
    plot_network_gv("out.png", network)

    # Want function _find_data_node with options:
    # ignore_defaults
    # allow_more_general
    #node = _make_goal(
    #    SignalFile, dict(format="pcl", preprocess="rma"))
    #node_id = _find_data_node(network.nodes, node)
    #print "NODE:", node_id

    #node = _make_goal(
    #    SignalFile, dict(format=["tdf", "pcl", "gct"], preprocess="rma"))
    #node_id = _find_data_node(network.nodes, node)
    #print "NODE:", node_id



if __name__ == '__main__':
    test_bie()
    #import cProfile; cProfile.run("test_bie()")
