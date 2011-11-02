"""

Functions:
fix_environ_path        PATH messed up within GenePattern.  Fix this.

format_parameter_info   Write a list of ParameterInfo objects to GP format.
parse_parameter_info
format_rtbparameter     Write a list of RTBParameter objects in JSON format.
parse_rtbparameter

param2elem
elem2param
param2rtbparam
elem2rtbparam

Classes:
ParameterInfo
RTBParameter
Element
EFileChooser
EDropdown
ETextBox

"""

# _parse_formats
# _format_formats
#
# _parse_choices
# _format_choices

# These are listed in the order that they appear in the XML file.
ATTRIBUTES = [
    "fileFormat", "TYPE", "MODE", "optional", "type",
    "domain", "default_value", "prefix_when_specified"]
    
class ParameterInfo:
    def __init__(self, name, value, description, attributes):
        self.name = name
        self.value = value
        self.description = description
        self.attributes = attributes

class RTBParameter:
    # Parameter class from RunTaskBean.
    def __init__(
        self, name, displayName, displayDesc, optional, inputType,
        defaultValue, choices):
        # name          Name of the parameter.  "rma_expression_file"
        # displayName   Usually the same as <name>.
        # displayDesc   Description for the user.
        # optional      Is this optional?  True/False.
        # inputType     "password", "file", "select" or "text".
        # defaultValue  Not necessary for "select" inputs.
        # choices       List of dictionaries.  Each dict describes one option.
        #
        # In choices, each option is a dictionary of:
        #   value           Passed to module.            "HIERARCHICAL"
        #   label           Pretty label shown to user.  "Hierarchical"
        #   defaultOption   True/False

        # Do some basic checking.
        if inputType == "select":
            assert choices
        self.name = name
        self.displayName = displayName
        self.displayDesc = displayDesc
        self.optional = optional
        self.inputType = inputType
        self.defaultValue = defaultValue
        self.choices = choices

# For Element.
TEXT, INTEGER, FLOAT, PASSWORD = range(4)

class Element:
    def __init__(self, name, description, default=None, optional=None):
        optional = optional or False
        
        self.name = name
        self.description = description
        self.default = default
        self.optional = optional

class EFileChooser(Element):
    def __init__(
        self, name, description, formats, default=None, optional=None,
        prefix_when_specified=None):
        Element.__init__(
            self, name, description, default=default, optional=optional)
        assert type(formats) is type([])
        self.formats = formats
        self.prefix_when_specified = prefix_when_specified

class EDropdown(Element):
    def __init__(
        self, name, description, choices, formats=None, default=None,
        optional=None, input_type=None):
        Element.__init__(
            self, name, description, default=default, optional=optional)
        assert formats is None or type(formats) is type([])
        self.formats = formats
        if input_type is None:
            input_type = TEXT
        assert input_type in range(4)
        assert type(choices) is type([])
        
        # Make sure choices is a list of (value, pretty_value).  If an
        # item is a singleton, then make value and pretty_value the
        # same.
        for i in range(len(choices)):
            if type(choices[i]) is type(()):
                continue
            assert type(choices[i]) is type(""), "choice must be a string."
            value = choices[i]
            choices[i] = value, value
        self.choices = choices
        self.input_type = input_type

        # Make sure default is one of the choices.
        if default:
            for value, pretty_value in choices:
                if value == default:
                    break
            else:
                raise AssertionError, (
                    "default value [%s] is not a choice" % default)
        
        
class ETextBox(Element):
    def __init__(
        self, name, description, formats=None, default=None, optional=None,
        input_type=None, prefix_when_specified=None):
        Element.__init__(
            self, name, description, default=default, optional=optional)
        formats = formats or []
        assert type(formats) is type([])
        self.formats = formats
        if input_type is None:
            input_type = TEXT
        assert input_type in range(4)
        self.input_type = input_type
        self.prefix_when_specified = prefix_when_specified

def fix_environ_path(binpath=[]):
    # For some reason, GenePattern screws up the PATH environment
    # variable.  It's automatically set to:
    # '/sbin:/usr/sbin:/bin:/usr/bi /opt/GenePatternServer/taskLib/  \
    #  BinReg.20.204'
    # Detect this and fix it.
    import os
    
    if os.environ["PATH"].find("GenePatternServer") < 0 and not binpath:
        return

    # Make sure binpath is specified correctly.
    for x in binpath:
        assert ":" not in x

    # BUG: Should allow configuration of path, java path.
    PATH = [
        #"/usr/java/jdk1.5.0_14/bin",   # Need Java for matlab.
        "/usr/local/bin",
        "/usr/bin",
        "/bin",
        #"/sbin",
        ]
    PATH = PATH + binpath
    
    # /usr/java/jdk1.5.0_14/bin:/usr/local/bin:/bin:/usr/bin:
    os.environ["PATH"] = ":".join(PATH)

def format_parameter_info(params):
    import StringIO

    handle = StringIO.StringIO()
    w = handle.write
    w('<?xml version="1.0" encoding="UTF-8"?>\n')
    w("\n")
    w("<ANALYSISPARAMETERS>\n")
    need_newline1 = False
    for param in params:
        if need_newline1:
            w("\n")
        need_newline1 = True
        w('  <PARAMETER name="%s" value="%s">\n' % (param.name, param.value))
        # Should do a better job of this.
        x = param.description
        x = x.replace("&", "&amp;")
        x = x.replace(">", "&gt;")
        x = x.replace("<", "&lt;")
        w('    <DESCRIPTION>%s</DESCRIPTION>\n' % x)
        # Make sure all the attributes are known.
        for key in param.attributes:
            assert key in ATTRIBUTES
        need_newline2 = False
        for key in ATTRIBUTES:
            if key not in param.attributes:
                continue
            if need_newline2:
                w("\n")
            need_newline2 = True
            value = param.attributes[key]
            w('    <ATTRIBUTE key="%s">%s</ATTRIBUTE>' % (key, value))
        w('</PARAMETER>')
    w("</ANALYSISPARAMETERS>\n")
    handle.seek(0)
    return handle.read()

def parse_parameter_info(parameter_info):
    # Parse an XML formatted string into a list of ParameterInfo
    # objects.
    from xml.parsers import expat

    class Parser:
        def __init__(self):
            self._reset_parameters()
            self.params = []
        def _reset_parameters(self):
            self.name = self.value = None
            self.description = ""
            self.attrs = {}   # dict of key -> value
            self._tag = None
            self._attr_key = ""
            self._attr_value = ""
        def start_element(self, name, attrs):
            # Wrapped around all the parameters.
            if name == "ANALYSISPARAMETERS":
                return
            elif name == "PARAMETER":
                assert self.name is None
                assert self.value is None
                self.name = attrs["name"]
                self.value = attrs["value"]
            elif name == "DESCRIPTION":
                assert not self.description
                assert self._tag is None
                self._tag = name
            elif name == "ATTRIBUTE":
                assert self._tag is None
                assert self._attr_key is ""
                self._tag = name
                self._attr_key = attrs["key"]
            #print 'Start element:', name, attrs
        def end_element(self, name):
            if self._tag == "DESCRIPTION":
                pass
            elif self._tag == "ATTRIBUTE":
                assert self._attr_key
                # _attr_value might be blank.
                self.attrs[self._attr_key] = self._attr_value
                self._attr_key = ""
                self._attr_value = ""
            if name == "PARAMETER":
                x = ParameterInfo(
                    name=self.name, value=self.value,
                    description=self.description, attributes=self.attrs)
                self.params.append(x)
                self._reset_parameters()
            self._tag = None
            #print 'End element:', name
        def char_data(self, data):
            if self._tag == "DESCRIPTION":
                self.description = self.description + str(data)
            elif self._tag == "ATTRIBUTE":
                self._attr_value = self._attr_value + str(data)
            #print 'Character data:', repr(data)

    parser = Parser()
    p = expat.ParserCreate()
    p.StartElementHandler = parser.start_element
    p.EndElementHandler = parser.end_element
    p.CharacterDataHandler = parser.char_data
    p.Parse(parameter_info, True)
    return parser.params

def format_rtbparameter(params):
    import json
    
    x = json.dumps([p.__dict__ for p in params])
    return x

def parse_rtbparameter(rtbparameter_str):
    # Parse a JSON formatted string into a list of RTBParameter objects.
    import json
 
    x = json.loads(rtbparameter_str)
    params = [RTBParameter(**x) for x in x]
    return params

def param2elem(param):
    # Figure out which type of parameter this is.
    # CONDITION                           TYPE
    # <value> given                       EDropdown
    # attributes["type"] == java.io.File  EFileChooser
    # otherwise                           ETextBox
    attrs = param.attributes

    # Make sure I know how to interpret each attribute in this
    # parameter.
    for k in attrs:
        assert k in ATTRIBUTES, "Unknown attribute: %s" % k

    if param.value:
        elem_type = EDropdown
    elif attrs["type"] == "java.io.File":
        elem_type = EFileChooser
    else:
        elem_type = ETextBox

    attr_type2input_type = {
        "java.lang.String" : TEXT,
        "java.lang.Integer" : INTEGER,
        "java.lang.Float" : FLOAT,
        "PASSWORD" : PASSWORD,
        }

    elem = None
    if elem_type == EDropdown:
        assert param.value
        choices = _parse_choices(param.value)
        # Sometimes fileFormat is missing, sometimes it's empty.
        formats = _parse_formats(attrs.get("fileFormat", ""))
        assert "TYPE" not in attrs
        # Sometimes MODE is missing, sometimes it's empty or "IN".
        assert attrs.get("MODE", "") in ["", "IN"]
        optional = attrs["optional"]
        assert attrs["type"] in attr_type2input_type, \
               "Unknown attribute type: %s" % attrs["type"]
        input_type = attr_type2input_type[attrs["type"]]
        assert attrs.get("domain", "") == ""
        default = attrs["default_value"]
        assert attrs["prefix_when_specified"] == ""
        elem = EDropdown(
            param.name, param.description, choices, default=default,
            optional=optional, input_type=input_type)
    elif elem_type == EFileChooser:
        assert not param.value
        # Sometimes fileFormat is missing.
        formats = _parse_formats(attrs.get("fileFormat", ""))
        assert attrs["TYPE"] == "FILE"
        assert attrs["MODE"] == "IN"
        optional = attrs["optional"]
        assert attrs["type"] == "java.io.File"
        assert attrs.get("domain", "") == ""
        # Sometimes default_value is missing.
        default = attrs.get("default_value", "")
        prefix_when_specified = attrs["prefix_when_specified"]
        elem = EFileChooser(
            param.name, param.description, formats, default=default,
            optional=optional, prefix_when_specified=prefix_when_specified)
    elif elem_type == ETextBox:
        assert not param.value
        # Sometimes fileFormat is missing, sometimes it's empty.
        formats = _parse_formats(attrs.get("fileFormat", ""))
        assert "TYPE" not in attrs
        # Sometimes MODE is missing, sometimes it's empty or "IN".
        assert attrs.get("MODE", "") in ["", "IN"]
        optional = attrs["optional"]
        # type might be missing.  If it's missing, then assume String.
        attrs_type = attrs["type"]
        if not attrs_type:
            attrs_type = "java.lang.String"
        assert attrs_type in attr_type2input_type, \
               "Unknown attribute type: %s" % attrs_type
        input_type = attr_type2input_type[attrs_type]
        assert attrs.get("domain", "") == ""
        default = attrs["default_value"]
        prefix = attrs["prefix_when_specified"]
        elem = ETextBox(
            param.name, param.description, formats=formats, default=default,
            optional=optional, input_type=input_type, 
            prefix_when_specified=prefix)
    else:
        raise AssertionError, "Unknown elem_type."
    return elem

def elem2param(elem):
    name = elem.name
    value = ""
    description = elem.description
    attributes = {}
    attributes["default_value"] = elem.default or ""
    attributes["optional"] = ""
    if elem.optional:
        attributes["optional"] = "on"
    
    input_type2attr_type = {
        TEXT : "java.lang.String",
        INTEGER : "java.lang.Integer",
        FLOAT : "java.lang.Float",
        PASSWORD : "PASSWORD",
        }

    if elem.__class__ is EFileChooser:
        attributes["fileFormat"] = _format_formats(elem.formats)
        attributes["TYPE"] = "FILE"
        attributes["MODE"] = "IN"
        attributes["type"] = "java.io.File"
        attributes["prefix_when_specified"] = elem.prefix_when_specified or ""
    elif elem.__class__ is EDropdown:
        value = _format_choices(elem.choices)
        # Leave fileFormat empty unless it's explicitly specified.
        if elem.formats:
            attributes["fileFormat"] = _format_formats(elem.formats)
        attributes["type"] = input_type2attr_type[elem.input_type]
        attributes["prefix_when_specified"] = ""
    elif elem.__class__ is ETextBox:
        # Leave fileFormat empty unless it's explicitly specified.
        if elem.formats:
            attributes["fileFormat"] = _format_formats(elem.formats)
        attributes["type"] = input_type2attr_type[elem.input_type]
        attributes["prefix_when_specified"] = elem.prefix_when_specified
    else:
        raise AssertionError, "Unknown element: %s." % elem.__class__
    param = ParameterInfo(name, value, description, attributes)
    return param

def param2rtbparam(param):
    name = param.name
    
    x = param.name
    x.replace(".", " ")
    displayName = x
    
    displayDesc = param.description
    optional = param.attributes.get("optional", False)
    defaultValue = param.attributes.get("default_value", "")
    
    choices = []
    if param.value:
        for x in _parse_choices(param.value):
            value, pretty_value = x
            x = {}
            x["value"] = value
            x["label"] = pretty_value
            d = defaultValue in [value, pretty_value]
            x["defaultOption"] = d
            choices.append(x)
    
    # Should be: "password", "file", "select" or "text"
    if choices:
        inputType = "select"
    elif param.attributes.get("type") == "java.io.File":
        inputType = "file"
    elif param.attributes.get("type") == "PASSWORD":
        inputType = "password"
    else:
        inputType = "text"

    x = RTBParameter(
        name, displayName, displayDesc, optional, inputType, defaultValue,
        choices)
    return x

def elem2rtbparam(elem):
    name = elem.name
    displayName = name
    
    displayDesc = elem.description
    optional = elem.optional
    defaultValue = elem.default or ""

    choices = []
    
    if elem.__class__ is ETextBox:
        inputType = "text"
        if elem.input_type == PASSWORD:
            inputType = "password"
    elif elem.__class__ is EFileChooser:
        inputType = "file"
    elif elem.__class__ is EDropdown:
        inputType = "select"
        for value, pretty_value in elem.choices:
            x = {}
            x["value"] = value
            x["label"] = pretty_value
            x["defaultOption"] = elem.default in [value, pretty_value]
            choices.append(x)
    else:
        raise AssertionError, "Unknown element: %s." % elem.__class__

    x = RTBParameter(
        name, displayName, displayDesc, optional, inputType, defaultValue,
        choices)
    return x

def _parse_formats(formats_str):
    # Semicolon-separated list of choices.
    # Examples:
    # gct;res;Dataset
    import urllib
    x = urllib.unquote(formats_str)
    return x.split(";")

def _format_formats(formats):
    import urllib
    x = ";".join(formats)
    return urllib.quote(x)

def _parse_choices(choices_str):
    # Parse the choices string and return a list of (value,
    # pretty_value).
    # 
    # Semicolon-separated list of choices.
    # Examples:
    # <value>=<pretty_value>;<value>=<pretty_value>
    # HIERARCHICAL=Hierarchical;SOM;NMF;KMEANS=KMeans
    # -q=Yes;=No
    #
    # <value> can be blank, indicating a blank parameter will be
    # passed to the module if this option is chosen.
    # =<pretty_value> is optional.  If not specified, then
    # <pretty_value> is the same as <value>.
    import urllib
    
    x = urllib.unquote(choices_str)
    assert ";" in x, "no choices"
    choices = []
    for x in x.split(";"):
        x = str(x)   # no unicode.
        value = pretty = x
        if "=" in x:
            value, pretty = x.split("=", 1)
        choices.append((value, pretty))
    return choices

def _format_choices(choices):
    import urllib
    choices_str = []
    for value, pretty in choices:
        x = value
        if value != pretty:
            x = "%s=%s" % (value, pretty)
        choices_str.append(x)
    x = ";".join(choices_str)
    x = urllib.quote(x)
    return x
