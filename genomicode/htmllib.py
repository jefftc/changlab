"""

Functions:
HEAD
TITLE

P
BR
HR

B
H1
H2
H3
EM
SPAN
FONT
CENTER

UL
LI

TABLE
TR
TH
TD

A
IMG

"""
# _format_attrs
# _format_se
# _format_s
# _join_sme

def HEAD(s):
    return _join_sme("HEAD", s, join="\n")

def TITLE(s):
    return _join_sme("TITLE", s)

def P():
    return "<P>"

def BR():
    return "<BR>"

def HR():
    return "<HR>"

def B(s):
    return _join_sme("B", s)

def H1(s):
    return _join_sme("H1", s)

def H2(s):
    return _join_sme("H2", s)

def H3(s):
    return _join_sme("H3", s)

def EM(s):
    return _join_sme("EM", s)

def SPAN(s, style=None):
    attrs = [
        ("STYLE", style),
        ]
    return _join_sme("SPAN", s, attrs=attrs)

def FONT(s, size=None):
    attrs = [
        ("SIZE", size),
        ]
    return _join_sme("SIZE", s, attrs=attrs)

def CENTER(s):
    return _join_sme("CENTER", s)

def UL(s):
    return _join_sme("UL", s, join="\n")

def LI():
    return "<LI>"

def TABLE(s, border=None, cellpadding=None, cellspacing=None, width=None):
    attrs = [
        ("BORDER", border),
        ("CELLPADDING", cellpadding),
        ("CELLSPACING", cellspacing),
        ("WIDTH", width),
        ]
    return _join_sme("TABLE", s, attrs=attrs, join="\n")

def TR(s):
    return _join_sme("TR", s, join="\n")

def TH(s, colspan=None, align=None, valign=None):
    attrs = [
        ("COLSPAN", colspan),
        ("ALIGN", align),
        ("VALIGN", valign),
        ]
    return _join_sme("TH", s, attrs=attrs, join="\n")

def TD(s, colspan=None, align=None, valign=None):
    attrs = [
        ("COLSPAN", colspan),
        ("ALIGN", align),
        ("VALIGN", valign),
        ]
    return _join_sme("TD", s, attrs=attrs, join="\n")

def A(s, href=None):
    attrs = [
        ("HREF", href),
        ]
    return _join_sme("A", s, attrs=attrs, join="\n")

def IMG(height=None, src=None):
    attrs = [
        ("HEIGHT", height),
        ("SRC", src),
        ]
    return _format_s("IMG", attrs=attrs)

def _format_attrs(attrs):
    if not attrs:
        return ""

    attrs_list = []
    for k, v in attrs:
        if " " in str(v):
            v = "'%s'" % v
        x = "%s=%s" % (k, v)
        attrs_list.append(x)
    return " ".join(attrs_list)

def _format_se(tag, attrs=[]):
    # attrs should be a list of (key, value).  If value is None, then
    # skip it.
    attrs = [x for x in attrs if x[1] is not None]
    
    START = "<%s>" % tag
    if attrs:
        START = "<%s %s>" % (tag, _format_attrs(attrs))
    END = "</%s>" % tag
    return START, END

def _format_s(tag, attrs=[]):
    attrs = [x for x in attrs if x[1] is not None]
    
    START = "<%s>" % tag
    if attrs:
        START = "<%s %s>" % (tag, _format_attrs(attrs))
    return START

def _join_sme(tag, m, attrs=[], join=""):
    START, END = _format_se(tag, attrs=attrs)
    return "%s%s%s%s%s" % (START, join, m, join, END)

