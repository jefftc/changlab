"""
Parse strings that come from the user.  Also does some formatting.

Functions:
parse_ranges

pretty_int
pretty_range
pretty_pvalue

HTML PARSING
remove                 Remove text from a string.
remove_tags            Remove HTML tags from a string.
remove_tags_and_contents  Remove HTML tags and the text between them.
remove_all_tags        Remove all HTML tags from a string.
remove_comments        Remove comments from HTML.

get_tags_and_contents  Retrieve HTML tags and the text between them.
get_contents           Get just the text between HTML tags.
get_first_contents     Get the text between the first HTML tag found.

find_start_end         Find the indexes between two pieces of text.
get_href               Extract the URL from an HREF.
get_src                Extract the SRC from an IMG.
resolve_relative_links Make relative links into absolute links.

collapse               Remove extraneous whitespace from a string.
clean                  Remove tags and extraneous whitespace.

"""
import re

def parse_ranges(s):
    """Return a list of (start, end) indexes, according to what is
    specified in s.  It is up to the client to interpret the semantics
    of the ranges.

    Samples:
    16           [(16, 16)]
    1-6          [(1, 6)]
    1,3,5-8,99   [(1, 1), (3, 3), (5, 8), (99, 99)]

    """
    import re
    assert type(s) is type(""), "s must be string"
    assert not re.search(r"[\r\n]", s), "linebreaks not allowed"
    s = re.sub(r"\s+", "", s)    # no spaces
    assert not re.search(r"[^0-9,-]", s), "invalid characters"
    
    ranges = []
    for r in s.split(","):
        # Case 1:  1       ['1']
        # Case 2:  1-2     ['1', '2']
        # Case 3:  -1      ['', '1']
        # Case 4:  -1-5    ['', '1', '5']
        # Case 5:  -5-1    ['', '5', '1']
        # Case 6:  -5--1   ['', '5', '', '1']
        # Case 7:  --1     ['', '', '5']        ERROR
        # Case 8:  1--5    ['1', '', '5']       ERROR
        # Case 9:  5-      ['5', '']            ERROR
        parts = r.split("-")

        # Blanks should be interpreted as negatives for the next
        # digit.
        i = 0
        while i < len(parts):
            if parts[i] == '':
                assert i < len(parts)-1, "Invalid range: %s" % r
                assert parts[i+1] != "", "Invalid range: %s" % r
                parts[i+1] = "-" + parts[i+1]
                del parts[i]
            else:
                i += 1
        
        assert len(parts) <= 2, "I do not understand ranges %s" % r
        if len(parts) == 1:
            start = int(parts[0])
            stop = start
        else:
            start, stop = map(int, parts)
            stop = stop
        assert start <= stop, "Invalid range %s" % r
        ranges.append((start, stop))
    return ranges

def test_parse_ranges():
    tests = [
        ("16", [(16, 17)]),
        ("1-6", [(1, 7)]),
        ("2-3, 4-5", [(2, 4), (4, 6)]),
        ("-1", [(-1, 0)]),
        ("-5--1,-6", [(-5, 0), (-6, -5)]),
        ("-5-", "exception"),
        ("1--5", "exception"),
        ("1-2;3-4", "exception"),
        ("12-4\r34", "exception"),
        ("15-3", "exception"),
        ]

    for s, gold_standard in tests:
        try:
            output = parse_ranges(s)
        except Exception, x:
            output = "exception"
            #raise
        status = "PASSED"
        if str(output) != str(gold_standard):
            status = "FAILED"
        x = s, output, gold_standard, status
        print "\t".join(map(str, x))

def pretty_int(num):
    num = int(num)
    num = str(num)
    # Split into groups of three characters.
    groups = []
    while num:
        groups.insert(0, num[-3:])
        num = num[:-3]
    # If the first character is a "-", make sure it's not alone.
    assert groups
    if groups[0] == "-":
        groups[1] = groups[0] + groups[1]
        del groups[0]
    pretty = ",".join(groups)
    return pretty

def pretty_range(start, stop):
    # Return a list of numbers from start to stop-1.  e.g.  "000",
    # "001", ... <stop-1>
    import math
    assert stop > start
    assert start >= 0, "Negative numbers not implemented."
    ndigits = 1   # Default handles stop==1 boundary case.
    if stop >= 10:
        # Doesn't work due to floating point representation.
        # e.g. math.log(100, 10) == 2, math.log(1000, 10) == 2.999999996
        # ndigits = int(math.floor(math.log(stop-1, 10)))+1
        ndigits = len(str(stop-1))
    x = ["%0*d" % (ndigits, i) for i in range(start, stop)]
    return x

def pretty_pvalue(pvalue, nsig=1):
    import math

    assert pvalue >= 0 and pvalue <= 1, pvalue
    assert nsig >= 1
    digits_after_decimal = int(math.ceil(-math.log(pvalue, 10)))
    digits_after_decimal += nsig - 1  # add more significant digits
    x = "%.*f" % (digits_after_decimal, pvalue)
    return x

# NOTE: The functions in this package use quick and dirty regular
# expressions to handle HTML.  There are cases where the code can get
# fooled by clever or obscure patterns.  For example, it looks for
# HTML tags using the expresion "<[^>]*>", which will fail if there
# are special characters, e.g. <A HREF="foo>bar">.  There are ways to
# handle this better (e.g. using a real HTML parser), but nearly all
# of them result in more code complexity or performance loss.

def _make_tag_patterns(tag):
    # Return patterns for "<tag ...>" and "</tag>"
    return r"<%s\b[^>]*>" % re.escape(tag), r"</%s>" % re.escape(tag)

def remove(html, *args):
    """remove(html[, str][...]) -> html

    Return the html with each string str removed (replaced with a
    blank space).

    """
    reobj = re.compile("|".join(args), re.IGNORECASE)
    return reobj.sub(" ", html)

def remove_tags(html, *tags):
    """remove_tags(html[, tag][...]) -> html

    Return the html with each tag removed (replaced with a blank
    space).

    """
    to_remove = []
    for t in tags:
        to_remove.extend(list(_make_tag_patterns(t)))
    return remove(html, *to_remove)

def remove_tags_and_contents(html, *tags):
    """remove_tags_and_contents(html[, tag][...]) -> html

    Return the html with the tags and the text between them removed
    (replaced with a blank space).

    """
    # remove *tags and the stuff between them.  tag is "table"
    for tag in tags:
        start, end = _make_tag_patterns(tag)
        reobj = re.compile(r"%s.*?%s" % (start, end), re.IGNORECASE|re.DOTALL)
        html = reobj.sub(" ", html)
    return html

def remove_all_tags(html):
    """remove_all_tags(html) -> html

    Remove all tags from an HTML string.

    """
    # remove all <...>
    reobj = re.compile(r"<[^>]*>", re.IGNORECASE|re.DOTALL)
    return reobj.sub(" ", html)

def remove_comments(html):
    """remove_comments(html) -> html, without comments"""
    return re.sub(r"<!--.*?-->", " ", html)

def get_contents(html, *tags):
    """get_contents(html[, tag][...]) -> list of strings

    Return the text between open and close tags.  For example:
    get_contents("<B>Hi!</B><B>Bye</B>", "b") -> ["Hi!", "Bye"]

    """
    # return the stuff in between the open and close tags.
    tags_and_contents = get_tags_and_contents(html, *tags)
    found = []
    for tac in tags_and_contents:
        tac = tac[tac.find(">")+1:]   # strip off the tags
        tac = tac[:tac.rfind("<")]
        found.append(tac)
    return found
    
def get_first_contents(html, *tags):
    """get_first_contents(html[, tag][...]) -> string or None

    Similar to get_contents, but only returns the contents of the
    first tag matched.

    """
    # return the stuff in between the first tag found or None
    x = get_contents(html, *tags)
    if not x:
        return None
    return x[0]

def get_tags_and_contents(html, *tags):
    """get_tags_and_contents(html[, tag][...]) -> list of strings

    Return the tags found and the text between them.  For example:
    get_contents("<B>Hi!</B><B>Bye</B>", "b") -> ["<B>Hi!</B>", "<B>Bye</B>"]

    """
    # return the tags found, plus the stuff between them
    found = []
    for tag in tags:
        start, end = _make_tag_patterns(tag)
        i = 0
        while 1:
            m = re.search(start, html[i:], re.IGNORECASE|re.DOTALL)
            if not m:
                break
            s = m.start()+i
            
            m = re.search(end, html[s:], re.IGNORECASE|re.DOTALL)
            if not m:
                break
            e = m.end()+s
            
            found.append(html[s:e])
            i = e
        # This segfaults sometimes on Python 2.2!
        #reobj = re.compile(r"%s.*?%s" % (start, end), re.IGNORECASE|re.DOTALL)
        #found.extend(reobj.findall(html))
    return found

def find_start_end(text, start_text, end_text, start=0):
    """find_start_end(text, start_text, end_text[, start]) -> (start, end) or None

    Find the start and end indexes of a block of text that starts with
    start_text and ends with end_text.  start_text and end_text are
    included in the indexes.  Return None if not found.

    """
    # return (s, e) or None
    s = text.find(start_text, start)
    if s < 0:
        return None
    e = text.find(end_text, s+1)
    if e < 0:
        return None
    e += len(end_text)
    return s, e

def collapse(s):
    """collapse(s) -> s, with runs of whitespace replaced with single spaces"""
    return ' '.join(s.split()).strip()

def clean(s):
    return collapse(remove_all_tags(s))

def get_href(text, base_url=None):
    """get_href(text[, base_url]) -> href or None

    Extract the URL out of an HREF tag.  If base_url is provided,
    will attempt to resolve relative links.

    """
    m = re.search(r'href\s*=\s*["\']?([^"\'> ]+)["\'> ]', text, re.I)
    if not m:
        return None
    link = m.group(1).strip()
    if base_url and not link.lower().startswith("http"):
        import urlparse
        link = urlparse.urljoin(base_url, link)
    return link

def get_src(text, base_url=None):
    """get_src(text[, base_url]) -> href or None

    Extract the SRC out of an IMG tag.  If base_url is provided,
    will attempt to resolve relative links.

    """
    # get the value of SRC in an IMG tag, or None if not found
    m = re.search(r'src\s*=\s*"?([^>" ]+)"?', text, re.I)
    if not m:
        return None
    link = m.group(1)
    if base_url and not link.lower().startswith("http"):
        import urlparse
        link = urlparse.urljoin(base_url, link)
    return link

def resolve_relative_links(html, base_url):
    """Return the html with all the relative links made into absolute
    links, using base_url as the source of the page."""
    import sgmllib
    import urlparse
    import StringIO

    class MyParser(sgmllib.SGMLParser):
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self.outhandle = StringIO.StringIO()
        def unknown_starttag(self, tag, attributes):
            contents = [tag]
            for name, value in attributes:
                if name == 'href':
                    value = urlparse.urljoin(base_url, value)
                contents.append("%s=%s" % (name, value))
            self.outhandle.write("<%s>" % ' '.join(contents))
        def unknown_endtag(self, tag):
            self.outhandle.write("</%s>" % tag)
        def handle_data(self, data):
            self.outhandle.write(data)
    parser = MyParser()
    parser.feed(html)
    parser.outhandle.seek(0)
    return parser.outhandle.read()
    
if __name__ == '__main__':
    test_parse_ranges()
