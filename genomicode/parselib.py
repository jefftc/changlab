"""
Parse strings that come from the user.  Also does some formatting.

Functions:
parse_ranges

pretty_int
pretty_float
pretty_ordinal
pretty_range
pretty_pvalue
pretty_date
pretty_list
pretty_time_delta
pretty_filesize

linesplit              Split one long string into lines.
print_split


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
    # Put commas every 3 digits.
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

def pretty_float(num, ndigits=None):
    import math
    
    left = int(math.floor(num))
    right = num - left
    assert right < 1 and right >= 0
    
    left = pretty_int(left)
    if ndigits:
        right = "%.*f" % (ndigits, right)
    else:
        right = str(right)
    assert right.startswith("0.")
    right = right[2:]
    
    x = "%s.%s" % (left, right)
    return x

def pretty_ordinal(num):
    # Negative numbers and 0 not handled.
    assert type(num) is type(0)
    assert num >= 1
    suffix = "th"
    if num == 1:
        suffix = "st"
    elif num == 2:
        suffix = "nd"
    elif num == 3:
        suffix = "rd"
    return "%d%s" % (num, suffix)
    

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

def pretty_date(ctime=None, format=None):
    # ctime should be seconds since epoch or struct_time.
    import time

    ctime = ctime or time.time()
    # e.g. Monday, 06 June 2011, 01:09 PM
    format = format or "%A, %d %B %Y, %I:%M %p"
    
    if type(ctime) is type(0.0):
        ctime = time.localtime(ctime)
    time_str = time.strftime(format, ctime)
    return time_str


def pretty_list(items, max_items=None, conjunction="and"):
    # conjunction is typically "and" or "or".
    assert max_items is None or max_items > 2
    assert type(conjunction) is type("")
    
    assert len(items) >= 0
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    if len(items) == 2:
        return "%s %s %s" % (items[0], conjunction, items[1])
    if max_items is None or len(items) <= max_items:
        x = ", ".join(items[:-1])
        return "%s, %s %s" % (x, conjunction, items[-1])
    x = items[:max_items] + ["..."]
    x = ", ".join(x)
    return x


def pretty_time_delta(delta):
    # Delta is difference in time in seconds.
    assert delta >= 0
    days, x = divmod(delta, 60*60*24)
    hours, x = divmod(x, 60*60)
    minutes, seconds = divmod(x, 60)

    if not days and not hours and not minutes and seconds < 1:
        return "instant"
    if not days and not hours and not minutes:
        if seconds == 1:
            return "1 sec"
        return "%d secs" % seconds
    if not days and not hours:
        x = minutes + seconds/60.
        return "%.1f mins" % x
    if not days:
        x = hours + minutes/60. + seconds/3600.
        return "%.1f hrs" % x

    day_or_days = "day"
    if days > 1:
        day_or_days = "days"
    x = hours + minutes/60. + seconds/3600.
    return "%d %s, %.1f hours" % (days, day_or_days, x)


def pretty_filesize(size):
    # size is the size in bytes.
    assert size >= 0
    bytes = size
    kbytes, bytes = divmod(bytes, 1024)
    mbytes, kbytes = divmod(kbytes, 1024)
    gbytes, mbytes = divmod(mbytes, 1024)
    tbytes, gbytes = divmod(gbytes, 1024)
    if tbytes:
        return "%.2f Tb" % (tbytes + gbytes/1024.)
    if gbytes:
        return "%.2f Gb" % (gbytes + mbytes/1024.)
    if mbytes:
        return "%.2f Mb" % (mbytes + kbytes/1024.)
    if kbytes:
        return "%.2f kb" % (kbytes + bytes/1024.)
    return "%d b" % bytes


def linesplit(one_long_line, prefix1=0, prefixn=4, width=72):
    # Return list of lines.
    lines = one_long_line.split("\n")
    all_lines = []
    for i in range(len(lines)):
        p1, pn = prefix1, prefixn
        if i > 0:
            p1 = pn
        x = _linesplit_h(lines[i], p1, pn, width)
        all_lines.extend(x)
    return all_lines

def _linesplit_h(one_long_line, prefix1, prefixn, width):
    # prefix1 and prefixn can be integer indicating the number of
    # spaces, or a string indicating the prefix to print.
    assert width > 0
    if type(prefix1) is type(0):
        assert prefix1 >= 0
        prefix1 = " " * prefix1
    if type(prefixn) is type(0):
        assert prefixn >= 0
        prefixn = " " * prefixn
    assert width > len(prefix1)
    assert width > len(prefixn)

    assert "\n" not in one_long_line
    assert "\r" not in one_long_line
    assert "\t" not in one_long_line
    
    lines = []
    while 1:
        #ind = " "*indent1
        ind = prefix1
        if lines:
            #ind = " "*indento
            ind = prefixn
        #if ind:
        #    one_long_line = one_long_line.lstrip()  # no leading spaces
        if lines:
            one_long_line = one_long_line.lstrip()  # no leading spaces
        one_long_line = ind + one_long_line

        if len(one_long_line) < width:
            lines.append(one_long_line)
            break

        # Try to split on a space.
        w = width
        i = one_long_line.rfind(" ", len(ind), w)
        if i > 0:
            w = i
        x = one_long_line[:w]
        one_long_line = one_long_line[w:]
        lines.append(x)
    return lines


def print_split(one_long_line, prefix1=0, prefixn=4, width=72, outhandle=None):
    import sys
    outhandle = outhandle or sys.stdout

    lines = linesplit(
        one_long_line, prefix1=prefix1, prefixn=prefixn, width=width)
    for line in lines:
        print >>outhandle, line

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
