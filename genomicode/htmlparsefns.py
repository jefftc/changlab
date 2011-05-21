"""Miscellaneous functions to help parse HTML pages.

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

"""
# NOTE: The functions in this package use quick and dirty regular
# expressions to handle HTML.  There are cases where the code can get
# fooled by clever or obscure patterns.  For example, it looks for
# HTML tags using the expresion "<[^>]*>", which will fail if there
# are special characters, e.g. <A HREF="foo>bar">.  There are ways to
# handle this better (e.g. using a real HTML parser), but nearly all
# of them result in more code complexity or performance loss.

import re

def _make_tag_patterns(tag):
    # Return patterns for "<tag>" and "</tag>"
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
    get_tags_and_contents("<B>Hi!</B><B>Bye</B>", "b") ->
      ["<B>Hi!</B>", "<B>Bye</B>"]

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
