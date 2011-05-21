"""

Classes:
Canvas    Draw stuff on here

Functions:
plot_matrix

"""
import sys, os

import SVGdraw

def draw_svg(matrix, gene_names, array_names, boxwidth, boxheight):
    import svgplot
    canvas = svgplot.Canvas()
    svgplot.plot_matrix(canvas, matrix, 20, 20,
                        boxwidth, boxheight,
                        row_names=gene_names, row_name_width=100,
                        col_names=array_names, col_name_height=100)
    return canvas.format()


SVG_ELEMENTS = {
    "rect" : (SVGdraw.rect,
              ("x", "y", "width", "height"),
              ("fill", "stroke", "stroke_width"),
              ),

    "text" : (SVGdraw.text,
              ("x", "y", "text"),
              ("font_size", "font_family", "fill"),
              ),
    
    "group" : (SVGdraw.group,
              (),
              ("transform",),
              ),
    }

def make(name, *args, **keywds):
    if name not in SVG_ELEMENTS:
        raise ValueError, "I don't know how to deal with %s" % name
    fn, reqd_args, opt_args = SVG_ELEMENTS[name]
    all_args = reqd_args + opt_args
    assert len(args) <= len(all_args)
    params = {}
    for i in range(len(args)):
        params[all_args[i]] = args[i]
    params.update(keywds)
    for x in params:
        params[x] = str(params[x])
    return fn(**params)

class Canvas:
    def __init__(self, desc=None):
        self.sd = SVGdraw.svg()

    def add(self, name, *args, **keywds):
        elem = make(name, *args, **keywds)
        self.sd.addElement(elem)

    def addat(self, index, name, *args, **keywds):
        elem = make(name, *args, **keywds)
        self.sd.elements.insert(index, elem)

    def add_element(self, elem):
        self.sd.addElement(elem)

    def format(self):
        d = SVGdraw.drawing()
        d.setSVG(self.sd)
        return d.toXml()

def format_color(r, g, b):
    return "#%02x%02x%02x" % (r, g, b)

def color_red_green(v):
    # Assume values from -1 to 1.
    v = (v+1)/2.        # scale to 0-1
    c = int(256 * 2 * abs(v - 0.5))
    c = int(c * .75)    # make it a little darkerr
    if v >= 0.5:
        r, g, b = 0, c, 0
    else:
        r, g, b = c, 0, 0
    return format_color(r, g, b)

class Unit:
    def __init__(self, value, unit=None):
        self.value = value
        self.unit = unit or ""
    def __add__(self, v):
        unit = self.unit
        if isinstance(v, Unit):
            assert not self.unit or not v.unit or self.unit == v.unit
            if not unit:
                unit = v.unit
            v = getattr(v, "value", v)
        return Unit(self.value+v, unit)
    def __sub__(self, v):
        return self.__add__(-v)
    def __neg__(self):
        return Unit(-self.value, self.unit)
    def __str__(self):
        if type(self.value) is type(0):
            return "%d%s" % (self.value, self.unit)
        return "%g%s" % (self.value, self.unit)

def as_unit(x):
    try:
        x = float(x)
    except ValueError:
        pass
    if type(x) is type(0.0):
        v, u = x, ""
    else:
        v, u = float(x[:-2]), x[-2:]
    assert u in ["", "em", "ex", "px", "pt", "pc", "cm", "mm", "in"]
    return Unit(v, u)

def plot_matrix(canvas, data, x, y, boxwidth, boxheight,
                row_names=None, row_name_width=None,
                col_names=None, col_name_height=None,
                textcolor=None, data2color_fn=None):
    # Initialize all the coordinates.
    orig_x, orig_y = as_unit(x), as_unit(y)
    boxwidth, boxheight = as_unit(boxwidth), as_unit(boxheight)
    row_name_width = row_name_width or 0
    col_name_height = col_name_height or 0

    # Check input data.
    assert len(data)
    x = map(len, data)
    assert min(x) == max(x)
    num_rows, num_cols = len(data), len(data[0])
    assert not row_names or len(row_names) == num_rows
    assert not col_names or len(col_names) == num_cols
    # Choose a good default
    data2color_fn = data2color_fn or color_red_green
    textcolor = textcolor or "black"

    # Draw the row names.
    if row_names and row_name_width:
        x, y = orig_x+row_name_width-2, orig_y+col_name_height+boxheight
        for name in (row_names or []):
            canvas.add("text", x, y, name, **{"text-anchor" : "end"})
            y += boxheight

    # Draw the col names.
    # The text is going to be rotated 90, so the axes will be different.
    if col_names and col_name_height:
        x, y = orig_y+col_name_height-2, -orig_x-row_name_width
        for name in (col_names or []):
            canvas.add("text", x, y, name, transform="rotate(90)",
                       fill=textcolor,
                       **{"text-anchor" : "end"})
            y -= boxwidth

    # Draw the actual matrix.
    matrix_x, matrix_y = orig_x+row_name_width, orig_y+col_name_height
    g = make("group", transform="translate(%s, %s)" % (matrix_x, matrix_y))
    y = Unit(0, y.unit)
    for i in range(len(data)):
        x = Unit(0, x.unit)
        for j in range(len(data[i])):
            d = data[i][j]
            c = data2color_fn(d)
            g.addElement(make("rect", x, y, boxwidth, boxheight, fill=c))
            x += boxwidth
        y += boxheight
    canvas.add_element(g)
