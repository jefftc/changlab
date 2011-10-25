"""

povray

include
declare
local
global_settings
background
camera
light_source

Objects:
object_
sphere
cylinder
box
text

Transformations:
translate
translate_by_extent
rotate
scale

Object modifiers:
texture
pigment
finish

Simple items:
ambient_light    for global_settings
projection       for camera
location
up
right
angle
sky
look_at

light_type       for light_source
area_args
radius
falloff
tightness
point_at

ambient          for finish
diffuse
brilliance
crand
phong
phong_size
specular
roughness
reflection
refraction
ior

no_shadow

blob_threshold
blob_sphere
blob_cylinder

Miscellaneous:
color
vector

"""
# _fmt_simple
# _fmt_complex
#
# _fmt_nolabel_fn
# _fmt_only_label_fn
# _fmt_amount_fn
# _fmt_color_fn
# _fmt_vector_fn

def _fmt_simple(name, arg):
    """Format a simple item that consists of a name and argument."""
    if arg is None:
        return None
    return "%s %s" % (name, arg)

def _fmt_complex(name, *items, **params):
    """Format a complex item.  A complex item takes a set of arguments
    as well as nested items.

    *items contains the items to include in this item.  They can be
    simple items or other complex items.

    "args" is a keyword parameter that contains the arguments for this
    complex item.  They are unlabeled and rendered as a
    comma-separated at the top of the item.

    """
    # If there are arguments, then format them and add them to the
    # list of items.
    if "args" in params:
        x = ", ".join(map(str, params["args"]))
        items = [x] + list(items)

    # Do the formatting.
    items = [x for x in items if x]
    if not items:
        return ""
    if len(items) == 1:
        return "%s { %s }" % (name, items[0])
    s = "%s {\n" % name
    for x in items:
        x = str(x).strip()
        x = x.replace("\n", "\n  ")  # indent each line in x.
        x = "  %s\n" % x
        s += x
    s += "}\n"
    return s

class _fmt_nolabel_fn:
    def __init__(self, name):
        self.name = name
    def __call__(self, value):
        return value

class _fmt_only_label_fn:
    def __init__(self, name):
        self.name = name
    def __call__(self):
        return self.name

class _fmt_amount_fn:
    def __init__(self, name):
        self.name = name
    def __call__(self, amount):
        return _fmt_simple(self.name, amount)

class _fmt_color_fn:
    def __init__(self, name):
        self.name = name
    def __call__(self, *args, **keywds):
        return _fmt_simple(self.name, color(*args, **keywds))

class _fmt_vector_fn:
    def __init__(self, name):
        self.name = name
    def __call__(self, *args):
        return _fmt_simple(self.name, vector(*args))

def povray(
    filename, outfile=None, height=None, width=None, antialias=None,
    quality=None, povray_bin=None):
    # Return a handle to the results.
    # antialias  Which colors to anti-alias (0-3.0).
    # quality    0-9.
    import os
    import subprocess
    
    assert os.path.exists(filename)
    povray_bin = povray_bin or "povray"
    
    args = []
    args.append("-D")   # don't show output.
    args.append("-J")   # turn off jitter.
    if outfile:
        args.append("+O%s" % outfile)
    if height is not None:
        args.append("-H%d" % height)
    if width is not None:
        args.append("-W%d" % width)
    if antialias is not None:
        assert antialias >= 0 and antialias <= 3
        args.append("+A%g" % antialias)
    if quality is not None:
        assert quality >= 0 and quality <= 9
        args.append("-Q%d" % quality)
    args.append(filename)

    cmd = "%s %s" % (povray_bin, " ".join(args))
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    w, r = p.stdin, p.stdout
    #w, r = os.popen4(cmd)
    w.close()
    return r

def include(file_inc):
    return '#include "%s"' % file_inc

def declare(name, value):
    return "#declare %s = %s" % (name, value)

def local(name, value):
    return "#local %s = %s;" % (name, value)

def global_settings(*items):
    # ambient_light COLOR    def "rgb <1,1,1>"
    return _fmt_complex("global_settings", *items)

def background(color):
    return _fmt_complex("background", color)

def camera(*items):
    # projection        perspective, orthographic, fisheye, etc...
    # location VECTOR   coords of the camera.
    # up VECTOR
    # right VECTOR
    # angle AMOUNT      viewing angle of camera used.
    # rotate VECTOR
    # translate VECTOR
    # sky VECTOR        direction of sky (def <0, 1, 0>).
    # look_at VECTOR    where the camera is pointed.
    #
    # Notes:
    # o Projections:
    #   perspective   classic pinhole camera
    #   orthographic  uses parallel camera rays
    # o up and right tell POV-Ray relative height (up) and width
    #   (right) of the view screen.  By default:
    #     right 4/3*x
    #     up y
    #   In the default perspective camera, also defines initial plane
    #   of view screen before moving it with the look_at or rotate
    #   vectors.  In the orthographic projection, lengths of up and
    #   right vectors set the size of the viewing window regardless of
    #   the direction vector.
    # o look_at should almost always be last item in camera statement.
    #    If other items placed after this, then camera may not
    #    continue to look at the specified point.
    return _fmt_complex("camera", *items)

def light_source(location, color, *items):
    # location VECTOR
    # color COLOR       2,2,2 will be 2x as bright as 1,1,1.
    # light_type VALUE  spotlight, shadowless, cylinder, area_light
    # light_type VALUE  parallel
    # area_args ARGS    area_light.  REQUIRED.
    # radius AMOUNT     spotlight, cylinder  (def 30, 0.75)
    # falloff AMOUNT    spotlight, cylinder  (def 45, 1)
    # tightness AMOUNT  spotlight, cylinder  (def 0, 0)
    # point_at VECTOR   spotlight, parallel  (def <0, 0, 0>)
    # see docs for more...
    #
    # NOTES:
    # o parallel can be used with any type of light source.  Makes
    #   light rays parallel to simulate very distant light sources.
    # o For normal point lights, point_at must come after parallel.
    return _fmt_complex("light_source", args=(location, color), *items)

def object_(name, *items):
    x = _fmt_complex("object", args=(name,), *items)
    if not x.endswith("\n"):
        x += "\n"
    return x

def sphere(center, radius, *items):
    # center VECTOR
    # radius AMOUNT
    # + modifiers
    return _fmt_complex("sphere", args=(center, radius), *items)

def blob_sphere(center, radius, strength, *items):
    # center VECTOR
    # radius AMOUNT
    # strength AMOUNT
    # + modifiers
    return _fmt_complex("sphere", args=(center, radius, strength), *items)

def superellipsoid(value_e, value_n, *items):
    # http://local.wasp.uwa.edu.au/~pbourke/geometry/superellipse/
    # http://www.f-lohmueller.de/pov_tut/all_shapes/shapes350e.htm
    x = "<%s, %s>" % (value_e, value_n)
    return _fmt_complex("superellipsoid", args=(x,), *items)

def cylinder(base_point, cap_point, radius, *items):
    # base_point VECTOR
    # cap_point VECTOR
    # radius AMOUNT
    # + modifiers
    return _fmt_complex(
        "cylinder", args=(base_point, cap_point, radius), *items)

def blob_cylinder(base_point, cap_point, radius, strength, *items):
    # base_point VECTOR
    # cap_point VECTOR
    # radius AMOUNT
    # strength AMOUNT
    # + modifiers
    return _fmt_complex(
        "cylinder", args=(base_point, cap_point, radius, strength), *items)

def box(corner_1, corner_2, *items):
    # corner_1 VECTOR
    # corner_2 VECTOR
    # + modifers
    return _fmt_complex("box", args=(corner_1, corner_2), *items)

def text(file, string, thickness, offset, *items):
    # If "file" is the name of a file, then you need to enclose it
    # within quotes.  Otherwise, it will be interpreted as a variable
    # name.  thickness is best from 0.5 to 2.  Offset will add kerning
    # distance.
    # Text starts with origin in the lower left.  Thickness goes in
    # the +z direction.
    x = 'ttf %s "%s" %g, %g' % (file, string, thickness, offset)
    items = (x,) + items
    return _fmt_complex("text", *items)

translate = _fmt_vector_fn("translate")

def translate_by_extent(name, x_trans, y_trans, z_trans):
    from StringIO import StringIO
    diff_name = "%s_DEL" % name
    handle = StringIO()
    w = handle.write
    w(local(diff_name, "max_extent(%s)-min_extent(%s)" % (name, name))+"\n")
    x = y = z = "0"
    if x_trans:
        x = "%s.x*%g" % (diff_name, x_trans)
    if y_trans:
        y = "%s.y*%g" % (diff_name, y_trans)
    if z_trans:
        z = "%s.z*%g" % (diff_name, z_trans)
    w(translate(x, y, z)+"\n")
    handle.seek(0)
    return handle.read()

# Rotates around the origin by x-, then y-, then z- degrees.
# Right-handed coordinate system, so positive z is counter-clockwise.
# Most of the time should rotate before translate.
rotate = _fmt_vector_fn("rotate")

# Should translate after scale.  Scaling after translate will multiply
# the translate amount, e.g. translate <5, 6, 7> scale 4 will
# translate to <20, 24, 28>.
scale = _fmt_vector_fn("scale")

def texture(*items):
    # pigment
    # normal
    # finish
    return _fmt_complex("texture", *items)

def pigment(color):
    return _fmt_complex("pigment", color)

def finish(*items):
    # ambient COLOR     def <0.1, 0.1, 0.1>
    # diffuse VALUE
    # brilliance VALUE
    # crand VALUE       OBSOLETE
    # phong VALUE
    # phong_size VALUE
    # specular VALUE
    # roughness VALUE
    # metallic
    # reflection VALUE
    # refraction VALUE    NOT NECESSARY ANYMORE
    # ior        VALUE    OBSOLETE.  NOW IN INTERIOR.
    # more...
    #
    # Notes:
    # o ambient controls the amount of ambient light that is produced
    #   by the object.  0.0 means object will be black if not directly
    #   lighted.  Default 0.1.
    # o diffuse controls how much the reflected light is scattered in
    #   a variety of directions.  Only very smooth polished surfaces
    #   have perfect specular (not diffuse) reflections.  Range from
    #   0 (dull) - 1 (shiny), default is 0.6.
    # o brilliance varies the amount that light falls off depending on
    #   the angle of incidence.  Higher values look more metallic,
    #   bolder.  Default is 1.0, and 5-10 causes light to fall off
    #   less at medium to low angles.
    # o crand causes minor random darkening, e.g. for very rough
    #   surfaces such as concrete or sand.  Default is 0, and
    #   typically ranges from 0.01 to 0.5.
    # o phong causes bright and shiny spots that are the color of the
    #   light source being reflected.  Typically from 0.0-1.0.  Bigger
    #   means brighter shiny spot.
    # o phong_size controls the size of the highlight spot.  Larger
    #   values generate tighter, smaller, shinier highlights.  Values
    #   range from 1 (very dull) to 250 (highly polished).  Default is
    #   40 (plastic).
    # o specular is similar to phong, but uses a different model.
    #   Values range from 0-1.
    # o roughness determines size of spot for specular.  Values range
    #   from 0.0005 (very smooth, small highlight) to 1.0 (very rough,
    #   large highlight).  Default is 0.05 (plastic).
    # o metallic can be used with phong or specular to model a
    #   metallic surface.
    # o reflection indicates how reflective the object is.  Ranges
    #   from 0-1.  Default 0.15.
    return _fmt_complex("finish", *items)

def interior(*items):
    # ior        VALUE
    #
    # Notes:
    # o refraction is either 0 (disable) or 1 (enable).  For
    #   transparent objects, will refract the light going through it.
    # o ior is the index of refraction.  Empty space is 1.0 (default).
    #   1.000292 is air, 1.33 water, and 1.5 glass.
    return _fmt_complex("interior", *items)

ambient_light = _fmt_color_fn("ambient_light")

projection = _fmt_nolabel_fn("projection")
location = _fmt_vector_fn("location")
up = _fmt_vector_fn("up")
right = _fmt_vector_fn("right")
angle = _fmt_amount_fn("angle")
sky = _fmt_vector_fn("sky")
look_at = _fmt_vector_fn("look_at")

light_type = _fmt_nolabel_fn("light_type")

def area_args(axis_1, axis_2, size_1, size_2):
    # axis_1 VECTOR  Defines the location of the area light.
    # axis_2 VECTOR
    # size_1 AMOUNT  Number of rows of lights.
    # size_2 AMOUNT  Number of columns of lights.
    #
    # NOTES:
    # o The location vector of the light source is the center of the
    #   rectangle created by the two axis vectors.
    x = axis_1, axis_2, size_1, size_2
    x = ", ".join(map(str, x))
    return x

radius = _fmt_amount_fn("radius")
falloff = _fmt_amount_fn("falloff")
tightness = _fmt_amount_fn("tightness")
point_at = _fmt_vector_fn("point_at")

ambient = _fmt_color_fn("ambient")
diffuse = _fmt_amount_fn("diffuse")
brilliance = _fmt_amount_fn("brilliance")
crand = _fmt_amount_fn("crand")
phong = _fmt_amount_fn("phong")
phong_size = _fmt_amount_fn("phong_size")
specular = _fmt_amount_fn("specular")
roughness = _fmt_amount_fn("roughness")
metallic = _fmt_only_label_fn("metallic")
reflection = _fmt_amount_fn("reflection")
refraction = _fmt_amount_fn("refraction")
ior = _fmt_amount_fn("ior")

no_shadow = _fmt_only_label_fn("no_shadow")

blob_threshold = _fmt_amount_fn("threshold")
    
def color(red, green=None, blue=None, filter=None, transmit=None):
    # filter    AMOUNT
    # transmit  AMOUNT  
    #
    # Notes:
    # o For filter, the light passing through is tinted by the
    #   appropriate color (e.g. stained glass window).
    #   1, 1, 1, filter=1 makes a clear object.
    # o For transmit, all frequencies of light pass through tiny holes
    #   in the surface (e.g. mesh netting).

    if green is None and blue is None:
        green = blue = red
    
    # Don't display integers as floating point, e.g. 1.0 -> 1.
    if abs(red-round(red, 0))<1E-3:
        red = int(round(red, 0))
    if abs(green-round(green, 0))<1E-3:
        green = int(round(green, 0))
    if abs(blue-round(blue, 0))<1E-3:
        blue = int(round(blue, 0))
        
    fn = "rgb"
    if filter is not None:
        fn += "f"
    if transmit is not None:
        fn += "t"
    if fn == "rgb" and red == green == blue:
        # Use simplified expression.
        return "rgb %g" % red

    x = [red, green, blue, filter, transmit]
    x = [x for x in x if x is not None]
    x = map(str, x)
    return "%s <%s>" % (fn, ", ".join(x))

def vector(x, y, z):
    # Need to be %s because vectors can contain expressions instead of
    # numbers.
    if x == y == z:
        return "%s" % x
    if x and not y and not z:
        if type(x) is type(""):
            x = '(%s)' % x
        return "x*%s" % x
    if y and not x and not z:
        if type(y) is type(""):
            y = '(%s)' % y
        return "y*%s" % y
    if z and not x and not y:
        if type(z) is type(""):
            z = '(%s)' % z
        return "z*%s" % z
    return "<%s, %s, %s>" % (x, y, z)
