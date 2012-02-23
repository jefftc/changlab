"""

Uses the 4th quadrant in the Cartesian coordinate system. (0, 0, 0) is
at the upper left corner of the screen, extending positive to the
right, down, and toward the viewer.
    
Color is tuple of (R, G, B) where RGB is from 0.0 to 1.0.


Functions:
image

box
sphere
cylinder
text
text90
get_text_size

write

"""
# _format_appearance
# _coord2pr
# _position_camera
# _declare_fontfile
# _draw_background

class Image:
    def __init__(self, width, height, depth):
        import StringIO
        self.width = width
        self.height = height
        self.depth = depth
        self.handle = StringIO.StringIO()
        self._object_id = 0
        self._objects = []    # list of graph objects.  subset of all objects
    def _make_object_name(self):
        name = "OBJECT%05d" % self._object_id
        self._object_id += 1
        self._objects.append(name)
        return name
    def format(self, outhandle):
        if type(outhandle) is type(""):
            outhandle = open(outhandle, 'w')
        outhandle.write(self.handle.getvalue())

def image(width, height, depth, make_3d):
    im = Image(width, height, depth)
    _declare_fontfile(im)
    _position_camera(im, width, height, depth, make_3d)
    #_declare_points(im)
    _draw_background(im)
    # Clear the canvas for the client to draw on.
    del im._objects[:]
    return im

def box(image, coord, extent, color, shadow=False, finish=None):
    # coord is (x, y, z).  extent is (width, height, depth).
    import povray as pr
    x = _format_appearance(color, finish)
    finish, pigment, interior = x
    w = image.handle.write
    coord1 = _coord2pr(image, coord)
    x, y, z = coord1
    width, height, depth = extent
    coord2 = x+width, y-height, z-depth
    name = image._make_object_name()
    #print coord1, coord2
    ns = None
    if not shadow:
        ns = pr.no_shadow()
    w(pr.declare(name, pr.box(
        pr.vector(*coord1), pr.vector(*coord2),
        pigment,
        finish,
        interior,
        ns,
        )))
    w(pr.object_(name))
    w("\n")
    
    #extent = coord2[0]-coord1[0], coord2[1]-coord1[1], coord2[2]-coord1[2]
    #image._add_object(name, coord1, extent)

        
def sphere(image, coord, radius, color, shadow=False, shape=None,
           finish=None):
    # coord is the center of the sphere, with radius.
    import povray as pr
    import graphconst as gc

    x = _format_appearance(color, finish)
    finish, pigment, interior = x
    coord = _coord2pr(image, coord)

    ns = None
    if not shadow:
        ns = pr.no_shadow()
        
    # DEFAULT, CIRCLE, SQUARE, DIAMOND
    if shape is None:
        shape = gc.DEFAULT
    if shape in [gc.DEFAULT, gc.CIRCLE]:
        point = pr.superellipsoid(
            1.0, 0.50, pr.scale(radius, radius, radius),
            pr.translate(*coord), pigment, ns, finish, interior)
    # Square
    elif shape == gc.SQUARE:
        point = pr.superellipsoid(
            0.25, 0.25, pr.scale(radius, radius, radius),
            pr.translate(*coord), pigment, ns, finish, interior)
    # Diamond
    elif shape == gc.DIAMOND:
        point = pr.superellipsoid(
            1.75, 0.50, pr.scale(radius, radius, radius),
            pr.translate(*coord), pigment, ns, finish, interior)
    else:
        raise AssertionError, "Unknown point: %s" % (shape)
    
    w = image.handle.write
    name = image._make_object_name()
    w(pr.declare(name, point))
    w(pr.object_(name))
    w("\n")

    #x, y, z = coord
    #coord = x-radius, y-radius, z-radius
    #extent = radius*2, radius*2, radius*2
    #image._add_object(name, coord, extent)

def cylinder(
    image, coord, extent, radius, color, shadow=False, finish=None, blob=None):
    # coord is the center of one end, with the center of the other end
    # in the direction determined by extent.  blob can be True or a
    # blob strength.
    import povray as pr
    x = _format_appearance(color, finish)
    finish, pigment, interior = x
    if blob is True:
        blob = 1.0
        
    w = image.handle.write
    coord1 = _coord2pr(image, coord)
    x, y, z = coord1
    width, height, depth = extent
    coord2 = x+width, y-height, z-depth
    name = image._make_object_name()
    ns = None
    if not shadow:
        ns = pr.no_shadow()
    w(pr.declare(name, pr.cylinder(
        pr.vector(*coord1), pr.vector(*coord2), radius,
        pigment, finish, interior, ns, 
        )))
    w(pr.object_(name))

    # HACK: min_extent and max_extent doesn't work for blobs for some
    # reason.  Solution is to draw the object so that the extent can
    # be calculated, and then draw a blob at the same place.
    if blob:
        w(pr._fmt_complex(
            "blob", pr.blob_threshold(0.01), pr.blob_cylinder(
            pr.vector(*coord1), pr.vector(*coord2), radius, blob),
            pigment, finish, interior, ns, 
            ))
    w("\n")

def text(
    image, text, coord, depth, fontsize, color, shadow=False,
    center_x=False, center_y=False, center_z=False,
    wrong_x=False, wrong_y=False, wrong_z=False,
    min_x=False, max_x=False, min_y=False, max_y=False,
    background=None, vertical=False, finish=None):
    # The coordinate should be at the upper left point of the text.
    # depth grows in the +z direction (toward the user).
    # center_x means the x-coord is at the middle, and wrong_x means
    # the x-coord is at the right.  They should not both be True.
    import operator
    import povray as pr
    import graphconst as gc
    
    assert not (center_x and wrong_x)
    assert not (center_y and wrong_y)
    assert not (center_z and wrong_z)
    assert not (min_x and max_x)
    assert not (min_y and max_y)
    # Background should be a true value, or the color of the background.
    if background and not operator.isSequenceType(background):
        background = (1, 1, 1)
    assert not background or len(background) == 3

    x = _format_appearance(color, finish)
    pr_finish, pr_pigment, pr_interior = x

    w = image.handle.write

    # Fix the coordinates.  If min_x, max_x, min_y, or max_y are
    # given, then x or y are relative coordinates.
    #old_coord = coord
    pr_coord = _coord2pr(image, coord)
    x, y, z = coord
    pr_x, pr_y, pr_z = pr_coord
    if min_x or max_x:
        pr_x = x
    if min_y or max_y:
        pr_y = -y
    pr_z = pr_z-depth   # povray z is flipped
    coord = pr_x, pr_y, pr_z

    # Declare the text object and apply translations.
    rotate = pr.rotate(0, 0, 0)
    if vertical:
        rotate = pr.rotate(0, 0, 90)
        
    name = image._make_object_name()
    ns = None
    if not shadow:
        ns = pr.no_shadow()
    w(pr.declare(name, pr.text(
        "FONTFILE", text, depth, 0,
        pr_pigment,
        pr_finish,
        pr_interior,
        ns,
        rotate,
        pr.scale(fontsize, fontsize, 1),
        pr.translate(*coord),
        )))
    
    # Translate by min or max of everything drawn.
    x_trans = y_trans = 0
    names = [x for x in image._objects if x != name]
    if min_x and image._objects:
        x = ["min_extent(%s).x" % x for x in names]
        x = "min(\n  %s)" % ",\n  ".join(x)
        x_trans = x
    elif max_x and image._objects:
        x = ["max_extent(%s).x" % x for x in names]
        x = "max(\n  %s)" % ",\n  ".join(x)
        x_trans = x
    if min_y and image._objects:
        # Convert coordinates to POV-RAY.
        x = ["max_extent(%s).y" % x for x in names]
        x = "max(\n  %s)" % ",\n  ".join(x)
        y_trans = x
    elif max_y and image._objects:
        # Hack: max_y can be an integer that indicates the number of
        # last objects to ignore.
        n = names
        if type(max_y) is type(0):
            max_y = min(max_y, len(names))
            n = names[:-max_y]
        if n:
            x = ["min_extent(%s).y" % x for x in n]
            x = "min(\n  %s)" % ",\n  ".join(x)
            y_trans = x
    if x_trans or y_trans:
        w(pr.declare(name, pr.object_(
            name, pr.translate(x_trans, y_trans, 0))))
        
    # Translate by extent.
    x_extent = 0
    y_extent = -1      # povray coordinates are from the lower left
    z_extent = 0
    if vertical:
        x_extent += 1  # Rotating will throw off the coordinates.
    if center_x:
        x_extent += -0.5
    elif wrong_x:
        x_extent += -1.0
    if center_y:
        y_extent += 0.5
    elif wrong_y:
        y_extent += 1.0
    if center_z:
        z_extent += 0.5
    elif wrong_z:
        z_extent += 1.0
    if x_extent or y_extent or z_extent:
        w(pr.declare(name, pr.object_(
            name, pr.translate_by_extent(name, x_extent, y_extent, z_extent))))
    w(pr.object_(name))

    # Coordinates are hard to calculate.  Hope this isn't necessary.
    #image._add_object(name, None, None)
    
    if not background:
        return
    #background = (1, 0, 0)
    x = _format_appearance(background, gc.TRANSPARENT)
    pr_finish, pr_pigment, pr_interior = x
    BACKGROUND_DEPTH = 0.2*depth
    coord1 = (
        "min_extent(%s).x" % name, "min_extent(%s).y" % name, pr_z)
    coord2 = (
        "max_extent(%s).x" % name, "max_extent(%s).y" % name,
        pr_z-BACKGROUND_DEPTH)
    w(pr.box(
        pr.vector(*coord1), pr.vector(*coord2),
        pr_pigment, pr_finish, pr_interior, pr.no_shadow()))
    w("\n")
    
    
def text90(*args, **keywds):
    keywds["vertical"] = True
    return text(*args, **keywds)

def get_text_size(text, fontsize):
    # This only works for Verdana Bold!.  Also, doesn't get the exact
    # size.  Just an estimate.

    REF_FONTSIZE = 100
    #DEFAULT_HEIGHT = 100  # includes drop letters
    DEFAULT_HEIGHT = 118  # includes drop letters
    # On average, uses 9.75 pixels between letters (see 110706).
    SPACING = 9.75

    # Sizes of letters in Verdana Bold at 100 fontsize.  The heights
    # are pretty much useless because some of the letters drop down.
    letter2size = {
        "A" : (75, 73), "B" : (63, 73), "C" : (63, 77), "D" : (69, 73),
        "E" : (53, 73), "F" : (52, 73), "G" : (69, 78), "H" : (66, 73),
        "I" : (43, 73), "J" : (45, 75), "K" : (68, 73), "L" : (53, 73),
        "M" : (77, 73), "N" : (66, 73), "O" : (75, 78), "P" : (60, 73),
        "Q" : (75, 116), "R" : (70, 73), "S" : (63, 77), "T" : (65, 73),
        "U" : (65, 76), "V" : (74, 73), "W" : (107, 73), "X" : (74, 73),
        "Y" : (73, 73), "Z" : (61, 73),
        "a" : (55, 60), "b" : (58, 78), "c" : (52, 60), "d" : (58, 79),
        "e" : (58, 60), "f" : (42, 77), "g" : (58, 98), "h" : (55, 76),
        "i" : (19, 76), "j" : (35, 118), "k" : (58, 76), "l" : (18, 76),
        "m" : (90, 57), "n" : (55, 57), "o" : (60, 60), "p" : (58, 97),
        "q" : (58, 97), "r" : (40, 55), "s" : (52, 60), "t" : (42, 73),
        "u" : (55, 58), "v" : (62, 55), "w" : (94, 55), "x" : (64, 55),
        "y" : (62, 95), "z" : (52, 55),
        "0" : (61, 78), "1" : (49, 73), "2" : (57, 74), "3" : (58, 78),
        "4" : (63, 73), "5" : (57, 76), "6" : (60, 77), "7" : (57, 73),
        "8" : (62, 78), "9" : (60, 77),
        "#" : (72, 73), "*" : (51, 50), "+" : (66, 65), "-" : (38, 20),
        "/" : (54, 108), "=" : (63, 50), "_" : (73, 28),
        }

    # Calculate a default width to use for symbols not in this table.
    total = sum([w for (w, h) in letter2size.itervalues()])
    DEFAULT_WIDTH = float(total)/len(letter2size)

    width = 0
    height = DEFAULT_HEIGHT   # Use a default height of 73.
    for l in text:
        w, h = letter2size.get(str(l), (DEFAULT_WIDTH, DEFAULT_HEIGHT))
        width += w
    width += SPACING * (len(text)-1)

    # Scale the height and width by fontsize.
    width = width * fontsize / float(REF_FONTSIZE)
    height = height * fontsize / float(REF_FONTSIZE)
    
    return width, height

## def line(image, x, y, width, height, color):
##     raise NotImplementedError

def write(image, handle, povray_bin=None):
    # Format the image with POV-RAY and write the results to handle.
    import os
    import tempfile
    import povray

    if type(handle) is type(""):
        handle = open(handle, 'w')

    pov_file = out_file = None
    try:
        x, pov_file = tempfile.mkstemp(dir=".", suffix=".pov"); os.close(x)
        x, out_file = tempfile.mkstemp(dir=".", suffix=".png"); os.close(x)
        if os.path.exists(out_file):
            os.unlink(out_file)

        image.format(open(pov_file, 'w'))
        r = povray.povray(
            pov_file, outfile=out_file, 
            height=image.height, width=image.width, antialias=0.5, quality=9,
            povray_bin=povray_bin)
        output = r.read()
        assert os.path.exists(out_file), "POV-RAY failed.\n%s" % output
        handle.write(open(out_file).read())
    finally:
        if pov_file and os.path.exists(pov_file):
            os.unlink(pov_file)
        if out_file and os.path.exists(out_file):
            os.unlink(out_file)
    output = output.replace("\r", "\n")
    return output

def _format_appearance(color, finish):
    import povray as pr
    import graphconst as gc

    pr_finish = pr_pigment = pr_interior = ""

    pr_pigment = pr.pigment(pr.color(*color))
    
    if finish is None:
        finish = gc.SIMPLE
    if finish == gc.SIMPLE:
        # Ambient of 0.8 and diffuse 0.3 is too dark.  At default
        # camera position, background is gray.
        pr_finish = pr.finish(
            pr.ambient(0.95),
            pr.diffuse(0.4),
            )
    elif finish == gc.METALLIC:
        ## x = pr.finish(
        ##     pr.ambient(0.4, 0.4, 0.4),
        ##     pr.diffuse(0.3),
        ##     pr.phong(0.3),
        ##     pr.phong_size(5),
        ##     )
        #
        # F_MetalA.  Bold.  Probably best one.
        # F_MetalB.  Looks pretty good.
        # F_MetalC.  Colors washed out.  Bright.
        # F_MetalE.  Too washed out.
        pr_finish = pr.finish(
            pr.ambient(0.45),
            #pr.ambient(0.30),
            pr.brilliance(2),
            pr.diffuse(0.3),
            pr.metallic(),
            pr.specular(0.80),
            pr.roughness(1./20),
            pr.reflection(0.1),
            )
    elif finish == gc.ROUGH:
        pr_finish = pr.finish(
            pr.ambient(0.7),
            pr.diffuse(0.3),
            pr.crand(0.1),
            )
    elif finish == gc.SHINY:
        # Doesn't work very well for graphs.
        pr_finish = pr.finish(
            pr.ambient(0.7),
            pr.diffuse(0.1),
            pr.reflection(0.3),
            #pr.ior(1.5),
            )
        pr_interior = pr.interior(pr.ior(1.5))
    elif finish == gc.TRANSPARENT:
        pr_finish = pr.finish(
            #pr.refraction(1),
            #pr.ior(1.3),
            pr.ambient(0.8),
            pr.diffuse(0.3),
            )
        pr_pigment = pr.pigment(pr.color(*color, filter=0.9))
        pr_interior = pr.interior(pr.ior(1.3))
    else:
        raise AssertionError, "Unknown finish: %s" % finish
    return pr_finish, pr_pigment, pr_interior


# Finishes:
# DIFF  CRAND  AMB  REFR  IOR  PHO  PHS  DESCRIPTION        USES
#              1.0   1    1.5            Shines through     onpoint label
#                                        pigment 1, 1, 1, 0.25
# 0.1          0.5                                          Error bar
# 
# ROUGH
# 0.4    0.2   0.5                       Rough, textured    Bar graph
# 
# METALLIC
# 0.6          0.6             0.2    5                     axes
# 0.4          0.3             0.4    5                     Points
# 0.4          0.5             0.3   10  metallic           line
# 0.6          0.6             0.3   10                     Tick
# 0.6          0.6             0.3   10                     axis cap
# 
# SIMPLE
#             0.75                                          background
# 0.3          0.8

def _coord2pr(image, coord):
    # Convert the coordinate system for this library to one that is
    # compatible with POV-RAY.  In POV-RAY, (0, 0, 0) is at the bottom
    # left, extending positive to the right, up, and away from the
    # viewer.  The Y and Z coordinates need to be fixed.
    x, y, z = coord
    y = image.height - y
    z = -z
    return x, y, z

def _position_camera(image, width, height, depth, make_3d):
    import math
    import povray as pr
    
    CAMERA_HEIGHT = depth * 1.2
    LIGHT_ANGLE = 70        # lower means longer shadows, darker colors
    LIGHT_COLOR = 1, 1, 1

    w = image.handle.write

    # Set the camera to the middle of the plot, looking down.
    # Have to look at the middle of the entire plot, or else the
    # borders will be off.
    x_mid, y_mid, z_mid = width*0.5, height*0.5, depth*0.5
    if not make_3d:
        # Default camera location.
        camera = (x_mid, y_mid, -CAMERA_HEIGHT)
        #look = (x_mid, y_mid, z_mid)   # should be -z_mid?
        look = (x_mid, y_mid, 0)
        dist_scale = 1.0
        # Make this a bit dimmer to compensate for the fact that the
        # 3D projection is further away.
        w(pr.global_settings(pr.ambient_light(0.67, 0.67, 0.67)))
        w("\n")
    else:
        camera = (width*0.60, height*0.50, -CAMERA_HEIGHT)
        look = (x_mid*0.25, y_mid*0.35, z_mid)
        dist_scale = 1.5
        #dist_scale = 1.0
        #camera = (width*0.70, height*0.50, -depth*1.2)
        #look = (x_mid*0.25, y_mid*0.35, z_mid)
        #dist_scale = 1.35
    
    w(pr.camera(
        pr.projection("orthographic"),
        pr.location(*camera),
        pr.right(width*dist_scale, 0, 0),
        pr.up(0, height*dist_scale, 0),
        pr.look_at(*look),
        ))
        
    # The light source is at the upper right of the plot, shining
    # toward the lower left.  Calculate the height such that the
    # light is shining at a specific angle relative to the plane
    # of the plot.  Lower LIGHT_ANGLE means light source is closer
    # to the plane (longer shadows).
    light_x = width * 0.95
    light_y = height * 0.95
    target_x = width * 0.15
    target_y = height * 0.15
    side = math.sqrt((light_x-target_x)**2 + (light_y-target_y)**2)
    z = side * math.tan(math.pi*LIGHT_ANGLE/180.0)
    #print z
    w(pr.light_source(
        pr.vector(light_x, light_y, -z),
        pr.color(*LIGHT_COLOR),
        pr.light_type("parallel"),
        ))
    w("\n")

def _declare_fontfile(image, fontfile=None):
    import os
    import povray as pr

    if not fontfile:
        x = os.path.join(os.path.split(__file__)[0], "Verdana Bold.ttf")
        fontfile = os.path.realpath(x)
    fontfile = fontfile or default_font
    assert os.path.exists(fontfile)

    w = image.handle.write
    w(pr.declare("FONTFILE", '"%s"' % fontfile)+"\n")
    w("\n")

## def _declare_points(image):
##     import povray as pr
##     w = image.handle.write
    
##     # Circle
##     w(pr.declare("POINT_CIRCLE", pr.superellipsoid(
##         1.0, 0.50, pr.scale(RADIUS, RADIUS, RADIUS))))
##     # Square
##     w(pr.declare("POINT_SQUARE", pr.superellipsoid(
##         0.25, 0.25, pr.scale(RADIUS, RADIUS, RADIUS))))
##     # Diamond
##     w(pr.declare("POINT_DIAMOND", pr.superellipsoid(
##         1.75, 0.50,
##         pr.scale(RADIUS, RADIUS, RADIUS),
##         #pr.scale(RADIUS*1.25, RADIUS*1.25, RADIUS*1.25),
##         #2.50, 2.00, pr.scale(RADIUS*1.25, RADIUS*1.25, RADIUS*1.25),
##         )))

def _draw_background(image):
    import povray as pr
    import graphconst as gc
    
    BACKGROUND_DEPTH = 10   # extends into the Z axis.
    BACKGROUND_COLOR = 1, 1, 1

    w = image.handle.write
    w(pr.background(pr.color(*BACKGROUND_COLOR)))
    w("\n")

    # Draw a box for the background.  But make sure it doesn't appear
    # in the list of graph objects.

    SIZE = 100
    coord = -image.width*SIZE/2, -image.height*SIZE/2, 0
    extent = image.width*SIZE, image.height*SIZE, -BACKGROUND_DEPTH
    objs = image._objects[:]
    box(image, coord, extent, BACKGROUND_COLOR, finish=gc.SIMPLE)
    #box(image, coord, extent, BACKGROUND_COLOR, finish=gc.ROUGH)
    image._objects = objs

    #coord = -image.width*SIZE/2, image.height, -image.depth*SIZE/2
    #extent = image.width*SIZE, BACKGROUND_DEPTH, image.depth*SIZE
    #box(image, coord, extent, BACKGROUND_COLOR, finish=gc.SIMPLE)
