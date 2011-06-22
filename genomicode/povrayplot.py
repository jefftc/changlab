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

write

"""
# _make_finish
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
        self._objects = []  # list of object names
    def _make_object_name(self):
        name = "OBJECT%05d" % len(self._objects)
        self._objects.append(name)
        return name
    def format(self, outhandle):
        outhandle.write(self.handle.getvalue())

def image(width, height, depth):
    im = Image(width, height, depth)
    _declare_fontfile(im)
    _position_camera(im, width, height, depth)
    #_declare_points(im)
    _draw_background(im)
    # Clear the canvas for the client to draw on.
    del im._objects[:]
    return im

def box(image, coord, extent, color, shadow=False, finish=None):
    # coord is (x, y, z).  extent is (width, height, depth).
    import povray as pr
    finish = _make_finish(finish)
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
        pr.pigment(pr.color(*color)),
        ns,
        finish,
        )))
    w(pr.object_(name))

    #extent = coord2[0]-coord1[0], coord2[1]-coord1[1], coord2[2]-coord1[2]
    #image._add_object(name, coord1, extent)

        
def sphere(image, coord, radius, color, shadow=False, shape=None,
           finish=None):
    # coord is the center of the sphere, with radius.
    import povray as pr
    import graphconst as gc

    finish = _make_finish(finish)
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
            pr.translate(*coord), pr.pigment(pr.color(*color)), ns, finish)
    # Square
    elif shape == gc.SQUARE:
        point = pr.superellipsoid(
            0.25, 0.25, pr.scale(radius, radius, radius),
            pr.translate(*coord), pr.pigment(pr.color(*color)), ns, finish)
    # Diamond
    elif shape == gc.DIAMOND:
        point = pr.superellipsoid(
            1.75, 0.50, pr.scale(radius, radius, radius),
            pr.translate(*coord), pr.pigment(pr.color(*color)), ns, finish)
    else:
        raise AssertionError, "Unknown point: %s" % (shape)
    
    w = image.handle.write
    name = image._make_object_name()
    w(pr.declare(name, point))
    w(pr.object_(name))

    #x, y, z = coord
    #coord = x-radius, y-radius, z-radius
    #extent = radius*2, radius*2, radius*2
    #image._add_object(name, coord, extent)

def cylinder(
    image, coord, extent, radius, color, shadow=False, finish=None):
    # coord is the center of one end, with the center of the other end
    # in the direction determined by extent.
    import povray as pr
    finish = _make_finish(finish)
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
        pr.pigment(pr.color(*color)),
        ns, 
        finish,
        )))
    w(pr.object_(name))

def text(
    image, text, coord, depth, fontsize, color, shadow=False,
    center_x=False, center_y=False, center_z=False,
    wrong_x=False, wrong_y=False, wrong_z=False,
    min_x=False, max_x=False, min_y=False, max_y=False,
    vertical=False, finish=None):
    # The coordinate should be at the upper left point of the text.
    # depth grows in the +z direction (toward the user).
    # center_x means the x-coord is at the middle, and wrong_x means
    # the x-coord is at the right.  They should not both be True.
    import povray as pr
    
    assert not (center_x and wrong_x)
    assert not (center_y and wrong_y)
    assert not (center_z and wrong_z)
    assert not (min_x and max_x)
    assert not (min_y and max_y)

    finish = _make_finish(finish)

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
        pr.pigment(pr.color(*color)),
        finish,
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
            n = names[:-max_y]
        x = ["min_extent(%s).y" % x for x in n]
        x = "min(\n  %s)" % ",\n  ".join(x)
        y_trans = x
    if x_trans or y_trans:
        w(pr.declare(name, pr.object_(
            name, pr.translate(x_trans, y_trans, 0))))
        
    # Translate by extent.
    x_extent = 0
    y_extent = -1   # povray coordinates are from the lower left
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
    
    
def text90(*args, **keywds):
    keywds["vertical"] = True
    return text(*args, **keywds)

## def line(image, x, y, width, height, color):
##     raise NotImplementedError

def write(image, handle):
    # Format the image with POV-RAY and write the results to handle.
    import os
    import tempfile
    import povray

    pov_file = out_file = None
    try:
        x, pov_file = tempfile.mkstemp(dir=".", suffix=".pov"); os.close(x)
        x, out_file = tempfile.mkstemp(dir=".", suffix=".png"); os.close(x)
        if os.path.exists(out_file):
            os.unlink(out_file)
            
        image.format(open(pov_file, 'w'))
        r = povray.povray(
            pov_file, outfile=out_file, 
            height=image.height, width=image.width, antialias=0.5, quality=9)
        output = r.read()
        assert os.path.exists(out_file), "POV-RAY failed.\n%s" % output
        handle.write(open(out_file).read())
    finally:
        if pov_file and os.path.exists(pov_file):
            os.unlink(pov_file)
        if out_file and os.path.exists(out_file):
            os.unlink(out_file)
    return output

def _make_finish(finish):
    import povray as pr
    import graphconst as gc

    if finish is None:
        finish = gc.SIMPLE
        
    if finish == gc.SIMPLE:
        x = pr.finish(
            pr.ambient(0.8, 0.8, 0.8),
            pr.diffuse(0.3),
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
        x = pr.finish(
            pr.ambient(0.45),
            pr.brilliance(2),
            pr.diffuse(0.3),
            pr.metallic(),
            pr.specular(0.80),
            pr.roughness(1./20),
            pr.reflection(0.1),
            )
    elif finish == gc.ROUGH:
        x = pr.finish(
            pr.ambient(0.8, 0.8, 0.8),
            pr.diffuse(0.4),
            pr.crand(0.2),
            )
    elif finish == gc.SHINY:
        # Doesn't work very well for graphs.
        x = pr.finish(
            pr.ambient(0.7),
            pr.diffuse(0.1),
            pr.reflection(0.3),
            pr.ior(1.5),
            )
    else:
        raise AssertionError, "Unknown finish: %s" % finish
    return x

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

def _position_camera(image, width, height, depth):
    import math
    import povray as pr
    
    CAMERA_HEIGHT = depth * 1.1
    LIGHT_ANGLE = 78        # lower means longer shadows, darker colors
    LIGHT_COLOR = 1, 1, 1

    # Set the camera to the middle of the plot, looking down.
    # Have to look at the middle of the entire plot, or else the
    # borders will be off.
    x_mid, y_mid, z_mid = width*0.5, height*0.5, depth*0.5
    if 1:
        # Default camera location.
        camera = (x_mid, y_mid, -CAMERA_HEIGHT)
        look = (x_mid, y_mid, z_mid)
        dist_scale = 1.0
    elif 0:
        camera = (width*0.60, height*0.50, -depth*1.2)
        look = (x_mid*0.25, y_mid*0.35, z_mid)
        dist_scale = 1.50
    elif 0:
        camera = (width*0.70, height*0.50, -depth*1.2)
        look = (x_mid*0.25, y_mid*0.35, z_mid)
        dist_scale = 1.35
    
    w = image.handle.write
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

    SIZE = 100
    coord = -image.width*SIZE/2, -image.height*SIZE/2, 0
    extent = image.width*SIZE, image.height*SIZE, -BACKGROUND_DEPTH
    box(image, coord, extent, BACKGROUND_COLOR, finish=gc.SIMPLE)

    #coord = -image.width*SIZE/2, image.height, -image.depth*SIZE/2
    #extent = image.width*SIZE, BACKGROUND_DEPTH, image.depth*SIZE
    #box(image, coord, extent, BACKGROUND_COLOR, finish=gc.SIMPLE)
