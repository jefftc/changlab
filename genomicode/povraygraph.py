"""

Functions:
povray
scatter
scatter_3d

place_ticks

"""
import os, sys

def _choose_limits(X, ticks):
    import math
    x_len = max(X) - min(X)
    x_min = min(X) - x_len*0.1
    x_max = max(X) + x_len*0.1
    if ticks:
        x_min = min(x_min, *ticks)
        x_max = max(x_max, *ticks)
    return x_min, x_max

def _choose_tick_delta(v_min, v_max, num_ticks=10):
    # Figure out the right delta to make at most num_ticks tick marks
    # between v_min and v_max.
    import math

    assert v_min < v_max
    assert num_ticks > 0 and num_ticks <= 100

    # Generate a list of the allowable DELTAs.  The resolution depends
    # on the number of tickmarks the user wants.
    DELTAS = [0.5, 0.25, 0.20, 0.1]
    if num_ticks > 10:
        DELTAS = DELTAS + [x/10 for x in DELTAS]
    DELTAS = [5.0, 2.5, 2.0, 1.0] + DELTAS
    
    # Normalize the range to 1.0 or above.
    range = v_max - v_min
    x = math.log(range, 10)
    multiplier = 10**-int(x)

    # Calculate the ideal delta, and choose the smallest DELTAS that
    # is >= the ideal one.  Assume DELTAS sorted from largest to
    # smallest.
    delta_ideal = float(range) * multiplier / num_ticks
    #print range, multiplier, delta_ideal
    x = [x for x in DELTAS if x >= delta_ideal]
    delta = x[-1]
    delta = delta / multiplier

    return delta

def test_choose_tick_delta():
    print choose_tick_delta(0, 1)
    print choose_tick_delta(0, 100)
    print choose_tick_delta(0, 21)
    print choose_tick_delta(0, 0.8, 2)
    print choose_tick_delta(0, 1, 40)
    print choose_tick_delta(-12, 15, 2)
    print choose_tick_delta(0, 99, 2)

def place_ticks(v_min, v_max, num_ticks=10, delta=None):
    import math
    
    assert v_min < v_max
    assert num_ticks > 0 and num_ticks <= 100

    if delta is None:
        delta = _choose_tick_delta(v_min, v_max, num_ticks=num_ticks)

    # Do calculation in integers to 4 decimal places.
    x = math.log(delta, 10)
    multiplier = 10**(-int(x)+4)

    delta = int(delta * multiplier)
    tick_min, tick_max = int(v_min*multiplier), int(v_max*multiplier)

    # Round tick_min down to the nearest delta and tick_max up to the
    # nearest delta.
    tick_min = tick_min - tick_min % delta
    if tick_max % delta != 0:
        tick_max = tick_max + (delta - tick_max % delta)

    ticks = [float(i)/multiplier for i in range(tick_min, tick_max+1, delta)]
    return ticks

def test_place_ticks():
    print place_ticks(0.00, 1.00, 5)
    print place_ticks(0.03, 0.94, 3)
    print place_ticks(-12, 15, 10)

def povray(
    filename, outfile=None, height=None, width=None, antialias=None,
    quality=None, povray_bin=None):
    # Return a handle to the results.
    # antialias  Which colors to anti-alias (0-3.0).
    # quality    0-9.
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

DEFAULT, CIRCLE, SQUARE, DIAMOND = range(4)
def scatter(
    X, Y, error_bar=None, color=None, xlim=None, ylim=None, pch=None,
    point_size=1.0,
    xtick=None, xtick_label=None, ytick=None, ytick_label=None, tick_size=1.0,
    xlabel=None, ylabel=None, label_size=1.0, font=None,
    plot_width=None, plot_height=None):
    # X, Y        Parallel lists that indicate the coordinates of the points.
    # error_bar   List of sizes of the error bar for each point.
    # color       List of the colors for each point (<r>,<g>,<b>) from 0-1.
    # xlim        Tuple of (minimum x, maximum x).
    # pch         DEFAULT, CIRCLE, SQUARE, or DIAMOND
    # point_size  Scales the size of each point.
    # xtick       Coordinates of the tickmarks on the x-axis (or True).
    # xtick_label Labels for the tickmarks.  Should be parallel to xtick.
    # tick_size   Scales the size of the tick mark labels.
    # xlabel      Label for the X-axes.
    # label_size  Scales the size of the labels.
    import math
    import operator
    from StringIO import StringIO
    import povray as pr

    # Constants.
    TOTAL_WIDTH = plot_width or 1024   # Size of entire drawing area.
    TOTAL_HEIGHT = plot_height or 768
    TOTAL_X_MID = TOTAL_WIDTH/2.0
    TOTAL_Y_MID = TOTAL_HEIGHT/2.0
    UNIT = TOTAL_WIDTH/100.0           # Scale to this unit.

    CAMERA_HEIGHT = UNIT*20
    LIGHT_ANGLE = 55
    LIGHT_COLOR = 1, 1, 1
    BACKGROUND_DEPTH = UNIT
    BACKGROUND_COLOR = 1, 1, 1
    AXIS_THICKNESS = UNIT*0.2          # Thickness of the axes.
    AXIS_COLOR = 0.4, 0.4, 0.4
    TICK_THICKNESS = UNIT*0.1          # Thickness of the ticks.
    TICK_COLOR = AXIS_COLOR
    TICK_WIDTH = AXIS_THICKNESS*5
    TICK_SIZE = UNIT*2*tick_size, UNIT*2*tick_size, 1
    #RADIUS = UNIT*0.75*point_size      # Radius of each point.
    RADIUS = UNIT*0.65*point_size      # Radius of each point.
    ERROR_THICKNESS = UNIT*0.1         # Error bars.
    ERROR_WIDTH = RADIUS*2*1.2
    LABEL_SIZE = UNIT*3*label_size, UNIT*3*label_size, 1
    LABEL_COLOR = AXIS_COLOR

    BORDER = 0.25*min(TOTAL_WIDTH, TOTAL_HEIGHT)
    X1_BORDER = X2_BORDER = Y1_BORDER = Y2_BORDER = BORDER/2.0
    if xlabel:
        X1_BORDER += LABEL_SIZE[0]
    if xtick:
        X1_BORDER += TICK_SIZE[0]
    if ylabel:
        Y1_BORDER += LABEL_SIZE[1]
    if ytick:
        Y1_BORDER += TICK_SIZE[1]
    
    PLOT_WIDTH = TOTAL_WIDTH-(X1_BORDER+X2_BORDER)
    PLOT_HEIGHT = TOTAL_HEIGHT-(Y1_BORDER+Y2_BORDER)
    PLOT_X_MID = X1_BORDER+PLOT_WIDTH/2.0
    PLOT_Y_MID = Y1_BORDER+PLOT_HEIGHT/2.0
    PLOT_X_MIN = X1_BORDER
    PLOT_X_MAX = PLOT_X_MIN + PLOT_WIDTH
    PLOT_Y_MIN = Y1_BORDER
    PLOT_Y_MAX = PLOT_Y_MIN + PLOT_HEIGHT

    # Check the inputs
    assert len(X) == len(Y)
    assert not error_bar or len(X) == len(error_bar)
    color = color or [(0, 0, 0)] * len(X)
    assert len(color) == len(X)
    assert len(color) == len(Y)
    pch = pch or [DEFAULT] * len(X)
    assert len(pch) == len(X)
    if xtick and not operator.isSequenceType(xtick):
        xtick = place_ticks(min(X), max(X), num_ticks=5)
    if ytick and not operator.isSequenceType(ytick):
        ytick = place_ticks(min(Y), max(Y), num_ticks=5)
    xtick = xtick or []
    ytick = ytick or []
    xtick_label = xtick_label or []
    ytick_label = ytick_label or []
    if xtick and not xtick_label:
        xtick_label = [str(x) for x in xtick]
    if ytick and not ytick_label:
        ytick_label = [str(x) for x in ytick]
    assert len(xtick) == len(xtick_label)
    assert len(ytick) == len(ytick_label)
    xlim = xlim or _choose_limits(X, xtick)   # set after xtick
    ylim = ylim or _choose_limits(Y, ytick)
    default_font = os.path.join(os.path.split(__file__)[0], "Verdana.ttf")
    font = font or default_font

    # Each error can be a single number or a tuple of (lower, upper).
    # Normalize so each of them are tuples of (lower, upper).
    if error_bar:
        for i in range(len(error_bar)):
            x = error_bar[i]
            if operator.isNumberType(x):
                err_l = err_u = x
            elif operator.isSequenceType(x):
                assert len(x) == 2
                err_l, err_u = x
            else:
                raise AssertionError, "Unknown error: %s" % repr(x)
            assert err_l >= 0 and err_u >= 0, "Errors: %s, %s" % (err_l, err_u)
            error_bar[i] = err_l, err_u

    # Rescale X, Y, and error_bar to fit the coordinates of the plot.
    (x_min, x_max), (y_min, y_max) = xlim, ylim
    x_scale = PLOT_WIDTH / (x_max - x_min)
    y_scale = PLOT_HEIGHT / (y_max - y_min)
    x_off, y_off = PLOT_X_MIN, PLOT_Y_MIN
    X_plot, Y_plot = [None]*len(X), [None]*len(Y)
    error_bar_plot = None
    if error_bar:
        error_bar_plot = [None] * len(error_bar)
    for i in range(len(X)):
        x, y = X[i], Y[i]
        X_plot[i] = (x-x_min)*x_scale + x_off
        Y_plot[i] = (y-y_min)*y_scale + y_off
        if error_bar:
            x1, x2 = error_bar[i]
            error_bar_plot[i] = x1*y_scale, x2*y_scale

    # Scale the axes to the plot's coordinate system.
    #x_axis_at, y_axis_at = 0, 0
    x_axis_at, y_axis_at = y_min, x_min
    x_axis_1 = (x_min-x_min)*x_scale+x_off, (x_axis_at-y_min)*y_scale+y_off, 0
    x_axis_2 = (x_max-x_min)*x_scale+x_off, (x_axis_at-y_min)*y_scale+y_off, 0
    y_axis_1 = (y_axis_at-x_min)*x_scale+x_off, (y_min-y_min)*y_scale+y_off, 0
    y_axis_2 = (y_axis_at-x_min)*x_scale+x_off, (y_max-y_min)*y_scale+y_off, 0

    # Extend the X axis out to the left so there's no gap at the
    # intersection.  Subtract AXIS_THICKNESS instead of
    # AXIS_THICKNESS/2.0 because Y axis is not really in the right
    # place.  We rotated it from the AXIS object around one corner of
    # the axis instead of around the middle, so it's actually too far
    # to the left by 1/2 AXIS_THICKNESS.
    x = x_axis_1[0] - AXIS_THICKNESS
    x_axis_1 = x, x_axis_1[1], x_axis_1[2]

    # Scale xtick and ytick.
    xtick_plot = [None] * len(xtick)
    for i in range(len(xtick)):
        x = xtick[i]
        xtick_plot[i] = (x-x_min)*x_scale+x_off
    ytick_plot = [None] * len(ytick)
    for i in range(len(ytick)):
        y = ytick[i]
        ytick_plot[i] = (y-y_min)*y_scale+y_off


    handle = StringIO()
    w = handle.write

    w(pr.declare("FONTFILE", '"%s"' % font)+"\n")
    w("\n")
    
    w(pr.camera(
        pr.projection("orthographic"),
        pr.location(TOTAL_X_MID, TOTAL_Y_MID, -CAMERA_HEIGHT),
        pr.right(TOTAL_WIDTH, 0, 0), pr.up(0, TOTAL_HEIGHT, 0),
        pr.look_at(TOTAL_X_MID, TOTAL_Y_MID, 0),
        ))
    x, y = PLOT_X_MAX, PLOT_Y_MAX
    hyp = math.sqrt(PLOT_WIDTH**2+PLOT_HEIGHT**2)
    z = hyp * math.tan(math.pi*LIGHT_ANGLE/180.0)
    w(pr.light_source(
        pr.vector(x, y, -z), pr.color(*LIGHT_COLOR), pr.light_type("parallel"),
        ))
    w("\n")

    w("// Draw the background plane.\n")
    w(pr.box(
        pr.vector(0, 0, 0),
        pr.vector(TOTAL_WIDTH, TOTAL_HEIGHT, BACKGROUND_DEPTH),
        pr.pigment(pr.color(*BACKGROUND_COLOR)),
        pr.finish(pr.ambient(0.75, 0.75, 0.75))))
    w("\n")
    #w(pr.background(pr.color(1, 1, 1)))
    #w("\n")

    w("// MARKER tracks the bottom right coordinates.\n")
    w(pr.declare("MARKER", pr.sphere(
        pr.vector(0, 0, 0), 10.0, pr.pigment(pr.color(0, 0, 0)))))
    

    w("// Draw axes.\n")
    w(pr.declare("AXIS", pr.cylinder(
        pr.vector(0, 0, 0), pr.vector(1, 0, 0), AXIS_THICKNESS,
        pr.pigment(pr.color(*AXIS_COLOR)),
        pr.finish(
        #pr.phong(0.8), pr.phong_size(5), 
        pr.diffuse(0.6), pr.ambient(0.6, 0.6, 0.6)),
        pr.no_shadow())))
    w(pr.object_(
        "AXIS", pr.scale(x_axis_2[0]-x_axis_1[0], 1, 1),
        pr.translate(*x_axis_1)))
    w(pr.object_(
        "AXIS", pr.rotate(0, 0, 90), pr.scale(1, y_axis_2[1]-y_axis_1[1], 1),
        pr.translate(*y_axis_1)))
    w(pr.declare("MARKER", pr.object_(
        "MARKER", pr.translate(x_axis_1[0], y_axis_1[1], 0))))
    w("\n")

    w("// Draw the tick marks.\n")
    #w(pr.object_("MARKER"))
    w(pr.declare("TICK", pr.cylinder(
        pr.vector(0, 0, 0), pr.vector(1, 0, 0), TICK_THICKNESS,
        pr.pigment(pr.color(*TICK_COLOR)),
        pr.finish(
        pr.phong(0.8), pr.phong_size(5), pr.diffuse(0.6),
        pr.ambient(0.6, 0.6, 0.6)),
        pr.no_shadow())))

    # Bug: Places marker based on number of characters in the name.
    # Does not take the rendered size into account.
    longest_i = None
    for i in range(len(ytick_label)):
        if longest_i is None or (
            len(ytick_label[i]) > len(ytick_label[longest_i])):
            longest_i = i
    for i in range(len(ytick_plot)):
        # Draw the tickmark.
        w(pr.object_(
            "TICK", pr.scale(TICK_WIDTH, 1, 1),
            pr.translate(y_axis_1[0]-TICK_WIDTH/2.0, ytick_plot[i], 0)))
        # Label it.
        name = "YTICK_%d" % i
        w(pr.declare(name, pr.text(
            "FONTFILE", ytick_label[i], RADIUS, 0,
            pr.pigment(pr.color(*LABEL_COLOR)),
            pr.scale(*TICK_SIZE),
            pr.translate(y_axis_1[0]-TICK_WIDTH/2.0, ytick_plot[i], -RADIUS),
            pr.no_shadow())))
        w(pr.object_(name, pr.translate_by_extent(name, -1.2, -0.5, 0)))
        if i == longest_i:
            # Shift for the tick mark.
            x = "-min_extent(MARKER).x+%g" % (y_axis_1[0]-TICK_WIDTH/2.0)
            w(pr.declare("MARKER", pr.object_("MARKER", pr.translate(x,0,0))))
            # Shift for the tick label.
            w(pr.declare("MARKER", pr.object_(
                "MARKER", pr.translate_by_extent(name, -1.2, 0, 0))))
    #w(pr.object_("MARKER"))
            
    for i in range(len(xtick_plot)):
        w(pr.object_(
            "TICK", pr.scale(TICK_WIDTH, 1, 1),
            pr.rotate(0, 0, 90),
            pr.translate(xtick_plot[i], x_axis_1[1]-TICK_WIDTH/2.0, 0)))
        name = "XTICK_%d" % i
        w(pr.declare(name, pr.text(
            "FONTFILE", xtick_label[i], RADIUS, 0,
            pr.pigment(pr.color(*LABEL_COLOR)),
            pr.scale(*TICK_SIZE),
            pr.translate(xtick_plot[i], x_axis_1[1]-TICK_WIDTH/2.0, -RADIUS),
            pr.no_shadow())))
        w(pr.object_(name, pr.translate_by_extent(name, -0.5, -1.2, 0)))
        if i == 0:
            # Shift for the tick mark.
            x = "-min_extent(MARKER).y+%g" % (x_axis_1[1]-TICK_WIDTH/2.0)
            w(pr.declare("MARKER", pr.object_("MARKER", pr.translate(0,x,0))))
            # Shift for the tick label.
            w(pr.declare("MARKER", pr.object_(
                "MARKER", pr.translate_by_extent(name, 0, -1.2, 0))))
    w("\n")
    #w(pr.object_("MARKER"))

    w("// Plot each point.\n")
    finish = pr.finish(
        #pr.phong(0.8), pr.phong_size(10), pr.diffuse(0.4),
        pr.phong(0.4), pr.phong_size(5), pr.diffuse(0.4),
        pr.ambient(0.3, 0.3, 0.3))
        #pr.ambient(0.4, 0.4, 0.4))
    #w(pr.declare("POINT", pr.sphere(pr.vector(0, 0, 0), RADIUS, finish)))
    # Circle
    w(pr.declare("POINT_CIRCLE", pr.superellipsoid(
        1.0, 0.50, pr.scale(RADIUS, RADIUS, RADIUS), finish)))
    # Square
    w(pr.declare("POINT_SQUARE", pr.superellipsoid(
        0.25, 0.25, pr.scale(RADIUS, RADIUS, RADIUS), finish)))
    # Diamond
    w(pr.declare("POINT_DIAMOND", pr.superellipsoid(
        2.50, 2.50, pr.scale(RADIUS*1.25, RADIUS*1.25, RADIUS*1.25), finish)))
    
    finish = pr.finish(pr.diffuse(0.1), pr.ambient(0.5, 0.5, 0.5))
    w(pr.declare("ERRBAR", pr.cylinder(
        pr.vector(0, 0, 0), pr.vector(1, 0, 0), ERROR_THICKNESS, finish)))
        
    for i in range(len(X)):
        x, y, col = X_plot[i], Y_plot[i], color[i]
        pigment = pr.pigment(pr.color(*col))
        w("// Point %d of %d.\n" % (i+1, len(X)))
        # DEFAULT, CIRCLE, SQUARE, DIAMOND
        if pch[i] in [DEFAULT, CIRCLE]:
            point = "POINT_CIRCLE"
        elif pch[i] == SQUARE:
            point = "POINT_SQUARE"
        elif pch[i] == DIAMOND:
            point = "POINT_DIAMOND"
        else:
            raise AssertionError, "Unknown point [%d]: %s" % (i, pch[i])
        w(pr.object_(point, pr.translate(x, y, 0), pigment))

        if not error_bar_plot:
            continue
        err_l, err_u = error_bar_plot[i]
        err_total = err_l + err_u
        if err_total:
            w(pr.object_(
                "ERRBAR", pr.rotate(0, 0, 90), pr.scale(1, err_total, 1),
                pr.translate(x, y-err_l, 0), pigment))
            w(pr.object_(
                "ERRBAR", pr.scale(ERROR_WIDTH, 1, 1),
                pr.translate(x-ERROR_WIDTH/2, y+err_u, 0), pigment))
            
        w(pr.object_(
            "ERRBAR", pr.scale(ERROR_WIDTH, 1, 1),
            pr.translate(x-ERROR_WIDTH/2, y-err_l, 0), pigment))
    w("\n")

    if xlabel:
        # Nudges label down to add spacing.
        #XLABEL_NUDGE = 0.8
        XLABEL_NUDGE = 1.2
        w(pr.declare("XLABEL", pr.text(
            "FONTFILE", xlabel, RADIUS, 0,
            pr.pigment(pr.color(*LABEL_COLOR)),
            pr.scale(*LABEL_SIZE),
            pr.translate(PLOT_X_MID, "min_extent(MARKER).y", -RADIUS),
            pr.no_shadow())))
        w(pr.object_(
            "XLABEL", pr.translate_by_extent(
            "XLABEL", -0.5, -XLABEL_NUDGE, 0)))
        w(pr.declare("MARKER", pr.object_(
            "MARKER", pr.translate_by_extent("XLABEL", 0, -XLABEL_NUDGE, 0))))
        w("\n")
    if ylabel:
        # Nudges label left to add spacing.
        YLABEL_NUDGE = 1.4
        w(pr.declare("YLABEL", pr.text(
            "FONTFILE", ylabel, RADIUS, 0,
            pr.pigment(pr.color(*LABEL_COLOR)),
            pr.rotate(0, 0, 90),
            pr.scale(*LABEL_SIZE),
            pr.translate("min_extent(MARKER).x", PLOT_Y_MID, -RADIUS),
            pr.no_shadow()
            )))
        w(pr.object_(
            "YLABEL", pr.translate_by_extent(
            "YLABEL", -0.5, -(YLABEL_NUDGE-1), 0)))
        w(pr.declare("MARKER", pr.object_(
            "MARKER", pr.translate_by_extent("YLABEL", -YLABEL_NUDGE, 0, 0))))
        w("\n")

    #w(pr.object_("MARKER"))
    
    handle.seek(0)
    return handle.read()

def scatter_3d(X, Y, Z):
    # X, Y, Z     Parallel lists that indicate the coordinates of the points.
    raise NotImplementedError

    point_size = 1.0
    
    import math
    import operator
    from StringIO import StringIO
    import povray as pr

    # Constants.
    TOTAL_WIDTH = 1024   # Size of entire drawing area.
    TOTAL_HEIGHT = 768
    UNIT = TOTAL_WIDTH/100.0           # Scale to this unit.

    X_MID = (min(X)+max(X))/2.0
    Y_MID = (min(Y)+max(Y))/2.0
    Z_MID = (min(Z)+max(Z))/2.0
    RADIUS = UNIT*0.65*point_size      # Radius of each point.
    
    # Check the inputs
    assert len(X) == len(Y)
    assert len(X) == len(Z)

    handle = StringIO()
    w = handle.write

    w(pr.background(pr.color(1, 1, 1)))
    w("\n")

    CAMERA_X = min(X)
    CAMERA_Y = min(Y)
    CAMERA_Z = min(Z)

    w(pr.camera(
        pr.projection("orthographic"),
        pr.location(CAMERA_X, CAMERA_Y, CAMERA_Z),
        pr.look_at(X_MID, Y_MID, Z_MID),
        ))

    w("// Plot each point.\n")
    finish = pr.finish(
        pr.phong(0.4), pr.phong_size(5), pr.diffuse(0.4),
        pr.ambient(0.3, 0.3, 0.3))
    # Circle
    w(pr.declare("POINT_CIRCLE", pr.superellipsoid(
        1.0, 0.50, pr.scale(RADIUS, RADIUS, RADIUS), finish)))
    
    for i in range(len(X)):
        x, y, z = X[i], Y[i], Z[i]
        w("// Point %d of %d.\n" % (i+1, len(X)))
        point = "POINT_CIRCLE"
        w(pr.object_(point, pr.translate(x, y, z)))
    w("\n")
    
    handle.seek(0)
    return handle.read()
