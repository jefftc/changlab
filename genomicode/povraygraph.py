"""

Functions:
povray

scatter
line
bar

"""
# _make_graph
# _set_default_color
# _set_default_shape
# _set_default_label
# _set_default_error_bar
# _set_default_axes    lim, tick, tick_label
# _set_default_lim
# _set_default_tick
# _set_default_tick_label
# _coord2pixel


import os, sys

DEFAULT, CIRCLE, SQUARE, DIAMOND = range(4)


class Graph:
    # TOTAL_WIDTH   Number of pixels for the whole canvas.
    # TOTAL_HEIGHT
    # PLOT_X        The pixel coordinate where the drawable area starts.
    # PLOT_Y
    # PLOT_WIDTH    Number of pixels for the drawable area.
    # PLOT_HEIGHT
    # 
    # UNIT          A standard unit for setting relative sizes of things.
    #
    # _xlim
    # _ylim
    #
    # Methods:
    # declare_fontfile
    # position_camera
    # draw_background
    # draw_axes
    # draw_tick_marks
    # draw_tick_labels
    # draw_labels
    # 
    # draw_points
    # draw_error_bars
    # draw_onpoint_labels
    # draw_overpoint_labels
    # draw_line
    # draw_bar
    # 
    # draw
    # 
    # _get_axis_coords
    # _z_height
    def __init__(self, handle, xlim, ylim, width=None, height=None):
        self.handle = handle

        # Constants.
        self.TOTAL_WIDTH = width or 1024
        self.TOTAL_HEIGHT = height or 768
        self.UNIT = self.TOTAL_WIDTH/100.0
        self.POINT_RADIUS = 0.6 * self.UNIT  # Radius of each point.

        # Calculate the borders around the plot.
        BORDER = 0.25*min(self.TOTAL_WIDTH, self.TOTAL_HEIGHT)
        X1_BORDER = BORDER * 0.65   # Make a little larger for the labels.
        X2_BORDER = BORDER - X1_BORDER
        Y1_BORDER = BORDER * 0.6
        Y2_BORDER = BORDER - Y1_BORDER

        self.PLOT_WIDTH = self.TOTAL_WIDTH-(X1_BORDER+X2_BORDER)
        self.PLOT_HEIGHT = self.TOTAL_HEIGHT-(Y1_BORDER+Y2_BORDER)
        self.PLOT_X = X1_BORDER
        self.PLOT_Y = Y1_BORDER
        self.DRAWN = {}
        
        self._xlim = xlim
        self._ylim = ylim

    def declare_fontfile(self, fontfile=None):
        import povray as pr

        if not fontfile:
            x = os.path.join(os.path.split(__file__)[0], "Verdana Bold.ttf")
            fontfile = os.path.realpath(x)
        fontfile = fontfile or default_font
        assert os.path.exists(fontfile)

        w = self.handle.write
        w(pr.declare("FONTFILE", '"%s"' % fontfile)+"\n")
        w("\n")
        self.DRAWN["fontfile"] = 1

    def position_camera(self):
        import math
        import povray as pr

        CAMERA_HEIGHT = 10*self.UNIT
        LIGHT_ANGLE = 70   # lower means longer shadows, darker colors
        LIGHT_COLOR = 1, 1, 1

        w = self.handle.write
        # Set the camera to the middle of the plot, looking down.
        x_mid = self.TOTAL_WIDTH/2.0
        y_mid = self.TOTAL_HEIGHT/2.0
        w(pr.camera(
            pr.projection("orthographic"),
            pr.location(x_mid, y_mid, -CAMERA_HEIGHT),
            pr.right(self.TOTAL_WIDTH, 0, 0),
            pr.up(0, self.TOTAL_HEIGHT, 0),
            pr.look_at(x_mid, y_mid, 0),
            ))
        
        # The light source is at the upper right of the plot, shining
        # toward the lower left.  Calculate the height such that the
        # light is shining at a specific angle relative to the plane
        # of the plot.  Lower LIGHT_ANGLE means light source is closer
        # to the plane (longer shadows).
        x = self.PLOT_X + self.PLOT_WIDTH
        y = self.PLOT_Y + self.PLOT_HEIGHT
        side = math.sqrt(self.PLOT_WIDTH**2+self.PLOT_HEIGHT**2)
        z = side * math.tan(math.pi*LIGHT_ANGLE/180.0)
        w(pr.light_source(
            pr.vector(x, y, -z),
            pr.color(*LIGHT_COLOR),
            pr.light_type("parallel"),
            ))
        w("\n")
        self.DRAWN["camera"] = 1

    def draw_background(self):
        import povray as pr
        
        BACKGROUND_DEPTH = self.UNIT
        BACKGROUND_COLOR = 1, 1, 1

        w = self.handle.write
        w("// Draw the background plane.\n")
        w(pr.box(
            pr.vector(0, 0, 0),
            pr.vector(self.TOTAL_WIDTH, self.TOTAL_HEIGHT, BACKGROUND_DEPTH),
            pr.pigment(pr.color(*BACKGROUND_COLOR)),
            pr.finish(pr.ambient(0.75, 0.75, 0.75))))
        w("\n")
        
        w("// MARKER tracks the bottom left coordinates.\n")
        w(pr.declare("MARKER", pr.sphere(
            pr.vector(2.5, 2.5, 2.5), 5.0,
            pr.finish(pr.ambient(0.6, 0.6, 0.6)))))
        w("\n")
        self.DRAWN["background"] = 1

    def draw_axes(self):
        import povray as pr

        AXIS_THICKNESS = 0.15*self.UNIT     # Thickness of the axes.
        AXIS_COLOR = 0.6, 0.6, 0.6

        x = self._get_axis_coords()
        (x_axis_1, x_axis_2), (y_axis_1, y_axis_2) = x

        # Extend the X axis out to the left so there's no gap at the
        # intersection.  Subtract AXIS_THICKNESS instead of
        # AXIS_THICKNESS/2.0 because Y axis is not really in the right
        # place.  We rotated it from the AXIS object around one corner of
        # the axis instead of around the middle, so it's actually too far
        # to the left by 1/2 AXIS_THICKNESS.
        x = x_axis_1[0] - AXIS_THICKNESS
        x_axis_1 = x, x_axis_1[1], x_axis_1[2]
        
        w = self.handle.write
        w("// Draw axes.\n")
        w(pr.declare("AXIS", pr.cylinder(
            pr.vector(0, 0, 0), pr.vector(1, 0, 0), AXIS_THICKNESS,
            pr.pigment(pr.color(*AXIS_COLOR)),
            pr.finish(
                pr.diffuse(0.6),
                pr.ambient(0.6, 0.6, 0.6),
                pr.phong(0.2),
                pr.phong_size(5),
                #pr.phong(0.8), pr.phong_size(5),   # too shiny
                ),
            pr.no_shadow(),
            )))
        # Give the axis a rounded cap so it looks smoother.
        w(pr.declare("AXIS_CAP", pr.sphere(
            pr.vector(0, 0, 0), AXIS_THICKNESS,
            pr.pigment(pr.color(*AXIS_COLOR)),
            pr.finish(
                pr.diffuse(0.6),
                pr.ambient(0.6, 0.6, 0.6),
                pr.phong(0.3),
                pr.phong_size(10),
                #pr.phong(0.8), pr.phong_size(5),   # too shiny
                ),
            pr.no_shadow(),
            )))
        # X-axis.
        w(pr.object_(
            "AXIS",
            pr.scale(x_axis_2[0]-x_axis_1[0], 1, 1),
            pr.translate(*x_axis_1)))
        w(pr.object_(
            "AXIS_CAP",
            pr.translate(*x_axis_2)))
        # Y-axis.
        w(pr.object_(
            "AXIS",
            pr.rotate(0, 0, 90),
            pr.scale(1, y_axis_2[1]-y_axis_1[1], 1),
            pr.translate(*y_axis_1)))
        w(pr.object_(
            "AXIS_CAP",
            pr.translate(*y_axis_2)))

        w(pr.declare("MARKER", pr.object_(
            "MARKER", pr.translate(x_axis_1[0], y_axis_1[1], 0))))
        #w(pr.object_("MARKER"))
        w("\n")
        self.DRAWN["axes"] = 1

    def _get_axis_coords(self):
        # Return (x_axis_1, x_axis_2), (y_axis_1, y_axis_2).
        x_min = self.PLOT_X
        x_max = self.PLOT_X + self.PLOT_WIDTH
        y_min = self.PLOT_Y
        y_max = self.PLOT_Y + self.PLOT_HEIGHT
        x_axis_at = self.PLOT_Y
        y_axis_at = self.PLOT_X
        x_axis_1 = (x_min, x_axis_at, 0)
        x_axis_2 = (x_max, x_axis_at, 0)
        y_axis_1 = (y_axis_at, y_min, 0)
        y_axis_2 = (y_axis_at, y_max, 0)
        return (x_axis_1, x_axis_2), (y_axis_1, y_axis_2)

    def draw_tick_marks(self, xtick, ytick, tick_size=1.0):
        import povray as pr

        TICK_THICKNESS = 0.1*self.UNIT*tick_size # Thickness of the tick marks
        TICK_COLOR = 0.6, 0.6, 0.6
        TICK_WIDTH = 1.2 * self.UNIT

        x = self._get_axis_coords()
        (x_axis_1, x_axis_2), (y_axis_1, y_axis_2) = x

        # Scale xtick and ytick.
        xlim, ylim = self._xlim, self._ylim
        x_min, x_len = xlim[0], xlim[1]-xlim[0]
        y_min, y_len = ylim[0], ylim[1]-ylim[0]
        xtick_pixel = [
            _coord2pixel(x, x_min, x_len, self.PLOT_X, self.PLOT_WIDTH)
            for x in xtick]
        ytick_pixel = [
            _coord2pixel(y, y_min, y_len, self.PLOT_Y, self.PLOT_HEIGHT)
            for y in ytick]

        w = self.handle.write
        w("// Draw the tick marks.\n")
        w(pr.declare("TICK", pr.cylinder(
            pr.vector(0, 0, 0), pr.vector(1, 0, 0), TICK_THICKNESS,
            pr.pigment(pr.color(*TICK_COLOR)),
            pr.finish(
                #pr.phong(0.8), pr.phong_size(5),
                pr.diffuse(0.6),
                pr.ambient(0.6, 0.6, 0.6),
                pr.phong(0.3),
                pr.phong_size(10),
                ),
            pr.no_shadow())))

        for i in range(len(xtick_pixel)):
            name = "XTICK_%d" % i
            w(pr.declare(name, pr.object_(
                "TICK", pr.scale(TICK_WIDTH, 1, 1),
                pr.rotate(0, 0, 90),
                pr.translate(xtick_pixel[i], x_axis_1[1]-TICK_WIDTH/2.0, 0))))
            w(pr.object_(name))
        for i in range(len(ytick_pixel)):
            name = "YTICK_%d" % i
            w(pr.declare(name, pr.object_(
                "TICK", pr.scale(TICK_WIDTH, 1, 1),
                pr.translate(y_axis_1[0]-TICK_WIDTH/2.0, ytick_pixel[i], 0))))
            w(pr.object_(name))

        # Shift the marker for the ticks.
        #w(pr.object_("MARKER", pr.pigment(pr.color(1.0, 0, 0)),))
        x_trans = "-min_extent(MARKER).x+%g" % (y_axis_1[0]-TICK_WIDTH/2.0)
        y_trans = "-min_extent(MARKER).y+%g" % (x_axis_1[1]-TICK_WIDTH/2.0)
        w(pr.declare(
            "MARKER", pr.object_("MARKER", pr.translate(x_trans, y_trans, 0))))
        #w(pr.object_("MARKER", pr.pigment(pr.color(0, 1.0, 0)),))
        w("\n")
        self.DRAWN["tick_marks"] = 1

    def draw_tick_labels(self, xtick_label, ytick_label, label_size=1.0):
        import povray as pr

        LABEL_SIZE = 1.6*self.UNIT*label_size, 1.6*self.UNIT*label_size, 1
        LABEL_COLOR = 0.6, 0.6, 0.6
        THICKNESS = 0.5*self.UNIT
        SPACE = 0.5*self.UNIT

        if not xtick_label and not ytick_label:
            return
        assert "tick_marks" in self.DRAWN, "Cannot draw labels with tick marks"

        x = self._get_axis_coords()
        (x_axis_1, x_axis_2), (y_axis_1, y_axis_2) = x

        w = self.handle.write
        w("// Draw the tick labels.\n")

        for i in range(len(xtick_label)):
            name = "XTLABEL_%d" % i
            w(pr.declare("X", pr.text(
                "FONTFILE", xtick_label[i], THICKNESS, 0,
                pr.pigment(pr.color(*LABEL_COLOR)),
                pr.no_shadow(),
                pr.scale(*LABEL_SIZE),
                pr.translate(
                    "min_extent(XTICK_%d).x" % i,
                    "min_extent(XTICK_%d).y" % i,
                    -THICKNESS),
                pr.translate(0, -SPACE, 0),
                )))
            w(pr.declare("X", pr.object_(
                "X", pr.translate_by_extent("XTICK_%d"%i, 0.5, 0, 0))))
            w(pr.declare(name, pr.object_(
                "X", pr.translate_by_extent("X", -0.5, -1, 0))))
            w(pr.object_(name))

        for i in range(len(ytick_label)):
            name = "YTLABEL_%d" % i
            w(pr.declare("X", pr.text(
                "FONTFILE", ytick_label[i], THICKNESS, 0,
                pr.pigment(pr.color(*LABEL_COLOR)),
                pr.no_shadow(),
                pr.scale(*LABEL_SIZE),
                pr.translate(
                    "min_extent(YTICK_%d).x" % i,
                    "min_extent(YTICK_%d).y" % i,
                    -THICKNESS),
                pr.translate(-SPACE, 0, 0),
                )))
            w(pr.declare("X", pr.object_(
                "X", pr.translate_by_extent("YTICK_%d"%i, 0, 0.5, 0))))
            w(pr.declare(name, pr.object_(
                "X", pr.translate_by_extent("X", -1, -0.5, 0))))
            w(pr.object_(name))

        # DEBUG:
        #w(pr.declare("Min", "min_extent(YTICK_2);"))
        #w(pr.declare("Max", "max_extent(YTICK_2);"))
        #w(pr.box("Min", "Max", pr.pigment(pr.color(1.0, 0.0, 1.0))))

        #w(pr.object_("MARKER", pr.pigment(pr.color(1.0, 0.0, 0)),))
        x_trans = y_trans = 0
        if xtick_label:
            x = ["min_extent(XTLABEL_%d).y" % i
                 for i in range(len(xtick_label))]
            x = "min(\n  %s)" % ",\n  ".join(x)
            y_trans = "-min_extent(MARKER).y + %s" % x
        if ytick_label:
            x = ["min_extent(YTLABEL_%d).x" % i
                 for i in range(len(ytick_label))]
            x = "min(\n  %s)" % ",\n  ".join(x)
            x_trans = "-min_extent(MARKER).x + %s" % x
        w(pr.declare(
            "MARKER", pr.object_("MARKER", pr.translate(x_trans, y_trans, 0))))
        #w(pr.object_("MARKER", pr.pigment(pr.color(1.0, 1.0, 0)),))
        
        w("\n")
        self.DRAWN["tick_labels"] = 1

    def draw_labels(self, xlabel, ylabel, label_size=1.0):
        import povray as pr
        
        THICKNESS = 0.5*self.UNIT
        LABEL_SIZE = (
            2.0*self.UNIT*label_size, 2.0*self.UNIT*label_size, 1)
        # Color (0.6, 0.6, 0.6) is too light for pybinreg predictions.
        LABEL_COLOR = 0.2, 0.2, 0.2
        SPACE = 1.5*self.UNIT
        
        PLOT_X_MID = self.PLOT_X+self.PLOT_WIDTH/2.0
        PLOT_Y_MID = self.PLOT_Y+self.PLOT_HEIGHT/2.0
        
        w = self.handle.write

        if xlabel:
            w(pr.declare("XLABEL", pr.text(
                "FONTFILE", xlabel, 1.0, 0,
                pr.pigment(pr.color(*LABEL_COLOR)),
                pr.no_shadow(),
                pr.scale(*LABEL_SIZE),
                pr.translate(PLOT_X_MID, "min_extent(MARKER).y", -THICKNESS),
                pr.translate(0, -SPACE, 0),
                )))
            w(pr.declare("XLABEL", pr.object_(
                "XLABEL", pr.translate_by_extent("XLABEL", -0.5, -1, 0))))
            w(pr.object_("XLABEL"))

        if ylabel:
            w(pr.declare("YLABEL", pr.text(
                "FONTFILE", ylabel, 1.0, 0,
                pr.pigment(pr.color(*LABEL_COLOR)),
                pr.no_shadow(),
                pr.scale(*LABEL_SIZE),
                pr.rotate(0, 0, 90),
                pr.translate("min_extent(MARKER).x", PLOT_Y_MID, -THICKNESS),
                pr.translate(-SPACE, 0, 0),
                )))
            w(pr.declare("YLABEL", pr.object_(
                "YLABEL", pr.translate_by_extent("YLABEL", 0, -0.5, 0))))
            w(pr.object_("YLABEL"))

        x_trans = y_trans = 0
        if xlabel:
            y_trans = "-min_extent(MARKER).y+min_extent(XLABEL).y"
        if ylabel:
            x_trans = "-min_extent(MARKER).x+min_extent(YLABEL).x"
        w(pr.declare("MARKER", pr.object_(
            "MARKER", pr.translate(x_trans, y_trans, 0))))
        #w(pr.object_("MARKER", pr.pigment(pr.color(0, 1.0, 0)),))
        w("\n")

    def draw_points(
        self, points, color, shape, point_size=1.0, height=0):
        import povray as pr
        
        RADIUS = self.POINT_RADIUS*point_size      # Radius of each point.
        HEIGHT = self._z_height(height)

        xlim, ylim = self._xlim, self._ylim
        x_min, x_len = xlim[0], xlim[1]-xlim[0]
        y_min, y_len = ylim[0], ylim[1]-ylim[0]

        w = self.handle.write
        w("// Plot each point.\n")
        finish = pr.finish(
            #pr.phong(0.8), pr.phong_size(10), pr.diffuse(0.4),
            pr.phong(0.4), pr.phong_size(5), pr.diffuse(0.4),
            pr.ambient(0.3, 0.3, 0.3))
            #pr.ambient(0.4, 0.4, 0.4))
        # Circle
        w(pr.declare("POINT_CIRCLE", pr.superellipsoid(
            1.0, 0.50, pr.scale(RADIUS, RADIUS, RADIUS), finish)))
        # Square
        w(pr.declare("POINT_SQUARE", pr.superellipsoid(
            0.25, 0.25, pr.scale(RADIUS, RADIUS, RADIUS), finish)))
        # Diamond
        w(pr.declare("POINT_DIAMOND", pr.superellipsoid(
            1.75, 0.50,
            pr.scale(RADIUS, RADIUS, RADIUS),
            #pr.scale(RADIUS*1.25, RADIUS*1.25, RADIUS*1.25),
            #2.50, 2.00, pr.scale(RADIUS*1.25, RADIUS*1.25, RADIUS*1.25),
            finish)))
        # WEIRD SHAPE
        #w(pr.declare("POINT_OTHER", pr.superellipsoid(
        #    1.00, 2.00, pr.scale(RADIUS, RADIUS, RADIUS),
        #    finish)))
        for i in range(len(points)):
            x_coord, y_coord = points[i]

            x_pixel = _coord2pixel(
                x_coord, x_min, x_len, self.PLOT_X, self.PLOT_WIDTH)
            y_pixel = _coord2pixel(
                y_coord, y_min, y_len, self.PLOT_Y, self.PLOT_HEIGHT)

            # Don't draw points that are off the drawing area.
            if x_pixel<self.PLOT_X or x_pixel >= self.PLOT_X+self.PLOT_WIDTH:
                continue
            if y_pixel<self.PLOT_Y or y_pixel >= self.PLOT_Y+self.PLOT_HEIGHT:
                continue
            
            pigment = pr.pigment(pr.color(*color[i]))
            w("// Point %d of %d.\n" % (i+1, len(points)))
            # DEFAULT, CIRCLE, SQUARE, DIAMOND
            if shape[i] in [DEFAULT, CIRCLE]:
                point = "POINT_CIRCLE"
            elif shape[i] == SQUARE:
                point = "POINT_SQUARE"
            elif shape[i] == DIAMOND:
                point = "POINT_DIAMOND"
            #elif shape[i] == OTHER:
            #    point = "POINT_OTHER"
            else:
                raise AssertionError, "Unknown point [%d]: %s" % (i, shape[i])
            z_pixel = -HEIGHT
            w(pr.object_(
                point, pr.translate(x_pixel, y_pixel, z_pixel), pigment))
                
        w("\n")

    def draw_error_bars(self, points, error_bar, color, height=0):
        import povray as pr
        
        ERROR_THICKNESS = 0.1*self.UNIT     # Error bars.
        ERROR_WIDTH = 1.5*self.UNIT
        HEIGHT = self._z_height(height)

        if not error_bar:
            return
        assert len(points) == len(error_bar)

        # Convert points and error bar coordinates to pixels.
        xlim, ylim = self._xlim, self._ylim
        x_min, x_max = xlim
        x_len = x_max - x_min
        y_min, y_max = ylim
        y_len = y_max - y_min
        x_scale = self.PLOT_WIDTH / float(x_len)
        y_scale = self.PLOT_HEIGHT / float(y_len)
    
        points_pixel = [None] * len(points)
        error_bar_pixel = [None] * len(error_bar)
        for i in range(len(points)):
            x, y = points[i]
            x_pix = _coord2pixel(x, x_min, x_len, self.PLOT_X, self.PLOT_WIDTH)
            y_pix = _coord2pixel(y, y_min, y_len, self.PLOT_Y,self.PLOT_HEIGHT)
            points_pixel[i] = x_pix, y_pix
            
            x1, x2 = error_bar[i]
            x1_pix, x2_pix = x1*y_scale, x2*y_scale
            error_bar_pixel[i] = x1_pix, x2_pix

        
        w = self.handle.write
        finish = pr.finish(pr.diffuse(0.1), pr.ambient(0.5, 0.5, 0.5))
        w(pr.declare("ERRBAR", pr.cylinder(
            pr.vector(0, 0, 0), pr.vector(1, 0, 0), ERROR_THICKNESS, finish)))

        for i in range(len(points_pixel)):
            pigment = pr.pigment(pr.color(*color[i]))
            
            x, y = points_pixel[i]
            err_l, err_u = error_bar_pixel[i]
            err_total = err_l + err_u
            if not err_total:
                continue
            w(pr.object_(
                "ERRBAR", pr.rotate(0, 0, 90), pr.scale(1, err_total, 1),
                pr.translate(x, y-err_l, -HEIGHT), pigment))
            w(pr.object_(
                "ERRBAR", pr.scale(ERROR_WIDTH, 1, 1),
                pr.translate(x-ERROR_WIDTH/2, y+err_u, -HEIGHT), pigment))
            w(pr.object_(
                "ERRBAR", pr.scale(ERROR_WIDTH, 1, 1),
                pr.translate(x-ERROR_WIDTH/2, y-err_l, -HEIGHT), pigment))

    def draw_onpoint_labels(self, points, label, point_size=1.0, height=0):
        import povray as pr

        RADIUS = self.POINT_RADIUS*point_size      # Radius of each point.
        HEIGHT = self._z_height(height)
        SCALE = RADIUS * 1.25
        THICKNESS = RADIUS/2

        xlim, ylim = self._xlim, self._ylim
        x_min, x_len = xlim[0], xlim[1]-xlim[0]
        y_min, y_len = ylim[0], ylim[1]-ylim[0]

        w = self.handle.write
        
        # lower filter is brighter
        pigment = pr.pigment(pr.color(1, 1, 1, 0.25))
        finish = pr.finish(
            pr.ambient(1, 1, 1),
            pr.refraction(1),
            pr.ior(1.5),
            )

        w("// Draw labels on each point.\n")
        for i in range(len(label)):
            if label[i] is None:
                continue
            
            x_coord, y_coord = points[i]

            x_pixel = _coord2pixel(
                x_coord, x_min, x_len, self.PLOT_X, self.PLOT_WIDTH)
            y_pixel = _coord2pixel(
                y_coord, y_min, y_len, self.PLOT_Y, self.PLOT_HEIGHT)
            z_pixel = -HEIGHT

            # Don't draw points that are off the drawing area.
            if x_pixel<self.PLOT_X or x_pixel >= self.PLOT_X+self.PLOT_WIDTH:
                continue
            if y_pixel<self.PLOT_Y or y_pixel >= self.PLOT_Y+self.PLOT_HEIGHT:
                continue
            
            w(pr.declare("X", pr.text(
                "FONTFILE", label[i], THICKNESS, 0,
                pigment, finish, 
                pr.no_shadow(),
                pr.scale(SCALE, SCALE, SCALE),
                pr.translate(x_pixel, y_pixel, z_pixel-RADIUS),
                )))
            # Center the label.
            w(pr.declare("X", pr.object_(
                "X", pr.translate_by_extent("X", -0.5, -0.5, 0))))
            w(pr.object_("X"))
                
        w("\n")

    def draw_overpoint_labels(
        self, points, error_bar, label, point_size=1.0, height=0):
        import povray as pr

        RADIUS = self.POINT_RADIUS*point_size      # Radius of each point.
        HEIGHT = self._z_height(height)
        THICKNESS = RADIUS/2.0
        SCALE = RADIUS * 1.75    # how big to make the label
        BUFFER = RADIUS * 0.50   # how much space between point and label.
        COLOR = 0.5, 0.5, 0.5

        xlim, ylim = self._xlim, self._ylim
        x_min, x_len = xlim[0], xlim[1]-xlim[0]
        y_min, y_len = ylim[0], ylim[1]-ylim[0]
        y_scale = self.PLOT_HEIGHT / float(y_len)
        
        w = self.handle.write
        
        pigment = pr.pigment(pr.color(*COLOR))
        finish = pr.finish()

        w("// Draw labels over each point.\n")
        for i in range(len(label)):
            if label[i] is None:
                continue

            err_l = err_u = 0
            if error_bar[i]:
                x1, x2 = error_bar[i]
                err_l, err_u = x1*y_scale, x2*y_scale

            x_coord, y_coord = points[i]

            x_pixel = _coord2pixel(
                x_coord, x_min, x_len, self.PLOT_X, self.PLOT_WIDTH)
            y_pixel = _coord2pixel(
                y_coord, y_min, y_len, self.PLOT_Y, self.PLOT_HEIGHT)
            z_pixel = -HEIGHT

            draw_above = y_coord < y_min + y_len/2

            if draw_above:
                # Draw above the point.
                y_pixel += RADIUS + BUFFER + err_u
            else:
                y_pixel -= RADIUS + BUFFER + err_l
            z_pixel -= RADIUS

            # Don't draw points that are off the drawing area.
            if x_pixel<self.PLOT_X or x_pixel >= self.PLOT_X+self.PLOT_WIDTH:
                continue
            if y_pixel<self.PLOT_Y or y_pixel >= self.PLOT_Y+self.PLOT_HEIGHT:
                continue
            
            w(pr.declare("X", pr.text(
                "FONTFILE", label[i], THICKNESS, 0,
                pigment, finish, 
                pr.no_shadow(),
                pr.rotate(0, 0, 90),
                pr.scale(SCALE, SCALE, SCALE),
                pr.translate(x_pixel, y_pixel, z_pixel),
                )))
            # Center the label.
            w(pr.declare("X", pr.object_(
                "X", pr.translate_by_extent("X", 0.5, 0, 0))))
            if not draw_above:
                w(pr.declare("X", pr.object_(
                    "X", pr.translate_by_extent("X", 0, -1, 0))))
            w(pr.object_("X"))
                
        w("\n")
            
    def draw_line(self, points, color, line_size=1.0, height=0):
        # height should be an integer.  0 is lowest.
        import math
        import povray as pr

        THICKNESS = 0.40 * self.UNIT * line_size
        # BUG: HEIGHT should be scaled based on line_size.
        HEIGHT = self._z_height(height)

        xlim, ylim = self._xlim, self._ylim
        x_min, x_len = xlim[0], xlim[1]-xlim[0]
        y_min, y_len = ylim[0], ylim[1]-ylim[0]

        w = self.handle.write
        w("// Plot each segment of the line.\n")
        finish = pr.finish(
            # Metallic finish.
            pr.diffuse(0.4),
            pr.ambient(0.5, 0.5, 0.5),
            pr.phong(0.3),              # bigger means brighter shiny spot
            pr.phong_size(10),          # bigger means smaller shiny spot
            # Clean and simple.
            #pr.diffuse(0.3),             # lower numbers are darker
            #pr.ambient(0.8, 0.8, 0.8),   # lower numbers are darker
            #pr.brilliance(5),            # higher numbers look bolder
            # Extra stuff.
            #pr.reflection(0.8),
            #pr.specular(0.3),           # Makes things shinier.
            #pr.roughness(0.05),
            )
        w(pr.declare("LINESEG", pr.cylinder(
            pr.vector(0, 0, -HEIGHT), pr.vector(1, 0, -HEIGHT), THICKNESS,
            finish,
            #pr.no_shadow(),
            )))
        w(pr.declare("LINEEND", pr.sphere(
            pr.vector(0, 0, -HEIGHT), THICKNESS, finish,
            #pr.no_shadow(),
            )))

        # For convenience.
        PLOT_X, PLOT_Y = self.PLOT_X, self.PLOT_Y
        PLOT_WIDTH, PLOT_HEIGHT = self.PLOT_WIDTH, self.PLOT_HEIGHT
                
        for i in range(len(points)-1):
            xc1, yc1 = points[i]
            xc2, yc2 = points[i+1]

            xp1 = _coord2pixel(xc1, x_min, x_len, PLOT_X, PLOT_WIDTH)
            yp1 = _coord2pixel(yc1, y_min, y_len, PLOT_Y, PLOT_HEIGHT)
            xp2 = _coord2pixel(xc2, x_min, x_len, PLOT_X, PLOT_WIDTH)
            yp2 = _coord2pixel(yc2, y_min, y_len, PLOT_Y, PLOT_HEIGHT)

            # Make sure all the points are within the drawing area.
            # Both X-coordinates are off the plot.
            if max(xp1, xp2) < PLOT_X or min(xp1, xp2) >= PLOT_X+PLOT_WIDTH:
                continue
            # Both Y-coordinates are off the plot.
            if max(yp1, yp2) < PLOT_Y or min(yp1, yp2) >= PLOT_Y + PLOT_HEIGHT:
                continue
            # At this point, at least one of the points are within the
            # drawing area.
            EPS = 1E-8
            if abs(xp2-xp1) < EPS and abs(yp2-yp1) < EPS:
                # Just a point.  Ignore.
                continue
            elif abs(xp2-xp1) < EPS:
                # Slope is infinity, straight line up and down.
                if xp1 < PLOT_X or xp1 > PLOT_X+PLOT_WIDTH:
                    continue
                yp1 = min(max(yp1, PLOT_Y), PLOT_Y+PLOT_HEIGHT)
                yp2 = min(max(yp2, PLOT_Y), PLOT_Y+PLOT_HEIGHT)
            elif abs(yp2-yp1) < EPS:
                # Slope is 0, straight line left and right.
                if yp1 < PLOT_Y or yp1 > PLOT_Y+PLOT_HEIGHT:
                    continue
                xp1 = min(max(xp1, PLOT_X), PLOT_X+PLOT_WIDTH)
                xp2 = min(max(xp2, PLOT_X), PLOT_X+PLOT_WIDTH)
            else:
                m = float((yp2-yp1)) / (xp2-xp1)
                b = yp1 - m * xp1
                # Make sure the x-coords are in the plot.
                min_x, max_x = PLOT_X, PLOT_X+PLOT_WIDTH
                min_y, max_y = m*min_x + b, m*max_x + b
                if xp1 < min_x:
                    xp1, yp1 = min_x, min_y
                if xp1 > max_x:
                    xp1, yp1 = max_x, max_y
                if xp2 < min_x:
                    xp2, yp2 = min_x, min_y
                if xp2 > max_x:
                    xp2, yp2 = max_x, max_y
                # Make sure the y-coords are in the plot.
                min_y, max_y = PLOT_Y, PLOT_Y+PLOT_HEIGHT
                min_x, max_x = (min_y-b)/m, (max_y-b)/m
                if yp1 < min_y:
                    xp1, yp1 = min_x, min_y
                if yp1 > max_y:
                    xp1, yp1 = max_x, max_y
                if yp2 < min_y:
                    xp2, yp2 = min_x, min_y
                if yp2 > max_y:
                    xp2, yp2 = max_x, max_y
                
            # Draw a line.
            pigment = pr.pigment(pr.color(*color))
            linelen = math.sqrt((xp2-xp1)**2 + (yp2-yp1)**2)
            angle = 90
            if abs(xp2-xp1) > 1E-8:
                angle = math.degrees(math.atan(float(yp2-yp1)/(xp2-xp1)))
            w(pr.object_(
                "LINESEG", pr.scale(linelen, 1, 1), pr.rotate(0, 0, angle),
                pr.translate(xp1, yp1, 0), pigment))
            w(pr.object_(
                "LINEEND", pr.translate(xp2, yp2, 0),
                #pr.pigment(pr.color(1, 0, 1)),
                pigment,
                ))
            if i == 0:
                w(pr.object_(
                    "LINEEND", pr.translate(xp1, yp1, 0),
                    #pr.pigment(pr.color(1, 0, 1)),
                    pigment,
                    ))
        w("\n")

    def draw_bars(self, bars, barwidth, color, width_size=1.0):
        # height should be an integer.  0 is lowest.
        import math
        import povray as pr

        THICKNESS = 0.6 * self.UNIT * 0.1

        xlim, ylim = self._xlim, self._ylim
        x_min, x_len = xlim[0], xlim[1]-xlim[0]
        y_min, y_len = ylim[0], ylim[1]-ylim[0]

        w = self.handle.write
        w("// Plot each bar.\n")
        finish = pr.finish(
            # Rough, textured finish.
            pr.diffuse(0.4),
            pr.crand(0.2),
            pr.ambient(0.5, 0.5, 0.5),
            )
        w(pr.declare("BAR", pr.box(
            pr.vector(0, 0, -THICKNESS), pr.vector(1, 1, -THICKNESS), 
            finish,
            )))

        # For convenience.
        PLOT_X, PLOT_Y = self.PLOT_X, self.PLOT_Y
        PLOT_WIDTH, PLOT_HEIGHT = self.PLOT_WIDTH, self.PLOT_HEIGHT

        for i in range(len(bars)):
            if len(bars) == 3:
                xc, yminc, ymaxc = bars[i]
            else:
                xc, ymaxc = bars[i]
                yminc = 0.0
            assert yminc <= ymaxc
                
            xc -= barwidth/2.0

            xp = _coord2pixel(xc, x_min, x_len, self.PLOT_X, self.PLOT_WIDTH)
            yp = _coord2pixel(
                yminc, y_min, y_len, self.PLOT_Y, self.PLOT_HEIGHT)
            ymaxp = _coord2pixel(
                ymaxc, y_min, y_len, self.PLOT_Y, self.PLOT_HEIGHT)
            width = barwidth / float(x_len) * self.PLOT_WIDTH
            height = ymaxp - yp

            width = width * width_size

            # Make sure the points are within the drawing area.
            if xp < PLOT_X:
                width = width - (PLOT_X-xp)
                xp = PLOT_X
            if xp+width > PLOT_X+PLOT_WIDTH:
                width = width - ((xp+width)-(PLOT_X+PLOT_WIDTH))
            if yp < PLOT_Y:
                height = height - (PLOT_Y-yp)
                yp = PLOT_Y
            if yp+height > PLOT_Y+PLOT_HEIGHT:
                height = height - ((yp+height)-(PLOT_Y+PLOT_HEIGHT))
            if width <= 0 or height <= 0:
                continue
            
            #if max(xp1, xp2) < self.PLOT_X or \
            #       min(xp1, xp2) >= self.PLOT_X + self.PLOT_WIDTH:
            #    continue
                
            # Draw the box.
            pigment = pr.pigment(pr.color(*color))
            w(pr.object_(
                "BAR", pr.scale(width, height, 1), pr.translate(xp, yp, 0),
                pigment))
        w("\n")

    def draw(self):
        return self.handle.getvalue()

    def _z_height(self, height, scale=1.0):
        HEIGHT_UNIT = 0.50 * self.UNIT * scale
        z = (0.5+height)*HEIGHT_UNIT
        return z


def _make_graph(
    X_min, X_max, Y_min, Y_max, xlim=None, ylim=None, 
    xtick=None, xtick_label=None, ytick=None, ytick_label=None, tick_size=1.0,
    xlabel=None, ylabel=None, label_size=1.0, font=None,
    width=None, height=None, **keywds):
    # xlim         Tuple of (minimum x, maximum x).
    # xtick        Coordinates of the tickmarks on the x-axis (or True).
    # xtick_label  Labels for the tickmarks.  Should be parallel to xtick.
    # tick_size    Scales the size of the tick mark labels.
    #
    # xlabel       Label for the X-axis.
    # ylabel
    # label_size   Scales the size of the labels.
    # font
    # 
    # width        Number of pixels wide for the plot.
    # height
    from StringIO import StringIO

    # Check the inputs
    x = _set_default_axes(
        X_min, X_max, Y_min, Y_max,
        xlim, ylim, xtick, ytick, xtick_label, ytick_label)
    xlim, ylim, xtick, ytick, xtick_label, ytick_label = x

    handle = StringIO()
    graph = Graph(handle, xlim, ylim, width=width, height=height)
    graph.declare_fontfile(fontfile=font)
    graph.position_camera()
    graph.draw_background()
    graph.draw_axes()
    graph.draw_tick_marks(xtick, ytick, tick_size=tick_size)
    graph.draw_tick_labels(xtick_label, ytick_label, label_size=label_size)
    graph.draw_labels(xlabel, ylabel, label_size=label_size)
    return graph

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

def scatter(
    X, Y, color=None, shape=None, onpoint_label=None, overpoint_label=None, 
    error_bar=None, point_size=1.0, graph=None, **keywds):
    """Return the povray-formatted script to generate a scatter plot.
    
    X, Y             Parallel lists of the coordinates of the points.
    color            List of the color for each point (<r>,<g>,<b>) from 0-1.
    shape            List of DEFAULT, CIRCLE, SQUARE, or DIAMOND.
    onpoint_label    What label to put on each point.
    overpoint_label  What label to put over each point.
    error_bar        List of sizes of the error bar for each point.
    point_size       Scales the size of each point.

    More from _make_graph.
    
    """
    from StringIO import StringIO

    # Check the inputs
    color = _set_default_color(color, len(X))
    shape = _set_default_shape(shape, len(X))
    onpoint_label = _set_default_label(onpoint_label, len(X))
    overpoint_label = _set_default_label(overpoint_label, len(X))
    error_bar = _set_default_error_bar(error_bar)
    assert len(X) == len(Y)
    assert len(color) == len(X)
    assert len(shape) == len(X)
    assert len(onpoint_label) == len(X)
    assert len(overpoint_label) == len(X)
    assert not error_bar or len(error_bar) == len(X)
    points = zip(X, Y)

    if graph is None:
        graph = _make_graph(min(X), max(X), min(Y), max(Y), **keywds)

    graph.draw_points(points, color, shape, point_size=point_size)
    graph.draw_onpoint_labels(points, onpoint_label, point_size=point_size)
    graph.draw_overpoint_labels(
        points, error_bar, overpoint_label, point_size=point_size)
    graph.draw_error_bars(points, error_bar, color)
    return graph

def line(*args, **keywds):
    """Return the povray-formatted script to generate a line graph.
    
    line1        List of the (x, y) coordinates of the points.
    line2        ...
    color        List of colors, 1 for each line (<r>,<g>,<b>) from 0-1.
    shape        List of DEFAULT, CIRCLE, SQUARE, or DIAMOND.
    line_size    Scales the size of the line.
    draw_points  Whether to draw the points or not.
    point_size   Scales the size of each point.
    
    """
    import itertools

    #KNOWN_KEYWDS = [
    #    "color", "shape", "draw_points", "point_size",
    #    "xlim", "ylim",
    #    "xtick", "xtick_label", "ytick", "ytick_label", "tick_size",
    #    "xlabel", "ylabel", "label_size", "font",
    #    "width", "height"
    #    ]
    #for word in keywds:
    #    assert word in KNOWN_KEYWDS, "I don't know the keyword: %s" % word
    assert args, "No lines given."

    color = keywds.get("color", None)
    shape = keywds.get("shape", None)
    line_size = keywds.get("line_size", 1.0)
    draw_points = keywds.get("draw_points", False)
    point_size = keywds.get("point_size", 1.0)
    graph = keywds.get("graph", None)
    
    # Check the inputs
    lines = args
    for line in lines:
        for x in line:
            assert len(x) == 2
    X = list(itertools.chain.from_iterable(
        [[x[0] for x in line] for line in lines]))
    Y = list(itertools.chain.from_iterable(
        [[x[1] for x in line] for line in lines]))
    color = _set_default_color(color, len(lines))
    shape = _set_default_shape(shape, len(lines))
    assert len(color) == len(lines)
    assert len(shape) == len(lines)

    if not graph:
        graph = _make_graph(min(X), max(X), min(Y), max(Y), **keywds)
        
    for i, line in enumerate(lines):
        graph.draw_line(line, color[i], line_size=line_size, height=i)
        if not draw_points:
            continue
        # Should let the user specify which points to plot.
        col = [color[i]] * len(line)
        shp = [shape[i]] * len(line)
        graph.draw_points(line, col, shp, None, height=i)

    return graph

def bar(*args, **keywds):
    """Return the povray-formatted script to generate a bar graph.

    series1      List of the (x, y) coordinates of this set of bars.
                 Can also be (x, y_min, y_max).
    series2      ...
    color        List of colors, 1 for each bar (<r>,<g>,<b>) from 0-1.
    width_size   Width of bar.  1.0 means bars do not overlap.

    """
    import itertools
    
    # TODO: ERROR BAR
    color = keywds.get("color", None)
    width_size = keywds.get("width_size", 1.0)
    graph = keywds.get("graph", None)

    # Check the inputs
    assert args, "No bars given."
    all_series = args
    X, Y = [], []
    for series in all_series:
        for bar in series:
            assert len(bar) in [2, 3]
            X.append(bar[0])
            Y.append(bar[1])
            if len(bar) == 3:
                Y.append(bar[2])
            else:
                Y.append(0.0)
    color = _set_default_color(color, len(all_series))
    assert len(color) == len(all_series)

    barwidth = _set_bar_width(X)

    if graph is None:
        X_min = min(X) - barwidth/2.0
        X_max = max(X) + barwidth/2.0
        graph = _make_graph(X_min, X_max, min(Y), max(Y), **keywds)
        
    for i, series in enumerate(all_series):
        graph.draw_bars(series, barwidth, color[i], width_size=width_size)

    return graph

def _set_bar_width(X):
    # Width should be the minimum distance from one bar to the next.
    # Can draw the bars with no overlap.
    X = sorted(X)
    dists = [X[i]-X[i-1] for i in range(1, len(X))]
    width = min(dists)
    return width

def _set_default_tick(tick, coord_min, coord_max):
    import operator
    import plotlib
    
    if not tick:
        tick = []
    elif operator.isSequenceType(tick):
        pass
    elif tick == True:
        tick = plotlib.place_ticks(coord_min, coord_max, num_ticks=5)
    elif type(tick) is type(0):
        tick = plotlib.place_ticks(coord_min, coord_max, num_ticks=tick)
    else:
        raise AssertionError, "I don't understand: %s" % str(tick)
    return tick

def _set_default_tick_label(label, tick):
    import operator
    if not tick:
        label = []
    elif operator.isSequenceType(label):
        assert len(label) == len(tick)
    elif label:
        # Figure out the number of digits after the decimal place.
        x = [max(len(str(abs(x)%1))-2, 0) for x in tick]
        digits = max(x)
        label = ["%.*f" % (digits, x) for x in tick]
    else:
        label = []
    return label

def _set_default_axes(
    min_x, max_x, min_y, max_y,
    xlim, ylim, xtick, ytick, xtick_label, ytick_label):
    if xtick_label and not xtick:
        xtick = True
    if ytick_label and not ytick:
        ytick = True
    if xlim:
        assert len(xlim) == 2
        xtick = _set_default_tick(xtick, xlim[0], xlim[1])
    else:
        xtick = _set_default_tick(xtick, min_x, max_x)
    if ylim:
        assert len(ylim) == 2
        ytick = _set_default_tick(ytick, ylim[0], ylim[1])
    else:
        ytick = _set_default_tick(ytick, min_y, max_y)
    xtick_label = _set_default_tick_label(xtick_label, xtick)
    ytick_label = _set_default_tick_label(ytick_label, ytick)
    if xtick:
        min_x, max_x = min(min_x, min(xtick)), max(max_x, max(xtick))
    if ytick:
        min_y, max_y = min(min_y, min(ytick)), max(max_y, max(ytick))
    xlim = _set_default_lim(xlim, min_x, max_x)
    ylim = _set_default_lim(ylim, min_y, max_y)
    assert not xtick_label or len(xtick) == len(xtick_label)
    assert not ytick_label or len(ytick) == len(ytick_label)
    return xlim, ylim, xtick, ytick, xtick_label, ytick_label

def _set_default_lim(lim, coord_min, coord_max):
    # Return (min, max).
    if lim:
        assert len(lim) == 2
        return lim
    coord_min, coord_max = float(coord_min), float(coord_max)
    x_len = coord_max - coord_min
    x_min = coord_min - x_len*0.05
    x_max = coord_max + x_len*0.05
    return x_min, x_max

def _set_default_color(color, n):
    # Return a list of colors.
    import operator
    import colorlib

    assert n > 0
    if color == True:
        if n == 1:
            color = (0, 0, 0)
        else:
            color = colorlib.bild_colors(n)
    elif type(color) is type(_set_default_color):
        color = color(n)
    if not color:
        color = (0, 0, 0)
    assert operator.isSequenceType(color)
    # Distinguish between (0, 0, 0) and [(0, 0, 0), (0, 0, 0), (0, 0, 0)]
    if not operator.isSequenceType(color[0]):
        color = [color] * n
    for c in color:
        assert len(c) == 3, "Does not look like color: %s" % str(c)
    assert len(color) == n, "You specified %d colors, but I expected %d." % (
        len(color), n)
    return color

def _set_default_shape(shape, n):
    # Return a list of shapes.
    import operator
    if shape is None:
        shape = DEFAULT
    if type(shape) is type(DEFAULT):
        shape = [shape]
    assert operator.isSequenceType(shape)
    if len(shape) == 1:
        shape = shape * n
    assert len(shape) == n
    return shape

def _set_default_label(label, n):
    # Return a list of labels.  label can be:
    # 1.  List of strings.  None means don't label this point.
    # 2.  None
    # 3.  []
    import operator
    if label is None or label == []:
        label = [None]
    assert operator.isSequenceType(label)
    if len(label) == 1:
        label = label * n
    assert len(label) == n
    return label

def _set_default_error_bar(error_bar):
    import operator
    
    if not error_bar:
        return None
    # Each error can be a single number or a tuple of (lower, upper).
    # Normalize so each of them are tuples of (lower, upper).
    error_bar = error_bar[:]
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
    return error_bar
    
def _coord2pixel(x_coord, coord_min, coord_length, pixel_min, pixel_length):
    # The user provides data in coordinates of their own scale.
    # Convert those coordinates to pixels on the plot.
    scale = float(pixel_length) / coord_length
    x_pixel = (x_coord-coord_min)*scale + pixel_min
    return x_pixel
